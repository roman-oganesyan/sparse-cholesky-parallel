#include "CSR_Matrix.h"
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <algorithm>

Matrix load_mtx(const std::string& filename) {
    std::ifstream file(filename);
    if (!file)
        throw std::runtime_error("Couldn't open file " + filename);

    std::string line;
    bool symmetric = false;
    bool pattern = false;

    while (std::getline(file, line)) {
        if (line[0] != '%') break;
        else {
            if (line.find("symmetric") != std::string::npos) symmetric = true;
            if (line.find("pattern") != std::string::npos) pattern = true;
        }
    }

    std::istringstream iss(line);
    int rows, cols, nnz;
    iss >> rows >> cols >> nnz;

    std::vector<Matrix::Entry> entries;
    entries.reserve(symmetric ? nnz * 2 : nnz);

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '%') continue;

        std::istringstream line_stream(line);
        int row, col;
        double value = 1.0;

        if (pattern) {
            line_stream >> row >> col;
        }
        else {
            line_stream >> row >> col >> value;
        }

        if (row > rows || col > cols)
            throw std::out_of_range("Invalid file data: index out of bounds");

        row--; col--;
        entries.emplace_back(row, col, value);
        if (symmetric && row != col) {
            entries.emplace_back(col, row, value);
        }
    }

    std::sort(entries.begin(), entries.end());
    return Matrix(rows, cols, entries);
}

Matrix::Matrix(int rows, int cols, int nnz)
    : m_rows(rows), m_cols(cols), m_nnz(nnz) {
    if (rows <= 0 || cols <= 0 || nnz < 0)
        throw std::runtime_error("Invalid matrix dimensions");
    m_elems = new double[nnz];
    m_col_inds = new int[nnz];
    m_row_ptrs = new int[rows + 1];
}

Matrix::Matrix(int rows, int cols, const std::vector<Entry>& entries)
    : m_rows(rows), m_cols(cols), m_nnz(entries.size()) {
    if (rows <= 0 || cols <= 0)
        throw std::runtime_error("Invalid matrix dimensions");
    if (entries.empty())
        throw std::runtime_error("Matrix must have at least one non-zero element");

    m_elems = new double[m_nnz];
    m_col_inds = new int[m_nnz];
    m_row_ptrs = new int[rows + 1]();

    for (const auto& entry : entries) {
        int row = std::get<0>(entry);
        m_row_ptrs[row + 1]++;
    }

    for (int i = 0; i < rows; i++) {
        m_row_ptrs[i + 1] += m_row_ptrs[i];
    }

    std::vector<int> offsets(rows, 0);
    for (const auto& entry : entries) {
        int row = std::get<0>(entry);
        int col = std::get<1>(entry);
        double val = std::get<2>(entry);
        int idx = m_row_ptrs[row] + offsets[row];
        m_elems[idx] = val;
        m_col_inds[idx] = col;
        offsets[row]++;
    }
}

Matrix::~Matrix() {
    delete[] m_elems;
    delete[] m_col_inds;
    delete[] m_row_ptrs;
}

Matrix::Matrix(Matrix&& other) noexcept
    : m_rows(other.m_rows), m_cols(other.m_cols), m_nnz(other.m_nnz),
    m_elems(other.m_elems), m_col_inds(other.m_col_inds), m_row_ptrs(other.m_row_ptrs) {
    other.m_rows = 0;
    other.m_cols = 0;
    other.m_nnz = 0;
    other.m_elems = nullptr;
    other.m_col_inds = nullptr;
    other.m_row_ptrs = nullptr;
}

Matrix& Matrix::operator=(Matrix&& other) noexcept {
    if (this != &other) {
        delete[] m_elems;
        delete[] m_col_inds;
        delete[] m_row_ptrs;

        m_rows = other.m_rows;
        m_cols = other.m_cols;
        m_nnz = other.m_nnz;
        m_elems = other.m_elems;
        m_col_inds = other.m_col_inds;
        m_row_ptrs = other.m_row_ptrs;

        other.m_rows = 0;
        other.m_cols = 0;
        other.m_nnz = 0;
        other.m_elems = nullptr;
        other.m_col_inds = nullptr;
        other.m_row_ptrs = nullptr;
    }
    return *this;
}

int Matrix::rows() const { return m_rows; }
int Matrix::cols() const { return m_cols; }
int Matrix::nnz() const { return m_nnz; }
const double* Matrix::values() const { return m_elems; }
const int* Matrix::col_indices() const { return m_col_inds; }
const int* Matrix::row_pointers() const { return m_row_ptrs; }

bool Matrix::is_symmetric() const {
    if (m_rows != m_cols) return false;
    for (int i = 0; i < m_rows; ++i) {
        for (int idx = m_row_ptrs[i]; idx < m_row_ptrs[i + 1]; ++idx) {
            int j = m_col_inds[idx];
            double val_ij = m_elems[idx];
            double val_ji = get_element(j, i);
            if (std::abs(val_ij - val_ji) > 1e-9) return false;
        }
    }
    return true;
}

double Matrix::get_element(int i, int j) const {
    if (i < 0 || i >= m_rows || j < 0 || j >= m_cols) return 0.0;
    for (int idx = m_row_ptrs[i]; idx < m_row_ptrs[i + 1]; ++idx) {
        if (m_col_inds[idx] == j) return m_elems[idx];
    }
    return 0.0;
}

std::vector<double> Matrix::multiply_by_vector(const std::vector<double>& vec) const {
    if (static_cast<int>(vec.size()) != m_cols) {
        throw std::invalid_argument("Vector size must match matrix columns");
    }
    std::vector<double> result(m_rows, 0.0);
    for (int i = 0; i < m_rows; ++i) {
        double sum = 0.0;
        for (int idx = m_row_ptrs[i]; idx < m_row_ptrs[i + 1]; ++idx) {
            int col = m_col_inds[idx];
            sum += m_elems[idx] * vec[col];
        }
        result[i] = sum;
    }
    return result;
}

std::vector<double> Matrix::transpose_multiply_by_vector(const std::vector<double>& vec) const {
    if (static_cast<int>(vec.size()) != m_rows) {
        throw std::invalid_argument("Vector size must match matrix rows");
    }
    std::vector<double> result(m_cols, 0.0);
    for (int i = 0; i < m_rows; ++i) {
        for (int idx = m_row_ptrs[i]; idx < m_row_ptrs[i + 1]; ++idx) {
            int col = m_col_inds[idx];
            result[col] += m_elems[idx] * vec[i];
        }
    }
    return result;
}