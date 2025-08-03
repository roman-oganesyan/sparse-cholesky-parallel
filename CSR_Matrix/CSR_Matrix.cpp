#include "CSR_Matrix.h"
#include <stdexcept>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <algorithm>


// .mtx CSR matrix loader
Matrix mtx_loader(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) 
        throw std::runtime_error("Couldn't open file " + filename);

    std::string line;
    bool symmetric = false;
    bool pattern = false;

    // Reading the header of the .mtx file
    while (std::getline(file, line)) {
        if (line[0] != '%') break;
        else {
            if (line.find("symmetric") != std::string::npos) symmetric = true;
            if (line.find("pattern") != std::string::npos) pattern = true;
        }
    }

    // Getting the numbers of rows, columns and non-zeros in the matrix
    std::istringstream iss(line);
    int rows, cols, nnz;
    iss >> rows >> cols >> nnz;

    // Getting the non-zero elements
    std::vector<Matrix::Entry> entries;
    entries.reserve(symmetric ? nnz * 2 : nnz);

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '%') continue;

        std::istringstream line_stream(line);
        int row, col;
        double value = 1.0; // For pattern matrices

        if (pattern) {
            line_stream >> row >> col;
        }
        else {
            line_stream >> row >> col >> value;
        }

        if (row > rows)
            throw std::out_of_range("Invalid file data: row index out of bounds");
        if (col > cols)
            throw std::out_of_range("Invalid file data: column index out of bounds");

        row--; col--; // converting to the 0-indexation
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
	if (rows <= 0 || cols <= 0 || nnz <= 0)
		throw std::runtime_error("Invalid matrix dimensions");
	m_elems.reserve(nnz);
	m_col_inds.reserve(nnz);
	m_row_ptrs.resize(rows + 1, 0);
}

Matrix::Matrix(int rows, int cols, const std::vector<Entry>& entries) :
	m_rows(rows), m_cols(cols), m_nnz(entries.size()), m_row_ptrs(rows + 1, 0) {
	if (rows <= 0 || cols <= 0) 
		throw std::runtime_error("Invalid matrix dimensions");
    if (entries.size() == 0)
        throw std::runtime_error("Matrix must have at least one non-zero element");

    m_elems.reserve(m_nnz);
    m_col_inds.reserve(m_nnz);
    int counter = 0;
    m_row_ptrs[0] = 0;

    for (const auto& entry : entries) {
        const int row = std::get<0>(entry);
        const int col = std::get<1>(entry);
        const double val = std::get<2>(entry);
        m_row_ptrs[row + 1]++; // row ptr calculation 1 (number of previous row elements)
        m_col_inds.push_back(col);
        m_elems.push_back(val);
    }

    for (int r = 0; r < rows; ++r) { //row ptr calculation 2 (adding the previous row starting index)
        m_row_ptrs[r + 1] += m_row_ptrs[r];
    }
}

std::vector<double> Matrix::multiply_by_vector(const std::vector<double>& vec) const {
    if (static_cast<int>(vec.size()) != m_cols) {
        throw std::invalid_argument("Vector size must match matrix columns");
    }

    std::vector<double> result(m_rows, 0.0);

    for (int i = 0; i < m_rows; ++i) {
        const int row_start = m_row_ptrs[i];
        const int row_end = m_row_ptrs[i + 1];
        double sum = 0.0;

        for (int j = row_start; j < row_end; ++j) {
            const int col_index = m_col_inds[j];
            sum += m_elems[j] * vec[col_index];
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
        const int row_start = m_row_ptrs[i];
        const int row_end = m_row_ptrs[i + 1];

        for (int j = row_start; j < row_end; ++j) {
            const int col_index = m_col_inds[j];
            result[col_index] += m_elems[j] * vec[i];
        }
    }

    return result;
}

/*
std::vector<Matrix::Entry> CSR_to_mtx(Matrix& csr) {
    
}
*/