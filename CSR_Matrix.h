#pragma once
#include <vector>
#include <tuple>
#include <string>

// CSR Matrix class 
class Matrix {
private:
	int m_rows;
	int m_cols;
	int m_nnz;
	std::vector<double> m_elems;
	std::vector<int> m_col_inds;
	std::vector<int> m_row_ptrs;
public:
	using Entry = std::tuple<int, int, double>;
	Matrix() = delete;
	Matrix(int rows, int cols, int nnz);
	Matrix(int rows, int cols, const std::vector<Entry>& entries);

	int rows() const;
	int cols() const;
	int nnz() const;
	const std::vector<double>& values() const;
	const std::vector<int>& col_indices() const;
	const std::vector<int>& row_pointers() const;

	bool is_symmetric() const;
	double get_element(int i, int j) const;
	std::vector<double> multiply_by_vector(const std::vector<double>& vec) const;
	std::vector<double> transpose_multiply_by_vector(const std::vector<double>& vec) const;
};

Matrix load_mtx(const std::string& filename);