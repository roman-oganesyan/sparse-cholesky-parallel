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
	//Matrix(int rows, int cols, const std::vector<int>& col_ptrs, const std::vector<int>& row_ptrs, const std::vector<double>& elems);

	std::vector<double> multiply_by_vector(const std::vector<double>& vec) const;
	std::vector<double> transpose_multiply_by_vector(const std::vector<double>& vec) const;

	int rows() const { return m_rows; }
	int cols() const { return m_cols; }
	int nnz() const { return m_nnz; }
};

Matrix mtx_loader(const std::string& filename);