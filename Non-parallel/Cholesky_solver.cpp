#include "Cholesky_solver.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <unordered_map>

std::vector<double> CholeskySolver::decompose(const Matrix& A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("Matrix must be square");
    }
    if (!A.is_symmetric()) {
        throw std::invalid_argument("Matrix must be symmetric");
    }

    const int n = A.rows();
    const auto& row_ptrs = A.row_pointers();
    const auto& col_inds = A.col_indices();
    const auto& values = A.values();

    std::vector<double> L_vals(values); // The L matrix non-xero values

    // Easing access to matrix by columns
    std::vector<std::unordered_map<int, double>> col_map(n);

    for (int i = 0; i < n; ++i) {
        for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; ++idx) {
            int j = col_inds[idx];
            if (i >= j) { 
                col_map[j][i] = L_vals[idx];
            }
        }
    }

    for (int j = 0; j < n; ++j) {
        // Processing diagonal elements of the matrix
        double diag_val = col_map[j][j];
        for (int k = 0; k < j; ++k) {
            if (col_map[k].find(j) != col_map[k].end()) {
                diag_val -= col_map[k][j] * col_map[k][j];
            }
        }
        // Checking suitability for the Cholesky decomposition
        if (diag_val <= 0) {
            throw std::runtime_error("Matrix is not positive definite");
        }

        diag_val = std::sqrt(diag_val);
        col_map[j][j] = diag_val;

        // Processing the rest of the matrix
        for (int i = j + 1; i < n; ++i) {
            if (col_map[j].find(i) != col_map[j].end()) {
                double val = col_map[j][i];
                for (int k = 0; k < j; ++k) {
                    if (col_map[k].find(i) != col_map[k].end() &&
                        col_map[k].find(j) != col_map[k].end()) {
                        val -= col_map[k][i] * col_map[k][j];
                    }
                }
                col_map[j][i] = val / diag_val;
            }
        }
    }

    std::vector<double> result(values.size(), 0.0);
    for (int i = 0; i < n; ++i) {
        for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; ++idx) {
            int j = col_inds[idx];
            if (i >= j) {
                result[idx] = col_map[j][i];
            }
        }
    }

    return result;
}

std::vector<double> CholeskySolver::solve(
    const Matrix& A,
    const std::vector<double>& b,
    const std::vector<double>& L_vals
) {
    const int n = A.rows();
    const auto& row_ptrs = A.row_pointers();
    const auto& col_inds = A.col_indices();

    // Solution of the Ly = b equation
    std::vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = b[i];
        for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; ++idx) {
            int j = col_inds[idx];
            if (j < i) { 
                sum -= L_vals[idx] * y[j];
            }
        }
        for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; ++idx) {
            if (col_inds[idx] == i) {
                y[i] = sum / L_vals[idx];
                break;
            }
        }
    }

    // Inverse substitution (L^Tx = y)
    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double sum = y[i];
        for (int j = i + 1; j < n; ++j) {
            for (int idx = row_ptrs[j]; idx < row_ptrs[j + 1]; ++idx) {
                if (col_inds[idx] == i) {
                    sum -= L_vals[idx] * x[j];
                    break;
                }
            }
        }
        for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; ++idx) {
            if (col_inds[idx] == i) {
                x[i] = sum / L_vals[idx];
                break;
            }
        }
    }

    return x;
}

std::pair<long, long> CholeskySolver::test(const Matrix& A, const std::vector<double>& b) {

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> L_vals = CholeskySolver::decompose(A);
    auto end = std::chrono::high_resolution_clock::now();
    auto time_decompose = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    std::vector<double> x = CholeskySolver::solve(A, b, L_vals);
    end = std::chrono::high_resolution_clock::now();
    auto time_solve = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    return std::make_pair(time_decompose, time_solve);
}