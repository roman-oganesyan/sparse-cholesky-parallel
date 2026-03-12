#include "Cholesky_solver.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <unordered_map>
#include <omp.h>

std::vector<double> CholeskySolver::decompose(const Matrix& A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("Matrix must be square");
    }
    if (!A.is_symmetric()) {
        throw std::invalid_argument("Matrix must be symmetric");
    }

    const int n = A.rows();
    const int* row_ptrs = A.row_pointers();
    const int* col_inds = A.col_indices();
    const double* values = A.values();
    const int nnz = A.nnz();

    std::vector<double> L_vals(values, values + nnz);
    std::vector<std::unordered_map<int, double>> col_map(n);

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; ++idx) {
            int j = col_inds[idx];
            if (i >= j) {
#pragma omp critical
                col_map[j][i] = L_vals[idx];
            }
        }
    }

    for (int j = 0; j < n; ++j) {
        double diag_val = col_map[j][j];
#pragma omp parallel for
        for (int k = 0; k < j; ++k) {
            if (col_map[k].count(j)) {
                diag_val -= col_map[k][j] * col_map[k][j];
            }
        }
        if (diag_val <= 0) {
            throw std::runtime_error("Matrix is not positive definite");
        }
        diag_val = std::sqrt(diag_val);
        col_map[j][j] = diag_val;

#pragma omp parallel for schedule(dynamic)
        for (int i = j + 1; i < n; ++i) {
            if (col_map[j].count(i)) {
                double val = col_map[j][i];
                for (int k = 0; k < j; ++k) {
                    if (col_map[k].count(i) && col_map[k].count(j)) {
                        val -= col_map[k][i] * col_map[k][j];
                    }
                }
                col_map[j][i] = val / diag_val;
            }
        }
    }

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; ++idx) {
            int j = col_inds[idx];
            if (i >= j) {
                L_vals[idx] = col_map[j][i];
            }
        }
    }

    return L_vals;
}

std::vector<double> CholeskySolver::solve(
    const Matrix& A,
    const std::vector<double>& b,
    const std::vector<double>& L_vals
) {
    const int n = A.rows();
    const int* row_ptrs = A.row_pointers();
    const int* col_inds = A.col_indices();

    // Прямая подстановка (Ly = b)
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

    // Обратная подстановка (L^Tx = y)
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

    // Precision check
    std::vector<double> Ax = A.multiply_by_vector(x);
    double discrepancy = 0.0;
    int error_count = 0;
    const double tolerance = 1e-6;

    double b_norm = 0.0;
    for (int i = 0; i < n; ++i) {
        b_norm += b[i] * b[i];
    }
    b_norm = std::sqrt(b_norm);

    double Axb_norm = 0.0;
    for (int i = 0; i < n; ++i) {
        Axb_norm += (Ax[i] - b[i]) * (Ax[i] - b[i]);
    }
    Axb_norm = std::sqrt(Axb_norm);

    // dim > 10^5 -> /b_norm

    std::cout << "Norm difference: " << Axb_norm / b_norm << std::endl;
    std::cout << "First 100 elements:" << std::endl;
    for (int i = 0; i < 100; ++i) {
        std::cout << i << ' ' << x[i] << std::endl;
    }
    return x;
}

std::pair<long, long> CholeskySolver::test(int thread_num, const Matrix& A, const std::vector<double>& b) {
    omp_set_num_threads(thread_num);
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> L_vals = decompose(A);
    auto end = std::chrono::high_resolution_clock::now();
    auto time_decompose = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    std::vector<double> x = solve(A, b, L_vals);
    end = std::chrono::high_resolution_clock::now();
    auto time_solve = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    return std::make_pair(time_decompose, time_solve);
}