#include "CSR_Matrix.h"
#include "Cholesky_solver.h"
#include <iostream>
#include <chrono>

int main() {
    Matrix A = load_mtx("5k.mtx");

    if (!A.is_symmetric()) {
        throw std::runtime_error("Matrix is not symmetric");
    }

    std::vector<double> b(A.rows(), 1.0);

    std::cout << "Cholesky decomposition started..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<double> L_vals = CholeskySolver::decompose(A);

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Decomposition time: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " ms" << std::endl;

    std::cout << "Solver started..." << std::endl;
    start = std::chrono::high_resolution_clock::now();

    std::vector<double> x = CholeskySolver::solve(A, b, L_vals);

    end = std::chrono::high_resolution_clock::now();
    std::cout << "Solving time: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " ms" << std::endl;
    return 0;
}