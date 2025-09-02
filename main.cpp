#include "CSR_Matrix.h"
#include "Cholesky_solver.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>

int main() {
    using Log = std::tuple<int, std::pair<long, long>>;
    std::vector<Log> tests_res(21);
    std::ofstream mem("mem.txt");

    Matrix A = load_mtx("5k.mtx");
    std::vector<double> b(A.rows(), 1.0);

    std::cout << "Number of threads | Decomposition time | Solving time | Total time" << std::endl;
    for (int thread_num = 1; thread_num <= 12; ++thread_num) {
        tests_res.emplace_back(thread_num, CholeskySolver::test(thread_num, A, b));
        mem << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
        std::cout << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
    }
 
    for (int thread_num = 20; thread_num <= 100; thread_num += 10) {
        tests_res.emplace_back(thread_num, CholeskySolver::test(thread_num, A, b));
        mem << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
        std::cout << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
    }

    return 0;
}
