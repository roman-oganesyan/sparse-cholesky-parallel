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

    Matrix A = load_mtx("25k.mtx");
    std::vector<double> b(A.rows(), 1.0);

    for (int thread_num = 1; thread_num <= 100; ++thread_num) {
        tests_res.emplace_back(thread_num, CholeskySolver::test(thread_num, A, b));
        mem << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
        std::cout << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
    }
    /*
    std::cout << "Number of threads | Decomposition time | Solving time | Total time" << std::endl;
    tests_res.emplace_back(7, CholeskySolver::test(7, A, b));
    mem << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
    std::cout << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
    tests_res.emplace_back(8, CholeskySolver::test(8, A, b));
    mem << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
    std::cout << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
    tests_res.emplace_back(20, CholeskySolver::test(20, A, b));
    mem << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
    std::cout << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
    tests_res.emplace_back(40, CholeskySolver::test(40, A, b));
    mem << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
    std::cout << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
    tests_res.emplace_back(50, CholeskySolver::test(50, A, b));
    mem << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
    std::cout << std::get<0>(tests_res.back()) << '\t' << std::get<1>(tests_res.back()).first << '\t' << std::get<1>(tests_res.back()).second << std::endl;
    */
   
    return 0;
}
