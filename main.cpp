#include "CSR_Matrix.h"
#include "Cholesky_solver.h"
#include <iostream>
#include <vector>
#include <chrono>

int main() {
    using Log = std::tuple<int, std::pair<long, long>>;
    std::vector<Log> tests_res(21);

    Matrix A = load_mtx("5k.mtx");
    std::vector<double> b(A.rows(), 1.0);

    for (int thread_num = 1; thread_num <= 12; ++thread_num) {
        tests_res.emplace_back(thread_num, CholeskySolver::test(thread_num, A, b));
    }
    for (int thread_num = 20; thread_num <= 100; thread_num += 10) {
        tests_res.emplace_back(thread_num, CholeskySolver::test(thread_num, A, b));
    }

    std::cout << "Number of threads | Decomposition time | Solving time | Total time" << std::endl;
    for (auto test : tests_res) {
        std::cout << std::get<0>(test) << '\t' << std::get<1>(test).first << '\t' << std::get<1>(test).second << std::endl;
    }
    return 0;
}
