#include "CSR_Matrix.h"
#include <vector>
#include <chrono>

class CholeskySolver {
public:
    static std::vector<double> decompose(const Matrix& A);

    static std::vector<double> solve(
        const Matrix& A,
        const std::vector<double>& b,
        const std::vector<double>& L_vals
    );

    static std::pair<long, long> test(const Matrix& A, const std::vector<double>& b);
};