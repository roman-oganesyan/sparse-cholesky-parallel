#include "CSR_Matrix.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <unordered_map>

class CholeskySolver {
public:
    static std::vector<double> decompose(const Matrix& A);

    static std::vector<double> solve(
        const Matrix& A,
        const std::vector<double>& b,
        const std::vector<double>& L_vals
    );
};