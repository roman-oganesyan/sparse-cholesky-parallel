#include "Cholesky_solver.h"

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

    // Копируем значения матрицы для модификации
    std::vector<double> L_vals(values);

    // Структура для быстрого доступа к элементам по столбцам
    std::vector<std::unordered_map<int, double>> col_map(n);
    for (int i = 0; i < n; ++i) {
        for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; ++idx) {
            int j = col_inds[idx];
            if (i >= j) { // Нижний треугольник
                col_map[j][i] = L_vals[idx];
            }
        }
    }

    // Основной цикл разложения
    for (int j = 0; j < n; ++j) {
        // Обработка диагонального элемента
        double diag_val = col_map[j][j];
        for (int k = 0; k < j; ++k) {
            if (col_map[k].find(j) != col_map[k].end()) {
                diag_val -= col_map[k][j] * col_map[k][j];
            }
        }

        if (diag_val <= 0) {
            throw std::runtime_error("Matrix is not positive definite");
        }

        diag_val = std::sqrt(diag_val);
        col_map[j][j] = diag_val;

        // Обновление поддиагональных элементов
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

    // Конвертируем результат обратно в плоский массив
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

    // Прямая подстановка (Ly = b)
    std::vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = b[i];
        for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; ++idx) {
            int j = col_inds[idx];
            if (j < i) { // Элементы левее диагонали
                sum -= L_vals[idx] * y[j];
            }
        }
        // Диагональный элемент
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
        // Используем доступ по столбцам через предварительное построение карты
        for (int j = i + 1; j < n; ++j) {
            // Поиск элемента L(j,i)
            for (int idx = row_ptrs[j]; idx < row_ptrs[j + 1]; ++idx) {
                if (col_inds[idx] == i) {
                    sum -= L_vals[idx] * x[j];
                    break;
                }
            }
        }
        // Диагональный элемент
        for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; ++idx) {
            if (col_inds[idx] == i) {
                x[i] = sum / L_vals[idx];
                break;
            }
        }
    }

    return x;
}