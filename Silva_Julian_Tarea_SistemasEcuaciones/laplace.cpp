#include "laplace.h"
#include <stdexcept>

double determinanteLaplace(const Matrix& A) {
    size_t n = A.size();
    if (n == 0 || A[0].size() != n)
        throw std::invalid_argument("La matriz debe ser cuadrada.");

    if (n == 1) return A[0][0];
    if (n == 2) return A[0][0] * A[1][1] - A[0][1] * A[1][0];

    double det = 0.0;
    for (size_t j = 0; j < n; j++) {
        Matrix sub(n - 1, std::vector<double>(n - 1));
        for (size_t r = 1; r < n; r++) {
            size_t colSub = 0;
            for (size_t c = 0; c < n; c++) {
                if (c == j) continue;
                sub[r - 1][colSub++] = A[r][c];
            }
        }
        double cofactor = ((j % 2 == 0) ? 1 : -1) * A[0][j] * determinanteLaplace(sub);
        det += cofactor;
    }
    return det;
}
