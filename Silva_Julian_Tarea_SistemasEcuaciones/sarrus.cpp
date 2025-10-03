#include "sarrus.h"
#include <stdexcept>

double determinanteSarrus(const Matrix& A) {
    if (A.size() != 3 || A[0].size() != 3)
        throw std::invalid_argument("La matriz debe ser 3x3 para aplicar la regla de Sarrus.");

    return A[0][0] * A[1][1] * A[2][2]
        + A[0][1] * A[1][2] * A[2][0]
        + A[0][2] * A[1][0] * A[2][1]
        - A[0][2] * A[1][1] * A[2][0]
        - A[0][0] * A[1][2] * A[2][1]
        - A[0][1] * A[1][0] * A[2][2];
}
