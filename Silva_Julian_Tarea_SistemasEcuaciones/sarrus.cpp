#include "sarrus.h"
#include <iostream>
#include <stdexcept>
using namespace std;
double determinanteSarrus(const Matrix& A) {
    if (A.size() != 3 || A[0].size() != 3) {
        throw invalid_argument("La matriz debe ser 3x3 para aplicar la regla de Sarrus.");
    }

    double det = 0.0;
    det = A[0][0] * A[1][1] * A[2][2]
        + A[0][1] * A[1][2] * A[2][0]
        + A[0][2] * A[1][0] * A[2][1]
        - A[0][2] * A[1][1] * A[2][0]
        - A[0][0] * A[1][2] * A[2][1]
        - A[0][1] * A[1][0] * A[2][2];

    return det;
}