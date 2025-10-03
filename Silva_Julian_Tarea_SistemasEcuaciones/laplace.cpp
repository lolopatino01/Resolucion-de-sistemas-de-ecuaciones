#include "laplace.h"
#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;
using Matrix = vector<vector<double>>;
// Función recursiva para calcular el determinante por expansión de Laplace
double determinanteLaplace(const Matrix& A) {
    size_t n = A.size();
    if (n == 0 || A[0].size() != n) {
        throw invalid_argument("La matriz debe ser cuadrada.");
    }
    // Casos base
    if (n == 1) return A[0][0];
    if (n == 2) return A[0][0] * A[1][1] - A[0][1] * A[1][0];
    // Caso general: expansión por la primera fila
    double det = 0.0;
    for (size_t j = 0; j < n; j++) {
        // Construir submatriz excluyendo fila 0 y columna j
        Matrix sub(n - 1, vector<double>(n - 1));
        for (size_t r = 1; r < n; r++) {
            size_t colSub = 0;
            for (size_t c = 0; c < n; c++) {
                if (c == j) continue;
                sub[r - 1][colSub++] = A[r][c];
            }
        }
        // Cofactor con signo alternante
        double cofactor = ((j % 2 == 0) ? 1 : -1) * A[0][j] * determinanteLaplace(sub);
        det += cofactor;
    }
    return det;
}
// Función que pide datos al usuario y ejecuta el cálculo