#include "gauss.h"
#include <cmath>
#include <iostream> 
using namespace std;    
static constexpr double EPS = 1e-12;
Vector resolverGauss(Matrix A, Vector b) {
    size_t n = A.size();
    if (A.empty() || A[0].size() != n || b.size() != n) {
        throw std::invalid_argument("La matriz no es cuadrada o el vector no coincide en tamaño.");
    }
    // Eliminación hacia adelante
    for (size_t k = 0; k < n; ++k) {
        if (std::fabs(A[k][k]) < EPS) {
            throw std::runtime_error("Pivote nulo: el sistema es singular.");
        }
        for (size_t i = k + 1; i < n; ++i) {
            double m = A[i][k] / A[k][k];
            for (size_t j = k; j < n; ++j) {
                A[i][j] -= m * A[k][j];
            }
            b[i] -= m * b[k];
        }
    }
    // Sustitución hacia atrás
    Vector x(n, 0.0);
    for (int i = int(n) - 1; i >= 0; --i) {
        double sum = 0.0;
        for (size_t j = i + 1; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        if (std::fabs(A[i][i]) < EPS) {
            throw std::runtime_error("División por cero en la diagonal.");
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
    return x;
}
Vector resolverGaussJordan(Matrix A, Vector b) {
    size_t n = A.size();
    if (A.empty() || A[0].size() != n || b.size() != n) {
        throw std::invalid_argument("La matriz no es cuadrada o el vector no coincide en tamaño.");
    }

    for (size_t k = 0; k < n; ++k) {
        if (std::fabs(A[k][k]) < EPS) {
            throw std::runtime_error("Pivote nulo: el sistema es singular.");
        }

        // Normalizar fila pivote
        double piv = A[k][k];
        for (size_t j = 0; j < n; ++j) A[k][j] /= piv;
        b[k] /= piv;

        // Eliminar en todas las demás filas
        for (size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            double m = A[i][k];
            for (size_t j = 0; j < n; ++j) A[i][j] -= m * A[k][j];
            b[i] -= m * b[k];
        }
    }
    return b; // b ya contiene la solución
}