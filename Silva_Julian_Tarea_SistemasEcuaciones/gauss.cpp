#include "gauss.h"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <sstream>

static constexpr double EPS = 1e-12;

static void imprimirPaso(const Matrix& A, const Vector& b, const std::string& desc, bool mostrar) {
    if (!mostrar) return;
    std::cout << desc << "\n";
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            std::cout << std::setw(8) << A[i][j] << " ";
        }
        std::cout << "| " << std::setw(8) << b[i] << "\n";
    }
    std::cout << "-------------------------\n";
}

Vector resolverGauss(Matrix A, Vector b, bool mostrarPasos) {
    size_t n = A.size();
    if (A.empty() || A[0].size() != n || b.size() != n)
        throw std::invalid_argument("La matriz no es cuadrada o el vector no coincide en tamaño.");

    imprimirPaso(A, b, "Matriz aumentada inicial", mostrarPasos);

    for (size_t k = 0; k < n; ++k) {
        if (std::fabs(A[k][k]) < EPS)
            throw std::runtime_error("Pivote nulo: el sistema es singular.");

        for (size_t i = k + 1; i < n; ++i) {
            double m = A[i][k] / A[k][k];
            for (size_t j = k; j < n; ++j) A[i][j] -= m * A[k][j];
            b[i] -= m * b[k];

            std::ostringstream oss;
            oss << "R" << i << " = R" << i << " - " << m << "*R" << k;
            imprimirPaso(A, b, oss.str(), mostrarPasos);
        }
    }

    Vector x(n, 0.0);
    for (int i = int(n) - 1; i >= 0; --i) {
        double sum = 0.0;
        for (size_t j = i + 1; j < n; ++j) sum += A[i][j] * x[j];
        if (std::fabs(A[i][i]) < EPS)
            throw std::runtime_error("División por cero en la diagonal.");
        x[i] = (b[i] - sum) / A[i][i];
    }

    return x; // ? devolvemos la solución, no la imprimimos
}

Vector resolverGaussJordan(Matrix A, Vector b, bool mostrarPasos) {
    size_t n = A.size();
    if (A.empty() || A[0].size() != n || b.size() != n)
        throw std::invalid_argument("La matriz no es cuadrada o el vector no coincide en tamaño.");

    imprimirPaso(A, b, "Matriz aumentada inicial", mostrarPasos);

    for (size_t k = 0; k < n; ++k) {
        if (std::fabs(A[k][k]) < EPS)
            throw std::runtime_error("Pivote nulo: el sistema es singular.");

        double piv = A[k][k];
        for (size_t j = 0; j < n; ++j) A[k][j] /= piv;
        b[k] /= piv;

        {
            std::ostringstream oss;
            oss << "Normalizando fila " << k;
            imprimirPaso(A, b, oss.str(), mostrarPasos);
        }

        for (size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            double m = A[i][k];
            for (size_t j = 0; j < n; ++j) A[i][j] -= m * A[k][j];
            b[i] -= m * b[k];

            std::ostringstream oss;
            oss << "R" << i << " = R" << i << " - " << m << "*R" << k;
            imprimirPaso(A, b, oss.str(), mostrarPasos);
        }
    }

    return b; // 
}
