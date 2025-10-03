#pragma once
#include <vector>
#include <stdexcept>  // Para usar runtime_error, invalid_argument, etc.
using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;
// Resolver sistema por Eliminación de Gauss
Vector resolverGauss(Matrix A, Vector b);
// Resolver sistema por Gauss-Jordan
Vector resolverGaussJordan(Matrix A, Vector b);
// ejecutable

