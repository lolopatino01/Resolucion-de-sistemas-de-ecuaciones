#pragma once
#include <vector>
using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;
Vector resolverGauss(Matrix A, Vector b, bool mostrarPasos = true);
Vector resolverGaussJordan(Matrix A, Vector b, bool mostrarPasos = true);
