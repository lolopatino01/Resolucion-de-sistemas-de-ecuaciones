#pragma once
#include <vector>
// Alias de tipos para todo el proyecto
using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;
double determinanteSarrus(const Matrix& A);