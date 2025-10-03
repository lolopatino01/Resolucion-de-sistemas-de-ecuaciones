#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <fstream>
#include "gauss.h"
#include "sarrus.h"
#include "laplace.h"
#include "ejecutarAutomatico.h"
using namespace std;
// Alias locales (sin matrix.h)
using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;
// Función auxiliar para leer enteros positivos con validación
int leerEnteroPositivo(const std::string& mensaje) {
    std::string linea;
    int valor;
    while (true) {
        std::cout << mensaje;
        std::getline(std::cin, linea);

        if (linea.empty()) {
            std::cerr << "No puede dejar vacio, ingrese un numero.\n";
            continue;
        }
        std::istringstream iss(linea);
        if (iss >> valor && valor > 0) {
            return valor;
        }
        else {
            std::cerr << "Entrada invalida, debe ser un entero positivo.\n";
        }
    }
}
void ejecutarAutomatico() {
    try {
        // Lectura robusta de dimensiones
        int nEq_in = leerEnteroPositivo("Ingrese cantidad de ecuaciones: ");
        int nInc_in = leerEnteroPositivo("Ingrese cantidad de incognitas: ");

        size_t nEq = static_cast<size_t>(nEq_in);
        size_t nInc = static_cast<size_t>(nInc_in);

        if (nEq != nInc) {
            throw invalid_argument("El sistema debe ser cuadrado (mismo numero de ecuaciones e incognitas).");
        }

        Matrix A(nEq, Vector(nInc));
        Vector b(nEq);

        cout << "Ingrese los coeficientes de la matriz A y el vector b:\n";
        for (size_t i = 0; i < nEq; i++) {
            for (size_t j = 0; j < nInc; j++) {
                cout << "A[" << i << "][" << j << "]: ";
                if (!(cin >> A[i][j])) throw invalid_argument("Entrada invalida en coeficiente.");
            }
            cout << "b[" << i << "]: ";
            if (!(cin >> b[i])) throw invalid_argument("Entrada invalida en termino independiente.");
        }
        // Abrimos archivo para guardar resultados
        ofstream archivo("D:/julian/resultados.txt", ios::app);
        if (!archivo) throw runtime_error("No se pudo abrir el archivo de resultados.");
        // Selección automática
        if (nEq == 3) {
            cout << "\n=== Metodo seleccionado: Regla de Sarrus (3x3) ===\n";
            double det = determinanteSarrus(A);
            cout << "Determinante = " << det << "\n";
            cout << "Resolviendo sistema con Gauss-Jordan...\n";
            Vector sol = resolverGaussJordan(A, b);
            for (size_t i = 0; i < sol.size(); i++) {
                cout << "x" << i + 1 << " = " << sol[i] << "\n";
            }
            // Guardar en archivo
            archivo << "=== Metodo seleccionado: Regla de Sarrus (3x3) ===\n";
            archivo << "Matriz A:\n";
            for (size_t i = 0; i < nEq; i++) {
                for (size_t j = 0; j < nInc; j++) archivo << A[i][j] << " ";
                archivo << "\n";
            }
            archivo << "Vector b:\n";
            for (size_t i = 0; i < nEq; i++) archivo << b[i] << "\n";
            archivo << "Determinante = " << det << "\n";
            for (size_t i = 0; i < sol.size(); i++) archivo << "x" << i + 1 << " = " << sol[i] << "\n";
            archivo << "----------------------------------------\n";
        }
        else if (nEq == 2) {
            cout << "\n=== Metodo seleccionado: Determinante 2x2 + Gauss ===\n";
            double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
            cout << "Determinante = " << det << "\n";
            Vector sol = resolverGauss(A, b);
            for (size_t i = 0; i < sol.size(); i++) {
                cout << "x" << i + 1 << " = " << sol[i] << "\n";
            }
            archivo << "=== Metodo seleccionado: Determinante 2x2 + Gauss ===\n";
            archivo << "Matriz A:\n";
            for (size_t i = 0; i < nEq; i++) {
                for (size_t j = 0; j < nInc; j++) archivo << A[i][j] << " ";
                archivo << "\n";
            }
            archivo << "Vector b:\n";
            for (size_t i = 0; i < nEq; i++) archivo << b[i] << "\n";
            archivo << "Determinante = " << det << "\n";
            for (size_t i = 0; i < sol.size(); i++) archivo << "x" << i + 1 << " = " << sol[i] << "\n";
            archivo << "----------------------------------------\n";
        }
        else if (nEq > 3) {
            cout << "\n=== Metodo seleccionado: Expansion de Laplace (n x n) ===\n";
            double det = determinanteLaplace(A);
            cout << "Determinante = " << det << "\n";
            cout << "Resolviendo sistema con Gauss...\n";
            Vector sol = resolverGauss(A, b);
            for (size_t i = 0; i < sol.size(); i++) {
                cout << "x" << i + 1 << " = " << sol[i] << "\n";
            }

            archivo << "=== Metodo seleccionado: Expansion de Laplace (n x n) ===\n";
            archivo << "Matriz A:\n";
            for (size_t i = 0; i < nEq; i++) {
                for (size_t j = 0; j < nInc; j++) archivo << A[i][j] << " ";
                archivo << "\n";
            }
            archivo << "Vector b:\n";
            for (size_t i = 0; i < nEq; i++) archivo << b[i] << "\n";
            archivo << "Determinante = " << det << "\n";
            for (size_t i = 0; i < sol.size(); i++) archivo << "x" << i + 1 << " = " << sol[i] << "\n";
            archivo << "----------------------------------------\n";
        }
        else { // nEq == 1
            cout << "\n=== Metodo seleccionado: Gauss (n=1) ===\n";
            if (A[0][0] == 0) throw runtime_error("El sistema no tiene solucion unica (division por cero).");
            double sol = b[0] / A[0][0];
            cout << "x1 = " << sol << "\n";

            archivo << "=== Metodo seleccionado: Gauss (n=1) ===\n";
            archivo << "A[0][0] = " << A[0][0] << "\n";
            archivo << "b[0] = " << b[0] << "\n";
            archivo << "x1 = " << sol << "\n";
            archivo << "----------------------------------------\n";
        }
    }
    catch (const invalid_argument& e) {
        cerr << "Excepcion de argumento invalido: " << e.what() << "\n";
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }
    catch (const runtime_error& e) {
        cerr << "Excepcion en tiempo de ejecucion: " << e.what() << "\n";
    }
    catch (const bad_alloc& e) {
        cerr << "Excepcion de memoria: " << e.what() << "\n";
    }
    catch (...) {
        cerr << "Se produjo una excepcion desconocida.\n";
    }
}