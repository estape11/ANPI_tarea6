#include "include/Tarea6.hpp"

/**
 * Programa de prueba Tarea 6
 * @return
 */
int main() {

    std::cout << "--- Test Tarea 6 Eigensistemas ---" << std::endl << std::endl;
    anpi::Matrix<double> matrizA = randomSymmetricSqr<double>(10);

    std::cout << "-> Matriz cuadrada 10x10 generada con numeros aleatorios de -100 a 100" << std::endl;
    matrizA.printm();

    anpi::Matrix<double> matrizE(10, 10, anpi::Matrix<double>::DoNotInitialize);
    std::vector<double> testVal;

    const clock_t begin_time = clock();
    eig<double>(matrizA, testVal, matrizE);

    std::cout << "-> Eigenvectores utilizando LAPACK: " << std::endl;
    matrizE.printm();

    std::cout << "-> Eigenvalores utilizando LAPACK: " << std::endl;
    printVector(testVal);

    std::cout << std::endl << "-> Tiempo tomado por LAPACK: " << float(clock() - begin_time) / CLOCKS_PER_SEC
              << std::endl << std::endl;

    anpi::Matrix<double> matrizRecons = reconstruirMatrix<double>(matrizE, testVal);

    std::cout << "-> Matriz reconstruida: " << std::endl;
    matrizRecons.printm();

    std::cout << std::endl << "-> Norma de la matriz reconstruida (LAPACK): "
              << normaMatrix<double>(matrizA, matrizRecons) << std::endl << std::endl;

    anpi::Matrix<double> matrizB;
    std::vector<double> v;

    const clock_t begin_time2 = clock();

    jacobi<double>(matrizA, v, matrizB);

    std::cout << "-> Eigenvectores utilizando Jacobi:" << std::endl;
    matrizB.printm();

    std::cout << "-> Eigenvalores utilizando Jacobi:" << std::endl;
    printVector(v);

    std::cout << std::endl << "-> Tiempo tomado por Jacobi: " << float(clock() - begin_time2) / CLOCKS_PER_SEC
              << std::endl << std::endl;

    anpi::Matrix<double> matrizRecons2 = reconstruirMatrix<double>(matrizE, v);

    std::cout << "-> Matriz reconstruida a partir de los resultados de Jacobi:" << std::endl;
    matrizRecons2.printm();
    std::cout << std::endl << "-> Norma de la matriz reconstruida (Jacobi):"
              << normaMatrix<double>(matrizA, matrizRecons2) << std::endl << std::endl;

    sort<double>(v, matrizE);

    std::cout << "-> Eigenvectores producidos por LAPACK ordenados:" << std::endl;
    matrizE.printm();

    std::cout << "-> Eigenvalores producidos por Jacobi ordenados:" << std::endl;
    printVector(v);

    std::cout << "\n--- FIN ---" << std::endl;

    return true;
}