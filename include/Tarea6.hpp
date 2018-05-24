#include "Matrix.hpp"

/**
 * Util funcion para visualizar matrices
 * @tparam T
 * @param matrix
 */
template<typename T>
void printMatrix(anpi::Matrix<T> matrix) {
    for (unsigned int i = 0; i < matrix.rows(); i++) {
        std::cout<<"{ ";
        for (unsigned int j = 0; j < matrix.cols(); j++) {
            if(j==matrix.cols()-1){
                std::cout<<matrix[i][j];
            } else std::cout<<matrix[i][j]<<", ";
        }
        std::cout<<" }"<<std::endl;
    }
}

/**
 * Funcion que genera numeros aleatorios entre un rango de -100 a 100
 * @tparam T
 * @return
 */
template<typename T>
T getRandom() {
    T temp = T(rand() % 200 - 100); // rango -100 a  100
    return temp;
}

/**
 * Genera una matriz cuadrada llena de numeros aleatorios
 * @tparam T
 * @param N
 * @return
 */
template<typename T>
anpi::Matrix<T> randomSymmetricSqr(const size_t N) {
    anpi::Matrix<T> temp = anpi::Matrix<T>(N, N);
    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < N; j++) {
            temp[i][j] = getRandom<T>();
        }
    }
    return temp;
}
