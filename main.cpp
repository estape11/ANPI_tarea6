#include <iostream>
#include "include/Matrix.hpp"
#include "include/Tarea6.hpp"
int main() {
    std::cout << "Hello, World!" << std::endl;
    anpi::Matrix<int> jej=randomSymmetricSqr<int>(10);
    printMatrix(jej);
    return 0;
}