#include <iostream>
#include "include/Tarea6.hpp"


int main (){
    
    std::cout<<"--- Test Tarea 6 Eigensistemas ---"<<std::endl;
    anpi::Matrix<double> testA = randomSymmetricSqr<double>(10);

    std::cout<<"Matriz de 10x10 generada con numeros aleatorios de -100 a 100"<<std::endl;
    testA.printm();

    anpi::Matrix<double> testE(10, 10, anpi::Matrix<double>::DoNotInitialize);
    std::vector<double> testVal;

    const clock_t begin_time = clock();
    eig<double>(testA, testVal, testE);

    std::cout<<"Eigenvectores producidos usando LAPACK:"<<std::endl;
    testE.printm();

    std::cout<<"Eigenvalores producidos usando LAPACK:"<<std::endl;
    for (int i = 0; i<10; i++){
        std::cout<<testVal[i]<<std::endl;
    }

    std::cout<<std::endl<<"Tiempo que tomo el algortimo de LAPACK: "<< float(clock() - begin_time)/CLOCKS_PER_SEC<<std::endl<<std::endl;

    anpi::Matrix<double> X = reconstruirMatrix<double>(testE, testVal);

    std::cout<<"Matriz reconstruida a partir de los resultados de LAPACK:"<<std::endl;
    X.printm();
    std::cout<<std::endl<<"Norma de la matriz reconstruida a partir de los resultados de LAPACK:"<<std::endl;
    std::cout<<normaMatrix<double>(testA, X)<<std::endl<<std::endl;

    anpi::Matrix<double> testB;
    std::vector<double> v;

    const clock_t begin_time2 = clock();

    jacobi<double>(testA, v, testB);

    std::cout<<"Eigenvectores producidos usando Jacobi:"<<std::endl;
    testB.printm();


    std::cout<<"Eigenvalores producidos usando Jacobi:"<<std::endl;
    for (std::vector<double>::const_iterator i = v.begin(); i != v.end(); ++i)
        std::cout << *i << ' '<<std::endl;

    std::cout<<std::endl<< "Tiempo que tomo el algortimo de Jacobi: "<< float(clock() - begin_time2)/CLOCKS_PER_SEC<<std::endl<<std::endl;

    anpi::Matrix<double> Y = reconstruirMatrix<double>(testE, testVal);

    std::cout<<"Matriz reconstruida a partir de los resultados de Jacobi:"<<std::endl;
    Y.printm();
    std::cout<<std::endl<<"Norma de la matriz reconstruida a partir de los resultados de Jacobi:"<<std::endl;
    std::cout<<normaMatrix<double>(testA, Y)<<std::endl<<std::endl;

    sort<double>(testVal, testE);

    std::cout<<"Comprobacion de la funcion sort() (comparando con los resultados de LAPACK)"<<std::endl<<"Eigenvectores producidos por LAPACK ordenados por sort():"<<std::endl;
    testE.printm();

    std::cout<<"Eigenvalores producidos por Jacobi ordenados por sort():"<<std::endl;
    for (int i = 0; i<10; i++){
        std::cout<<testVal[i]<<std::endl;
    }


    return(1);
}