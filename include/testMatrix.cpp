/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 */

/**
 * Unit tests for the matrix class
 */

#include <iostream>
#include <exception>

#include <cstdlib>

#include "Matrix.hpp"

// TODO: Use a proper unit-testing framework, like boost's.

// Explicit instantiation of all methods of Matrix
template class anpi::Matrix<double>;
template class anpi::Matrix<float>;

/// Minimalistic testing entry
int testsOk=0;
int testsFailed=0;
void test(const bool state,const std::string& msg="") {
  if (!state) {
    ++testsFailed;
    std::cout << "Failed: " << msg << std::endl;
  } else {
    ++testsOk;
  }
  
}


/**
 * Instantiate and test the methods of the Matrix class
 */
void testMatrix() {
  // Constructors
  { // default
    anpi::Matrix<float> a;
    test( a.rows() == 0 ,"Wrong rows");
    test( a.cols() == 0 ,"Wrong cols");
  }
  { // unitilialized
    anpi::Matrix<float> a(2,3,anpi::Matrix<float>::DoNotInitialize);
    test( a.rows() == 2,"Wrong rows");
    test( a.cols() == 3,"Wrong cols");
  }
  { // default initialized
    anpi::Matrix<float> a(3,2);
    test( a.rows() == 3,"Wrong rows");
    test( a.cols() == 2,"Wrong cols");
    test( a(0,0) == 0.f,"Wrong value");
  }
  { // default initialized
    anpi::Matrix<double> a(3,2,4.);
    test( a.rows() == 3,"Wrong rows");
    test( a.cols() == 2,"Wrong cols");
    test( a(0,0) == 4.,"Wrong value");
  }
  { // initializer_list
    anpi::Matrix<float> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    test( a.rows() == 3,"Wrong rows");
    test( a.cols() == 4,"Wrong cols");
    
    test( a(0,0) == 1.f,"Wrong value" );
    test( a(1,2) == 7.f,"Wrong value" );
    test( a(2,3) == 12.f,"Wrong value" );
  }
  { // Copy constructor
    anpi::Matrix<float> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    anpi::Matrix<float> b(a);

    test( a==b                ,"Not equal");
    test( b.rows() == 3       ,"Wrong rows");
    test( b.cols() == 4       ,"Wrong cols");
    test( b.data() != a.data(),"Wrong data");
  }

  { // Move constructor
    anpi::Matrix<float> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    anpi::Matrix<float> b(std::move(a));

    test( b.rows() == 3,"Wrong rows");
    test( b.cols() == 4,"Wrong cols");

    test( a.empty(), "Not moved");
  }
  { // Mem constructor
    anpi::Matrix<float> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    anpi::Matrix<float> b(a.rows(),a.cols(),a.data());

    test( a==b                ,"Not equal");
    test( b.rows() == 3       ,"Wrong rows");
    test( b.cols() == 4       ,"Wrong cols");
    test( b.data() != a.data(),"Wrong data");
  }
  { // == and !=
    anpi::Matrix<int> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    anpi::Matrix<int> b = { {1,2,3,4},{5,6,8,8},{9,10,11,12} };

    test( (a!=b), "Equal");

    b(1,2)=7;

    test( (a==b), "Not equal");
  }
  { // Move assignment
    anpi::Matrix<int> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    anpi::Matrix<int> c(a);
    anpi::Matrix<int> b;
    b=std::move(a);
    test(a.empty(), "Not moved");
    test(!b.empty(), "Not moved");
    test(b.rows()==3,"Wrong rows");
    test(b.cols()==4,"Wrong cols");
    test(b==c,"Wrong content");
  }
  { // assignment
    anpi::Matrix<int> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    anpi::Matrix<int> b;
    b=a;
    test(a==b,"Wrong content");
  }
  { // swap
    anpi::Matrix<int> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    anpi::Matrix<int> b = { {13,14},{15,16} };

    anpi::Matrix<int> c(a);
    anpi::Matrix<int> d(b);

    test(a==c,"Not correctly copy constructed");
    test(d==b,"Not correctly copy constructed");
    
    c.swap(d);
    test(a==d,"Not correctly swapped");
    test(b==c,"Not correctly swapped");
  }
  {
    anpi::Matrix<int> a = { {1,2,3},{ 4, 5, 6} };
    anpi::Matrix<int> b = { {7,8,9},{10,11,12} };
    anpi::Matrix<int> r = { {8,10,12},{14,16,18} };
    
    anpi::Matrix<int> c(a);
    c+=b;
    test(c==r,"Wrong inplace sum");
    c=a+b;
    test(c==r,"Wrong external sum");


    c=anpi::Matrix<int>{ {1,2,3},{ 4, 5, 6} } + b;
    test(c==r,"Wrong cast sum");

    c=a+anpi::Matrix<int>{ {7,8,9},{10,11,12} };
    test(c==r,"Wrong cast sum");
  }

  {
    anpi::Matrix<int> a = { {1,2,3},{ 4, 5, 6} };
    anpi::Matrix<int> b = { {7,8,9},{10,11,12} };
    anpi::Matrix<int> r = { {-6,-6,-6},{-6,-6,-6} };
    
    anpi::Matrix<int> c(a);
    c-=b;
    test(c==r,"Wrong inplace subtract");
    c=a-b;
    test(c==r,"Wrong external subtract");


    c=anpi::Matrix<int>{ {1,2,3},{ 4, 5, 6} } - b;
    test(c==r,"Wrong cast subtract");

    c=a-anpi::Matrix<int>{ {7,8,9},{10,11,12} };
    test(c==r,"Wrong cast subtract");
  }
  
}

int main() {
  try {
    testMatrix();
    std::cout << "Total tests: " << testsOk + testsFailed << std::endl;
    std::cout << "  Succeeded: " << testsOk               << std::endl;
    std::cout << "     Failed: " << testsFailed           << std::endl;
    
  }
  catch(const std::exception& exc) {
    std::cout << "Error: " << exc.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
  
