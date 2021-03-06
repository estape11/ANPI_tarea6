/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 */

#ifndef _ANPI_MATRIX_HPP
#define _ANPI_MATRIX_HPP

#include <cstddef>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <initializer_list>

namespace anpi
{
  /**
   * Row-major simple matrix class
   */
  template<typename T>
  class Matrix {
  private:
    /// All matrix data
    T* _data;
    /// Number of rows
    size_t _rows;
    /// Number of columns
    size_t _cols;

  public:
    /**
     * Type to make explicit a desired unitialized matrix
     */
    enum InitializationType {
      DoNotInitialize
    };


    /// Construct an empty matrix
    Matrix();
   
    /**
     * Construct a matrix rows x cols and do not initialize its
     * elements
     *
     * Call as Matrix<T> m(5,5,anpi::Matrix<T>::DoNotInitialize)
     */
    explicit Matrix(const size_t rows,
		    const size_t cols,
		    const InitializationType);

    /**
     * Construct a matrix rows x cols and initialize all
     * elements with the given value     
     */
    explicit Matrix(const size_t rows,
		    const size_t cols,
		    const T initVal=T());

    /**
     * Construct a matrix rows x cols and initialize all
     * elements with the memory content at the given pointer.
     *
     * Note that the content in the given block is copied
     */
    explicit Matrix(const size_t rows,
		    const size_t cols,
		    const T *const initMem);

    /**
     * Constructs a matrix from a std::initializer_list
     *
     * This allows to construct and initialize a 2x3 matrix in this way:
     *
     * \code
     * anpi::Matrix<int> m={ {1,2,3}, {4,5,6} };
     * \endcode
     */
    Matrix(std::initializer_list< std::initializer_list<T> > lst);
    
    
    /**
     * Copy constructor will do a deep copy on the given matrix
     */
    Matrix(const Matrix<T>& other);

    /**
     * The move constructor will transfer the contents of the other matrix
     * into the newly created one
     * 
     * Note that this is possible only for C++11 or later
     */
    Matrix(Matrix<T>&& other);
    
    /** 
     * Release all memory
     */
    ~Matrix();

    /**
     * Deep copy another matrix
     */
    Matrix<T>& operator=(const Matrix<T>& other);

    /**
     * Move assignment operator
     */
    Matrix<T>& operator=(Matrix<T>&& other);

    /**
     * Compare two matrices for equality
     *
     * This is slow, as all componentes are elementwise compared
     */
    bool operator==(const Matrix<T>& other) const;

    /**
     * Compare two matrices for equality
     *
     * This is slow, as all componentes are elementwise compared
     */
    bool operator!=(const Matrix<T>& other) const;
    
    /// Return pointer to a given row
    inline T* operator[](const size_t row) {return _data + row*_cols;}

    /// Return read-only pointer to a given row
    const T* operator[](const size_t row) const {return _data + row*_cols;}

    /// Return reference to the element at the r row and c column
    T& operator()(const size_t row,const size_t col) {
      return *(_data + (row*_cols + col));
    }

    /// Return const reference to the element at the r row and c column
    const T& operator()(const size_t row,const size_t col) const {
      return *(_data + (row*_cols + col));
    }

    /**
     * Swap the contents of the other matrix with this one
     */
    void swap(Matrix<T>& other);

    //Prints the whole matrix in the console
    void printm();

    //Encuantra la norma de la matriz con el metodo de Frobenius
    T norma();
    
    /**
     * Allocate memory for the given number of rows and cols
     */
    void allocate(const size_t row,
		  const size_t col);

    /**
     * Fill all elements of the matrix with the given value
     */
    void fill(const T val);

    /**
     * Fill all elements of the matrix with the given memory block
     *
     * The user must ensure that the given memory block has enough elements
     */
    void fill(const T* mem);

    /**
     * Check if the matrix is empty (zero rows or columns)
     */
    inline bool empty() const { return (_rows==0) || (_cols==0); }
    
    /**
     * Number of rows
     */
    inline size_t rows() const { return _rows; }

    /**
     * Number of columns
     */
    inline size_t cols() const { return _cols; }

    /**
     * Total number of entries (rows x cols)
     */
    inline size_t entries() const { return _rows*_cols; }

    /**
     * Pointer to data block
     */
    inline T* data() { return _data; }

    /**
     * Pointer to data block
     */
    inline const T* data() const { return _data; }

    /**
     * @name Arithmetic operators
     */
    //@{

    /// Sum this and another matrix, and leave the result in here
    Matrix& operator+=(const Matrix& other);

    /// Subtract another matrix to this one, and leave the result in here
    Matrix& operator-=(const Matrix& other);
    
    //@}

    
  }; // class Matrix


  // External arithmetic operators
  template<class T>
  Matrix<T> operator+(const Matrix<T>& a,
		      const Matrix<T>& b);

  template<class T>
  Matrix<T> operator-(const Matrix<T>& a,
		      const Matrix<T>& b);


  template<class T>
  Matrix<T> operator*(const Matrix<T>& a,
              const Matrix<T>& b);

} // namespace ANPI

// include the template implementations
#include "Matrix.tpp"

#endif
