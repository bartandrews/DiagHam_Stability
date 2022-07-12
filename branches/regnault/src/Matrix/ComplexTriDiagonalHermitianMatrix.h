////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of complex tridiagonal hermitian matrix              //
//                                                                            //
//                        last modification : 18/01/2001                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef COMPLEXTRIDIAGONALHERMITIANMATRIX_H
#define COMPLEXTRIDIAGONALHERMITIANMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Polynomial/Polynomial.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class ComplexTriDiagonalHermitianMatrix : protected Matrix
{

  friend class ComplexVector;
  friend class HermitianMatrix;

 protected:

  double* DiagonalElements;

  double* RealUpperDiagonalElements;
  double* ImaginaryUpperDiagonalElements;

  GarbageFlag Flag;

 public:

  // default constructor
  //
  ComplexTriDiagonalHermitianMatrix();

  // constructor from matrix elements (without duplicating datas)
  //
  // diagonal = pointer to diagonal element array
  // realUpperDiagonal = pointer to real part of upper diagonal element
  // imaginaryUpperDiagonal = pointer to imaginary part of upper diagonal element
  // dimension = matrix dimension
  ComplexTriDiagonalHermitianMatrix(double* diagonal, double* realUpperDiagonal, double* imaginaryUpperDiagonal, int dimension);

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  ComplexTriDiagonalHermitianMatrix(const ComplexTriDiagonalHermitianMatrix& M);

  // copy constructor from a real tridiagonal symmetric matrix (without duplicating diagonal elements)
  //
  // M = matrix to copy
  ComplexTriDiagonalHermitianMatrix(const RealTriDiagonalSymmetricMatrix& M);

  // destructor
  //
  ~ComplexTriDiagonalHermitianMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  ComplexTriDiagonalHermitianMatrix& operator = (const ComplexTriDiagonalHermitianMatrix& M);

  // assignement from  a real tridiagonal symmetric matrix (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix  
  ComplexTriDiagonalHermitianMatrix& operator = (const RealTriDiagonalSymmetricMatrix& M);

  // return pointer on a clone matrix (without duplicating datas)
  //
  // retrun value = pointer on new matrix 
  Matrix* Clone ();

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  void SetMatrixElement(int i, int j, double x);

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  void SetMatrixElement(int i, int j, const Complex& x);

  // add a value to a matrix element
  //
  // i = line position
  // j = column position
  // x = value to add to matrix element
  void AddToMatrixElement(int i, int j, double x);

  // add a value  a matrix element
  //
  // i = line position
  // j = column position
  // x = value to add to matrix element
  void AddToMatrixElement(int i, int j, const Complex& x);

  // add two matrices
  //
  // M1 = first matrix
  // M2 = second matrix
  // return value = sum of the two matrices
  friend ComplexTriDiagonalHermitianMatrix operator + (const ComplexTriDiagonalHermitianMatrix& M1, 
						       const ComplexTriDiagonalHermitianMatrix& M2);

  // add two matrices where the left one is a real tridiagonal symmetric matrix
  //
  // M1 = left matrix
  // M2 = right matrix
  // return value = sum of the two matrices
  friend ComplexTriDiagonalHermitianMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, 
						       const ComplexTriDiagonalHermitianMatrix& M2);

  // add two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M1 = left matrix
  // M2 = right matrix
  // return value = sum of the two matrices
  friend ComplexTriDiagonalHermitianMatrix operator + (const ComplexTriDiagonalHermitianMatrix& M1, 
						       const RealTriDiagonalSymmetricMatrix& M2);

  // substract two matrices
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend ComplexTriDiagonalHermitianMatrix operator - (const ComplexTriDiagonalHermitianMatrix& M1, 
						       const ComplexTriDiagonalHermitianMatrix& M2);

  // substract two matrices where the left one is a real tridiagonal symmetric matrix
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend ComplexTriDiagonalHermitianMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, 
						       const ComplexTriDiagonalHermitianMatrix& M2);

  // substract two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend ComplexTriDiagonalHermitianMatrix operator - (const ComplexTriDiagonalHermitianMatrix& M1,
						       const RealTriDiagonalSymmetricMatrix& M2);

  // multiply a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend ComplexTriDiagonalHermitianMatrix operator * (const ComplexTriDiagonalHermitianMatrix& M, double x);

  // multiply a matrix by a real number (left multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend ComplexTriDiagonalHermitianMatrix operator * (double x, const ComplexTriDiagonalHermitianMatrix& M);

  // divide a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = division result
  friend ComplexTriDiagonalHermitianMatrix operator / (const ComplexTriDiagonalHermitianMatrix& M, double x);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  ComplexTriDiagonalHermitianMatrix& operator += (const ComplexTriDiagonalHermitianMatrix& M);

  // add two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  ComplexTriDiagonalHermitianMatrix& operator += (const RealTriDiagonalSymmetricMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  ComplexTriDiagonalHermitianMatrix& operator -= (const ComplexTriDiagonalHermitianMatrix& M);

  // substract two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  ComplexTriDiagonalHermitianMatrix& operator -= (const RealTriDiagonalSymmetricMatrix& M);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  ComplexTriDiagonalHermitianMatrix& operator *= (double x);

  // divide a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  ComplexTriDiagonalHermitianMatrix& operator /= (double x) ;

  // evaluate matrix trace
  //
  // return value = matrix trace 
  double Tr ();

  // evaluate matrix determinant
  //
  // return value = matrix determinant 
  double Det ();

  // return matrix characteritic equation
  //
  // return value =  reference one polynomial corresponding to matrix characteritic equation  
  Polynomial& CharacteristicEquation();

  // evaluate a normalized eigenvector for a given eigenvalue (supposing the eigenvalue is non-degenerate)
  //
  // eigenvalue = eigenvalue to use
  // eigenvector = vector where the eigenvector has to be stored
  // return value = reference on eigenvector
  ComplexVector& Eigenvector(double eigenvalue, ComplexVector& eigenvector);

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const ComplexTriDiagonalHermitianMatrix& P);

};

#endif
