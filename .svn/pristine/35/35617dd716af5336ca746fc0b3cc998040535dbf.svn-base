////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class of complex matrix                         //
//                                                                            //
//                        last modification : 05/02/2001                      //
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


#ifndef COMPLEXMATRIX_H
#define COMPLEXMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#include "Output/MathematicaOutput.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include <iostream>


using std::ostream;


class ComplexVector;


class ComplexMatrix : protected Matrix
{

  friend class HermitianMatrix;
  friend class RealVector;
  friend class ComplexVector;
  friend class RealTriDiagonalSymmetricMatrix;
  friend class SingleParticle;

 protected:

  ComplexVector* Columns; 
  int* ColumnGarbageFlag;

 public:

  // default constructor
  //
  ComplexMatrix();

  // constructor from matrix elements (without duplicating datas)
  //
  // columns = pointer an array of vector
  // nbrColumn = number of columns
  ComplexMatrix(ComplexVector* columns, int nbrColumn);

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  ComplexMatrix(const ComplexMatrix& M);

  // destructor
  //
  ~ComplexMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  ComplexMatrix& operator = (const ComplexMatrix& M);

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

  // Resize matrix
  //
  // nbrRow = new number of rows
  // nbrColumn = new number of columns
  void Resize (int nbrRow, int nbrColumn);

  // Resize matrix and set to zero all elements that have been added
  //
  // nbrRow = new number of rows
  // nbrColumn = new number of columns
  void ResizeAndClean (int nbrRow, int nbrColumn);

  // add two matrices
  //
  // M1 = first matrix
  // M2 = second matrix
  // return value = sum of the two matrices
  friend ComplexMatrix operator + (const ComplexMatrix& M1, const ComplexMatrix& M2);

  // add two matrices where the left one is a real tridiagonal symmetric matrix
  //
  // M1 = left matrix
  // M2 = right matrix
  // return value = sum of the two matrices
  friend ComplexMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const ComplexMatrix& M2);

  // add two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M1 = left matrix
  // M2 = right matrix
  // return value = sum of the two matrices
  friend ComplexMatrix operator + (const ComplexMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);

  // substract two matrices
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend ComplexMatrix operator - (const ComplexMatrix& M1, const ComplexMatrix& M2);

  // substract two matrices where the left one is a real tridiagonal symmetric matrix
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend ComplexMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, const ComplexMatrix& M2);

  // substract two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend ComplexMatrix operator - (const ComplexMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);

  // multiply two matrices
  //
  // M1 = first matrix
  // M2 = matrix to multiply to M1
  // return value = product of the two matrices
  friend ComplexMatrix operator * (const ComplexMatrix& M1, const ComplexMatrix& M2);

  // multiply a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend ComplexMatrix operator * (const ComplexMatrix& M, double x);

  // multiply a matrix by a real number (left multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend ComplexMatrix operator * (double x, const ComplexMatrix& M);

  // divide a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = division result
  friend ComplexMatrix operator / (const ComplexMatrix& M, double x);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  ComplexMatrix& operator += (const ComplexMatrix& M);

  // add two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  ComplexMatrix& operator += (const RealTriDiagonalSymmetricMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  ComplexMatrix& operator -= (const ComplexMatrix& M);

  // substract two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  ComplexMatrix& operator -= (const RealTriDiagonalSymmetricMatrix& M);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  ComplexMatrix& operator *= (double x);

  // divide a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  ComplexMatrix& operator /= (double x);

  // normalize matrix column vectors
  //
  // return value = reference on current matrix
  ComplexMatrix& NormalizeColumns ();

  // orthonormalize matrix column vectors
  //
  // return value = reference on current matrix
  ComplexMatrix& OrthoNormalizeColumns ();

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const ComplexMatrix& P);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexMatrix& P);

};

#endif
