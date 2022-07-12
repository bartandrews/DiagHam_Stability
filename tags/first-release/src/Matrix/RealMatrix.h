////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                             class of real matrix                           //
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


#ifndef REALMATRIX_H
#define REALMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#include "Output/MathematicaOutput.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"

#include <iostream>


using std::ostream;


class RealVector;


class RealMatrix : public Matrix
{

  friend class HermitianMatrix;
  friend class RealSymmetricMatrix;
  friend class RealAntisymmetricMatrix;
  friend class RealTriDiagonalSymmetricMatrix;
  friend class RealVector;
  friend class ComplexVector;

 protected:

  RealVector* Columns; 
  int* ColumnGarbageFlag;

 public:

  // default constructor
  //
  RealMatrix();

  // constructor for an empty matrix
  //
  // nbrRow = number of rows
  // nbrColumn = number of columns
  // zero = tue if matrix elements have to be set to zero
  RealMatrix(int nbrRow, int nbrColumn, bool zero = false);

  // constructor from matrix elements (without duplicating datas)
  //
  // columns = pointer an array of vector
  // nbrColumn = number of columns
  RealMatrix(RealVector* columns, int nbrColumn);

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  RealMatrix(const RealMatrix& M);

  // copy constructor (duplicating all datas)
  //
  // M = matrix to copy
  RealMatrix(Matrix& M);

  // destructor
  //
  ~RealMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  RealMatrix& operator = (const RealMatrix& M);

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

  // get reference of a given matrix element
  //
  // i = line position
  // j = column position
  // return value = reference om matrix elememt
  double& operator () (int i, int j);

  // get reference to a given column
  //
  // i = column position
  // return value = column reference 
  RealVector& operator [] (int i);

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
  friend RealMatrix operator + (const RealMatrix& M1, const RealMatrix& M2);

  // add two matrices where the left one is a real tridiagonal symmetric matrix
  //
  // M1 = left matrix
  // M2 = right matrix
  // return value = sum of the two matrices
  friend RealMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const RealMatrix& M2);

  // add two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M1 = left matrix
  // M2 = right matrix
  // return value = sum of the two matrices
  friend RealMatrix operator + (const RealMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);

  // substract two matrices
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend RealMatrix operator - (const RealMatrix& M1, const RealMatrix& M2);

  // substract two matrices where the left one is a real tridiagonal symmetric matrix
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend RealMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, const RealMatrix& M2);

  // substract two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend RealMatrix operator - (const RealMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);

  // multiply two matrices
  //
  // M1 = first matrix
  // M2 = matrix to multiply to M1
  // return value = product of the two matrices
  friend RealMatrix operator * (const RealMatrix& M1, const RealMatrix& M2);

  // multiply a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend RealMatrix operator * (const RealMatrix& M, double x);

  // multiply a matrix by a real number (left multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend RealMatrix operator * (double x, const RealMatrix& M);

  // multiply a matrix to the right by another matrix without using temporary matrix
  //
  // M = matrix used as multiplicator
  // return value = reference on current matrix
  RealMatrix& Multiply (const RealMatrix& M);

  // multiply a matrix to the right by another matrix without using temporary matrix and in a given range of indices
  // beware the matrix is not resized after multiplication in order the operation to be thread safe
  //
  // M = matrix used as multiplicator
  // startLine = starting line in destination matrix
  // nbrLine = number of lines to multiply
  // return value = reference on current matrix
  RealMatrix& Multiply (const RealMatrix& M, int startLine, int nbrLine);

  // divide a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = division result
  friend RealMatrix operator / (const RealMatrix& M, double x);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  RealMatrix& operator += (const RealMatrix& M);

  // add two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  RealMatrix& operator += (const RealTriDiagonalSymmetricMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  RealMatrix& operator -= (const RealMatrix& M);

  // substract two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  RealMatrix& operator -= (const RealTriDiagonalSymmetricMatrix& M);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  RealMatrix& operator *= (double x);

  // divide a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  RealMatrix& operator /= (double x) ;

  // normalize matrix column vectors
  //
  // return value = reference on current matrix
  RealMatrix& NormalizeColumns ();

  // transpose matrix
  //
  // return value = reference on current matrix
  RealMatrix& Transpose ();

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const RealMatrix& P);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const RealMatrix& P);

};

// get reference of a given matrix element
//
// i = line position
// j = column position
// return value = reference on matrix elememt

inline double& RealMatrix::operator () (int i, int j)
{
  return  this->Columns[j].Components[i];
}

// get reference to a given column
//
// i = column position
// return value = column reference 

inline RealVector& RealMatrix::operator [] (int i)
{
  return this->Columns[i];
}

#endif
