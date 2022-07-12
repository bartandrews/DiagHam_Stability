////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of real upper triangular matrix                  //
//                                                                            //
//                        last modification : 07/01/2003                      //
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


#ifndef REALUPPERTRIANGULARMATRIX_H
#define REALUPPERTRIANGULARMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Output/MathematicaOutput.h"

#include <iostream>
#include <fstream>


using std::ostream;
using std::ofstream;
using std::ifstream;


class RealMatrix;
class BlockDiagonalMatrix;


class RealUpperTriangularMatrix : public Matrix
{

  friend class RealVector;
  friend class ComplexVector;

 protected:

  double* DiagonalElements;
  int* DiagonalGarbageFlag;

  double* OffDiagonalElements;
  int* OffDiagonalGarbageFlag;

  int Increment;

  // dummy variable whose reference is send when an element of the lower part of the matrix is asked (initialize to 0)
  double Dummy;

 public:

  // default constructor
  //
  RealUpperTriangularMatrix();

  // constructor for an empty matrix
  //
  // dimension = matrix dimension
  // zero = true if matrix has to be filled with zeros
  RealUpperTriangularMatrix(int dimension, bool zero = false);

  // constructor from matrix elements (without duplicating datas)
  //
  // diagonal = pointer to diagonal element array
  // offDiagonal = pointer to off-diagonal element array
  // dimension = matrix dimension
  RealUpperTriangularMatrix(double* diagonal, double* offDiagonal, int dimension);

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  RealUpperTriangularMatrix(const RealUpperTriangularMatrix& M);

  // copy constructor from a real tridiagonal symmetric matrix (without duplicating diagonal elements)
  //
  // M = matrix to copy
  RealUpperTriangularMatrix(const RealTriDiagonalSymmetricMatrix& M);

  // destructor
  //
  ~RealUpperTriangularMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  RealUpperTriangularMatrix& operator = (const RealUpperTriangularMatrix& M);

  // assignement from  a real tridiagonal symmetric matrix (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix  
  RealUpperTriangularMatrix& operator = (const RealTriDiagonalSymmetricMatrix& M);

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

  // get reference of a given diagonal matrix element
  //
  // i = line position
  // return value = reference on matrix diagonal elememt
  double& operator [] (int i);

  // get reference of a given matrix element
  //
  // i = line position
  // j = column position
  // return value = reference om matrix elememt
  double& operator () (int i, int j);

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
  friend RealUpperTriangularMatrix operator + (const RealUpperTriangularMatrix& M1, 
					       const RealUpperTriangularMatrix& M2);

  // substract two matrices
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend RealUpperTriangularMatrix operator - (const RealUpperTriangularMatrix& M1, 
					       const RealUpperTriangularMatrix& M2);

  // multiply a real matrix with a real upper triangular matrix
  //
  // m1 = real matrix
  // m2 = real upper triangular matrix
  // return value = product result
  friend RealMatrix operator * (RealMatrix& m1, const RealUpperTriangularMatrix& m2);

  // multiply a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend RealUpperTriangularMatrix operator * (const RealUpperTriangularMatrix& M, double x);

  // multiply a matrix by a real number (left multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend RealUpperTriangularMatrix operator * (double x, const RealUpperTriangularMatrix& M);

  // divide a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = division result
  friend RealUpperTriangularMatrix operator / (const RealUpperTriangularMatrix& M, double x);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  RealUpperTriangularMatrix& operator += (const RealUpperTriangularMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  RealUpperTriangularMatrix& operator -= (const RealUpperTriangularMatrix& M);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  RealUpperTriangularMatrix& operator *= (double x);

  // divide a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  RealUpperTriangularMatrix& operator /= (double x) ;

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  double MatrixElement (RealVector& V1, RealVector& V2);

  // evaluate matrix trace
  //
  // return value = matrix trace 
  double Tr ();

  // evaluate matrix determinant
  //
  // return value = matrix determinant 
  double Det ();

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const RealUpperTriangularMatrix& P);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const RealUpperTriangularMatrix& P);

  // output file stream overload
  //
  // file = reference on output file stream
  // matrix = reference on matrix to save
  // return value = reference on output file stream
//  friend ofstream& operator << (ofstream& file, const RealUpperTriangularMatrix& matrix);

  // input file stream overload
  //
  // file = reference on output file stream
  // matrix = reference on matrix to load
  // return value = reference on output file stream
  friend ifstream& operator >> (ifstream& file, RealUpperTriangularMatrix& matrix);

};

// get reference of a given diagonal matrix element
//
// i = line position
// return value = reference on matrix diagonal elememt

inline double& RealUpperTriangularMatrix::operator [] (int i)
{
  return this->DiagonalElements[i];
}

#endif
