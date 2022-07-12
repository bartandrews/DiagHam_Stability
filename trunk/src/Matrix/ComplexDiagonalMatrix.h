////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of complex diagonal matrix                     //
//                                                                            //
//                        last modification : 09/03/2001                      //
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


#ifndef COMPLEXDIAGONALMATRIX_H
#define COMPLEXDIAGONALMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif

#include <iostream>


using std::ostream;


class RealVector;
class ComplexVector;
class RealMatrix;
class ComplexMatrix;
class RealDiagonalMatrix;

class ComplexDiagonalMatrix : public Matrix
{

  friend class RealVector;
  friend class ComplexVector;
  friend class RealSymmetricMatrix;
  friend class HermitianMatrix;

 protected:

  Complex* DiagonalElements;
  int* DiagonalGarbageFlag;

 public:

  // default constructor
  //
  ComplexDiagonalMatrix();

  // constructor for an empty matrix
  //
  // dimension = matrix dimension
  // zero = true if matrix has to be filled with zeros
  ComplexDiagonalMatrix(int dimension, bool zero = false);

  // constructor from matrix elements (without duplicating datas)
  //
  // diagonal = pointer to diagonal element array
  // dimension = matrix dimension
  ComplexDiagonalMatrix(Complex* diagonal, int dimension);

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  ComplexDiagonalMatrix(const ComplexDiagonalMatrix& M);

  
  // constructor from a full complex matrix
  // M = matrix to copy
  // isDiagonal = returns whether off-diagonal elements were non-zero
  // tolerance = maximal value of off-diagonal elements before isDiagonal is set false
  //
  ComplexDiagonalMatrix(const ComplexMatrix& M, bool &isDiagonal, double tolerance=1e-12);
  

  // destructor
  //
  ~ComplexDiagonalMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  ComplexDiagonalMatrix& operator = (const ComplexDiagonalMatrix& M);

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

  // set a matrix element
  //
  // i = line position
  // j = column position
  // real = real part of new value for matrix element
  // imag = imaginary part of new value for matrix element
  void SetMatrixElement(int i, int j, double real, double imag);


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
  // return value = reference on real part of matrix elememt
  double& operator () (int i, int j);

  // get reference of a given matrix diagonal element
  //
  // i = line position
  // return value = reference om matrix elememt
  Complex& operator [] (int i);


  // get the diagonal elements as an array
  //
  // return value = pointer to the array of diagonal elements
  Complex* GetDiagonalElements();

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

#ifdef USE_HILBERT_SPACE
  // project matrix into a given subspace
  //
  // subspace = reference on subspace structure
  // return value = pointer to projected matrix
  Matrix* Project (SubspaceSpaceConverter& subspace);  
#endif

  // add two matrices
  //
  // M1 = first matrix
  // M2 = second matrix
  // return value = sum of the two matrices
  friend ComplexDiagonalMatrix operator + (const ComplexDiagonalMatrix& M1, 
					const ComplexDiagonalMatrix& M2);

  // substract two matrices
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend ComplexDiagonalMatrix operator - (const ComplexDiagonalMatrix& M1, 
					 const ComplexDiagonalMatrix& M2);


  // multiply a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend ComplexDiagonalMatrix operator * (const ComplexDiagonalMatrix& M, double x);

  // multiply a matrix by a real number (left multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend ComplexDiagonalMatrix operator * (double x, const ComplexDiagonalMatrix& M);

  // multiply two matrices
  //
  // M1 = first matrix
  // M2 = matrix to multiply to M1
  // return value = product of the two matrices
  friend ComplexMatrix operator * (const ComplexMatrix  & M1, const ComplexDiagonalMatrix & M2);
  
  // multiply two matrices
  //
  // M1 = first matrix
  // M2 = matrix to multiply to M1
  // return value = product of the two matrices
  friend ComplexMatrix operator * (const  ComplexDiagonalMatrix & M1, const ComplexMatrix & M2);

  // divide a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = division result
  friend ComplexDiagonalMatrix operator / (const ComplexDiagonalMatrix& M, double x);

  // multiply a matrix by a complex number (right multiplication)
  //
  // M = source matrix
  // x = complex number to use
  // return value = product result
  friend ComplexDiagonalMatrix operator * (const ComplexDiagonalMatrix& M, const Complex &x);

  // multiply a matrix by a real number (left multiplication)
  //
  // M = source matrix
  // x = Complex number to use
  // return value = product result
  friend ComplexDiagonalMatrix operator * (const Complex &x, const ComplexDiagonalMatrix& M);

  // divide a matrix by a Complex number (right multiplication)
  //
  // M = source matrix
  // x = Complex number to use
  // return value = division result
  friend ComplexDiagonalMatrix operator / (const ComplexDiagonalMatrix& M, const Complex &x);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  ComplexDiagonalMatrix& operator += (const RealDiagonalMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  ComplexDiagonalMatrix& operator -= (const RealDiagonalMatrix& M);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  ComplexDiagonalMatrix& operator += (const ComplexDiagonalMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  ComplexDiagonalMatrix& operator -= (const ComplexDiagonalMatrix& M);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  ComplexDiagonalMatrix& operator *= (double x);

  // divide a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  ComplexDiagonalMatrix& operator /= (double x) ;

    // multiply a matrix by a complex number
  //
  // x = real number to use
  // return value = reference on current matrix
  ComplexDiagonalMatrix& operator *= (Complex &x);

  // divide a matrix by a complex number
  //
  // x = real number to use
  // return value = reference on current matrix
  ComplexDiagonalMatrix& operator /= (Complex &x) ;


  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // evaluate matrix trace
  //
  // return value = matrix trace 
  Complex Trace ();

  // evaluate matrix determinant
  //
  // return value = matrix determinant 
  Complex Determinant ();

  // Sort Matrix such that diagnonal elements are sort in decreasing order
  //
  // normSort = sort with respect to the norm instead of the real part
  // normError = when sort with respect to the norm, consider two norms to be identical when their difference is lower than normError, 
  //             then sort values with the same norm with respect to their phase (0 if no sort has to be applied on the phase)
  // return value = reference on current Matrix
  ComplexDiagonalMatrix& SortMatrixDownOrder(bool normSort = false, double normError = 0.0);
  
  // Sort Matrix such that diagnonal elements are sort in decreasing order
  // and apply corresponding transformation to column of a given complex matrix 
  //
  // matrix = matrix on which transformation has to be applied
  // normSort = sort with respect to the norm instead of the real part
  // normError = when sort with respect to the norm, consider two norms to be identical when their difference is lower than normError, 
  //             then sort values with the same norm with respect to their phase (0 if no sort has to be applied on the phase)
  // return value = reference on current Matrix
  ComplexDiagonalMatrix& SortMatrixDownOrder(ComplexMatrix& matrix, bool normSort = false, double normError = 0.0);

  // Sort Matrix such that diagnonal elements are sort in increasing order
  //
  // normSort = sort with respect to the norm instead of the real part
  // normError = when sort with respect to the norm, consider two norms to be identical when their difference is lower than normError, 
  //             then sort values with the same norm with respect to their phase (0 if no sort has to be applied on the phase)
  // return value = reference on current Matrix
  ComplexDiagonalMatrix& SortMatrixUpOrder(bool normSort = false, double normError = 0.0);
  
  // Sort Matrix such that diagnonal elements are sort in increasing order
  // and apply corresponding transformation to column of a given complex matrix 
  //
  // matrix = matrix on which transformation has to be applied
  // normSort = sort with respect to the norm instead of the real part
  // normError = when sort with respect to the norm, consider two norms to be identical when their difference is lower than normError, 
  //             then sort values with the same norm with respect to their phase (0 if no sort has to be applied on the phase)
  // return value = reference on current Matrix
  ComplexDiagonalMatrix& SortMatrixUpOrder(ComplexMatrix& matrix, bool normSort = false, double normError = 0.0);

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const ComplexDiagonalMatrix& P);

#ifdef USE_OUTPUT

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexDiagonalMatrix& P);

#endif

};

// get the diagonal elements as an array
//
// return value = pointer to the array of diagonal elements

inline Complex* ComplexDiagonalMatrix::GetDiagonalElements()
{
  return this->DiagonalElements;
}

#endif
