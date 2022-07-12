////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of real band-diagonal symmetric matrix                //
//                                                                            //
//                        last modification : 16/03/2005                      //
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


#ifndef REALBANDDIAGONALSYMMETRICMATRIX_H
#define REALBANDDIAGONALSYMMETRICMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif
#include "Vector/RealVector.h"
#ifdef USE_POLYNOMIAL
#include "Polynomial/Polynomial.h"
#endif
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class Complex;
class RealTriDiagonalSymmetricMatrix;


class RealBandDiagonalSymmetricMatrix : public Matrix
{

 protected:

  // array which contains diagonal elements
  double* DiagonalElements;
  // array which contains upper off-diagonal elements (first index is used as row index)
  double** UpperOffDiagonalElements;

  // number of bands in used in the upper part of the matrix
  int NbrBands;
  // total number of bands in the upper part of the matrix
  int TrueNbrBands;

  // garbage flag used to avoid data duplication
  GarbageFlag Flag;

  // dummy variable used if an element outside diagonal or upper(lower) diagonal is requested
  double Dummy;

 public:

  // default constructor
  //
  RealBandDiagonalSymmetricMatrix();

  // constructor for an empty matrix
  //
  // dimension = matrix dimension
  // nbrBands = number of bands in the upper part of the matrix
  // zero = true if matrix has to be filled with zeros
  RealBandDiagonalSymmetricMatrix(int dimension, int nbrBands, bool zero = false);

  // constructor from matrix elements (without duplicating datas)
  //
  // diagonal = pointer to diagonal element array
  // upperOffDiagonal = pointer to the array which contains upper off-diagonal elements (second index is used as row index)
  // dimension = matrix dimension
  // nbrBands = number of bands in the upper part of the matrix
  RealBandDiagonalSymmetricMatrix(double* diagonal, double** upperOffDiagonal, int dimension, int nbrBands);

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  RealBandDiagonalSymmetricMatrix(const RealBandDiagonalSymmetricMatrix& M);

  // destructor
  //
  ~RealBandDiagonalSymmetricMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  RealBandDiagonalSymmetricMatrix& operator = (const RealBandDiagonalSymmetricMatrix& M);

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

  // get a matrix element (real part if complex)
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
  void GetMatrixElement(int i, int j, double& x) const;

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

  // return refernce on real part of a given matrix element
  //
  // i = line position
  // j = column position
  // return value = reference on real part
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

  // Resize matrix and change the number of bands
  //
  // nbrRow = new number of rows
  // nbrColumn = new number of columns
  // nbrBands = new number of bands
  void Resize (int nbrRow, int nbrColumn, int nbrBands);

  // copy matrix
  //
  // M = matrix to copy
  // return value = refence on current matrix
  RealBandDiagonalSymmetricMatrix& Copy (RealBandDiagonalSymmetricMatrix& M);

  // add two matrices
  //
  // M1 = first matrix
  // M2 = second matrix
  // return value = sum of the two matrices
  friend RealBandDiagonalSymmetricMatrix operator + (const RealBandDiagonalSymmetricMatrix& M1, const RealBandDiagonalSymmetricMatrix& M2);

  // substract two matrices
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend RealBandDiagonalSymmetricMatrix operator - (const RealBandDiagonalSymmetricMatrix& M1, const RealBandDiagonalSymmetricMatrix& M2);

  // multiply a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend RealBandDiagonalSymmetricMatrix operator * (const RealBandDiagonalSymmetricMatrix& M, double x);

  // multiply a matrix by a real number (left multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend RealBandDiagonalSymmetricMatrix operator * (double x, const RealBandDiagonalSymmetricMatrix& M);

  // divide a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = division result
  friend RealBandDiagonalSymmetricMatrix operator / (const RealBandDiagonalSymmetricMatrix& M, double x);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  RealBandDiagonalSymmetricMatrix& operator += (const RealBandDiagonalSymmetricMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  RealBandDiagonalSymmetricMatrix& operator -= (const RealBandDiagonalSymmetricMatrix& M);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  RealBandDiagonalSymmetricMatrix& operator *= (double x);

  // divide a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  RealBandDiagonalSymmetricMatrix& operator /= (double x) ;

  // get a matrix element 
  // 
  // i = Row number
  // j = Column number
  // return value = matrix element M_(i,j)
  double GetElement(int i, int j);

  // access to i-th diagonal element
  // 
  // i = position 
  // return value = reference on i-th diagonal element
  double& DiagonalElement(int i);

  // evaluate matrix trace
  //
  // return value = matrix trace 
  double Tr ();

  // evaluate matrix determinant
  //
  // return value = matrix determinant 
  double Det ();

  // Tridiagonalize a real band symmetric matrix using Rutishauer-Schwarz  (modifying current matrix)
  //
  // M = reference on real tridiagonal symmetric matrix where result has to be stored
  // err = absolute error on matrix element
  // return value = reference on real tridiagonal symmetric matrix
  RealTriDiagonalSymmetricMatrix& Tridiagonalize (RealTriDiagonalSymmetricMatrix& M, double err);

  // Tridiagonalize a real band symmetric matrix using Rutishauer-Schwarz and evaluate transformation matrix  (modifying current matrix)
  //
  // M = reference on real tridiagonal symmetric matrix where result has to be stored
  // err = absolute error on matrix element
  // Q = matrix where transformation matrix has to be stored
  // return value = reference on real tridiagonal symmetric matrix
  RealTriDiagonalSymmetricMatrix& Tridiagonalize (RealTriDiagonalSymmetricMatrix& M, double err, RealMatrix& Q);

#ifdef __LAPACK__

  // Diagonalize a real symmetric matrix using the LAPACK library (modifying current matrix)
  //
  // M = reference on real diagonal matrix where result has to be stored
  // err = absolute error on matrix element
  // maxIter = maximum number of iteration to fund an eigenvalue
  // return value = reference on real tridiagonal symmetric matrix
  RealDiagonalMatrix& LapackDiagonalize (RealDiagonalMatrix& M, double err = 1e-7, int maxIter = 50);

  // Diagonalize a real symmetric matrix and evaluate transformation matrix using the LAPACK library (modifying current matrix)
  //
  // M = reference on real diagonal matrix where result has to be stored
  // Q = matrix where transformation matrix has to be stored
  // err = absolute error on matrix element
  // maxIter = maximum number of iteration to fund an eigenvalue
  // return value = reference on real tridiagonal symmetric matrix
  RealDiagonalMatrix& LapackDiagonalize (RealDiagonalMatrix& M, RealMatrix& Q, double err = 1e-7, int maxIter = 50);

#endif

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const RealBandDiagonalSymmetricMatrix& P);

#ifdef USE_OUTPUT

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const RealBandDiagonalSymmetricMatrix& P);

#endif

};

#endif
