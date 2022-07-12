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
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Vector/ComplexVector.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class ComplexUpperTriangularMatrix;


class ComplexMatrix : public Matrix
{

  friend class ComplexSkewSymmetricMatrix;
  friend class ComplexUpperTriangularMatrix;
  friend class HermitianMatrix;
  friend class RealVector;
  friend class ComplexVector;
  friend class RealTriDiagonalSymmetricMatrix;
  friend class SingleParticle;

 protected:

  ComplexVector* Columns; 
  GarbageFlag Flag;

 public:

  // default constructor
  //
  ComplexMatrix();

  // constructor for an empty matrix
  //
  // nbrRow = number of rows
  // nbrColumn = number of columns
  // zero = tue if matrix elements have to be set to zero
  ComplexMatrix(int nbrRow, int nbrColumn, bool zero = false);

  // constructor from matrix elements (without duplicating datas)
  //
  // columns = pointer an array of vector
  // nbrColumn = number of columns
  ComplexMatrix(ComplexVector* columns, int nbrColumn);

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  ComplexMatrix(const ComplexMatrix& M);

  // copy constructor (duplicating all datas)
  //
  // M = matrix to copy
  ComplexMatrix(Matrix& M);

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

  // set a matrix element
  //
  // i = line position
  // j = column position
  // real = new real value for matrix element
  // imag = new imaginary value for matrix element
  void SetMatrixElement(int i, int j, double real, double imag);

  // get a matrix element (real part if complex)
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
  void GetMatrixElement(int i, int j, double& x) const;

  // get a matrix element
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
  void GetMatrixElement(int i, int j, Complex& x) const;

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

  // get reference to a given column
  //
  // i = column position
  // return value = column reference 
  ComplexVector& operator [] (int i);

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

  // multiply a complex matrix with a complex upper triangular matrix
  //
  // m1 = complex matrix
  // m2 = complex upper triangular matrix
  // return value = product result
  friend ComplexMatrix operator * (ComplexMatrix& m1, const ComplexUpperTriangularMatrix& m2);

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

  // get adjoint (hermitian conjugate) matrix 
  //
  // return value = reference on modified matrix
  ComplexMatrix GetAdjoint();

  // evaluate matrix determinant (skrewing up matrix elements)
  //
  // return value = matrix determinant 
  Complex Determinant ();

  // evaluate permanent associated to the (square) matrix using Ryser algorithm
  //
  // return value = permanent associated to the matrix
  Complex Permanent();                                                                                                                                     

  // evaluate minor develomment of permanent associated to the (square) matrix using Ryser algorithm
  //
  // column = index of the column from which permnanent will developped
  // minors = reference on an array where minors will be stored
  void PermanentMinorDevelopment(int column, Complex*& minors);
  
  // evaluate permanent associated to the (square) matrix using Ryser algorithm and precalculation array (faster)
  //
  // changeBit = array indicating which bit is changed at the i-th iteration of the Gray code
  // changeBitSign = array with -1 if the changed bit is from 1 to 0, +1 either
  // return value = permanent associated to the matrix
  Complex FastPermanent(int* changeBit, int* changeBitSign);

  // evaluate minor develomment of permanent associated to the (square) matrix using Ryser algorithm and precalculation array (faster)
  //
  // changeBit = array indicating which bit is changed at the i-th iteration of the Gray code
  // changeBitSign = array with -1 if the changed bit is from 1 to 0, +1 either
  // column = index of the column from which permnanent will developped
  // minors = reference on an array where minors will be stored
  void FastPermanentMinorDevelopment(int* changeBit, int* changeBitSign, int column, Complex*& minors);
  
  // evaluate precalculation array needed for the fast permanent calculation
  //
  // changeBit = reference on the array indicating which bit is changed at the i-th iteration of the Gray code
  // changeBitSign = reference on array with -1 if the changed bit is from 1 to 0, +1 either
  // minor = flag that indicated if precalculation will be used for minor development
  void EvaluateFastPermanentPrecalculationArray(int*& changeBit, int*& changeBitSign, bool minor = false);

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const ComplexMatrix& P);

#ifdef USE_OUTPUT

 // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexMatrix& P);

#endif

  // Diagonalize an hermitian matrix (modifying current matrix)
  //
  // M = reference on real diagonal matrix where result has to be stored
  // err = absolute error on matrix element
  // maxIter = maximum number of iteration to fund an eigenvalue
  // return value = reference on real tridiagonal symmetric matrix
  ComplexDiagonalMatrix& Diagonalize (ComplexDiagonalMatrix& M);

  // Diagonalize an hermitian matrix and evaluate transformation matrix (modifying current matrix)
  //
  // M = reference on real diagonal matrix where result has to be stored
  // Q = matrix where transformation matrix has to be stored
  // err = absolute error on matrix element
  // maxIter = maximum number of iteration to fund an eigenvalue
  // return value = reference on real tridiagonal symmetric matrix
  ComplexDiagonalMatrix& Diagonalize (ComplexDiagonalMatrix& M, ComplexMatrix& Q);

#ifdef __LAPACK__

  // calculate a determinant using the LAPACK library (conserving current matrix)
  //
  Complex LapackDeterminant ();

    // Diagonalize a hermitian matrix using the LAPACK library (modifying current matrix)
  //
  // M = reference on real diagonal matrix of eigenvalues
  // err = absolute error on matrix element
  // maxIter = maximum number of iteration to fund an eigenvalue
  // return value = reference on real matrix consisting of eigenvalues
  ComplexDiagonalMatrix& LapackDiagonalize (ComplexDiagonalMatrix& M);

  // Diagonalize a hermitian matrix and evaluate transformation matrix using the LAPACK library (modifying current matrix)
  //
  // M = reference on real diagonal matrix of eigenvalues
  // Q = matrix where transformation matrix has to be stored
  // err = absolute error on matrix element
  // maxIter = maximum number of iteration to fund an eigenvalue
  // return value = reference on real matrix consisting of eigenvalues
  ComplexDiagonalMatrix& LapackDiagonalize (ComplexDiagonalMatrix& M, ComplexMatrix& Q);

  // reduce a complex matrix to its Schur form S
  //
  // M = reference on real diagonal matrix of eigenvalues
  // Q = matrix where transformation matrix has to be stored
  // S = matrix where Schur form of matrix has to be stored
  // return value = reference o1n real matrix consisting of eigenvalues
  ComplexDiagonalMatrix& LapackSchurForm (ComplexDiagonalMatrix& M, ComplexMatrix& Q, ComplexMatrix &S);

 private:
  int LapackWorkAreaDimension;
  doublecomplex *LapackMatrix;
  doublecomplex *LapackEVMatrix;
  doublecomplex *LapackWorkingArea;
  double *LapackRealWorkingArea;  

#endif

};

// get reference to a given column
//
// i = column position
// return value = column reference 

inline ComplexVector& ComplexMatrix::operator [] (int i)
{
  return this->Columns[i];
}

// get a matrix element (real part if complex)
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

inline void ComplexMatrix::GetMatrixElement(int i, int j, double& x) const
{
  x = this->Columns[j].Re(i);
}

// get a matrix element
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

inline void ComplexMatrix::GetMatrixElement(int i, int j, Complex& x) const
{
  x = this->Columns[j][i];
}

#endif
