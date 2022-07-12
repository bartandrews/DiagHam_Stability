////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of real symmetric matrix                      //
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


#ifndef REALSYMMETRICMATRIX_H
#define REALSYMMETRICMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif

#include <iostream>
#include <fstream>


using std::ostream;
using std::ofstream;
using std::ifstream;


class RealMatrix;
class RealDiagonalMatrix;
class BlockDiagonalMatrix;


class RealSymmetricMatrix : public Matrix
{

  friend class RealVector;
  friend class ComplexVector;

 protected:

  double* DiagonalElements;
  int* DiagonalGarbageFlag;

  double* OffDiagonalElements;
  int* OffDiagonalGarbageFlag;

  int Increment;

 public:

  // default constructor
  //
  RealSymmetricMatrix();

  // constructor for an empty matrix
  //
  // dimension = matrix dimension
  // zero = true if matrix has to be filled with zeros
  RealSymmetricMatrix(int dimension, bool zero = false);

  // constructor from matrix elements (without duplicating datas)
  //
  // diagonal = pointer to diagonal element array
  // offDiagonal = pointer to off-diagonal element array
  // dimension = matrix dimension
  RealSymmetricMatrix(double* diagonal, double* offDiagonal, int dimension);

  // constructor from a real matrix Q (new matrix = Qt * Q)
  //
  RealSymmetricMatrix(const RealMatrix& Q);

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  RealSymmetricMatrix(const RealSymmetricMatrix& M);

  // copy constructor from any matrix (only keeping real part of elements of and above the diagonal, duplicating datas)
  //
  // M = matrix to copy
  RealSymmetricMatrix(Matrix& M);

  // destructor
  //
  ~RealSymmetricMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  RealSymmetricMatrix& operator = (const RealSymmetricMatrix& M);

  // assignement from  a real tridiagonal symmetric matrix (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix  
  RealSymmetricMatrix& operator = (const RealTriDiagonalSymmetricMatrix& M);

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

  // copy matrix
  //
  // M = matrix to copy
  // return value = refence on current matrix
  RealSymmetricMatrix& Copy (RealSymmetricMatrix& M);

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
  friend RealSymmetricMatrix operator + (const RealSymmetricMatrix& M1, 
					 const RealSymmetricMatrix& M2);

  // add two matrices where the left one is a real tridiagonal symmetric matrix
  //
  // M1 = left matrix
  // M2 = right matrix
  // return value = sum of the two matrices
  friend RealSymmetricMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, 
					 const RealSymmetricMatrix& M2);

  // add two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M1 = left matrix
  // M2 = right matrix
  // return value = sum of the two matrices
  friend RealSymmetricMatrix operator + (const RealSymmetricMatrix& M1, 
					 const RealTriDiagonalSymmetricMatrix& M2);

  // substract two matrices
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend RealSymmetricMatrix operator - (const RealSymmetricMatrix& M1, 
					 const RealSymmetricMatrix& M2);

  // substract two matrices where the left one is a real tridiagonal symmetric matrix
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend RealSymmetricMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, 
					 const RealSymmetricMatrix& M2);

  // substract two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend RealSymmetricMatrix operator - (const RealSymmetricMatrix& M1,
					 const RealTriDiagonalSymmetricMatrix& M2);

  // multiply a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend RealSymmetricMatrix operator * (const RealSymmetricMatrix& M, double x);

  // multiply a matrix by a real number (left multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend RealSymmetricMatrix operator * (double x, const RealSymmetricMatrix& M);

  // divide a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = division result
  friend RealSymmetricMatrix operator / (const RealSymmetricMatrix& M, double x);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  RealSymmetricMatrix& operator += (const RealSymmetricMatrix& M);

  // add two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  RealSymmetricMatrix& operator += (const RealTriDiagonalSymmetricMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  RealSymmetricMatrix& operator -= (const RealSymmetricMatrix& M);

  // substract two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  RealSymmetricMatrix& operator -= (const RealTriDiagonalSymmetricMatrix& M);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  RealSymmetricMatrix& operator *= (double x);

  // divide a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  RealSymmetricMatrix& operator /= (double x) ;

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  double MatrixElement (RealVector& V1, RealVector& V2);

  // conjugate an hermitian matrix with an unitary matrix (Ut M U)
  //
  // UnitaryM = unitary matrix to use
  // return value = pointer to conjugated matrix
  Matrix* Conjugate(RealMatrix& UnitaryM);

  // conjugate an hermitian matrix with an unitary matrix (Ut M U)
  //
  // UnitaryM = unitary matrix to use
  // return value = pointer to conjugated matrix
  Matrix* Conjugate(BlockDiagonalMatrix& UnitaryM);

  // conjugate a block of the matrix with an unitary matrix (Ut M U)
  //
  // UnitaryM = unitary matrix to use
  // sourcePosition = index of the row where the block to conjugate starts
  // destinationPosition = index of the row where the conjugated block has to be stored
  // matrix = matrix where result has to be stored
  void Conjugate(RealMatrix& UnitaryM, int sourcePosition, int destinationPosition, 
		 RealSymmetricMatrix& matrix);

  // conjugate a block of the matrix (in the upper diagonal part) with two matrix matrix (Vt M U)
  //
  // UnitaryMl = unitary matrix to use at the left hand side
  // UnitaryMr = unitary matrix to use at the right hand side
  // sourceRowIndex = index of the row where the block to conjugate starts
  // sourceColumnIndex = index of the column where the block to conjugate starts
  // destinationRowIndex = index of the row where the conjugated block has to be stored
  // destinationColumnIndex = index of the column where the conjugated block has to be stored
  // matrix = matrix where result has to be stored
  void Conjugate(RealMatrix& UnitaryMl, RealMatrix& UnitaryMr, int sourceRowIndex, int sourceColumnIndex, 
		 int destinationRowIndex, int destinationColumnIndex, RealSymmetricMatrix& matrix);

  // add the symmetric matrix m1.m2.m3^t.m4^t + m4.m3.m2^t.m1^t to current symmetric matrix
  //
  // m1 = first matrix
  // m2 = second matrix
  // m3 = third matrix
  // m4 = fourth matrix
  // coefficient = optional global multiplicative factor in front of m1.m2.m3^t.m4^t + m4.m3.m2^t.m1^t
  RealSymmetricMatrix& AddAAAtAt(RealMatrix& m1, RealMatrix& m2, RealMatrix& m3, RealMatrix& m4, double coefficient = 1.0);

  // add the symmetric matrix m1.m2.m3^t.m4^t + m4.m3.m2^t.m1^t to current symmetric matrix within a given range of indices
  //
  // m1 = first matrix
  // m2 = second matrix
  // m3 = third matrix
  // m4 = fourth matrix
  // rowIndex = row index of the first element to add
  // columnIndex = column index of the first element to add
  // nbrElement = number of element to add (starting from first element knowing that line i contains n - i + 1)
  // coefficient = optional global multiplicative factor in front of m1.m2.m3^t.m4^t + m4.m3.m2^t.m1^t
  RealSymmetricMatrix& AddAAAtAt(RealMatrix& m1, RealMatrix& m2, RealMatrix& m3, RealMatrix& m4, int rowIndex, int columnIndex,
				 int nbrElement, double coefficient = 1.0);

  // evaluate matrix trace
  //
  // return value = matrix trace 
  double Tr ();

  // evaluate matrix determinant
  //
  // return value = matrix determinant 
  double Det ();

  // Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step
  //
  // dimension = maximum iteration number
  // M = reference on complex tridiagonal hermitian matrix where result has to be stored
  // V1 = reference on complex vector used as first vector (will contain last produced vector at the end)
  // return value = reference on complex tridiagonal hermitian matrix
  RealTriDiagonalSymmetricMatrix& Lanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, RealVector& V1);

  // Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step
  //
  // dimension = maximum iteration number
  // M = reference on real tridiagonal symmetric matrix where result has to be stored
  // Q = matrix where new orthonormalized base will be stored (first column is used as first vector)
  // return value = reference on complex tridiagonal hermitian matrix
  RealTriDiagonalSymmetricMatrix& Lanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, RealMatrix& Q);

  // Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step, if during process a 
  // null vector appears then new random vector is evaluated
  //
  // dimension = maximum iteration number
  // M = reference on real tridiagonal symmetric matrix where result has to be stored
  // Q = matrix where new orthonormalized base will be stored (first column is used as first vector)
  // err = absolute error on vector norm
  // return value = reference on complex tridiagonal hermitian matrix
  RealTriDiagonalSymmetricMatrix& OrthoLanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, RealMatrix& Q, 
						double err = 0.00000001);

  // Tridiagonalize a real symmetric matrix using Householder algorithm (modifying current matrix)
  //
  // M = reference on real tridiagonal symmetric matrix where result has to be stored
  // err = absolute error on matrix element
  // return value = reference on real tridiagonal symmetric matrix
  RealTriDiagonalSymmetricMatrix& Householder (RealTriDiagonalSymmetricMatrix& M, double err);

  // Tridiagonalize a real symmetric matrix using Householder algorithm and evaluate transformation matrix  (modifying current matrix)
  //
  // M = reference on real tridiagonal symmetric matrix where result has to be stored
  // err = absolute error on matrix element
  // Q = matrix where transformation matrix has to be stored
  // return value = reference on real tridiagonal symmetric matrix
  RealTriDiagonalSymmetricMatrix& Householder (RealTriDiagonalSymmetricMatrix& M, double err, RealMatrix& Q);

  // Diagonalize a real symmetric matrix (modifying current matrix)
  //
  // M = reference on real diagonal matrix where result has to be stored
  // err = absolute error on matrix element
  // maxIter = maximum number of iteration to fund an eigenvalue
  // return value = reference on real tridiagonal symmetric matrix
  RealDiagonalMatrix& Diagonalize (RealDiagonalMatrix& M, double err = 1e-7, int maxIter = 50);

  // Diagonalize a real symmetric matrix and evaluate transformation matrix (modifying current matrix)
  //
  // M = reference on real diagonal matrix where result has to be stored
  // Q = matrix where transformation matrix has to be stored
  // err = absolute error on matrix element
  // maxIter = maximum number of iteration to fund an eigenvalue
  // return value = reference on real tridiagonal symmetric matrix
  RealDiagonalMatrix& Diagonalize (RealDiagonalMatrix& M, RealMatrix& Q, double err = 1e-7, int maxIter = 50);

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
  friend ostream& operator << (ostream& Str, const RealSymmetricMatrix& P);

#ifdef USE_OUTPUT

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const RealSymmetricMatrix& P);

#endif

  // output file stream overload
  //
  // file = reference on output file stream
  // matrix = reference on matrix to save
  // return value = reference on output file stream
//  friend ofstream& operator << (ofstream& file, const RealSymmetricMatrix& matrix);

  // input file stream overload
  //
  // file = reference on output file stream
  // matrix = reference on matrix to load
  // return value = reference on output file stream
  friend ifstream& operator >> (ifstream& file, RealSymmetricMatrix& matrix);

};

#endif
