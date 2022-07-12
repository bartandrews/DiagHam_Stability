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
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "MathTools/Complex.h"

#include <iostream>
#include <cassert>

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
  friend class RealDiagonalMatrix;
  friend class ComplexMatrix;

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

  // constructor for one dimensional array
  //
  // array = one dimensional array where the matrix elements are stored
  // nbrRow = number of rows
  // nbrColumn = number of columns
  // columnOrder = elements in array are ordered column-wise  (all components of the first column, then all components of the second column,...)
  RealMatrix(double* array, int nbrRow, int nbrColumn, bool columnOrder = true);

#ifdef __MPI__
  // constructor from informations sent using MPI
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts or sends the vector
  // broadcast = true if the vector is broadcasted
  RealMatrix(MPI::Intracomm& communicator, int id, bool broadcast = true);
#endif

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

  // copy a matrix into another (duplicating data)
  //
  // matrix = matrix to copy
  // return value = reference on current matrix
  RealMatrix& Copy (RealMatrix& matrix);

  // get a matrix element (real part if complex)
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
  void GetMatrixElement(int i, int j, double& x) const;

  // get a matrix element (real part if complex)
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
  void GetMatrixElement(int i, int j, Complex& x) const;

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

  // Set all entries in matrix to zero
  //
  void ClearMatrix ();

  // set matrix to identity 
  //
  void SetToIdentity();

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

  // multiply two matrices
  //
  // M1 = first matrix
  // M2 = matrix to multiply to M1
  // return value = product of the two matrices
  friend RealMatrix operator * (const RealMatrix  & M1, const RealDiagonalMatrix & M2);

  // multiply two matrices
  //
  // M1 = first matrix
  // M2 = matrix to multiply to M1
  // return value = product of the two matrices
  friend RealMatrix operator * (const  RealDiagonalMatrix & M1, const RealMatrix & M2);

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

  // multiply two matrices and add the result to the current one
  //
  // M1 = left matrix 
  // M2 = right matrix 
  // return value = reference on current matrix
  RealMatrix& AddMultiply (const RealMatrix  & M1, const RealMatrix & M2);

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

  // divide a matrix by a real diagonal matrix (all entries must be non zero)
  //
  // M = source matrix
  // x = real number to use
  // return value = division result
  friend RealMatrix operator / (const RealMatrix& M1, const RealDiagonalMatrix& M2);



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

  // orthonormalize matrix column vectors
  //
  // return value = reference on current matrix
  RealMatrix& OrthoNormalizeColumns ();

  // orthonormalize matrix column vectors, computing the transformation matrix to the new orthonormal basis
  //
  // transformation= reference on the transformation matrix
  // return value = reference on current matrix
  RealMatrix& OrthoNormalizeColumns (RealMatrix& transformation);

  // transpose matrix
  //
  // return value = reference on current matrix
  RealMatrix& Transpose ();

  // duplicate and transpose a matrix
  //
  // return value = transposed matrix
  RealMatrix DuplicateAndTranspose();
    
  // evaluate matrix trace
  //
  // return value = matrix trace 
  virtual double Tr ();

  // evaluate matrix determinant (skrewing up matrix elements)
  //
  // return value = matrix determinant 
  double Determinant ();

  // evaluate permanent associated to the (square) matrix using Ryser algorithm
  //
  // return value = permanent associated to the matrix
  double Permanent();

  // compute singular value decomposition U D V^t
  // 
  // uMatrix = reference on the U matrix
  // vMatrix = reference on the V matrix
  // truncatedUVFlag = if false, set JOBZ = 'A' (returns full U, V matrices)
  // return value = pointer on the diagonal elements of D
  double* SingularValueDecomposition(RealMatrix& uMatrix, RealMatrix& vMatrix, bool truncatedUVFlag = true);

  // compute singular value decomposition U D V^t
  // 
  // uMatrix = reference on the U matrix
  // diagonal = reference on the diagonal D matrix
  // vMatrix = reference on the V matrix
  void SingularValueDecomposition(RealMatrix& uMatrix, RealDiagonalMatrix& diagonal, RealMatrix& vMatrix, bool truncatedUVFlag = true);

  // compute the diagonal part of the singular value decomposition U D V^t
  // 
  // return value = pointer on the diagonal elements of D
  double* SingularValueDecomposition();

  // Diagonalize a real matrix using the LAPACK library
  //
  // M = reference on complex diagonal matrix where result has to be stored
  // leftFlag = compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
  // return value = reference on complex diagonal matrix
  ComplexDiagonalMatrix& LapackDiagonalize (ComplexDiagonalMatrix& M, bool leftFlag = false);

  // Diagonalize a real matrix and evaluate the left eigenstates using the LAPACK library
  //
  // M = reference on complex diagonal matrix where result has to be stored
  // Q = matrix where transformation matrix has to be stored
  // leftFlag = compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
  // return value = reference on complex diagonal matrix
  ComplexDiagonalMatrix& LapackDiagonalize (ComplexDiagonalMatrix& M, ComplexMatrix& Q, bool leftFlag = false);

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const RealMatrix& P);

#ifdef USE_OUTPUT

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const RealMatrix& P);

#endif

#ifdef __MPI__

  // send a matrix to a given MPI process
  // 
  // communicator = reference on the communicator to use
  // id = id of the destination MPI process
  // return value = reference on the current matrix
  virtual Matrix& SendMatrix(MPI::Intracomm& communicator, int id);

  // broadcast a matrix to all MPI processes associated to the same communicator
  // 
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the matrix
  // return value = reference on the current matrix
  virtual Matrix& BroadcastMatrix(MPI::Intracomm& communicator,  int id);

  // receive a matrix from a MPI process
  // 
  // communicator = reference on the communicator to use 
  // id = id of the source MPI process
  // return value = reference on the current matrix
  virtual Matrix& ReceiveMatrix(MPI::Intracomm& communicator, int id);

  // add current matrix to the current matrix of a given MPI process
  // 
  // communicator = reference on the communicator to use 
  // id = id of the destination MPI process
  // return value = reference on the current matrix
  virtual Matrix& SumMatrix(MPI::Intracomm& communicator, int id);

  // reassemble matrix from a scattered one
  // 
  // communicator = reference on the communicator to use 
  // id = id of the destination MPI process
  // return value = reference on the current matrix
  virtual Matrix& ReassembleMatrix(MPI::Intracomm& communicator, int id);

  // create a new matrix on each MPI node which is an exact clone of the broadcasted one
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the matrix
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new matrix 
  virtual Matrix* BroadcastClone(MPI::Intracomm& communicator, int id);

  // create a new matrix on each MPI node with same size and same type but non-initialized components
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the matrix
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new matrix 
  virtual Matrix* BroadcastEmptyClone(MPI::Intracomm& communicator, int id, bool zeroFlag = false);

#endif

};

// get a matrix element (real part if complex)
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

inline void RealMatrix::GetMatrixElement(int i, int j, double& x) const
{
  x = this->Columns[j].Components[i];
}

// get a matrix element (real part if complex)
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

inline void RealMatrix::GetMatrixElement(int i, int j, Complex& x) const
{
  x = this->Columns[j].Components[i];
}

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
  assert(i<this->NbrColumn);
  return this->Columns[i];
}

#endif
