////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of hermitian matrix                        //
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


#ifndef HERMITIANMATRIX_H
#define HERMITIANMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/ComplexTriDiagonalHermitianMatrix.h"
#include "Output/MathematicaOutput.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class ComplexVector;
class RealMatrix;


class HermitianMatrix : public Matrix
{

  friend class RealVector;
  friend class ComplexVector;
  friend class SingleParticle;

 protected:

  double* DiagonalElements;

  double* RealOffDiagonalElements;
  double* ImaginaryOffDiagonalElements;

  GarbageFlag Flag;

  int Increment;

 public:

  // default constructor
  //
  HermitianMatrix();

  // constructor from matrix elements (without duplicating datas)
  //
  // diagonal = pointer to diagonal element array
  // realOffDiagonal = pointer to real part of off-diagonal elements
  // imaginaryOffDiagonal = pointer to imaginary part of off-diagonal elements
  // dimension = matrix dimension
  HermitianMatrix(double* diagonal, double* realOffDiagonal, double* imaginaryOffDiagonal, int dimension);

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  HermitianMatrix(const HermitianMatrix& M);

  // copy constructor from a real tridiagonal symmetric matrix (without duplicating diagonal elements)
  //
  // M = matrix to copy
  HermitianMatrix(const RealTriDiagonalSymmetricMatrix& M);

  // destructor
  //
  ~HermitianMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  HermitianMatrix& operator = (const HermitianMatrix& M);

  // assignement from  a real tridiagonal symmetric matrix (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix  
  HermitianMatrix& operator = (const RealTriDiagonalSymmetricMatrix& M);

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
  friend HermitianMatrix operator + (const HermitianMatrix& M1, 
						       const HermitianMatrix& M2);

  // add two matrices where the left one is a real tridiagonal symmetric matrix
  //
  // M1 = left matrix
  // M2 = right matrix
  // return value = sum of the two matrices
  friend HermitianMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, 
						       const HermitianMatrix& M2);

  // add two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M1 = left matrix
  // M2 = right matrix
  // return value = sum of the two matrices
  friend HermitianMatrix operator + (const HermitianMatrix& M1, 
						       const RealTriDiagonalSymmetricMatrix& M2);

  // substract two matrices
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend HermitianMatrix operator - (const HermitianMatrix& M1, 
						       const HermitianMatrix& M2);

  // substract two matrices where the left one is a real tridiagonal symmetric matrix
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend HermitianMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, 
						       const HermitianMatrix& M2);

  // substract two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend HermitianMatrix operator - (const HermitianMatrix& M1,
						       const RealTriDiagonalSymmetricMatrix& M2);

  // multiply a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend HermitianMatrix operator * (const HermitianMatrix& M, double x);

  // multiply a matrix by a real number (left multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend HermitianMatrix operator * (double x, const HermitianMatrix& M);

  // divide a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = division result
  friend HermitianMatrix operator / (const HermitianMatrix& M, double x);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  HermitianMatrix& operator += (const HermitianMatrix& M);

  // add two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  HermitianMatrix& operator += (const RealTriDiagonalSymmetricMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  HermitianMatrix& operator -= (const HermitianMatrix& M);

  // substract two matrices where the right one is a real tridiagonal symmetric matrix
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  HermitianMatrix& operator -= (const RealTriDiagonalSymmetricMatrix& M);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  HermitianMatrix& operator *= (double x);

  // divide a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  HermitianMatrix& operator /= (double x) ;

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // conjugate a matrix with an unitary matrix (Ut M U)
  //
  // UnitaryM = unitary matrix to use
  // return value = pointer to conjugated matrix
  Matrix* Conjugate(RealMatrix& UnitaryM);

  // conjugate an hermitian matrix with an unitary matrix (Ut M U)
  //
  // UnitaryM = unitary matrix to use
  // return value = conjugated matrix
  HermitianMatrix Conjugate(ComplexMatrix& UnitaryM);

  // conjugate an hermitian matrix with an hermitian transposed unitary matrix (U M Ut)
  //
  // UnitaryM = unitary matrix to use
  // return value = conjugated matrix  
  HermitianMatrix InvConjugate(ComplexMatrix& UnitaryM);

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
  RealTriDiagonalSymmetricMatrix& Lanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, ComplexVector& V1);

  // Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step
  //
  // dimension = maximum iteration number
  // M = reference on real tridiagonal symmetric matrix where result has to be stored
  // Q = matrix where new orthonormalized base will be stored (first column is used as first vector)
  // return value = reference on complex tridiagonal hermitian matrix
  RealTriDiagonalSymmetricMatrix& Lanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, ComplexMatrix& Q);

  // Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step, if during process a 
  // null vector appears then new random vector is evaluated
  //
  // dimension = maximum iteration number
  // M = reference on real tridiagonal symmetric matrix where result has to be stored
  // Q = matrix where new orthonormalized base will be stored (first column is used as first vector)
  // err = absolute error on vector norm
  // return value = reference on complex tridiagonal hermitian matrix
  RealTriDiagonalSymmetricMatrix& OrthoLanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, ComplexMatrix& Q, 
						double err = 0.00000001);

  // Tridiagonalize a hermitian matrix using Householder algorithm and evaluate transformation matrix
  //
  // M = reference on real tridiagonal symmetric matrix where result has to be stored
  // err = absolute error on matrix element
  // Q = matrix where transformation matrix has to be stored
  // return value = reference on real tridiagonal symmetric matrix
  RealTriDiagonalSymmetricMatrix& Householder (RealTriDiagonalSymmetricMatrix& M, double err, ComplexMatrix& Q);

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const HermitianMatrix& P);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const HermitianMatrix& P);

};

#endif
