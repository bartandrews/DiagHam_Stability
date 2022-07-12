////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of complex upper Hessenberg matrix                 //
//                                                                            //
//                        last modification : 27/11/2012                      //
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


#ifndef COMPLEXUPPERHESSENBERGMATRIX_H
#define COMPLEXUPPERHESSENBERGMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif
#include "GeneralTools/GarbageFlag.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Matrix/ComplexUpperTriangularMatrix.h"
#include "MathTools/Complex.h"

#include <iostream>
#include <fstream>


using std::ostream;
using std::ofstream;
using std::ifstream;


class ComplexMatrix;


class ComplexUpperHessenbergMatrix : public Matrix
{

  friend class RealVector;
  friend class ComplexVector;
  friend class ComplexMatrix;

 protected:

  // upper off-diagonal elements
  Complex* UpperOffDiagonalElements;
  // diagonal elements
  Complex* DiagonalElements;
  // lower diagonal elements (the first element is unused, row index has to be used as the array index)
  Complex* LowerDiagonalElements;

  // garbage flag 
  GarbageFlag Flag;

  // dummy variable whose reference is send when an element of the lower part of the matrix is asked (initialize to 0)
  double Dummy;

 public:

  // default constructor
  //
  ComplexUpperHessenbergMatrix();

  // constructor for an empty matrix
  //
  // dimension = matrix dimension
  // zero = true if matrix has to be filled with zeros
  ComplexUpperHessenbergMatrix(int dimension, bool zero = false);

  // constructor from matrix elements (without duplicating datas)
  //
  // diagonalElements = diagonal elements
  // offDiagonalElements = upper off-diagonal elements
  // lowerDiagonalElements = lower diagonal elements
  // dimension = matrix dimension
  ComplexUpperHessenbergMatrix(Complex* diagonalElements, Complex* offDiagonalElements, 
			       Complex* lowerDiagonalElements, int dimension);

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  ComplexUpperHessenbergMatrix(const ComplexUpperHessenbergMatrix& M);

  // constructor from a matrix, copying the data and discarding all elements below the lower diagonal
  //
  // M = reference on the matrix
  ComplexUpperHessenbergMatrix(Matrix& M);

  // destructor
  //
  ~ComplexUpperHessenbergMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  ComplexUpperHessenbergMatrix& operator = (const ComplexUpperHessenbergMatrix& M);

  // return pointer on a clone matrix (without duplicating datas)
  //
  // retrun value = pointer on new matrix 
  Matrix* Clone ();

  // copy a matrix into another (duplicating data)
  //
  // matrix = matrix to copy
  // return value = reference on current matrix
  ComplexUpperHessenbergMatrix& Copy (ComplexUpperHessenbergMatrix& matrix);

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
  friend ComplexUpperHessenbergMatrix operator + (const ComplexUpperHessenbergMatrix& M1, 
						  const ComplexUpperHessenbergMatrix& M2);

  // substract two matrices
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend ComplexUpperHessenbergMatrix operator - (const ComplexUpperHessenbergMatrix& M1, 
						  const ComplexUpperHessenbergMatrix& M2);

  // multiply a matrix by a complex number (right multiplication)
  //
  // M = source matrix
  // x = complex number to use
  // return value = product result
  friend ComplexUpperHessenbergMatrix operator * (const ComplexUpperHessenbergMatrix& M, double x);

  // multiply a matrix by a complex number (left multiplication)
  //
  // M = source matrix
  // x = complex number to use
  // return value = product result
  friend ComplexUpperHessenbergMatrix operator * (double x, const ComplexUpperHessenbergMatrix& M);

  // divide a matrix by a complex number (right multiplication)
  //
  // M = source matrix
  // x = complex number to use
  // return value = division result
  friend ComplexUpperHessenbergMatrix operator / (const ComplexUpperHessenbergMatrix& M, double x);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  ComplexUpperHessenbergMatrix& operator += (const ComplexUpperHessenbergMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  ComplexUpperHessenbergMatrix& operator -= (const ComplexUpperHessenbergMatrix& M);

  // multiply a matrix by a complex number
  //
  // x = complex number to use
  // return value = reference on current matrix
  ComplexUpperHessenbergMatrix& operator *= (double x);

  // divide a matrix by a complex number
  //
  // x = complex number to use
  // return value = reference on current matrix
  ComplexUpperHessenbergMatrix& operator /= (double x) ;

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // evaluate matrix trace
  //
  // return value = matrix trace 
  double Tr ();

  // evaluate matrix determinant
  //
  // return value = matrix determinant 
  double Det ();

  // shift all diagonal elements 
  //
  // shift = shift to apply
  // return value = reference on current matrix
  ComplexUpperHessenbergMatrix& ShiftDiagonal(const Complex& shift);

  // conjugate matrix with an unitary matrix (Uh M U), assuming the Hessenberg from will be preserved
  //
  // unitaryM = unitary matrix to use
  // conjugatedMatrix = reference on the matrix where conjugate matrix will be stored
  // return value = pointer to conjugated matrix
  ComplexUpperHessenbergMatrix& Conjugate(ComplexMatrix& unitaryM, ComplexUpperHessenbergMatrix& conjugatedMatrix);

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

  // find QR factorization using the LAPACK library
  //
  // R = reference on the triangular matrix
  // Q = reference on the transformation matrix
  // return value = reference on upper triangular matrix
  ComplexUpperTriangularMatrix& LapackQRFactorization (ComplexUpperTriangularMatrix& R, ComplexMatrix& Q);

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const ComplexUpperHessenbergMatrix& P);

#ifdef USE_OUTPUT

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexUpperHessenbergMatrix& P);

#endif

};

// get a matrix element (real part if complex)
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

inline void ComplexUpperHessenbergMatrix::GetMatrixElement(int i, int j, double& x) const
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      x= this->DiagonalElements[i].Re;
    }
  else
    {
      if (i < j)
	{
	  x = this->UpperOffDiagonalElements[i + (j * (j - 1l)) / 2l].Re;
	}
      else
	{
	  if ((j + 1) == i)
	    {
	      x = this->LowerDiagonalElements[i].Re;
	    }
	  else
	    x = 0.0;
	}
    }      
  return;
}

// get a matrix element (real part if complex)
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

inline void ComplexUpperHessenbergMatrix::GetMatrixElement(int i, int j, Complex& x) const
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      x= this->DiagonalElements[i];
    }
  else
    {
      if (i < j)
	{
	  x = this->UpperOffDiagonalElements[i + (j * (j - 1l)) / 2l];
	}
      else
	{
	  if ((j + 1) == i)
	    {
	      x = this->LowerDiagonalElements[i];
	    }
	  else
	    x = 0.0;
	}
    }      
  return;
}

#endif
