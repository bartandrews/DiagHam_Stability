////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of complex upper triangular matrix                 //
//                                                                            //
//                        last modification : 20/08/2004                      //
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


#ifndef COMPLEXUPPERTRIANGULARMATRIX_H
#define COMPLEXUPPERTRIANGULARMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"

#include <iostream>
#include <fstream>


using std::ostream;
using std::ofstream;
using std::ifstream;


class ComplexMatrix;
class ComplexLowerTriangularMatrix;


class ComplexUpperTriangularMatrix : public Matrix
{

  friend class RealVector;
  friend class ComplexVector;
  friend class ComplexMatrix;

 protected:

  // real part of the upper off diagonal elements
  Complex* OffDiagonalElements;
  // garbage flag used for the matrix upper off diagonal elements
  GarbageFlag OffDiagonalFlag;

  // real part of the diagonal elements
  Complex* DiagonalElements;
  // garbage flag used for the matrix diagonal elements
  GarbageFlag DiagonalFlag;

  // dummy variable whose reference is send when an element of the lower part of the matrix is asked (initialize to 0)
  Complex Dummy;

 public:

  // default constructor
  //
  ComplexUpperTriangularMatrix();

  // constructor for an empty matrix
  //
  // dimension = matrix dimension
  // zero = true if matrix has to be filled with zeros
  ComplexUpperTriangularMatrix(int dimension, bool zero = false);

  // constructor from matrix elements (without duplicating datas)
  //
  // diagonal = pointer to the diagonal elements
  // offDiagonal = pointer to the off-diagonal elements
  // dimension = matrix dimension
  ComplexUpperTriangularMatrix(Complex* diagonal, Complex* offDiagonal, int dimension);

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  ComplexUpperTriangularMatrix(const ComplexUpperTriangularMatrix& M);

  // destructor
  //
  ~ComplexUpperTriangularMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  ComplexUpperTriangularMatrix& operator = (const ComplexUpperTriangularMatrix& M);

  // return pointer on a clone matrix (without duplicating datas)
  //
  // retrun value = pointer on new matrix 
  Matrix* Clone ();

  // copy a matrix into another (duplicating data)
  //
  // matrix = matrix to copy
  // return value = reference on current matrix
  ComplexUpperTriangularMatrix& Copy (ComplexUpperTriangularMatrix& matrix);

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
  friend ComplexUpperTriangularMatrix operator + (const ComplexUpperTriangularMatrix& M1, 
						  const ComplexUpperTriangularMatrix& M2);

  // substract two matrices
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend ComplexUpperTriangularMatrix operator - (const ComplexUpperTriangularMatrix& M1, 
						  const ComplexUpperTriangularMatrix& M2);

  // multiply a complex lower triangular matrix with a complex upper triangular matrix
  //
  // m1 = complex lower triangular matrix
  // m2 = complex upper triangular matrix
  // return value = product result
  friend ComplexMatrix operator * (ComplexLowerTriangularMatrix& m1, ComplexUpperTriangularMatrix& m2);

  // multiply a complex upper triangular matrix with a complex lower triangular matrix
  //
  // m1 = complex upper triangular matrix
  // m2 = complex lower triangular matrix
  // return value = product result
  friend ComplexMatrix operator * (ComplexUpperTriangularMatrix& m1, ComplexLowerTriangularMatrix& m2);

  // multiply a complex matrix with a complex upper triangular matrix
  //
  // m1 = complex matrix
  // m2 = complex upper triangular matrix
  // return value = product result
  friend ComplexMatrix operator * (ComplexMatrix& m1, const ComplexUpperTriangularMatrix& m2);

  // multiply a matrix by a complex number (right multiplication)
  //
  // M = source matrix
  // x = complex number to use
  // return value = product result
  friend ComplexUpperTriangularMatrix operator * (const ComplexUpperTriangularMatrix& M, double x);

  // multiply a matrix by a complex number (left multiplication)
  //
  // M = source matrix
  // x = complex number to use
  // return value = product result
  friend ComplexUpperTriangularMatrix operator * (double x, const ComplexUpperTriangularMatrix& M);

  // divide a matrix by a complex number (right multiplication)
  //
  // M = source matrix
  // x = complex number to use
  // return value = division result
  friend ComplexUpperTriangularMatrix operator / (const ComplexUpperTriangularMatrix& M, double x);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  ComplexUpperTriangularMatrix& operator += (const ComplexUpperTriangularMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  ComplexUpperTriangularMatrix& operator -= (const ComplexUpperTriangularMatrix& M);

  // multiply a matrix by a complex number
  //
  // x = complex number to use
  // return value = reference on current matrix
  ComplexUpperTriangularMatrix& operator *= (double x);

  // divide a matrix by a complex number
  //
  // x = complex number to use
  // return value = reference on current matrix
  ComplexUpperTriangularMatrix& operator /= (double x) ;

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // Solve the linear equation M x = y
  //
  // x = vector where the solution will be stored
  // y = vector that gives the right hand side of the equation
  // return value = true if no error occured
  bool SolveLinearEquation (ComplexVector& x, ComplexVector& y);

  // invert the current matrix
  //
  // return value = true if no error occured
  bool Invert ();

  // compute the invert of a matrix from its PLU decomposition
  // 
  // lowerMatrix = reference on the matrix where the lower triangular matrix
  // upperMatrix = reference on the matrix where the upper triangular matrix
  // permutations = array that list all the permutations defining P. Each permutation is given at a pair corresponding to an index i and 
  //                the i-th entry in the array (i.e. i <-> permutations[i]). The sequence is performed from the latest entry of permutations to the first one
  // return value = inverted matrix
  friend ComplexMatrix InvertMatrixFromLUDecomposition(ComplexLowerTriangularMatrix& lowerMatrix, ComplexUpperTriangularMatrix& upperMatrix, int* permutations);

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
  friend ostream& operator << (ostream& Str, const ComplexUpperTriangularMatrix& P);

#ifdef USE_OUTPUT

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexUpperTriangularMatrix& P);

#endif

 protected:

  // get the linearized index to access an off-diagonal matrix element
  //
  // i = row index
  // j = column index
  // return value = linearized index
  long GetLinearizedOffDiagonalIndex(int i, int j) const;


};

// get the linearized index to access an off-diagonal matrix element
//
// i = row index
// j = column index
// return value = linearized index

inline long ComplexUpperTriangularMatrix::GetLinearizedOffDiagonalIndex(int i, int j) const 
{
  return (i + (j * (j - 1l)) / 2l);
}

// get a matrix element (real part if complex)
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

inline void ComplexUpperTriangularMatrix::GetMatrixElement(int i, int j, double& x) const
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
	  x = this->OffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, j)].Re;
	}
      else
	{
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

inline void ComplexUpperTriangularMatrix::GetMatrixElement(int i, int j, Complex& x) const
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      x = this->DiagonalElements[i];
    }
  else
    {
      if (i < j)
	{
	  x = this->OffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, j)];
	}
      else
	{
	  x = 0.0;
	}
    }      
  return;
}

#endif
