////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of real antisymmetric matrix                    //
//                                                                            //
//                        last modification : 03/04/2001                      //
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


#ifndef REALANTISYMMETRICMATRIX_H
#define REALANTISYMMETRICMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#include "Output/MathematicaOutput.h"

#include <iostream>


using std::ostream;


class RealMatrix;
class BlockDiagonalMatrix;


class RealAntisymmetricMatrix : public Matrix
{

  friend class RealVector;
  friend class ComplexVector;

 private:
  
  double Dummy;

 protected:

  double* OffDiagonalElements;
  int* OffDiagonalGarbageFlag;

  int Increment;

 public:

  // default constructor
  //
  RealAntisymmetricMatrix();

  // constructor for an empty matrix
  //
  // dimension = matrix dimension
  // zero = true if matrix has to be filled with zeros
  RealAntisymmetricMatrix(int dimension, bool zero = false);

  // constructor from matrix elements (without duplicating datas)
  //
  // upperDiagonal = pointer to upper-diagonal element array (with real part in even position and imaginary part in odd position)
  // dimension = matrix dimension
  RealAntisymmetricMatrix(double* upperDiagonal, int dimension) ;

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  RealAntisymmetricMatrix(const RealAntisymmetricMatrix& M);

  // destructor
  //
  ~RealAntisymmetricMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  RealAntisymmetricMatrix& operator = (const RealAntisymmetricMatrix& M);

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

  // get reference of a given matrix element supposing i < j
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

  // project matrix into a given subspace
  //
  // subspace = reference on subspace structure
  // return value = pointer to projected matrix
  Matrix* Project (SubspaceSpaceConverter& subspace);  

  // add two matrices
  //
  // M1 = first matrix
  // M2 = second matrix
  // return value = sum of the two matrices
  friend RealAntisymmetricMatrix operator + (const RealAntisymmetricMatrix& M1, 
					     const RealAntisymmetricMatrix& M2);

  // substract two matrices
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend RealAntisymmetricMatrix operator - (const RealAntisymmetricMatrix& M1, 
					     const RealAntisymmetricMatrix& M2);

  // multiply a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend RealAntisymmetricMatrix operator * (const RealAntisymmetricMatrix& M, double x);

  // multiply a matrix by a real number (left multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend RealAntisymmetricMatrix operator * (double x, const RealAntisymmetricMatrix& M);

  // divide a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = division result
  friend RealAntisymmetricMatrix operator / (const RealAntisymmetricMatrix& M, double x);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  RealAntisymmetricMatrix& operator += (const RealAntisymmetricMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  RealAntisymmetricMatrix& operator -= (const RealAntisymmetricMatrix& M);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  RealAntisymmetricMatrix& operator *= (double x);

  // divide a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  RealAntisymmetricMatrix& operator /= (double x) ;

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  double MatrixElement (RealVector& V1, RealVector& V2);

  // conjugate a matrix with an unitary real matrix (Ut M U)
  //
  // UnitaryM = unitary matrix to use
  // return value = pointer to conjugated matrix
  Matrix* Conjugate(RealMatrix& UnitaryM);

  // conjugate a matrix with an unitary block-diagonal matrix (Ut M U)
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
		 RealAntisymmetricMatrix& matrix);

  // conjugate a block of the matrix (in the upper diagonal part) with two matrix matrix (Vt M U)
  //
  // UnitaryMl = unitary matrix to use at the left hand side
  // UnitaryMr = unitary matrix to use at the right hand side
  // sourceRowIndex = index of the row where the block to conjugate starts
  // sourceColumnIndex = index of the column where the block to conjugate starts
  // destinationRowIndex = index of the row where the conjugated block has to be stored
  // destinationColumnIndex = index of the column where the conjugated block has to be stored
  // matrix = matrix where result has to be stored
  void Conjugate(RealMatrix& UnitaryMl, RealMatrix& UnitaryMr, int sourceRowIndex, 
		 int sourceColumnIndex, int destinationRowIndex,
		 int destinationColumnIndex, RealAntisymmetricMatrix& matrix);

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
  friend ostream& operator << (ostream& Str, const RealAntisymmetricMatrix& P);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const RealAntisymmetricMatrix& P);

};

#endif
