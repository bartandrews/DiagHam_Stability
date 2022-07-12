////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of complex skew symmetric matrix                  //
//                                                                            //
//                        last modification : 19/08/2004                      //
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


#ifndef COMPLEXSKEWSYMMETRICMATRIX_H
#define COMPLEXSKEWSYMMETRICMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class ComplexMatrix;
class BlockDiagonalMatrix;
class ComplexVector;
class ComplexVector;



class ComplexSkewSymmetricMatrix : public Matrix
{

  friend class RealVector;
  friend class ComplexVector;

 private:
  
 // dummy variable whose reference is send when an element of the lower part of the matrix is asked (initialize to 0)
  double Dummy;

 protected:

  // real part of the upper off diagonal elements
  double* RealOffDiagonalElements;
  // imaginary part of the upper off diagonal elements
  double* ImaginaryOffDiagonalElements;

  // garbage flag used for the matrix elements
  GarbageFlag Flag;

  // increment to add to the end of each line to go to the next line minus 1
  int Increment;

 public:

  // default constructor
  //
  ComplexSkewSymmetricMatrix();

  // constructor for an empty matrix
  //
  // dimension = matrix dimension
  // zero = true if matrix has to be filled with zeros
  ComplexSkewSymmetricMatrix(int dimension, bool zero = false);

  // constructor from matrix elements (without duplicating datas)
  //
  // realUpperDiagonal = pointer to real part of the upper-diagonal elements
  // imaginaryUpperDiagonal = pointer to imaginary part of the upper-diagonal elements
  // dimension = matrix dimension
  ComplexSkewSymmetricMatrix(double* realUpperDiagonal, double* imaginaryUpperDiagonal, int dimension) ;

  // copy constructor
  //
  // M = matrix to copy
  // duplicateFlag = true if datas have to be duplicated
  ComplexSkewSymmetricMatrix(const ComplexSkewSymmetricMatrix& M, bool duplicateFlag);

  // destructor
  //
  ~ComplexSkewSymmetricMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  ComplexSkewSymmetricMatrix& operator = (const ComplexSkewSymmetricMatrix& M);

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

  // add two matrices
  //
  // M1 = first matrix
  // M2 = second matrix
  // return value = sum of the two matrices
  friend ComplexSkewSymmetricMatrix operator + (const ComplexSkewSymmetricMatrix& M1, 
						const ComplexSkewSymmetricMatrix& M2);

  // substract two matrices
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend ComplexSkewSymmetricMatrix operator - (const ComplexSkewSymmetricMatrix& M1, 
						const ComplexSkewSymmetricMatrix& M2);

  // multiply a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend ComplexSkewSymmetricMatrix operator * (const ComplexSkewSymmetricMatrix& M, double x);

  // multiply a matrix by a real number (left multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend ComplexSkewSymmetricMatrix operator * (double x, const ComplexSkewSymmetricMatrix& M);

  // divide a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = division result
  friend ComplexSkewSymmetricMatrix operator / (const ComplexSkewSymmetricMatrix& M, double x);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  ComplexSkewSymmetricMatrix& operator += (const ComplexSkewSymmetricMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  ComplexSkewSymmetricMatrix& operator -= (const ComplexSkewSymmetricMatrix& M);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  ComplexSkewSymmetricMatrix& operator *= (double x);

  // divide a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  ComplexSkewSymmetricMatrix& operator /= (double x) ;

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // conjugate a matrix with an unitary complex matrix (Ut M U)
  //
  // UnitaryM = unitary matrix to use
  // return value = pointer to conjugated matrix
  Matrix* Conjugate(ComplexMatrix& UnitaryM);

  // conjugate a block of the matrix with an unitary matrix (Ut M U)
  //
  // UnitaryM = unitary matrix to use
  // sourcePosition = index of the row where the block to conjugate starts
  // destinationPosition = index of the row where the conjugated block has to be stored
  // matrix = matrix where result has to be stored
  void Conjugate(ComplexMatrix& UnitaryM, int sourcePosition, int destinationPosition,
		 ComplexSkewSymmetricMatrix& matrix);

  // conjugate a block of the matrix (in the upper diagonal part) with two matrix matrix (Vt M U)
  //
  // UnitaryMl = unitary matrix to use at the left hand side
  // UnitaryMr = unitary matrix to use at the right hand side
  // sourceRowIndex = index of the row where the block to conjugate starts
  // sourceColumnIndex = index of the column where the block to conjugate starts
  // destinationRowIndex = index of the row where the conjugated block has to be stored
  // destinationColumnIndex = index of the column where the conjugated block has to be stored
  // matrix = matrix where result has to be stored
  void Conjugate(ComplexMatrix& UnitaryMl, ComplexMatrix& UnitaryMr, int sourceRowIndex, 
		 int sourceColumnIndex, int destinationRowIndex,
		 int destinationColumnIndex, ComplexSkewSymmetricMatrix& matrix);

  // swap the i-th row/column with the j-th row/column (thus preserving the skew symmetric form)
  //
  // i = index of the first the row/column
  // j = index of the second the row/column
  // return value = reference on the current matrix
  ComplexSkewSymmetricMatrix& SwapRowColumn (int i, int j);

  // evaluate matrix trace
  //
  // return value = matrix trace 
  double Tr ();

  // evaluate matrix determinant
  //
  // return value = matrix determinant 
  double Det ();

  // evaluate matrix pfaffian
  //
  // return value = matrix pfaffian 
  Complex Pfaffian();

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const ComplexSkewSymmetricMatrix& P);

#ifdef USE_OUTPUT

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexSkewSymmetricMatrix& P);

#endif

};

#endif
