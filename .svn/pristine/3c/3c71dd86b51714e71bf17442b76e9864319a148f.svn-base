////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          base class of for matrix                          //
//                                                                            //
//                        last modification : 05/01/2001                      //
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


#ifndef MATRIX_H
#define MATRIX_H


#include "config.h"

#include <iostream>


using std::ostream;


class Complex;
class RealMatrix;
class BlockDiagonalMatrix;
class SubspaceSpaceConverter;


class Matrix
{

  friend class RealVector;

 private:
  
  double Dummy;

 protected:

  int NbrRow;
  int NbrColumn;

  int TrueNbrRow;
  int TrueNbrColumn;

  int MatrixType;

 public:

  enum MatrixProperties
  {
    RealElements = 0x0001,
    ComplexElements = 0x0002,
    Diagonal = 0x0010,
    TriDiagonal = 0x0020,
    Triangular = 0x0040,
    BlockDiagonal = 0x0080,
    Symmetric = 0x0100,
    Antisymmetric = 0x0200,
    Hermitian = 0x0400,
    AntiHermitian = 0x0800,    
    Block = 0x1000
  };

  // default constructor
  //
  Matrix();

  // virtual destructor
  //
  virtual ~Matrix();

  // return pointer on a clone matrix (without duplicating datas)
  //
  // retrun value = pointer on new matrix 
  virtual Matrix* Clone ();

  // project matrix into a given subspace
  //
  // subspace = reference on subspace structure
  // return value = pointer to projected matrix
  virtual Matrix* Project (SubspaceSpaceConverter& subspace);  

  // get matrix type
  //
  // return value = matrix type 
  virtual int GetMatrixType ();

  // get number of row
  //
  // return value = number of row
  virtual int GetNbrRow ();

  // get number of column
  //
  // return value = number of column
  virtual int GetNbrColumn ();

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  virtual void SetMatrixElement(int i, int j, double x);

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  virtual void SetMatrixElement(int i, int j, const Complex& x);

  // add a value to a matrix element
  //
  // i = line position
  // j = column position
  // x = value to add to matrix element
  virtual void AddToMatrixElement(int i, int j, double x);

  // add a value  a matrix element
  //
  // i = line position
  // j = column position
  // x = value to add to matrix element
  virtual void AddToMatrixElement(int i, int j, const Complex& x);

  // Resize matrix
  //
  // nbrRow = new number of rows
  // nbrColumn = new number of columns
  virtual void Resize (int nbrRow, int nbrColumn);

  // Resize matrix and set to zero all elements that have been added
  //
  // nbrRow = new number of rows
  // nbrColumn = new number of columns
  virtual void ResizeAndClean (int nbrRow, int nbrColumn);

  // return refernce on real part of a given matrix element
  //
  // i = line position
  // j = column position
  // return value = reference on real part
  virtual double& operator () (int i, int j);

  // conjugate matrix with an unitary real matrix (Ut M U)
  //
  // UnitaryM = unitary matrix to use
  // return value = pointer to conjugated matrix
  Matrix* Conjugate (RealMatrix& UnitaryM);

  // conjugate matrix with an unitary block diagonal matrix (Ut M U)
  //
  // UnitaryM = unitary matrix to use
  // return value = pointer to conjugated matrix
  virtual Matrix* Conjugate (BlockDiagonalMatrix& UnitaryM);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  //  virtual Matrix& operator *= (double x);

  // evaluate matrix trace
  //
  // return value = matrix trace 
  virtual double Tr ();

  // evaluate matrix determinant
  //
  // return value = matrix determinant 
  virtual double Det ();

  // Output Stream overload
  //
  // str = reference on output stream
  // matrix = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, const Matrix& matrix);

};

// get matrix type
//
// return value = matrix type 

inline int Matrix::GetMatrixType ()
{
  return this->MatrixType;
}

// get number of row
//
// return value = number of row

inline int Matrix::GetNbrRow ()
{
  return this->NbrRow;
}

// get number of column
//
// return value = number of column

inline int Matrix::GetNbrColumn ()
{
  return this->NbrColumn;
}

#endif
