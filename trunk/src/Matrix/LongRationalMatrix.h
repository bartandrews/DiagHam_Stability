////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of long rational matrix                      //
//                                                                            //
//                        last modification : 24/10/2012                      //
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


#ifndef LONGRATIONALMATRIX_H
#define LONGRATIONALMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif
#include "MathTools/LongRational.h"
#include "Vector/LongRationalVector.h"

#include <iostream>


using std::ostream;




class LongRationalMatrix : public Matrix
{

  friend class LongRationalVector;

 protected:

  
  LongRationalVector* Columns; 
  int* ColumnGarbageFlag;

  // dummy temporary double
  double Dummy;

 public:

  // default constructor
  //
  LongRationalMatrix();

  // constructor for an empty matrix
  //
  // nbrRow = number of rows
  // nbrColumn = number of columns
  // zero = tue if matrix elements have to be set to zero
  LongRationalMatrix(int nbrRow, int nbrColumn, bool zero = false);

  // constructor from matrix elements (without duplicating datas)
  //
  // columns = pointer an array of vector
  // nbrColumn = number of columns
  LongRationalMatrix(LongRationalVector* columns, int nbrColumn);

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  LongRationalMatrix(const LongRationalMatrix& M);

  // destructor
  //
  ~LongRationalMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  LongRationalMatrix& operator = (const LongRationalMatrix& M);

  // return pointer on a clone matrix (without duplicating datas)
  //
  // retrun value = pointer on new matrix 
  Matrix* Clone ();

  // copy a matrix into another (duplicating data)
  //
  // matrix = matrix to copy
  // return value = reference on current matrix
  LongRationalMatrix& Copy (LongRationalMatrix& matrix);

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
  void GetMatrixElement(int i, int j, LongRational& x) const;

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  void SetMatrixElement(int i, int j, const LongRational& x);

  // add a value  a matrix element
  //
  // i = line position
  // j = column position
  // x = value to add to matrix element
  void AddToMatrixElement(int i, int j, const LongRational& x);

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
  LongRationalVector& operator [] (int i);

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
  friend LongRationalMatrix operator + (const LongRationalMatrix& M1, const LongRationalMatrix& M2);

  // substract two matrices
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend LongRationalMatrix operator - (const LongRationalMatrix& M1, const LongRationalMatrix& M2);

  // multiply two matrices
  //
  // M1 = first matrix
  // M2 = matrix to multiply to M1
  // return value = product of the two matrices
  friend LongRationalMatrix operator * (const LongRationalMatrix& M1, const LongRationalMatrix& M2);

  // multiply a matrix by a long rational number (right multiplication)
  //
  // M = source matrix
  // x = long rational number to use
  // return value = product result
  friend LongRationalMatrix operator * (const LongRationalMatrix& M, const LongRational& x);

  // multiply a matrix by a long rational number (left multiplication)
  //
  // M = source matrix
  // x = long rational number to use
  // return value = product result
  friend LongRationalMatrix operator * (const LongRational& x, const LongRationalMatrix& M);

  // multiply a matrix to the right by another matrix without using temporary matrix
  //
  // M = matrix used as multiplicator
  // return value = reference on current matrix
  LongRationalMatrix& Multiply (const LongRationalMatrix& M);

  // multiply a matrix to the right by another matrix without using temporary matrix and in a given range of indices
  // beware the matrix is not resized after multiplication in order the operation to be thread safe
  //
  // M = matrix used as multiplicator
  // startLine = starting line in destination matrix
  // nbrLine = number of lines to multiply
  // return value = reference on current matrix
  LongRationalMatrix& Multiply (const LongRationalMatrix& M, int startLine, int nbrLine);

  // divide a matrix by a long rational number (right multiplication)
  //
  // M = source matrix
  // x = long rational number to use
  // return value = division result
  friend LongRationalMatrix operator / (const LongRationalMatrix& M, const LongRationalMatrix x);

  // divide a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = division result
  friend LongRationalMatrix operator / (const LongRationalMatrix& M, const LongRational& x);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  LongRationalMatrix& operator += (const LongRationalMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  LongRationalMatrix& operator -= (const LongRationalMatrix& M);

  // multiply a matrix by a long rational number
  //
  // x = long rational number to use
  // return value = reference on current matrix
  LongRationalMatrix& operator *= (const LongRational& x);

  // divide a matrix by a long rational number
  //
  // x = long rational number to use
  // return value = reference on current matrix
  LongRationalMatrix& operator /= (const LongRational& x) ;

  // transpose matrix
  //
  // return value = reference on current matrix
  LongRationalMatrix& Transpose ();

  // duplicate and transpose a matrix
  //
  // return value = transposed matrix
  LongRationalMatrix DuplicateAndTranspose();
    
  // evaluate matrix determinant (skrewing up matrix elements)
  //
  // return value = matrix determinant 
  LongRational Determinant ();

  // evaluate permanent associated to the (square) matrix using Ryser algorithm
  //
  // return value = permanent associated to the matrix
  LongRational Permanent();

  // evaluate matrix rank
  //
  // accuracy = numerical accuracy used to define linearly dependence 
  // return value = rank
  virtual int Rank(double accuracy = MACHINE_PRECISION);

  // write matrix in a file 
  //
  // file = reference on the output file stream
  // return value = true if no error occurs
  virtual bool WriteMatrix (ofstream& file);

  // write matrix in a file 
  //
  // fileName = name of the file where the matrix has to be stored
  // return value = true if no error occurs
  virtual bool WriteMatrix (char* fileName);

  // read matrix from a file 
  //
  // file = reference  on the input file stream
  // return value = true if no error occurs
  virtual bool ReadMatrix (ifstream& file);

  // read matrix from a file 
  //
  // fileName = name of the file where the matrix has to be read
  // return value = true if no error occurs
  virtual bool ReadMatrix (char* fileName);

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const LongRationalMatrix& P);

};

// get a matrix element (long rational part if complex)
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

inline void LongRationalMatrix::GetMatrixElement(int i, int j, double& x) const
{
  x = this->Columns[j][i].GetNumericalValue();
}

// get a matrix element
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

inline void LongRationalMatrix::GetMatrixElement(int i, int j, LongRational& x) const
{ 
  x = this->Columns[j].Components[i];
}

// get reference of a given matrix element
//
// i = line position
// j = column position
// return value = reference on matrix elememt

inline double& LongRationalMatrix::operator () (int i, int j)
{
  this->Dummy = this->Columns[j].Components[i].GetNumericalValue();
  return  this->Dummy;
}

// get reference to a given column
//
// i = column position
// return value = column reference 

inline LongRationalVector& LongRationalMatrix::operator [] (int i)
{
  return this->Columns[i];
}

#endif
