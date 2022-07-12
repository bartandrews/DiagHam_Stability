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


#include "Matrix/ComplexUpperTriangularMatrix.h"
#include "Matrix/ComplexLowerTriangularMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "GeneralTools/ListIterator.h"
#include "MathTools/Complex.h"

#include <stdlib.h>
#include <fstream>


using std::cout;
using std::endl;


// default constructor
//

ComplexUpperTriangularMatrix::ComplexUpperTriangularMatrix() 
{
  this->DiagonalElements = 0;
  this->OffDiagonalElements = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Triangular | Matrix::Upper;
  this->Dummy = 0.0;
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

ComplexUpperTriangularMatrix::ComplexUpperTriangularMatrix(int dimension, bool zero) 
{
  this->DiagonalFlag.Initialize();
  this->OffDiagonalFlag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Triangular | Matrix::Upper;
  this->DiagonalElements = new Complex [this->NbrRow];
  long TmpNbrOffDiagonalElements = (((long) this->NbrRow) * (((long) this->NbrRow) - 1l)) / 2l;
  this->OffDiagonalElements = new Complex [TmpNbrOffDiagonalElements];
  if (zero == true)
    {
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  this->DiagonalElements[i] = 0.0;
	}
      for (long i = 0l; i < TmpNbrOffDiagonalElements; ++i)
	{
	  this->OffDiagonalElements[i] = 0.0;
	}
    }
  this->Dummy = 0.0;
}

// constructor from matrix elements (without duplicating datas)
//
// diagonal = pointer to the diagonal elements
// offDiagonal = pointer to the off-diagonal elements
// dimension = matrix dimension

ComplexUpperTriangularMatrix::ComplexUpperTriangularMatrix(Complex* diagonal, Complex* offDiagonal, int dimension) 
{
  this->DiagonalElements = diagonal;
  this->OffDiagonalElements = offDiagonal;
  this->DiagonalFlag.Initialize();
  this->OffDiagonalFlag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Triangular | Matrix::Upper;
  this->Dummy = 0.0;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

ComplexUpperTriangularMatrix::ComplexUpperTriangularMatrix(const ComplexUpperTriangularMatrix& M) 
{
  this->DiagonalElements = M.DiagonalElements;
  this->DiagonalFlag = M.DiagonalFlag;
  this->OffDiagonalElements = M.OffDiagonalElements;
  this->OffDiagonalFlag = M.OffDiagonalFlag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Triangular | Matrix::Upper;
  this->Dummy = 0.0;
}

// destructor
//

ComplexUpperTriangularMatrix::~ComplexUpperTriangularMatrix() 
{
  if ((this->DiagonalElements != 0) && 
      (this->DiagonalFlag.Used() == true) && (this->DiagonalFlag.Shared() == false))
    {
      delete[] this->DiagonalElements;
    }
  if ((this->OffDiagonalElements != 0) && 
      (this->OffDiagonalFlag.Used() == true) && (this->OffDiagonalFlag.Shared() == false))
    {
      delete[] this->OffDiagonalElements;
    }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

ComplexUpperTriangularMatrix& ComplexUpperTriangularMatrix::operator = (const ComplexUpperTriangularMatrix& M) 
{
  if ((this->DiagonalElements != 0) && 
      (this->DiagonalFlag.Used() == true) && (this->DiagonalFlag.Shared() == false))
    {
      delete[] this->DiagonalElements;
    }
  if ((this->OffDiagonalElements != 0) && 
      (this->OffDiagonalFlag.Used() == true) && (this->OffDiagonalFlag.Shared() == false))
    {
      delete[] this->OffDiagonalElements;
    }
  this->DiagonalElements = M.DiagonalElements;
  this->DiagonalFlag = M.DiagonalFlag;
  this->OffDiagonalElements = M.OffDiagonalElements;
  this->OffDiagonalFlag = M.OffDiagonalFlag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Triangular | Matrix::Upper;
  this->Dummy = 0.0;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* ComplexUpperTriangularMatrix::Clone ()
{
  return ((Matrix*) new ComplexUpperTriangularMatrix (*this));
}

// copy a matrix into another (duplicating data)
//
// matrix = matrix to copy
// return value = reference on current matrix

ComplexUpperTriangularMatrix& ComplexUpperTriangularMatrix::Copy (ComplexUpperTriangularMatrix& matrix)
{
  this->Resize(matrix.NbrRow, matrix.NbrColumn);
  this->DiagonalElements = new Complex [this->NbrRow];
  for (int i = 0; i < this->NbrColumn; ++i)
    this->DiagonalElements[i] = matrix.DiagonalElements[i];
  long TmpNbrOffDiagonalElements = (((long) this->NbrRow) * (((long) this->NbrRow) - 1l)) / 2l;
  for (long i = 0l; i < TmpNbrOffDiagonalElements; ++i)
    this->OffDiagonalElements[i] = matrix.OffDiagonalElements[i];
  return *this;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexUpperTriangularMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn)  || (i > j))
    return;
  if (i == j)
    {
      this->DiagonalElements[i] = x;
    }
  else
    {
      this->OffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, j)] = x;
    }
  return;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexUpperTriangularMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn) || (i > j))
    return;
  if (i == j)
    {
      this->DiagonalElements[i] = x;
    }
  else
    {
      this->OffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, j)] = x;
    }
  return;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void ComplexUpperTriangularMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn) || (i > j))
    return;
  if (i == j)
    {
      this->DiagonalElements[i] += x;
    }
  else
    {
      this->OffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, j)] += x;
    }
  return;
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void ComplexUpperTriangularMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn) || (i > j))
    return;
  if (i == j)
    {
      this->DiagonalElements[i] += x;
    }
  else
    {
      this->OffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, j)] += x;
    }
  return;
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexUpperTriangularMatrix::Resize (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  Complex* TmpDiagonalElements = new Complex [nbrRow];
  long TmpNbrOffDiagonalElements = (((long) nbrRow) * (((long) nbrRow) - 1l)) / 2l;
  Complex* TmpOffDiagonalElements = new Complex [TmpNbrOffDiagonalElements];
  int MinNbrRow = nbrRow;
  if (nbrRow > this->NbrRow)
    MinNbrRow = this->NbrRow;
  for (int i = 0; i < MinNbrRow; ++i)
    {
      TmpDiagonalElements[i] = this->DiagonalElements[i];
      for (int j = i + 1; j < MinNbrRow; ++j)
	{
	  TmpOffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, j)] = this->OffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, j)];
	}
    }
  if ((this->DiagonalElements != 0) && 
      (this->DiagonalFlag.Used() == true) && (this->DiagonalFlag.Shared() == false))
    {
      delete[] this->DiagonalElements;
    }
  if ((this->OffDiagonalElements != 0) && 
      (this->OffDiagonalFlag.Used() == true) && (this->OffDiagonalFlag.Shared() == false))
    {
      delete[] this->OffDiagonalElements;
    }
  if ((this->DiagonalElements != 0) && 
      (this->DiagonalFlag.Used() == true) && (this->DiagonalFlag.Shared() == false))
    {
      delete[] this->DiagonalElements;
    }
  if ((this->OffDiagonalElements != 0) && 
      (this->OffDiagonalFlag.Used() == true) && (this->OffDiagonalFlag.Shared() == false))
    {
      delete[] this->OffDiagonalElements;
    }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->DiagonalElements = TmpDiagonalElements;
  this->OffDiagonalElements = TmpOffDiagonalElements;
  this->DiagonalFlag.Initialize();
  this->OffDiagonalFlag.Initialize();
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexUpperTriangularMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  Complex* TmpDiagonalElements = new Complex [nbrRow];
  long TmpNbrOffDiagonalElements = (((long) nbrRow) * (((long) nbrRow) - 1l)) / 2l;
  Complex* TmpOffDiagonalElements = new Complex [TmpNbrOffDiagonalElements];
  int MinNbrRow = nbrRow;
  if (nbrRow > this->NbrRow)
    MinNbrRow = this->NbrRow;
  for (int i = 0; i < MinNbrRow; ++i)
    {
      TmpDiagonalElements[i] = this->DiagonalElements[i];
      for (int j = i + 1; j < MinNbrRow; ++j)
	{
	  TmpOffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, j)] = this->OffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, j)];
	}
      for (int j = MinNbrRow; j < nbrRow; ++j)
	{
	  TmpOffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, j)] = 0.0;
	}
    }
  for (int i = this->NbrRow; i < nbrRow; ++i)
    {
      TmpDiagonalElements[i] = 0.0;
      for (int j = i + 1; j < nbrRow; ++j)
	{
	  TmpOffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, j)] = 0.0;
	}
    }
  if ((this->DiagonalElements != 0) && 
      (this->DiagonalFlag.Used() == true) && (this->DiagonalFlag.Shared() == false))
    {
      delete[] this->DiagonalElements;
    }
  if ((this->OffDiagonalElements != 0) && 
      (this->OffDiagonalFlag.Used() == true) && (this->OffDiagonalFlag.Shared() == false))
    {
      delete[] this->OffDiagonalElements;
    }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->DiagonalElements = TmpDiagonalElements;
  this->OffDiagonalElements = TmpOffDiagonalElements;
  this->DiagonalFlag.Initialize();
  this->OffDiagonalFlag.Initialize();
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

ComplexUpperTriangularMatrix operator + (const ComplexUpperTriangularMatrix& M1, const ComplexUpperTriangularMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexUpperTriangularMatrix();
  Complex* Diagonal = new Complex [M1.NbrRow];
  long TmpNbrOffDiagonalElements = (((long) M1.NbrRow) * (((long) M1.NbrRow) - 1l)) / 2l;
  Complex* OffDiagonal = new Complex [TmpNbrOffDiagonalElements];
  for (int i = 0; i < M1.NbrRow; i++)
    {
      Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
    }
  for (long i = 0l; i < TmpNbrOffDiagonalElements; ++i)
    {
      OffDiagonal[i] = M1.OffDiagonalElements[i] + M2.OffDiagonalElements[i];      
    }
  return ComplexUpperTriangularMatrix(Diagonal, OffDiagonal, M1.NbrRow);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexUpperTriangularMatrix operator - (const ComplexUpperTriangularMatrix& M1, const ComplexUpperTriangularMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexUpperTriangularMatrix();
  Complex* Diagonal = new Complex [M1.NbrRow];
  long TmpNbrOffDiagonalElements = (((long) M1.NbrRow) * (((long) M1.NbrRow) - 1l)) / 2l;
  Complex* OffDiagonal = new Complex [TmpNbrOffDiagonalElements];
  for (int i = 0; i < M1.NbrRow; i++)
    {
      Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
    }
  for (long i = 0l; i < TmpNbrOffDiagonalElements; ++i)
    {
      OffDiagonal[i] = M1.OffDiagonalElements[i] - M2.OffDiagonalElements[i];      
    }
  return ComplexUpperTriangularMatrix(Diagonal, OffDiagonal, M1.NbrRow);
}

// multiply a real matrix with a real upper triangular matrix
//
// m1 = real matrix
// m2 = real upper triangular matrix
// return value = product result

ComplexMatrix operator * (ComplexMatrix& m1, const ComplexUpperTriangularMatrix& m2)
{
  ComplexMatrix TmpM(m1.GetNbrRow(), m2.NbrRow);
  long Pos = 0l;
  for (int i = 0; i < m1.GetNbrRow(); ++i)
    {
      for (int j = 0; j < m2.NbrRow; ++j)
	{
	  Complex& Tmp = TmpM[j][i];
	  Tmp = m1[j][i] * m2.DiagonalElements[j];
	  Pos = ((j - 1l) * j) / 2l;
	  for (int k = 0; k < j; ++k)
	    {
	      Tmp += m1[k][i] * m2.OffDiagonalElements[Pos];
	      ++Pos;
	    }
	}
    }
  return TmpM;
}

// multiply a complex lower triangular matrix with a complex upper triangular matrix
//
// m1 = complex lower triangular matrix
// m2 = complex upper triangular matrix
// return value = product result

ComplexMatrix operator * (ComplexLowerTriangularMatrix& m1, ComplexUpperTriangularMatrix& m2)
{
  if (m1.NbrColumn != m2.NbrRow)
    {
      return ComplexMatrix();
    }
  ComplexMatrix TmpM(m1.NbrRow, m2.NbrColumn);
  for (int i = 0; i < m1.NbrRow; ++i)
     {
       for (int j = 0; j < i; ++j)
 	{
 	  Complex& Tmp = TmpM[j][i];
	  Tmp = m1.OffDiagonalElements[m1.GetLinearizedOffDiagonalIndex(i, j)] * m2.DiagonalElements[j];
	  for (int k = 0; k < j; ++k)
	    {
	      Tmp += m1.OffDiagonalElements[m1.GetLinearizedOffDiagonalIndex(i, k)] * m2.OffDiagonalElements[m2.GetLinearizedOffDiagonalIndex(k, j)];
	    }
	}
       {
	 Complex& Tmp = TmpM[i][i];
	 Tmp = m1.DiagonalElements[i] * m2.DiagonalElements[i];
	 for (int k = 0; k < i; ++k)
	   {
	     Tmp += m1.OffDiagonalElements[m1.GetLinearizedOffDiagonalIndex(i, k)] * m2.OffDiagonalElements[m2.GetLinearizedOffDiagonalIndex(k, i)];
	   }
       }
       for (int j = i + 1; j < m2.NbrColumn; ++j)
	 {
	   Complex& Tmp = TmpM[j][i];
	   Tmp = m1.DiagonalElements[i] * m2.OffDiagonalElements[m2.GetLinearizedOffDiagonalIndex(i, j)];
	   for (int k = 0; k < i; ++k)
	     {
	       Tmp += m1.OffDiagonalElements[m1.GetLinearizedOffDiagonalIndex(i, k)] * m2.OffDiagonalElements[m2.GetLinearizedOffDiagonalIndex(k, j)];
	     }
	 }
     }
  return TmpM;  
}

// multiply a complex upper triangular matrix with a complex lower triangular matrix
//
// m1 = complex upper triangular matrix
// m2 = complex lower triangular matrix
// return value = product result

ComplexMatrix operator * (ComplexUpperTriangularMatrix& m1, ComplexLowerTriangularMatrix& m2)
{
  if (m1.NbrColumn != m2.NbrRow)
    {
      return ComplexMatrix();
    }
  ComplexMatrix TmpM(m1.NbrRow, m2.NbrColumn);
  for (int i = 0; i < m1.NbrRow; ++i)
     {
       for (int j = 0; j < i; ++j)
 	{
 	  Complex& Tmp = TmpM[j][i];
	  Tmp = m1.OffDiagonalElements[m1.GetLinearizedOffDiagonalIndex(i, j)] * m2.DiagonalElements[j];
	  for (int k = 0; k < j; ++k)
	    {
	      Tmp += m1.OffDiagonalElements[m1.GetLinearizedOffDiagonalIndex(i, k)] * m2.OffDiagonalElements[m2.GetLinearizedOffDiagonalIndex(k, j)];
	    }
	}
       {
	 Complex& Tmp = TmpM[i][i];
	 Tmp = m1.DiagonalElements[i] * m2.DiagonalElements[i];
	 for (int k = 0; k < i; ++k)
	   {
	     Tmp += m1.OffDiagonalElements[m1.GetLinearizedOffDiagonalIndex(i, k)] * m2.OffDiagonalElements[m2.GetLinearizedOffDiagonalIndex(k, i)];
	   }
       }
       for (int j = i + 1; j < m2.NbrColumn; ++j)
	 {
	   Complex& Tmp = TmpM[j][i];
	   Tmp = m1.OffDiagonalElements[i] * m2.DiagonalElements[m2.GetLinearizedOffDiagonalIndex(i, j)];
	   for (int k = 0; k < i; ++k)
	     {
	       Tmp += m1.OffDiagonalElements[m1.GetLinearizedOffDiagonalIndex(i, k)] * m2.OffDiagonalElements[m2.GetLinearizedOffDiagonalIndex(k, j)];
	     }
	 }
     }
  return TmpM;  
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexUpperTriangularMatrix operator * (const ComplexUpperTriangularMatrix& M, double x) 
{
  Complex* Diagonal = new Complex [M.NbrRow];
  long TmpNbrOffDiagonalElements = (((long) M.NbrRow) * (((long) M.NbrRow) - 1l)) / 2l;
  Complex* OffDiagonal = new Complex [TmpNbrOffDiagonalElements];
  for (int i = 0; i < M.NbrRow; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
    }
  for (long i = 0l; i < TmpNbrOffDiagonalElements; ++i)
    {
      OffDiagonal[i] = M.OffDiagonalElements[i] * x;
    }
  return ComplexUpperTriangularMatrix(Diagonal, OffDiagonal, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexUpperTriangularMatrix operator * (double x, const ComplexUpperTriangularMatrix& M) 
{
  return (M * x);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

ComplexUpperTriangularMatrix operator / (const ComplexUpperTriangularMatrix& M, double x) 
{
  x = 1.0 / x;
  return (M * x);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

ComplexUpperTriangularMatrix& ComplexUpperTriangularMatrix::operator += (const ComplexUpperTriangularMatrix& M) 
{
  if ((this->NbrRow == 0) || (this->NbrRow != M.NbrRow))
    return *this;
  for (int i = 0; i < M.NbrRow; i++)
    {
      this->DiagonalElements[i] += M.DiagonalElements[i];
    }
  long TmpNbrOffDiagonalElements = (((long) M.NbrRow) * (((long) M.NbrRow) - 1l)) / 2l;
  for (long i = 0l; i < TmpNbrOffDiagonalElements; i++)
    {
      this->OffDiagonalElements[i] += M.OffDiagonalElements[i];
    }
  return *this;
}

// add two matrices where the right one is a real tridiagonal symmetric matrix
//
// M = matrix to add to current matrix
// return value = reference on current matrix

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

ComplexUpperTriangularMatrix& ComplexUpperTriangularMatrix::operator -= (const ComplexUpperTriangularMatrix& M) 
{
  if ((this->NbrRow == 0) || (this->NbrRow != M.NbrRow))
    return *this;
  for (int i = 0; i < M.NbrRow; i++)
    {
      this->DiagonalElements[i] -= M.DiagonalElements[i];
    }
  long TmpNbrOffDiagonalElements = (((long) M.NbrRow) * (((long) M.NbrRow) - 1l)) / 2l;
  for (long i = 0l; i < TmpNbrOffDiagonalElements; i++)
    {
      this->OffDiagonalElements[i] -= M.OffDiagonalElements[i];
    }
  return *this;
}  

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexUpperTriangularMatrix& ComplexUpperTriangularMatrix::operator *= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] *= x;
    }
  long TmpNbrOffDiagonalElements = (((long) this->NbrRow) * (((long) this->NbrRow) - 1l)) / 2l;
  for (long i = 0l; i < TmpNbrOffDiagonalElements; i++)
    {
      this->OffDiagonalElements[i] *= x;
    }
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexUpperTriangularMatrix& ComplexUpperTriangularMatrix::operator /= (double x)
{
  x = 1.0 / x;
  (*this) *= x;
  return *this;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ComplexUpperTriangularMatrix::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  Complex x = 0.0;
  if ((V1.Dimension != this->NbrRow) || (V2.Dimension != this->NbrColumn))
    return x;
  for (int i = 0; i < this->NbrRow ; i++)
    {
      Complex x2 = this->DiagonalElements[i] * V2.Components[i];
      for (int k = i + 1; k < this->NbrColumn; k++)
	{
	  x2 += this->OffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, k)] * V2.Components[k];
	}
      x += x2 * V1.Components[i];
    }
  return x;
}

// Solve the linear equation M x = y
//
// x = vector where the solution will be stored
// y = vector that gives the right hand side of the equation
// return value = true if no error occured

bool ComplexUpperTriangularMatrix::SolveLinearEquation (ComplexVector& x, ComplexVector& y)
{
  if ((this->NbrRow == 0) || (this->NbrRow != x.GetVectorDimension()) || (this->NbrRow != y.GetVectorDimension()))
    return false;
  if ((this->DiagonalElements[this->NbrRow - 1].Re == 0.0) && (this->DiagonalElements[this->NbrRow - 1].Im == 0.0))
    return false;
  x[this->NbrRow - 1] = y[this->NbrRow - 1] / this->DiagonalElements[this->NbrRow - 1];
  for (int i = this->NbrRow - 2; i >= 0; --i)
    {
      if ((this->DiagonalElements[i].Re == 0.0) && (this->DiagonalElements[i].Im == 0.0))
	return false;
      Complex Tmp = y[i];
      for (int j = i + 1 ; j < this->NbrColumn; ++j)
	{
	  Tmp -= this->OffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, j)] * x[j];	  
	}
      x[i] = Tmp / this->DiagonalElements[i];
    }
  return true;
}

// invert the current matrix
//
// return value = true if no error occured

bool ComplexUpperTriangularMatrix::Invert ()
{
  if (this->NbrRow == 0)
    return false;
  for (int i = this->NbrRow - 1; i >= 0; --i)
    {
      if ((this->DiagonalElements[i].Re == 0.0) && (this->DiagonalElements[i].Im == 0.0))
	return false;
      this->DiagonalElements[i] = 1.0 / this->DiagonalElements[i];
      for (int j = i - 1; j >=0 ; --j)
	{
	  Complex Tmp = 0.0;
	  for (int k = j + 1 ; k < this->NbrColumn; ++k)
	    {
	      Tmp += this->OffDiagonalElements[this->GetLinearizedOffDiagonalIndex(j, k)] * this->OffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, k)];	  
	    }
	  this->OffDiagonalElements[this->GetLinearizedOffDiagonalIndex(i, j)] = -Tmp / this->DiagonalElements[j];
	}
    }
  return true;
}

// evaluate matrix trace
//
// return value = matrix trace 

double ComplexUpperTriangularMatrix::Tr () 
{
  return 0.0;
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double ComplexUpperTriangularMatrix::Det () 
{
  return 0.0;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const ComplexUpperTriangularMatrix& P)
{
  Complex TmpZero = 0.0;
  for (int i = 0; i < P.NbrRow; i++)
    {
      for (int j = 0; j < i; j ++)
	{
	  Str << TmpZero << "    ";
	}
      Str << P.DiagonalElements[i] << "    ";
      for (int j = i + 1; j < P.NbrRow; j++)
	{
	  Str << P.OffDiagonalElements[P.GetLinearizedOffDiagonalIndex(i, j)] << "    ";
	}
      Str << endl;
    }
  return Str;
}

#ifdef USE_OUTPUT

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// P = matrix to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexUpperTriangularMatrix& P)
{
  Str << "{";
  for (int i = 0; i < P.NbrRow; ++i)
    {
      Str << "{";
      for (int j = 0; j < i; ++j)
	{
	  Str << "0,";
	}
      Str << P.DiagonalElements[i];
      if (i != (P.NbrRow - 1))
	{
	  Str << ",";	  
	  for (int j = i + 1; j < (P.NbrRow - 1); ++j)
	    {
	      Str << P.OffDiagonalElements[P.GetLinearizedOffDiagonalIndex(i, j)];
	      Str << ",";
	    }
	  Str << P.OffDiagonalElements[P.GetLinearizedOffDiagonalIndex(i, P.NbrRow - 1)];
	  Str << "},";
	}
      else
	Str << "}";
    }
  Str << "}";
  return Str;
}

#endif
