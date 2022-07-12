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


#include "Matrix/RealMatrix.h"
#include "Vector/RealVector.h"

#include <math.h>


using std::endl;
using std::cout;


// default constructor
//

RealMatrix::RealMatrix() 
{
  this->Columns = 0;
  this->ColumnGarbageFlag = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = 0;
  this->TrueNbrColumn = 0;
  this->MatrixType = Matrix::RealElements;
}

// constructor for an empty matrix
//
// nbrRow = number of rows
// nbrColumn = number of columns
// zero = tue if matrix elements have to be set to zero

RealMatrix::RealMatrix(int nbrRow, int nbrColumn, bool zero)
{
  this->ColumnGarbageFlag = new int;
  *(this->ColumnGarbageFlag) = 1;
  this->NbrColumn = nbrColumn;
  this->NbrRow = nbrRow;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Columns = new RealVector [this->NbrColumn];
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] = RealVector (this->NbrRow, zero);
  this->MatrixType = Matrix::RealElements;
}

// constructor from matrix elements (without duplicating datas)
//
// columns = pointer an array of vector
// nbrColumn = number of columns

RealMatrix::RealMatrix(RealVector* columns, int nbrColumn) 
{
  this->Columns = columns;
  this->ColumnGarbageFlag = new int;
  *(this->ColumnGarbageFlag) = 1;
  this->NbrRow = columns[0].GetVectorDimension();
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::RealElements;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

RealMatrix::RealMatrix(const RealMatrix& M) 
{
  if (M.ColumnGarbageFlag != 0)
    {
      this->Columns = M.Columns;
      this->ColumnGarbageFlag = M.ColumnGarbageFlag;
      (*(this->ColumnGarbageFlag))++;
      this->NbrRow = M.NbrRow;
      this->NbrColumn = M.NbrColumn;
      this->TrueNbrRow = M.TrueNbrRow;
      this->TrueNbrColumn = M.TrueNbrColumn;
      this->MatrixType = Matrix::RealElements;
    }
  else
    {
      this->Columns = 0;
      this->ColumnGarbageFlag = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::RealElements;
    }
}

// copy constructor (duplicating all datas)
//
// M = matrix to copy

RealMatrix::RealMatrix(Matrix& M)
{
  if ((M.GetNbrRow() == 0) || (M.GetNbrColumn() == 0))
    {
      this->Columns = 0;
      this->ColumnGarbageFlag = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::RealElements;
    }
  else
    {
      this->ColumnGarbageFlag = new int;
      *(this->ColumnGarbageFlag) = 1;
      this->NbrColumn = M.GetNbrColumn();
      this->NbrRow = M.GetNbrRow();
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;
      this->Columns = new RealVector [this->NbrColumn];
      double Tmp;
      for (int i = 0; i < this->NbrColumn; i++)
	{
	  this->Columns[i] = RealVector (this->NbrRow);
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      M.GetMatrixElement(j, i, Tmp);
	      this->Columns[i][j] = Tmp;
	    }
	}
      this->MatrixType = Matrix::RealElements;
    }
}

// destructor
//

RealMatrix::~RealMatrix() 
{
  if (this->ColumnGarbageFlag != 0)
    {
      if ((*(this->ColumnGarbageFlag)) == 1)
	{
	  delete[] this->Columns;
	  delete this->ColumnGarbageFlag;
	}
      else
	(*(this->ColumnGarbageFlag))--;
    }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

RealMatrix& RealMatrix::operator = (const RealMatrix& M) 
{
  if (this->ColumnGarbageFlag != 0)
    {
      if ((*(this->ColumnGarbageFlag)) == 1)
	{
	  delete[] this->Columns;
	  delete this->ColumnGarbageFlag;
	}
      else
	(*(this->ColumnGarbageFlag))--;
    }
  if (M.ColumnGarbageFlag != 0)
    {
      this->Columns = M.Columns;
      this->ColumnGarbageFlag = M.ColumnGarbageFlag;
      (*(this->ColumnGarbageFlag))++;
      this->NbrRow = M.NbrRow;
      this->NbrColumn = M.NbrColumn;
      this->TrueNbrRow = M.TrueNbrRow;
      this->TrueNbrColumn = M.TrueNbrColumn;
      this->MatrixType = Matrix::RealElements;
    }
  else
    {
      this->Columns = 0;
      this->ColumnGarbageFlag = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::RealElements;
    }
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* RealMatrix::Clone ()
{
  return ((Matrix*) new RealMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void RealMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  this->Columns[j].Components[i] = x;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element
void RealMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void RealMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  this->Columns[j].Components[i] += x;
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void RealMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void RealMatrix::Resize (int nbrRow, int nbrColumn)
{
  if (this->NbrRow != nbrRow)
    {
      for (int i = 0; i < this->NbrColumn; i++)
	{
	  this->Columns[i].Resize(nbrRow);
	}
      if (this->TrueNbrRow >= nbrRow)
	{
	  this->NbrRow = nbrRow;
	}
      else
	{
	  this->NbrRow = nbrRow;
	  this->TrueNbrRow = nbrRow;
	}
    }
  if (this->TrueNbrColumn >= nbrColumn)
    {
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	this->Columns[i].Resize(nbrRow);
      this->NbrColumn = nbrColumn;
    }
  else
    {
      RealVector* Tmp = new RealVector[nbrColumn];
      for (int i = 0; i < this->NbrColumn; i++)
	Tmp[i] = this->Columns[i];      
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	Tmp[i] = RealVector(nbrRow);
      delete[] this->Columns;
      this->Columns = Tmp;
      this->TrueNbrColumn = nbrColumn;
      this->NbrColumn = nbrColumn;
    }
  return;
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void RealMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (this->NbrRow != nbrRow)
    {
      for (int i = 0; i < this->NbrColumn; i++)
	this->Columns[i].ResizeAndClean(nbrRow);
      if (this->TrueNbrRow >= nbrRow)
	{
	  this->NbrRow = nbrRow;
	}
      else
	{
	  this->NbrRow = nbrRow;
	  this->TrueNbrRow = nbrRow;
	}
    }
  if (this->TrueNbrColumn >= nbrColumn)
    {
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	this->Columns[i].ResizeAndClean(nbrRow);
      this->TrueNbrColumn = nbrColumn;
    }
  else
    {
      RealVector* Tmp = new RealVector[nbrColumn];
      for (int i = 0; i < this->NbrColumn; i++)
	Tmp[i] = this->Columns[i];      
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	Tmp[i] = RealVector(nbrRow, true);
      delete[] this->Columns;
      this->Columns = Tmp;
      this->TrueNbrColumn = nbrColumn;
      this->NbrColumn = nbrColumn;
    }
  return;
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

RealMatrix operator + (const RealMatrix& M1, const RealMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return RealMatrix();
  RealVector* TmpColumns = new RealVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; i++)
    {
      TmpColumns[i] = RealVector (M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	TmpColumns[i].Components[j] = M1.Columns[i].Components[j] + M2.Columns[i].Components[j];
    }
  return RealMatrix(TmpColumns, M1.NbrColumn);
}

// add two matrices where the left one is a real tridiagonal symmetric matrix
//
// M1 = left matrix
// M2 = right matrix
// return value = sum of the two matrices

RealMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const RealMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return RealMatrix();
  RealVector* TmpColumns = new RealVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; i++)
    {
      TmpColumns[i] = RealVector(M1.NbrRow);
      int j = 0;
      for (; j < (i - 1); j++)
	TmpColumns[i].Components[j] = M2.Columns[i].Components[j];
      if (i > 0)
	{
	  TmpColumns[i].Components[j] = M1.UpperDiagonalElements[i - 1] + M2.Columns[i].Components[j];
	  j++;
	}
      TmpColumns[i].Components[j] = M1.DiagonalElements[i] + M2.Columns[i].Components[j];
      j++;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].Components[j] = M1.UpperDiagonalElements[i + 1] + M2.Columns[i].Components[j];
	  j++;
	}
      j++;
      for (; j < M1.NbrColumn; j++)
	TmpColumns[i].Components[j] = M2.Columns[i].Components[j];	
    }
  return RealMatrix(TmpColumns, M1.NbrColumn);
}

// add two matrices where the right one is a real tridiagonal symmetric matrix
//
// M1 = left matrix
// M2 = right matrix
// return value = sum of the two matrices

RealMatrix operator + (const RealMatrix& M1, 
				const RealTriDiagonalSymmetricMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return RealMatrix();
  RealVector* TmpColumns = new RealVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; i++)
    {
      TmpColumns[i] = RealVector(M1.NbrRow);
      int j = 0;
      for (; j < (i - 1); j++)
	TmpColumns[i].Components[j] = M1.Columns[i].Components[j];
      if (i > 0)
	{
	  TmpColumns[i].Components[j] = M2.UpperDiagonalElements[i - 1] + M1.Columns[i].Components[j];
	  j++;
	}
      TmpColumns[i].Components[j] = M2.DiagonalElements[i] + M1.Columns[i].Components[j];
      j++;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].Components[j] = M2.UpperDiagonalElements[i + 1] + M1.Columns[i].Components[j];
	  j++;
	}
      j++;
      for (; j < M1.NbrColumn; j++)
	TmpColumns[i].Components[j] = M1.Columns[i].Components[j];	
    }
  return RealMatrix(TmpColumns, M1.NbrColumn);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

RealMatrix operator - (const RealMatrix& M1, const RealMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return RealMatrix();
  RealVector* TmpColumns = new RealVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; i++)
    {
      TmpColumns[i] = RealVector(M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	TmpColumns[i].Components[j] = M1.Columns[i].Components[j] - M2.Columns[i].Components[j];
    }
  return RealMatrix(TmpColumns, M1.NbrColumn);
}

// substract two matrices where the left one is a real tridiagonal symmetric matrix
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

RealMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, 
				const RealMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return RealMatrix();
  RealVector* TmpColumns = new RealVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; i++)
    {
      TmpColumns[i] = RealVector(M1.NbrRow);
      int j = 0;
      for (; j < (i - 1); j++)
	TmpColumns[i].Components[j] = -M2.Columns[i].Components[j];
      if (i > 0)
	{
	  TmpColumns[i].Components[j] = M1.UpperDiagonalElements[i - 1] - M2.Columns[i].Components[j];
	  j++;
	}
      TmpColumns[i].Components[j] = M1.DiagonalElements[i] - M2.Columns[i].Components[j];
      j++;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].Components[j] = M1.UpperDiagonalElements[i + 1] - M2.Columns[i].Components[j];
	  j++;
	}
      j++;
      for (; j < M1.NbrColumn; j++)
	TmpColumns[i].Components[j] = -M2.Columns[i].Components[j];	
    }
  return RealMatrix(TmpColumns, M1.NbrColumn);
}

// substract two matrices where the right one is a real tridiagonal symmetric matrix
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

RealMatrix operator - (const RealMatrix& M1, 
				const RealTriDiagonalSymmetricMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return RealMatrix();
  RealVector* TmpColumns = new RealVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; i++)
    {
      TmpColumns[i] = RealVector(M1.NbrRow);
      int j = 0;
      for (; j < (i - 1); j++)
	TmpColumns[i].Components[j] = M1.Columns[i].Components[j];
      if (i > 0)
	{
	  TmpColumns[i].Components[j] = M1.Columns[i].Components[j] - M2.UpperDiagonalElements[i - 1];
	  j++;
	}
      TmpColumns[i].Components[j] = M1.Columns[i].Components[j] - M2.DiagonalElements[i];
      j++;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].Components[j] = M1.Columns[i].Components[j] - M2.UpperDiagonalElements[i + 1];
	  j++;
	}
      j++;
      for (; j < M1.NbrColumn; j++)
	TmpColumns[i].Components[j] = M1.Columns[i].Components[j];	
    }
  return RealMatrix(TmpColumns, M1.NbrColumn);
}

// multiply two matrices
//
// M1 = first matrix
// M2 = matrix to multiply to M1
// return value = product of the two matrices

RealMatrix operator * (const RealMatrix& M1, const RealMatrix& M2)
{
  if (M1.NbrColumn != M2.NbrRow)
    return RealMatrix();
  RealVector* TmpColumns = new RealVector [M2.NbrColumn];
  for (int i = 0; i < M2.NbrColumn; i++)
    {
      TmpColumns[i] = RealVector(M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	{
	  TmpColumns[i].Components[j] = 0.0;
	  for (int k = 0; k < M2.NbrRow; k++)	
	    TmpColumns[i].Components[j] += M1.Columns[k].Components[j] * M2.Columns[i].Components[k];
	}
    }
  return RealMatrix(TmpColumns, M2.NbrColumn);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealMatrix operator * (const RealMatrix& M, double x) 
{
  RealVector* TmpColumns = new RealVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = RealVector(M.NbrRow);
      for (int j = 0; j < M.NbrRow; j++)
	TmpColumns[i].Components[j] = M.Columns[i].Components[j] * x;
    }
  return RealMatrix(TmpColumns, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealMatrix operator * (double x, const RealMatrix& M) 
{
  RealVector* TmpColumns = new RealVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = RealVector(M.NbrRow);
      for (int j = 0; j < M.NbrRow; j++)
	TmpColumns[i].Components[j] = M.Columns[i].Components[j] * x;
    }
  return RealMatrix(TmpColumns, M.NbrRow);
}

// multiply a matrix to the right by another matrix without using temporary matrix
//
// M = matrix used as multiplicator
// return value = reference on current matrix

RealMatrix& RealMatrix::Multiply (const RealMatrix& M)
{
  if (M.NbrColumn >  this->TrueNbrColumn)
    {
      int OldNbrColumn = this->NbrColumn;
      this->Resize(this->NbrRow, M.NbrColumn);
      this->Resize(this->NbrRow, OldNbrColumn);
    }
  this->Multiply(M, 0, this->NbrRow);
  this->Resize(this->NbrRow, M.NbrColumn);
  return *this;
}

// multiply a matrix to the right by another matrix without using temporary matrix and in a given range of indices
// beware the matrix is not resized after multiplication in order the operation to be thread safe
//
// M = matrix used as multiplicator
// startLine = starting line in destination matrix
// nbrLine = number of lines to multiply
// return value = reference on current matrix

RealMatrix& RealMatrix::Multiply (const RealMatrix& M, int startLine, int nbrLine)
{
  if ((M.NbrRow != this->NbrColumn) || (M.NbrColumn >  this->TrueNbrColumn))
    return *this;
  int EndLine  = nbrLine + startLine;
  double* TmpElements = new double [this->NbrColumn];
  double Tmp;
  for (int i = startLine; i < EndLine; ++i)
    {
      for (int k = 0; k < this->NbrColumn; ++k)
	TmpElements[k] = this->Columns[k].Components[i];
      for (int j = 0; j < M.NbrColumn; ++j)
	{
	  Tmp = TmpElements[0] * M.Columns[j].Components[0];
	  for (int k = 1; k < this->NbrColumn; ++k)
	    {
	      Tmp += TmpElements[k] * M.Columns[j].Components[k];
	    }    
	  this->Columns[j].Components[i] = Tmp;
	}  
    }
  delete[] TmpElements;
  return *this;
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

RealMatrix operator / (const RealMatrix& M, double x) 
{
  RealVector* TmpColumns = new RealVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = RealVector(M.NbrRow);
      for (int j = 0; j < M.NbrRow; j++)
	TmpColumns[i].Components[j] = M.Columns[i].Components[j] * x;
    }
  return RealMatrix(TmpColumns, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

RealMatrix& RealMatrix::operator += (const RealMatrix& M) 
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow))
    return *this;
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] += M.Columns[i];
  return *this;
}

// add two matrices where the right one is a real tridiagonal symmetric matrix
//
// M = matrix to add to current matrix
// return value = reference on current matrix

RealMatrix& RealMatrix::operator += (const RealTriDiagonalSymmetricMatrix& M) 
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow) || (this->ColumnGarbageFlag == 0))
    return *this;  
  this->Columns[0].Components[0] += M.DiagonalElements[0];
  for (int i = 1; i < this->NbrColumn; i++)
    {
      this->Columns[i].Components[i] += M.DiagonalElements[i];
      this->Columns[i].Components[i - 1] += M.UpperDiagonalElements[i - 1];
      this->Columns[i - 1].Components[i] += M.UpperDiagonalElements[i - 1];
    }
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

RealMatrix& RealMatrix::operator -= (const RealMatrix& M) 
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow))
    return *this;
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] -= M.Columns[i];
  return *this;
}

// substract two matrices where the right one is a real tridiagonal symmetric matrix
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

RealMatrix& RealMatrix::operator -= (const RealTriDiagonalSymmetricMatrix& M) 
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow) || (this->ColumnGarbageFlag == 0))
    return *this;  
  this->Columns[0].Components[0] -= M.DiagonalElements[0];
  for (int i = 1; i < this->NbrColumn; i++)
    {
      this->Columns[i].Components[i] -= M.DiagonalElements[i];
      this->Columns[i].Components[i - 1] -= M.UpperDiagonalElements[i - 1];
      this->Columns[i - 1].Components[i] -= M.UpperDiagonalElements[i - 1];
    }
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

RealMatrix& RealMatrix::operator *= (double x) 
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] *= x;
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

RealMatrix& RealMatrix::operator /= (double x)
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] /= x;
  return *this;
}

// normalize matrix column vectors
//
// return value = reference on current matrix

RealMatrix& RealMatrix::NormalizeColumns ()
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i].Normalize();
  return *this;
}

// transpose matrix
//
// return value = reference on current matrix

RealMatrix& RealMatrix::Transpose ()
{
  if (this->NbrRow != this->NbrColumn)
    return *this; 
  double tmp;
  for (int i = 0; i < this->NbrColumn; i++)
    for (int j = i + 1; j < this->NbrColumn; j++)
      {
	tmp = this->Columns[i].Components[j];
	this->Columns[i].Components[j] = this->Columns[j].Components[i];
	this->Columns[j].Components[i] = tmp;;
      }
  return *this;
}

// evaluate matrix determinant (skrewing up matrix elements)
//
// return value = matrix determinant 

double RealMatrix::Determinant () 
{
  if (this->NbrColumn != this->NbrRow)
    return 0.0;
  double TmpDet = 1.0;
  int ReducedNbrRow = this->NbrRow - 1;
  double Pivot;
  double Factor;
  int PivotPos = 0;
  for (int k = 0; k < ReducedNbrRow; ++k)
    {
      Pivot = fabs(this->Columns[k][k]);
      PivotPos = k + 1;
      while ((PivotPos < this->NbrRow) && (fabs(this->Columns[PivotPos][k]) < Pivot))
	{
	  ++PivotPos;
	}
      if (PivotPos == this->NbrRow)
	{
	  Pivot = this->Columns[k][k];	  
	  if (Pivot == 0.0)
	    return 0.0;
	}
      else
	{
	  Pivot = this->Columns[PivotPos][k];	  
	  RealVector TmpColumn3(this->Columns[k]);
	  this->Columns[k] = this->Columns[PivotPos];
	  this->Columns[PivotPos] = TmpColumn3;	  
	  TmpDet *= -1.0;
	}
      TmpDet *= Pivot;
      Pivot = 1.0 / Pivot;       
      for (int i = k + 1; i < this->NbrRow; ++i)
	{
	  RealVector& TmpColumn = this->Columns[i];
	  RealVector& TmpColumn2 = this->Columns[k];
	  Factor = Pivot * TmpColumn[k];
	  for (int j = k + 1; j < this->NbrRow; ++j)
	    {
	      TmpColumn[j] -= TmpColumn2[j] * Factor;
	    }
	}
    }
  TmpDet *= this->Columns[ReducedNbrRow][ReducedNbrRow];
  return TmpDet;
}

// evaluate permanent associated to the (square) matrix using Ryser algorithm
//
// return value = permanent associated to the matrix

double RealMatrix::Permanent()
{
  if (this->NbrColumn != this->NbrRow)
    return 0.0;
  double Perm = 0.0;
  double Sign = 1.0;
  if ((this->NbrColumn & 1) == 0)
    Sign = -1.0;
  double* Tmp = new double [this->NbrColumn];
  double Tmp2;
  int Lim = 1 << this->NbrColumn;
  for (int i = 0; i < this->NbrColumn; ++i)
    Tmp[i] = 0.0;
  int GrayCode = 0;
  int ChangedBit;
  int Index;
  for (int k = 1; k < Lim; ++k)
    {
      ChangedBit = (k ^ (k >> 1)) ^ GrayCode;
      GrayCode = k ^ (k >> 1);
      if ((GrayCode & ChangedBit) == 0)
	{
	  Index = 0;
	  while (ChangedBit != 1)
	    {
	      ChangedBit >>= 1;
	      ++Index;
	    }
	  for (int i = 0; i < this->NbrColumn; ++i)
	    Tmp[i] -= this->Columns[Index].Components[i];
	}
      else
	{
	  Index = 0;
	  while (ChangedBit != 1)
	    {
	      ChangedBit >>= 1;
	      ++Index;
	    }
	  for (int i = 0; i < this->NbrColumn; ++i)
	    Tmp[i] += this->Columns[Index].Components[i];
	}
      Tmp2 = Tmp[0];
      for (int i = 1; i < this->NbrColumn; ++i)
	Tmp2 *= Tmp[i];
      Perm += Sign * Tmp2;
      Sign *= -1.0;
    }  
  delete[] Tmp;
  return Perm;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const RealMatrix& P) 
{
  for (int i = 0; i < (P.NbrRow - 1); i++)
    {
      for (int j = 0; j < (P.NbrColumn - 1); j ++)
	Str << P.Columns[j].Components[i] << "    ";      
      Str << P.Columns[P.NbrColumn - 1].Components[i] << endl;      
    }
  for (int j = 0; j < (P.NbrColumn - 1); j ++)
    Str << P.Columns[j].Components[P.NbrRow - 1] << "    ";      
  Str << P.Columns[P.NbrColumn - 1].Components[P.NbrRow - 1] << endl;      
  return Str;
}

#ifdef USE_OUTPUT

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// P = matrix to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, const RealMatrix& P) 
{
  Str << "[";
  for (int i = 0; i < (P.NbrRow - 1); i++)
    {
      Str << "[";
      for (int j = 0; j < (P.NbrColumn - 1); j ++)
	Str << P.Columns[j].Components[i] << ",";      
      Str << P.Columns[P.NbrColumn - 1].Components[i] << "],";      
    }
  Str << "[";
  for (int j = 0; j < (P.NbrColumn - 1); j ++)
    Str << P.Columns[j].Components[P.NbrRow - 1] << ",";      
  Str << P.Columns[P.NbrColumn - 1].Components[P.NbrRow - 1] << "]";      
  Str << "]";
  return Str;
}

#endif
