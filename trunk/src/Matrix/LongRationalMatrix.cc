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


#include "Matrix/LongRationalMatrix.h"
#include "Vector/LongRationalVector.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>


using std::endl;
using std::cout;
using std::ios;


// default constructor
//

LongRationalMatrix::LongRationalMatrix() 
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

LongRationalMatrix::LongRationalMatrix(int nbrRow, int nbrColumn, bool zero)
{
  this->ColumnGarbageFlag = new int;
  *(this->ColumnGarbageFlag) = 1;
  this->NbrColumn = nbrColumn;
  this->NbrRow = nbrRow;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Columns = new LongRationalVector [this->NbrColumn];
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] = LongRationalVector (this->NbrRow, zero);
  this->MatrixType = Matrix::LongRationalElements;
}

// constructor from matrix elements (without duplicating datas)
//
// columns = pointer an array of vector
// nbrColumn = number of columns

LongRationalMatrix::LongRationalMatrix(LongRationalVector* columns, int nbrColumn) 
{
  this->Columns = columns;
  this->ColumnGarbageFlag = new int;
  *(this->ColumnGarbageFlag) = 1;
  this->NbrRow = columns[0].GetVectorDimension();
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::LongRationalElements;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

LongRationalMatrix::LongRationalMatrix(const LongRationalMatrix& M) 
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
      this->MatrixType = Matrix::LongRationalElements;
    }
  else
    {
      this->Columns = 0;
      this->ColumnGarbageFlag = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::LongRationalElements;
    }
}

// destructor
//

LongRationalMatrix::~LongRationalMatrix() 
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

LongRationalMatrix& LongRationalMatrix::operator = (const LongRationalMatrix& M) 
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
      this->MatrixType = Matrix::LongRationalElements;
    }
  else
    {
      this->Columns = 0;
      this->ColumnGarbageFlag = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::LongRationalElements;
    }
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* LongRationalMatrix::Clone ()
{
  return ((Matrix*) new LongRationalMatrix (*this));
}

// copy a matrix into another (duplicating data)
//
// matrix = matrix to copy
// return value = reference on current matrix

LongRationalMatrix& LongRationalMatrix::Copy (LongRationalMatrix& matrix)
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
  if (matrix.ColumnGarbageFlag != 0)
    {
      this->ColumnGarbageFlag = new int;
      *(this->ColumnGarbageFlag) = 1;
      this->NbrRow = matrix.NbrRow;
      this->NbrColumn = matrix.NbrColumn;
      this->TrueNbrRow = matrix.TrueNbrRow;
      this->TrueNbrColumn = matrix.TrueNbrColumn;
      this->MatrixType = Matrix::RealElements;
      this->Columns = new LongRationalVector[this->NbrColumn];
      for (int i = 0; i < this->NbrColumn; i++)
	this->Columns[i].Copy(matrix.Columns[i]);
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

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void LongRationalMatrix::SetMatrixElement(int i, int j, const LongRational& x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  this->Columns[j].Components[i] = x;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void LongRationalMatrix::AddToMatrixElement(int i, int j, const LongRational& x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  this->Columns[j].Components[i] += x;
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void LongRationalMatrix::Resize (int nbrRow, int nbrColumn)
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
      LongRationalVector* Tmp = new LongRationalVector[nbrColumn];
      for (int i = 0; i < this->NbrColumn; i++)
	Tmp[i] = this->Columns[i];      
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	Tmp[i] = LongRationalVector(nbrRow);
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
      this->ColumnGarbageFlag = new int;
      *(this->ColumnGarbageFlag) = 1;
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

void LongRationalMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
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
      LongRationalVector* Tmp = new LongRationalVector[nbrColumn];
      for (int i = 0; i < this->NbrColumn; i++)
	Tmp[i] = this->Columns[i];      
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	Tmp[i] = LongRationalVector(nbrRow, true);
      delete[] this->Columns;
      this->Columns = Tmp;
      this->TrueNbrColumn = nbrColumn;
      this->NbrColumn = nbrColumn;
    }
  return;
}

// Set all entries in matrix to zero
//

void LongRationalMatrix::ClearMatrix ()
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i].ClearVector();
  return;
}

// set matrix to identity 
//

void LongRationalMatrix::SetToIdentity()
{
  this->ClearMatrix();
  if (this->NbrColumn <= this->NbrRow)
    {
      for (int i = 0; i < this->NbrColumn; i++)
	this->Columns[i][i] = 1l;
    }
  else
    {
      for (int i = 0; i < this->NbrRow; i++)
	this->Columns[i][i] = 1l;
    }
}


// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

LongRationalMatrix operator + (const LongRationalMatrix& M1, const LongRationalMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return LongRationalMatrix();
  LongRationalVector* TmpColumns = new LongRationalVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; i++)
    {
      TmpColumns[i] = LongRationalVector (M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	TmpColumns[i][j] = M1.Columns[i][j] + M2.Columns[i][j];
    }
  return LongRationalMatrix(TmpColumns, M1.NbrColumn);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

LongRationalMatrix operator - (const LongRationalMatrix& M1, const LongRationalMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return LongRationalMatrix();
  LongRationalVector* TmpColumns = new LongRationalVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; i++)
    {
      TmpColumns[i] = LongRationalVector(M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	TmpColumns[i][j] = M1.Columns[i][j] - M2.Columns[i][j];
    }
  return LongRationalMatrix(TmpColumns, M1.NbrColumn);
}

// multiply two matrices
//
// M1 = first matrix
// M2 = matrix to multiply to M1
// return value = product of the two matrices

LongRationalMatrix operator * (const LongRationalMatrix& M1, const LongRationalMatrix& M2)
{
  if (M1.NbrColumn != M2.NbrRow)
    return LongRationalMatrix();
  LongRationalVector* TmpColumns = new LongRationalVector [M2.NbrColumn];
  for (int i = 0; i < M2.NbrColumn; i++)
    {
      TmpColumns[i] = LongRationalVector(M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	{
	  TmpColumns[i][j] = 0l;
	  for (int k = 0; k < M2.NbrRow; k++)	
	    TmpColumns[i][j] += M1.Columns[k][j] * M2.Columns[i][k];
	}
    }
  return LongRationalMatrix(TmpColumns, M2.NbrColumn);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

LongRationalMatrix operator * (const LongRationalMatrix& M, const LongRational& x) 
{
  LongRationalVector* TmpColumns = new LongRationalVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = LongRationalVector(M.NbrRow);
      for (int j = 0; j < M.NbrRow; j++)
	TmpColumns[i][j] = M.Columns[i][j] * x;
    }
  return LongRationalMatrix(TmpColumns, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

LongRationalMatrix operator * (const LongRational& x, const LongRationalMatrix& M) 
{
  LongRationalVector* TmpColumns = new LongRationalVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = LongRationalVector(M.NbrRow);
      for (int j = 0; j < M.NbrRow; j++)
	TmpColumns[i][j] = M.Columns[i][j] * x;
    }
  return LongRationalMatrix(TmpColumns, M.NbrRow);
}

// multiply a matrix to the right by another matrix without using temporary matrix
//
// M = matrix used as multiplicator
// return value = reference on current matrix

LongRationalMatrix& LongRationalMatrix::Multiply (const LongRationalMatrix& M)
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

LongRationalMatrix& LongRationalMatrix::Multiply (const LongRationalMatrix& M, int startLine, int nbrLine)
{
  if ((M.NbrRow != this->NbrColumn) || (M.NbrColumn >  this->TrueNbrColumn))
    return *this;
  int EndLine  = nbrLine + startLine;
  LongRational* TmpElements = new LongRational [this->NbrColumn];
  LongRational Tmp;
  for (int i = startLine; i < EndLine; ++i)
    {
      for (int k = 0; k < this->NbrColumn; ++k)
	TmpElements[k] = this->Columns[k][i];
      for (int j = 0; j < M.NbrColumn; ++j)
	{
	  Tmp = TmpElements[0] * M.Columns[j][0];
	  for (int k = 1; k < this->NbrColumn; ++k)
	    {
	      Tmp += TmpElements[k] * M.Columns[j][k];
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

LongRationalMatrix operator / (const LongRationalMatrix& M, const LongRational& x) 
{
  LongRationalVector* TmpColumns = new LongRationalVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = LongRationalVector(M.NbrRow);
      for (int j = 0; j < M.NbrRow; j++)
	TmpColumns[i][j] = M.Columns[i][j] * x;
    }
  return LongRationalMatrix(TmpColumns, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

LongRationalMatrix& LongRationalMatrix::operator += (const LongRationalMatrix& M) 
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow))
    return *this;
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] += M.Columns[i];
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

LongRationalMatrix& LongRationalMatrix::operator -= (const LongRationalMatrix& M) 
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow))
    return *this;
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] -= M.Columns[i];
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

LongRationalMatrix& LongRationalMatrix::operator *= (const LongRational& x) 
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] *= x;
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

LongRationalMatrix& LongRationalMatrix::operator /= (const LongRational& x)
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] /= x;
  return *this;
}

// transpose matrix
//
// return value = reference on current matrix

LongRationalMatrix& LongRationalMatrix::Transpose ()
{
  if (this->NbrRow != this->NbrColumn)
    return *this; 
  LongRational tmp;
  for (int i = 0; i < this->NbrColumn; i++)
    for (int j = i + 1; j < this->NbrColumn; j++)
      {
	Swap(this->Columns[i].Components[j], this->Columns[j].Components[i]);
      }
  return *this;
}

// duplicate and transpose a matrix
//
// return value = transposed matrix

LongRationalMatrix LongRationalMatrix::DuplicateAndTranspose ()
{
  LongRationalMatrix TmpMatrix(this->NbrColumn, this->NbrRow);
  for (int i = 0; i < this->NbrRow; ++i)
    for (int j = 0; j < this->NbrColumn; ++j)
      {
	TmpMatrix.Columns[i][j] = this->Columns[j][i];
      }
  return TmpMatrix;
}

// evaluate matrix determinant (skrewing up matrix elements)
//
// return value = matrix determinant 

LongRational LongRationalMatrix::Determinant () 
{
  cout << "determinant is untested for LongRationalMatrix" << endl;
  if (this->NbrColumn != this->NbrRow)
    return 0l;
  LongRational TmpDet = 1.0;
  int ReducedNbrRow = this->NbrRow - 1;
  LongRational Pivot;
  LongRational Factor;
  int PivotPos = 0;
  for (int k = 0; k < ReducedNbrRow; ++k)
    {
      PivotPos = k + 1;
      while ((PivotPos < this->NbrRow) && (this->Columns[PivotPos][k].IsZero()))
	{
	  ++PivotPos;
	}
      if (PivotPos == this->NbrRow)
	{
	  return 0l;
	}
      else
	{
	  Pivot = this->Columns[PivotPos][k];	  
	  LongRationalVector TmpColumn3(this->Columns[k]);
	  this->Columns[k] = this->Columns[PivotPos];
	  this->Columns[PivotPos] = TmpColumn3;	  
	  TmpDet *= -1l;
	}
      TmpDet *= Pivot;
      Pivot = 1l / Pivot;       
      for (int i = k + 1; i < this->NbrRow; ++i)
	{
	  LongRationalVector& TmpColumn = this->Columns[i];
	  LongRationalVector& TmpColumn2 = this->Columns[k];
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

// evaluate matrix rank
//
// accuracy = numerical accuracy used to define linearly dependence 
// return value = rank

int LongRationalMatrix::Rank(double accuracy)
{
  cout << "dim = " << this->NbrRow << " " <<  this->NbrColumn << endl;
  cout << (*this) << endl;
  int ReducedDim = this->NbrColumn;
  if (ReducedDim > this->NbrRow)
    ReducedDim = this->NbrRow;
  --ReducedDim;
  LongRational Pivot;
  LongRational Factor;
  int PivotPos = 0;
  for (int k = 0; k < ReducedDim; ++k)
    {
      PivotPos = k;
      while ((PivotPos < this->NbrColumn) && (this->Columns[PivotPos][k].IsZero()))
	{
	  ++PivotPos;
	}
      if (PivotPos < this->NbrColumn)
	{
	  if (PivotPos != k)
	    {
	      LongRationalVector TmpColumn3(this->Columns[k]);
	      this->Columns[k] = this->Columns[PivotPos];
	      this->Columns[PivotPos] = TmpColumn3;	  
	    }
	  Pivot = 1l / this->Columns[k][k];       
	  for (int i = k + 1; i < this->NbrColumn; ++i)
	    {
	      LongRationalVector& TmpColumn = this->Columns[i];
	      LongRationalVector& TmpColumn2 = this->Columns[k];
	      if (TmpColumn[k].IsZero() == false)
		{
		  Factor = Pivot * TmpColumn[k];
		  for (int j = k; j < this->NbrRow; ++j)
		    {
		      TmpColumn[j] -= TmpColumn2[j] * Factor;
		    }
		}
	    }
	}
    }
  int Rank = 0;
  ++ReducedDim;
  for (int k = 0; k < ReducedDim; ++k)
    {
      bool Flag = true;
      for (int i = k; (i < this->NbrRow) && (Flag == true); ++i)
	Flag = this->Columns[k][i].IsZero();
      if (Flag == false)
	++Rank;
    }
//  cout << (*this) << endl;
//  cout << "rank = " << Rank << endl;
  return Rank;
}

// evaluate permanent associated to the (square) matrix using Ryser algorithm
//
// return value = permanent associated to the matrix

LongRational LongRationalMatrix::Permanent()
{
  cout << "permanent is untested for LongRationalMatrix" << endl;
  if (this->NbrColumn != this->NbrRow)
    return 0l;
  LongRational Perm(0l);
  long Sign = 1l;
  if ((this->NbrColumn & 1) == 0)
    Sign = -1l;
  LongRational* Tmp = new LongRational [this->NbrColumn];
  LongRational Tmp2;
  int Lim = 1 << this->NbrColumn;
  for (int i = 0; i < this->NbrColumn; ++i)
    Tmp[i] = 0l;
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
      Tmp2 *= Sign;
      Perm += Tmp2;
      Sign *= -1l;
    }  
  delete[] Tmp;
  return Perm;
}

// write matrix in a file 
//
// fileName = name of the file where the matrix has to be stored
// return value = true if no error occurs

bool LongRationalMatrix::WriteMatrix (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  this->WriteMatrix(File);
  File.close();
  return true;
}

// write matrix in a file 
//
// file = reference on the output file stream
// return value = true if no error occurs

bool LongRationalMatrix::WriteMatrix (ofstream& file)
{
  file << "# LongRationalMatrix" << endl;
  file << "# Nbr rows = " << this->NbrRow << endl;
  file << "# Nbr columns = " << this->NbrColumn << endl;
  for (int i = 0; i < this->NbrRow; ++i)
    for (int j = 0; j < this->NbrColumn; ++j)
      file << i << " " << j << " " << this->Columns[j][i] << endl;  
  return true;
}

// read matrix from a file 
//
// fileName = name of the file where the matrix has to be read
// return value = true if no error occurs

bool LongRationalMatrix::ReadMatrix (char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return false;
    }
  this->ReadMatrix(File);
  File.close();
  return true;
}

// read matrix from a file 
//
// file = reference  on the input file stream
// return value = true if no error occurs

bool LongRationalMatrix::ReadMatrix (ifstream& file)
{
  char* TmpBuffer = new char [256];  
  int Count = 0;
  for (int i = 0; (i < 256) && (Count < 3); ++i)
    {  
      file.read(TmpBuffer + i, 1);
      if (TmpBuffer[i] == '\n')
	++Count;
    }
  if (strcasestr(TmpBuffer, "# LongRationalMatrix") == 0)
    {
      return false;
    }
  char* TmpPos = strcasestr(TmpBuffer, "# Nbr rows = ");
  if (TmpPos == 0)
    return false;
  this->NbrRow = atoi (TmpPos + 13);
  TmpPos = strcasestr(TmpBuffer, "# Nbr columns = ");
  if (TmpPos == 0)
    return false;
  this->NbrColumn = atoi (TmpPos + 16);
  TmpPos = strcasestr(TmpPos, "\n");
  if (TmpPos == 0)
    return false;
  this->ColumnGarbageFlag = new int;
  *(this->ColumnGarbageFlag) = 1;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Columns = new LongRationalVector [this->NbrColumn];
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] = LongRationalVector (this->NbrRow, true);
  int TmpColumnIndex;
  int TmpRowIndex;
  for (int i = 0; i < this->NbrRow; ++i)
    for (int j = 0; j < this->NbrColumn; ++j)
      {
	file >> TmpRowIndex;
	file >> TmpColumnIndex;
	this->Columns[TmpColumnIndex][TmpRowIndex].Read(file);
      }  
  return true;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const LongRationalMatrix& P) 
{
  for (int i = 0; i < (P.NbrRow - 1); i++)
    {
      for (int j = 0; j < (P.NbrColumn - 1); j ++)
	Str << P.Columns[j][i] << "    ";      
      Str << P.Columns[P.NbrColumn - 1][i] << endl;      
    }
  for (int j = 0; j < (P.NbrColumn - 1); j ++)
    Str << P.Columns[j][P.NbrRow - 1] << "    ";      
  Str << P.Columns[P.NbrColumn - 1][P.NbrRow - 1] << endl;      
  return Str;
}


