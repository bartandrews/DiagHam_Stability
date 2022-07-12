////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class of integer matrix                         //
//                                                                            //
//                        last modification : 03/01/2022                      //
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


#include "Matrix/IntegerMatrix.h"
#include "Vector/IntegerVector.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>


using std::endl;
using std::cout;
using std::ios;


// default constructor
//

IntegerMatrix::IntegerMatrix() 
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

IntegerMatrix::IntegerMatrix(int nbrRow, int nbrColumn, bool zero)
{
  this->ColumnGarbageFlag = new int;
  *(this->ColumnGarbageFlag) = 1;
  this->NbrColumn = nbrColumn;
  this->NbrRow = nbrRow;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Columns = new IntegerVector [this->NbrColumn];
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] = IntegerVector (this->NbrRow, zero);
  this->MatrixType = Matrix::IntegerElements;
}

// constructor from matrix elements (without duplicating datas)
//
// columns = pointer an array of vector
// nbrColumn = number of columns

IntegerMatrix::IntegerMatrix(IntegerVector* columns, int nbrColumn) 
{
  this->Columns = columns;
  this->ColumnGarbageFlag = new int;
  *(this->ColumnGarbageFlag) = 1;
  this->NbrRow = columns[0].GetVectorDimension();
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::IntegerElements;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

IntegerMatrix::IntegerMatrix(const IntegerMatrix& M) 
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
      this->MatrixType = Matrix::IntegerElements;
    }
  else
    {
      this->Columns = 0;
      this->ColumnGarbageFlag = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::IntegerElements;
    }
}

// destructor
//

IntegerMatrix::~IntegerMatrix() 
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

IntegerMatrix& IntegerMatrix::operator = (const IntegerMatrix& M) 
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
      this->MatrixType = Matrix::IntegerElements;
    }
  else
    {
      this->Columns = 0;
      this->ColumnGarbageFlag = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::IntegerElements;
    }
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* IntegerMatrix::Clone ()
{
  return ((Matrix*) new IntegerMatrix (*this));
}

// copy a matrix into another (duplicating data)
//
// matrix = matrix to copy
// return value = reference on current matrix

IntegerMatrix& IntegerMatrix::Copy (IntegerMatrix& matrix)
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
      this->Columns = new IntegerVector[this->NbrColumn];
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

void IntegerMatrix::SetMatrixElement(int i, int j, const long& x)
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

void IntegerMatrix::AddToMatrixElement(int i, int j, const long& x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  this->Columns[j].Components[i] += x;
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void IntegerMatrix::Resize (int nbrRow, int nbrColumn)
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
      IntegerVector* Tmp = new IntegerVector[nbrColumn];
      for (int i = 0; i < this->NbrColumn; i++)
	Tmp[i] = this->Columns[i];      
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	Tmp[i] = IntegerVector(nbrRow);
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

void IntegerMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
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
      IntegerVector* Tmp = new IntegerVector[nbrColumn];
      for (int i = 0; i < this->NbrColumn; i++)
	Tmp[i] = this->Columns[i];      
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	Tmp[i] = IntegerVector(nbrRow, true);
      delete[] this->Columns;
      this->Columns = Tmp;
      this->TrueNbrColumn = nbrColumn;
      this->NbrColumn = nbrColumn;
    }
  return;
}

// Set all entries in matrix to zero
//

void IntegerMatrix::ClearMatrix ()
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i].ClearVector();
  return;
}

// set matrix to identity 
//

void IntegerMatrix::SetToIdentity()
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

IntegerMatrix operator + (const IntegerMatrix& M1, const IntegerMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return IntegerMatrix();
  IntegerVector* TmpColumns = new IntegerVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; i++)
    {
      TmpColumns[i] = IntegerVector (M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	TmpColumns[i][j] = M1.Columns[i][j] + M2.Columns[i][j];
    }
  return IntegerMatrix(TmpColumns, M1.NbrColumn);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

IntegerMatrix operator - (const IntegerMatrix& M1, const IntegerMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return IntegerMatrix();
  IntegerVector* TmpColumns = new IntegerVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; i++)
    {
      TmpColumns[i] = IntegerVector(M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	TmpColumns[i][j] = M1.Columns[i][j] - M2.Columns[i][j];
    }
  return IntegerMatrix(TmpColumns, M1.NbrColumn);
}

// multiply two matrices
//
// M1 = first matrix
// M2 = matrix to multiply to M1
// return value = product of the two matrices

IntegerMatrix operator * (const IntegerMatrix& M1, const IntegerMatrix& M2)
{
  if (M1.NbrColumn != M2.NbrRow)
    return IntegerMatrix();
  IntegerVector* TmpColumns = new IntegerVector [M2.NbrColumn];
  for (int i = 0; i < M2.NbrColumn; i++)
    {
      TmpColumns[i] = IntegerVector(M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	{
	  TmpColumns[i][j] = 0l;
	  for (int k = 0; k < M2.NbrRow; k++)	
	    TmpColumns[i][j] += M1.Columns[k][j] * M2.Columns[i][k];
	}
    }
  return IntegerMatrix(TmpColumns, M2.NbrColumn);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

IntegerMatrix operator * (const IntegerMatrix& M, const long& x) 
{
  IntegerVector* TmpColumns = new IntegerVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = IntegerVector(M.NbrRow);
      for (int j = 0; j < M.NbrRow; j++)
	TmpColumns[i][j] = M.Columns[i][j] * x;
    }
  return IntegerMatrix(TmpColumns, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

IntegerMatrix operator * (const long& x, const IntegerMatrix& M) 
{
  IntegerVector* TmpColumns = new IntegerVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = IntegerVector(M.NbrRow);
      for (int j = 0; j < M.NbrRow; j++)
	TmpColumns[i][j] = M.Columns[i][j] * x;
    }
  return IntegerMatrix(TmpColumns, M.NbrRow);
}

// multiply a matrix to the right by another matrix without using temporary matrix
//
// M = matrix used as multiplicator
// return value = reference on current matrix

IntegerMatrix& IntegerMatrix::Multiply (const IntegerMatrix& M)
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

IntegerMatrix& IntegerMatrix::Multiply (const IntegerMatrix& M, int startLine, int nbrLine)
{
  if ((M.NbrRow != this->NbrColumn) || (M.NbrColumn >  this->TrueNbrColumn))
    return *this;
  int EndLine  = nbrLine + startLine;
  long* TmpElements = new long [this->NbrColumn];
  long Tmp;
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

IntegerMatrix operator / (const IntegerMatrix& M, const long& x) 
{
  IntegerVector* TmpColumns = new IntegerVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = IntegerVector(M.NbrRow);
      for (int j = 0; j < M.NbrRow; j++)
	TmpColumns[i][j] = M.Columns[i][j] * x;
    }
  return IntegerMatrix(TmpColumns, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

IntegerMatrix& IntegerMatrix::operator += (const IntegerMatrix& M) 
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

IntegerMatrix& IntegerMatrix::operator -= (const IntegerMatrix& M) 
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

IntegerMatrix& IntegerMatrix::operator *= (const long& x) 
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] *= x;
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

IntegerMatrix& IntegerMatrix::operator /= (const long& x)
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] /= x;
  return *this;
}

// transpose matrix
//
// return value = reference on current matrix

IntegerMatrix& IntegerMatrix::Transpose ()
{
  if (this->NbrRow != this->NbrColumn)
    return *this; 
  long Tmp;
  for (int i = 0; i < this->NbrColumn; i++)
    {
      for (int j = i + 1; j < this->NbrColumn; j++)
	{
	  Tmp = this->Columns[i].Components[j];
	  this->Columns[i].Components[j] = this->Columns[j].Components[i];
	  this->Columns[j].Components[i] = Tmp;
	}
    }
  return *this;
}

// duplicate and transpose a matrix
//
// return value = transposed matrix

IntegerMatrix IntegerMatrix::DuplicateAndTranspose ()
{
  IntegerMatrix TmpMatrix(this->NbrColumn, this->NbrRow);
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

long IntegerMatrix::Determinant () 
{
  cout << "IntegerMatrix::Determinant is untested and not functional" << endl;
  if (this->NbrColumn != this->NbrRow)
    return 0l;
  long TmpDet = 1.0;
  int ReducedNbrRow = this->NbrRow - 1;
  long Pivot;
  long Factor;
  int PivotPos = 0;
  for (int k = 0; k < ReducedNbrRow; ++k)
    {
      PivotPos = k + 1;
      while ((PivotPos < this->NbrRow) && (this->Columns[PivotPos][k] == 0l))
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
	  IntegerVector TmpColumn3(this->Columns[k]);
	  this->Columns[k] = this->Columns[PivotPos];
	  this->Columns[PivotPos] = TmpColumn3;	  
	  TmpDet *= -1l;
	}
      TmpDet *= Pivot;
      Pivot = 1l / Pivot;       
      for (int i = k + 1; i < this->NbrRow; ++i)
	{
	  IntegerVector& TmpColumn = this->Columns[i];
	  IntegerVector& TmpColumn2 = this->Columns[k];
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

int IntegerMatrix::Rank(double accuracy)
{
  cout << "IntegerMatrix::Rank is untested and not functional" << endl;
  cout << "dim = " << this->NbrRow << " " <<  this->NbrColumn << endl;
  cout << (*this) << endl;
  int ReducedDim = this->NbrColumn;
  if (ReducedDim > this->NbrRow)
    ReducedDim = this->NbrRow;
  --ReducedDim;
  long Pivot;
  long Factor;
  int PivotPos = 0;
  for (int k = 0; k < ReducedDim; ++k)
    {
      PivotPos = k;
      while ((PivotPos < this->NbrColumn) && (this->Columns[PivotPos][k] == 0l))
	{
	  ++PivotPos;
	}
      if (PivotPos < this->NbrColumn)
	{
	  if (PivotPos != k)
	    {
	      IntegerVector TmpColumn3(this->Columns[k]);
	      this->Columns[k] = this->Columns[PivotPos];
	      this->Columns[PivotPos] = TmpColumn3;	  
	    }
	  Pivot = 1l / this->Columns[k][k];       
	  for (int i = k + 1; i < this->NbrColumn; ++i)
	    {
	      IntegerVector& TmpColumn = this->Columns[i];
	      IntegerVector& TmpColumn2 = this->Columns[k];
	      if (TmpColumn[k] != 0l)
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
	Flag = (this->Columns[k][i] == 0l);
      if (Flag == false)
	++Rank;
    }
  return Rank;
}

// evaluate permanent associated to the (square) matrix using Ryser algorithm
//
// return value = permanent associated to the matrix

long IntegerMatrix::Permanent()
{
  cout << "permanent is untested for IntegerMatrix" << endl;
  if (this->NbrColumn != this->NbrRow)
    return 0l;
  long Perm(0l);
  long Sign = 1l;
  if ((this->NbrColumn & 1) == 0)
    Sign = -1l;
  long* Tmp = new long [this->NbrColumn];
  long Tmp2;
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

bool IntegerMatrix::WriteMatrix (char* fileName)
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

bool IntegerMatrix::WriteMatrix (ofstream& file)
{
  file << "# IntegerMatrix" << endl;
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

bool IntegerMatrix::ReadMatrix (char* fileName)
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

bool IntegerMatrix::ReadMatrix (ifstream& file)
{
  char* TmpBuffer = new char [256];  
  int Count = 0;
  for (int i = 0; (i < 256) && (Count < 3); ++i)
    {  
      file.read(TmpBuffer + i, 1);
      if (TmpBuffer[i] == '\n')
	++Count;
    }
  if (strcasestr(TmpBuffer, "# IntegerMatrix") == 0)
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
  this->Columns = new IntegerVector [this->NbrColumn];
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] = IntegerVector (this->NbrRow, true);
  int TmpColumnIndex;
  int TmpRowIndex;
  for (int i = 0; i < this->NbrRow; ++i)
    for (int j = 0; j < this->NbrColumn; ++j)
      {
	file >> TmpRowIndex;
	file >> TmpColumnIndex;
	file >> this->Columns[TmpColumnIndex][TmpRowIndex];
      }  
  return true;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const IntegerMatrix& P) 
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


// compute the characteristic polynomial using the Faddeev–Le Verrier algorith and assuming a symmetric matrix
//

long* IntegerMatrix::CharacteristicPolynomialAssumingSymmetric()
{
  long* PolynomialCoefficients = new long [this->NbrRow + 1];

  PolynomialCoefficients[this->NbrRow] = 1l;
  
  IntegerMatrix TmpMatrix (this->NbrRow, this->NbrColumn);
  TmpMatrix.Copy(*this);
  IntegerMatrix TmpMatrix2 (this->NbrRow, this->NbrColumn, true);

  double TmpCoefficient = -this->Trace();
  PolynomialCoefficients[this->NbrRow - 1] = TmpCoefficient;

  for (int k = this->NbrRow - 2; k >= 0; --k)
    {
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  TmpMatrix.Columns[i][i] += PolynomialCoefficients[k + 1];
	}      
      //     cout << TmpMatrix << endl;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      TmpMatrix2.Columns[i][j] = this->Columns[i] * TmpMatrix.Columns[j];
	    }	  
	}
      IntegerMatrix TmpMatrix3 = TmpMatrix2;
      TmpMatrix2 = TmpMatrix;
      TmpMatrix = TmpMatrix3;
      // cout << TmpMatrix << endl;
      // cout << "trace = " << TmpMatrix.Trace() << endl;
      TmpCoefficient = -(TmpMatrix.Trace() / (this->NbrRow - k));
      PolynomialCoefficients[k] = TmpCoefficient;
      
    }
  //  for (int k = this->NbrRow; k >= 0; --k)
  // for (int k = 0; k <= this->NbrRow; ++k)
  //   {
  //     cout << PolynomialCoefficients[k] << endl;
  //   }  
  return PolynomialCoefficients;
}
