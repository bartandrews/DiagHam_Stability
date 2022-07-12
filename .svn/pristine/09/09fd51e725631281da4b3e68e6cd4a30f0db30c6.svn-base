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


#include "Matrix/LongIntegerMatrix.h"
#include "Vector/LongIntegerVector.h"
#include "Architecture/ArchitectureOperation/LongIntegerMatrixCharacteristicPolynomialOperation.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>


using std::endl;
using std::cout;
using std::ios;


// default constructor
//

LongIntegerMatrix::LongIntegerMatrix() 
{
  this->Columns = 0;
  this->ColumnGarbageFlag = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = 0;
  this->TrueNbrColumn = 0;
  this->MatrixType = Matrix::LongIntegerElements;
}

// constructor for an empty matrix
//
// nbrRow = number of rows
// nbrColumn = number of columns
// zero = tue if matrix elements have to be set to zero

LongIntegerMatrix::LongIntegerMatrix(int nbrRow, int nbrColumn, bool zero)
{
  this->ColumnGarbageFlag = new int;
  *(this->ColumnGarbageFlag) = 1;
  this->NbrColumn = nbrColumn;
  this->NbrRow = nbrRow;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Columns = new LongIntegerVector [this->NbrColumn];
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] = LongIntegerVector (this->NbrRow, zero);
  this->MatrixType = Matrix::LongIntegerElements;
}

// constructor from matrix elements (without duplicating datas)
//
// columns = pointer an array of vector
// nbrColumn = number of columns

LongIntegerMatrix::LongIntegerMatrix(LongIntegerVector* columns, int nbrColumn) 
{
  this->Columns = columns;
  this->ColumnGarbageFlag = new int;
  *(this->ColumnGarbageFlag) = 1;
  this->NbrRow = columns[0].GetVectorDimension();
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::LongIntegerElements;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

LongIntegerMatrix::LongIntegerMatrix(const LongIntegerMatrix& M) 
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
      this->MatrixType = Matrix::LongIntegerElements;
    }
  else
    {
      this->Columns = 0;
      this->ColumnGarbageFlag = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::LongIntegerElements;
    }
}

// copy constructor from a real matrrix (duplicating data)
//
// M = matrix to copy
// scalingFactor = scaling factor to apply to M before casting it into an integer matrix

LongIntegerMatrix::LongIntegerMatrix(const Matrix& M, double scalingFactor)
{
  this->MatrixType = Matrix::LongIntegerElements;
  this->NbrColumn = M.GetNbrColumn();
  this->NbrRow = M.GetNbrRow();
  if ((this->NbrRow != 0) && ( this->NbrColumn != 0))
    {
      this->ColumnGarbageFlag = new int;
      *(this->ColumnGarbageFlag) = 1;
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;
      this->Columns = new LongIntegerVector [this->NbrColumn];
      double Tmp;
      for (int i = 0; i < this->NbrColumn; i++)
	{
	  this->Columns[i] = LongIntegerVector (this->NbrRow);
	  for (int j = 0; j < this->NbrRow; ++j)
	    {
	      M.GetMatrixElement(j, i, Tmp);
	      Tmp *= scalingFactor;
	      if (fabs(nearbyint(Tmp) - Tmp) > 1e-10)
		{
		  cout << "error when converting matrix to int (" << j << "," << i << "): " << Tmp << endl;
		}
	      else
		{
		  this->SetMatrixElement(j, i, lrint(Tmp));
		}
	    }
	}
    }
  else
    {
      this->Columns = 0;
      this->ColumnGarbageFlag = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::LongIntegerElements;
    }
}

// destructor
//

LongIntegerMatrix::~LongIntegerMatrix() 
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

LongIntegerMatrix& LongIntegerMatrix::operator = (const LongIntegerMatrix& M) 
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
      this->MatrixType = Matrix::LongIntegerElements;
    }
  else
    {
      this->Columns = 0;
      this->ColumnGarbageFlag = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::LongIntegerElements;
    }
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* LongIntegerMatrix::Clone ()
{
  return ((Matrix*) new LongIntegerMatrix (*this));
}

// copy a matrix into another (duplicating data)
//
// matrix = matrix to copy
// return value = reference on current matrix

LongIntegerMatrix& LongIntegerMatrix::Copy (LongIntegerMatrix& matrix)
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
      this->MatrixType = Matrix::LongIntegerElements;
      this->Columns = new LongIntegerVector[this->NbrColumn];
      for (int i = 0; i < this->NbrColumn; i++)
	{
	  this->Columns[i].Copy(matrix.Columns[i]);
	}
    }
  else
    {
      this->Columns = 0;
      this->ColumnGarbageFlag = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::LongIntegerElements;
    }
  return *this;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void LongIntegerMatrix::SetMatrixElement(int i, int j, const long& x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
#ifdef __GMP__
  mpz_set_si(this->Columns[j].Components[i], x);
#else
  this->Columns[j].Components[i] = (LONGLONG) x;
#endif  
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void LongIntegerMatrix::AddToMatrixElement(int i, int j, const long& x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
#ifdef __GMP__
  if (x >= 0)
    {
      mpz_add_ui(this->Columns[j].Components[i], this->Columns[j].Components[i], (unsigned long) x);
    }
  else
    {
      unsigned long Tmp = (unsigned long) (-x);
      mpz_sub_ui(this->Columns[j].Components[i], this->Columns[j].Components[i], (unsigned long) Tmp);
    }
#else  
  this->Columns[j].Components[i] += x;
#endif
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void LongIntegerMatrix::Resize (int nbrRow, int nbrColumn)
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
      LongIntegerVector* Tmp = new LongIntegerVector[nbrColumn];
      for (int i = 0; i < this->NbrColumn; i++)
	Tmp[i] = this->Columns[i];      
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	Tmp[i] = LongIntegerVector(nbrRow);
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

void LongIntegerMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
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
      LongIntegerVector* Tmp = new LongIntegerVector[nbrColumn];
      for (int i = 0; i < this->NbrColumn; i++)
	Tmp[i] = this->Columns[i];      
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	Tmp[i] = LongIntegerVector(nbrRow, true);
      delete[] this->Columns;
      this->Columns = Tmp;
      this->TrueNbrColumn = nbrColumn;
      this->NbrColumn = nbrColumn;
    }
  return;
}

// Set all entries in matrix to zero
//

void LongIntegerMatrix::ClearMatrix ()
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i].ClearVector();
  return;
}

// set matrix to identity 
//

void LongIntegerMatrix::SetToIdentity()
{
  this->ClearMatrix();
  if (this->NbrColumn <= this->NbrRow)
    {
      for (int i = 0; i < this->NbrColumn; i++)
	{
#ifdef __GMP__
	    mpz_set_ui(this->Columns[i][i], 1ul);
#else
	    this->Columns[i][i] = (LONGLONG) 1l;
#endif
	}
    }
  else
    {
      for (int i = 0; i < this->NbrRow; i++)
	{
#ifdef __GMP__
	    mpz_set_ui(this->Columns[i][i], 1ul);
#else
	    this->Columns[i][i] = (LONGLONG) 1l;
#endif
	}
    }
}


// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

LongIntegerMatrix operator + (const LongIntegerMatrix& M1, const LongIntegerMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return LongIntegerMatrix();
  LongIntegerVector* TmpColumns = new LongIntegerVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; i++)
    {
      TmpColumns[i] = LongIntegerVector (M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	{
#ifdef __GMP__
	  mpz_add(TmpColumns[i][j], M1.Columns[i][j], M2.Columns[i][j]);
#else
	  TmpColumns[i][j] = M1.Columns[i][j] + M2.Columns[i][j];
#endif
	}
    }
  return LongIntegerMatrix(TmpColumns, M1.NbrColumn);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

LongIntegerMatrix operator - (const LongIntegerMatrix& M1, const LongIntegerMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return LongIntegerMatrix();
  LongIntegerVector* TmpColumns = new LongIntegerVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; i++)
    {
      TmpColumns[i] = LongIntegerVector(M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	{
#ifdef __GMP__
	  mpz_sub(TmpColumns[i][j], M1.Columns[i][j], M2.Columns[i][j]);
#else
	  TmpColumns[i][j] = M1.Columns[i][j] - M2.Columns[i][j];
#endif
	}
    }
  return LongIntegerMatrix(TmpColumns, M1.NbrColumn);
}

// multiply two matrices
//
// M1 = first matrix
// M2 = matrix to multiply to M1
// return value = product of the two matrices

LongIntegerMatrix operator * (const LongIntegerMatrix& M1, const LongIntegerMatrix& M2)
{
  if (M1.NbrColumn != M2.NbrRow)
    return LongIntegerMatrix();
  LongIntegerVector* TmpColumns = new LongIntegerVector [M2.NbrColumn];
  for (int i = 0; i < M2.NbrColumn; i++)
    {
      TmpColumns[i] = LongIntegerVector(M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	{
#ifdef __GMP__
	  for (int k = 0; k < M2.NbrRow; k++)
	    {
	      mpz_addmul(TmpColumns[i][j], M1.Columns[k][j], M2.Columns[i][k]);
	    }
#else	  
	  TmpColumns[i][j] = (LONGLONG) 0l;
	  for (int k = 0; k < M2.NbrRow; k++)
	    {
	      TmpColumns[i][j] += M1.Columns[k][j] * M2.Columns[i][k];
	    }
#endif	  
	}
    }
  return LongIntegerMatrix(TmpColumns, M2.NbrColumn);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

LongIntegerMatrix operator * (const LongIntegerMatrix& M, const long& x) 
{
  LongIntegerVector* TmpColumns = new LongIntegerVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = LongIntegerVector(M.NbrRow);
      for (int j = 0; j < M.NbrRow; j++)
	{
#ifdef __GMP__
	  mpz_mul_si(TmpColumns[i][j], M.Columns[i][j], x);
#else  
	  TmpColumns[i][j] = M.Columns[i][j] * ((LONGLONG) x);
#endif      
	}
    }
  return LongIntegerMatrix(TmpColumns, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

LongIntegerMatrix operator * (const long& x, const LongIntegerMatrix& M) 
{
  LongIntegerVector* TmpColumns = new LongIntegerVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = LongIntegerVector(M.NbrRow);
      for (int j = 0; j < M.NbrRow; j++)
	{
#ifdef __GMP__
	  mpz_mul_si(TmpColumns[i][j], M.Columns[i][j], x);
#else  
	  TmpColumns[i][j] = M.Columns[i][j] * ((LONGLONG) x);
#endif      
	}
    }
  return LongIntegerMatrix(TmpColumns, M.NbrRow);
}

// multiply a matrix to the right by another matrix without using temporary matrix
//
// M = matrix used as multiplicator
// return value = reference on current matrix

LongIntegerMatrix& LongIntegerMatrix::Multiply (const LongIntegerMatrix& M)
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

LongIntegerMatrix& LongIntegerMatrix::Multiply (const LongIntegerMatrix& M, int startLine, int nbrLine)
{
  if ((M.NbrRow != this->NbrColumn) || (M.NbrColumn >  this->TrueNbrColumn))
    return *this;
  int EndLine  = nbrLine + startLine;
#ifdef __GMP__
  mpz_t* TmpElements = new mpz_t [this->NbrColumn];
  mpz_t Tmp;
  for (int k = 0; k < this->NbrColumn; ++k)
    {
      mpz_init(TmpElements[k]);
    }
  mpz_init(Tmp);
  for (int i = startLine; i < EndLine; ++i)
    {
      for (int k = 0; k < this->NbrColumn; ++k)
	{
	  mpz_set(TmpElements[k], this->Columns[k][i]);
	}
      for (int j = 0; j < M.NbrColumn; ++j)
	{
	  mpz_mul(Tmp, TmpElements[0], M.Columns[j][0]);
	  for (int k = 1; k < this->NbrColumn; ++k)
	    {
	      mpz_addmul(Tmp, TmpElements[k], M.Columns[j][k]);
	    }    
	  mpz_set(this->Columns[j].Components[i], Tmp);
	}  
    }
  for (int k = 0; k < this->NbrColumn; ++k)
    {
      mpz_clear(TmpElements[k]);
    }
  mpz_clear(Tmp);
  delete[] TmpElements;
#else  
  LONGLONG* TmpElements = new LONGLONG  [this->NbrColumn];
  LONGLONG Tmp;
  for (int i = startLine; i < EndLine; ++i)
    {
      for (int k = 0; k < this->NbrColumn; ++k)
	{
	  TmpElements[k] = this->Columns[k][i];
	}
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
#endif  
  return *this;
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

LongIntegerMatrix operator / (const LongIntegerMatrix& M, const long& x) 
{
  LongIntegerVector* TmpColumns = new LongIntegerVector [M.NbrColumn];
#ifdef __GMP__
  mpz_t Tmp;
  mpz_init_set_si(Tmp, x);
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = LongIntegerVector(M.NbrRow);
      for (int j = 0; j < M.NbrRow; j++)
	{
	  mpz_divexact(TmpColumns[i][j], M.Columns[i][j], Tmp);
	}
    }
  mpz_clear(Tmp);
#else
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = LongIntegerVector(M.NbrRow);
      for (int j = 0; j < M.NbrRow; j++)
	{
	  TmpColumns[i][j] = M.Columns[i][j] / ((LONGLONG) x);
	}
    }
#endif
  return LongIntegerMatrix(TmpColumns, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

LongIntegerMatrix& LongIntegerMatrix::operator += (const LongIntegerMatrix& M) 
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

LongIntegerMatrix& LongIntegerMatrix::operator -= (const LongIntegerMatrix& M) 
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

LongIntegerMatrix& LongIntegerMatrix::operator *= (const long& x) 
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] *= x;
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

LongIntegerMatrix& LongIntegerMatrix::operator /= (const long& x)
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] /= x;
  return *this;
}

// transpose matrix
//
// return value = reference on current matrix

LongIntegerMatrix& LongIntegerMatrix::Transpose ()
{
  if (this->NbrRow != this->NbrColumn)
    return *this; 
  long Tmp;
  for (int i = 0; i < this->NbrColumn; i++)
    {
      for (int j = i + 1; j < this->NbrColumn; j++)
	{
#ifdef __GMP__
	  mpz_swap(this->Columns[i].Components[j], this->Columns[j].Components[i]);
#else	  
	  Tmp = this->Columns[i].Components[j];
	  this->Columns[i].Components[j] = this->Columns[j].Components[i];
	  this->Columns[j].Components[i] = Tmp;
#endif
	}
    }
  return *this;
}

// duplicate and transpose a matrix
//
// return value = transposed matrix

LongIntegerMatrix LongIntegerMatrix::DuplicateAndTranspose ()
{
  LongIntegerMatrix TmpMatrix(this->NbrColumn, this->NbrRow);
  for (int i = 0; i < this->NbrRow; ++i)
    for (int j = 0; j < this->NbrColumn; ++j)
      {
#ifdef __GMP__
	mpz_set(TmpMatrix.Columns[i][j], this->Columns[j][i]);
#else
	TmpMatrix.Columns[i][j] = this->Columns[j][i];
#endif
      }
  return TmpMatrix;
}

// evaluate matrix determinant (skrewing up matrix elements)
//
// return value = matrix determinant 

// long LongIntegerMatrix::Determinant () 
// {
//   cout << "LongIntegerMatrix::Determinant is untested and not functional" << endl;
//   if (this->NbrColumn != this->NbrRow)
//     return 0l;
//   long TmpDet = 1.0;
//   int ReducedNbrRow = this->NbrRow - 1;
//   long Pivot;
//   long Factor;
//   int PivotPos = 0;
//   for (int k = 0; k < ReducedNbrRow; ++k)
//     {
//       PivotPos = k + 1;
//       while ((PivotPos < this->NbrRow) && (this->Columns[PivotPos][k] == 0l))
// 	{
// 	  ++PivotPos;
// 	}
//       if (PivotPos == this->NbrRow)
// 	{
// 	  return 0l;
// 	}
//       else
// 	{
// 	  Pivot = this->Columns[PivotPos][k];	  
// 	  LongIntegerVector TmpColumn3(this->Columns[k]);
// 	  this->Columns[k] = this->Columns[PivotPos];
// 	  this->Columns[PivotPos] = TmpColumn3;	  
// 	  TmpDet *= -1l;
// 	}
//       TmpDet *= Pivot;
//       Pivot = 1l / Pivot;       
//       for (int i = k + 1; i < this->NbrRow; ++i)
// 	{
// 	  LongIntegerVector& TmpColumn = this->Columns[i];
// 	  LongIntegerVector& TmpColumn2 = this->Columns[k];
// 	  Factor = Pivot * TmpColumn[k];
// 	  for (int j = k + 1; j < this->NbrRow; ++j)
// 	    {
// 	      TmpColumn[j] -= TmpColumn2[j] * Factor;
// 	    }
// 	}
//     }
//   TmpDet *= this->Columns[ReducedNbrRow][ReducedNbrRow];
//   return TmpDet;
// }

// evaluate matrix rank
//
// accuracy = numerical accuracy used to define linearly dependence 
// return value = rank

// int LongIntegerMatrix::Rank(double accuracy)
// {
//   cout << "LongIntegerMatrix::Rank is untested and not functional" << endl;
//   cout << "dim = " << this->NbrRow << " " <<  this->NbrColumn << endl;
//   cout << (*this) << endl;
//   int ReducedDim = this->NbrColumn;
//   if (ReducedDim > this->NbrRow)
//     ReducedDim = this->NbrRow;
//   --ReducedDim;
//   long Pivot;
//   long Factor;
//   int PivotPos = 0;
//   for (int k = 0; k < ReducedDim; ++k)
//     {
//       PivotPos = k;
//       while ((PivotPos < this->NbrColumn) && (this->Columns[PivotPos][k] == 0l))
// 	{
// 	  ++PivotPos;
// 	}
//       if (PivotPos < this->NbrColumn)
// 	{
// 	  if (PivotPos != k)
// 	    {
// 	      LongIntegerVector TmpColumn3(this->Columns[k]);
// 	      this->Columns[k] = this->Columns[PivotPos];
// 	      this->Columns[PivotPos] = TmpColumn3;	  
// 	    }
// 	  Pivot = 1l / this->Columns[k][k];       
// 	  for (int i = k + 1; i < this->NbrColumn; ++i)
// 	    {
// 	      LongIntegerVector& TmpColumn = this->Columns[i];
// 	      LongIntegerVector& TmpColumn2 = this->Columns[k];
// 	      if (TmpColumn[k] != 0l)
// 		{
// 		  Factor = Pivot * TmpColumn[k];
// 		  for (int j = k; j < this->NbrRow; ++j)
// 		    {
// 		      TmpColumn[j] -= TmpColumn2[j] * Factor;
// 		    }
// 		}
// 	    }
// 	}
//     }
//   int Rank = 0;
//   ++ReducedDim;
//   for (int k = 0; k < ReducedDim; ++k)
//     {
//       bool Flag = true;
//       for (int i = k; (i < this->NbrRow) && (Flag == true); ++i)
// 	Flag = (this->Columns[k][i] == 0l);
//       if (Flag == false)
// 	++Rank;
//     }
//   return Rank;
// }

// evaluate permanent associated to the (square) matrix using Ryser algorithm
//
// return value = permanent associated to the matrix

// long LongIntegerMatrix::Permanent()
// {
//   cout << "permanent is untested for LongIntegerMatrix" << endl;
//   if (this->NbrColumn != this->NbrRow)
//     return 0l;
//   long Perm(0l);
//   long Sign = 1l;
//   if ((this->NbrColumn & 1) == 0)
//     Sign = -1l;
//   long* Tmp = new long [this->NbrColumn];
//   long Tmp2;
//   int Lim = 1 << this->NbrColumn;
//   for (int i = 0; i < this->NbrColumn; ++i)
//     Tmp[i] = 0l;
//   int GrayCode = 0;
//   int ChangedBit;
//   int Index;
//   for (int k = 1; k < Lim; ++k)
//     {
//       ChangedBit = (k ^ (k >> 1)) ^ GrayCode;
//       GrayCode = k ^ (k >> 1);
//       if ((GrayCode & ChangedBit) == 0)
// 	{
// 	  Index = 0;
// 	  while (ChangedBit != 1)
// 	    {
// 	      ChangedBit >>= 1;
// 	      ++Index;
// 	    }
// 	  for (int i = 0; i < this->NbrColumn; ++i)
// 	    Tmp[i] -= this->Columns[Index].Components[i];
// 	}
//       else
// 	{
// 	  Index = 0;
// 	  while (ChangedBit != 1)
// 	    {
// 	      ChangedBit >>= 1;
// 	      ++Index;
// 	    }
// 	  for (int i = 0; i < this->NbrColumn; ++i)
// 	    Tmp[i] += this->Columns[Index].Components[i];
// 	}
//       Tmp2 = Tmp[0];
//       for (int i = 1; i < this->NbrColumn; ++i)
// 	Tmp2 *= Tmp[i];
//       Tmp2 *= Sign;
//       Perm += Tmp2;
//       Sign *= -1l;
//     }  
//   delete[] Tmp;
//   return Perm;
// }

// write matrix in a file 
//
// fileName = name of the file where the matrix has to be stored
// return value = true if no error occurs

bool LongIntegerMatrix::WriteMatrix (char* fileName)
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

bool LongIntegerMatrix::WriteMatrix (ofstream& file)
{
  file << "# LongIntegerMatrix" << endl;
  file << "# Nbr rows = " << this->NbrRow << endl;
  file << "# Nbr columns = " << this->NbrColumn << endl;
  for (int i = 0; i < this->NbrRow; ++i)
    {
      for (int j = 0; j < this->NbrColumn; ++j)
	{
#ifdef __GMP__
	  file << i << " " << j << " " << this->Columns[j][i] << endl;
#else
	  file << i << " " << j << " " << ((long) this->Columns[j][i]) << endl;
#endif	  
	}
    }
  return true;
}

// read matrix from a file 
//
// fileName = name of the file where the matrix has to be read
// return value = true if no error occurs

bool LongIntegerMatrix::ReadMatrix (char* fileName)
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

bool LongIntegerMatrix::ReadMatrix (ifstream& file)
{
  char* TmpBuffer = new char [256];  
  int Count = 0;
  for (int i = 0; (i < 256) && (Count < 3); ++i)
    {  
      file.read(TmpBuffer + i, 1);
      if (TmpBuffer[i] == '\n')
	++Count;
    }
  if (strcasestr(TmpBuffer, "# LongIntegerMatrix") == 0)
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
  this->Columns = new LongIntegerVector [this->NbrColumn];
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] = LongIntegerVector (this->NbrRow, true);
  int TmpColumnIndex;
  int TmpRowIndex;
#ifdef __GMP__
  for (int i = 0; i < this->NbrRow; ++i)
    {
      for (int j = 0; j < this->NbrColumn; ++j)
	{
	  file >> TmpRowIndex;
	  file >> TmpColumnIndex;
	  file >> this->Columns[TmpColumnIndex][TmpRowIndex];
	}
    }
#else
  long Tmp;
  for (int i = 0; i < this->NbrRow; ++i)
    {
      for (int j = 0; j < this->NbrColumn; ++j)
	{
	  file >> TmpRowIndex;
	  file >> TmpColumnIndex;
	  file >> Tmp;
	  this->Columns[TmpColumnIndex][TmpRowIndex] = (LONGLONG) Tmp;
	}
    }
#endif
  return true;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const LongIntegerMatrix& P) 
{
#ifdef __GMP__  
  for (int i = 0; i < (P.NbrRow - 1); i++)
    {
      for (int j = 0; j < (P.NbrColumn - 1); j ++)
	{
	  Str << P.Columns[j][i] << "    ";
	}
      Str << P.Columns[P.NbrColumn - 1][i] << endl;      
    }
  for (int j = 0; j < (P.NbrColumn - 1); j ++)
    {
      Str << P.Columns[j][P.NbrRow - 1] << "    ";
    }  
  Str << P.Columns[P.NbrColumn - 1][P.NbrRow - 1] << endl;
  
  // Str << "[";
  // for (int i = 0; i < (P.NbrRow - 1); i++)
  //   {
  //     Str << "[";
  //     for (int j = 0; j < (P.NbrColumn - 1); j ++)
  // 	Str << P.Columns[j][i] << ",";      
  //     Str << P.Columns[P.NbrColumn - 1][i] << "],";      
  //   }
  // Str << "[";
  // for (int j = 0; j < (P.NbrColumn - 1); j ++)
  //   Str << P.Columns[j][P.NbrRow - 1] << ",";      
  // Str << P.Columns[P.NbrColumn - 1][P.NbrRow - 1] << "]";      
  // Str << "]";
  return Str;

#else
  for (int i = 0; i < (P.NbrRow - 1); i++)
    {
      for (int j = 0; j < (P.NbrColumn - 1); j ++)
	{
	  Str << (long) P.Columns[j][i] << "    ";
	}
      Str << (long) P.Columns[P.NbrColumn - 1][i] << endl;      
    }
  for (int j = 0; j < (P.NbrColumn - 1); j ++)
    {
      Str << (long) P.Columns[j][P.NbrRow - 1] << "    ";
    }
  Str << (long) P.Columns[P.NbrColumn - 1][P.NbrRow - 1] << endl;
#endif  
  return Str;
}


#ifdef __GMP__
// compute the characteristic polynomial using the Faddeev–Le Verrier algorithm
//
// architecture = pointer to the architecture

mpz_t* LongIntegerMatrix::CharacteristicPolynomial(AbstractArchitecture* architecture)
{
  int* TmpNbrMatrixElements = new int[this->NbrRow];
  int** TmpMatrixElementPositions = new int*[this->NbrRow];
  int* TmpMatrixElementPositions2 = new int[this->NbrRow];  
  for (int i = 0; i < this->NbrRow; ++i)
    {
      TmpNbrMatrixElements[i] = 0;
      for (int j = 0; j < this->NbrColumn; ++j)
	{
	  if (mpz_sgn(this->Columns[j][i]) != 0)
	    {
	      TmpMatrixElementPositions2[TmpNbrMatrixElements[i]] = j;
	      TmpNbrMatrixElements[i]++;
	    }	  
	}
      if (TmpNbrMatrixElements[i] > 0)
	{
	  TmpMatrixElementPositions[i] = new int[TmpNbrMatrixElements[i]];
	  for (int j = 0; j < TmpNbrMatrixElements[i]; ++j)
	    {
	      TmpMatrixElementPositions[i][j] =  TmpMatrixElementPositions2[j];
	    }
	}
      else
	{
	  TmpMatrixElementPositions[i] = 0;
	}
    }     
  delete[]  TmpMatrixElementPositions2;
  
  mpz_t* PolynomialCoefficients = 0;
  if (architecture != 0)
    {
      LongIntegerMatrixCharacteristicPolynomialOperation TmpOperation(this, TmpNbrMatrixElements, TmpMatrixElementPositions);
      TmpOperation.ApplyOperation(architecture);
      PolynomialCoefficients = TmpOperation.GetCharacteristicPolynomial();
    }
  else
    {
      PolynomialCoefficients = new mpz_t [this->NbrRow + 1];
      mpz_t TmpTrace;
      
      for (int i = 0; i <= this->NbrRow; ++i)
	{
	  mpz_init(PolynomialCoefficients[i]);
	}
      mpz_init(TmpTrace);
      mpz_set_ui(PolynomialCoefficients[this->NbrRow], 1ul);
      
      
      LongIntegerMatrix TmpMatrix (this->NbrRow, this->NbrColumn);
      TmpMatrix.Copy(*this);
      LongIntegerMatrix TmpMatrix2 (this->NbrRow, this->NbrColumn, true);
      
      this->Trace(TmpTrace);
      mpz_neg(TmpTrace, TmpTrace);
      mpz_set(PolynomialCoefficients[this->NbrRow - 1], TmpTrace);
      
      for (int k = this->NbrRow - 2; k >= 0; --k)
	{
	  for (int i = 0; i < this->NbrRow; ++i)
	    {
	      mpz_add (TmpMatrix.Columns[i][i], TmpMatrix.Columns[i][i], PolynomialCoefficients[k + 1]);
	    }      
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      for (int i = 0; i < this->NbrRow; ++i)
		{
		  mpz_set_ui(TmpMatrix2.Columns[j][i], 0ul);
		  for (int l = 0; l < TmpNbrMatrixElements[i]; ++l)
		    {
		      mpz_addmul(TmpMatrix2.Columns[j][i], this->Columns[TmpMatrixElementPositions[i][l]][i], TmpMatrix.Columns[j][TmpMatrixElementPositions[i][l]]);
		    }
		}	  
	    }
	  LongIntegerMatrix TmpMatrix3 = TmpMatrix2;
	  TmpMatrix2 = TmpMatrix;
	  TmpMatrix = TmpMatrix3;
	  TmpMatrix.Trace(TmpTrace);
	  mpz_divexact_ui(TmpTrace, TmpTrace, (unsigned long) (this->NbrRow - k));
	  mpz_neg(TmpTrace, TmpTrace);
	  mpz_set(PolynomialCoefficients[k], TmpTrace);      
	}
      mpz_clear(TmpTrace);
    }
  for (int i = 0; i < this->NbrRow; ++i)
    {
      if (TmpNbrMatrixElements[i] > 0)
	{
	  delete[] TmpMatrixElementPositions[i];
	}
    }
  delete[] TmpNbrMatrixElements;
  delete[] TmpMatrixElementPositions;
  return PolynomialCoefficients;
}
#else
// compute the characteristic polynomial using the Faddeev–Le Verrier algorithm
//
// architecture = pointer to the architecture

LONGLONG* LongIntegerMatrix::CharacteristicPolynomial(AbstractArchitecture* architecture)
{
  LONGLONG* PolynomialCoefficients = new LONGLONG [this->NbrRow + 1];

  PolynomialCoefficients[this->NbrRow] = 1l;
  
  LongIntegerMatrix TmpMatrix (this->NbrRow, this->NbrColumn);
  TmpMatrix.Copy(*this);
  LongIntegerMatrix TmpMatrix2 (this->NbrRow, this->NbrColumn, true);

  LONGLONG TmpTrace;
  this->Trace(TmpTrace);
  TmpTrace *= (LONGLONG) -1l;
  
  PolynomialCoefficients[this->NbrRow - 1] = TmpTrace;

  for (int k = this->NbrRow - 2; k >= 0; --k)
    {
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  TmpMatrix.Columns[i][i] += PolynomialCoefficients[k + 1];
	}      
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      TmpMatrix2.Columns[j][i] = (LONGLONG) 0l;
	      for (int l = 0; l < this->NbrColumn; ++l)
		{
		  TmpMatrix2.Columns[j][i] = this->Columns[l][i] * TmpMatrix.Columns[j][l];
		}
	    }	  
	}
      LongIntegerMatrix TmpMatrix3 = TmpMatrix2;
      TmpMatrix2 = TmpMatrix;
      TmpMatrix = TmpMatrix3;
      TmpMatrix.Trace(TmpTrace);
      TmpTrace /= -((LONGLONG) (this->NbrRow - k));
      PolynomialCoefficients[k] = TmpTrace;
      
    }
  return PolynomialCoefficients;
}
#endif


#ifdef __GMP__
// compute the characteristic polynomial using the Faddeev–Le Verrier algorithm and assuming a symmetric matrix
//

mpz_t* LongIntegerMatrix::CharacteristicPolynomialAssumingSymmetric()
{
  mpz_t* PolynomialCoefficients = new mpz_t [this->NbrRow + 1];
  mpz_t TmpTrace;

  for (int i = 0; i <= this->NbrRow; ++i)
    {
      mpz_init(PolynomialCoefficients[i]);
    }
  mpz_init(TmpTrace);
  mpz_set_ui(PolynomialCoefficients[this->NbrRow], 1ul);
  
  LongIntegerMatrix TmpMatrix (this->NbrRow, this->NbrColumn);
  TmpMatrix.Copy(*this);
  LongIntegerMatrix TmpMatrix2 (this->NbrRow, this->NbrColumn, true);

  this->Trace(TmpTrace);
  mpz_neg(TmpTrace, TmpTrace);
  mpz_set(PolynomialCoefficients[this->NbrRow - 1], TmpTrace);

  int* TmpNbrMatrixElements = new int[this->NbrRow];
  int** TmpMatrixElementPositions = new int*[this->NbrRow];
  int* TmpMatrixElementPositions2 = new int[this->NbrRow];  
  for (int i = 0; i < this->NbrRow; ++i)
    {
      TmpNbrMatrixElements[i] = 0;
      for (int j = 0; j < this->NbrColumn; ++j)
	{
	  if (mpz_sgn(this->Columns[i][j]) != 0)
	    {
	      TmpMatrixElementPositions2[TmpNbrMatrixElements[i]] = j;
	      TmpNbrMatrixElements[i]++;
	    }	  
	}
      if (TmpNbrMatrixElements[i] > 0)
	{
	  TmpMatrixElementPositions[i] = new int[TmpNbrMatrixElements[i]];
	  for (int j = 0; j < TmpNbrMatrixElements[i]; ++j)
	    {
	      TmpMatrixElementPositions[i][j] =  TmpMatrixElementPositions2[j];
	    }
	}
      else
	{
	  TmpMatrixElementPositions[i] = 0;
	}
    }     
  delete[]  TmpMatrixElementPositions2;
  
  for (int k = this->NbrRow - 2; k >= 0; --k)
    {
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  mpz_add (TmpMatrix.Columns[i][i], TmpMatrix.Columns[i][i], PolynomialCoefficients[k + 1]);
	}      
      // for (int i = 0; i < this->NbrRow; ++i)
      // 	{
      // 	  for (int j = 0; j < this->NbrColumn; ++j)
      // 	    {
      // 	      ScalarProduct(TmpMatrix2.Columns[i][j], this->Columns[i], TmpMatrix.Columns[j]);
      // 	    }	  
      // 	}
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = i; j < this->NbrColumn; ++j)
	    {
	      mpz_set_ui(TmpMatrix2.Columns[i][j], 0ul);
	      for (int l = 0; l < TmpNbrMatrixElements[i]; ++l)
		{
		  mpz_addmul(TmpMatrix2.Columns[i][j], this->Columns[i][TmpMatrixElementPositions[i][l]], TmpMatrix.Columns[j][TmpMatrixElementPositions[i][l]]);
		}
	    }	  
	}
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < i; ++j)
	    {
	      mpz_set(TmpMatrix2.Columns[i][j], TmpMatrix2.Columns[j][i]);
	    }
	}
      LongIntegerMatrix TmpMatrix3 = TmpMatrix2;
      TmpMatrix2 = TmpMatrix;
      TmpMatrix = TmpMatrix3;
      TmpMatrix.Trace(TmpTrace);
      mpz_divexact_ui(TmpTrace, TmpTrace, (unsigned long) (this->NbrRow - k));
      mpz_neg(TmpTrace, TmpTrace);
      mpz_set(PolynomialCoefficients[k], TmpTrace);      
    }
  mpz_clear(TmpTrace);
  for (int i = 0; i < this->NbrRow; ++i)
    {
      if (TmpNbrMatrixElements[i] > 0)
	{
	  delete[] TmpMatrixElementPositions[i];
	}
    }
  delete[] TmpNbrMatrixElements;
  delete[] TmpMatrixElementPositions;
  return PolynomialCoefficients;
}
#else
// compute the characteristic polynomial using the Faddeev–Le Verrier algorithm and assuming a symmetric matrix
//

LONGLONG* LongIntegerMatrix::CharacteristicPolynomialAssumingSymmetric()
{
  LONGLONG* PolynomialCoefficients = new LONGLONG [this->NbrRow + 1];

  PolynomialCoefficients[this->NbrRow] = 1l;
  
  LongIntegerMatrix TmpMatrix (this->NbrRow, this->NbrColumn);
  TmpMatrix.Copy(*this);
  LongIntegerMatrix TmpMatrix2 (this->NbrRow, this->NbrColumn, true);

  LONGLONG TmpTrace;
  this->Trace(TmpTrace);
  TmpTrace *= (LONGLONG) -1l;
  PolynomialCoefficients[this->NbrRow - 1] = TmpTrace;

  for (int k = this->NbrRow - 2; k >= 0; --k)
    {
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  TmpMatrix.Columns[i][i] += PolynomialCoefficients[k + 1];
	}      
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
		{
		  TmpMatrix2.Columns[i][j] = this->Columns[i] * TmpMatrix.Columns[j];
		}
	    }	  
	}
      LongIntegerMatrix TmpMatrix3 = TmpMatrix2;
      TmpMatrix2 = TmpMatrix;
      TmpMatrix = TmpMatrix3;
      TmpMatrix.Trace(TmpTrace);
      TmpTrace /= -((LONGLONG) (this->NbrRow - k));
      PolynomialCoefficients[k] = TmpTrace;      
    }
  return PolynomialCoefficients;
}
#endif

// compute the number of columns equal to a zero vector
//
// return value = number of null columns 

int LongIntegerMatrix::NbrNullColumns()
{
  int TmpNbrZeroColumns = 0;
  for (int i = 0; i < this->NbrColumn; ++i)
    {
      if (this->Columns[i].IsNullVector() == true)
	{
	  TmpNbrZeroColumns++;
	}
    }
  return TmpNbrZeroColumns;
}

