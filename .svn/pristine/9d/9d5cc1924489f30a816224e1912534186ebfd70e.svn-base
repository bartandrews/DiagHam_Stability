////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of real upper triangular matrix                  //
//                                                                            //
//                        last modification : 07/01/2003                      //
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


#include "Matrix/RealUpperTriangularMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "GeneralTools/ListIterator.h"

#include <stdlib.h>
#include <fstream>


using std::endl;


// default constructor
//

RealUpperTriangularMatrix::RealUpperTriangularMatrix() 
{
  this->DiagonalElements = 0;
  this->OffDiagonalElements = 0;
  this->DiagonalGarbageFlag =  0;
  this->OffDiagonalGarbageFlag =  0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = 0;
  this->MatrixType = Matrix::RealElements | Matrix::Symmetric;
  this->Dummy = 0.0;
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

RealUpperTriangularMatrix::RealUpperTriangularMatrix(int dimension, bool zero) 
{
  this->DiagonalGarbageFlag =  new int;
  *(this->DiagonalGarbageFlag) = 1;
  this->OffDiagonalGarbageFlag =  new int;
  *(this->OffDiagonalGarbageFlag) = 1;
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::RealElements | Matrix::Symmetric;
  this->DiagonalElements = new double [this->NbrRow];
  this->OffDiagonalElements = new double [(this->NbrRow * (this->NbrRow - 1)) / 2];
  if (zero == true)
    {
      int pos = 0;
      for (int i = 0; i < this->NbrRow; i++)
	{
	  this->DiagonalElements[i] = 0.0;
	  for (int j = i + 1; j < this->NbrRow; j++)
	    {
	      this->OffDiagonalElements[pos] = 0.0;
	      pos++;
	    }
	}
    }
  this->Dummy = 0.0;
}

// constructor from matrix elements (without duplicating datas)
//
// diagonal = pointer to diagonal element array
// offDiagonal = pointer to off-diagonal element array (with real part in even position and imaginary part in odd position)
// dimension = matrix dimension

RealUpperTriangularMatrix::RealUpperTriangularMatrix(double* diagonal, double* offDiagonal, int dimension) 
{
  this->DiagonalElements = diagonal;
  this->OffDiagonalElements = offDiagonal;
  this->DiagonalGarbageFlag =  new int;
  *(this->DiagonalGarbageFlag) = 1;
  this->OffDiagonalGarbageFlag =  new int;
  *(this->OffDiagonalGarbageFlag) = 1;
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::RealElements | Matrix::Symmetric;
  this->Dummy = 0.0;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

RealUpperTriangularMatrix::RealUpperTriangularMatrix(const RealUpperTriangularMatrix& M) 
{
  this->DiagonalElements = M.DiagonalElements;
  this->DiagonalGarbageFlag = M.DiagonalGarbageFlag;
  if (this->DiagonalGarbageFlag != 0)
    (*(this->DiagonalGarbageFlag))++;
  this->OffDiagonalElements = M.OffDiagonalElements;
  this->OffDiagonalGarbageFlag = M.OffDiagonalGarbageFlag;
  if (this->OffDiagonalGarbageFlag != 0)
    (*(this->OffDiagonalGarbageFlag))++;  
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::RealElements | Matrix::Symmetric;
  this->Dummy = 0.0;
}

// destructor
//

RealUpperTriangularMatrix::~RealUpperTriangularMatrix() 
{
  if (this->OffDiagonalGarbageFlag != 0)
    {
      if ((*(this->OffDiagonalGarbageFlag)) == 1)
	{
	  delete this->OffDiagonalGarbageFlag;
	  delete[] this->OffDiagonalElements;
	}
      else
	(*(this->OffDiagonalGarbageFlag))--;
    }
  if (this->DiagonalGarbageFlag != 0)
    {
      if ((*(this->DiagonalGarbageFlag)) == 1)
	{
	  delete this->DiagonalGarbageFlag;
	  delete[] this->DiagonalElements;
	}
      else
	(*(this->DiagonalGarbageFlag))--;
    }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

RealUpperTriangularMatrix& RealUpperTriangularMatrix::operator = (const RealUpperTriangularMatrix& M) 
{
  if (this->OffDiagonalGarbageFlag != 0)
    {
      if ((*(this->OffDiagonalGarbageFlag)) == 1)
	{
	  delete this->OffDiagonalGarbageFlag;
	  delete[] this->OffDiagonalElements;
	}
      else
	(*(this->OffDiagonalGarbageFlag))--;
    }
  if (this->DiagonalGarbageFlag != 0)
    {
      if ((*(this->DiagonalGarbageFlag)) == 1)
	{
	  delete this->DiagonalGarbageFlag;
	  delete[] this->DiagonalElements;
	}
      else
	(*(this->DiagonalGarbageFlag))--;
    }
  this->DiagonalElements = M.DiagonalElements;
  this->DiagonalGarbageFlag = M.DiagonalGarbageFlag;
  if (this->DiagonalGarbageFlag != 0)
    (*(this->DiagonalGarbageFlag))++;
  this->OffDiagonalElements = M.OffDiagonalElements;
  this->OffDiagonalGarbageFlag = M.OffDiagonalGarbageFlag;
  if (this->OffDiagonalGarbageFlag != 0)
    (*(this->OffDiagonalGarbageFlag))++;  
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->Increment = this->TrueNbrRow - this->NbrRow;
  this->Dummy = 0.0;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* RealUpperTriangularMatrix::Clone ()
{
  return ((Matrix*) new RealUpperTriangularMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void RealUpperTriangularMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      this->DiagonalElements[i] = x;
    }
  else
    if (i < j)
      {
	i += (j * (j - 1)) / 2;
	this->OffDiagonalElements[i] = x;
      }
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element
void RealUpperTriangularMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  return;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void RealUpperTriangularMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      this->DiagonalElements[i] += x;
    }
  else
    if (i < j)
      {
	i += (j * (j - 1)) / 2;
	this->OffDiagonalElements[i] += x;
      }
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void RealUpperTriangularMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  return;
}

// get reference of a given matrix element
//
// i = line position
// j = column position
// return value = reference om matrix elememt

double& RealUpperTriangularMatrix::operator () (int i, int j)
{
  if (i == j)
    {
      return this->DiagonalElements[i];
    }
  else
    {
      if (i < j)
	{
	  i += (j * (j - 1)) / 2;
	  return this->OffDiagonalElements[i];
	}
      else
	{
	  return this->Dummy;
	}
    }
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void RealUpperTriangularMatrix::Resize (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  if (nbrRow <= this->TrueNbrRow)
    {
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      this->Increment = this->TrueNbrRow - this->NbrRow;
      return;
    }
  double* TmpDiag = new double [nbrRow];
  int Tot = (nbrRow * (nbrRow - 1)) / 2;
  double* TmpOffDiag = new double [Tot];
  for (int i = 0; i < this->NbrRow; i++)
    TmpDiag [i] = this->DiagonalElements[i];
  for (int i = this->NbrRow; i < nbrRow; i++)
    TmpDiag [i]  = 0.0;
  int k = 0;
  int l = 0;
  for (int i = 0; i < (this->NbrRow - 1); i++)
    {
      for (int j = i + 1; j < this->NbrRow; j++)
	TmpOffDiag[k++] = this->OffDiagonalElements[l++];
      l += this->Increment;
      for (int j = this->NbrRow; j < nbrRow; j++)
	{
	  TmpOffDiag[k++] = 0.0;
	}      
    }
  for (int i = this->NbrRow * (this->NbrRow - 1); i < Tot; i++)
    TmpOffDiag[i] = 0.0;
  if (this->OffDiagonalGarbageFlag != 0)
    {
      if ((*(this->OffDiagonalGarbageFlag)) == 1)
	{
	  delete this->OffDiagonalGarbageFlag;
	  delete[] this->OffDiagonalElements;
	}
      else
	(*(this->OffDiagonalGarbageFlag))--;
    }
  if (this->DiagonalGarbageFlag != 0)
    {
      if ((*(this->DiagonalGarbageFlag)) == 1)
	{
	  delete this->DiagonalGarbageFlag;
	  delete[] this->DiagonalElements;
	}
      else
	(*(this->DiagonalGarbageFlag))--;
    }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = this->TrueNbrRow - this->NbrRow;
  this->DiagonalElements = TmpDiag;
  this->OffDiagonalElements = TmpOffDiag;
  this->DiagonalGarbageFlag =  new int;
  *(this->DiagonalGarbageFlag) = 1;
  this->OffDiagonalGarbageFlag =  new int;
  *(this->OffDiagonalGarbageFlag) = 1;
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void RealUpperTriangularMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  if (nbrRow <= this->TrueNbrRow)
    {
      if (this->NbrRow < nbrRow)
	{
	  int Tot = (nbrRow * (nbrRow - 1));
	  for (int i = this->NbrRow; i < nbrRow; i++)
	    this->DiagonalElements[i] = 0.0;
	  int k = (this->NbrRow - 1);
	  for (int i = 0; i < (this->NbrRow - 1); i++)
	    {
	      for (int j = this->NbrRow; j < nbrRow; j++)
		{
		  this->OffDiagonalElements[k++] = 0.0;
		}
	      k += (this->NbrRow - 2 - i);
	    }
	  for (int i = this->NbrRow * (this->NbrRow - 1); i < Tot; i++)
	    this->OffDiagonalElements[i] = 0.0;
	}
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      this->Increment = (this->TrueNbrRow - this->NbrRow);
      return;
    }
  double* TmpDiag = new double [nbrRow];
  int Tot = (nbrRow * (nbrRow - 1)) / 2;
  double* TmpOffDiag = new double [Tot];
  for (int i = 0; i < this->NbrRow; i++)
    TmpDiag [i] = this->DiagonalElements[i];
  for (int i = this->NbrRow; i < nbrRow; i++)
    TmpDiag [i]  = 0.0;
  int k = 0;
  int l = 0;
  for (int i = 0; i < (this->NbrRow - 1); i++)
    {
      for (int j = i + 1; j < this->NbrRow; j++)
	TmpOffDiag[k++] = this->OffDiagonalElements[l++];
      l += this->Increment;
      for (int j = this->NbrRow; j < nbrRow; j++)
	{
	  TmpOffDiag[k++] = 0.0;
	}      
    }
  for (int i = (this->NbrRow * (this->NbrRow - 1)) / 2; i < Tot; i++)
    TmpOffDiag[i] = 0.0;
  if (this->OffDiagonalGarbageFlag != 0)
    {
      if ((*(this->OffDiagonalGarbageFlag)) == 1)
	{
	  delete this->OffDiagonalGarbageFlag;
	  delete[] this->OffDiagonalElements;
	}
      else
	(*(this->OffDiagonalGarbageFlag))--;
    }
  if (this->DiagonalGarbageFlag != 0)
    {
      if ((*(this->DiagonalGarbageFlag)) == 1)
	{
	  delete this->DiagonalGarbageFlag;
	  delete[] this->DiagonalElements;
	}
      else
	(*(this->DiagonalGarbageFlag))--;
    }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = 2 * (this->TrueNbrRow - this->NbrRow);
  this->DiagonalElements = TmpDiag;
  this->OffDiagonalElements = TmpOffDiag;
  this->DiagonalGarbageFlag =  new int;
  *(this->DiagonalGarbageFlag) = 1;
  this->OffDiagonalGarbageFlag =  new int;
  *(this->OffDiagonalGarbageFlag) = 1;
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

RealUpperTriangularMatrix operator + (const RealUpperTriangularMatrix& M1, const RealUpperTriangularMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return RealUpperTriangularMatrix();
  double* Diagonal = new double [M1.NbrRow];
  int ReducedNbr = M1.NbrRow - 1;
  double* OffDiagonal = new double [M1.NbrRow * ReducedNbr];
  for (int i = 0; i < M1.NbrRow; i++)
    {
      Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
    }
  int k = 0;
  int l1 = 0;
  int l2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = 0; j < i; j++)
	OffDiagonal[k++] = M1.OffDiagonalElements[l1++] + M2.OffDiagonalElements[l2++];      
      l1 += M2.Increment;
      l2 += M2.Increment;
    }
  return RealUpperTriangularMatrix(Diagonal, OffDiagonal, M1.NbrRow);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

RealUpperTriangularMatrix operator - (const RealUpperTriangularMatrix& M1, const RealUpperTriangularMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return RealUpperTriangularMatrix();
  double* Diagonal = new double [M1.NbrRow];
  int ReducedNbr = M1.NbrRow - 1;
  double* OffDiagonal = new double [M1.NbrRow * ReducedNbr];
  for (int i = 0; i < M1.NbrRow; i++)
    {
      Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
    }
  int k = 0;
  int l1 = 0;
  int l2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = 0; j < i; j++)
	OffDiagonal[k++] = M1.OffDiagonalElements[l1++] - M2.OffDiagonalElements[l2++];      
      l1 += M2.Increment;
      l2 += M2.Increment;
    }
  return RealUpperTriangularMatrix(Diagonal, OffDiagonal, M1.NbrRow);
}

// multiply a real matrix with a real upper triangular matrix
//
// m1 = real matrix
// m2 = real upper triangular matrix
// return value = product result

RealMatrix operator * (RealMatrix& m1, const RealUpperTriangularMatrix& m2)
{
  RealMatrix TmpM(m1.GetNbrRow(), m2.NbrRow);
  double Tmp = 0.0;
  int Pos = 0;
  for (int i = 0; i < m1.GetNbrRow(); ++i)
    {
      for (int j = 0; j < m2.NbrRow; ++j)
	{
	  Tmp = m1(i, j) * m2.DiagonalElements[j];
	  Pos = ((j - 1) * j) / 2;
	  for (int k = 0; k < j; ++k)
	    {
	      Tmp += m1(i, k) * m2.OffDiagonalElements[Pos];
	      ++Pos;
	    }
	  TmpM(i, j) = Tmp;	
	}
    }
  return TmpM;
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealUpperTriangularMatrix operator * (const RealUpperTriangularMatrix& M, double x) 
{
  double* Diagonal = new double [M.NbrRow];
  int ReducedNbr = M.NbrRow - 1;
  double* OffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  for (int i = 0; i < M.NbrRow; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
    }
  int k = 0;
  int k2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	OffDiagonal[k++] = M.OffDiagonalElements[k2++] * x;
      k2 += M.Increment;
    }
  return RealUpperTriangularMatrix(Diagonal, OffDiagonal, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealUpperTriangularMatrix operator * (double x, const RealUpperTriangularMatrix& M) 
{
  double* Diagonal = new double [M.NbrRow];
  int ReducedNbr = M.NbrRow - 1;
  double* OffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  for (int i = 0; i < M.NbrRow; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
    }
  int k = 0;
  int k2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	OffDiagonal[k++] = M.OffDiagonalElements[k2++] * x;
      k2 += M.Increment;
    }
  return RealUpperTriangularMatrix(Diagonal, OffDiagonal, M.NbrRow);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

RealUpperTriangularMatrix operator / (const RealUpperTriangularMatrix& M, double x) 
{
  x = 1.0 / x;
  double* Diagonal = new double [M.NbrRow];
  int ReducedNbr = M.NbrRow - 1;
  double* OffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  for (int i = 0; i < M.NbrRow; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
    }
  int k = 0;
  int k2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	OffDiagonal[k++] = M.OffDiagonalElements[k2++] * x;
      k2 += M.Increment;
    }
  return RealUpperTriangularMatrix(Diagonal, OffDiagonal, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

RealUpperTriangularMatrix& RealUpperTriangularMatrix::operator += (const RealUpperTriangularMatrix& M) 
{
  if (this->NbrRow == 0)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  for (int i = 0; i < M.NbrRow; i++)
    {
      this->DiagonalElements[i] += M.DiagonalElements[i];
    }
  int k = 0;
  int k2 = 0;  
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	this->OffDiagonalElements[k++] += M.OffDiagonalElements[k2++];
      k += this->Increment;
      k2 += M.Increment;
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

RealUpperTriangularMatrix& RealUpperTriangularMatrix::operator -= (const RealUpperTriangularMatrix& M) 
{
  if (this->NbrRow == 0)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  for (int i = 0; i < M.NbrRow; i++)
    {
      this->DiagonalElements[i] -= M.DiagonalElements[i];
    }
  int k = 0;
  int k2 = 0;  
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	this->OffDiagonalElements[k++] -= M.OffDiagonalElements[k2++];
      k += this->Increment;
      k2 += M.Increment;
    }
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

RealUpperTriangularMatrix& RealUpperTriangularMatrix::operator *= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  int ReducedNbr = this->NbrRow - 1;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] *= x;
    }
  int k = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	this->OffDiagonalElements[k++] *= x;
      k += this->Increment;
    }
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

RealUpperTriangularMatrix& RealUpperTriangularMatrix::operator /= (double x)
{
  if (this->NbrRow == 0)
    return *this;
  x = 1.0 / x;
  int ReducedNbr = this->NbrRow - 1;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] *= x;
    }
  int k = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	this->OffDiagonalElements[k++] *= x;
      k += this->Increment;
    }
  return *this;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

double RealUpperTriangularMatrix::MatrixElement (RealVector& V1, RealVector& V2)
{
  double x = 0.0;
  if ((V1.Dimension != this->NbrRow) || (V2.Dimension != this->NbrColumn))
    return x;
  for (int i = 0; i < this->NbrRow ; i++)
    {
      double x2 = this->DiagonalElements[i] * V2.Components[i];
      int l = (i - 1);
      for (int k = 0; k < i; k++)
	{
	  x2 += this->OffDiagonalElements[l] * V2.Components[k];
	  l += (this->NbrColumn - 2 - k) + this->Increment;
	}
      x += x2 * V1.Components[i];
    }
  return x;
}

// evaluate matrix trace
//
// return value = matrix trace 

double RealUpperTriangularMatrix::Tr () 
{
  if (this->NbrRow == 0)
    return 0.0;
  double x = this->DiagonalElements[0];
  for (int i = 1; i < this->NbrRow; i++)
    {
      x += this->DiagonalElements[i];
    }
  return x;
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double RealUpperTriangularMatrix::Det () 
{
  if (this->NbrRow == 0)
    return 0.0;
  double x = this->DiagonalElements[0];
  for (int i = 1; i < this->NbrRow; i++)
    {
      x *= this->DiagonalElements[i];
    }
  return x;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const RealUpperTriangularMatrix& P)
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      for (int j = 0; j < i; j ++)
	{
	  Str << 0.0 << "    ";
	}
      Str << P.DiagonalElements[i] << "    ";
      for (int j = i + 1; j < P.NbrRow; j++)
	{
	  Str << P.OffDiagonalElements[i + ((j - 1) * j) / 2] << "    ";
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

MathematicaOutput& operator << (MathematicaOutput& Str, const RealUpperTriangularMatrix& P)
{
  Str << "{";
  for (int i = 0; i < P.NbrRow; i++)
    {
      Str << "{";
      int pos = (i - 1);
      for (int j = 0; j < i; j ++)
	{
	  Str << P.OffDiagonalElements[pos] << ",";
	  pos += (P.NbrRow - j - 2) + P.Increment;
	}
      Str << P.DiagonalElements[i];
      if (i != (P.NbrRow - 1))
	{
	  Str << ",";	  
	  pos++;
	  for (int j = i + 1; j < (P.NbrRow - 1); j++)
	    {
	      Str << 0.0 << ",";
	    }
	  Str << 0.0;
	  Str << "},";
	}
      else
	Str << "}";
    }
  Str << "}";
  return Str;
}

#endif

// output file stream overload
//
// file = reference on output file stream
// matrix = reference on matrix to save
// return value = reference on output file stream

/*ofstream& operator << (ofstream& file, const RealUpperTriangularMatrix& matrix)
{
  file.write ((char*) &(matrix.NbrRow), sizeof(int));
  file.write ((char*) &(matrix.NbrColumn), sizeof(int));
  file.write ((char*) (matrix.DiagonalElements), sizeof(double) * matrix.NbrRow);
  int Pos = 0;
  int ReducedNbrRow = matrix.NbrColumn - 1;
  for (int i = 0; i < ReducedNbrRow; i++)
    {
      file.write ((char*) &(matrix.OffDiagonalElements[Pos]), sizeof(double) * (ReducedNbrRow - i));
      Pos += ReducedNbrRow - i + matrix.Increment;
    }
  return file;
}*/

// input file stream overload
//
// file = reference on output file stream
// matrix = reference on matrix to load
// return value = reference on output file stream

ifstream& operator >> (ifstream& file, RealUpperTriangularMatrix& matrix)
{
  return file;
}

