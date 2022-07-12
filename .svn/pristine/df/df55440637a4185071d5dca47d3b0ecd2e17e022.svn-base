////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of real symmetric matrix                      //
//                                                                            //
//                        last modification : 09/03/2001                      //
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


#include "config.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "GeneralTools/ListIterator.h"
#ifdef USE_HILBERT_SPACE
#include "HilbertSpace/SubspaceSpaceConverter.h"
#endif

#include <stdlib.h>
#include <math.h>
#include <fstream>


#ifdef HAVE_LAPACK

// binding to the LAPACK dsyev function
//
extern "C" void FORTRAN_NAME(dsyev)(const char* jobz, const char* uplo, const int* dimension, const double* matrix, const int* leadingDimension,
				    const double* eigenvalues, const double* workingArea, const int* workingAreaSize, const int* information);

// binding to the LAPACK dsyevd function
//
extern "C" void FORTRAN_NAME(dsyevd)(const char* jobz, const char* uplo, const int* dimension, const double* matrix, const int* leadingDimension,
                                     const double* eigenvalues, const double* workingArea, const int* workingAreaSize, 
				     const int* integerWorkingArea, const int* integerWorkingAreaSize, const int* information);

#endif



using std::endl;


// default constructor
//

RealSymmetricMatrix::RealSymmetricMatrix() 
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
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

RealSymmetricMatrix::RealSymmetricMatrix(int dimension, bool zero) 
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
}

// constructor from matrix elements (without duplicating datas)
//
// diagonal = pointer to diagonal element array
// offDiagonal = pointer to off-diagonal element array (with real part in even position and imaginary part in odd position)
// dimension = matrix dimension

RealSymmetricMatrix::RealSymmetricMatrix(double* diagonal, double* offDiagonal, int dimension) 
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
}

// constructor from a real matrix Q (new matrix = Qt * Q)
//

RealSymmetricMatrix::RealSymmetricMatrix(const RealMatrix& Q)
{
  this->NbrRow = Q.NbrColumn;
  this->NbrColumn = Q.NbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::RealElements | Matrix::Symmetric;
  this->DiagonalGarbageFlag =  new int;
  *(this->DiagonalGarbageFlag) = 1;
  this->OffDiagonalGarbageFlag =  new int;
  *(this->OffDiagonalGarbageFlag) = 1;
  this->DiagonalElements = new double [this->NbrRow];
  this->OffDiagonalElements = new double [(this->NbrRow * (this->NbrRow - 1)) / 2];
  int pos = 0;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] = Q.Columns[i][0] * Q.Columns[0][i];
      for (int k = 1; k < Q.NbrColumn; k++)
	this->DiagonalElements[i] += Q.Columns[i][k] * Q.Columns[k][i];
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  this->OffDiagonalElements[pos] = Q.Columns[j][0] * Q.Columns[0][i];
	  for (int k = 1; k < Q.NbrColumn; k++)
	    this->OffDiagonalElements[pos] += Q.Columns[j][k] * Q.Columns[k][i];
	  pos++;
	}
    }
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

RealSymmetricMatrix::RealSymmetricMatrix(const RealSymmetricMatrix& M) 
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
}

// copy constructor from any matrix (only keeping real part of elements of and above the diagonal, duplicating datas)
//
// M = matrix to copy

RealSymmetricMatrix::RealSymmetricMatrix(Matrix& M)
{
  int Min = M.GetNbrRow();
  if (Min > M.GetNbrColumn())
    {
      Min = M.GetNbrColumn();
    }
  this->NbrRow = Min;
  this->NbrColumn = Min;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::RealElements | Matrix::Symmetric;
  this->DiagonalGarbageFlag =  new int;
  *(this->DiagonalGarbageFlag) = 1;
  this->OffDiagonalGarbageFlag =  new int;
  *(this->OffDiagonalGarbageFlag) = 1;
  this->DiagonalElements = new double [this->NbrRow];
  this->OffDiagonalElements = new double [(this->NbrRow * (this->NbrRow - 1)) / 2];
  int pos = 0;
  double Tmp;
  for (int i = 0; i < this->NbrRow; i++)
    {
      M.GetMatrixElement(i, i, Tmp);
      this->DiagonalElements[i] = Tmp;
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  M.GetMatrixElement(i, j, Tmp);	  
	  this->OffDiagonalElements[pos] = Tmp;
	  pos++;
	}
    }
}



// destructor
//

RealSymmetricMatrix::~RealSymmetricMatrix() 
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

RealSymmetricMatrix& RealSymmetricMatrix::operator = (const RealSymmetricMatrix& M) 
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
  return *this;
}

// assignement from  a real tridiagonal symmetric matrix (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix  

RealSymmetricMatrix& RealSymmetricMatrix::operator = (const RealTriDiagonalSymmetricMatrix& M) 
{
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* RealSymmetricMatrix::Clone ()
{
  return ((Matrix*) new RealSymmetricMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void RealSymmetricMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      this->DiagonalElements[i] = x;
    }
  else
    {
      if (i > j)
	{
	  int tmp = j;
	  j = i;
	  i = tmp;
	}
      j -= (i * (i + 1)) / 2 - i * (this->NbrRow + this->Increment - 1) + 1;
      this->OffDiagonalElements[j] = x;
    }
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element
void RealSymmetricMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  return;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void RealSymmetricMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      this->DiagonalElements[i] += x;
    }
  else
    {
      if (i > j)
	{
	  int tmp = j;
	  j = i;
	  i = tmp;
	}
      j -= (i * (i + 1)) / 2 - i * (this->NbrRow + this->Increment - 1) + 1;
      this->OffDiagonalElements[j] += x;
    }
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void RealSymmetricMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  return;
}

// get reference of a given matrix element
//
// i = line position
// j = column position
// return value = reference om matrix elememt

double& RealSymmetricMatrix::operator () (int i, int j)
{
  if (i == j)
    {
      return this->DiagonalElements[i];
    }
  else
    {
      if (i > j)
	{
	  int tmp = j;
	  j = i;
	  i = tmp;
	}
      j -= (i * (i + 1)) / 2 - i * (this->NbrRow + this->Increment - 1) + 1;
      return this->OffDiagonalElements[j];
    }
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void RealSymmetricMatrix::Resize (int nbrRow, int nbrColumn)
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

void RealSymmetricMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
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

// copy matrix
//
// M = matrix to copy
// return value = refence on current matrix

RealSymmetricMatrix& RealSymmetricMatrix::Copy (RealSymmetricMatrix& M)
{
  if (this->NbrRow != M.NbrRow)
    this->Resize(M.NbrRow, M.NbrColumn);
  int Pos1 = 0;
  int Pos2 = 0;
  for (int i = 0; i < M.NbrColumn; i++)
    {
      for (int j = i + 1; j < M.NbrColumn; ++j)
	{
	  this->OffDiagonalElements[Pos1] = M.OffDiagonalElements[Pos2];
	  ++Pos1;
	  ++Pos2;
	}
      Pos1 += this->Increment;
      Pos2 += M.Increment;
      this->DiagonalElements[i] = M.DiagonalElements[i];
    }
  return *this;
}

#ifdef USE_HILBERT_SPACE
// project matrix into a given subspace
//
// subspace = reference on subspace structure
// return value = pointer to projected matrix

Matrix* RealSymmetricMatrix::Project (SubspaceSpaceConverter& subspace)
{
  RealSymmetricMatrix* TmpM = new RealSymmetricMatrix (subspace.SubspaceDimension);
  int Pos;
  int Pos2;
  int Pos3 = 0;
  for (int i = 0; i < subspace.SubspaceDimension; i++)
    {
      Pos = subspace.SubspaceSpaceConverterArray[i];
      TmpM->DiagonalElements[i] = this->DiagonalElements[Pos];
      for (int j = i + 1; j < subspace.SubspaceDimension; j++)
	{
	  Pos2 = subspace.SubspaceSpaceConverterArray[j];
	  if (Pos2 > Pos)
	    TmpM->OffDiagonalElements[Pos3++] = this->OffDiagonalElements[Pos2 - 1 - (Pos * (Pos + 3)) / 2 
									 + Pos * this->NbrRow];	
	  else
	    TmpM->OffDiagonalElements[Pos3++] = this->OffDiagonalElements[Pos - 1 - (Pos2 * (Pos2 + 3)) / 2
									 + Pos2 * this->NbrRow];	
	}
    }
  return TmpM;
}
#endif

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

RealSymmetricMatrix operator + (const RealSymmetricMatrix& M1, const RealSymmetricMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return RealSymmetricMatrix();
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
  return RealSymmetricMatrix(Diagonal, OffDiagonal, M1.NbrRow);
}

// add two matrices where the left one is a real tridiagonal symmetric matrix
//
// M1 = left matrix
// M2 = right matrix
// return value = sum of the two matrices

RealSymmetricMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const RealSymmetricMatrix& M2)
{
  return RealSymmetricMatrix();
}

// add two matrices where the right one is a real tridiagonal symmetric matrix
//
// M1 = left matrix
// M2 = right matrix
// return value = sum of the two matrices

RealSymmetricMatrix operator + (const RealSymmetricMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2)
{
  return RealSymmetricMatrix();
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

RealSymmetricMatrix operator - (const RealSymmetricMatrix& M1, const RealSymmetricMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return RealSymmetricMatrix();
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
  return RealSymmetricMatrix(Diagonal, OffDiagonal, M1.NbrRow);
}

// substract two matrices where the left one is a real tridiagonal symmetric matrix
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

RealSymmetricMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, const RealSymmetricMatrix& M2)
{
  return RealSymmetricMatrix();
}

// substract two matrices where the right one is a real tridiagonal symmetric matrix
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

RealSymmetricMatrix operator - (const RealSymmetricMatrix& M1,  const RealTriDiagonalSymmetricMatrix& M2)
{
  return RealSymmetricMatrix();
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealSymmetricMatrix operator * (const RealSymmetricMatrix& M, double x) 
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
  return RealSymmetricMatrix(Diagonal, OffDiagonal, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealSymmetricMatrix operator * (double x, const RealSymmetricMatrix& M) 
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
  return RealSymmetricMatrix(Diagonal, OffDiagonal, M.NbrRow);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

RealSymmetricMatrix operator / (const RealSymmetricMatrix& M, double x) 
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
  return RealSymmetricMatrix(Diagonal, OffDiagonal, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

RealSymmetricMatrix& RealSymmetricMatrix::operator += (const RealSymmetricMatrix& M) 
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

RealSymmetricMatrix& RealSymmetricMatrix::operator += (const RealTriDiagonalSymmetricMatrix& M) 
{
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

RealSymmetricMatrix& RealSymmetricMatrix::operator -= (const RealSymmetricMatrix& M) 
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

// substract two matrices where the right one is a real tridiagonal symmetric matrix
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

RealSymmetricMatrix& RealSymmetricMatrix::operator -= (const RealTriDiagonalSymmetricMatrix& M) 
{
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

RealSymmetricMatrix& RealSymmetricMatrix::operator *= (double x) 
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

RealSymmetricMatrix& RealSymmetricMatrix::operator /= (double x)
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

double RealSymmetricMatrix::MatrixElement (RealVector& V1, RealVector& V2)
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
      l++;
      for (int k = i + 1; k < this->NbrColumn; k++)
	{
	  x2 += this->OffDiagonalElements[l++] * V2.Components[k];
	}      
      x += x2 * V1.Components[i];
    }
  return x;
}

// conjugate matrix with an unitary matrix (Ut M U)
//
// UnitaryM = unitary matrix to use
// return value = pointer to conjugated matrix

Matrix* RealSymmetricMatrix::Conjugate(RealMatrix& UnitaryM)
{
  if (UnitaryM.NbrRow != this->NbrColumn)
    return 0;
  double* TmpDiag = new double [UnitaryM.NbrColumn];
  int NbrOffDiag = (UnitaryM.NbrColumn * (UnitaryM.NbrColumn - 1)) / 2;
  double* TmpOffDiag = new double [NbrOffDiag];
  int Inc = this->NbrColumn - 2 + this->Increment;
  int k;
  int l;
  double tmp;
  for (int i = 0; i < UnitaryM.NbrColumn; i++)
    {
      TmpDiag[i] = 0.0;
      for (int j = 0; j < this->NbrColumn; j++)
	{
	  tmp = 0.0;
	  k = 0;
	  l = (j - 1);
	  for (; k < j; k++)
	    {
	      tmp += this->OffDiagonalElements[l] * UnitaryM.Columns[i].Components[k];
	      l += Inc - k;
	    }
	  l++;
	  tmp += this->DiagonalElements[j] * UnitaryM.Columns[i].Components[j];
	  k++;
	  for (; k < this->NbrColumn; k++)
	    {
	      tmp += this->OffDiagonalElements[l] * UnitaryM.Columns[i].Components[k];
	      l++;
	    }
	  TmpDiag[i] += tmp * UnitaryM.Columns[i].Components[j];
	}
    }    
  int i2 = 0;
  Inc--;
  int ReducedNbrColumn = UnitaryM.NbrColumn - 1;
  for (int i = 0; i < ReducedNbrColumn; i++)
    {
      for (int m = i + 1; m < UnitaryM.NbrColumn; m++)
	{    
	  TmpOffDiag[i2] = 0.0;
	  for (int j = 0; j < this->NbrColumn; j++)
	    {
	      tmp = 0.0;
	      k = 0;
	      l = (j - 1);
	      for (; k < j; k++)
		{
		  tmp += this->OffDiagonalElements[l++] * UnitaryM.Columns[m].Components[k];
		  l += Inc  - k;
		}
	      l++;
	      tmp += this->DiagonalElements[j] * UnitaryM.Columns[m].Components[j];
	      k++;
	      for (; k < this->NbrColumn; k++)
		{
		  tmp += this->OffDiagonalElements[l++] * UnitaryM.Columns[m].Components[k];
		}
	      TmpOffDiag[i2] += tmp * UnitaryM.Columns[i].Components[j];
	    }
	  i2++;
	}    
    }
  return new RealSymmetricMatrix(TmpDiag, TmpOffDiag, UnitaryM.NbrColumn);
}

// conjugate an hermitian matrix with an unitary matrix (Ut M U)
//
// UnitaryM = unitary matrix to use
// return value = pointer to conjugated matrix

Matrix* RealSymmetricMatrix::Conjugate(BlockDiagonalMatrix& UnitaryM)
{
  if (UnitaryM.NbrRow != this->NbrColumn)
    return 0;
  RealSymmetricMatrix* Result = new RealSymmetricMatrix(UnitaryM.NbrColumn);
  int NbrBlock = UnitaryM.Blocks.GetNbrElement();
  Matrix** Blocks = new Matrix* [NbrBlock];
  Matrix** TmpM;
  ListIterator<Matrix*> IterBlocks (UnitaryM.Blocks);
  int i = 0;
  while ((TmpM = IterBlocks()))
    {
      Blocks[i++] = *TmpM;
    }
  int SrcRowIndex = 0;
  int DestRowIndex = 0;
  int SrcColumnIndex = 0;
  int DestColumnIndex = 0;
  for (int i = 0; i < NbrBlock; i++)
    {
//      SrcColumnIndex = SrcRowIndex;
//      DestColumnIndex = DestRowIndex;
      SrcRowIndex = UnitaryM.BlockRowPosition[i];
      DestRowIndex = UnitaryM.BlockColumnPosition[i];
      SrcColumnIndex =  UnitaryM.BlockRowPosition[i];
      DestColumnIndex = UnitaryM.BlockColumnPosition[i];
      if (Blocks[i]->GetMatrixType() == Matrix::RealElements)
	{
	  this->Conjugate(*((RealMatrix*) Blocks[i]), SrcRowIndex, DestRowIndex, *Result);
	}
//      SrcColumnIndex += Blocks[i]->GetNbrRow();
//      DestColumnIndex += Blocks[i]->GetNbrColumn();
      for (int j = i + 1; j < NbrBlock; j++)
	{
	  SrcColumnIndex = UnitaryM.BlockRowPosition[j];
	  DestColumnIndex = UnitaryM.BlockColumnPosition[j];
	  if ((Blocks[i]->GetMatrixType() == Matrix::RealElements) && 
	      (Blocks[j]->GetMatrixType() == Matrix::RealElements))
	    {
	      this->Conjugate(*((RealMatrix*) Blocks[i]), *((RealMatrix*) Blocks[j]), 
			      SrcRowIndex, SrcColumnIndex, DestRowIndex, DestColumnIndex, 
			      *Result);
	    }
//	  SrcColumnIndex += Blocks[j]->GetNbrRow();
//	  DestColumnIndex += Blocks[j]->GetNbrColumn();
	}
      SrcRowIndex += Blocks[i]->GetNbrRow();
      DestRowIndex += Blocks[i]->GetNbrColumn();
    }
  delete[] Blocks;
  return Result;
}

// conjugate a block of the matrix with an unitary matrix (Ut M U)
//
// UnitaryM = unitary matrix to use
// sourcePosition = index of the row where the block to conjugate starts
// destinationPosition = index of the row where the conjugated block has to be stored
// matrix = matrix where result has to be stored

void RealSymmetricMatrix::Conjugate(RealMatrix& UnitaryM, int sourcePosition, int destinationPosition,
				    RealSymmetricMatrix& matrix)
{
  if (((UnitaryM.NbrRow + sourcePosition) > this->NbrColumn) || 
      ((UnitaryM.NbrColumn + destinationPosition) > matrix.NbrColumn))
    return;
  for (int i = 0; i < UnitaryM.NbrColumn; i++)
    {
      matrix.DiagonalElements[i + destinationPosition] = 0.0;
      for (int j = 0; j < UnitaryM.NbrRow; j++)
	{
	  double tmp = 0.0;
	  int k = 0;
	  int l = ((j + sourcePosition) - 1 - (sourcePosition * (sourcePosition + 1)) / 2 +
		   sourcePosition * (this->NbrRow + this->Increment - 1));
	  for (; k < j; k++)
	    {
	      tmp += this->OffDiagonalElements[l] * UnitaryM.Columns[i].Components[k];
	      l += (this->NbrColumn - 2 - k - sourcePosition) + this->Increment;
	    }
	  l++;
	  tmp += this->DiagonalElements[j + sourcePosition] * UnitaryM.Columns[i].Components[j];
	  k++;
	  for (; k < UnitaryM.NbrRow; k++)
	    {
	      tmp += this->OffDiagonalElements[l] * UnitaryM.Columns[i].Components[k];
	      l++;
	    }
	  matrix.DiagonalElements[i + destinationPosition] += tmp * UnitaryM.Columns[i].Components[j];
	}
    }    
  int i2 = (destinationPosition - (destinationPosition * (destinationPosition + 1)) / 2 +
	    destinationPosition * (matrix.NbrRow + matrix.Increment - 1));
  int ReducedNbrColumn = UnitaryM.NbrColumn - 1;
  for (int i = 0; i < ReducedNbrColumn; i++)
    {
      for (int m = i + 1; m < UnitaryM.NbrColumn; m++)
	{    
	  matrix.OffDiagonalElements[i2] = 0.0;
	  for (int j = 0; j < UnitaryM.NbrRow; j++)
	    {
	      double tmp1 = 0.0;
	      int k = 0;
	      int l = ((j + sourcePosition) - 1 - (sourcePosition * (sourcePosition + 1)) / 2 +
		       sourcePosition * (this->NbrRow + this->Increment - 1));
	      for (; k < j; k++)
		{
		  tmp1 += this->OffDiagonalElements[l++] * UnitaryM.Columns[m].Components[k];
		  l += (this->NbrColumn - 3 - k - sourcePosition) + this->Increment;
		}
	      l++;
	      tmp1 += this->DiagonalElements[j + sourcePosition] * UnitaryM.Columns[m].Components[j];
	      k++;
	      for (; k < UnitaryM.NbrRow; k++)
		{
		  tmp1 += this->OffDiagonalElements[l++] * UnitaryM.Columns[m].Components[k];
		}
	      matrix.OffDiagonalElements[i2] += tmp1 * UnitaryM.Columns[i].Components[j];
	    }
	  i2++;
	}
      i2 += matrix.NbrColumn - destinationPosition - UnitaryM.NbrColumn + matrix.Increment;
    }
  return;
}

// conjugate a block of the matrix (in the upper diagonal part) with two matrix matrix (Vt M U)
//
// UnitaryMl = unitary matrix to use at the left hand side
// UnitaryMr = unitary matrix to use at the right hand side
// sourceRowIndex = index of the row where the block to conjugate starts
// sourceColumnIndex = index of the column where the block to conjugate starts
// destinationRowIndex = index of the row where the conjugated block has to be stored
// destinationColumnIndex = index of the column where the conjugated block has to be stored
// matrix = matrix where result has to be stored

void RealSymmetricMatrix::Conjugate(RealMatrix& UnitaryMl, RealMatrix& UnitaryMr, int sourceRowIndex, 
				    int sourceColumnIndex, int destinationRowIndex,
				    int destinationColumnIndex, RealSymmetricMatrix& matrix)
{
  if (((UnitaryMr.NbrRow + sourceColumnIndex) > this->NbrColumn) || 
      ((UnitaryMl.NbrRow + sourceRowIndex) > this->NbrRow) || 
      ((UnitaryMr.NbrColumn + destinationColumnIndex) > matrix.NbrColumn) || 
      ((UnitaryMl.NbrColumn + destinationRowIndex) > matrix.NbrRow))
    return;
  double tmp1;
  int l;
  int i2;
  for (int i = 0; i < UnitaryMl.NbrColumn; i++)
    {
      i2 = (destinationColumnIndex - 1 - ((i + destinationRowIndex) * ((i + destinationRowIndex) + 1)) / 2 +
	    (i + destinationRowIndex) * (matrix.NbrRow + matrix.Increment - 1));
      RealVector& CurrentColumn1 = UnitaryMl.Columns[i];
      for (int m = 0; m < UnitaryMr.NbrColumn; ++m)
	{    
	  RealVector& CurrentColumn2 = UnitaryMr.Columns[m];
	  matrix.OffDiagonalElements[i2] = 0.0;
	  for (int j = 0; j < UnitaryMl.NbrRow; ++j)
	    {
	      tmp1 = 0.0;
	      l = (sourceColumnIndex - 1 - 
		   ((j + sourceRowIndex) * ((sourceRowIndex + j) + 1)) / 2 +
		   (sourceRowIndex + j) * (this->NbrRow + this->Increment - 1));
	      for (int k = 0; k < UnitaryMr.NbrRow; ++k)
		{
		  tmp1 += this->OffDiagonalElements[l] * CurrentColumn2.Components[k];
		  ++l;
		}
	      matrix.OffDiagonalElements[i2] += tmp1 * CurrentColumn1.Components[j];
	    }
//	  matrix.OffDiagonalElements[i2] = tmp2;
	  ++i2;
	}
    }
  return;
}

// evaluate matrix trace
//
// return value = matrix trace 

double RealSymmetricMatrix::Tr () 
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

double RealSymmetricMatrix::Det () 
{
  return 1.0;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const RealSymmetricMatrix& P)
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      int pos = (i - 1);
      for (int j = 0; j < i; j ++)
	{
	  Str << P.OffDiagonalElements[pos] << "    ";
	  pos += (P.NbrRow - j - 2) + P.Increment;
	}
      Str << P.DiagonalElements[i] << "    ";
      pos++;
      for (int j = i + 1; j < P.NbrRow; j++)
	{
	  Str << P.OffDiagonalElements[pos++] << "    ";
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

MathematicaOutput& operator << (MathematicaOutput& Str, const RealSymmetricMatrix& P)
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
	      Str << P.OffDiagonalElements[pos++] << ",";
	    }
	  Str << P.OffDiagonalElements[pos];
	  Str << "},";
	}
      else
	Str << "}";
    }
  Str << "}";
  return Str;
}

#endif

// Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step
//
// dimension = maximum iteration number
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// V1 = reference on complex vector used as first vector (will contain last produced vector at the end)
// return value = reference on complex tridiagonal hermitian matrix

RealTriDiagonalSymmetricMatrix& RealSymmetricMatrix::Lanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, RealVector& V1)
{  
  int Index = 0;
  V1 /= V1.Norm();
  RealVector V2(this->NbrRow);
  RealVector V3(this->NbrRow);
  V2.Multiply(*this, V1);
  M.DiagonalElements[Index] = (V1 * V2);
  V2.AddLinearCombination(-M.DiagonalElements[Index], V1);
  V2 /= V2.Norm();
  dimension -= 2;
  for (int i = 0; i < dimension; i++)
    {
      V3.Multiply(*this, V2);
      M.UpperDiagonalElements[Index] = (V1 * V3);
      M.DiagonalElements[Index + 1] = (V2 * V3);
      V3.AddLinearCombination(-M.DiagonalElements[Index + 1], V2);
      V3.AddLinearCombination(-M.UpperDiagonalElements[Index], V1);
      V3 /= V3.Norm();
      Index++;
      RealVector TmpV = V1;
      V1 = V2;
      V2 = V3;
      V3 = TmpV;
    }  
  V3.Multiply(*this, V2);
  M.UpperDiagonalElements[Index] = (V1 * V3);
  M.DiagonalElements[Index + 1] = (V2 * V3);
  RealVector TmpV = V1;
  V1 = V2;
  V2 = TmpV;
  return M;
}

// Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step
//
// dimension = maximum iteration number
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// Q = matrix where new orthonormalized base will be stored (first column is used as first vector)
// return value = reference on complex tridiagonal hermitian matrix

RealTriDiagonalSymmetricMatrix& RealSymmetricMatrix::Lanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, RealMatrix& Q)
{  
  int Index = 0;
  if ((Q.NbrColumn != dimension) || (Q.NbrRow != this->NbrRow))
    Q.Resize(this->NbrRow, dimension);
  if (M.NbrRow != dimension)
    M.Resize(dimension, dimension);
  Q.Columns[0] /= Q.Columns[0].Norm(); 
  Q.Columns[1].Multiply(*this, Q.Columns[0]);
  M.DiagonalElements[Index] = (Q.Columns[0] * Q.Columns[1]);
  Q.Columns[1].AddLinearCombination(-M.DiagonalElements[Index], Q.Columns[0]);
  Q.Columns[1] /= Q.Columns[1].Norm(); 
  for (int i = 2; i < dimension; i++)
    {
      Q.Columns[i].Multiply(*this, Q.Columns[i - 1]);
      M.UpperDiagonalElements[Index] = (Q.Columns[i - 2] * Q.Columns[i]);
      M.DiagonalElements[Index + 1] = (Q.Columns[i - 1] * Q.Columns[i]);
      Q.Columns[i].AddLinearCombination(-M.DiagonalElements[Index + 1], Q.Columns[i - 1]);
      Q.Columns[i].AddLinearCombination(-M.UpperDiagonalElements[Index], Q.Columns[i - 2]);
      Q.Columns[i] /= Q.Columns[i].Norm();
      Index++;
    }  
  M.UpperDiagonalElements[Index] = this->MatrixElement(Q.Columns[dimension - 2], Q.Columns[dimension - 1]);
  M.DiagonalElements[Index + 1] = this->MatrixElement(Q.Columns[dimension - 1], Q.Columns[dimension - 1]);
  return M;
}

// Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step, if during process a 
// null vector appears then new random vector is evaluated
//
// dimension = maximum iteration number
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// Q = matrix where new orthonormalized base will be stored (first column is used as first vector)
// err = absolute error on vector norm
// return value = reference on complex tridiagonal hermitian matrix

RealTriDiagonalSymmetricMatrix& RealSymmetricMatrix::OrthoLanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, RealMatrix& Q,
								   double err)
{  
  int Index = 0;
  if ((Q.NbrColumn != dimension) || (Q.NbrRow != this->NbrRow))
    Q.Resize(this->NbrRow, dimension);
  if (M.NbrRow != dimension)
    M.Resize(dimension, dimension);
  Q.Columns[0] /= Q.Columns[0].Norm(); 
  Q.Columns[1].Multiply(*this, Q.Columns[0]);
  M.DiagonalElements[Index] = (Q.Columns[0] * Q.Columns[1]);
  Q.Columns[1].AddLinearCombination(-M.DiagonalElements[Index], Q.Columns[0]);
  Q.Columns[1] /= Q.Columns[1].Norm(); 
  for (int i = 2; i < dimension; i++)
    {
      Q.Columns[i].Multiply(*this, Q.Columns[i - 1]);
      M.UpperDiagonalElements[Index] = (Q.Columns[i - 2] * Q.Columns[i]);
      M.DiagonalElements[Index + 1] = (Q.Columns[i - 1] * Q.Columns[i]);
      Q.Columns[i].AddLinearCombination(-M.DiagonalElements[Index + 1], Q.Columns[i - 1]);
      Q.Columns[i].AddLinearCombination(-M.UpperDiagonalElements[Index], Q.Columns[i - 2]);
      double VectorNorm = Q.Columns[i].Norm();
      while (VectorNorm < err)
	{
	  double tmp = 0;
	  for (int j = 0; j < Q.NbrRow; j++)
	    {
	      Q.Columns[i].Components[2 * j] = rand ();
	      tmp += Q.Columns[i].Components[2 * j] * Q.Columns[i].Components[2 * j];
	    }
	  tmp = sqrt(tmp);
	  Q.Columns[i] /= tmp;
	  Q.Columns[i] *= *this;
	  for (int j = 0; j < i; j++)
	    {
	      Q.Columns[i].AddLinearCombination(- (Q.Columns[j] * Q.Columns[i]), Q.Columns[j]);
	    }
	  VectorNorm = Q.Columns[i].Norm();
	}
      Q.Columns[i] /= VectorNorm;
      Index++;
    }  
  M.UpperDiagonalElements[Index] = this->MatrixElement(Q.Columns[dimension - 2], Q.Columns[dimension - 1]);
  M.DiagonalElements[Index + 1] = this->MatrixElement(Q.Columns[dimension - 1], Q.Columns[dimension - 1]);
  return M;
}

// Tridiagonalize a real symmetric matrix using Householder algorithm  (modifying current matrix)
//
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// err = absolute error on matrix element
// return value = reference on real tridiagonal symmetric matrix

RealTriDiagonalSymmetricMatrix& RealSymmetricMatrix::Householder (RealTriDiagonalSymmetricMatrix& M, double err)
{
  if (M.NbrRow != this->NbrRow)
    M.Resize(this->NbrRow, this->NbrColumn);
  int ReducedNbrRow = this->NbrRow - 1;
  int ReducedNbrRow2 ;
  double* TmpV = new double [ReducedNbrRow];
  double TmpNorm;
  double Coef;
  int Pos = 0;
  int TmpPos;
  int TmpPos2;
  int TmpPos3;
  int TmpPos4;
  int Inc;
  M.DiagonalElement(0) = this->DiagonalElements[0];
  double SquareErr = err * err;
  if (SquareErr < MACHINE_PRECISION)
    SquareErr = MACHINE_PRECISION;
  for (int i = 1; i < ReducedNbrRow; i++)
    {
      ReducedNbrRow2 = this->NbrRow - i;
      TmpNorm = 0.0;
      TmpPos = Pos;
      // construct vector for Householder transformation
      for (int j = 0; j < ReducedNbrRow2; j++)
	{
	  TmpNorm += this->OffDiagonalElements[TmpPos] * this->OffDiagonalElements[TmpPos];
	  TmpV[j] = 0.0;
	  TmpPos++;
	}
      TmpPos = Pos;
      Coef = sqrt(TmpNorm);
      M.UpperDiagonalElement(i- 1) = Coef;
      TmpNorm -= Coef * this->OffDiagonalElements[Pos];
      if  (TmpNorm > SquareErr)
	{
	  TmpNorm = 1.0 / sqrt(TmpNorm);
	  this->OffDiagonalElements[Pos] -= Coef;
	  for (int j = 0; j < ReducedNbrRow2; j++)
	    {
	      this->OffDiagonalElements[TmpPos++] *= TmpNorm;
	    }
	  
	  // store result of Hamiltonian applied to Householder vector
	  Coef = 0.0;
	  TmpPos2 = 0;
	  TmpPos4 = i * (this->NbrRow + this->Increment - 1) - (i * (i + 1)) / 2 - 1;
	  Inc = this->NbrRow - 2 + this->Increment;
	  for (int j = i; j < this->NbrRow; j++)
	    {
	      TmpPos = Pos;
	      TmpPos3 = j + TmpPos4;
	      int k = i;
	      for (; k < j; k++)
		{
		  TmpV[TmpPos2] +=  this->OffDiagonalElements[TmpPos++] * this->OffDiagonalElements[TmpPos3];
		  TmpPos3 += Inc - k;
		}
	      TmpV[TmpPos2] +=  this->DiagonalElements[k] * this->OffDiagonalElements[TmpPos++];
	      k++;
	      TmpPos3++;
	      for (; k < this->NbrRow; k++)
		{
		  TmpV[TmpPos2] +=  this->OffDiagonalElements[TmpPos++] * this->OffDiagonalElements[TmpPos3++];	      
		}
	      Coef += TmpV[TmpPos2] * this->OffDiagonalElements[Pos + TmpPos2];
	      TmpPos2++;
	    }
	  TmpPos = Pos;
	  Coef *= 0.5;
	  for (int j = 0; j < ReducedNbrRow2; j++)
	    {
	      TmpV[j] -=  this->OffDiagonalElements[TmpPos++] * Coef;
	    }
	  
	  // conjugate Hamiltonian
	  TmpPos2 = i * (this->NbrRow + this->Increment) - (i * (i + 1)) / 2;
	  TmpPos = Pos + 1;
	  for (int j = i; j < this->NbrRow; j++)
	    {	  
	      TmpPos4 = j - i;
	      this->DiagonalElements[j] -=  this->OffDiagonalElements[Pos] * 2.0 * TmpV[TmpPos4];
	      TmpPos3 = TmpPos;
	      for (int k = j + 1; k < this->NbrRow; k++)
		{
		  this->OffDiagonalElements[TmpPos2++] -=  (this->OffDiagonalElements[TmpPos3] * TmpV[TmpPos4] 
							    + this->OffDiagonalElements[Pos] *  TmpV[k - i]);
		  TmpPos3++;
		} 
	      TmpPos2 += this->Increment;
	      TmpPos++;
	      Pos++;
	    }
	}
      else
	{
	  Pos += (this->NbrRow - i);
	}
      M.DiagonalElement(i) = this->DiagonalElements[i];
      Pos += this->Increment;
    }
  M.UpperDiagonalElement(ReducedNbrRow - 1) = this->OffDiagonalElements[Pos - this->Increment];  
  M.DiagonalElement(ReducedNbrRow) = this->DiagonalElements[ReducedNbrRow];
  delete[] TmpV;
  return M;
}

// Tridiagonalize a real symmetric matrix using Householder algorithm and evaluate transformation matrix  (modifying current matrix)
//
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// err = absolute error on matrix element
// Q = matrix where transformation matrix has to be stored
// return value = reference on real tridiagonal symmetric matrix

RealTriDiagonalSymmetricMatrix& RealSymmetricMatrix::Householder (RealTriDiagonalSymmetricMatrix& M, double err, RealMatrix& Q)
{
  if (M.NbrRow != this->NbrRow)
    M.Resize(this->NbrRow, this->NbrColumn);
  if ((Q.NbrRow != M.NbrRow) || (Q.NbrColumn != M.NbrColumn))
    Q.Resize(this->NbrRow, this->NbrColumn);
   for (int i = 0; i < Q.NbrRow; i++)
     {
       for (int j = 0; j < i; j++)
	 Q.Columns[j].Components[i] = 0.0;       
       Q.Columns[i].Components[i] = 1.0;
       for (int j = i + 1; j < Q.NbrColumn; j++)
	 Q.Columns[j].Components[i] = 0.0;       
     }
  int ReducedNbrRow = this->NbrRow -1;
  int ReducedNbrRow2 ;
  double* TmpV = new double [ReducedNbrRow];
  double* TmpCoef = new double [ReducedNbrRow + 1];
  double TmpNorm;
  double Coef;
  int Pos = 0;
  int TmpPos;
  int TmpPos2;
  int TmpPos3;
  int TmpPos4;
  int TmpPos5;
  int TmpPos6;
  int Inc;
  M.DiagonalElement(0) = this->DiagonalElements[0];
  double SquareErr = err * err;
  if (SquareErr < MACHINE_PRECISION)
    SquareErr = MACHINE_PRECISION;
  for (int i = 1; i < ReducedNbrRow; i++)
    {
      ReducedNbrRow2 = this->NbrRow - i;
      TmpNorm = 0.0;
      TmpPos = Pos;
      // construct vector for Householder transformation
      for (int j = 0; j < ReducedNbrRow2; j++)
	{
	  TmpNorm += this->OffDiagonalElements[TmpPos] * this->OffDiagonalElements[TmpPos];
	  TmpV[j] = 0.0;
	  TmpPos++;
	}
      TmpPos = Pos;
      Coef = sqrt(TmpNorm);
      M.UpperDiagonalElement(i- 1) = Coef;
      TmpNorm -= Coef * this->OffDiagonalElements[Pos];
      if  (TmpNorm > SquareErr)
	{
	  TmpNorm = 1.0 / sqrt(TmpNorm);
	  this->OffDiagonalElements[Pos] -= Coef;
	  for (int j = 0; j < ReducedNbrRow2; j++)
	    {
	      this->OffDiagonalElements[TmpPos++] *= TmpNorm;
	    }
	  
	  // store result of Hamiltonian applied to Householder vector, evaluate coefficients for Q
	  Coef = 0.0;
	  TmpPos2 = 0;
	  TmpPos5 = 1;
	  TmpPos4 = i * (this->NbrRow + this->Increment - 1) - (i * (i + 1)) / 2 - 1;
	  Inc = this->NbrRow - 2 + this->Increment;
	  for (int j = 1; j < i; j++)
	    {
	      TmpPos = Pos;
	      TmpCoef[TmpPos5] = 0.0;
	      for (int k = i; k < this->NbrRow; k++)
		{
		  TmpCoef[TmpPos5] += Q.Columns[k].Components[j] * this->OffDiagonalElements[TmpPos++];
		}
	      TmpPos5++;
	    }
	  for (int j = i; j < this->NbrRow; j++)
	    {
	      TmpPos = Pos;
	      TmpPos3 = j + TmpPos4;
	      int k = i;
	      TmpCoef[TmpPos5] = 0.0;
	      for (; k < j; k++)
		{
		  TmpCoef[TmpPos5] += Q.Columns[k].Components[j] * this->OffDiagonalElements[TmpPos];
		  TmpV[TmpPos2] +=  this->OffDiagonalElements[TmpPos++] * this->OffDiagonalElements[TmpPos3];
		  TmpPos3 += Inc - k;
		}
	      TmpCoef[TmpPos5] += Q.Columns[j].Components[j] * this->OffDiagonalElements[TmpPos];
	      TmpV[TmpPos2] +=  this->DiagonalElements[k] * this->OffDiagonalElements[TmpPos++];
	      k++;
	      TmpPos3++;
	      for (; k < this->NbrRow; k++)
		{
		  TmpCoef[TmpPos5] += Q.Columns[k].Components[j] * this->OffDiagonalElements[TmpPos];
		  TmpV[TmpPos2] +=  this->OffDiagonalElements[TmpPos++] * this->OffDiagonalElements[TmpPos3++];	      
		}
	      Coef += TmpV[TmpPos2] * this->OffDiagonalElements[Pos + TmpPos2];
	      TmpPos2++;
	      TmpPos5++;
	    }
	  TmpPos = Pos;
	  Coef *= 0.5;
	  for (int j = 0; j < ReducedNbrRow2; j++)
	    {
	      TmpV[j] -=  this->OffDiagonalElements[TmpPos++] * Coef;
	    }
	  
	  // conjugate Hamiltonian
	  TmpPos2 = i * (this->NbrRow + this->Increment) - (i * (i + 1)) / 2;
	  TmpPos = Pos + 1;
	  TmpPos6 = Pos;
	  for (int j = 1; j < i; j++)
	    {
	      TmpPos5 = TmpPos6;
	      Coef = TmpCoef[j];
	      for (int k = i; k < this->NbrRow; k++)
		{
		  Q.Columns[k].Components[j] -=  this->OffDiagonalElements[TmpPos5++] * Coef;
		}	      
	    }
	  for (int j = i; j < this->NbrRow; j++)
	    {	  
	      TmpPos4 = j - i;
	      TmpPos5 = TmpPos6;
	      Coef = TmpCoef[j];
	      for (int k = i; k < j; k++)
		{
		  Q.Columns[k].Components[j] -=  this->OffDiagonalElements[TmpPos5++] * Coef;
		}	      
	      this->DiagonalElements[j] -=  this->OffDiagonalElements[Pos] * 2.0 * TmpV[TmpPos4];
	      Q.Columns[j].Components[j] -=  this->OffDiagonalElements[TmpPos5++] * Coef;
	      TmpPos3 = TmpPos;
	      for (int k = j + 1; k < this->NbrRow; k++)
		{
		  this->OffDiagonalElements[TmpPos2++] -=  (this->OffDiagonalElements[TmpPos3] * TmpV[TmpPos4] 
							    + this->OffDiagonalElements[Pos] *  TmpV[k - i]);
		  Q.Columns[k].Components[j] -=  this->OffDiagonalElements[TmpPos5++] * Coef;
		  TmpPos3++;
		} 
	      TmpPos2 += this->Increment;
	      TmpPos++;
	      Pos++;
	    }
	}
      else
	{
	  Pos += (this->NbrRow - i);
	}
      M.DiagonalElement(i) = this->DiagonalElements[i];
      Pos += this->Increment;
    }
  M.UpperDiagonalElement(ReducedNbrRow - 1) = this->OffDiagonalElements[Pos - this->Increment];  
  M.DiagonalElement(ReducedNbrRow) = this->DiagonalElements[ReducedNbrRow];
  delete[] TmpV;
  delete[] TmpCoef;
  return M;
}

// Diagonalize a real symmetric matrix (modifying current matrix)
//
// M = reference on real diagonal matrix where result has to be stored
// err = absolute error on matrix element
// maxIter = maximum number of iteration to fund an eigenvalue
// return value = reference on real tridiagonal symmetric matrix

RealDiagonalMatrix& RealSymmetricMatrix::Diagonalize (RealDiagonalMatrix& M, double err, int maxIter)
{
#ifdef __LAPACKONLY__
  return this->LapackDiagonalize(M, err, maxIter);
#endif
  if (M.GetNbrRow() != this->NbrRow)
    M.Resize(this->NbrRow, this->NbrColumn);
  RealTriDiagonalSymmetricMatrix TmpMatrix(this->NbrRow);
  this->Householder(TmpMatrix, err);
  TmpMatrix.Diagonalize(maxIter);
  for (int i = 0; i < this->NbrRow; ++i)
    M[i] = TmpMatrix.DiagonalElement(i);
  return M;
}

// Diagonalize a real symmetric matrix and evaluate transformation matrix (modifying current matrix)
//
// M = reference on real diagonal matrix where result has to be stored
// Q = matrix where transformation matrix has to be stored
// err = absolute error on matrix element
// maxIter = maximum number of iteration to fund an eigenvalue
// return value = reference on real tridiagonal symmetric matrix

RealDiagonalMatrix& RealSymmetricMatrix::Diagonalize (RealDiagonalMatrix& M, RealMatrix& Q, double err, int maxIter)
{
#ifdef __LAPACKONLY__
  return this->LapackDiagonalize(M, Q, err, maxIter);
#endif
  if (M.GetNbrRow() != this->NbrRow)
    M.Resize(this->NbrRow, this->NbrColumn);
  RealTriDiagonalSymmetricMatrix TmpMatrix(this->NbrRow);
  this->Householder(TmpMatrix, err, Q);
  TmpMatrix.Diagonalize(Q, maxIter);
  for (int i = 0; i < this->NbrRow; ++i)
    M[i] = TmpMatrix.DiagonalElement(i);
  return M;
}

#ifdef __LAPACK__

// Diagonalize a real symmetric matrix using the LAPACK library (modifying current matrix)
//
// M = reference on real diagonal matrix where result has to be stored
// err = absolute error on matrix element
// maxIter = maximum number of iteration to fund an eigenvalue
// return value = reference on real tridiagonal symmetric matrix

RealDiagonalMatrix& RealSymmetricMatrix::LapackDiagonalize (RealDiagonalMatrix& M, double err, int maxIter)
{
  if (M.GetNbrRow() != this->NbrRow)
    M.Resize(this->NbrRow, this->NbrColumn);
  double* TmpMatrix = new double [this->NbrRow * this->NbrRow];
  int Information = 0;
  int WorkingAreaSize = -1;
  double TmpWorkingArea;
  char Jobz = 'N';
  char UpperLower = 'L';
  int TotalIndex = 0;
  int Index2 = 0;
  for (int i = 0; i < this->NbrRow; ++i)
    {
      for (int j = 0; j < i; ++j)
	{
	  TmpMatrix[TotalIndex] = 0.0;
	  ++TotalIndex;
	}
      TmpMatrix[TotalIndex] = this->DiagonalElements[i];
      ++TotalIndex;
      for (int j = i + 1; j < this->NbrRow; ++j)
	{
	  TmpMatrix[TotalIndex] = this->OffDiagonalElements[Index2];
	  ++TotalIndex;
	  ++Index2;
	}
      Index2 += this->Increment;
    }
  FORTRAN_NAME(dsyev)(&Jobz, &UpperLower, &this->NbrRow, TmpMatrix, &this->NbrRow, M.DiagonalElements, &TmpWorkingArea, &WorkingAreaSize, &Information);
  WorkingAreaSize = (int) TmpWorkingArea;
  double* WorkingArea = new double [WorkingAreaSize];
  FORTRAN_NAME(dsyev)(&Jobz, &UpperLower, &this->NbrRow, TmpMatrix, &this->NbrRow, M.DiagonalElements, WorkingArea, &WorkingAreaSize, &Information);  
  delete[] WorkingArea;
  delete[] TmpMatrix;
  return M;
}

// Diagonalize a real symmetric matrix and evaluate transformation matrix using the LAPACK library (modifying current matrix)
//
// M = reference on real diagonal matrix where result has to be stored
// Q = matrix where transformation matrix has to be stored
// err = absolute error on matrix element
// maxIter = maximum number of iteration to fund an eigenvalue
// return value = reference on real tridiagonal symmetric matrix

RealDiagonalMatrix& RealSymmetricMatrix::LapackDiagonalize (RealDiagonalMatrix& M, RealMatrix& Q, double err, int maxIter)
{
  if (M.GetNbrRow() != this->NbrRow)
    M.Resize(this->NbrRow, this->NbrColumn);
  double* TmpMatrix = new double [this->NbrRow * this->NbrRow];
  int Information = 0;
  int WorkingAreaSize = -1;
  int IntegerWorkingAreaSize = -1;
  double TmpWorkingArea;
  int TmpIntegerWorkingArea;
  char Jobz = 'V';
  char UpperLower = 'L';
  int TotalIndex = 0;
  int Index2 = 0;
  for (int i = 0; i < this->NbrRow; ++i)
    {
      for (int j = 0; j < i; ++j)
	{
	  TmpMatrix[TotalIndex] = 0.0;
	  ++TotalIndex;
	}
      TmpMatrix[TotalIndex] = this->DiagonalElements[i];
      ++TotalIndex;
      for (int j = i + 1; j < this->NbrRow; ++j)
	{
	  TmpMatrix[TotalIndex] = this->OffDiagonalElements[Index2];
	  ++TotalIndex;
	  ++Index2;
	}
      Index2 += this->Increment;
     }
  FORTRAN_NAME(dsyevd)(&Jobz, &UpperLower, &this->NbrRow, TmpMatrix, &this->NbrRow, M.DiagonalElements, &TmpWorkingArea, &WorkingAreaSize, &TmpIntegerWorkingArea, &IntegerWorkingAreaSize, &Information);
  WorkingAreaSize = (int) TmpWorkingArea;
  double* WorkingArea = new double [WorkingAreaSize];
  IntegerWorkingAreaSize = TmpIntegerWorkingArea;
  int* IntegerWorkingArea = new int [IntegerWorkingAreaSize];
  FORTRAN_NAME(dsyevd)(&Jobz, &UpperLower, &this->NbrRow, TmpMatrix, &this->NbrRow, M.DiagonalElements, WorkingArea, &WorkingAreaSize, IntegerWorkingArea, &IntegerWorkingAreaSize, &Information);  
  TotalIndex = 0;
  for (int i = 0; i < this->NbrRow; ++i)
    for (int j = 0; j < this->NbrRow; ++j)
      {
	Q(j, i) = TmpMatrix[TotalIndex];
	++TotalIndex;
      }
  delete[] WorkingArea;
  delete[] IntegerWorkingArea;
  delete[] TmpMatrix;
  return M;
}

#endif


// output file stream overload
//
// file = reference on output file stream
// matrix = reference on matrix to save
// return value = reference on output file stream

/*ofstream& operator << (ofstream& file, const RealSymmetricMatrix& matrix)
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

ifstream& operator >> (ifstream& file, RealSymmetricMatrix& matrix)
{
  return file;
}

// add the symmetric matrix m1.m2.m3^t.m4^t + m4.m3.m2^t.m1^t to current symmetric matrix
//
// m1 = first matrix
// m2 = second matrix
// m3 = third matrix
// m4 = fourth matrix
// coefficient = optional global multiplicative factor in front of m1.m2.m3^t.m4^t + m4.m3.m2^t.m1^t

RealSymmetricMatrix& RealSymmetricMatrix::AddAAAtAt(RealMatrix& m1, RealMatrix& m2, RealMatrix& m3, RealMatrix& m4, double coefficient)
{
  return this->AddAAAtAt(m1, m2, m3, m4, 0, 0, (this->NbrRow * (this->NbrRow + 1)) >> 1, coefficient);
}

// add the symmetric matrix m1.m2.m3^t.m4^t + m4.m3.m2^t.m1^t to current symmetric matrix within a given range of indices
//
// m1 = first matrix
// m2 = second matrix
// m3 = third matrix
// m4 = fourth matrix
// rowIndex = row index of the first element to add
// columnIndex = column index of the first element to add
// nbrElement = number of element to add (starting from first element knowing that line i contains n - i + 1)
// coefficient = optional global multiplicative factor in front of m1.m2.m3^t.m4^t + m4.m3.m2^t.m1^t

RealSymmetricMatrix& RealSymmetricMatrix::AddAAAtAt(RealMatrix& m1, RealMatrix& m2, RealMatrix& m3, RealMatrix& m4, int rowIndex, int columnIndex,
						    int nbrElement, double coefficient)
{
  RealMatrix TmpMatrix (this->NbrRow, this->NbrRow);
  double TmpVal = 0.0;
  for (int i = 0; i < this->NbrRow; ++i)
    for (int j = 0; j < this->NbrRow; ++j)
      {
	TmpVal = 0.0;
	for (int k = 0; k < this->NbrRow; ++k)
	  TmpVal += m2(i, k) * m3(j, k);
	TmpMatrix(i,j) = TmpVal * coefficient;
      }
  if (rowIndex > columnIndex)
    {
      int tmp = rowIndex;
      rowIndex = columnIndex;
      columnIndex = tmp;
    }
  int pos;
  if (rowIndex == columnIndex)
    {
      pos = columnIndex - ((rowIndex * (rowIndex + 1)) >> 1) + rowIndex * (this->NbrRow + this->Increment - 1);
    }
  else
    {
      pos = columnIndex - ((rowIndex * (rowIndex + 1)) >> 1) + rowIndex * (this->NbrRow + this->Increment - 1) - 1;
    }
  while (nbrElement > 0)
    {
      TmpVal = 0.0;
      for (int i = 0; i < this->NbrRow; ++i)
	for (int j = 0; j < this->NbrRow; ++j)
	  TmpVal += m1(rowIndex, i) * TmpMatrix(i , j) * m4(columnIndex, j);
      if (rowIndex == columnIndex)
	this->DiagonalElements[rowIndex] += TmpVal;
      else
	{
	  this->OffDiagonalElements[pos] += TmpVal;
	  ++pos;
	}      
      ++columnIndex;
      if (columnIndex == this->NbrRow)
	{
	  ++rowIndex;
	  columnIndex = rowIndex;
	  pos += this->Increment;
	}
      --nbrElement;
    }
  return *this;
}


