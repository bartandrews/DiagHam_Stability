////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of real upper hessenberg matrix                  //
//                                                                            //
//                        last modification : 12/12/2012                      //
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


#include "Matrix/RealUpperHessenbergMatrix.h"
#include "Matrix/RealUpperTriangularMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "GeneralTools/ListIterator.h"
#include "MathTools/Complex.h"
#include "GeneralTools/Endian.h"

#include <stdlib.h>
#include <fstream>


using std::cout;
using std::endl;


#ifdef HAVE_LAPACK

// binding to the LAPACK function DHSEQR

extern "C" void FORTRAN_NAME(dhseqr)(const char* job, const char* computeZFlag, const int* nbrColumn, 
				     const int* triangularLowerIndex, const int* triangularHigherIndex, 
				     const double* matrix, const int* leadingDimension,
				     const double* eigenvaluesRealPart, const double* eigenvaluesImaginaryPart, 
				     const double* zMatrix, const int* leadingDimensionZMatrix,
				     const double* workingArea, const int* workingAreaSize, const int* information);

// binding to the LAPACK function DHSEIN
extern "C" void FORTRAN_NAME(dhsein)(const char* side, const char* eigenstateSource, const char* initalVectors, const int* selectEigenstates,
				     const int* nbrColumn, const double* matrix, const int* leadingDimension,
				     const double* eigenvaluesRealPart, const double* eigenvaluesImaginaryPart, 
				     const double* eigenvectorLeftMatrix, const int* leadingDimensionEigenvectorLeftMatrix,
				     const double* eigenvectorRightMatrix, const int* leadingDimensionEigenvectorRightMatrix,
				     const int* nbrColumEigenvectorMatrix, const int* nbrUsedColumEigenvectorMatrix,
				     const double* workingArea, const int* failedLeftEigenvectors, const int* failedRightEigenvectors,
				     const int* information);

// binding to the LAPACK function DGEQRF
extern "C" void FORTRAN_NAME(dgeqrf)(const int* nbrRow, const int* nbrColumn, const double* matrix, const int* leadingDimension,
				     const double* tau, const double* workingArea, const int* workingAreaSize,
				     const int* information);

// binding to the LAPACK function DORGQR
extern "C" void FORTRAN_NAME(dorgqr)(const int* nbrRow, const int* nbrColumn, const int* nbrReflectors, const double* matrix, const int* leadingDimension,
				     const double* tau, const double* workingArea, const int* workingAreaSize,
				     const int* information);

#endif


// default constructor
//

RealUpperHessenbergMatrix::RealUpperHessenbergMatrix() 
{
  this->UpperOffDiagonalElements = 0;
  this->DiagonalElements = 0;
  this->LowerDiagonalElements = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::Hessenberg | Matrix::Upper;
  this->Dummy = 0.0;
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

RealUpperHessenbergMatrix::RealUpperHessenbergMatrix(int dimension, bool zero) 
{
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::Hessenberg | Matrix::Upper;
  this->DiagonalElements = new double[this->NbrRow];
  this->LowerDiagonalElements = new double[this->NbrRow];
  long TmpNbrOffDiagonalElements = (((long) this->NbrRow) * (((long) this->NbrRow) - 1l)) / 2l;
  this->UpperOffDiagonalElements = new double[TmpNbrOffDiagonalElements];
  if (zero == true)
    {
      long pos = 0;
      for (int i = 0; i < this->NbrRow; i++)
	{
	  this->DiagonalElements[i] = 0.0;
	  this->LowerDiagonalElements[i] = 0.0;
	  for (int j = i + 1; j < this->NbrRow; j++)
	    {
	      this->UpperOffDiagonalElements[pos] = 0.0;
	      pos++;
	    }
	}
    }
  this->Dummy = 0.0;
}

// constructor from matrix elements (without duplicating datas)
//
// diagonalElements = diagonal elements
// offDiagonalElements = upper off-diagonal elements
// lowerDiagonalElements = lower diagonal elements
// dimension = matrix dimension

RealUpperHessenbergMatrix::RealUpperHessenbergMatrix(double* diagonalElements, double* offDiagonalElements, 
						     double* lowerDiagonalElements, int dimension) 
{
  this->DiagonalElements = diagonalElements;
  this->LowerDiagonalElements = lowerDiagonalElements;
  this->UpperOffDiagonalElements = offDiagonalElements;
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::Hessenberg | Matrix::Upper;
  this->Dummy = 0.0;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

RealUpperHessenbergMatrix::RealUpperHessenbergMatrix(const RealUpperHessenbergMatrix& M) 
{
  this->DiagonalElements = M.DiagonalElements;
  this->LowerDiagonalElements = M.LowerDiagonalElements;
  this->UpperOffDiagonalElements = M.UpperOffDiagonalElements;
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::Hessenberg | Matrix::Upper;
  this->Dummy = 0.0;
}

// constructor from a matrix, copying the data and discarding all elements below the lower diagonal
//
// M = reference on the matrix

RealUpperHessenbergMatrix::RealUpperHessenbergMatrix(Matrix& M) 
{
  this->Flag.Initialize();  
  if (M.GetNbrColumn() >= M.GetNbrRow())
    {
      this->NbrRow = M.GetNbrRow();
     this->NbrColumn = M.GetNbrRow();
     }
  else
    {    
      this->NbrRow = M.GetNbrColumn();
      this->NbrColumn = M.GetNbrColumn();
    }
  this->TrueNbrRow = this->TrueNbrRow;
  this->TrueNbrColumn = this->TrueNbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::Hessenberg | Matrix::Upper;
  this->Dummy = 0.0;
  this->DiagonalElements = new double[this->NbrRow];
  this->LowerDiagonalElements = new double[this->NbrRow];
  long TmpNbrOffDiagonalElements = (((long) this->NbrRow) * (((long) this->NbrRow) - 1l)) / 2l;
  this->UpperOffDiagonalElements = new double[TmpNbrOffDiagonalElements];
  double Tmp;
  for (int i = 0; i < this->NbrRow; ++i)
    {
      M.GetMatrixElement(i, i, Tmp);
      this->DiagonalElements[i] = Tmp;
      if (i > 0)
	M.GetMatrixElement(i, i - 1, Tmp);
      this->LowerDiagonalElements[i] = Tmp;
      for (int j = i + 1; j < this->NbrRow; ++j)
	{
	  M.GetMatrixElement(i, j, Tmp);
 	  this->UpperOffDiagonalElements[i + ((j * (j - 1l)) >> 1)] = Tmp;
	}
    }  
}

// destructor
//

RealUpperHessenbergMatrix::~RealUpperHessenbergMatrix() 
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      if (this->DiagonalElements != 0)
	{
	  delete[] this->DiagonalElements;
	}
      if (this->LowerDiagonalElements != 0)
	{
	  delete[] this->LowerDiagonalElements;
	}
      if (this->UpperOffDiagonalElements != 0)
	{
	  delete[] this->UpperOffDiagonalElements;
	}
    }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

RealUpperHessenbergMatrix& RealUpperHessenbergMatrix::operator = (const RealUpperHessenbergMatrix& M) 
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      if (this->DiagonalElements != 0)
	{
	  delete[] this->DiagonalElements;
	}
      if (this->LowerDiagonalElements != 0)
	{
	  delete[] this->LowerDiagonalElements;
	}
      if (this->UpperOffDiagonalElements != 0)
	{
	  delete[] this->UpperOffDiagonalElements;
	}
    }
  this->DiagonalElements = M.DiagonalElements;
  this->LowerDiagonalElements = M.LowerDiagonalElements;
  this->UpperOffDiagonalElements = M.UpperOffDiagonalElements;
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::Hessenberg | Matrix::Upper;
  this->Dummy = 0.0;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* RealUpperHessenbergMatrix::Clone ()
{
  return ((Matrix*) new RealUpperHessenbergMatrix (*this));
}

// copy a matrix into another (duplicating data)
//
// matrix = matrix to copy
// return value = reference on current matrix

RealUpperHessenbergMatrix& RealUpperHessenbergMatrix::Copy (RealUpperHessenbergMatrix& matrix)
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      if (this->DiagonalElements != 0)
	{
	  delete[] this->DiagonalElements;
	}
      if (this->LowerDiagonalElements != 0)
	{
	  delete[] this->LowerDiagonalElements;
	}
      if (this->UpperOffDiagonalElements != 0)
	{
	  delete[] this->UpperOffDiagonalElements;
	}
    }
  if (matrix.Flag.Used() == true)
    {
      this->Flag.Initialize();
      this->NbrRow = matrix.NbrRow;
      this->NbrColumn = matrix.NbrColumn;
      this->TrueNbrRow = matrix.TrueNbrRow;
      this->TrueNbrColumn = matrix.TrueNbrColumn;
      this->MatrixType = Matrix::RealElements | Matrix::Hessenberg | Matrix::Upper;
      this->Dummy = 0.0;
      this->DiagonalElements = new double[this->NbrRow];
      this->LowerDiagonalElements = new double[this->NbrRow];
      long TmpNbrOffDiagonalElements = (((long) this->NbrRow) * (((long) this->NbrRow) - 1l)) / 2l;
      this->UpperOffDiagonalElements = new double[TmpNbrOffDiagonalElements];
      long pos = 0;
      for (int i = 0; i < this->NbrRow; i++)
	{
	  this->DiagonalElements[i] = matrix.DiagonalElements[i];
	  this->LowerDiagonalElements[i] = matrix.LowerDiagonalElements[i];
	  for (int j = i + 1; j < this->NbrRow; j++)
	    {
	      this->UpperOffDiagonalElements[pos] = matrix.UpperOffDiagonalElements[pos];
	      pos++;
	    }
	}
    }
  else
    {
      this->UpperOffDiagonalElements = 0;
      this->DiagonalElements = 0;
      this->LowerDiagonalElements = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;
      this->MatrixType = Matrix::RealElements | Matrix::Hessenberg | Matrix::Upper;
      this->Dummy = 0.0;
   }
  return *this;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void RealUpperHessenbergMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      this->DiagonalElements[i] = x;
    }
  else
    {
      if (i < j)
	{
	  this->UpperOffDiagonalElements[i + (j * (j - 1l)) / 2l] = x;
	}
      else
	{
	  if ((j + 1) == i)
	    {
	      this->LowerDiagonalElements[i] = x;
	    }
	}
    }      
  return;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void RealUpperHessenbergMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  if (i == j)
    {
      this->DiagonalElements[i] = x.Re;
    }
  else
    {
      if (i < j)
	{
	  this->UpperOffDiagonalElements[i + (j * (j - 1l)) / 2l] = x.Re;
	}
      else
	{
	  if ((j + 1) == i)
	    {
	      this->LowerDiagonalElements[i] = x.Re;
	    }
	}
    }
  return;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void RealUpperHessenbergMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      this->DiagonalElements[i] = x;
    }
  else
    {
      if (i < j)
	{
	  this->UpperOffDiagonalElements[i + (j * (j - 1l)) / 2l] += x;
	}
      else
	{
	  if (i == (j + 1))
	    {
	      this->LowerDiagonalElements[i] += x;
	    }
	}
    }      
  return;
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void RealUpperHessenbergMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      this->DiagonalElements[i] = x.Re;
    }
  else
    {
      if (i < j)
	{
	  this->UpperOffDiagonalElements[i + (j * (j - 1l)) / 2l] += x.Re;
	}
      else
	{
	  if ((j + 1) == i)
	    {
	      this->LowerDiagonalElements[i] += x.Re;
	    }
	}
    }      
  return;
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void RealUpperHessenbergMatrix::Resize (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  double* TmpDiagonalElements = new double[nbrRow];
  double* TmpLowerDiagonalElements = new double[nbrRow];
  long TmpNbrOffDiagonalElements = (((long) nbrRow) * (((long) nbrRow) - 1l)) / 2l;
  double* TmpUpperOffDiagonalElements = new double[TmpNbrOffDiagonalElements];
  int MinNbrRow = nbrRow;
  if (nbrRow > this->NbrRow)
    MinNbrRow = this->NbrRow;
  for (int i = 0; i < MinNbrRow; ++i)
    {
      TmpDiagonalElements[i] = this->DiagonalElements[i];
      TmpLowerDiagonalElements[i] = this->LowerDiagonalElements[i];
      for (int j = i + 1; j < MinNbrRow; ++j)
	{
	  long Index = i + (j * (j - 1l)) / 2l;
	  TmpUpperOffDiagonalElements[Index] = this->UpperOffDiagonalElements[Index];
	}
    }
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      if (this->DiagonalElements != 0)
	{
	  delete[] this->DiagonalElements;
	}
      if (this->LowerDiagonalElements != 0)
	{
	  delete[] this->LowerDiagonalElements;
	}
      if (this->UpperOffDiagonalElements != 0)
	{
	  delete[] this->UpperOffDiagonalElements;
	}
    }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->DiagonalElements = TmpDiagonalElements;
  this->LowerDiagonalElements = TmpLowerDiagonalElements;
  this->UpperOffDiagonalElements = TmpUpperOffDiagonalElements;
  this->Flag.Initialize();
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void RealUpperHessenbergMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  double* TmpDiagonalElements = new double[nbrRow];
  double* TmpLowerDiagonalElements = new double[nbrRow];
  long TmpNbrOffDiagonalElements = (((long) nbrRow) * (((long) nbrRow) - 1l)) / 2l;
  double* TmpUpperOffDiagonalElements = new double[TmpNbrOffDiagonalElements];
  int MinNbrRow = nbrRow;
  if (nbrRow > this->NbrRow)
    MinNbrRow = this->NbrRow;
  for (int i = 0; i < MinNbrRow; ++i)
    {
      TmpDiagonalElements[i] = this->DiagonalElements[i];
      TmpLowerDiagonalElements[i] = this->LowerDiagonalElements[i];
      for (int j = i + 1; j < MinNbrRow; ++j)
	{
	  long Index = i + (j * (j - 1l)) / 2l;
	  TmpUpperOffDiagonalElements[Index] = this->UpperOffDiagonalElements[Index];
	}
      for (int j = MinNbrRow; j < nbrRow; ++j)
	{
	  long Index = i + (j * (j - 1l)) / 2l;
	  TmpUpperOffDiagonalElements[Index] = 0.0;
	}
    }
  for (int i = this->NbrRow; i < nbrRow; ++i)
    {
      TmpDiagonalElements[i] = 0.0;
      TmpLowerDiagonalElements[i] = 0.0;
      for (int j = i + 1; j < nbrRow; ++j)
	{
	  TmpUpperOffDiagonalElements[i + (j * (j - 1l)) / 2l] = 0.0;
	}
    }
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      if (this->DiagonalElements != 0)
	{
	  delete[] this->DiagonalElements;
	}
      if (this->LowerDiagonalElements != 0)
	{
	  delete[] this->LowerDiagonalElements;
	}
      if (this->UpperOffDiagonalElements != 0)
	{
	  delete[] this->UpperOffDiagonalElements;
	}
    }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->DiagonalElements = TmpDiagonalElements;
  this->LowerDiagonalElements = TmpLowerDiagonalElements;
  this->UpperOffDiagonalElements = TmpUpperOffDiagonalElements;
  this->Flag.Initialize();
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

RealUpperHessenbergMatrix operator + (const RealUpperHessenbergMatrix& M1, const RealUpperHessenbergMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return RealUpperHessenbergMatrix();
  double* TmpDiagonalElements = new double[M1.NbrRow];
  double* TmpLowerDiagonalElements = new double[M1.NbrRow];
  long TmpNbrOffDiagonalElements = (((long) M1.NbrRow) * (((long) M1.NbrRow) - 1l)) / 2l;
  double* TmpUpperOffDiagonalElements = new double[TmpNbrOffDiagonalElements];
  long pos = 0;
  for (int i = 0; i < M1.NbrRow; i++)
    {
      TmpDiagonalElements[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
      TmpLowerDiagonalElements[i] = M1.LowerDiagonalElements[i] + M2.LowerDiagonalElements[i];
      for (int j = i + 1; j < M1.NbrRow; j++)
	{
	  TmpUpperOffDiagonalElements[pos] = M1.UpperOffDiagonalElements[pos] + M2.UpperOffDiagonalElements[pos];
	  ++pos;
	}
    }
  return RealUpperHessenbergMatrix(TmpDiagonalElements, TmpUpperOffDiagonalElements, TmpLowerDiagonalElements, M1.NbrRow);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

RealUpperHessenbergMatrix operator - (const RealUpperHessenbergMatrix& M1, const RealUpperHessenbergMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return RealUpperHessenbergMatrix();
  double* TmpDiagonalElements = new double[M1.NbrRow];
  double* TmpLowerDiagonalElements = new double[M1.NbrRow];
  long TmpNbrOffDiagonalElements = (((long) M1.NbrRow) * (((long) M1.NbrRow) - 1l)) / 2l;
  double* TmpUpperOffDiagonalElements = new double[TmpNbrOffDiagonalElements];
  long pos = 0;
  for (int i = 0; i < M1.NbrRow; i++)
    {
      TmpDiagonalElements[i] = M1.DiagonalElements[i] - M2.DiagonalElements[i];
      TmpLowerDiagonalElements[i] = M1.LowerDiagonalElements[i]  - M2.LowerDiagonalElements[i];
      for (int j = i + 1; j < M1.NbrRow; j++)
	{
	  TmpUpperOffDiagonalElements[pos] = M1.UpperOffDiagonalElements[pos] - M2.UpperOffDiagonalElements[pos];
	  ++pos;
	}
    }
  return RealUpperHessenbergMatrix(TmpDiagonalElements, TmpUpperOffDiagonalElements, TmpLowerDiagonalElements, M1.NbrRow);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealUpperHessenbergMatrix operator * (const RealUpperHessenbergMatrix& M, double x) 
{
  double* TmpDiagonalElements = new double[M.NbrRow];
  double* TmpLowerDiagonalElements = new double[M.NbrRow];
  long TmpNbrOffDiagonalElements = (((long) M.NbrRow) * (((long) M.NbrRow) - 1l)) / 2l;
  double* TmpUpperOffDiagonalElements = new double[TmpNbrOffDiagonalElements];
  long pos = 0;
  for (int i = 0; i < M.NbrRow; i++)
    {
      TmpDiagonalElements[i] = M.DiagonalElements[i] * x;
      TmpLowerDiagonalElements[i] = M.LowerDiagonalElements[i] * x;
      for (int j = i + 1; j < M.NbrRow; j++)
	{
	  TmpUpperOffDiagonalElements[pos] = M.UpperOffDiagonalElements[pos] * x;
	  ++pos;
	}
    }
  return RealUpperHessenbergMatrix(TmpDiagonalElements, TmpUpperOffDiagonalElements, TmpLowerDiagonalElements, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealUpperHessenbergMatrix operator * (double x, const RealUpperHessenbergMatrix& M) 
{
  return (M * x);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

RealUpperHessenbergMatrix operator / (const RealUpperHessenbergMatrix& M, double x) 
{
  x = 1.0 / x;
  return (M * x);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

RealUpperHessenbergMatrix& RealUpperHessenbergMatrix::operator += (const RealUpperHessenbergMatrix& M) 
{
  if (this->NbrRow == 0)
    return *this;
  long pos = 0;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] += M.DiagonalElements[i];
      this->LowerDiagonalElements[i] += M.LowerDiagonalElements[i];
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  this->UpperOffDiagonalElements[pos] += M.UpperOffDiagonalElements[pos];
	  ++pos;
	}
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

RealUpperHessenbergMatrix& RealUpperHessenbergMatrix::operator -= (const RealUpperHessenbergMatrix& M) 
{
  if (this->NbrRow == 0)
    return *this;
  long pos = 0;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] -= M.DiagonalElements[i];
      this->LowerDiagonalElements[i] -= M.LowerDiagonalElements[i];
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  this->UpperOffDiagonalElements[pos] -= M.UpperOffDiagonalElements[pos];
	  ++pos;
	}
    }
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

RealUpperHessenbergMatrix& RealUpperHessenbergMatrix::operator *= (double x) 
{
  long pos = 0;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] *= x;
      this->LowerDiagonalElements[i] *= x;
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  this->UpperOffDiagonalElements[pos] *= x;
	  ++pos;
	}
    }
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

RealUpperHessenbergMatrix& RealUpperHessenbergMatrix::operator /= (double x)
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

double RealUpperHessenbergMatrix::MatrixElement (RealVector& V1, RealVector& V2)
{
  double x = 0.0;
  if ((V1.GetVectorDimension() != this->NbrRow) || (V2.GetVectorDimension() != this->NbrColumn))
    return x;
  for (int i = 0; i < this->NbrRow ; i++)
    {
      double x2 = this->DiagonalElements[i] * V2[i];
      x2 += this->LowerDiagonalElements[i] * V2[i - 1];
      for (int j = i + 1; j < this->NbrColumn; ++j)
	{
	  x2 +=  this->UpperOffDiagonalElements[i + (j * (j - 1l)) / 2l] * V2[j];
	}
      x += V1[i] * x2;
    }
  return x;
}

// evaluate matrix trace
//
// return value = matrix trace 

double RealUpperHessenbergMatrix::Tr () 
{
  double Trace = 0.0;
  for (int i = 0; i <  this->NbrRow; ++i)
    Trace += this->DiagonalElements[i];
  return Trace;
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double RealUpperHessenbergMatrix::Det () 
{
  return 0.0;
}

// shift all diagonal elements 
//
// shift = shift to apply
// return value = reference on current matrix

RealUpperHessenbergMatrix& RealUpperHessenbergMatrix::ShiftDiagonal(double shift)
{
  for (int i = 0; i <  this->NbrRow; ++i)
    this->DiagonalElements[i] += shift;
  return *this;
}

// conjugate matrix with an unitary matrix (Ut M U), assuming the Hessenberg from will be preserved
//
// unitaryM = unitary matrix to use
// conjugatedMatrix = reference on the matrix where conjugate matrix will be stored
// return value = pointer to conjugated matrix

RealUpperHessenbergMatrix& RealUpperHessenbergMatrix::Conjugate(RealMatrix& unitaryM, RealUpperHessenbergMatrix& conjugatedMatrix)
{
  if ((unitaryM.GetNbrRow() != this->NbrColumn) || (conjugatedMatrix.NbrColumn != this->NbrColumn))
    return conjugatedMatrix;
  for (int i = 0; i < this->NbrRow; ++i)
    {
      int j = i - 1;
      if (j < 0)
	j = 0;
      for (; j < this->NbrRow; ++j)
	{
	  double Tmp = 0.0;
	  for (int k = 0; k < this->NbrRow; ++k)
	    {
	      double Tmp2 = unitaryM[j][k] * this->DiagonalElements[k];
	      if (k > 0)
		Tmp2 += unitaryM[j][k - 1] * this->LowerDiagonalElements[k];	      
	      for (int l = k + 1; l < this->NbrColumn; ++l)
		{
		  Tmp2 += unitaryM[j][l] * this->UpperOffDiagonalElements[k + (l * (l - 1l)) / 2l];
		}
	      Tmp += Tmp2 * unitaryM[i][k];
	    }
	  conjugatedMatrix.SetMatrixElement(i, j , Tmp);
	}
    }    
  return conjugatedMatrix;
}

// Diagonalize a real matrix using the LAPACK library
//
// M = reference on complex diagonal matrix where result has to be stored
// leftFlag = compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
// return value = reference on complex diagonal matrix

ComplexDiagonalMatrix& RealUpperHessenbergMatrix::LapackDiagonalize (ComplexDiagonalMatrix& M, bool leftFlag)
{
#ifdef HAVE_LAPACK
   int Information = 0;
   int WorkingAreaSize = -1;
   char Job = 'E';
   char computeZFlag = 'N';
   int TriangularLowerIndex = 1;
   int TriangularHigherIndex = this->NbrColumn;
   double* TmpMatrix = new double [((long) this->NbrColumn) * this->NbrRow];
   long TotalIndex = 0l;
   for (int j = 0; j < this->NbrColumn; ++j)
     {
       for (int i = 0; i < j; ++i)
 	{
 	  TmpMatrix[TotalIndex] = this->UpperOffDiagonalElements[i + ((j * (j - 1l)) >> 1)];
 	  ++TotalIndex;
 	}
       TmpMatrix[TotalIndex] = this->DiagonalElements[j];
       ++TotalIndex;
       if ((j + 1) < this->NbrColumn)
	 {
	   TmpMatrix[TotalIndex] = this->LowerDiagonalElements[j + 1];
	   ++TotalIndex;
	 }      
       for (int i = j + 2; i < this->NbrRow; ++i)
	 {
	   TmpMatrix[TotalIndex] = 0.0;
	   ++TotalIndex;
	 }
     }
   double* TmpEigenvalueReal = new double[this->NbrColumn];
   double* TmpEigenvalueImaginary = new double[this->NbrColumn];
   int TmpLeadingDimension = 1;
   double* Dummy = 0;
   double TmpWorkingArea;
   FORTRAN_NAME(dhseqr)(&Job, &computeZFlag, &this->NbrRow, &TriangularLowerIndex, &TriangularHigherIndex, TmpMatrix, &this->NbrRow, 
			TmpEigenvalueReal, TmpEigenvalueImaginary,
			Dummy, &TmpLeadingDimension, &TmpWorkingArea, &WorkingAreaSize, &Information);
   WorkingAreaSize = (int) TmpWorkingArea;
   double* WorkingArea = new double [WorkingAreaSize];
   FORTRAN_NAME(dhseqr)(&Job, &computeZFlag, &this->NbrRow, &TriangularLowerIndex, &TriangularHigherIndex, TmpMatrix, &this->NbrRow, 
			TmpEigenvalueReal, TmpEigenvalueImaginary,
			Dummy, &TmpLeadingDimension, WorkingArea, &WorkingAreaSize, &Information);
   for (int i = 0; i < this->NbrRow; ++i)
     {
       M[i].Re = TmpEigenvalueReal[i];
       M[i].Im = TmpEigenvalueImaginary[i];
     }
   delete[] WorkingArea;
   delete[] TmpEigenvalueReal;
   delete[] TmpEigenvalueImaginary;
   delete[] TmpMatrix;
#endif
  return M;
}

// Diagonalize a real matrix and evaluate the left eigenstates using the LAPACK library
//
// M = reference on complex diagonal matrix where result has to be stored
// Q = matrix where transformation matrix has to be stored
// leftFlag = compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
// return value = reference on complex diagonal matrix

ComplexDiagonalMatrix& RealUpperHessenbergMatrix::LapackDiagonalize (ComplexDiagonalMatrix& M, ComplexMatrix& Q, bool leftFlag)
{
#ifdef HAVE_LAPACK
  int Information = 0;
  int WorkingAreaSize = -1;
  char Job = 'E';
  char computeZFlag = 'N';
  int TriangularLowerIndex = 1;
  int TriangularHigherIndex = this->NbrColumn;
  double* TmpMatrix = new double [((long) this->NbrColumn) * this->NbrRow];
  long TotalIndex = 0l;
  for (int j = 0; j < this->NbrColumn; ++j)
    {
      for (int i = 0; i < j; ++i)
 	{
 	  TmpMatrix[TotalIndex] = this->UpperOffDiagonalElements[i + ((j * (j - 1l)) >> 1)];
 	  ++TotalIndex;
 	}
      TmpMatrix[TotalIndex] = this->DiagonalElements[j];
      ++TotalIndex;
      if ((j + 1) < this->NbrColumn)
	{
	  TmpMatrix[TotalIndex] = this->LowerDiagonalElements[j + 1];
	  ++TotalIndex;
	}      
      for (int i = j + 2; i < this->NbrRow; ++i)
	{
	  TmpMatrix[TotalIndex] = 0.0;
	  ++TotalIndex;
	}
    }
  double* TmpEigenvalueReal = new double[this->NbrColumn];
  double* TmpEigenvalueImaginary = new double[this->NbrColumn];
  int TmpLeadingDimension = 1;
  double* Dummy = 0;
  double TmpWorkingArea;
  FORTRAN_NAME(dhseqr)(&Job, &computeZFlag, &this->NbrRow, &TriangularLowerIndex, &TriangularHigherIndex, TmpMatrix, &this->NbrRow, 
		       TmpEigenvalueReal, TmpEigenvalueImaginary,
		       Dummy, &TmpLeadingDimension, &TmpWorkingArea, &WorkingAreaSize, &Information);
  WorkingAreaSize = (int) TmpWorkingArea;
  double* WorkingArea = new double [WorkingAreaSize];
  FORTRAN_NAME(dhseqr)(&Job, &computeZFlag, &this->NbrRow, &TriangularLowerIndex, &TriangularHigherIndex, TmpMatrix, &this->NbrRow, 
		       TmpEigenvalueReal, TmpEigenvalueImaginary,
		       Dummy, &TmpLeadingDimension, WorkingArea, &WorkingAreaSize, &Information);
  delete[] WorkingArea;
  
  
  char Side = 'R';
  if (leftFlag == true)
    Side = 'L';
  char EigenvalueSource = 'Q';
  char InitialEigenvectors = 'N';
  int* SelectEigenvectors = new int[this->NbrColumn];
  int* FailedLeftEigenvectors = new int [this->NbrColumn];
  double* TmpLeftEigenstates = new double [((long) this->NbrColumn) * this->NbrRow];
  WorkingArea = new double [this->NbrColumn * (this->NbrColumn + 2l)];
  for (int j = 0; j < this->NbrColumn; ++j)
    SelectEigenvectors[j] = 1;
  int TmpLeadingLeftDimension;
  int TmpLeadingRightDimension;
  int NbrRequiredColumns = 0;
  if (leftFlag == true)
    {
      TmpLeadingLeftDimension = this->NbrColumn;
      TmpLeadingRightDimension = 1;
    }
  else
    {
      TmpLeadingRightDimension = this->NbrColumn;
      TmpLeadingLeftDimension = 1;
    }
  
  TotalIndex = 0l;
  for (int j = 0; j < this->NbrColumn; ++j)
    {
      for (int i = 0; i < j; ++i)
 	{
 	  TmpMatrix[TotalIndex] = this->UpperOffDiagonalElements[i + ((j * (j - 1l)) >> 1)];
 	  ++TotalIndex;
 	}
      TmpMatrix[TotalIndex] = this->DiagonalElements[j];
      ++TotalIndex;
      if ((j + 1) < this->NbrColumn)
	{
	  TmpMatrix[TotalIndex] = this->LowerDiagonalElements[j + 1];
	  ++TotalIndex;
	}      
      for (int i = j + 2; i < this->NbrRow; ++i)
	{
	  TmpMatrix[TotalIndex] = 0.0;
	  ++TotalIndex;
	}
    }
  
  if (leftFlag == true)
    {
      cout << "left " << this->NbrRow << endl;
      FORTRAN_NAME(dhsein)(&Side, &EigenvalueSource, &InitialEigenvectors, SelectEigenvectors,
			   &this->NbrRow, TmpMatrix, &this->NbrColumn, 
			   TmpEigenvalueReal, TmpEigenvalueImaginary,
			   TmpLeftEigenstates, &TmpLeadingLeftDimension, Dummy, &TmpLeadingRightDimension, 
			   &this->NbrRow, &NbrRequiredColumns, WorkingArea,
			   FailedLeftEigenvectors, FailedLeftEigenvectors, &Information);
    }
  else
    {
      FORTRAN_NAME(dhsein)(&Side, &EigenvalueSource, &InitialEigenvectors, SelectEigenvectors,
			   &this->NbrRow, TmpMatrix, &this->NbrColumn, 
			   TmpEigenvalueReal, TmpEigenvalueImaginary,
			   Dummy, &TmpLeadingLeftDimension, TmpLeftEigenstates, &TmpLeadingRightDimension, 
			   &this->NbrRow, &NbrRequiredColumns, WorkingArea,
			   FailedLeftEigenvectors, FailedLeftEigenvectors, &Information);
    }
  for (int i = 0; i < this->NbrRow; ++i)
    {
      M[i].Re = TmpEigenvalueReal[i];
      M[i].Im = TmpEigenvalueImaginary[i];
    }
  TotalIndex = 0l;
  for (int i = 0; i < this->NbrRow;)
    {
      if ((i == (this->NbrRow - 1)) || (TmpEigenvalueImaginary[i] != -TmpEigenvalueImaginary[i + 1]) 
	  || (TmpEigenvalueImaginary[i] == 0.0))
	{
	  Complex Tmp;
	  for (int j = 0; j < this->NbrRow; ++j)
	    {
	      Tmp.Re = TmpLeftEigenstates[TotalIndex];
	      Q.SetMatrixElement(j, i, Tmp);
	      ++TotalIndex;
	    }
	  ++i;
	}
      else
	{
	  Complex Tmp;
	  for (int j = 0; j < this->NbrRow; ++j)
	    {
	      Tmp.Re = TmpLeftEigenstates[TotalIndex];
	      Q.SetMatrixElement(j, i, Tmp);
	      Q.SetMatrixElement(j, i + 1, Tmp);
	      ++TotalIndex;
	    }
	  for (int j = 0; j < this->NbrRow; ++j)
	    {
	      Q.GetMatrixElement(j, i, Tmp);
	      Tmp.Im = TmpLeftEigenstates[TotalIndex];
	      Q.SetMatrixElement(j, i, Tmp);
	      Tmp.Im = -TmpLeftEigenstates[TotalIndex];
	      Q.SetMatrixElement(j, i + 1, Tmp);
	      ++TotalIndex;
	    }
	  i += 2;
	}
    }
  delete[] TmpEigenvalueReal;
  delete[] TmpEigenvalueImaginary;
  delete[] TmpMatrix;
  delete[] SelectEigenvectors;
  delete[] FailedLeftEigenvectors;
  delete[] TmpLeftEigenstates;
#endif
  return M;
}

// find QR factorization using the LAPACK library
//
// R = reference on the triangular matrix
// Q = reference on the transformation matrix
// return value = reference on upper triangular matrix

RealUpperTriangularMatrix& RealUpperHessenbergMatrix::LapackQRFactorization (RealUpperTriangularMatrix& R, RealMatrix& Q)
{
#ifdef HAVE_LAPACK
  int Information = 0;
  int WorkingAreaSize = -1;
  double* TmpMatrix = new double [((long) this->NbrColumn) * this->NbrRow];
  double* TmpTau = new double [this->NbrColumn];
  long TotalIndex = 0l;
  for (int j = 0; j < this->NbrColumn; ++j)
    {
      for (int i = 0; i < j; ++i)
 	{
 	  TmpMatrix[TotalIndex] = this->UpperOffDiagonalElements[i + ((j * (j - 1l)) >> 1)];
 	  ++TotalIndex;
 	}
      TmpMatrix[TotalIndex] = this->DiagonalElements[j];
      ++TotalIndex;
      if ((j + 1) < this->NbrColumn)
	{
	  TmpMatrix[TotalIndex] = this->LowerDiagonalElements[j + 1];
	  ++TotalIndex;
	}      
      for (int i = j + 2; i < this->NbrRow; ++i)
	{
	  TmpMatrix[TotalIndex] = 0.0;
	  ++TotalIndex;
	}
    }
  double TmpWorkingArea;
  FORTRAN_NAME(dgeqrf)(&this->NbrRow, &this->NbrColumn, TmpMatrix, &this->NbrColumn,
		       TmpTau, &TmpWorkingArea, &WorkingAreaSize,
		       &Information);
  WorkingAreaSize = (int) TmpWorkingArea;
  double* WorkingArea = new double [WorkingAreaSize];
  FORTRAN_NAME(dgeqrf)(&this->NbrRow, &this->NbrColumn, TmpMatrix, &this->NbrColumn,
		       TmpTau, WorkingArea, &WorkingAreaSize,
		       &Information);
  delete[] WorkingArea;
  TotalIndex = 0l;
  for (int i = 0; i < this->NbrColumn; ++i)
    {
      for (int j = 0; j < i; ++j)
	{
	  R.SetMatrixElement(j, i, TmpMatrix[TotalIndex]);
	  ++TotalIndex;
	}      
      R.SetMatrixElement(i, i, TmpMatrix[TotalIndex]);
      TotalIndex += (long) (this->NbrRow - i);
   }
  WorkingAreaSize = -1;
  FORTRAN_NAME(dorgqr)(&this->NbrRow, &this->NbrColumn, &this->NbrColumn, TmpMatrix, &this->NbrColumn,
		       TmpTau, &TmpWorkingArea, &WorkingAreaSize,
		       &Information);
  WorkingAreaSize = (int) TmpWorkingArea;
  WorkingArea = new double [WorkingAreaSize];
  FORTRAN_NAME(dorgqr)(&this->NbrRow, &this->NbrColumn, &this->NbrColumn, TmpMatrix, &this->NbrColumn,
		       TmpTau, WorkingArea, &WorkingAreaSize,
		       &Information);

  TotalIndex = 0l;
  for (int i = 0; i < this->NbrColumn; ++i)
    {
      for (int j = 0; j < this->NbrRow; ++j)
	{
	  Q.SetMatrixElement(j, i, TmpMatrix[TotalIndex]);
	  ++TotalIndex;
	}      
    }
  delete[] WorkingArea;
  delete[] TmpMatrix;
  delete[] TmpTau;
#endif
  return R;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const RealUpperHessenbergMatrix& P)
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      for (int j = 0; j < (i - 1); ++j)
	{
	  Str << 0.0 << "    ";
	}
      if (i > 0)
	{
	  Str << P.LowerDiagonalElements[i] << "    ";
	}
      Str << P.DiagonalElements[i] << "    ";
      for (int j = i + 1; j < P.NbrRow; ++j)
	{
	  Str << P.UpperOffDiagonalElements[i + (((j - 1l) * j) >> 1)] << "    ";
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

MathematicaOutput& operator << (MathematicaOutput& Str, const RealUpperHessenbergMatrix& P)
{
  Str << "{";
  for (int i = 0; i < P.NbrRow; ++i)
    {
      Str << "{";
      int pos = i - 1;
      for (int j = 0; j < (i - 1); ++j)
	{
	  Str << "0,";
	}
      if (i > 0)
	{
	  Str << P.LowerDiagonalElements[i];
	}
      for (int j = i + 1; j < P.NbrRow; j++)
	{
	  Str << "," << P.UpperOffDiagonalElements[i + ((j - 1l) * j) / 2l];
	}
      if (i != (P.NbrRow - 1))
	Str << "},";
      else
	Str << "}";
    }
  Str << "}";
  return Str;
}

// write matrix in a file 
//
// file = reference on the output file stream
// return value = true if no error occurs

bool RealUpperHessenbergMatrix::WriteMatrix (ofstream& file)
{
  WriteLittleEndian(file, this->MatrixType);
  WriteLittleEndian(file, this->NbrRow);
  WriteLittleEndian(file, this->NbrColumn);
  for (int i = 0; i < this->NbrRow; ++i)
    {
      WriteLittleEndian(file, this->LowerDiagonalElements[i]);
      WriteLittleEndian(file, this->DiagonalElements[i]);
    }
  long TmpNbrOffDiagonalElements = (((long) this->NbrRow) * (((long) this->NbrRow) - 1l)) / 2l;
  for (long i = 0l; i < TmpNbrOffDiagonalElements; ++i)
    WriteLittleEndian(file, this->UpperOffDiagonalElements[i]);
  return true;
}

// read matrix from a file 
//
// file = reference  on the input file stream
// return value = true if no error occurs

bool RealUpperHessenbergMatrix::ReadMatrix (ifstream& file)
{
  int TmpType = Matrix::RealElements;
  ReadLittleEndian(file, TmpType);
  if (TmpType != (Matrix::RealElements | Matrix::Hessenberg | Matrix::Upper))
    {
      return false;
    }
  int TmpNbrRow;
  int TmpNbrColumn;
  ReadLittleEndian(file, TmpNbrRow);
  ReadLittleEndian(file, TmpNbrColumn);
  this->Resize(TmpNbrRow, TmpNbrColumn);
  for (int i = 0; i < this->NbrRow; ++i)
    {
      ReadLittleEndian(file, this->LowerDiagonalElements[i]);
      ReadLittleEndian(file, this->DiagonalElements[i]);
    }
  long TmpNbrOffDiagonalElements = (((long) this->NbrRow) * (((long) this->NbrRow) - 1l)) / 2l;
  for (long i = 0l; i < TmpNbrOffDiagonalElements; ++i)
    ReadLittleEndian(file, this->UpperOffDiagonalElements[i]);
  return true;
}

#endif
