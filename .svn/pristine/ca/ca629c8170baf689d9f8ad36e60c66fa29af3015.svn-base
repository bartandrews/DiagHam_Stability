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


#include "Matrix/ComplexMatrix.h"
#include "Vector/ComplexVector.h"


using std::endl;
using std::cout;


#ifdef HAVE_LAPACK

// binding to the LAPACK zgetrf routine for LU decomposition and back-substitution
//
extern "C" void FORTRAN_NAME(zgetrf)(const int* dimensionM, const int* dimensionN, const doublecomplex* matrixA,
				     const int* leadingDimensionA, const int *ipiv, const int *info);

extern "C" void FORTRAN_NAME(zgetrs)(const char* transpose, const int* dimensionN, const int* numRHS,
				     const doublecomplex* matrixA, const int* leadingDimensionA, const int *ipiv,
				     const doublecomplex* matrixB, const int* leadingDimensionB, const int *info);

// sequence of routines to be called for diagonalization via Hessenberg matrix QR decomposition
// balance (optional)
extern "C" void FORTRAN_NAME(zgebal)(const char* jobz, const int* dimensionN, const doublecomplex* matrixA,
				     const int* leadingDimensionA, const int *iLow, const int *iHigh,
				     const double *scale, const int *info);
// calculate hessenberg form
extern "C" void FORTRAN_NAME(zgehrd)(const int* dimensionN, const int *iLow, const int *iHigh,
				     const doublecomplex* matrixA, const int* leadingDimensionA,
				     const doublecomplex* tau, const doublecomplex* complexWork,
				     const int *lComplexWork, const int *info);
// extract eigenvalue
extern "C" void FORTRAN_NAME(zhseqr)(const char* jobZ, const char* compZ, const int* dimensionN,
				     const int *iLow, const int *iHigh, const doublecomplex* matrixH,
				     const int* leadingDimensionH, const doublecomplex* eigenValues,
				     const doublecomplex* schurZ,  const int* leadingDimensionZ,
				     const doublecomplex* complexWork, const int *lComplexWork, const int *info);
// extract eigenvectors
extern "C" void FORTRAN_NAME(zhsein)(const char* side, const char* eigsrc, const char* initV, const int *select,
				     const int* dimensionN, const doublecomplex* matrixH,
				     const int* leadingDimensionH, const doublecomplex* eigenValues,
				     const doublecomplex* matrixVL, const int* leadingDimensionVL,
				     const doublecomplex* matrixVR, const int* leadingDimensionVR,
				     const int* columnsVecRL, const int* columnsVecRLrequired,
				     const doublecomplex* complexWork, const double * realWork,
				     const int * columnIFailL, const int * columnIFailR, const int *info);
// extract eigenvectors and optionally use transformation matrix to get eigenvectors of initial matrix
extern "C" void FORTRAN_NAME(ztrevc)(const char* side,const char* howMNY, const int *select,
				     const int* dimensionN, const doublecomplex* matrixT,
				     const int* leadingDimensionT, const doublecomplex* matrixVL,
				     const int* leadingDimensionVL, const doublecomplex* matrixVR,
				     const int* leadingDimensionVR, const int* columnsVecRL,
				     const int* columnsVecRLrequired, const doublecomplex* complexWork,
				     const double * realWork, const int *info);
// extract unitary matrix from short for returned by zgehrd
extern "C" void FORTRAN_NAME(zunghr)( const int* dimensionN, const int *iLow, const int *iHigh,
				      const doublecomplex *matrixA, const int* leadingDimensionA,
				      const doublecomplex* tau, const doublecomplex* complexWork,
				      const int *lComplexWork, const int *info);
// multiply other matrix with unitary matrix from short for returned by zgehrd
extern "C" void FORTRAN_NAME(zunmhr)(const char* side, const char* trans, const int *numRows, const int *numCols,
				     const int *iLow, const int *iHigh, const doublecomplex *matrixA,
				     const int* leadingDimensionA, const doublecomplex* tau,
				     const doublecomplex *outputMatrixC,  const int* leadingDimensionC,
				     const doublecomplex* complexWork, const int *lComplexWork, const int *info);
// after all this pain: the one function which does it all...
extern "C" void FORTRAN_NAME(zgeevx)(const char* balanc,const char* jobVL,const char* jobVR, const char* sense,
				     const int* dimensionN, const doublecomplex *matrixA,
				     const int* leadingDimensionA, const doublecomplex* eigenValues,
				     const doublecomplex* matrixVL, const int* leadingDimensionVL,
				     const doublecomplex* matrixVR, const int* leadingDimensionVR,
				     const int *iLow, const int *iHigh, const double *scale,
				     const double *absNorm, const double *reciCondVal, const double *reciCondVec, 
				     const doublecomplex* complexWork, const int *lComplexWork,
				     const double * realWork, const int *info);

#endif



// default constructor
//

ComplexMatrix::ComplexMatrix() 
{
  this->Columns = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;  
  this->MatrixType = Matrix::ComplexElements;
}

// constructor for an empty matrix
//
// nbrRow = number of rows
// nbrColumn = number of columns
// zero = tue if matrix elements have to be set to zero

ComplexMatrix::ComplexMatrix(int nbrRow, int nbrColumn, bool zero)
{
  this->Flag.Initialize();
  this->NbrColumn = nbrColumn;
  this->NbrRow = nbrRow;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Columns = new ComplexVector [this->NbrColumn];
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] = ComplexVector (this->NbrRow, zero);
  this->MatrixType = Matrix::ComplexElements;
}

// constructor from matrix elements (without duplicating datas)
//
// columns = pointer an array of vector
// nbrColumn = number of columns

ComplexMatrix::ComplexMatrix(ComplexVector* columns, int nbrColumn) 
{
  this->Columns = columns;
  this->Flag.Initialize();
  this->NbrRow = columns[0].GetVectorDimension();
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;  
  this->MatrixType = Matrix::ComplexElements;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

ComplexMatrix::ComplexMatrix(const ComplexMatrix& M) 
{
  this->Columns = M.Columns;
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;  
  this->MatrixType = Matrix::ComplexElements;
}

// copy constructor (duplicating all datas)
//
// M = matrix to copy

ComplexMatrix::ComplexMatrix(Matrix& M)
{
  if ((M.GetNbrRow() == 0) || (M.GetNbrColumn() == 0))
    {
      this->Columns = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::ComplexElements;
    }
  else
    {
      this->Flag.Initialize();
      this->NbrColumn = M.GetNbrColumn();
      this->NbrRow = M.GetNbrRow();
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;
      this->Columns = new ComplexVector [this->NbrColumn];
      Complex Tmp;
      for (int i = 0; i < this->NbrColumn; i++)
	{
	  this->Columns[i] = ComplexVector (this->NbrRow);
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      M.GetMatrixElement(j, i, Tmp);
	      this->Columns[i].Re(j) = Tmp.Re;
	      this->Columns[i].Im(j) = Tmp.Im;
	    }
	}
      this->MatrixType = Matrix::ComplexElements;
    }
}

// destructor
//

ComplexMatrix::~ComplexMatrix() 
{
  if ((this->Columns != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->Columns;
      }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

ComplexMatrix& ComplexMatrix::operator = (const ComplexMatrix& M) 
{
  if ((this->Columns != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->Columns;
      }
  this->Columns = M.Columns;
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;  
  this->MatrixType = Matrix::ComplexElements;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* ComplexMatrix::Clone ()
{
  return ((Matrix*) new ComplexMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  this->Columns[j].Components[i].Re = x;
  this->Columns[j].Components[i].Im = 0.0;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  this->Columns[j].Components[i].Re = x.Re;
  this->Columns[j].Components[i].Im = x.Im;
}

// set a matrix element
//
// i = line position
// j = column position
// real = new real value for matrix element
// imag = new imaginary value for matrix element
void ComplexMatrix::SetMatrixElement(int i, int j, double real, double imag)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  this->Columns[j].Components[i].Re = real;
  this->Columns[j].Components[i].Im = imag;
}


// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void ComplexMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  this->Columns[j].Components[i].Re += x;
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void ComplexMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  this->Columns[j].Components[i].Re += x.Re;
  this->Columns[j].Components[i].Im += x.Im;
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexMatrix::Resize (int nbrRow, int nbrColumn)
{
  if (this->NbrRow != nbrRow)
    {
      for (int i = 0; i < this->NbrColumn; i++)
	this->Columns[i].Resize(nbrRow);
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
      ComplexVector* Tmp = new ComplexVector[nbrColumn];
      for (int i = 0; i < this->NbrColumn; i++)
	Tmp[i] = this->Columns[i];      
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	Tmp[i] = ComplexVector(nbrRow);
      if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
	{
	  delete[] this->Columns;
	}
      this->Columns = Tmp;
      this->TrueNbrColumn = nbrColumn;
      this->NbrColumn = nbrColumn;
      this->Flag = GarbageFlag();
      this->Flag.Initialize();
    }
  return;
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
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
      ComplexVector* Tmp = new ComplexVector[nbrColumn];
      for (int i = 0; i < this->NbrColumn; i++)
	Tmp[i] = this->Columns[i];      
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	Tmp[i] = ComplexVector(nbrRow, true);
      if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
	{
	  delete[] this->Columns;
	}
      this->Columns = Tmp;
      this->Flag = GarbageFlag();
      this->Flag.Initialize();
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

ComplexMatrix operator + (const ComplexMatrix& M1, const ComplexMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; ++j)
	{
	  TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re + M2.Columns[i].Components[j].Re;
	  TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im + M2.Columns[i].Components[j].Im;
	}
    }
  return ComplexMatrix(TmpColumns, M1.NbrColumn);
}

// add two matrices where the left one is a real tridiagonal symmetric matrix
//
// M1 = left matrix
// M2 = right matrix
// return value = sum of the two matrices

ComplexMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const ComplexMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      int j = 0;
      for (; j < (i - 1); ++j)
	{
	  TmpColumns[i].Components[j].Re = M2.Columns[i].Components[j].Re;
	  TmpColumns[i].Components[j].Im = M2.Columns[i].Components[j].Im;
	}
      if (i > 0)
	{
	  TmpColumns[i].Components[j].Re = M1.UpperDiagonalElements[i - 1] + M2.Columns[i].Components[j].Re;
	  TmpColumns[i].Components[j].Im = M2.Columns[i].Components[j].Im;
	  ++j;
	}
      TmpColumns[i].Components[j].Re = M1.DiagonalElements[i] + M2.Columns[i].Components[j].Re;
      TmpColumns[i].Components[j].Im = M2.Columns[i].Components[j].Im;
      ++j;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].Components[j].Re = M1.UpperDiagonalElements[i + 1] + M2.Columns[i].Components[j].Re;
	  TmpColumns[i].Components[j].Im = M2.Columns[i].Components[j].Im;
	  ++j;
	}
      ++j;
      for (; j < M1.NbrColumn; ++j)
	{
	  TmpColumns[i].Components[j].Re = M2.Columns[i].Components[j].Re;	
	  TmpColumns[i].Components[j].Im = M2.Columns[i].Components[j].Im;	
	}
    }
  return ComplexMatrix(TmpColumns, M1.NbrColumn);
}

// add two matrices where the right one is a real tridiagonal symmetric matrix
//
// M1 = left matrix
// M2 = right matrix
// return value = sum of the two matrices

ComplexMatrix operator + (const ComplexMatrix& M1, 
			  const RealTriDiagonalSymmetricMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      int j = 0;
      for (; j < (i - 1); ++j)
	{
	  TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re;
	  TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im;
	}
      if (i > 0)
	{
	  TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re + M2.UpperDiagonalElements[i - 1];
	  TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im;
	  ++j;
	}
      TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re + M2.DiagonalElements[i];
      TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im;
      ++j;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re + M2.UpperDiagonalElements[i + 1];
	  TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im;
	  ++j;
	}
      ++j;
      for (; j < M1.NbrColumn; ++j)
	{
	  TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re;	
	  TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im;	
	}
    }
  return ComplexMatrix(TmpColumns, M1.NbrColumn);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexMatrix operator - (const ComplexMatrix& M1, const ComplexMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; ++j)
	{
	  TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re - M2.Columns[i].Components[j].Re;
	  TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im - M2.Columns[i].Components[j].Im;	  
	}
    }
  return ComplexMatrix(TmpColumns, M1.NbrColumn);
}

// substract two matrices where the left one is a real tridiagonal symmetric matrix
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, const ComplexMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      int j = 0;
      for (; j < (i - 1); ++j)
	{
	  TmpColumns[i].Components[j].Re = -M2.Columns[i].Components[j].Re;
	  TmpColumns[i].Components[j].Im = -M2.Columns[i].Components[j].Im;
	}
      if (i > 0)
	{
	  TmpColumns[i].Components[j].Re = M1.UpperDiagonalElements[i - 1] - M2.Columns[i].Components[j].Re;
	  TmpColumns[i].Components[j].Im = -M2.Columns[i].Components[j].Im;
	  ++j;
	}
      TmpColumns[i].Components[j].Re = M1.DiagonalElements[i] - M2.Columns[i].Components[j].Re;
      TmpColumns[i].Components[j].Im = -M2.Columns[i].Components[j].Im;
      ++j;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].Components[j].Re = M1.UpperDiagonalElements[i + 1] - M2.Columns[i].Components[j].Re;
	  TmpColumns[i].Components[j].Im = -M2.Columns[i].Components[j].Im;
	  ++j;
	}
      ++j;
      for (; j < M1.NbrColumn; ++j)
	{
	  TmpColumns[i].Components[j].Re = -M2.Columns[i].Components[j].Re;	
	  TmpColumns[i].Components[j].Im = -M2.Columns[i].Components[j].Im;	
	}
    }
  return ComplexMatrix(TmpColumns, M1.NbrColumn);
}

// substract two matrices where the right one is a real tridiagonal symmetric matrix
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexMatrix operator - (const ComplexMatrix& M1, 
			  const RealTriDiagonalSymmetricMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      int j = 0;
      for (; j < (i - 1); ++j)
	{
	  TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re;
	  TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im;
	}
      if (i > 0)
	{
	  TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re - M2.UpperDiagonalElements[i - 1];
	  TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im;
	  ++j;
	}
      TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re - M2.DiagonalElements[i];
      TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im;
      ++j;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re - M2.UpperDiagonalElements[i + 1];
	  TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im;
	  ++j;
	}
      ++j;
      for (; j < M1.NbrColumn; ++j)
	{
	  TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re;	
	  TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im;	
	}
    }
  return ComplexMatrix(TmpColumns, M1.NbrColumn);
}

// multiply two matrices
//
// M1 = first matrix
// M2 = matrix to multiply to M1
// return value = product of the two matrices

ComplexMatrix operator * (const ComplexMatrix& M1, const ComplexMatrix& M2)
{
  if (M1.NbrColumn != M2.NbrRow)
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M2.NbrColumn];
  for (int i = 0; i < M2.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; ++j)
	{
	  TmpColumns[i].Components[j].Re = 0.0;
	  TmpColumns[i].Components[j].Im = 0.0;
	  for (int k = 0; k < M2.NbrRow; ++k)	
	    {
	      TmpColumns[i].Components[j].Re += (M1.Columns[k].Components[j].Re * M2.Columns[i].Components[k].Re - 
						  M1.Columns[k].Components[j].Im * M2.Columns[i].Components[k].Im);
	      TmpColumns[i].Components[j].Im += (M1.Columns[k].Components[j].Re * M2.Columns[i].Components[k].Im + 
						       M1.Columns[k].Components[j].Im * M2.Columns[i].Components[k].Re);
	    }
	}
    }
  return ComplexMatrix(TmpColumns, M2.NbrColumn);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexMatrix operator * (const ComplexMatrix& M, double x) 
{
  ComplexVector* TmpColumns = new ComplexVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M.NbrRow);
      for (int j = 0; j < M.NbrRow; ++j)
	{
	  TmpColumns[i].Components[j].Re = M.Columns[i].Components[j].Re * x;
	  TmpColumns[i].Components[j].Im = M.Columns[i].Components[j].Im * x;
	}
    }
  return ComplexMatrix(TmpColumns, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexMatrix operator * (double x, const ComplexMatrix& M) 
{
  ComplexVector* TmpColumns = new ComplexVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M.NbrRow);
      for (int j = 0; j < M.NbrRow; ++j)
	{
	  TmpColumns[i].Components[j].Re = M.Columns[i].Components[j].Re * x;
	  TmpColumns[i].Components[j].Im = M.Columns[i].Components[j].Im * x;
	}
    }
  return ComplexMatrix(TmpColumns, M.NbrRow);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

ComplexMatrix operator / (const ComplexMatrix& M, double x) 
{
  ComplexVector* TmpColumns = new ComplexVector [M.NbrColumn];
  x = 1.0 / x;
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = ComplexVector (M.NbrRow);
      for (int j = 0; j < M.NbrRow; ++j)
	{
	  TmpColumns[i].Components[j].Re = M.Columns[i].Components[j].Re * x;
	  TmpColumns[i].Components[j].Im = M.Columns[i].Components[j].Im * x;
	}
    }
  return ComplexMatrix(TmpColumns, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::operator += (const ComplexMatrix& M) 
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

ComplexMatrix& ComplexMatrix::operator += (const RealTriDiagonalSymmetricMatrix& M) 
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow) || (this->Columns == 0))
    return *this;  
  this->Columns[0].Components[0].Re += M.DiagonalElements[0];
  for (int i = 1; i < this->NbrColumn; i++)
    {
      this->Columns[i].Components[i].Re += M.DiagonalElements[i];
      this->Columns[i].Components[i - 1].Re += M.UpperDiagonalElements[i - 1];
      this->Columns[i - 1].Components[i].Re += M.UpperDiagonalElements[i - 1];
    }
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::operator -= (const ComplexMatrix& M) 
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

ComplexMatrix& ComplexMatrix::operator -= (const RealTriDiagonalSymmetricMatrix& M) 
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow) || (this->Columns == 0))
    return *this;  
  this->Columns[0].Components[0].Re -= M.DiagonalElements[0];
  for (int i = 1; i < this->NbrColumn; i++)
    {
      this->Columns[i].Components[i].Re -= M.DiagonalElements[i];
      this->Columns[i].Components[i - 1].Re -= M.UpperDiagonalElements[i - 1];
      this->Columns[i - 1].Components[i].Re -= M.UpperDiagonalElements[i - 1];
    }
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::operator *= (double x) 
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] *= x;
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::operator /= (double x)
{
  x = 1.0 / x;;
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] *= x;
  return *this;
}

// normalize matrix column vectors
//
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::NormalizeColumns ()
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i].Normalize();
  return *this;
}

// orthonormalize matrix column vectors
//
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::OrthoNormalizeColumns ()
{
  Complex* tmp = new Complex [this->NbrColumn];
  for (int i = 0; i < this->NbrColumn; i++)
    {
      for (int j = 0; j < i; j++)
	{
	  tmp[j] = this->Columns[i] * this->Columns[j];
	}
      for (int j = 0; j < i; j++)
	{
	  this->Columns[i].AddLinearCombination(tmp[j], this->Columns[j]);
	}
      this->Columns[i].Normalize();
    }      
  delete[] tmp;
  return *this;
}


// get adjoint (hermitian conjugate) matrix 
//
// return value = reference on modified matrix

ComplexMatrix ComplexMatrix::GetAdjoint()
{
  ComplexMatrix rst(this->NbrColumn, this->NbrRow);
  for (int j=0; j<this->NbrColumn; ++j)
    for (int i=0; i<this->NbrRow; ++i)
      rst.SetMatrixElement(j,i,Conj(this->Columns[j][i]));
  return rst;
}

// evaluate matrix determinant (screwing up matrix elements)
//
// return value = matrix determinant 

Complex ComplexMatrix::Determinant () 
{
#ifdef __LAPACKONLY__
  return this->LapackDeterminant();
#endif
  if (this->NbrColumn != this->NbrRow)
    return 0.0;
  Complex TmpDet (1.0);
  int ReducedNbrRow = this->NbrRow - 1;
  Complex Pivot;
  Complex Factor;
  int PivotPos = 0;
  double PivotNorm;
  for (int k = 0; k < ReducedNbrRow; ++k)
    {
      Pivot.Re = this->Columns[k].Re(k);
      Pivot.Im = this->Columns[k].Im(k);
      PivotNorm = (Pivot.Re * Pivot.Re) + (Pivot.Im * Pivot.Im);
      PivotPos = k + 1;
      while ((PivotPos < this->NbrRow) && 
	     (((this->Columns[PivotPos].Re(k) * this->Columns[PivotPos].Re(k)) + (this->Columns[PivotPos].Im(k) * this->Columns[PivotPos].Im(k)))< PivotNorm))
	{
	  ++PivotPos;
	}
      if (PivotPos == this->NbrRow)
	{
	  if (PivotNorm == 0.0)
	    return Complex(0.0);
	}
      else
	{
	  Pivot.Re = this->Columns[PivotPos].Re(k);
	  Pivot.Im = this->Columns[PivotPos].Im(k);
	  ComplexVector TmpColumn3(this->Columns[k]);
	  this->Columns[k] = this->Columns[PivotPos];
	  this->Columns[PivotPos] = TmpColumn3;
	  TmpDet *= -1.0;
	}
      if (PivotNorm == 0.0)
	return Complex(0.0);
      TmpDet *= Pivot;
      Pivot = 1.0 / Pivot;       
      for (int i = k + 1; i < this->NbrRow; ++i)
	{
	  ComplexVector& TmpColumn = this->Columns[i];
	  ComplexVector& TmpColumn2 = this->Columns[k];
	  Factor.Re = ((Pivot.Re * TmpColumn.Re(k)) - (Pivot.Im * TmpColumn.Im(k)));
	  Factor.Im = ((Pivot.Im * TmpColumn.Re(k)) + (Pivot.Re * TmpColumn.Im(k)));
	  for (int j = k + 1; j < this->NbrRow; ++j)
	    {
	      TmpColumn.Re(j) -= ((TmpColumn2.Re(j) * Factor.Re) - (TmpColumn2.Im(j) * Factor.Im)); 
	      TmpColumn.Im(j) -= ((TmpColumn2.Re(j) * Factor.Im) + (TmpColumn2.Im(j) * Factor.Re)); 
	    }
	}
    } 
  Pivot.Re = this->Columns[ReducedNbrRow].Re(ReducedNbrRow);
  Pivot.Im = this->Columns[ReducedNbrRow].Im(ReducedNbrRow);
  TmpDet *= Pivot;
  return TmpDet;
}


// evaluate permanent associated to the (square) matrix using Ryser algorithm
//
// return value = permanent associated to the matrix
                                                                                                                                          
Complex ComplexMatrix::Permanent()
{
  if (this->NbrColumn != this->NbrRow)
    return 0.0;
  Complex Perm;
  double Sign = 1.0;
  if ((this->NbrColumn & 1) == 0)
    Sign = -1.0;
  Complex* Tmp = new Complex [this->NbrColumn];
  Complex Tmp2;
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
	  ComplexVector& TmpColumn = this->Columns[Index];
	  for (int i = 0; i < this->NbrColumn; ++i)
	    {
	      Tmp[i].Re -= TmpColumn.Components[i].Re;
	      Tmp[i].Im -= TmpColumn.Components[i].Im;
	    }
	}
      else
	{
	  Index = 0;
	  while (ChangedBit != 1)
	    {
	      ChangedBit >>= 1;
	      ++Index;
	    }
	  ComplexVector& TmpColumn = this->Columns[Index];
	  for (int i = 0; i < this->NbrColumn; ++i)
	    {
	      Tmp[i].Re += TmpColumn.Components[i].Re;
	      Tmp[i].Im += TmpColumn.Components[i].Im;
	    }
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

// evaluate minor develomment of permanent associated to the (square) matrix using Ryser algorithm
//
// column = index of the column from which permnanent will developped
// minors = reference on an array where minors will be stored
                                                                                                                                          
void ComplexMatrix::PermanentMinorDevelopment(int column, Complex*& minors)
{
  if (this->NbrColumn != this->NbrRow)
    return;
  int ReducedNbrColumn = this->NbrColumn - 1;
  Complex* Tmp = new Complex  [ReducedNbrColumn];
  Complex Tmp2;
  int Lim = 1 << ReducedNbrColumn;
  for (int l = 0; l < this->NbrColumn; ++l)
    {
      minors[l].Re = 0.0;
      minors[l].Im = 0.0;
      double Sign = 1.0;
      if ((ReducedNbrColumn & 1) == 0)
	Sign = -1.0;
      for (int i = 0; i < ReducedNbrColumn; ++i)
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
	      if (Index >= column)
		++Index;
	      int i = 0;
	      ComplexVector& TmpColumn = this->Columns[Index];
	      for (; i < l; ++i)
		{
		  Tmp[i].Re -= TmpColumn.Components[i].Re;
		  Tmp[i].Im -= TmpColumn.Components[i].Im;
		}
	      for (; i < ReducedNbrColumn; ++i)
		{
		  Tmp[i].Re -= TmpColumn.Components[i + 1].Re;
		  Tmp[i].Im -= TmpColumn.Components[i + 1].Im;
		}
	    }
	  else
	    {
	      Index = 0;
	      while (ChangedBit != 1)
		{
		  ChangedBit >>= 1;
		  ++Index;
		}
	      if (Index >= column)
		++Index;
	      int i = 0;
	      ComplexVector& TmpColumn = this->Columns[Index];
	      for (; i < l; ++i)
		{
		  Tmp[i].Re += TmpColumn.Components[i].Re;
		  Tmp[i].Im += TmpColumn.Components[i].Im;
		}
	      for (; i < ReducedNbrColumn; ++i)
		{
		  Tmp[i].Re += TmpColumn.Components[i + 1].Re;
		  Tmp[i].Im += TmpColumn.Components[i + 1].Im;
		}
	    }
 	  Tmp2 = Tmp[0];
	  for (int i = 1; i < ReducedNbrColumn; ++i)
	    Tmp2 *= Tmp[i];
	  minors[l] += Sign * Tmp2;
	  Sign *= -1.0;
	}
    }
  delete[] Tmp;
}

// evaluate permanent associated to the (square) matrix using Ryser algorithm using precalculation array (faster)
//
// changeBit = array indicating which bit is changed at the i-th iteration of the Gray code
// changeBitSign = array with 0 if the changed bit is from 1 to 0, +1 either
// return value = permanent associated to the matrix
                                                                                                                                          
Complex ComplexMatrix::FastPermanent(int* changeBit, int* changeBitSign)
{
  Complex Perm;
  double Sign = 1.0;
  if ((this->NbrColumn & 1) == 0)
    Sign = -1.0;
  Complex* Tmp = new Complex [this->NbrColumn];
  Complex Tmp2;
  int Lim = 1 << this->NbrColumn;
  for (int i = 0; i < this->NbrColumn; ++i)
    Tmp[i] = 0.0;
  for (int k = 1; k < Lim; ++k)
    {
      ComplexVector& TmpColumn = this->Columns[changeBit[k]];
      if (changeBitSign[k] == 0)
	{
	  for (int i = 0; i < this->NbrColumn; ++i)
	    {
	      Tmp[i].Re -= TmpColumn.Components[i].Re;
	      Tmp[i].Im -= TmpColumn.Components[i].Im;	  
	    }
	}
      else
	{
	  for (int i = 0; i < this->NbrColumn; ++i)
	    {
	      Tmp[i].Re += TmpColumn.Components[i].Re;
	      Tmp[i].Im += TmpColumn.Components[i].Im;	  
	    }
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


// evaluate minor develomment of permanent associated to the (square) matrix using Ryser algorithm and precalculation array (faster)
//
// changeBit = array indicating which bit is changed at the i-th iteration of the Gray code
// changeBitSign = array with -1 if the changed bit is from 1 to 0, +1 either
// column = index of the column from which permnanent will developped
// minors = reference on an array where minors will be stored

void ComplexMatrix::FastPermanentMinorDevelopment(int* changeBit, int* changeBitSign, int column, Complex*& minors)
{
  int ReducedNbrColumn = this->NbrColumn - 1;
  Complex* Tmp = new Complex [ReducedNbrColumn];
  Complex Tmp2;
  ComplexVector* TmpColumn;
  int Lim = 1 << ReducedNbrColumn;
  for (int l = 0; l < this->NbrColumn; ++l)
    {
      minors[l].Re = 0.0;
      minors[l].Im = 0.0;
      double Sign = 1.0;
      if ((ReducedNbrColumn & 1) == 0)
	Sign = -1.0;
      for (int i = 0; i < ReducedNbrColumn; ++i)
	Tmp[i] = 0.0;
      for (int k = 1; k < Lim; ++k)
	{
	  if (changeBit[k] < column)
	    {
	      TmpColumn = &(this->Columns[changeBit[k]]);
	    }
	  else
	    {
	      TmpColumn = &(this->Columns[changeBit[k] + 1]);
	    }
	  if (changeBitSign[k] == 0)
	    {
	      int i = 0;
	      for (; i < l; ++i)
		{
		  Tmp[i].Re -= TmpColumn->Components[i].Re;
		  Tmp[i].Im -= TmpColumn->Components[i].Im;	  
		}
	      for (; i < ReducedNbrColumn; ++i)
		{
		  Tmp[i].Re -= TmpColumn->Components[i + 1].Re;
		  Tmp[i].Im -= TmpColumn->Components[i + 1].Im;
		}
	    }
	  else
	    {
	      int i = 0;
	      for (; i < l; ++i)
		{
		  Tmp[i].Re += TmpColumn->Components[i].Re;
		  Tmp[i].Im += TmpColumn->Components[i].Im;	  
		}
	      for (; i < ReducedNbrColumn; ++i)
		{
		  Tmp[i].Re += TmpColumn->Components[i + 1].Re;
		  Tmp[i].Im += TmpColumn->Components[i + 1].Im;	  
		}
	    }
	  Tmp2 = Tmp[0];
	  for (int i = 1; i < ReducedNbrColumn; ++i)
	    Tmp2 *= Tmp[i];
	  minors[l] += Sign * Tmp2;
	  Sign *= -1.0;
	}
    }
  delete[] Tmp;
}
  
// evaluate precalculation array  neede for the fast permanent calculation
//
// changeBit = reference on the array indicating which bit is changed at the i-th iteration of the Gray code
// changeBitSign = reference on array with 0 if the changed bit is from 1 to 0, +1 either
// minor = flag that indicated if precalculation will be used for minor development

void ComplexMatrix::EvaluateFastPermanentPrecalculationArray(int*& changeBit, int*& changeBitSign, bool minor)
{
  int Lim ;
  if (minor == false)
    {
      Lim = 1 << this->NbrColumn;
    }
  else
    {
      Lim = 1 << (this->NbrColumn - 1);
    }
  int GrayCode = 0;
  int ChangedBit;
  int Index;
  changeBit = new int [Lim];
  changeBitSign = new int [Lim];
  for (int k = 1; k < Lim; ++k)
    {
      ChangedBit = (k ^ (k >> 1)) ^ GrayCode;
      GrayCode = k ^ (k >> 1);
      if ((GrayCode & ChangedBit) == 0)
	{
	  changeBitSign[k] = 0;
	}
      else
	{
	  changeBitSign[k] = 1;
	}
      Index = 0;
      while (ChangedBit != 1)
	{
	  ChangedBit >>= 1;
	  ++Index;
	}
      changeBit[k] = Index;
   }
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const ComplexMatrix& P) 
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      for (int j = 0; j < (P.NbrColumn - 1); j ++)
	{
	  Str << P.Columns[j].Components[i].Re;      
	  if (P.Columns[j].Components[i].Im < 0.0)
	    Str << P.Columns[j].Components[i].Im << "i    ";
	  else
	    if (P.Columns[j].Components[i].Im != 0.0)
	      Str << "+" << P.Columns[j].Components[i].Im << "i    ";
	    else
	      Str << "    ";
	}
      Str << P.Columns[P.NbrColumn - 1].Components[i].Re;      
      if (P.Columns[P.NbrColumn - 1].Components[i].Im < 0.0)
	Str << P.Columns[P.NbrColumn - 1].Components[i].Im << "i";
      else
	if (P.Columns[P.NbrColumn - 1].Components[i].Im != 0.0)
	  Str << "+" << P.Columns[P.NbrColumn - 1].Components[i].Im << "i";
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

MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexMatrix& P) 
{
  Str << "{";
  for (int i = 0; i < (P.NbrRow - 1); ++i)
    {
      Str << "{";
      for (int j = 0; j < (P.NbrColumn - 1); ++j)
	{
	  Str << P.Columns[j].Components[i].Re;      
	  if (P.Columns[j].Components[i].Im < 0.0)
	    Str << P.Columns[j].Components[i].Im << "I,";
	  else
	    if (P.Columns[j].Components[i].Im != 0.0)
	      Str << "+" << P.Columns[j].Components[i].Im << "I,";
	    else
	      Str << ",";
	}
      Str << P.Columns[P.NbrColumn - 1].Components[i].Re;      
      if (P.Columns[P.NbrColumn - 1].Components[i].Im < 0.0)
	Str << P.Columns[P.NbrColumn - 1].Components[i].Im << "I";
      else
	if (P.Columns[P.NbrColumn - 1].Components[i].Im != 0.0)
	  Str << "+" << P.Columns[P.NbrColumn - 1].Components[i].Im << "I";
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (P.NbrColumn - 1); ++j)
    {
      Str << P.Columns[j].Components[P.NbrRow - 1].Re;
      if (P.Columns[j].Components[P.NbrRow - 1].Im < 0.0)
	Str << P.Columns[j].Components[P.NbrRow - 1].Im << "I,";
      else
	if (P.Columns[j].Components[P.NbrRow - 1].Im != 0.0)
	  Str << "+" << P.Columns[j].Components[P.NbrRow - 1].Im << "I,";
	else
	  Str << ",";
    }
  Str << P.Columns[P.NbrColumn - 1].Components[P.NbrRow - 1].Re;      
  if (P.Columns[P.NbrColumn - 1].Components[P.NbrRow - 1].Im < 0.0)
    Str << P.Columns[P.NbrColumn - 1].Components[P.NbrRow - 1].Im << "I";
  else
    if (P.Columns[P.NbrColumn - 1].Components[P.NbrRow - 1].Im != 0.0)
      Str << "+" << P.Columns[P.NbrColumn - 1].Components[P.NbrRow - 1].Im << "I";
  Str << "}}";
  return Str;
}

#endif



#ifdef __LAPACK__

// calculate a determinant using the LAPACK library (conserving current matrix)
//
Complex ComplexMatrix::LapackDeterminant ()
{
  if (this->NbrColumn != this->NbrRow)
    return 0.0;
  doublecomplex* TmpMatrix = new doublecomplex [this->NbrRow * this->NbrRow];
  Complex *TmpColumn;
  for (int j=0;j<NbrRow;++j)
    {
      TmpColumn=this->Columns[j].Components;
      for (int i=0; i<NbrRow;++i)
	{
	  TmpMatrix[i+j*NbrRow].r=TmpColumn[i].Re;
	  TmpMatrix[i+j*NbrRow].i=TmpColumn[i].Im;
	}
    }
  int Information = 0;
  int DimensionM=NbrRow;
  int *Permutation = new int[NbrRow];
  FORTRAN_NAME(zgetrf)(&DimensionM, &DimensionM, TmpMatrix, &DimensionM , Permutation, &Information);

  if (Information < 0)
    {
      cout << "Illegal argument " << -Information << " in LAPACK function call in ComplexMatrix.cc, line "<< __LINE__<<endl;
      exit(1);
    }

  int sign=0;
  Complex Result(1.0,0.0);
  
  for (int i=0; i<DimensionM; ++i)
    {
      if (Permutation[i]!=i+1) sign ^= 1;
      Result *= Complex(TmpMatrix[i+DimensionM*i].r,TmpMatrix[i+DimensionM*i].i);
    }
  if (sign & 1)
      Result*=-1.0;

  delete [] TmpMatrix;
  delete [] Permutation;
  
  return Result;
}

#endif



// Diagonalize an hermitian matrix (modifying current matrix)
//
// M = reference on real diagonal matrix where result has to be stored
// return value = reference on real tridiagonal symmetric matrix

ComplexDiagonalMatrix& ComplexMatrix::Diagonalize (ComplexDiagonalMatrix& M)
{
#ifdef __LAPACK__
  return this->LapackDiagonalize(M);
#else
  cout << "ComplexMatrix::Diagonalize needs to be implemented!"<<endl;
  exit(1);
  return M;
#endif
}

// Diagonalize an hermitian matrix and evaluate transformation matrix (modifying current matrix)
//
// M = reference on real diagonal matrix where result has to be stored
// Q = matrix where transformation matrix has to be stored
// err = absolute error on matrix element
// maxIter = maximum number of iteration to fund an eigenvalue
// return value = reference on real tridiagonal symmetric matrix

ComplexDiagonalMatrix& ComplexMatrix::Diagonalize (ComplexDiagonalMatrix& M, ComplexMatrix& Q)
{
#ifdef __LAPACK__
  return this->LapackDiagonalize(M, Q);
#else
  cout << "ComplexMatrix::Diagonalize needs to be implemented!"<<endl;
  exit(1);
  return M;
#endif
}



#ifdef __LAPACK__
  
// Diagonalize a complex skew symmetric matrix using the LAPACK library (modifying current matrix)
//
// M = reference on real diagonal matrix of eigenvalues
// err = absolute error on matrix element
// maxIter = maximum number of iteration to fund an eigenvalue
// return value = reference on real matrix consisting of eigenvalues

ComplexDiagonalMatrix& ComplexMatrix::LapackDiagonalize (ComplexDiagonalMatrix& M)
{
  if (this->NbrColumn != this->NbrRow)
    return M;
  if (M.GetNbrColumn() != this->NbrColumn)
    M.Resize(this->NbrColumn, this->NbrColumn);
  doublecomplex* TmpMatrix = new doublecomplex [this->NbrRow * this->NbrRow];
  doublecomplex* Tau = new doublecomplex [this->NbrRow];
  doublecomplex* Eigenvalues = new doublecomplex [this->NbrRow];
  //double *Scale = new double[this->NbrRow];
  Complex *TmpColumn;
  for (int j=0;j<NbrRow;++j)
    {
      TmpColumn=this->Columns[j].Components;
      for (int i=0; i<NbrRow;++i)
	{
	  TmpMatrix[i+j*NbrRow].r=TmpColumn[i].Re;
	  TmpMatrix[i+j*NbrRow].i=TmpColumn[i].Im;
	}
    }
  int Information = 0;
  int DimensionM=NbrRow;
  int iLow=1;
  int iHigh=DimensionM;
  doublecomplex *complexWork=new doublecomplex[1];
  int lComplexWork=-1;
  // balancing omitted
  // workspace query
  FORTRAN_NAME(zgehrd)(&DimensionM, &iLow, &iHigh, TmpMatrix, &DimensionM, Tau,
		       complexWork, &lComplexWork, &Information);
  lComplexWork=(int)(complexWork[0].r);
  delete [] complexWork;
  complexWork=new doublecomplex[lComplexWork];
  FORTRAN_NAME(zgehrd)(&DimensionM, &iLow, &iHigh, TmpMatrix, &DimensionM, Tau,
		       complexWork, &lComplexWork, &Information);  
  if (Information < 0)
    {
      cout << "Illegal argument " << -Information << " in LAPACK function call to zgehrd in ComplexMatrix.cc in LapackDiagonalize, line "<< __LINE__<<endl;
      exit(1);
    }
  doublecomplex* SchurZ=NULL;
  int dimZ=1;  
  lComplexWork=-1;
  // workspace query
  FORTRAN_NAME(zhseqr)("E", "N", &DimensionM, &iLow, &iHigh, TmpMatrix, &DimensionM,
		       Eigenvalues, SchurZ, &dimZ, complexWork, &lComplexWork, &Information);
  lComplexWork=(int)(complexWork[0].r);
  delete [] complexWork;
  complexWork=new doublecomplex[lComplexWork];
  // perform eigenvalue calculation
  FORTRAN_NAME(zhseqr)("E", "N", &DimensionM, &iLow, &iHigh, TmpMatrix, &DimensionM,
		       Eigenvalues, SchurZ, &dimZ, complexWork, &lComplexWork, &Information);
  if (Information < 0)
    {
      cout << "Illegal argument " << -Information << " in LAPACK function call to zhseqr in ComplexMatrix.cc in LapackDiagonalize, line "<< __LINE__<<endl;
      exit(1);
    }
  if (Information > 0)
    {
      cout << "Only part of eigenvalues calculated in ComplexMatrix::LapackDiagonalize, line "<< __LINE__<<endl;
    }

  for (int i=0; i<DimensionM; ++i)
    M.SetMatrixElement(i,i,Eigenvalues[i].r, Eigenvalues[i].i);
  
  delete [] TmpMatrix;
  delete [] Tau;
  delete [] Eigenvalues;
  //delete [] Scale;
  delete [] complexWork;
  
  return M;
}

// Diagonalize a complex skew symmetric matrix and evaluate transformation matrix using the LAPACK library (modifying current matrix)
//
// M = reference on real diagonal matrix of eigenvalues
// Q = matrix where transformation matrix has to be stored
// err = absolute error on matrix element
// maxIter = maximum number of iteration to fund an eigenvalue
// return value = reference on real matrix consisting of eigenvalues

ComplexDiagonalMatrix& ComplexMatrix::LapackDiagonalize (ComplexDiagonalMatrix& M, ComplexMatrix& Q)
{
  if (M.GetNbrColumn() != this->NbrColumn)
    M.Resize(this->NbrColumn, this->NbrColumn);
  if (Q.GetNbrColumn() != this->NbrColumn)
    Q.Resize(this->NbrColumn, this->NbrColumn);
  doublecomplex* matrixA = new doublecomplex [this->NbrRow * this->NbrRow];
  double* scale = new double [this->NbrRow];
  double *reciCondVal = new double [this->NbrRow];
  double *reciCondVec = new double [this->NbrRow];
  doublecomplex* TransQ = new doublecomplex [this->NbrRow * this->NbrRow];
  doublecomplex* Eigenvalues = new doublecomplex [this->NbrRow];
  doublecomplex* EigenvectorsR = new doublecomplex [this->NbrRow * this->NbrRow];
  doublecomplex* EigenvectorsL = NULL;
  Complex *TmpColumn;
  for (int j=0;j<NbrRow;++j)
    {
      TmpColumn=this->Columns[j].Components;
      for (int i=0; i<NbrRow;++i)
	{
	  matrixA[i+j*NbrRow].r=TmpColumn[i].Re;
	  matrixA[i+j*NbrRow].i=TmpColumn[i].Im;
	}
    }  
  int Information = 0;
  int Dim = this->NbrRow;
  int iLow=1, iHigh=this->NbrRow;
  doublecomplex *complexWork = new doublecomplex[1];
  double *realWork = new double[2*NbrRow];
  double absNorm;
  int lComplexWork=-1;
  // workspace query
  FORTRAN_NAME(zgeevx)("Balance","No left eigenvectors","Vectors on right", "No condition numbers",
		       &Dim, matrixA, &Dim, Eigenvalues, EigenvectorsL, &Dim, EigenvectorsR, &Dim,
		       &iLow, &iHigh, scale, &absNorm, reciCondVal, reciCondVec, 
		       complexWork, &lComplexWork, realWork, &Information);  
  lComplexWork=(int)(complexWork[0].r);
  delete [] complexWork;
  complexWork = new doublecomplex[lComplexWork];
  FORTRAN_NAME(zgeevx)("Balance","No left eigenvectors","Vectors on right", "No condition numbers",
		       &Dim, matrixA, &Dim, Eigenvalues, EigenvectorsL, &Dim, EigenvectorsR, &Dim,
		       &iLow, &iHigh, scale, &absNorm, reciCondVal, reciCondVec, 
		       complexWork, &lComplexWork, realWork, &Information);
  if (Information < 0)
    {
      cout << "Illegal argument " << -Information << " in LAPACK function call to zgeevx in ComplexMatrix.cc in LapackDiagonalize, line "<< __LINE__<<endl;
      exit(1);
    }
  else if (Information > 0)
    {
      cout << "Attention: Only part of eigenvalues converged in LAPACK function call to zgeevx in ComplexMatrix.cc in LapackDiagonalize, line "<< __LINE__<<endl;
    }
    
  // recover values of eigenvalues and eigenvectors
  for (int j=0;j<NbrRow;++j)
    {
      M.SetMatrixElement(j, j, Eigenvalues[j].r, Eigenvalues[j].i);
      for (int i=0; i<NbrRow;++i)
	Q.SetMatrixElement(i, j, EigenvectorsR[i+j*NbrRow].r, EigenvectorsR[i+j*NbrRow].i);
    }
  
  delete [] matrixA;
  delete [] Eigenvalues;
  delete [] EigenvectorsR;
  delete [] TransQ;
  delete [] reciCondVal;
  delete [] reciCondVec;
  delete [] scale;
  delete [] complexWork;
  delete [] realWork;
    
  return M;
}

// reduce a complex matrix to its Schur form S
//
// M = reference on real diagonal matrix of eigenvalues
// Q = matrix where transformation matrix has to be stored
// S = matrix where Schur form of matrix has to be stored
// return value = reference on real matrix consisting of eigenvalues

ComplexDiagonalMatrix& ComplexMatrix::LapackSchurForm (ComplexDiagonalMatrix& M, ComplexMatrix& Q, ComplexMatrix &S)
{
  if (M.GetNbrColumn() != this->NbrColumn)
    M.Resize(this->NbrColumn, this->NbrColumn);
  if (Q.GetNbrColumn() != this->NbrColumn)
    Q.Resize(this->NbrColumn, this->NbrColumn);
  if (S.GetNbrColumn() != this->NbrColumn)
    S.Resize(this->NbrColumn, this->NbrColumn);
  doublecomplex* TmpMatrix = new doublecomplex [this->NbrRow * this->NbrRow];
  doublecomplex* TransQ = new doublecomplex [this->NbrRow * this->NbrRow];
  doublecomplex* Tau = new doublecomplex [this->NbrRow];
  doublecomplex* Eigenvalues = new doublecomplex [this->NbrRow];
  Complex *TmpColumn;
  for (int j=0;j<NbrRow;++j)
    {
      TmpColumn=this->Columns[j].Components;
      for (int i=0; i<NbrRow;++i)
	{
	  TmpMatrix[i+j*NbrRow].r=TmpColumn[i].Re;
	  TmpMatrix[i+j*NbrRow].i=TmpColumn[i].Im;
	}
    }
  int Information = 0;
  int DimensionM=NbrRow;
  int iLow=1;
  int iHigh=DimensionM;
  int dimZ=DimensionM;
  doublecomplex *complexWork=new doublecomplex[this->NbrRow * this->NbrRow];
  int lComplexWork=this->NbrRow * this->NbrRow;
  // balancing omitted
  // workspace query

//   FORTRAN_NAME(zgehrd)(&DimensionM, &iLow, &iHigh, TmpMatrix, &DimensionM, Tau,
// 		       complexWork, &lComplexWork, &Information);
//   lComplexWork=(int)(complexWork[0].r);
//   delete [] complexWork;
//   complexWork=new doublecomplex[lComplexWork];
  
  // reduce to hessenberg form

  FORTRAN_NAME(zgehrd)(&DimensionM, &iLow, &iHigh, TmpMatrix, &DimensionM, Tau,
		       complexWork, &lComplexWork, &Information);  
  if (Information < 0)
    {
      cout << "Illegal argument " << -Information << " in LAPACK function call to zgehrd in ComplexMatrix.cc in LapackDiagonalize, line "<< __LINE__<<endl;
      exit(1);
    }
  // store part of TmpMatrix associated with deflectors
  for (int j=0;j<NbrRow;++j)
    for (int i=0; i<NbrRow;++i)
      {
	TransQ[i+j*NbrRow].r = TmpMatrix[i+j*NbrRow].r;
	TransQ[i+j*NbrRow].i = TmpMatrix[i+j*NbrRow].i;
      }

  // extract unitary matrix from short for returned by zgehrd
  FORTRAN_NAME(zunghr)( &DimensionM, &iLow, &iHigh, TransQ, &DimensionM, Tau, complexWork,
			&lComplexWork, &Information);
  if (Information < 0)
    {
      cout << "Illegal argument " << -Information << " in LAPACK function call to zunghr in ComplexMatrix.cc in LapackDiagonalize, line "<< __LINE__<<endl;
      exit(1);
    }  
  //  lComplexWork=-1;
  // workspace query

//   FORTRAN_NAME(zhseqr)("S", "I", &DimensionM, &iLow, &iHigh, TmpMatrix, &DimensionM,
// 		       Eigenvalues, SchurZ, &dimZ, complexWork, &lComplexWork, &Information);
//   lComplexWork=(int)(complexWork[0].r);
//   delete [] complexWork;
//   complexWork=new doublecomplex[lComplexWork];
  
  // perform eigenvalue and schur decomposition calculation
  FORTRAN_NAME(zhseqr)("S", "V", &DimensionM, &iLow, &iHigh, TmpMatrix, &DimensionM,
		       Eigenvalues, TransQ, &dimZ, complexWork, &lComplexWork, &Information);    
  if (Information < 0)
    {
      cout << "Illegal argument " << -Information << " in LAPACK function call to zhseqr in ComplexMatrix.cc in LapackDiagonalize, line "<< __LINE__<<endl;
      exit(1);
    }
  if (Information > 0)
    {
      cout << "Only part of eigenvalues calculated in ComplexMatrix::LapackDiagonalize, line "<< __LINE__<<endl;
    }
  cout << "Schur matrix:"<<endl;
  for (int i=0; i<dimZ; ++i)
    {
      cout << TransQ[i*dimZ].r<<"+I*"<<TransQ[i*dimZ].i;
      for (int j=1; j<dimZ; ++j)
	cout << "  " <<TransQ[i*dimZ+j].r<<"+I*"<<TransQ[i*dimZ+j].i;      
      cout << endl;
    }

  cout << "Schur form of initial matrix:"<<endl;
  for (int i=0; i<dimZ; ++i)
    {
      cout << TmpMatrix[i*dimZ].r<<"+I*"<<TmpMatrix[i*dimZ].i;
      for (int j=1; j<dimZ; ++j)
	cout << "  " <<TmpMatrix[i*dimZ+j].r<<"+I*"<<TmpMatrix[i*dimZ+j].i;      
      cout << endl;
    }
  
  // store eigenvalues as return argument
  for (int i=0; i<DimensionM; ++i)
    M.SetMatrixElement(i,i,Eigenvalues[i].r, Eigenvalues[i].i);  
 
  // recover values of schur form and transition matrix
  for (int j=0;j<NbrRow;++j)
    for (int i=0; i<NbrRow;++i)
      {
	Q.SetMatrixElement(i, j, TransQ[i+j*NbrRow].r, TransQ[i+j*NbrRow].i);
	S.SetMatrixElement(i, j, TmpMatrix[i+j*NbrRow].r, TmpMatrix[i+j*NbrRow].i);
      }
  
  delete [] TmpMatrix;
  delete [] Tau;
  delete [] Eigenvalues;
  delete [] TransQ;
  //delete [] Scale;
  delete [] complexWork;
  
  return M;
}

#endif
