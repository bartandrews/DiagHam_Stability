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
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "MathTools/Complex.h"
#include "Matrix/RealDiagonalMatrix.h"

#include <math.h>


using std::endl;
using std::cout;

#ifdef HAVE_LAPACK

// TODO: review which of the lapack functions are used based on the selection of parallelized routines!

// binding to the LAPACK dgesvd function
//
extern "C" void FORTRAN_NAME(dgesvd)(const char* jobu, const char* jobv, const int* nbrRow, const int* nbrColumn, const double* matrix, const int* leadingDimension,
				     const double* singularValues, const double* uMatrix, const int* uLeadingDimension, const double* vMatrix, 
				     const int* vLeadingDimension, const double* workingArea, const int* workingAreaSize, const int* information);

// binding to the LAPACK dgesdd function
//
extern "C" void FORTRAN_NAME(dgesdd)(const char* jobz, const int* nbrRow, const int* nbrColumn, const double* matrix, const int* leadingDimension,
				     const double* singularValues, const double* uMatrix, const int* uLeadingDimension, const double* vMatrix, 
				     const int* vLeadingDimension, const double* workingArea, const int* workingAreaSize, const int* workingAreaInteger, const int* information);

// binding to the LAPACK function
//
extern "C" void FORTRAN_NAME(dgeev)(const char* jobVL, const char* jobVR, const int* nbrColumn, const double* matrix, const int* leadingDimension,
				    const double* eigenvaluesRealPart, const double* eigenvaluesImaginaryPart, 
				    const double* eigenvectorLeftMatrix, const int* leadingDimensionEigenvectorLeftMatrix,
				    const double* eigenvectorRightMatrix, const int* leadingDimensionEigenvectorRightMatrix,
				    const double* workingArea, const int* workingAreaSize, const int* information);

#endif

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

// constructor for one dimensional array
//
// array = one dimensional array where the matrix elements are stored (all components of the first column, then all components of the second column,...)
// nbrRow = number of rows
// nbrColumn = number of columns
// columnOrder = elements in array are ordered column-wise  (all components of the first column, then all components of the second column,...)

RealMatrix::RealMatrix(double* array, int nbrRow, int nbrColumn, bool columnOrder)
{
  if (columnOrder == true)
   {
     this->ColumnGarbageFlag = new int;
     *(this->ColumnGarbageFlag) = 1;
     this->NbrColumn = nbrColumn;
     this->NbrRow = nbrRow;
     this->TrueNbrRow = this->NbrRow;
     this->TrueNbrColumn = this->NbrColumn;
     this->Columns = new RealVector [this->NbrColumn];
     for (int i = 0; i < this->NbrColumn; i++)
      {
        this->Columns[i] = RealVector (this->NbrRow);
      }
  
     long Index = 0;
     for (int j = 0; j < this->NbrRow; j++)
       for (int i = 0; i < this->NbrColumn; i++)
         this->Columns[i][j] = array[Index++];

     this->MatrixType = Matrix::RealElements;
   }
  else //order by rows instead
   {
     this->ColumnGarbageFlag = new int;
     *(this->ColumnGarbageFlag) = 1;
     this->NbrColumn = nbrColumn;
     this->NbrRow = nbrRow;
     this->TrueNbrRow = this->NbrRow;
     this->TrueNbrColumn = this->NbrColumn;
     this->Columns = new RealVector [this->NbrColumn];
     for (int i = 0; i < this->NbrColumn; i++)
      {
        this->Columns[i] = RealVector (this->NbrRow);
      }
  
     long Index = 0;
     for (int i = 0; i < this->NbrRow; i++)
       for (int j = 0; j < this->NbrColumn; j++)
         this->Columns[i][j] = array[Index++];

     this->MatrixType = Matrix::RealElements;
   }
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

#ifdef __MPI__

// constructor from informations sent using MPI
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts or sends the vector
// broadcast = true if the vector is broadcasted

RealMatrix::RealMatrix(MPI::Intracomm& communicator, int id, bool broadcast)
{
  this->MatrixType = Matrix::RealElements;
  int TmpArray[4];
  if (broadcast == true)
    communicator.Bcast(TmpArray, 3, MPI::INT, id);      
  else
    communicator.Recv(TmpArray, 3, MPI::INT, id, 1);   
  this->NbrRow = TmpArray[0];
  this->NbrColumn = TmpArray[1];
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Columns = new RealVector [this->NbrColumn];
  if (TmpArray[2] == 1)
    {
      for (int i = 0; i < this->NbrColumn; i++)
	this->Columns[i] = RealVector (this->NbrRow, true);
    }
  else
    if (TmpArray[2] == 2)
      {
	for (int i = 0; i < this->NbrColumn; i++)
	  this->Columns[i] = RealVector (communicator, id, broadcast);
      }
  this->ColumnGarbageFlag = new int;
  *(this->ColumnGarbageFlag) = 1;
}

#endif

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
      for (int i = 0; i < this->NbrColumn; ++i)
	{
	  this->Columns[i] = RealVector (this->NbrRow);
	  for (int j = 0; j < this->NbrRow; ++j)
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

// copy a matrix into another (duplicating data)
//
// matrix = matrix to copy
// return value = reference on current matrix

RealMatrix& RealMatrix::Copy (RealMatrix& matrix)
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
      this->Columns = new RealVector[this->NbrColumn];
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

// Set all entries in matrix to zero
//

void RealMatrix::ClearMatrix ()
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i].ClearVector();
  return;
}

// set matrix to identity 
//

void RealMatrix::SetToIdentity()
{
  this->ClearMatrix();
  if (this->NbrColumn <= this->NbrRow)
    {
      for (int i = 0; i < this->NbrColumn; i++)
	this->Columns[i][i] = 1.0;
    }
  else
    {
      for (int i = 0; i < this->NbrRow; i++)
	this->Columns[i][i] = 1.0;
    }
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
    {
    return RealMatrix();
    }
  RealVector* TmpColumns = new RealVector [M2.NbrColumn];
  for (int i = 0; i < M2.NbrColumn; i++)
    {
      TmpColumns[i] = RealVector(M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	{
	  TmpColumns[i].Components[j] = 0.0;
	  for (int k = 0; k < M2.NbrRow; k++)
          {	
	    TmpColumns[i].Components[j] += M1.Columns[k].Components[j] * M2.Columns[i].Components[k];
          }
	}
    }
  return RealMatrix(TmpColumns, M2.NbrColumn);
}


// multiply two matrices
//
// M1 = first matrix
// M2 = matrix to multiply to M1
// return value = product of the two matrices

RealMatrix operator * (const RealMatrix & M1, const RealDiagonalMatrix & M2)
{
  if (M1.NbrColumn != M2.NbrRow)
    {
     return RealMatrix();
    }
  RealVector* TmpColumns = new RealVector [M2.NbrColumn];
  for (int i = 0; i < M2.NbrColumn; i++)
    {
      TmpColumns[i] = RealVector(M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	{
	  TmpColumns[i].Components[j] = M1.Columns[i].Components[j] * M2.DiagonalElements[i];	
	}
    }
  return RealMatrix(TmpColumns, M2.NbrColumn);
} 



// multiply two matrices
//
// M1 = first matrix
// M2 = matrix to multiply to M1
// return value = product of the two matrices

RealMatrix operator * (const  RealDiagonalMatrix & M1, const RealMatrix & M2) 
{
 if (M1.NbrColumn != M2.NbrRow)
   {
    return RealMatrix();
   }
  RealVector* TmpColumns = new RealVector [M2.NbrColumn];
  for (int i = 0; i < M2.NbrColumn; i++)
    {
      TmpColumns[i] = RealVector(M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	{
	  TmpColumns[i].Components[j] = M1.DiagonalElements[j] * M2.Columns[i].Components[j]; 
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


// multiply two matrices and add the result to the current one
//
// M1 = left matrix 
// M2 = right matrix 
// return value = reference on current matrix

RealMatrix& RealMatrix::AddMultiply (const RealMatrix  & M1, const RealMatrix & M2)
{
  if ((M1.NbrRow != this->NbrRow) || (M2.NbrColumn != this->NbrColumn))
    {
      cout << "incompatible matrix dimensions in RealMatrix::AddMultiply" << endl;
      return *this;
    }
  for (int i = 0; i < M2.NbrColumn; i++)
    {
      for (int j = 0; j < M1.NbrRow; j++)
	{
	  double Tmp = 0.0;
	  for (int k = 0; k < M2.NbrRow; k++)
	    {	
	      Tmp += M1.Columns[k].Components[j] * M2.Columns[i].Components[k];
	    }
	  this->Columns[i].Components[j] += Tmp;
	}
    }  
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

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

RealMatrix operator / (const RealMatrix& M1, const RealDiagonalMatrix& M2) 
{
  RealVector* TmpColumns = new RealVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; i++)
    {
      TmpColumns[i] = RealVector(M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	TmpColumns[i].Components[j] = M1.Columns[i].Components[j] / M2.DiagonalElements[i];
    }
  return RealMatrix(TmpColumns, M1.NbrRow);
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

// orthonormalize matrix column vectors
//
// return value = reference on current matrix

RealMatrix& RealMatrix::OrthoNormalizeColumns ()
{
  double* tmp = new double [this->NbrColumn];
  for (int i = 0; i < this->NbrColumn; i++)
    {
      for (int j = 0; j < i; j++)
	{
	  tmp[j] = this->Columns[j] * this->Columns[i];
	}
      for (int j = 0; j < i; j++)
	{
	  this->Columns[i].AddLinearCombination(-tmp[j], this->Columns[j]);
	}
      this->Columns[i].Normalize();
    }      
  delete[] tmp;
  return *this;
}

// orthonormalize matrix column vectors, computing the transformation matrix to the new orthonormal basis
//
// transformation= reference on the transformation matrix
// return value = reference on current matrix

RealMatrix& RealMatrix::OrthoNormalizeColumns (RealMatrix& transformation)
{
  double* tmp = new double [this->NbrColumn];
  transformation = RealMatrix (this->NbrColumn, this->NbrColumn, true);
  transformation.SetToIdentity();
  for (int i = 0; i < this->NbrColumn; i++)
    {
      for (int j = 0; j < i; j++)
	{
	  tmp[j] = this->Columns[j] * this->Columns[i];
	}
      for (int j = 0; j < i; j++)
	{
	  this->Columns[i].AddLinearCombination(-tmp[j], this->Columns[j]);
	  transformation[i].AddLinearCombination(-tmp[j], transformation[j]);
	}
      double TmpNorm = 1.0 / this->Columns[i].Norm();
      this->Columns[i] *= TmpNorm;
      transformation[i] *= TmpNorm;
    }
  delete[] tmp;
  return *this;
}

// transpose matrix
//
// return value = reference on current matrix

RealMatrix& RealMatrix::Transpose ()
{
  if (this->NbrRow == this->NbrColumn)
    {
      double tmp;
      for (int i = 0; i < this->NbrColumn; i++)
	for (int j = i + 1; j < this->NbrColumn; j++)
	  {
	    tmp = this->Columns[i].Components[j];
	    this->Columns[i].Components[j] = this->Columns[j].Components[i];
	    this->Columns[j].Components[i] = tmp;
	  }
    }
  else
    {
      RealVector* TmpColumns = new RealVector [this->NbrRow];
      for (int i = 0; i < this->NbrRow; i++)
	{
	  TmpColumns[i] = RealVector(this->NbrColumn);
	  for (int j = 0; j < this->NbrColumn; j++)
	    TmpColumns[i][j] = this->Columns[j][i];
	}
      if (this->ColumnGarbageFlag != 0)
	{
	  if ((*(this->ColumnGarbageFlag)) == 1)
	    {
	      delete[] this->Columns;
	    }
	  else
	    {
	      (*(this->ColumnGarbageFlag))--;
	    }
	}
      this->Columns = TmpColumns;
      int Tmp = this->NbrRow;
      this->NbrRow = this->NbrColumn;
      this->NbrColumn = Tmp;
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;      
    }
  return *this;
}

// duplicate and transpose a matrix
//
// return value = transposed matrix

RealMatrix RealMatrix::DuplicateAndTranspose ()
{
  RealMatrix TmpMatrix(this->NbrColumn, this->NbrRow);
  for (int i = 0; i < this->NbrRow; ++i)
    for (int j = 0; j < this->NbrColumn; ++j)
      {
	TmpMatrix.Columns[i][j] = this->Columns[j][i];
      }
  return TmpMatrix;
}

// evaluate matrix trace
//
// return value = matrix trace 

double RealMatrix::Tr ()
{
  double Tmp = 0.0;
  for (int i = 0; i < this->NbrRow; ++i)
    Tmp += this->Columns[i][i];
  return Tmp;
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

// compute singular value decomposition U D V^t
// 
// uMatrix = reference on the U matrix
// vMatrix = reference on the V matrix
// truncatedUVFlag = if false, set JOBZ = 'A' (returns full U, V matrices)
// return value = pointer on the diagonal elements of D

double* RealMatrix::SingularValueDecomposition(RealMatrix& uMatrix, RealMatrix& vMatrix, bool truncatedUVFlag)
{
#ifdef HAVE_LAPACK
  int MinDimension = this->NbrColumn;
  if (this->NbrColumn > this->NbrRow)
    MinDimension = this->NbrRow;
  double* SigmaMatrix = new double[MinDimension];
  int Information = 0;
  int WorkingAreaSize = -1;
  int IntegerWorkingAreaSize = -1;
  double TmpWorkingArea;
  int TmpIntegerWorkingArea;
  char Jobz;
  if (truncatedUVFlag == false)
     Jobz = 'A';
  else
     Jobz = 'S';
  double* TmpMatrix = new double [this->NbrRow * this->NbrColumn];
  int TotalIndex = 0;
  for (int j = 0; j < this->NbrColumn; ++j)
    {
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  TmpMatrix[TotalIndex] = this->Columns[j][i];
	  ++TotalIndex;
	}
    }
  double* TmpUMatrix = new double [this->NbrRow * this->NbrRow];
  for (int i = 0; i < (this->NbrRow * this->NbrRow); ++i)
    TmpUMatrix[i] = 0;
  double* TmpVMatrix = new double [this->NbrColumn * this->NbrColumn];
  for (int i = 0; i < (this->NbrColumn * this->NbrColumn); ++i)
    TmpVMatrix[i] = 0;

  int SizeLDU = this->NbrRow;
  int SizeLDVT;
  if (truncatedUVFlag == false)
    SizeLDVT = this->NbrColumn;
  else 
    SizeLDVT = MinDimension;


  FORTRAN_NAME(dgesdd)(&Jobz, &this->NbrRow, &this->NbrColumn, TmpMatrix, &this->NbrRow, SigmaMatrix, TmpUMatrix, &SizeLDU, TmpVMatrix, &SizeLDVT, &TmpWorkingArea, &WorkingAreaSize, &TmpIntegerWorkingArea, &Information); 
  WorkingAreaSize = (int) TmpWorkingArea;
  double* WorkingArea = new double [WorkingAreaSize];
  IntegerWorkingAreaSize = 8 * MinDimension;
  int* IntegerWorkingArea = new int [IntegerWorkingAreaSize];
  FORTRAN_NAME(dgesdd)(&Jobz, &this->NbrRow, &this->NbrColumn, TmpMatrix, &this->NbrRow, SigmaMatrix, TmpUMatrix, &SizeLDU, TmpVMatrix, &SizeLDVT, WorkingArea, &WorkingAreaSize, IntegerWorkingArea, &Information);
  uMatrix = RealMatrix(TmpUMatrix, this->NbrRow, this->NbrRow, false);
  vMatrix = RealMatrix(TmpVMatrix, this->NbrColumn, this->NbrColumn, false);
  delete[] TmpMatrix;
  delete[] TmpUMatrix;
  delete[] TmpVMatrix;
  delete[] WorkingArea;
  delete[] IntegerWorkingArea;
  return SigmaMatrix;
#else
  return 0;
#endif
}

// compute singular value decomposition U D V^t
// 
// uMatrix = reference on the U matrix
// diagonal = reference on the diagonal D matrix
// vMatrix = reference on the V matrix
// truncatedUVFlag = if false, set JOBZ = 'A' (returns full U, V matrices)

void RealMatrix::SingularValueDecomposition(RealMatrix& uMatrix, RealDiagonalMatrix& diagonal, RealMatrix& vMatrix, bool truncatedUVFlag)
{
  double* TmpDiag = this->SingularValueDecomposition(uMatrix, vMatrix, truncatedUVFlag);
  diagonal = RealDiagonalMatrix(TmpDiag, (uMatrix.NbrColumn>vMatrix.NbrRow)? vMatrix.NbrRow :uMatrix.NbrColumn );

}

// compute the diagonal part of the singular value decomposition U D V^t
// 
// return value = pointer on the diagonal elements of D

double* RealMatrix::SingularValueDecomposition()
{
#ifdef HAVE_LAPACK
  if ((this->NbrColumn == 1) || (this->NbrRow == 1))
    {
      double* SigmaMatrix = new double[1];
      SigmaMatrix[0] = 0.0;
      if (this->NbrColumn == 1)
	{
	  for (int i = 0; i < this->NbrRow; ++i)
	    SigmaMatrix[0] += this->Columns[0][i] * this->Columns[0][i];
	}
      else
	{
	  for (int i = 0; i < this->NbrColumn; ++i)
	    SigmaMatrix[0] += this->Columns[i][0] * this->Columns[i][0];
	}
      SigmaMatrix[0] = sqrt(SigmaMatrix[0]);
      return SigmaMatrix;
    }

  int MinDimension = this->NbrColumn;
  int MaxDimension = this->NbrRow;
  if (this->NbrColumn > this->NbrRow)
    {
      MinDimension = this->NbrRow;
      MaxDimension = this->NbrColumn;
    }

  // code for the slower (but less error proned) dgesvd
//   double* SigmaMatrix = new double[MinDimension];
//   int Information = 0;
//   int WorkingAreaSize = -1;
//   int IntegerWorkingAreaSize = 8 * MinDimension;
//   double TmpWorkingArea;
//   int TmpIntegerWorkingArea;
//   char Jobu = 'N';
//   char Jobvt = 'N';
//   double* TmpMatrix = new double [((long) this->NbrRow) * ((long) this->NbrColumn)];
//   long TotalIndex = 0l;
//   for (int j = 0; j < this->NbrColumn; ++j)
//     {
//       for (int i = 0; i < this->NbrRow; ++i)
// 	{
// 	  TmpMatrix[TotalIndex] = this->Columns[j][i];
// 	  ++TotalIndex;
// 	}
//     }
//   double* TmpUMatrix = new double [this->NbrColumn];
//   double* TmpVMatrix = new double [this->NbrRow];
//   int DummySize = 1;
//   FORTRAN_NAME(dgesvd)(&Jobu, &Jobvt, &this->NbrRow, &this->NbrColumn, TmpMatrix, &this->NbrRow, SigmaMatrix,
// 		       TmpUMatrix, &DummySize, TmpVMatrix, &DummySize, &TmpWorkingArea, &WorkingAreaSize, &Information);
//   WorkingAreaSize = (int) TmpWorkingArea;
//   double* WorkingArea = new double [WorkingAreaSize];
//   FORTRAN_NAME(dgesvd)(&Jobu, &Jobvt, &this->NbrRow, &this->NbrColumn, TmpMatrix, &this->NbrRow, SigmaMatrix,
// 		       TmpUMatrix, &DummySize, TmpVMatrix, &DummySize, WorkingArea, &WorkingAreaSize, &Information);

  double* SigmaMatrix = new double[MinDimension];
  int Information = 0;
  int WorkingAreaSize = -1;
  int IntegerWorkingAreaSize = 8 * MinDimension;
  double TmpWorkingArea;
  int TmpIntegerWorkingArea;
  char Jobz = 'N';
  double* TmpMatrix = new double [((long) this->NbrRow) * ((long) this->NbrColumn)];
  long TotalIndex = 0l;

  // it seems that highly non-square matrix (i.e this->NbrRow << this->NbrColumn) leads to
  // inconsistent SVD
  int LocalNbrColumn = this->NbrColumn;
  int LocalNbrRow = this->NbrRow;
  if (this->NbrRow >= this->NbrColumn)
    {
      for (int j = 0; j < this->NbrColumn; ++j)
	{
	  for (int i = 0; i < this->NbrRow; ++i)
	    {
	      TmpMatrix[TotalIndex] = this->Columns[j][i];
	      ++TotalIndex;
	    }
	}
    }
  else
    {
      LocalNbrColumn = this->NbrRow;
      LocalNbrRow = this->NbrColumn;
      for (int j = 0; j < LocalNbrColumn; ++j)
	{
	  for (int i = 0; i < LocalNbrRow; ++i)
	    {
	      TmpMatrix[TotalIndex] = this->Columns[i][j];
	      ++TotalIndex;
	    }
	}
    }

  double* TmpUMatrix = new double [LocalNbrColumn];
  double* TmpVMatrix = new double [LocalNbrRow];
  int DummySize = 1;
  FORTRAN_NAME(dgesdd)(&Jobz, &LocalNbrRow, &LocalNbrColumn, TmpMatrix, &LocalNbrRow, SigmaMatrix, TmpUMatrix, &DummySize, TmpVMatrix, &DummySize, &TmpWorkingArea, &WorkingAreaSize, &TmpIntegerWorkingArea, &Information);
  WorkingAreaSize = (int) TmpWorkingArea;
  int MinWorkingAreaSize = 3 * MinDimension;
  if (MaxDimension >= (7 * MinDimension))
    {
      MinWorkingAreaSize += MaxDimension;
    }
  else
    {
      MinWorkingAreaSize += 7 * MinDimension;
    }
  if (MinWorkingAreaSize > WorkingAreaSize)
    WorkingAreaSize = MinWorkingAreaSize;
  double* WorkingArea = new double [WorkingAreaSize];
  IntegerWorkingAreaSize = 8 * MinDimension;
  int* IntegerWorkingArea = new int [IntegerWorkingAreaSize];
  FORTRAN_NAME(dgesdd)(&Jobz, &LocalNbrRow, &LocalNbrColumn, TmpMatrix, &LocalNbrRow, SigmaMatrix, TmpUMatrix, &DummySize, TmpVMatrix, &DummySize, WorkingArea, &WorkingAreaSize, IntegerWorkingArea, &Information);
  delete[] IntegerWorkingArea;
  delete[] WorkingArea;
  delete[] TmpUMatrix;
  delete[] TmpVMatrix;
  delete[] TmpMatrix;
  return SigmaMatrix;
#else
  return 0;
#endif
}

// Diagonalize a real matrix using the LAPACK library
//
// M = reference on complex diagonal matrix where result has to be stored
// leftFlag = compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
// return value = reference on complex diagonal matrix

ComplexDiagonalMatrix& RealMatrix::LapackDiagonalize (ComplexDiagonalMatrix& M, bool leftFlag)
{
  if (this->NbrColumn != this->NbrRow)
    {
      cout << "only square matrices can be diagonalized" << endl;
      return M;
    }
#ifdef HAVE_LAPACK
  int Information = 0;
  int WorkingAreaSize = -1;
  char JobVL;
  char JobVR;
  if (leftFlag == true)
    {
      JobVL = 'N';
      JobVR = 'N';
    }
  else
    {
      JobVL = 'N';
      JobVR = 'N';
    }
  double* TmpMatrix = new double [((long) this->NbrColumn) * this->NbrRow];
  long TotalIndex = 0l;
  for (int j = 0; j < this->NbrColumn; ++j)
    {
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  TmpMatrix[TotalIndex] = this->Columns[j][i];
	  ++TotalIndex;
	}
    }
  double* TmpEigenvalueReal = new double[this->NbrColumn];
  double* TmpEigenvalueImaginary = new double[this->NbrColumn];
  int TmpLeadingDimension = 1;
  double* Dummy = 0;
  double TmpWorkingArea;
  FORTRAN_NAME(dgeev)(&JobVL, &JobVR, &this->NbrRow, TmpMatrix, &this->NbrRow, 
		      TmpEigenvalueReal, TmpEigenvalueImaginary,
		      Dummy, &TmpLeadingDimension, Dummy, &TmpLeadingDimension,
		      &TmpWorkingArea, &WorkingAreaSize, &Information);
  WorkingAreaSize = (int) TmpWorkingArea;
  double* WorkingArea = new double [WorkingAreaSize];
  FORTRAN_NAME(dgeev)(&JobVL, &JobVR, &this->NbrRow, TmpMatrix, &this->NbrRow, 
		      TmpEigenvalueReal, TmpEigenvalueImaginary,
		      Dummy, &TmpLeadingDimension, Dummy, &TmpLeadingDimension,
		      WorkingArea, &WorkingAreaSize, &Information);
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

ComplexDiagonalMatrix& RealMatrix::LapackDiagonalize (ComplexDiagonalMatrix& M, ComplexMatrix& Q, bool leftFlag)
{
  if (this->NbrColumn != this->NbrRow)
    {
      cout << "only square matrices can be diagonalized" << endl;
      return M;
    }
#ifdef HAVE_LAPACK
  int Information = 0;
  int WorkingAreaSize = -1;
  char JobVL;
  char JobVR;
  if (leftFlag == true)
    {
      JobVL = 'V';
      JobVR = 'N';
    }
  else
    {
      JobVL = 'N';
      JobVR = 'V';
    }
  double* TmpMatrix = new double [((long) this->NbrColumn) * this->NbrRow];
  double* TmpLeftEigenstates = new double [((long) this->NbrColumn) * this->NbrRow];
  long TotalIndex = 0l;
  for (int j = 0; j < this->NbrColumn; ++j)
    {
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  TmpMatrix[TotalIndex] = this->Columns[j][i];
	  ++TotalIndex;
	}
    }
  double* TmpEigenvalueReal = new double[this->NbrColumn];
  double* TmpEigenvalueImaginary = new double[this->NbrColumn];
  int TmpLeadingLeftDimension;
  int TmpLeadingRightDimension;
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
  double* Dummy = 0;
  double TmpWorkingArea;
  if (leftFlag == true)
    {
      FORTRAN_NAME(dgeev)(&JobVL, &JobVR, &this->NbrRow, TmpMatrix, &this->NbrRow, 
			  TmpEigenvalueReal, TmpEigenvalueImaginary,
			  TmpLeftEigenstates, &TmpLeadingLeftDimension, Dummy, &TmpLeadingRightDimension, 
			  &TmpWorkingArea, &WorkingAreaSize, &Information);
    }
  else
    {
      FORTRAN_NAME(dgeev)(&JobVL, &JobVR, &this->NbrRow, TmpMatrix, &this->NbrRow, 
			  TmpEigenvalueReal, TmpEigenvalueImaginary,
			  Dummy, &TmpLeadingLeftDimension, TmpLeftEigenstates, &TmpLeadingRightDimension, 
			  &TmpWorkingArea, &WorkingAreaSize, &Information);
    }
  WorkingAreaSize = (int) TmpWorkingArea;
  double* WorkingArea = new double [WorkingAreaSize];
  if (leftFlag == true)
    {
      FORTRAN_NAME(dgeev)(&JobVL, &JobVR, &this->NbrRow, TmpMatrix, &this->NbrRow, 
			  TmpEigenvalueReal, TmpEigenvalueImaginary,
			  TmpLeftEigenstates, &TmpLeadingLeftDimension, Dummy, &TmpLeadingRightDimension, 
			  WorkingArea, &WorkingAreaSize, &Information);
    }
  else
    {
      FORTRAN_NAME(dgeev)(&JobVL, &JobVR, &this->NbrRow, TmpMatrix, &this->NbrRow, 
			  TmpEigenvalueReal, TmpEigenvalueImaginary,
			  Dummy, &TmpLeadingLeftDimension, TmpLeftEigenstates, &TmpLeadingRightDimension, 
			  WorkingArea, &WorkingAreaSize, &Information);
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
  delete[] TmpLeftEigenstates;
  delete[] WorkingArea;
  delete[] TmpEigenvalueReal;
  delete[] TmpEigenvalueImaginary;
  delete[] TmpMatrix;
#endif
  return M;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const RealMatrix& P) 
{
  if ((P.NbrColumn > 0) && (P.NbrRow > 0))
    {
      for (int i = 0; i < (P.NbrRow - 1); ++i)
	{
	  for (int j = 0; j < (P.NbrColumn - 1); ++j)
	    Str << P.Columns[j].Components[i] << "    ";      
	  Str << P.Columns[P.NbrColumn - 1].Components[i] << endl;      
	}
      for (int j = 0; j < (P.NbrColumn - 1); j ++)
	Str << P.Columns[j].Components[P.NbrRow - 1] << "    ";      
      Str << P.Columns[P.NbrColumn - 1].Components[P.NbrRow - 1] << endl;      
    }
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

#ifdef __MPI__

// send a matrix to a given MPI process
// 
// communicator = reference on the communicator to use
// id = id of the destination MPI process
// return value = reference on the current matrix

Matrix& RealMatrix::SendMatrix(MPI::Intracomm& communicator, int id)
{
  communicator.Send(&this->MatrixType, 1, MPI::INT, id, 1);
  communicator.Send(&this->NbrRow, 1, MPI::INT, id, 1); 
  communicator.Send(&this->NbrColumn, 1, MPI::INT, id, 1); 
  int Acknowledge = 0;
  communicator.Recv(&Acknowledge, 1, MPI::INT, id, 1);
  if (Acknowledge != 0)
    return *this;
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i].SendVector(communicator, id);
  return *this;
}

// broadcast a matrix to all MPI processes associated to the same communicator
// 
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the matrix
// return value = reference on the current matrix

Matrix& RealMatrix::BroadcastMatrix(MPI::Intracomm& communicator,  int id)
{
  int TmpMatrixType = this->MatrixType;
  int TmpNbrRow = this->NbrRow;
  int TmpNbrColumn = this->NbrColumn;
  int Acknowledge = 0;
  communicator.Bcast(&TmpMatrixType, 1, MPI::INT, id);
  communicator.Bcast(&TmpNbrRow, 1, MPI::INT, id);
  communicator.Bcast(&TmpNbrColumn, 1, MPI::INT, id);
  if (this->MatrixType != TmpMatrixType)
    {
      Acknowledge = 1;
    }
  if (id != communicator.Get_rank())
    communicator.Send(&Acknowledge, 1, MPI::INT, id, 1);      
  else
    {
      int NbrMPINodes = communicator.Get_size();
      bool Flag = false;
      for (int i = 0; i < NbrMPINodes; ++i)
	if (id != i)
	  {
	    communicator.Recv(&Acknowledge, 1, MPI::INT, i, 1);      
	    if (Acknowledge == 1)
	      Flag = true;
	  }
      if (Flag == true)
	Acknowledge = 1;
    }
  communicator.Bcast(&Acknowledge, 1, MPI::INT, id);
  if (Acknowledge != 0)
    return *this;
  if ((TmpNbrRow != this->NbrRow) || (TmpNbrColumn != this->NbrColumn))
    {
      this->Resize(TmpNbrRow, TmpNbrColumn);      
    }
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i].BroadcastVector(communicator, id);
  return *this;
}

// receive a matrix from a MPI process
// 
// communicator = reference on the communicator to use 
// id = id of the source MPI process
// return value = reference on the current matrix

Matrix& RealMatrix::ReceiveMatrix(MPI::Intracomm& communicator, int id)
{
  int TmpMatrixType = 0;
  int TmpNbrRow = 0;
  int TmpNbrColumn = 0;
  communicator.Recv(&TmpMatrixType, 1, MPI::INT, id, 1);
  communicator.Recv(&TmpNbrRow, 1, MPI::INT, id, 1); 
  communicator.Recv(&TmpNbrColumn, 1, MPI::INT, id, 1); 
  if (TmpMatrixType != this->MatrixType)
    {
      TmpNbrRow = 1;
      communicator.Send(&TmpNbrRow, 1, MPI::INT, id, 1);
      return *this;
    }
  else
    {
      if ((TmpNbrRow != this->NbrRow) || (TmpNbrColumn != this->NbrColumn))
	{
	  this->Resize(TmpNbrRow, TmpNbrColumn);      
	}
      TmpNbrRow = 0;
      communicator.Send(&TmpNbrRow, 1, MPI::INT, id, 1);
    }
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i].ReceiveVector(communicator, id);
  return *this;
}

// add current matrix to the current matrix of a given MPI process
// 
// communicator = reference on the communicator to use 
// id = id of the destination MPI process
// return value = reference on the current matrix

Matrix& RealMatrix::SumMatrix(MPI::Intracomm& communicator, int id)
{
  int TmpMatrixType = this->MatrixType;
  int TmpNbrRow = this->NbrRow;
  int TmpNbrColumn = this->NbrColumn;
  int Acknowledge = 0;
  communicator.Bcast(&TmpMatrixType, 1, MPI::INT, id);
  communicator.Bcast(&TmpNbrRow, 1, MPI::INT, id);
  communicator.Bcast(&TmpNbrColumn, 1, MPI::INT, id);
  if ((this->MatrixType != TmpMatrixType) || (TmpNbrRow != this->NbrRow) || (TmpNbrColumn != this->NbrColumn))
    {
      Acknowledge = 1;
    }
  if (id != communicator.Get_rank())
    communicator.Send(&Acknowledge, 1, MPI::INT, id, 1);      
  else
    {
      int NbrMPINodes = communicator.Get_size();
      bool Flag = false;
      for (int i = 0; i < NbrMPINodes; ++i)
	if (id != i)
	  {
	    communicator.Recv(&Acknowledge, 1, MPI::INT, i, 1);      
	    if (Acknowledge == 1)
	      Flag = true;
	  }
      if (Flag == true)
	Acknowledge = 1;
    }
  communicator.Bcast(&Acknowledge, 1, MPI::INT, id);
  if (Acknowledge != 0)
    {
      return *this;
    }
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i].SumVector(communicator, id);
  return *this;
}

// reassemble matrix from a scattered one
// 
// communicator = reference on the communicator to use 
// id = id of the destination MPI process
// return value = reference on the current matrix

Matrix& RealMatrix::ReassembleMatrix(MPI::Intracomm& communicator, int id)
{
  if (id == communicator.Get_rank())
    {
      int NbrMPINodes = communicator.Get_size();
      int TmpArray[2];
      for (int i = 0; i < NbrMPINodes; ++i)
	if (id != i)
	  {
	    TmpArray[0] = 0;
	    TmpArray[1] = 0;
	    communicator.Recv(TmpArray, 2, MPI::INT, i, 1); 
	    int Lim = TmpArray[0] + TmpArray[1];
	    for (int i = TmpArray[0]; i < Lim; i++)
	      this->Columns[i].ReceiveVector(communicator, id);
	  }      
    }
  else
    {
      int TmpArray[2];
      TmpArray[0] = 0;
      TmpArray[1] = this->NbrColumn;
      communicator.Send(TmpArray, 2, MPI::INT, id, 1);
      int Lim = TmpArray[0] + TmpArray[1];
      for (int i = TmpArray[0]; i < Lim; i++)
	this->Columns[i].SendVector(communicator, id);
    }
  return *this;
}

// create a new matrix on each MPI node which is an exact clone of the broadcasted one
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the matrix
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new matrix 

Matrix* RealMatrix::BroadcastClone(MPI::Intracomm& communicator, int id)
{
  if (id == communicator.Get_rank())
    {
      communicator.Bcast(&this->MatrixType, 1, MPI::INT, id);
      int TmpArray[3];
      TmpArray[0] = this->NbrRow;
      TmpArray[1] = this->NbrColumn;
      TmpArray[2] = 2;
      communicator.Bcast(TmpArray, 3, MPI::INT, id);      
      for (int i = 0; i < this->NbrColumn; i++)
	this->Columns[i].BroadcastClone(communicator, id);
    }
  else
    {
      int Type = 0;
      communicator.Bcast(&Type, 1, MPI::INT, id);  
      return new RealMatrix(communicator, id);
    }
  return 0;
}

// create a new matrix on each MPI node with same size and same type but non-initialized components
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the matrix
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new matrix 

Matrix* RealMatrix::BroadcastEmptyClone(MPI::Intracomm& communicator, int id, bool zeroFlag)
{
  if (id == communicator.Get_rank())
    {
      communicator.Bcast(&this->MatrixType, 1, MPI::INT, id);
      int TmpArray[3];
      TmpArray[0] = this->NbrRow;
      TmpArray[1] = this->NbrColumn;
      TmpArray[2] = 0;
      if (zeroFlag == true)
	{
	  TmpArray[2] = 1;
	}
      communicator.Bcast(TmpArray, 3, MPI::INT, id);      
    }
  else
    {
      int Type = 0;
      communicator.Bcast(&Type, 1, MPI::INT, id);  
      return new RealMatrix(communicator, id);
    }
  return 0;
}

#endif

