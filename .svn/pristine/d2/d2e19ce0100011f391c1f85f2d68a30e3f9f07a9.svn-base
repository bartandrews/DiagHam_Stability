///////////////////////////////////////////////////////////////////////////////
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
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexLowerTriangularMatrix.h"
#include "Matrix/ComplexUpperTriangularMatrix.h"
#include <iostream>
#include <cstdlib>


using std::endl;
using std::cout;


#ifdef HAVE_LAPACK

typedef struct { float r, i; } LAcomplex;

// binding to the LAPACK zgesdd function
//
extern "C" void FORTRAN_NAME(zgesdd)(const char* jobz, const int* nbrRow, const int* nbrColumn, const doublecomplex* matrix, const int* leadingDimension,
				     const double* singularValues, const doublecomplex* uMatrix, const int* uLeadingDimension, const doublecomplex* vMatrix, 
				     const int* vLeadingDimension, const doublecomplex* workingArea, const int* workingAreaSize, const double* workingAreaR, const int* workingAreaInteger, const int* information);

// binding to the LAPACK zgetrf routine for LU decomposition and back-substitution
//
extern "C" void FORTRAN_NAME(zgetrf)(const int* dimensionM, const int* dimensionN, const doublecomplex* matrixA,
				     const int* leadingDimensionA, const int *ipiv, const int *info);

extern "C" void FORTRAN_NAME(zgetrs)(const char* transpose, const int* dimensionN, const int* numRHS,
				     const doublecomplex* matrixA, const int* leadingDimensionA, const int *ipiv,
				     const doublecomplex* matrixB, const int* leadingDimensionB, const int *info);

// binding to the LAPACK zgetri routine for matrix inversion using LU decomposition
//
extern "C" void FORTRAN_NAME(zgetri)(const int* dimensionM, const doublecomplex* matrixA, 
				     const int* leadingDimensionA, const int *ipiv, const doublecomplex* workingArea, 
				     const int* workingAreaSize, const int *info);

// binding to the LAPACK cgeqrf routine for matrix QR  decomposition
//
extern "C" void FORTRAN_NAME(zgeqrf) (const int* dimensionM, const int *  dimensionN,  const doublecomplex * matrixA,
				      const int* leadingDimensionA,  const doublecomplex * TAU,  const doublecomplex * WORK, const int * LWORK,
				      const int * INFO );
extern "C" void FORTRAN_NAME(zungqr) (const int* dimensionM, const int *  dimensionN, const int *  dimensionK ,  const doublecomplex * matrixA,
				      const int* leadingDimensionA,  const doublecomplex * TAU,  const doublecomplex * WORK, const int * LWORK,
				      const int * INFO );


// sequence of routines to be called for diagonalization via Hessenberg matrix QR decomposition
// balance (optional)
extern "C" void FORTRAN_NAME(zgebal)(const char* jobz, const int* dimensionN, const doublecomplex * matrixA,
				     const int* leadingDimensionA, const int *iLow, const int *iHigh,
				     const double *scale, const int *info);
// calculate hessenberg form
extern "C" void FORTRAN_NAME(zgehrd)(const int* dimensionN, const int *iLow, const int *iHigh,
				     const doublecomplex* matrixA, const int* leadingDimensionA,
				     const doublecomplex* tau, const doublecomplex* complexWork,
				     const int *lComplexWork, const int *info);
// calculate hessenberg form
extern "C" void FORTRAN_NAME(cgehrd)(const int* dimensionN, const int *iLow, const int *iHigh,
				     const LAcomplex* matrixA, const int* leadingDimensionA,
				     const LAcomplex* tau, const LAcomplex* complexWork,
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

extern "C" void FORTRAN_NAME(cgeevx)(const char* balanc,const char* jobVL,const char* jobVR, const char* sense,
				     const int* dimensionN, const LAcomplex *matrixA,
				     const int* leadingDimensionA, const LAcomplex* eigenValues,
				     const LAcomplex* matrixVL, const int* leadingDimensionVL,
				     const LAcomplex* matrixVR, const int* leadingDimensionVR,
				     const int *iLow, const int *iHigh, const double *scale,
				     const double *absNorm, const double *reciCondVal, const double *reciCondVec, 
				     const LAcomplex* complexWork, const int *lComplexWork,
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

#ifdef __LAPACK__
// constructor for one dimensional array
//
// array = one dimensional array where the matrix elements are stored (all components of the first column, then all components of the second column,...)
// nbrRow = number of rows
// nbrColumn = number of columns
// columnOrder = elements in array are ordered column-wise  (all components of the first column, then all components of the second column,...)

ComplexMatrix::ComplexMatrix(doublecomplex* array, int nbrRow, int nbrColumn, bool columnOrder)
{
  if (columnOrder == true)
   {
     this->Flag.Initialize();
     this->NbrColumn = nbrColumn;
     this->NbrRow = nbrRow;
     this->TrueNbrRow = this->NbrRow;
     this->TrueNbrColumn = this->NbrColumn;
     this->Columns = new ComplexVector [this->NbrColumn];
     for (int i = 0; i < this->NbrColumn; i++)
      {
        this->Columns[i] = ComplexVector (this->NbrRow);
      }
  
     long Index = 0;
     for (int j = 0; j < this->NbrRow; j++)
       for (int i = 0; i < this->NbrColumn; i++)
        {
          this->Columns[i][j].Re = array[Index].r;
          this->Columns[i][j].Im = array[Index].i;
          Index++; 
        }

     this->MatrixType = Matrix::ComplexElements;
   }
  else //order by rows instead
   {
     this->Flag.Initialize();
     this->NbrColumn = nbrColumn;
     this->NbrRow = nbrRow;
     this->TrueNbrRow = this->NbrRow;
     this->TrueNbrColumn = this->NbrColumn;
     this->Columns = new ComplexVector [this->NbrColumn];
     for (int i = 0; i < this->NbrColumn; i++)
      {
        this->Columns[i] = ComplexVector (this->NbrRow);
      }
  
     long Index = 0;
     for (int i = 0; i < this->NbrRow; i++)
       for (int j = 0; j < this->NbrColumn; j++)
        {
         this->Columns[i][j].Re = array[Index].r;
         this->Columns[i][j].Im = array[Index].i;
         Index++;
        }

     this->MatrixType = Matrix::ComplexElements;
   }
}
#endif

#ifdef __MPI__

// constructor from informations sent using MPI
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts or sends the vector
// broadcast = true if the vector is broadcasted

ComplexMatrix::ComplexMatrix(MPI::Intracomm& communicator, int id, bool broadcast)
{
  this->MatrixType = Matrix::ComplexElements;
  int TmpArray[4];
  if (broadcast == true)
    communicator.Bcast(TmpArray, 3, MPI::INT, id);      
  else
    communicator.Recv(TmpArray, 3, MPI::INT, id, 1);   
  this->NbrRow = TmpArray[0];
  this->NbrColumn = TmpArray[1];
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Columns = new ComplexVector [this->NbrColumn];
  if (TmpArray[2] == 1)
    {
      for (int i = 0; i < this->NbrColumn; i++)
	this->Columns[i] = ComplexVector (this->NbrRow, true);
    }
  else
    {
      if (TmpArray[2] == 2)
	{
	  for (int i = 0; i < this->NbrColumn; i++)
	    this->Columns[i] = ComplexVector (communicator, id, broadcast);
	}
    }
  this->Flag.Initialize();
}

#endif

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
	      this->Columns[i][j] = Tmp;
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

// copy a matrix into another (duplicating data)
//
// matrix = matrix to copy
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::Copy (ComplexMatrix& matrix)
{
  this->Resize(matrix.NbrRow, matrix.NbrColumn);
  for (int j = 0; j < this->NbrColumn; ++j)
    this->Columns[j].Copy(matrix.Columns[j]);
  return *this;
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

// Set all entries in matrix to zero
//

void ComplexMatrix::ClearMatrix ()
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i].ClearVector();
  return;
}

// set matrix to identity 
//

void ComplexMatrix::SetToIdentity()
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

// conjugate an complex matrix with an unitary matrix (Ut M U)
//
// UnitaryM = unitary matrix to use
// return value = conjugated matrix

ComplexMatrix ComplexMatrix::Conjugate(ComplexMatrix& UnitaryM)
{
  ComplexMatrix TmpMatrix (UnitaryM.NbrColumn, UnitaryM.NbrColumn, true);
  for (int i = 0; i < UnitaryM.NbrColumn; ++i)
    for (int j = 0; j < UnitaryM.NbrColumn; ++j)
      {
	Complex Tmp = 0.0;
	for (int k = 0; k < UnitaryM.NbrRow; ++k)
	  for (int l = 0; l < UnitaryM.NbrRow; ++l)
	    {
	      Tmp += Conj(UnitaryM.Columns[i][k]) * this->Columns[l][k] * UnitaryM.Columns[j][l];
	    }
	TmpMatrix.Columns[j][i] = Tmp;
      }
  return TmpMatrix;
}

// conjugate an complex matrix with a complex transposed unitary matrix (U M Ut)
//
// UnitaryM = unitary matrix to use
// return value = conjugated matrix  

ComplexMatrix ComplexMatrix::InvConjugate(ComplexMatrix& UnitaryM)
{
  ComplexMatrix TmpMatrix (UnitaryM.NbrRow, UnitaryM.NbrRow, true);
  for (int i = 0; i < UnitaryM.NbrRow; ++i)
    for (int j = 0; j < UnitaryM.NbrRow; ++j)
      {
	Complex Tmp = 0.0;
	for (int k = 0; k < UnitaryM.NbrColumn; ++k)
	  for (int l = 0; l < UnitaryM.NbrColumn; ++l)
	    {
	      Tmp += UnitaryM.Columns[k][i] * this->Columns[l][k] * Conj(UnitaryM.Columns[l][j]);
	    }
	TmpMatrix.Columns[j][i] = Tmp;
      }
  return TmpMatrix;
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


// multiply two matrices
//
// M1 = first matrix
// M2 = matrix to multiply to M1
// return value = product of the two matrices

ComplexMatrix operator * (const RealDiagonalMatrix & M1, const ComplexMatrix& M2)
{
  if (M1.NbrColumn != M2.NbrRow)
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M2.NbrColumn];
  for (int i = 0; i < M2.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; ++j)
	{
      TmpColumns[i].Components[j].Re = M1.DiagonalElements[j] *  M2.Columns[i].Components[j].Re;
      TmpColumns[i].Components[j].Im = M1.DiagonalElements[j] *  M2.Columns[i].Components[j].Im;
	   }
     }
  return ComplexMatrix(TmpColumns, M2.NbrColumn);
}

// multiply two matrices
//
// M1 = first matrix
// M2 = matrix to multiply to M1
// return value = product of the two matrices

ComplexMatrix operator * (const ComplexMatrix& M1, const RealDiagonalMatrix& M2)
{
  if (M1.NbrColumn != M2.NbrRow)
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M2.NbrColumn];
  for (int i = 0; i < M2.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; ++j)
	{
      TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re * M2.DiagonalElements[i];
      TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im * M2.DiagonalElements[i];
	 }
	}
  return ComplexMatrix(TmpColumns, M2.NbrColumn);
}


// multiply two matrices
//
// M1 = first matrix
// M2 = matrix to multiply to M1
// return value = product of the two matrices

ComplexMatrix operator * (const ComplexDiagonalMatrix & M1, const ComplexMatrix& M2)
{
  if (M1.NbrColumn != M2.NbrRow)
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M2.NbrColumn];
  for (int i = 0; i < M2.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; ++j)
	{
      TmpColumns[i].Components[j].Re = M1.DiagonalElements[j].Re *  M2.Columns[i].Components[j].Re - M1.DiagonalElements[j].Im *  M2.Columns[i].Components[j].Im;
      TmpColumns[i].Components[j].Im = M1.DiagonalElements[j].Re *  M2.Columns[i].Components[j].Im + M1.DiagonalElements[j].Im *  M2.Columns[i].Components[j].Re;
	   }
     }
  return ComplexMatrix(TmpColumns, M2.NbrColumn);
}

// multiply two matrices
//
// M1 = first matrix
// M2 = matrix to multiply to M1
// return value = product of the two matrices

ComplexMatrix operator * (const ComplexMatrix& M1, const ComplexDiagonalMatrix& M2)
{
  if (M1.NbrColumn != M2.NbrRow)
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M2.NbrColumn];
  for (int i = 0; i < M2.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; ++j)
	{
      TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re * M2.DiagonalElements[i].Re - M1.Columns[i].Components[j].Im * M2.DiagonalElements[i].Re;
      TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im * M2.DiagonalElements[i].Re +  M1.Columns[i].Components[j].Re * M2.DiagonalElements[i].Im;
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

// multiply a matrix to the right by another matrix without using temporary matrix
//
// M = matrix used as multiplicator
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::Multiply (const ComplexMatrix& M)
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

// multiply a matrix to the right by another matrix without using temporary matrix
//
// M = matrix used as multiplicator
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::Multiply (const RealMatrix& M)
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

ComplexMatrix& ComplexMatrix::Multiply (const ComplexMatrix& M, int startLine, int nbrLine)
{
  if ((M.NbrRow != this->NbrColumn) || (M.NbrColumn >  this->TrueNbrColumn))
    return *this;
  int EndLine  = nbrLine + startLine;
  Complex* TmpElements = new Complex [this->NbrColumn];
  Complex Tmp;
  for (int i = startLine; i < EndLine; ++i)
    {
      for (int k = 0; k < this->NbrColumn; ++k)
	TmpElements[k] = this->Columns[k].Components[i];
      for (int j = 0; j < M.NbrColumn; ++j)
	{
// 	  Tmp = TmpElements[0] * M.Columns[j].Components[0];
// 	  for (int k = 1; k < this->NbrColumn; ++k)
// 	    {
// 	      Tmp += TmpElements[k] * M.Columns[j].Components[k];
// 	    }    
	  Tmp = 0.0;
	  Complex* TmpColumn = M.Columns[j].Components;
	  for (int k = 0; k < this->NbrColumn; ++k)
	    {
	      Tmp += TmpElements[k] * TmpColumn[k];
	    }    
	  this->Columns[j].Components[i] = Tmp;
	}  
    }
  delete[] TmpElements;
  return *this;
}

// multiply a matrix to the right by another matrix without using temporary matrix and in a given range of indices
// beware the matrix is not resized after multiplication in order the operation to be thread safe
//
// M = matrix used as multiplicator
// startLine = starting line in destination matrix
// nbrLine = number of lines to multiply
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::Multiply (const RealMatrix& M, int startLine, int nbrLine)
{
  if ((M.NbrRow != this->NbrColumn) || (M.NbrColumn >  this->TrueNbrColumn))
    return *this;
  int EndLine  = nbrLine + startLine;
  Complex* TmpElements = new Complex [this->NbrColumn];
  Complex Tmp;
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

// multiply two matrices, taking the hermitian conjugate of the left nmatrix first (i.e. M1^+ M2)
//
// M1 = first matrix
// M2 = matrix to multiply to M1
// return value = product of the two matrices

ComplexMatrix HermitianMultiply (const ComplexMatrix& M1, const ComplexMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M2.NbrColumn];
  for (int i = 0; i < M2.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrColumn);
      for (int j = 0; j < M1.NbrColumn; ++j)
	{
	  TmpColumns[i][j] = M1.Columns[j] * M2.Columns[i];
	}
    }
  return ComplexMatrix(TmpColumns, M2.NbrColumn);
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


// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

ComplexMatrix operator / (const ComplexMatrix& M1, const RealDiagonalMatrix& M2) 
{
  ComplexVector* TmpColumns = new ComplexVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; i++)
    {
      double Tmp = 1.0 /  M2.DiagonalElements[i];
      TmpColumns[i] = ComplexVector(M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; j++)
	{
        TmpColumns[i].Components[j].Re = M1.Columns[i].Components[j].Re * Tmp;
        TmpColumns[i].Components[j].Im = M1.Columns[i].Components[j].Im * Tmp;
	}
    }
  return ComplexMatrix(TmpColumns, M2.NbrRow);
 }

// add another complex matrices
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

// add another hermitian matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::operator += (const HermitianMatrix& M) 
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow))
    return *this;
  Complex TmpC;
  for (int i = 0; i < this->NbrColumn; ++i)
    {
      this->Columns[i][i] += M.DiagonalElements[i];
      for (int j = i+1; j < this->NbrColumn; ++j)
	{
	  M.GetMatrixElement(i, j, TmpC);
	  this->Columns[i][j] += TmpC;
	  this->Columns[j][i] += Conj(TmpC);
	}
    }
  return *this;
}

// add a linear combination of another complex matrix
// x = prefactor for added terms
// M = added matrix
ComplexMatrix& ComplexMatrix::AddLinearCombination(double x, const ComplexMatrix &M)
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow))
    return *this;
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i].AddLinearCombination(x,M.Columns[i]);
  return *this;
}

// add a linear combination of another complex matrix
// x = prefactor for added terms
// M = added matrix
ComplexMatrix& ComplexMatrix::AddLinearCombination(double x, const HermitianMatrix &M)
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow))
    return *this;
  Complex TmpC;
  for (int i = 0; i < this->NbrColumn; ++i)
    {
      this->Columns[i][i] += x*M.DiagonalElements[i];
      for (int j = i+1; j < this->NbrColumn; ++j)
	{
	  M.GetMatrixElement(i, j, TmpC);
	  this->Columns[j][i] += x*TmpC;
	  this->Columns[i][j] += x*Conj(TmpC);
	}
    }
  return *this;
}

// add a linear combination of another complex matrix
// x = prefactor for added terms
// M = added matrix
ComplexMatrix& ComplexMatrix::AddLinearCombination(const Complex &x, const ComplexMatrix &M)
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow))
    return *this;
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i].AddLinearCombination(x,M.Columns[i]);
  return *this;
}

// add a linear combination of another complex matrix
// x = prefactor for added terms
// M = added matrix
ComplexMatrix& ComplexMatrix::AddLinearCombination(const Complex &x, const HermitianMatrix &M)
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow))
    return *this;
  Complex TmpC;
  for (int i = 0; i < this->NbrColumn; ++i)
    {
      this->Columns[i][i] += x*M.DiagonalElements[i];
      for (int j = i+1; j < this->NbrColumn; ++j)
	{
	  M.GetMatrixElement(i, j, TmpC);
	  this->Columns[j][i] += x*TmpC;
	  this->Columns[i][j] += x*Conj(TmpC);
	}
    }
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

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::operator *= (const Complex& x) 
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] *= x;
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::operator /= (const Complex& x)
{
  Complex InvX = 1.0 / x;
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] *= InvX;
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

ComplexMatrix& ComplexMatrix::OrthoNormalizeColumns (ComplexMatrix& transformation)
{
  Complex* tmp = new Complex [this->NbrColumn];
  transformation = ComplexMatrix (this->NbrColumn, this->NbrColumn, true);
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

// compute the hermitian transpose of the current matrix
//
// return value = reference on the current matrix

ComplexMatrix& ComplexMatrix::HermitianTranspose ()
{
  if (this->NbrRow == this->NbrColumn)
    {
      Complex tmp;
      for (int i = 0; i < this->NbrColumn; i++)
	{
	  for (int j = i + 1; j < this->NbrColumn; j++)
	    {
	      tmp = this->Columns[i].Components[j];
	      this->Columns[i].Components[j] = Conj(this->Columns[j].Components[i]);
	      this->Columns[j].Components[i] = Conj(tmp);
	    }
	  this->Columns[i].Components[i] = Conj(this->Columns[i].Components[i]);
	}
    }
  else
    {
      ComplexVector* TmpColumns = new ComplexVector [this->NbrRow];
      for (int i = 0; i < this->NbrRow; i++)
	{
	  TmpColumns[i] = ComplexVector(this->NbrColumn);
	  for (int j = 0; j < this->NbrColumn; j++)
	    TmpColumns[i][j] = Conj(this->Columns[j][i]);
	}
      if (this->Flag.Shared() == false)
	{
	  delete[] this->Columns;
	}
      this->Flag.Initialize();
      this->Columns = TmpColumns;
      int Tmp = this->NbrRow;
      this->NbrRow = this->NbrColumn;
      this->NbrColumn = Tmp;
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;
    }
  return *this;
}

// compute the transpose of the current matrix
//
// return value = reference on the current matrix

ComplexMatrix& ComplexMatrix::Transpose ()
{
  if (this->NbrRow == this->NbrColumn)
    {
      Complex tmp;
      for (int i = 0; i < this->NbrColumn; i++)
	{
	  for (int j = i + 1; j < this->NbrColumn; j++)
	    {
	      tmp = this->Columns[i].Components[j];
	      this->Columns[i].Components[j] = this->Columns[j].Components[i];
	      this->Columns[j].Components[i] = tmp;
	    }
	  this->Columns[i].Components[i] = this->Columns[i].Components[i];
	}
    }
  else
    {
      ComplexVector* TmpColumns = new ComplexVector [this->NbrRow];
      for (int i = 0; i < this->NbrRow; i++)
	{
	  TmpColumns[i] = ComplexVector(this->NbrColumn);
	  for (int j = 0; j < this->NbrColumn; j++)
	    TmpColumns[i][j] = this->Columns[j][i];
	}
      if (this->Flag.Shared() == false)
	{
	  delete[] this->Columns;
	}
      this->Flag.Initialize();
      this->Columns = TmpColumns;
      int Tmp = this->NbrRow;
      this->NbrRow = this->NbrColumn;
      this->NbrColumn = Tmp;
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;
    }
  return *this;
}

// compute the complex conjugate of the current matrix
//
// return value = reference on the current matrix

ComplexMatrix& ComplexMatrix::ComplexConjugate ()
{
  for (int i = 0; i < this->NbrColumn; i++)
    {
      for (int j = 0; j < this->NbrRow; j++)
	{
	  this->Columns[i].Components[j].Im *= -1.0;
	}
    }
  return *this;
}

// compute the number of non-zero matrix elements (zero having strictly zero square norm)
//
// return value = number of non-zero matrix elements

long ComplexMatrix::ComputeNbrNonZeroMatrixElements()
{
  long NbrNonZero = 0l;
  for (int j=0; j < this->NbrColumn; ++j)
    for (int i=0; i < this->NbrRow; ++i)
      if (SqrNorm(this->Columns[j][i]) != 0.0)
	++NbrNonZero;
  return NbrNonZero;
}

// check if a complex matrix is hermitian
//
// error = maximum relative error allowed
// return value = true if the matrix is hermitian

bool ComplexMatrix::TestHermitian(double error)
{
  cout << "check hermiticity" << endl;
  if (this->NbrRow != this->NbrColumn)
    {
      cout << "error, not a square matrix" << endl;
      return false;
    }
  Complex Tmp1;
  Complex Tmp2;
  double AverageNorm = 0.0;
  for (int i = 0; i < this->NbrRow; ++i)
    for (int j = i; j < this->NbrColumn; ++j)
      {
	AverageNorm += Norm(this->Columns[j][i]);
      }
  AverageNorm /= 0.5 * ((double) this->NbrRow) * ((double) (this->NbrRow + 1));
  for (int i = 0; i < this->NbrRow; ++i)
    for (int j = i; j < this->NbrColumn; ++j)
      {
	Tmp1 = this->Columns[j][i];
	Tmp2 = this->Columns[i][j];
	if (Norm(Tmp1 - Conj(Tmp2)) > (error * AverageNorm))
	  {
	    cout << "error at " << i << " " << j << " : " << Tmp1 << " " << Tmp2 << " " << Norm(Tmp1 - Conj(Tmp2)) << " (should be lower than " << (error * AverageNorm) << ")" << endl;
	  }
      }
  cout << "check done" << endl;
  return true;
}

// discard the columns that are strictly zero
//

void ComplexMatrix::RemoveZeroColumns()
{
  if ((this->NbrRow > 0) && (this->NbrColumn > 0))
    {
      bool* TmpTestZero = new bool[this->NbrColumn];
      int TmpNbrColumn = 0;
      for (int i = 0; i < this->NbrColumn; ++i)
	{
	  if (this->Columns[i].SqrNorm() == 0.0)
	    {
	      TmpTestZero[i] = true;
	    }
	  else
	    {
	      TmpTestZero[i] = false;
	      ++TmpNbrColumn;	      
	    }
	}
      if (TmpNbrColumn == 0)
	{
	  (*this) = ComplexMatrix();
	}
      else
	{
	  ComplexVector* TmpVectors = new ComplexVector[TmpNbrColumn];
	  TmpNbrColumn = 0;
	  for (int i = 0; i < this->NbrColumn; ++i)
	    {
	      if (TmpTestZero[i] == false)
		{
		  TmpVectors[TmpNbrColumn] = this->Columns[i];
		  ++TmpNbrColumn;
		}	      
	    }
	  (*this) = ComplexMatrix(TmpVectors, TmpNbrColumn);
	}
      delete[] TmpTestZero;
    }
}


// discard the rows that are strictly zero
//

void ComplexMatrix::RemoveZeroRows()
{
  if ((this->NbrRow > 0) && (this->NbrColumn > 0))
    {
      bool* TmpTestZero = new bool[this->NbrRow];
      int TmpNbrRow = 0;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  TmpTestZero[i] = true;
	  for (int j = 0; (j < this->NbrColumn) && (TmpTestZero[i] == true); ++j)
	    {
	      if ((this->Columns[j][i].Re != 0.0) || (this->Columns[j][i].Im != 0.0))
		{
		  TmpTestZero[i] = false;
		  ++TmpNbrRow;
		}	      
	    }
	}
      if (TmpNbrRow == 0)
	{
	  (*this) = ComplexMatrix();
	}
      else
	{
	  ComplexVector* TmpVectors = new ComplexVector[this->NbrColumn];
	  for (int j = 0; j < this->NbrColumn; ++j)
	    TmpVectors[j] = ComplexVector(TmpNbrRow);
	  TmpNbrRow = 0;
	  for (int i = 0; i < this->NbrRow; ++i)
	    {
	      if (TmpTestZero[i] == false)
		{
		  for (int j = 0; j < this->NbrColumn; ++j)
		    TmpVectors[j][TmpNbrRow] = this->Columns[j][i];
		  ++TmpNbrRow;
		}
	    }
	  (*this) = ComplexMatrix(TmpVectors, this->NbrColumn);
	}
      delete[] TmpTestZero;
    }
}


// apply a sequence of row permutations
//
// permutations = array that list all the permutations. Each permutation is given at a pair corresponding to an index i and the i-th entry in the array (i.e. i <-> permutations[i]). 
//                The sequence is performed from the latest entry of permutations to the first one
// nbrPermutations = number of permutations to apply

void ComplexMatrix::ApplyRowPermutations(int* permutations, int nbrPermutations)
{
  Complex Tmp;
  for (int i = nbrPermutations - 1; i >= 0; --i)
    {
      if (i != permutations[i])
	{
	  int TmpIndex = permutations[i];
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      Tmp = this->Columns[j][i];
	      this->Columns[j][i] = this->Columns[j][TmpIndex];
	      this->Columns[j][TmpIndex] = Tmp;
	    }
	}
    }
}

// evaluate the real part of the matrix trace
//
// return value = real part of the matrix trace 

double ComplexMatrix::Tr ()
{
  double Trace = 0.0;
  int Max = this->NbrRow;
  if (this->NbrColumn < Max)
    Max = this->NbrColumn;
  for (int i = 0; i < Max; ++i)
    Trace += this->Columns[i][i].Re;
  return Trace;
}

// evaluate the matrix trace
//
// return value = matrix trace 

Complex ComplexMatrix::ComplexTr ()
{
  Complex Trace = 0.0;
  int Max = this->NbrRow;
  if (this->NbrColumn < Max)
    Max = this->NbrColumn;
  for (int i = 0; i < Max; ++i)
    Trace += this->Columns[i][i].Re;
  return Trace;
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

// compute the LU decompostion of the matrix 
// 
// lowerMatrix = reference on the matrix where the lower triangular matrix will be stored
// upperMatrix = reference on the matrix where the upper triangular matrix will be stored
// return value = array that  describe the additional row permutation

int* ComplexMatrix::LUDecomposition(ComplexLowerTriangularMatrix& lowerMatrix, ComplexUpperTriangularMatrix& upperMatrix)
{
#ifdef __LAPACKONLY__
  return this->LapackLUDecomposition(lowerMatrix, upperMatrix);
#endif
  if (this->NbrRow != this->NbrColumn)
    {
      cout << "LU decomposition is only performed for square matrices" << endl;
      return 0; 
    }
  cout << "warning, ComplexMatrix::LUDecomposition requires Lapack" << endl;
  return 0;
}
  
// compute the invert of a matrix from its PLU decomposition
// 
// lowerMatrix = reference on the matrix where the lower triangular matrix
// upperMatrix = reference on the matrix where the upper triangular matrix
// permutations = array that list all the permutations defining P. Each permutation is given at a pair corresponding to an index i and 
//                the i-th entry in the array (i.e. i <-> permutations[i]). The sequence is performed from the latest entry of permutations to the first one
// return value = inverted matrix

ComplexMatrix InvertMatrixFromLUDecomposition(ComplexLowerTriangularMatrix& lowerMatrix, ComplexUpperTriangularMatrix& upperMatrix, int* permutations)
{
  ComplexMatrix TmpM (lowerMatrix.GetNbrRow(), upperMatrix.GetNbrColumn());
  ComplexVector TmpV (lowerMatrix.GetNbrRow(), true);
  Complex Tmp;
  int TmpIndex = 0;
  for (int i = 0; i < TmpM.NbrColumn; ++i)
    {
      TmpV.ClearVector();
      TmpIndex = i;
      for (int j = 0; j < TmpM.NbrRow; ++j)
	{
	  if (j != permutations[j])
	    {
	      if (j == TmpIndex)
		{
		  TmpIndex = permutations[j];
		}
	      else
		{
		  if (permutations[j] == TmpIndex)
		    {
		      TmpIndex = j;
		    }
		}
	    }
	}
      TmpV[TmpIndex] = 1.0;
      if (lowerMatrix.SolveLinearEquation(TmpM.Columns[i], TmpV) == false)
	{
	  return ComplexMatrix();
	}
      ComplexVector TmpV2 = TmpV;
      TmpV = TmpM.Columns[i];
      TmpM.Columns[i] = TmpV2;
      if (upperMatrix.SolveLinearEquation(TmpM.Columns[i], TmpV) == false)
	{
	  return ComplexMatrix();
	}      
    }
  return TmpM;
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




// calculate a determinant using the LAPACK library (conserving current matrix)
//
Complex ComplexMatrix::LapackDeterminant ()
{
  if (this->NbrColumn != this->NbrRow)
    return 0.0;
#ifdef __LAPACK__
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

#else
  cout << "Warning, using LapackDeterminant without the lapack library" << endl;
  Complex Result(0.0,0.0);
#endif
  return Result;
}


// compute singular value decomposition U D V^t
// 
// uMatrix = reference on the U matrix
// vMatrix = reference on the V matrix
// truncatedUVFlag = if false, set JOBZ = 'A' (returns full U, V matrices)
// return value = pointer on the diagonal elements of D

double* ComplexMatrix::SingularValueDecomposition(ComplexMatrix& uMatrix, ComplexMatrix& vMatrix, bool truncatedUVFlag)
{
#ifdef HAVE_LAPACK
  int MinDimension = this->NbrColumn;
  int MaxDimension = this->NbrRow;
  if (this->NbrColumn > this->NbrRow)
    {
      MinDimension = this->NbrRow;
      MaxDimension = this->NbrColumn;
    }
  double* SigmaMatrix = new double[MinDimension];
  int Information = 0;
  int WorkingAreaSize = -1;
  int IntegerWorkingAreaSize = -1;
  doublecomplex TmpWorkingArea;
  int TmpIntegerWorkingArea;
  char Jobz;
  if (truncatedUVFlag == false)
     Jobz = 'A';
  else
     Jobz = 'S';

  doublecomplex* TmpMatrix = new doublecomplex [this->NbrRow * this->NbrColumn];

  Complex *TmpColumn;
  for (int j = 0; j < this->NbrColumn; ++j)
    {
      TmpColumn=this->Columns[j].Components;
      for (int i = 0; i < this->NbrRow;++i)
	{
	  TmpMatrix[i+j*this->NbrRow].r=TmpColumn[i].Re;
	  TmpMatrix[i+j*this->NbrRow].i=TmpColumn[i].Im;
	}
    }

  doublecomplex* TmpUMatrix = new doublecomplex [this->NbrRow * this->NbrRow];
  for (int i = 0; i < (this->NbrRow * this->NbrRow); ++i)
   {
    TmpUMatrix[i].r = 0.0;
    TmpUMatrix[i].i = 0.0;
   }
  doublecomplex* TmpVMatrix = new doublecomplex [this->NbrColumn * this->NbrColumn];
  for (int i = 0; i < (this->NbrColumn * this->NbrColumn); ++i)
   {
    TmpVMatrix[i].r = 0.0;
    TmpVMatrix[i].i = 0.0;
   }

  int SizeLDU = this->NbrRow;
  int SizeLDVT;
  if (truncatedUVFlag == false)
    SizeLDVT = this->NbrColumn;
  else 
    SizeLDVT = MinDimension;

  int RWorkDim, LRWorkDim;

  LRWorkDim = MinDimension;
  if ((5*MinDimension+7) > (2*MaxDimension+2*MinDimension+1))
     LRWorkDim *= (5*MinDimension+7);
  else
     LRWorkDim *= (2*MaxDimension+2*MinDimension+1);

  if (LRWorkDim > 1)
     RWorkDim = LRWorkDim;
  else
     RWorkDim = 1;

  double* TmpRWork = new double[RWorkDim];

  FORTRAN_NAME(zgesdd)(&Jobz, &this->NbrRow, &this->NbrColumn, TmpMatrix, &this->NbrRow, SigmaMatrix, TmpUMatrix, &SizeLDU, TmpVMatrix, &SizeLDVT, &TmpWorkingArea, &WorkingAreaSize, TmpRWork, &TmpIntegerWorkingArea, &Information); 

  WorkingAreaSize = (int) TmpWorkingArea.r;
  doublecomplex* WorkingArea = new doublecomplex [WorkingAreaSize];
  IntegerWorkingAreaSize = 8 * MinDimension;
  int* IntegerWorkingArea = new int [IntegerWorkingAreaSize];

  FORTRAN_NAME(zgesdd)(&Jobz, &this->NbrRow, &this->NbrColumn, TmpMatrix, &this->NbrRow, SigmaMatrix, TmpUMatrix, &SizeLDU, TmpVMatrix, &SizeLDVT, WorkingArea, &WorkingAreaSize, TmpRWork, IntegerWorkingArea, &Information);

  uMatrix = ComplexMatrix(TmpUMatrix, this->NbrRow, this->NbrRow, false);
  vMatrix = ComplexMatrix(TmpVMatrix, this->NbrColumn, this->NbrColumn, false);
  delete[] TmpUMatrix;
  delete[] TmpVMatrix;
  delete[] WorkingArea;
  delete[] IntegerWorkingArea;
  delete[] TmpRWork;
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

void ComplexMatrix::SingularValueDecomposition(ComplexMatrix& uMatrix, RealDiagonalMatrix& diagonal, ComplexMatrix& vMatrix, bool truncatedUVFlag)
{
  double* TmpDiag = this->SingularValueDecomposition(uMatrix, vMatrix, truncatedUVFlag);
  diagonal = RealDiagonalMatrix(TmpDiag, (uMatrix.NbrColumn>vMatrix.NbrRow)? vMatrix.NbrRow :uMatrix.NbrColumn);
}

// compute the diagonal part of the singular value decomposition U D V^t
// 
// return value = pointer on the diagonal elements of D

double* ComplexMatrix::SingularValueDecomposition()
{
#ifdef HAVE_LAPACK
  if ((this->NbrColumn == 0) || (this->NbrRow == 0))
    {
      return 0;
    }
  if ((this->NbrColumn == 1) || (this->NbrRow == 1))
    {
      double* SigmaMatrix = new double[1];
      SigmaMatrix[0] = 0.0;
      if (this->NbrColumn == 1)
	{
	  for (int i = 0; i < this->NbrRow; ++i)
           { 
	     SigmaMatrix[0] += SqrNorm(this->Columns[0][i]);
           }
	}
      else
	{
	  for (int i = 0; i < this->NbrColumn; ++i)
           {
	     SigmaMatrix[0] += SqrNorm(this->Columns[i][0]);
           }
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

  double* SigmaMatrix = new double[MinDimension];
  int Information = 0;
  int WorkingAreaSize = -1;
  int IntegerWorkingAreaSize = -1;
  doublecomplex TmpWorkingArea;
  int TmpIntegerWorkingArea;
  char Jobz = 'N';

  doublecomplex* TmpMatrix = new doublecomplex [this->NbrRow * this->NbrColumn];

  Complex* TmpColumn;
  for (int j = 0; j < this->NbrColumn; ++j)
    {
      TmpColumn = this->Columns[j].Components;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  TmpMatrix[i + (j * this->NbrRow)].r = TmpColumn[i].Re;
	  TmpMatrix[i + (j * this->NbrRow)].i = TmpColumn[i].Im;
	}
    }

  doublecomplex* TmpUMatrix = new doublecomplex [this->NbrColumn];
  doublecomplex* TmpVMatrix = new doublecomplex [this->NbrRow];

  int RWorkDim, LRWorkDim;

  LRWorkDim = 5 * MinDimension;

  if (LRWorkDim > 1)
     RWorkDim = LRWorkDim;
  else
     RWorkDim = 1;

  double* TmpRWork = new double[10 * RWorkDim];

  int DummySize = 1;
  FORTRAN_NAME(zgesdd)(&Jobz, &this->NbrRow, &this->NbrColumn, TmpMatrix, &this->NbrRow, SigmaMatrix, TmpUMatrix, &DummySize, TmpVMatrix, &DummySize, &TmpWorkingArea, &WorkingAreaSize, TmpRWork, &TmpIntegerWorkingArea, &Information);
  if (Information != 0)
    {
      cout << "warning Lapack zgesdd potential error while requesting workspace " << Information << endl;
    }
  WorkingAreaSize = ((int) TmpWorkingArea.r);
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
  doublecomplex* WorkingArea = new doublecomplex [WorkingAreaSize];
  IntegerWorkingAreaSize = 8 * MinDimension;
  int* IntegerWorkingArea = new int [IntegerWorkingAreaSize];
  FORTRAN_NAME(zgesdd)(&Jobz, &this->NbrRow, &this->NbrColumn, TmpMatrix, &this->NbrRow, SigmaMatrix, TmpUMatrix, &DummySize, TmpVMatrix, &DummySize, WorkingArea, &WorkingAreaSize, TmpRWork, IntegerWorkingArea, &Information);
  delete[] WorkingArea;
  delete[] TmpMatrix;
  delete[] TmpUMatrix;
  delete[] TmpVMatrix;
  delete[] TmpRWork;
  return SigmaMatrix;
#else
  return 0;
#endif
}

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



// Diagonalize a complex skew symmetric matrix using the LAPACK library (modifying current matrix)
//
// M = reference on real diagonal matrix of eigenvalues
// leftFlag = compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
// return value = reference on real matrix consisting of eigenvalues

ComplexDiagonalMatrix& ComplexMatrix::LapackDiagonalize (ComplexDiagonalMatrix& M, bool leftFlag)
{
  if (this->NbrColumn != this->NbrRow)
    return M;
  if (M.GetNbrColumn() != this->NbrColumn)
    M.Resize(this->NbrColumn, this->NbrColumn);
#ifdef __LAPACK__
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
  doublecomplex *complexWork = new doublecomplex[1];
  int lComplexWork=-1;
  // balancing omitted
  // workspace query
  FORTRAN_NAME(zgehrd)(&DimensionM, &iLow, &iHigh, TmpMatrix, &DimensionM, Tau,
		       complexWork, &lComplexWork, &Information);
  lComplexWork=(int)(complexWork[0].r);
  delete [] complexWork;
  complexWork = new doublecomplex[lComplexWork];
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
#else
  cout << "Warning, using ComplexMatrix::LapackDiagonalize without the lapack library" << endl;
#endif  
  return M;
}

// Diagonalize a complex skew symmetric matrix and evaluate transformation matrix using the LAPACK library (modifying current matrix)
//
// M = reference on real diagonal matrix of eigenvalues
// Q = matrix where transformation matrix has to be stored
// leftFlag = compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
// return value = reference on real matrix consisting of eigenvalues

ComplexDiagonalMatrix& ComplexMatrix::LapackDiagonalize (ComplexDiagonalMatrix& M, ComplexMatrix& Q, bool leftFlag)
{
  if (M.GetNbrColumn() != this->NbrColumn)
    M.Resize(this->NbrColumn, this->NbrColumn);
  if (Q.GetNbrColumn() != this->NbrColumn)
    Q.Resize(this->NbrColumn, this->NbrColumn);
#ifdef __LAPACK__
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
  doublecomplex* matrixA = new doublecomplex [((long) this->NbrColumn) * this->NbrRow];
  double* scale = new double [this->NbrRow];
  double* reciCondVal = new double [this->NbrRow];
  double* reciCondVec = new double [this->NbrRow];
  doublecomplex* Eigenvalues = new doublecomplex [2 * this->NbrRow];
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
  doublecomplex* EigenvectorsR = new doublecomplex [((long) this->NbrColumn) * this->NbrRow];
  doublecomplex* EigenvectorsL = new doublecomplex [((long) this->NbrColumn) * this->NbrRow];
  Complex* TmpColumn;
  for (int j = 0; j < this->NbrRow; ++j)
    {
      TmpColumn = this->Columns[j].Components;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  matrixA[i + (j * this->NbrRow)].r = TmpColumn[i].Re;
	  matrixA[i + (j * this->NbrRow)].i = TmpColumn[i].Im;
	}
    }  
  int Information = 0;
  int Dim = this->NbrRow;
  int iLow = 1;
  int iHigh = this->NbrRow;
  doublecomplex* complexWork = new doublecomplex[1];
  double* realWork = new double[2 * this->NbrRow];
  double absNorm;
  int lComplexWork = -1;
  // workspace query
  FORTRAN_NAME(zgeevx)("Balance", &JobVL, &JobVR, "No condition numbers",
		       &Dim, matrixA, &Dim, Eigenvalues, EigenvectorsL, &TmpLeadingLeftDimension, 
		       EigenvectorsR, &TmpLeadingRightDimension,
		       &iLow, &iHigh, scale, &absNorm, reciCondVal, reciCondVec, 
		       complexWork, &lComplexWork, realWork, &Information); 
  if (Information < 0)
    {
      cout << "Illegal argument " << -Information << " in LAPACK function call to zgeevx in ComplexMatrix.cc in LapackDiagonalize, line "<< __LINE__<<endl;
      exit(1);
    }
  lComplexWork=(int)(complexWork[0].r);
  delete[] complexWork;
  complexWork = new doublecomplex[lComplexWork];
  FORTRAN_NAME(zgeevx)("Balance", &JobVL, &JobVR, "No condition numbers",
		       &Dim, matrixA, &Dim, Eigenvalues, EigenvectorsL, &TmpLeadingLeftDimension, 
		       EigenvectorsR, &TmpLeadingRightDimension,
		       &iLow, &iHigh, scale, &absNorm, reciCondVal, reciCondVec, 
		       complexWork, &lComplexWork, realWork, &Information);
  if (Information < 0)
    {
      cout << "Illegal argument " << -Information << " in LAPACK function call to zgeevx in ComplexMatrix.cc in LapackDiagonalize, line "<< __LINE__<<endl;
      exit(1);
    }
  if (Information > 0)
    {
      cout << "Attention: Only part of eigenvalues converged in LAPACK function call to zgeevx in ComplexMatrix.cc in LapackDiagonalize, line "<< __LINE__<<endl;
    }
    
  // recover values of eigenvalues and eigenvectors
  if (leftFlag == true)
    {
      for (int j = 0;j < this->NbrRow; ++j)
	{
	  M.SetMatrixElement(j, j, Eigenvalues[j].r, Eigenvalues[j].i);
	  for (int i = 0; i < this->NbrRow; ++i)
	    Q.SetMatrixElement(i, j, EigenvectorsL[i + (j * NbrRow)].r, EigenvectorsL[i + (j * NbrRow)].i);
	}
    }
  else
    {
      for (int j = 0; j < this->NbrRow; ++j)
	{
	  M.SetMatrixElement(j, j, Eigenvalues[j].r, Eigenvalues[j].i);
	  for (int i = 0; i < this->NbrRow; ++i)
	    Q.SetMatrixElement(i, j, EigenvectorsR[i + (j * NbrRow)].r, EigenvectorsR[i + (j * NbrRow)].i);
	}
    }
  
  delete [] matrixA;
  delete [] Eigenvalues;
  delete [] EigenvectorsR;
  delete [] reciCondVal;
  delete [] reciCondVec;
  delete [] scale;
  delete [] complexWork;
  delete [] realWork;
    
#else
  cout << "Warning, using ComplexMatrix::LapackDiagonalize without the lapack library" << endl;
#endif  
  return M;
}


// Diagonalize a complex skew symmetric matrix and evaluate transformation matrix using the LAPACK library, truncating to single precision (modifying current matrix)
//
// M = reference on real diagonal matrix of eigenvalues
// Q = matrix where transformation matrix has to be stored
// leftFlag = compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
// return value = reference on real matrix consisting of eigenvalues

ComplexDiagonalMatrix& ComplexMatrix::LapackDiagonalizeSinglePrecision (ComplexDiagonalMatrix& M, ComplexMatrix& Q, bool leftFlag)
{
  if (M.GetNbrColumn() != this->NbrColumn)
    M.Resize(this->NbrColumn, this->NbrColumn);
  if (Q.GetNbrColumn() != this->NbrColumn)
    Q.Resize(this->NbrColumn, this->NbrColumn);
#ifdef __LAPACK__
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
  LAcomplex* matrixA = new LAcomplex [((long) this->NbrColumn) * this->NbrRow];
  double* scale = new double [this->NbrRow];
  double* reciCondVal = new double [this->NbrRow];
  double* reciCondVec = new double [this->NbrRow];
  LAcomplex* Eigenvalues = new LAcomplex [2 * this->NbrRow];
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
  LAcomplex* EigenvectorsR = new LAcomplex [((long) this->NbrColumn) * this->NbrRow];
  LAcomplex* EigenvectorsL = new LAcomplex [((long) this->NbrColumn) * this->NbrRow];
  Complex* TmpColumn;
  for (int j = 0; j < this->NbrRow; ++j)
    {
      TmpColumn = this->Columns[j].Components;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  matrixA[i + (j * this->NbrRow)].r = TmpColumn[i].Re; // involves a typecast
	  matrixA[i + (j * this->NbrRow)].i = TmpColumn[i].Im;
	}
    }  
  int Information = 0;
  int Dim = this->NbrRow;
  int iLow = 1;
  int iHigh = this->NbrRow;
  LAcomplex* complexWork = new LAcomplex[1];
  double* realWork = new double[2 * this->NbrRow];
  double absNorm;
  int lComplexWork = -1;
  // workspace query
  FORTRAN_NAME(cgeevx)("Balance", &JobVL, &JobVR, "No condition numbers",
		       &Dim, matrixA, &Dim, Eigenvalues, EigenvectorsL, &TmpLeadingLeftDimension, 
		       EigenvectorsR, &TmpLeadingRightDimension,
		       &iLow, &iHigh, scale, &absNorm, reciCondVal, reciCondVec, 
		       complexWork, &lComplexWork, realWork, &Information); 
  if (Information < 0)
    {
      cout << "Illegal argument " << -Information << " in LAPACK function call to zgeevx in ComplexMatrix.cc in LapackDiagonalize, line "<< __LINE__<<endl;
      exit(1);
    }
  lComplexWork=(int)(complexWork[0].r);
  delete[] complexWork;
  complexWork = new LAcomplex[lComplexWork];
  FORTRAN_NAME(cgeevx)("Balance", &JobVL, &JobVR, "No condition numbers",
		       &Dim, matrixA, &Dim, Eigenvalues, EigenvectorsL, &TmpLeadingLeftDimension, 
		       EigenvectorsR, &TmpLeadingRightDimension,
		       &iLow, &iHigh, scale, &absNorm, reciCondVal, reciCondVec, 
		       complexWork, &lComplexWork, realWork, &Information);
  if (Information < 0)
    {
      cout << "Illegal argument " << -Information << " in LAPACK function call to zgeevx in ComplexMatrix.cc in LapackDiagonalize, line "<< __LINE__<<endl;
      exit(1);
    }
  if (Information > 0)
    {
      cout << "Attention: Only part of eigenvalues converged in LAPACK function call to zgeevx in ComplexMatrix.cc in LapackDiagonalize, line "<< __LINE__<<endl;
    }
    
  // recover values of eigenvalues and eigenvectors
  if (leftFlag == true)
    {
      for (int j = 0;j < this->NbrRow; ++j)
	{
	  M.SetMatrixElement(j, j, Eigenvalues[j].r, Eigenvalues[j].i);
	  for (int i = 0; i < this->NbrRow; ++i)
	    Q.SetMatrixElement(i, j, EigenvectorsL[i + (j * NbrRow)].r, EigenvectorsL[i + (j * NbrRow)].i);
	}
    }
  else
    {
      for (int j = 0; j < this->NbrRow; ++j)
	{
	  M.SetMatrixElement(j, j, Eigenvalues[j].r, Eigenvalues[j].i);
	  for (int i = 0; i < this->NbrRow; ++i)
	    Q.SetMatrixElement(i, j, EigenvectorsR[i + (j * NbrRow)].r, EigenvectorsR[i + (j * NbrRow)].i);
	}
    }
  
  delete [] matrixA;
  delete [] Eigenvalues;
  delete [] EigenvectorsR;
  delete [] reciCondVal;
  delete [] reciCondVec;
  delete [] scale;
  delete [] complexWork;
  delete [] realWork;
    
#else
  cout << "Warning, using ComplexMatrix::LapackDiagonalize without the lapack library" << endl;
#endif  
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
#ifdef __LAPACK__
  doublecomplex* TmpMatrix = new doublecomplex [this->NbrRow * this->NbrRow];
  doublecomplex* TransQ = new doublecomplex [this->NbrRow * this->NbrRow];
  doublecomplex* Tau = new doublecomplex [this->NbrRow];
  doublecomplex* Eigenvalues = new doublecomplex [this->NbrRow];
  Complex *TmpColumn;
  for (int j = 0; j < NbrRow; ++j)
    {
      TmpColumn = this->Columns[j].Components;
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
  
#else
  cout << "Warning, using ComplexMatrix::LapackSchurForm without the lapack library" << endl;
#endif  
  return M;
}

// compute the LU decompostion of the matrix using the LAPACK library (conserving current matrix)
// 
// lowerMatrix = reference on the matrix where the lower triangular matrix will be stored
// upperMatrix = reference on the matrix where the upper triangular matrix will be stored
// return value = array that describes the additional row permutation

int* ComplexMatrix::LapackLUDecomposition(ComplexLowerTriangularMatrix& lowerMatrix, ComplexUpperTriangularMatrix& upperMatrix)
{
  if (this->NbrRow != this->NbrColumn)
    {
      cout << "LU decomposition is only performed for square matrices" << endl;
      return 0; 
    }
#ifdef __LAPACK__
  int* PermutationArray = new int [this->NbrRow];
  doublecomplex* TmpMatrix = new doublecomplex [((long) this->NbrRow) * ((long) this->NbrColumn)];
  long Pos = 0l;
  for (int j = 0; j < this->NbrColumn;++j)
    {
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  TmpMatrix[Pos].r = this->Columns[j][i].Re;
	  TmpMatrix[Pos].i = this->Columns[j][i].Im;
	  ++Pos;
	}
    }
  int Information = 0;
  int DimensionM = NbrRow;
  FORTRAN_NAME(zgetrf)(&DimensionM, &DimensionM, TmpMatrix, &DimensionM , PermutationArray, &Information);

  if (Information < 0)
    {
      cout << "Illegal argument " << -Information << " in LAPACK function call in ComplexMatrix.cc, line "<< __LINE__<<endl;
      exit(1);
    }

  
  lowerMatrix = ComplexLowerTriangularMatrix(this->NbrColumn, true);
  upperMatrix = ComplexUpperTriangularMatrix(this->NbrColumn, true);

  Complex Tmp;
  Pos = 0;
  for (int j = 0; j < this->NbrColumn;++j)
    {
      for (int i = 0; i < j; ++i)
	{
	  Tmp.Re = TmpMatrix[Pos].r;
	  Tmp.Im = TmpMatrix[Pos].i;
	  upperMatrix.SetMatrixElement(i ,j, Tmp);
	  ++Pos;
 	}
      lowerMatrix.SetMatrixElement(j ,j, 1.0);
      Tmp.Re = TmpMatrix[Pos].r;
      Tmp.Im = TmpMatrix[Pos].i;
      upperMatrix.SetMatrixElement(j ,j, Tmp);
      ++Pos;
     for (int i = j + 1 ; i < this->NbrRow; ++i)
	{
	  Tmp.Re = TmpMatrix[Pos].r;
	  Tmp.Im = TmpMatrix[Pos].i;
	  lowerMatrix.SetMatrixElement(i ,j, Tmp);
	  ++Pos;
	}
    }
  for (int i = 0; i < this->NbrRow; ++i)
    {
      PermutationArray[i]--;
    }
  delete [] TmpMatrix;
  return PermutationArray;
#else
  cout << "Warning, using ComplexMatrix::LapackLUDecomposition without the lapack library" << endl;
  return 0;
#endif  
}


void ComplexMatrix::QRDecompositionFromLapack (ComplexMatrix & Q, ComplexMatrix & R)
{
#ifdef __LAPACK__
  doublecomplex* TmpMatrix = new doublecomplex [((long) this->NbrRow) * ((long) this->NbrColumn)];

  long Pos = 0l;
  for (int j = 0; j < this->NbrColumn;++j)
    {
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  TmpMatrix[Pos].r = this->Columns[j][i].Re;
	  TmpMatrix[Pos].i = this->Columns[j][i].Im;
	  ++Pos;
	}
    }
  int DimensionRow = this->NbrRow;
  int DimensionColumn =  this->NbrColumn;
  doublecomplex * Tau = new doublecomplex [DimensionRow];
  int Information = 0;

  int LWork=-1;
  doublecomplex * IWork = new doublecomplex[DimensionRow] ;
  
  doublecomplex * complexWork = new doublecomplex[this->NbrRow * this->NbrColumn];
  int lComplexWork = this->NbrRow * this->NbrColumn;
  
  
  FORTRAN_NAME(zgeqrf)(&DimensionRow, &DimensionColumn, TmpMatrix, &DimensionRow, Tau, complexWork, &lComplexWork, &Information);

  Complex Tmp;
  for (int j = 0; j < this->NbrColumn;++j)
    {
      for (int i = 0; i < j; ++i)
	{
	  Tmp.Re = TmpMatrix[i+this->NbrRow *j].r;
	  Tmp.Im = TmpMatrix[i+this->NbrRow *j].i;
	  R.SetMatrixElement(i ,j, Tmp);
 	}
      Tmp.Re = TmpMatrix[j*(1+this->NbrRow)].r;
      Tmp.Im = TmpMatrix[j*(1+this->NbrRow)].i;
      R.SetMatrixElement(j ,j, Tmp);
    }

  
  int MinimumRowColumn = DimensionRow ;
  if (DimensionColumn <  MinimumRowColumn)
    MinimumRowColumn = DimensionColumn;
					 
  FORTRAN_NAME(zungqr) (&DimensionRow, &DimensionColumn, &MinimumRowColumn, TmpMatrix, &DimensionRow, Tau, complexWork, &lComplexWork, &Information); 
  
  delete[] complexWork;
  
  Pos = 0;
  for (int j = 0; j < this->NbrColumn;++j)
    {
      for (int i = 0; i < j; ++i)
	{
	  Tmp.Re = TmpMatrix[Pos].r;
	  Tmp.Im = TmpMatrix[Pos].i;
	  Q.SetMatrixElement(i ,j, Tmp);
	  ++Pos;
 	}
      Tmp.Re = TmpMatrix[Pos].r;
      Tmp.Im = TmpMatrix[Pos].i;
      Q.SetMatrixElement(j ,j, Tmp);
      ++Pos;
      for (int i = j + 1 ; i < this->NbrRow; ++i)
	{
	  Tmp.Re = TmpMatrix[Pos].r;
	  Tmp.Im = TmpMatrix[Pos].i;
	  Q.SetMatrixElement(i ,j, Tmp);
	  ++Pos;
	}
    }



  cout <<Q<< " " <<R <<endl;

  cout << Q*R<<endl;
  
  for(int i = 0 ; i <  R.NbrColumn ; i++)
    {
      R.GetMatrixElement(i ,i,Tmp);
      cout <<Tmp<<endl;
      if ( Tmp.Re < 0 )
	{
	  for(int j = 0 ; j <  R.NbrRow ; j++)
	    {
	      Q.GetMatrixElement(j ,i, Tmp);
	      Q.SetMatrixElement(j ,i, -1.0*Tmp);
	      R.GetMatrixElement(i ,j, Tmp);
	      R.SetMatrixElement(i ,j, -1.0*Tmp);
	    }
	}
    }

  cout <<Q<< " " <<R <<endl;
  cout << Q*R<<endl;
  
#else
  cout << "Warning, using ComplexMatrix::QRDecompositionFromLapack without the lapack library" << endl;
#endif  
}
 

 
// invert the current matrix using the LAPACK library
// 

void ComplexMatrix::LapackInvert()
{
  if (this->NbrRow != this->NbrColumn)
    {
      cout << "Matrix inversion is only performed for square matrices" << endl;
      return ; 
    }
#ifdef __LAPACK__
  int* PermutationArray = new int [this->NbrRow];
  doublecomplex* TmpMatrix = new doublecomplex [((long) this->NbrRow) * ((long) this->NbrColumn)];
  long Pos = 0l;
  for (int j = 0; j < this->NbrColumn;++j)
    {
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  TmpMatrix[Pos].r = this->Columns[j][i].Re;
	  TmpMatrix[Pos].i = this->Columns[j][i].Im;
	  ++Pos;
	}
    }
  int Information = 0;
  int DimensionM = NbrRow;
  FORTRAN_NAME(zgetrf)(&DimensionM, &DimensionM, TmpMatrix, &DimensionM , PermutationArray, &Information);

  if (Information < 0)
    {
      cout << "Illegal argument " << -Information << " in LAPACK function call in ComplexMatrix.cc, line "<< __LINE__<<endl;
      exit(1);
    }

  int WorkingAreaSize = -1;
  doublecomplex TmpWorkingArea;
  FORTRAN_NAME(zgetri)(&DimensionM, TmpMatrix, &DimensionM , PermutationArray, &TmpWorkingArea, &WorkingAreaSize, &Information);

  WorkingAreaSize = (int) TmpWorkingArea.r;
  doublecomplex* WorkingArea = new doublecomplex [WorkingAreaSize];
  FORTRAN_NAME(zgetri)(&DimensionM, TmpMatrix, &DimensionM , PermutationArray, WorkingArea, &WorkingAreaSize, &Information); 

  Pos = 0l;
  for (int j = 0; j < this->NbrColumn;++j)
    {
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  this->Columns[j][i].Re = TmpMatrix[Pos].r;
	  this->Columns[j][i].Im = TmpMatrix[Pos].i;
	  ++Pos;
	}
    }
  delete [] TmpMatrix;
  delete [] PermutationArray;
#else
  cout << "Warning, using ComplexMatrix::LapackInvert without the lapack library" << endl;
#endif  
}

// build a random unitary matrix
//
// return value = reference on the current matrix

ComplexMatrix& ComplexMatrix::RandomUnitaryMatrix()
{

  for (int j = 0; j < this->NbrColumn; ++j)
    {
      for (int i = 0; i < this->NbrRow;++i)
	{
	  this->Columns[j][i] = drand48() * Phase (2.0 * M_PI * drand48());
	}
      double TmpNorm = this->Columns[j].Norm();
      this->Columns[j] /= TmpNorm;    
    }
  this->OrthoNormalizeColumns();
  return *this;
}

 // compute the Frobenius norm of the current matrix 
 //
 // return value = value of the norm

double ComplexMatrix::FrobeniusNorm()
{
  double Tmp = 0.0;
  for(int i =0; i < this->NbrColumn; i++)
    Tmp += this->Columns[i].SqrNorm();
  return sqrt(Tmp);
}

// compute the Frobenius scalar product of two matrices
//
// return value = value of the scalar product
Complex FrobeniusScalarProduct(ComplexMatrix & matrixA, ComplexMatrix & matrixB)
{
  Complex Tmp = 0.0;
  for(int i =0 ; i < matrixA.NbrRow; i++)
    {
      for(int j =0 ; j < matrixA.NbrRow; j++)
	{
	  Tmp+= Conj(matrixA.Columns[j][i]) *  matrixB.Columns[j][i];
	}
    }
  return Tmp;
}

#ifdef __MPI__

// send a matrix to a given MPI process
// 
// communicator = reference on the communicator to use
// id = id of the destination MPI process
// return value = reference on the current matrix

Matrix& ComplexMatrix::SendMatrix(MPI::Intracomm& communicator, int id)
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

Matrix& ComplexMatrix::BroadcastMatrix(MPI::Intracomm& communicator,  int id)
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

Matrix& ComplexMatrix::ReceiveMatrix(MPI::Intracomm& communicator, int id)
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

Matrix& ComplexMatrix::SumMatrix(MPI::Intracomm& communicator, int id)
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

Matrix& ComplexMatrix::ReassembleMatrix(MPI::Intracomm& communicator, int id)
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

Matrix* ComplexMatrix::BroadcastClone(MPI::Intracomm& communicator, int id)
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
      return new ComplexMatrix(communicator, id);
    }
  return 0;
}

// create a new matrix on each MPI node with same size and same type but non-initialized components
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the matrix
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new matrix 

Matrix* ComplexMatrix::BroadcastEmptyClone(MPI::Intracomm& communicator, int id, bool zeroFlag)
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
      return new ComplexMatrix(communicator, id);
    }
  return 0;
}

#endif

