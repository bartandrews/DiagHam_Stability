////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of real band-diagonal symmetric matrix                //
//                                                                            //
//                        last modification : 16/03/2005                      //
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


#include "Matrix/BandDiagonalHermitianMatrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealMatrix.h"

#include <cstdlib>
#include <cmath>

#ifdef __LAPACK__
// binding to the LAPACK zhpev function
extern "C" void FORTRAN_NAME(zhpev)(const char* jobz, const char* uplo, const int* dimension, const doublecomplex* matrix, const double *eigenvalues, const doublecomplex *eigenvectors, const int* leadingDimension, const doublecomplex *work, const doublereal *rwork, const int* information );
extern "C" void FORTRAN_NAME(zhbev)(const char* jobz, const char* uplo, const int* dimension, const int* numdiag, const doublecomplex* bandmatrixAD, const int* leadingDimensionAB, const double *eigenvalues, const doublecomplex *eigenvectorsZ, const int* leadingDimensionZ, const doublecomplex *work, const doublereal *rwork, const int* information);
#endif


using std::cout;
using std::endl;


// default constructor
//

BandDiagonalHermitianMatrix::BandDiagonalHermitianMatrix() 
{
  this->DiagonalElements = 0;
  this->RealUpperOffDiagonalElements = 0;
  this->ImaginaryUpperOffDiagonalElements = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->NbrBands = 0; 
  this->TrueNbrBands = this->NbrBands;
  this->TrueNbrRow = 0;
  this->TrueNbrColumn = 0;
  this->MatrixType = Matrix::ComplexElements | Matrix::Hermitian | Matrix::BandDiagonal;
  this->Dummy = 0.0;
#ifdef __LAPACK__
  this->LapackWorkAreaDimension=0;
#endif

}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

BandDiagonalHermitianMatrix::BandDiagonalHermitianMatrix(int dimension, int nbrBands, bool zero)
{
  this->NbrBands = nbrBands;
  this->TrueNbrBands = this->NbrBands;
  this->DiagonalElements = new double [dimension];
  this->RealUpperOffDiagonalElements = new double* [this->NbrBands];
  this->ImaginaryUpperOffDiagonalElements = new double* [this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->RealUpperOffDiagonalElements[i] = new double [dimension];
      this->ImaginaryUpperOffDiagonalElements[i] = new double [dimension];
    }
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = dimension;
  this->TrueNbrColumn = dimension;
  this->MatrixType = Matrix::ComplexElements | Matrix::Hermitian | Matrix::BandDiagonal;
  if (zero == true)
    {
      for (int i = 0; i < this->NbrRow; ++i)
	this->DiagonalElements[i] = 0.0;
      double* Tmp, *TmpI;
      for (int j = 0; j < this->NbrBands; ++j)
	{
	  Tmp = this->RealUpperOffDiagonalElements[j];
	  TmpI = this->ImaginaryUpperOffDiagonalElements[j];
	  for (int i = 0; i < this->NbrRow; ++i)
	    {
	      Tmp[i] = 0.0;
	      TmpI[i] = 0.0;
	    }
	}
    }
  this->Dummy = 0.0;
#ifdef __LAPACK__
  this->LapackWorkAreaDimension=0;
#endif
}

// constructor from matrix elements (without duplicating datas)
//
// diagonal = pointer to diagonal element array
// upperOffDiagonal = pointer to the array which contains upper off-diagonal elements (first index is used as row index)
// dimension = matrix dimension
// nbrBands = number of bands in the upper part of the matrix

BandDiagonalHermitianMatrix::BandDiagonalHermitianMatrix(double* diagonal, double** realUpperDiagonal, double ** imaginaryUpperDiagonal, int dimension, int nbrBands) 
{
  this->DiagonalElements = diagonal;
  this->RealUpperOffDiagonalElements = realUpperDiagonal;
  this->ImaginaryUpperOffDiagonalElements = imaginaryUpperDiagonal;
  this->NbrBands = nbrBands;
  this->TrueNbrBands = this->NbrBands;
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = dimension;
  this->TrueNbrColumn = dimension;
  this->MatrixType = Matrix::ComplexElements | Matrix::Hermitian | Matrix::BandDiagonal;
  this->Dummy = 0.0;
#ifdef __LAPACK__
  this->LapackWorkAreaDimension=0;
#endif

}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

BandDiagonalHermitianMatrix::BandDiagonalHermitianMatrix(const BandDiagonalHermitianMatrix& M) 
{  
  this->DiagonalElements = M.DiagonalElements;
  this->RealUpperOffDiagonalElements = M.RealUpperOffDiagonalElements;
  this->ImaginaryUpperOffDiagonalElements = M.ImaginaryUpperOffDiagonalElements;
  this->NbrBands = M.NbrBands;
  this->TrueNbrBands = M.NbrBands;
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Hermitian | Matrix::BandDiagonal;
  this->Dummy = 0.0;
#ifdef __LAPACK__
  this->LapackWorkAreaDimension=0;
#endif
}

// destructor
//

BandDiagonalHermitianMatrix::~BandDiagonalHermitianMatrix() 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < this->TrueNbrBands; ++i)
	{
	  delete[] this->RealUpperOffDiagonalElements[i];
	  delete[] this->ImaginaryUpperOffDiagonalElements[i];
	}
      delete[] this->RealUpperOffDiagonalElements;
      delete[] this->ImaginaryUpperOffDiagonalElements;
      delete[] this->DiagonalElements;

#ifdef __LAPACK__
      if (this->LapackWorkAreaDimension>0)
	{
	  delete[] this->LapackMatrix;
	  if (this->LapackEVMatrix!=NULL)
	    delete[] this->LapackEVMatrix;
	  delete[] this->LapackWorkingArea;
	  delete[] this->LapackRealWorkingArea;
	}
#endif      
    }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

BandDiagonalHermitianMatrix& BandDiagonalHermitianMatrix::operator = (const BandDiagonalHermitianMatrix& M) 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < this->TrueNbrBands; ++i)
	{
	  delete[] this->RealUpperOffDiagonalElements[i];
	  delete[] this->ImaginaryUpperOffDiagonalElements[i];
	}
      delete[] this->RealUpperOffDiagonalElements;
      delete[] this->ImaginaryUpperOffDiagonalElements;
      delete[] this->DiagonalElements;
    }
  this->DiagonalElements = M.DiagonalElements;
  this->RealUpperOffDiagonalElements = M.RealUpperOffDiagonalElements;
  this->ImaginaryUpperOffDiagonalElements = M.ImaginaryUpperOffDiagonalElements;
  this->NbrBands = M.NbrBands;
  this->TrueNbrBands = M.NbrBands;
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Hermitian | Matrix::BandDiagonal;
  this->Dummy = 0.0;
#ifdef __LAPACK__
  this->LapackWorkAreaDimension=0;
#endif
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* BandDiagonalHermitianMatrix::Clone ()
{
  return ((Matrix*) new BandDiagonalHermitianMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void BandDiagonalHermitianMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i == j) && (i < this->NbrRow))
    {
      this->DiagonalElements[i] = x;
    }
  else
    {
      if (j > i)
	{
	  j -= i;
	  if ((j <= this->NbrBands) && (i < (this->NbrRow - 1)))
	    {
	      this->RealUpperOffDiagonalElements[j - 1][i] = x;
	      this->ImaginaryUpperOffDiagonalElements[j - 1][i] = 0.0;
	    }
	}
      else
	{
	  i -= j;
	  if ((i <= this->NbrBands) && (j < (this->NbrRow - 1)))
	    {
	      this->RealUpperOffDiagonalElements[i - 1][j] = x;
	      this->ImaginaryUpperOffDiagonalElements[i - 1][j] = 0.0;
	    }	
	}
    }    
}

// get a matrix element (real part if complex)
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

void BandDiagonalHermitianMatrix::GetMatrixElement(int i, int j, double& x) const
{
  if ((i == j) && (i < this->NbrRow))
    {
      x = this->DiagonalElements[i];
      return;
    }
  else
    {
      if (j > i)
	{
	  j -= i;
	  if ((j <= this->NbrBands) && (i < (this->NbrRow - 1)))
	    {
	      x = this->RealUpperOffDiagonalElements[j - 1][i];
	      return;
	    }
	}
      else
	{
	  i -= j;
	  if ((i <= this->NbrBands) && (j < (this->NbrRow - 1)))
	    {
	      x = this->RealUpperOffDiagonalElements[i - 1][j];
	      return;
	    }	
	}
    }    
  x = 0.0;
}

// get a matrix element (real part if complex)
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

void BandDiagonalHermitianMatrix::GetMatrixElement(int i, int j, Complex& x) const
{
  if ((i == j) && (i < this->NbrRow))
    {
      x.Re = this->DiagonalElements[i];
      x.Im = 0.0;
      return;
    }
  else
    {
      if (j > i)
	{
	  j -= i;
	  if ((j <= this->NbrBands) && (i < (this->NbrRow - 1)))
	    {
	      x.Re = this->RealUpperOffDiagonalElements[j - 1][i];
	      x.Im = this->ImaginaryUpperOffDiagonalElements[j - 1][i];
	      return;
	    }
	}
      else
	{
	  i -= j;
	  if ((i <= this->NbrBands) && (j < (this->NbrRow - 1)))
	    {
	      x.Re = this->RealUpperOffDiagonalElements[i - 1][j];
	      x.Im = -this->ImaginaryUpperOffDiagonalElements[i - 1][j];
	      return;
	    }	
	}
    }    
  x = 0.0;
}

// return refernce on real part of a given matrix element
//
// i = line position
// j = column position
// return value = reference on real part

double& BandDiagonalHermitianMatrix::operator () (int i, int j)
{
  if ((i == j) && (i < this->NbrRow))
    {
      return this->DiagonalElements[i];
    }
  else
    {
      if (j > i)
	{
	  j -= i;
	  if ((j <= this->NbrBands) && (i < (this->NbrRow - 1)))
	    {
	      return this->RealUpperOffDiagonalElements[j - 1][i];
	    }
	}
      else
	{
	  i -= j;
	  if ((i <= this->NbrBands) && (j < (this->NbrRow - 1)))
	    {
	      return this->RealUpperOffDiagonalElements[i - 1][j];
	    }	
	}
    }    
  return this->Dummy;
}

// get a matrix element 
// 
// i = Row number
// j = Column number
// return value = matrix element M_(i,j)

Complex BandDiagonalHermitianMatrix::GetElement(int i, int j)
{
  if ((i == j) && (i < this->NbrRow))
    {
      return Complex(this->DiagonalElements[i]);
    }
  else
    {
      if (j > i)
	{
	  j -= i;
	  if ((j <= this->NbrBands) && (i < (this->NbrRow - 1)))
	    {
	      return Complex(this->RealUpperOffDiagonalElements[j - 1][i], this->ImaginaryUpperOffDiagonalElements[j - 1][i]);
	    }
	}
      else
	{
	  i -= j;
	  if ((i <= this->NbrBands) && (j < (this->NbrRow - 1)))
	    {
	      return Complex(this->RealUpperOffDiagonalElements[i - 1][j], -this->ImaginaryUpperOffDiagonalElements[i - 1][j]);
	    }	
	}
    }    
  return 0.0;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void BandDiagonalHermitianMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  if ((i == j) && (i < this->NbrRow))
    {
      this->DiagonalElements[i] = x.Re;
    }
  else
    {
      if (j > i)
	{
	  j -= i;
	  if ((j <= this->NbrBands) && (i < (this->NbrRow - 1)))
	    {
	      this->RealUpperOffDiagonalElements[j - 1][i] = x.Re;
	      this->ImaginaryUpperOffDiagonalElements[j - 1][i] = x.Im;
	    }
	}
      else
	{
	  i -= j;
	  if ((i <= this->NbrBands) && (j < (this->NbrRow - 1)))
	    {
	      this->RealUpperOffDiagonalElements[i - 1][j] = x.Re;
	      this->ImaginaryUpperOffDiagonalElements[i - 1][j] = -x.Im;
	    }	
	}
    }    

}

// access to i-th diagonal element
// 
// i = position 
// return value = reference on i-th diagonal element

double& BandDiagonalHermitianMatrix::DiagonalElement(int i)
{
  return this->DiagonalElements[i];
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void BandDiagonalHermitianMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i == j) && (i < this->NbrRow))
    {
      this->DiagonalElements[i] = x;
    }
  else
    {
      if (j > i)
	{
	  j -= i;
	  if ((j <= this->NbrBands) && (i < (this->NbrRow - 1)))
	    {
	      this->RealUpperOffDiagonalElements[j - 1][i] += x;
	    }
	}
      else
	{
	  i -= j;
	  if ((i <= this->NbrBands) && (j < (this->NbrRow - 1)))
	    {
	      this->RealUpperOffDiagonalElements[i - 1][j] += x;
	    }	
	}
    }    
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void BandDiagonalHermitianMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  if ((i == j) && (i < this->NbrRow))
    {
      this->DiagonalElements[i] = x.Re;
    }
  else
    {
      if (j > i)
	{
	  j -= i;
	  if ((j <= this->NbrBands) && (i < (this->NbrRow - 1)))
	    {
	      this->RealUpperOffDiagonalElements[j - 1][i] += x.Re;
	      this->ImaginaryUpperOffDiagonalElements[j - 1][i] += x.Im;
	    }
	}
      else
	{
	  i -= j;
	  if ((i <= this->NbrBands) && (j < (this->NbrRow - 1)))
	    {
	      this->RealUpperOffDiagonalElements[i - 1][j] += x.Re;
	      this->ImaginaryUpperOffDiagonalElements[j - 1][i] += -x.Im;
	    }	
	}
    }    

}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void BandDiagonalHermitianMatrix::Resize (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  if (nbrRow <= this->TrueNbrRow)
    {
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      return;
    }
  double* TmpDiag = new double [nbrRow];
  double** TmpUpperDiagRe = new double* [this->TrueNbrBands];
  double** TmpUpperDiagIm = new double* [this->TrueNbrBands];
  for (int i = 0; i < this->TrueNbrBands; ++i)
    {
      TmpUpperDiagRe[i] = new double [nbrRow];
      TmpUpperDiagIm[i] = new double [nbrRow];
    }
  if (this->Flag.Used() == true)
    {
      for (int i = 0; i < this->NbrRow; i++)
	TmpDiag[i] = this->DiagonalElements[i];
      for (int j = 0; j < this->NbrBands; ++j)
	for (int i = 0; i < this->NbrRow; ++i)
	  {
	    TmpUpperDiagRe[j][i] = this->RealUpperOffDiagonalElements[j][i];
	    TmpUpperDiagIm[j][i] = this->ImaginaryUpperOffDiagonalElements[j][i];
	  }
    }
   if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int j = 0; j < this->TrueNbrBands; ++j)
	{
	  delete[] this->RealUpperOffDiagonalElements[j];
	  delete[] this->ImaginaryUpperOffDiagonalElements[j];
	}
      delete[] this->RealUpperOffDiagonalElements;
      delete[] this->ImaginaryUpperOffDiagonalElements;
      delete[] this->DiagonalElements;
    }
  this->DiagonalElements = TmpDiag;
  this->RealUpperOffDiagonalElements = TmpUpperDiagRe;
  this->ImaginaryUpperOffDiagonalElements = TmpUpperDiagRe;
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = nbrRow;
  this->TrueNbrColumn = nbrColumn;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
  return;
}

// Resize matrix and change the number of bands
//
// nbrRow = new number of rows
// nbrColumn = new number of columns
// nbrBands = new number of bands

void BandDiagonalHermitianMatrix::Resize (int nbrRow, int nbrColumn, int nbrBands)
{
  if (nbrRow != nbrColumn)
    return;
  if ((nbrRow <= this->TrueNbrRow) && (nbrBands <= this->TrueNbrBands))
    {
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      this->NbrBands = nbrBands;
      return;
    }
  double* TmpDiag = new double [nbrRow];
  double** TmpUpperDiagRe = new double* [nbrBands];
  double** TmpUpperDiagIm = new double* [nbrBands];
  for (int i = 0; i < nbrBands; ++i)
    {
      TmpUpperDiagRe[i] = new double [nbrRow];
      TmpUpperDiagIm[i] = new double [nbrRow];
    }    
  if (this->Flag.Used() == true)
    {
      for (int i = 0; i < this->NbrRow; i++)
	TmpDiag[i] = this->DiagonalElements[i];
      for (int j = 0; j < this->NbrBands; ++j)
	for (int i = 0; i < this->NbrRow; ++i)
	  {
	    TmpUpperDiagRe[j][i] = this->RealUpperOffDiagonalElements[j][i];
	    TmpUpperDiagIm[j][i] = this->ImaginaryUpperOffDiagonalElements[j][i];
	  }
    }
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int j = 0; j < this->TrueNbrBands; ++j)
	{
	  delete[] this->RealUpperOffDiagonalElements[j];
	  delete[] this->ImaginaryUpperOffDiagonalElements[j];
	}
      delete[] this->RealUpperOffDiagonalElements;
      delete[] this->ImaginaryUpperOffDiagonalElements;
      delete[] this->DiagonalElements;
    }
  this->DiagonalElements = TmpDiag;
  this->RealUpperOffDiagonalElements = TmpUpperDiagRe;
  this->ImaginaryUpperOffDiagonalElements = TmpUpperDiagRe;
  this->NbrBands = nbrBands;
  this->TrueNbrBands = nbrBands;
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = nbrRow;
  this->TrueNbrColumn = nbrColumn;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
  return;
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void BandDiagonalHermitianMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  if (nbrRow <= this->TrueNbrRow)
    {
      for (int i = this->NbrRow; i < nbrRow; i++)
	{
	  this->DiagonalElements[i] = 0.0;
	}
      for (int j = 0; j < this->NbrBands; ++j)
	for (int i = this->NbrRow; i < nbrRow; ++i)
	  {
	    this->RealUpperOffDiagonalElements[j][i] = 0.0;
	    this->ImaginaryUpperOffDiagonalElements[j][i] = 0.0;
	  }
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      return;
    }
  double* TmpDiag = new double [nbrRow];
  double** TmpUpperDiagRe = new double* [this->NbrBands];
  double** TmpUpperDiagIm = new double* [this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      TmpUpperDiagRe[i] = new double [nbrRow];
      TmpUpperDiagIm[i] = new double [nbrRow];
    }    
  if (this->Flag.Used() == true)
    {
      int i = 0;
      for (; i < this->NbrRow; i++)
	TmpDiag[i] = this->DiagonalElements[i];
      for (; i < nbrRow; i++)
	TmpDiag[i] = 0.0;
      for (int j = 0; j < this->NbrBands; ++j)
	for (i = 0; i < this->NbrRow; ++i)
	  {
	    TmpUpperDiagRe[j][i] = this->RealUpperOffDiagonalElements[j][i];
	    TmpUpperDiagIm[j][i] = this->ImaginaryUpperOffDiagonalElements[j][i];
	  }
      for (int j = 0; j < this->NbrBands; ++j)
	for (i = nbrRow; i < this->NbrRow; ++i)
	  {
	    TmpUpperDiagRe[j][i] = 0.0;
	    TmpUpperDiagIm[j][i] = 0.0;
	  }
    }
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int j = 0; j < this->NbrBands; ++j)
	{
	  delete[] this->RealUpperOffDiagonalElements[j];
	  delete[] this->ImaginaryUpperOffDiagonalElements[j];
	}
      delete[] this->RealUpperOffDiagonalElements;
      delete[] this->ImaginaryUpperOffDiagonalElements;
      delete[] this->DiagonalElements;
    }
  this->DiagonalElements = TmpDiag;
  this->RealUpperOffDiagonalElements = TmpUpperDiagRe;
  this->ImaginaryUpperOffDiagonalElements = TmpUpperDiagRe;
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = nbrRow;
  this->TrueNbrColumn = nbrColumn;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
  return;
}

// copy matrix
//
// M = matrix to copy
// return value = refence on current matrix

BandDiagonalHermitianMatrix& BandDiagonalHermitianMatrix::Copy (BandDiagonalHermitianMatrix& M)
{
  if (this->NbrBands != M.NbrBands)
    this->Resize(M.NbrRow, M.NbrColumn, M.NbrBands);
  else
    if (this->NbrRow != M.NbrRow)
      this->Resize(M.NbrRow, M.NbrColumn);
  for (int i = 0; i < M.NbrColumn; i++)
    this->DiagonalElements[i] = M.DiagonalElements[i];
  double* Tmp1;
  double* Tmp2;
  for (int j = 0; j < this->NbrBands; ++j)
    {
      Tmp1 = this->RealUpperOffDiagonalElements[j];
      Tmp2 = M.RealUpperOffDiagonalElements[j];
      for (int i = 0; i < this->NbrRow; i++)
	Tmp1[i] = Tmp2[i];
      Tmp1 = this->ImaginaryUpperOffDiagonalElements[j];
      Tmp2 = M.ImaginaryUpperOffDiagonalElements[j];
      for (int i = 0; i < this->NbrRow; i++)
	Tmp1[i] = Tmp2[i];

    }
  return *this;
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

BandDiagonalHermitianMatrix operator + (const BandDiagonalHermitianMatrix& M1, const BandDiagonalHermitianMatrix& M2) 
{
  if (M1.NbrRow != M2.NbrRow)
    return BandDiagonalHermitianMatrix();
  int MaxNbrBands = M1.NbrBands;
  int MinNbrBands = M2.NbrBands;
  if (M1.NbrBands < M2.NbrBands)
    {
      MaxNbrBands = M2.NbrBands;
      MinNbrBands = M1.NbrBands;
    }
  double* Diagonal = new double [M1.NbrRow];
  double** UpperDiagonalRe = new double* [MaxNbrBands];
  double** UpperDiagonalIm = new double* [MaxNbrBands];
  for (int i = 0; i < MaxNbrBands; ++i)
    {
      UpperDiagonalRe[i] = new double [M1.NbrRow];
      UpperDiagonalIm[i] = new double [M1.NbrRow];
    }
  for (int i = 0; i < M1.NbrRow; i++)
    Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
  double* Tmp1;
  double* Tmp2;
  double* Tmp3;
  for (int j = 0; j < MinNbrBands; ++j)
    {      
      Tmp1 = M1.RealUpperOffDiagonalElements[j];
      Tmp2 = M2.RealUpperOffDiagonalElements[j];
      Tmp3 = UpperDiagonalRe[j];
      for (int i = 0; i < M1.NbrRow; i++)
	Tmp3[i] = Tmp1[i] + Tmp2[i];
      Tmp1 = M1.ImaginaryUpperOffDiagonalElements[j];
      Tmp2 = M2.ImaginaryUpperOffDiagonalElements[j];
      Tmp3 = UpperDiagonalIm[j];
      for (int i = 0; i < M1.NbrRow; i++)
	Tmp3[i] = Tmp1[i] + Tmp2[i];
    }
  for (int j = MinNbrBands; j < MaxNbrBands; ++j)
    {      
      if (M1.NbrBands == MaxNbrBands)
	Tmp1 = M1.RealUpperOffDiagonalElements[j];
      else
	Tmp1 = M2.RealUpperOffDiagonalElements[j];
      Tmp3 = UpperDiagonalRe[j];
      for (int i = 0; i < M1.NbrRow; i++)
	Tmp3[i] = Tmp1[i];
      if (M1.NbrBands == MaxNbrBands)
	Tmp1 = M1.ImaginaryUpperOffDiagonalElements[j];
      else
	Tmp1 = M2.ImaginaryUpperOffDiagonalElements[j];
      Tmp3 = UpperDiagonalIm[j];
      for (int i = 0; i < M1.NbrRow; i++)
	Tmp3[i] = Tmp1[i];

    }
  return BandDiagonalHermitianMatrix(Diagonal, UpperDiagonalRe, UpperDiagonalIm, M1.NbrRow, MaxNbrBands);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

BandDiagonalHermitianMatrix operator - (const BandDiagonalHermitianMatrix& M1, const BandDiagonalHermitianMatrix& M2) 
{
  if (M1.NbrRow != M2.NbrRow)
    return BandDiagonalHermitianMatrix();
  int MaxNbrBands = M1.NbrBands;
  int MinNbrBands = M2.NbrBands;
  if (M1.NbrBands < M2.NbrBands)
    {
      MaxNbrBands = M2.NbrBands;
      MinNbrBands = M1.NbrBands;
    }
  double* Diagonal = new double [M1.NbrRow];
  double** UpperDiagonalRe = new double* [MaxNbrBands];
  double** UpperDiagonalIm = new double* [MaxNbrBands];
  for (int i = 0; i < MaxNbrBands; ++i)
    {
      UpperDiagonalRe[i] = new double [M1.NbrRow];
      UpperDiagonalIm[i] = new double [M1.NbrRow];
    }
  for (int i = 0; i < M1.NbrRow; i++)
    Diagonal[i] = M1.DiagonalElements[i] - M2.DiagonalElements[i];
  double* Tmp1;
  double* Tmp2;
  double* Tmp3;
  for (int j = 0; j < MinNbrBands; ++j)
    {      
      Tmp1 = M1.RealUpperOffDiagonalElements[j];
      Tmp2 = M2.RealUpperOffDiagonalElements[j];
      Tmp3 = UpperDiagonalRe[j];
      for (int i = 0; i < M1.NbrRow; i++)
	Tmp3[i] = Tmp1[i] + Tmp2[i];
      Tmp1 = M1.ImaginaryUpperOffDiagonalElements[j];
      Tmp2 = M2.ImaginaryUpperOffDiagonalElements[j];
      Tmp3 = UpperDiagonalIm[j];
      for (int i = 0; i < M1.NbrRow; i++)
	Tmp3[i] = Tmp1[i] + Tmp2[i];

    }
  if (M1.NbrBands == MaxNbrBands)
    for (int j = MinNbrBands; j < MaxNbrBands; ++j)
      {      
	Tmp1 = M1.RealUpperOffDiagonalElements[j];
	Tmp3 = UpperDiagonalRe[j];
	for (int i = 0; i < M1.NbrRow; i++)
	  Tmp3[i] = Tmp1[i];
	Tmp1 = M1.ImaginaryUpperOffDiagonalElements[j];
	Tmp3 = UpperDiagonalIm[j];
	for (int i = 0; i < M1.NbrRow; i++)
	  Tmp3[i] = Tmp1[i];
      }
  else
    for (int j = MinNbrBands; j < MaxNbrBands; ++j)
      {      
	Tmp1 = M2.RealUpperOffDiagonalElements[j];
	Tmp3 = UpperDiagonalRe[j];
	for (int i = 0; i < M1.NbrRow; i++)
	  Tmp3[i] = -Tmp1[i];
	Tmp1 = M2.ImaginaryUpperOffDiagonalElements[j];
	Tmp3 = UpperDiagonalIm[j];
	for (int i = 0; i < M1.NbrRow; i++)
	  Tmp3[i] = -Tmp1[i];

      }
  return BandDiagonalHermitianMatrix(Diagonal, UpperDiagonalRe, UpperDiagonalIm, M1.NbrRow, MaxNbrBands);
}

// multiply a matrix with a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

BandDiagonalHermitianMatrix operator * (const BandDiagonalHermitianMatrix& M, double x) 
{
  double* Diagonal = new double [M.NbrRow];
  double** UpperDiagonalRe = new double* [M.NbrBands];
  double** UpperDiagonalIm = new double* [M.NbrBands];
  for (int i = 0; i < M.NbrRow; ++i)
    Diagonal[i] = M.DiagonalElements[i] * x;      
  for (int i = 0; i < M.NbrBands; ++i)
    {
      UpperDiagonalRe[i] = new double [M.NbrRow];
      for (int j = 0; j < M.NbrRow; ++j)
	UpperDiagonalRe[i][j] = M.RealUpperOffDiagonalElements[i][j] * x;
      UpperDiagonalIm[i] = new double [M.NbrRow];
      for (int j = 0; j < M.NbrRow; ++j)
	UpperDiagonalIm[i][j] = M.ImaginaryUpperOffDiagonalElements[i][j] * x;      	

    }
  return BandDiagonalHermitianMatrix(Diagonal, UpperDiagonalRe, UpperDiagonalIm, M.NbrRow, M.NbrBands);
}

// multiply a matrix with a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

BandDiagonalHermitianMatrix operator * (double x, const BandDiagonalHermitianMatrix& M) 
{
  double* Diagonal = new double [M.NbrRow];
  double** UpperDiagonalRe = new double* [M.NbrBands];
  double** UpperDiagonalIm = new double* [M.NbrBands];
  for (int i = 0; i < M.NbrRow; ++i)
    Diagonal[i] = M.DiagonalElements[i] * x;      
  for (int i = 0; i < M.NbrBands; ++i)
    {
      UpperDiagonalRe[i] = new double [M.NbrRow];
      for (int j = 0; j < M.NbrRow; ++j)
	UpperDiagonalRe[i][j] = M.RealUpperOffDiagonalElements[i][j] * x;
      UpperDiagonalIm[i] = new double [M.NbrRow];
      for (int j = 0; j < M.NbrRow; ++j)
	UpperDiagonalIm[i][j] = M.ImaginaryUpperOffDiagonalElements[i][j] * x;      	
    }
  return BandDiagonalHermitianMatrix(Diagonal, UpperDiagonalRe, UpperDiagonalIm, M.NbrRow, M.NbrBands);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

BandDiagonalHermitianMatrix operator / (const BandDiagonalHermitianMatrix& M, double x) 
{
  x = 1.0 / x;
  double* Diagonal = new double [M.NbrRow];
  double** UpperDiagonalRe = new double* [M.NbrBands];
  double** UpperDiagonalIm = new double* [M.NbrBands];
  for (int i = 0; i < M.NbrRow; ++i)
    Diagonal[i] = M.DiagonalElements[i] * x;      
  for (int i = 0; i < M.NbrBands; ++i)
    {
      UpperDiagonalRe[i] = new double [M.NbrRow];
      for (int j = 0; j < M.NbrRow; ++j)
	UpperDiagonalRe[i][j] = M.RealUpperOffDiagonalElements[i][j] * x;
      UpperDiagonalIm[i] = new double [M.NbrRow];
      for (int j = 0; j < M.NbrRow; ++j)
	UpperDiagonalIm[i][j] = M.ImaginaryUpperOffDiagonalElements[i][j] * x;
    }
  return BandDiagonalHermitianMatrix(Diagonal, UpperDiagonalRe, UpperDiagonalIm, M.NbrRow, M.NbrBands);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

BandDiagonalHermitianMatrix& BandDiagonalHermitianMatrix::operator += (const BandDiagonalHermitianMatrix& M) 
{
  if ((this->NbrRow != M.NbrRow) || (this->NbrBands < M.NbrBands))
    return *this;
  for (int i = 0; i < this->NbrRow; ++i)
    this->DiagonalElements[i] += M.DiagonalElements[i];
  for (int j = 0; j < M.NbrBands; ++j)
    for (int i = 0; i < this->NbrRow; ++i)
      {
	this->RealUpperOffDiagonalElements[j][i] += M.RealUpperOffDiagonalElements[j][i];
	this->ImaginaryUpperOffDiagonalElements[j][i] += M.ImaginaryUpperOffDiagonalElements[j][i];
      }
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

BandDiagonalHermitianMatrix& BandDiagonalHermitianMatrix::operator -= (const BandDiagonalHermitianMatrix& M) 
{
  if ((this->NbrRow != M.NbrRow) || (this->NbrBands < M.NbrBands))
    return *this;
  for (int i = 0; i < this->NbrRow; ++i)
    this->DiagonalElements[i] -= M.DiagonalElements[i];
  for (int j = 0; j < M.NbrBands; ++j)
    for (int i = 0; i < this->NbrRow; ++i)
      {
	this->RealUpperOffDiagonalElements[j][i] -= M.RealUpperOffDiagonalElements[j][i];
	this->ImaginaryUpperOffDiagonalElements[j][i] -= M.ImaginaryUpperOffDiagonalElements[j][i];
      }
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

BandDiagonalHermitianMatrix& BandDiagonalHermitianMatrix::operator *= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  for (int i = 0; i < this->NbrRow; i++)
    this->DiagonalElements[i] *= x;
  for (int j = 0; j < this->NbrBands; ++j)
    for (int i = 0; i < this->NbrRow; i++)
      {
	this->RealUpperOffDiagonalElements[i][j] *= x;
	this->ImaginaryUpperOffDiagonalElements[i][j] *= x;
      }
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

BandDiagonalHermitianMatrix& BandDiagonalHermitianMatrix::operator /= (double x) 
{  
  if (this->NbrRow == 0)
    return *this;
  x = 1.0 / x;
  for (int i = 0; i < this->NbrRow; i++)
    this->DiagonalElements[i] *= x;
  for (int j = 0; j < this->NbrBands; ++j)
    for (int i = 0; i < this->NbrRow; i++)    
      {
	this->RealUpperOffDiagonalElements[i][j] *= x;
	this->ImaginaryUpperOffDiagonalElements[i][j] *= x;
      }
  return *this;
}

// evaluate matrix trace
//
// return value = matrix trace 

double BandDiagonalHermitianMatrix::Tr () 
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

double BandDiagonalHermitianMatrix::Det () 
{
  return 1.0;
}

// Tridiagonalize a real band symmetric matrix using Rutishauser-Schwarz (modifying current matrix)
//
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// err = absolute error on matrix element
// return value = reference on real tridiagonal symmetric matrix

RealTriDiagonalSymmetricMatrix& BandDiagonalHermitianMatrix::Tridiagonalize (RealTriDiagonalSymmetricMatrix& M, double err)
{
  cout << "DiagHam native routine BandDiagonalHermitianMatrix::Tridiagonalize not defined, please use LAPACK version."<<endl;
  exit(1);
  /*
  if (M.NbrRow != this->NbrRow)
    M.Resize(this->NbrRow, this->NbrColumn);
  int ReducedNbrRow = this->NbrRow - 1;
  int ReducedNbrRow2 = this->NbrRow - 2;
  double TmpNorm;
  double Cosinus;
  double Sinus;
  double Tmp;
  double Tmp2;
  double FillInElement;
  int MinJ;
  int Pos;
  int Pos2;
  int Max;
  int GivenColumnPosition;
  int GivenRowPosition;
  double SquareErr = err * err;
  if (SquareErr < MACHINE_PRECISION)
    SquareErr = MACHINE_PRECISION;

  for (int i = 0; i < ReducedNbrRow2; ++i)
    {
      MinJ = ReducedNbrRow - i - 1;
      if (MinJ >= this->NbrBands)
	MinJ = this->NbrBands - 1;
      for (int j = MinJ; j >= 1; --j)
	{
	  GivenRowPosition = i;
	  GivenColumnPosition = j - 1;
	  FillInElement = this->UpperOffDiagonalElements[GivenColumnPosition + 1][GivenRowPosition]; 
	  Cosinus = this->UpperOffDiagonalElements[GivenColumnPosition][GivenRowPosition];
	  Sinus = -FillInElement;
	  TmpNorm =  ((Cosinus * Cosinus) + (Sinus * Sinus));
	  while ((GivenRowPosition < this->NbrRow) && (TmpNorm > SquareErr))
	    {
	      // zeroing outmost element of the i-th line using Given rotation and apllying Given rotations to chase out the fill-in element produced by the previous Given rotation
	      TmpNorm = sqrt (TmpNorm);
	      this->UpperOffDiagonalElements[GivenColumnPosition][GivenRowPosition] = TmpNorm;
	      TmpNorm = 1.0 / TmpNorm;
	      Cosinus *= TmpNorm;
	      Sinus *= TmpNorm;
	      Pos = GivenRowPosition + 1;
	      Pos2 = GivenColumnPosition;
	      while (Pos2 > 0)
		{
		  Tmp = this->UpperOffDiagonalElements[Pos2][Pos];
		  this->UpperOffDiagonalElements[Pos2][Pos] *= Cosinus;
		  this->UpperOffDiagonalElements[Pos2][Pos] += Sinus * this->UpperOffDiagonalElements[Pos2 - 1][Pos];
		  this->UpperOffDiagonalElements[Pos2 - 1][Pos] *= Cosinus;
		  this->UpperOffDiagonalElements[Pos2 - 1][Pos] -= Sinus * Tmp;
		  --Pos2;
		  ++Pos;
		}
	      

	      Tmp = this->UpperOffDiagonalElements[0][Pos];
	      Tmp2 = this->DiagonalElements[Pos];
	      this->DiagonalElements[Pos] *= Cosinus * Cosinus;
	      this->DiagonalElements[Pos] += Sinus * ((Sinus * this->DiagonalElements[Pos + 1]) - (2.0 * Cosinus * Tmp));
	      this->UpperOffDiagonalElements[0][Pos] *= ((Cosinus * Cosinus) - (Sinus * Sinus));
	      this->UpperOffDiagonalElements[0][Pos] += Cosinus * Sinus * (Tmp2 - this->DiagonalElements[Pos + 1]);
	      this->DiagonalElements[Pos + 1] *= Cosinus * Cosinus;
	      this->DiagonalElements[Pos + 1] += Sinus * ((Sinus * Tmp2) + (2.0 * Cosinus * Tmp));
	      
	      Pos2 = 1;
	      Max = ReducedNbrRow - Pos;
	      if (Max > this->NbrBands)
		Max = this->NbrBands;
	      while (Pos2 < Max)
		{
		  Tmp = this->UpperOffDiagonalElements[Pos2][Pos];
		  this->UpperOffDiagonalElements[Pos2][Pos] *= Cosinus;
		  this->UpperOffDiagonalElements[Pos2][Pos] -= Sinus * this->UpperOffDiagonalElements[Pos2 - 1][Pos + 1];
		  this->UpperOffDiagonalElements[Pos2 - 1][Pos + 1] *= Cosinus;
		  this->UpperOffDiagonalElements[Pos2 - 1][Pos + 1] += Sinus * Tmp;
		  ++Pos2;
		}
	      
	      if ((ReducedNbrRow - Pos) > this->NbrBands)
		{
		  FillInElement = -this->UpperOffDiagonalElements[Max - 1][Pos + 1] * Sinus;
		  this->UpperOffDiagonalElements[Max - 1][Pos + 1] *= Cosinus;
		  GivenRowPosition = Pos;
		  GivenColumnPosition = this->NbrBands - 1;
		  Cosinus = this->UpperOffDiagonalElements[GivenColumnPosition][GivenRowPosition];
		  Sinus = -FillInElement;
		  TmpNorm =  ((Cosinus * Cosinus) + (Sinus * Sinus));
		}
	      else
		GivenRowPosition = this->NbrRow;
	    }
	}
    }

  for (int i = 0; i < this->NbrRow; ++i)
    M.DiagonalElements[i] = this->DiagonalElements[i];
  ++ReducedNbrRow;
  double* TmpColumn = this->UpperOffDiagonalElements[0];
  for (int i = 0; i < ReducedNbrRow; ++i)
    M.UpperDiagonalElements[i] = TmpColumn[i];
  */
  return M;
}

// Tridiagonalize a real band symmetric matrix using Rutishauser-Schwarz and evaluate transformation matrix  (modifying current matrix)
//
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// err = absolute error on matrix element
// Q = matrix where transformation matrix has to be stored
// return value = reference on real tridiagonal symmetric matrix

RealTriDiagonalSymmetricMatrix& BandDiagonalHermitianMatrix::Tridiagonalize (RealTriDiagonalSymmetricMatrix& M, double err, ComplexMatrix& Q)
{
  cout << "DiagHam native routine BandDiagonalHermitianMatrix::Tridiagonalize not defined, please use LAPACK version."<<endl;
  exit(1);
  /*
  if (M.NbrRow != this->NbrRow)
    M.Resize(this->NbrRow, this->NbrColumn);
  if ((Q.GetNbrRow() != this->NbrRow) || (Q.GetNbrColumn() != this->NbrColumn))
    Q.Resize(this->NbrRow, this->NbrColumn);
  for (int i = 0; i < this->NbrRow; i++)
    {
      RealVector& Vector1 =  Q[i];
      for (int j = 0; j < i; j++)
	Vector1[j] = 0.0;       
      Vector1[i] = 1.0;
      for (int j = i + 1; j < this->NbrColumn; j++)
	Vector1[j] = 0.0;       
    }
  int ReducedNbrRow = this->NbrRow - 1;
  int ReducedNbrRow2 = this->NbrRow - 2;
  double TmpNorm;
  double Cosinus;
  double Sinus;
  double Tmp;
  double Tmp2;
  double FillInElement;
  int MinJ;
  int Pos;
  int Pos2;
  int Max;
  int GivenColumnPosition;
  int GivenRowPosition;
  double SquareErr = err * err;
  if (SquareErr < MACHINE_PRECISION)
    SquareErr = MACHINE_PRECISION;

  for (int i = 0; i < ReducedNbrRow2; ++i)
    {
      MinJ = ReducedNbrRow - i - 1;
      if (MinJ >= this->NbrBands)
        MinJ = this->NbrBands - 1;
      for (int j = MinJ; j >= 1; --j)
        {
          GivenRowPosition = i;
          GivenColumnPosition = j - 1;
          FillInElement = this->UpperOffDiagonalElements[GivenColumnPosition + 1][GivenRowPosition];
          Cosinus = this->UpperOffDiagonalElements[GivenColumnPosition][GivenRowPosition];
          Sinus = -FillInElement;
          TmpNorm =  ((Cosinus * Cosinus) + (Sinus * Sinus));
          while ((GivenRowPosition < this->NbrRow) && (TmpNorm > SquareErr))
            {
              // zeroing outmost element of the i-th line using Given rotation and apllying Given rotations to chase out the fill-in element produced by the previous Given rotation
              TmpNorm = sqrt (TmpNorm);
              this->UpperOffDiagonalElements[GivenColumnPosition][GivenRowPosition] = TmpNorm;
              TmpNorm = 1.0 / TmpNorm;
              Cosinus *= TmpNorm;
              Sinus *= TmpNorm;
              Pos = GivenRowPosition + 1;
              Pos2 = GivenColumnPosition;
              while (Pos2 > 0)
                {
                  Tmp = this->UpperOffDiagonalElements[Pos2][Pos];
                  this->UpperOffDiagonalElements[Pos2][Pos] *= Cosinus;
                  this->UpperOffDiagonalElements[Pos2][Pos] += Sinus * this->UpperOffDiagonalElements[Pos2 - 1][Pos];
                  this->UpperOffDiagonalElements[Pos2 - 1][Pos] *= Cosinus;
                  this->UpperOffDiagonalElements[Pos2 - 1][Pos] -= Sinus * Tmp;
                  --Pos2;
                  ++Pos;
                }

              Tmp = this->UpperOffDiagonalElements[0][Pos];
              Tmp2 = this->DiagonalElements[Pos];
              this->DiagonalElements[Pos] *= Cosinus * Cosinus;
              this->DiagonalElements[Pos] += Sinus * ((Sinus * this->DiagonalElements[Pos + 1]) - (2.0 * Cosinus * Tmp));
              this->UpperOffDiagonalElements[0][Pos] *= ((Cosinus * Cosinus) - (Sinus * Sinus));
              this->UpperOffDiagonalElements[0][Pos] += Cosinus * Sinus * (Tmp2 - this->DiagonalElements[Pos + 1]);
              this->DiagonalElements[Pos + 1] *= Cosinus * Cosinus;
              this->DiagonalElements[Pos + 1] += Sinus * ((Sinus * Tmp2) + (2.0 * Cosinus * Tmp));

              Pos2 = 1;
              Max = ReducedNbrRow - Pos;
              if (Max > this->NbrBands)
                Max = this->NbrBands;
              while (Pos2 < Max)
                {
                  Tmp = this->UpperOffDiagonalElements[Pos2][Pos];
                  this->UpperOffDiagonalElements[Pos2][Pos] *= Cosinus;
                  this->UpperOffDiagonalElements[Pos2][Pos] -= Sinus * this->UpperOffDiagonalElements[Pos2 - 1][Pos + 1];
                  this->UpperOffDiagonalElements[Pos2 - 1][Pos + 1] *= Cosinus;
                  this->UpperOffDiagonalElements[Pos2 - 1][Pos + 1] += Sinus * Tmp;
                  ++Pos2;
                }
	      
	      
	      RealVector& Vector1 =  Q[GivenRowPosition + GivenColumnPosition + 1];
	      RealVector& Vector2 =  Q[GivenRowPosition + GivenColumnPosition + 2];
	      for (int k = 0; k < this->NbrRow; ++k)
		{
		  Tmp = Vector2[k];
		  Vector2[k] *= Cosinus;
		  Vector2[k] += Sinus * Vector1[k];
		  Vector1[k] *= Cosinus;
		  Vector1[k] -= Sinus * Tmp;
		}
	      
             if ((ReducedNbrRow - Pos) > this->NbrBands)
                {
                  FillInElement = -this->UpperOffDiagonalElements[Max - 1][Pos + 1] * Sinus;
                  this->UpperOffDiagonalElements[Max - 1][Pos + 1] *= Cosinus;
                  GivenRowPosition = Pos;
                  GivenColumnPosition = this->NbrBands - 1;
                  Cosinus = this->UpperOffDiagonalElements[GivenColumnPosition][GivenRowPosition];
                  Sinus = -FillInElement;
                  TmpNorm =  ((Cosinus * Cosinus) + (Sinus * Sinus));
                }
	     else
                GivenRowPosition = this->NbrRow;
            }
        }
    }

  for (int i = 0; i < this->NbrRow; ++i)
    M.DiagonalElements[i] = this->DiagonalElements[i];
  ++ReducedNbrRow;
  double* TmpColumn = this->UpperOffDiagonalElements[0];
  for (int i = 0; i < ReducedNbrRow; ++i)
    M.UpperDiagonalElements[i] = TmpColumn[i];
  */
  return M;
}

#ifdef __LAPACK__

// Diagonalize a real symmetric matrix using the LAPACK library (modifying current matrix)
//
// M = reference on real diagonal matrix where result has to be stored
// err = absolute error on matrix element
// maxIter = maximum number of iteration to fund an eigenvalue
// return value = reference on real tridiagonal symmetric matrix

RealDiagonalMatrix& BandDiagonalHermitianMatrix::LapackDiagonalize (RealDiagonalMatrix& M, double err, int maxIter)
{
  // lazy butt method of doing it: full complex matrix (slow!) - should implement call to ZHBEV
  // to call: extern "C" void FORTRAN_NAME(zhbev)(const char* jobz, const char* uplo, const int* dimension, const int* numdiag, const doublecomplex* bandmatrixAD, const int* leadingDimensionAB, const double *eigenvalues, const doublecomplex *eigenvectorsZ, const int* leadingDimensionZ, const doublecomplex *work, const doublereal *rwork, const int* information);

  if (M.GetNbrRow() != this->NbrRow)
    M.Resize(this->NbrRow, this->NbrColumn);
  if (this->NbrRow>500)
    cout << "Attention, using slow diagonalization routine in BandDiagonalHermitianMatrix::LapackDiagonalize!"<<endl;
  Complex Tmp;
  if (this->LapackWorkAreaDimension<this->NbrRow)
    {
      if (this->LapackWorkAreaDimension>0)
	{
	  delete [] LapackMatrix;
	  delete [] LapackWorkingArea;
	  delete [] LapackRealWorkingArea;
	  if (LapackEVMatrix!=0) delete [] LapackEVMatrix;	  
	}      
      this->LapackMatrix = new doublecomplex [this->NbrRow * (this->NbrRow+1)/2];
      this->LapackEVMatrix = NULL;	  
      this->LapackWorkingArea = new doublecomplex [2*this->NbrRow-1];
      this->LapackRealWorkingArea = new double [3*this->NbrRow-2];
      this->LapackWorkAreaDimension=this->NbrRow;
    }
  
  int Information = 0;  
  const char* Jobz = "N";
  const char* UpperLower = "U";
  int TotalIndex = 0;
  for (int j = 0; j < this->NbrRow; ++j)
    {
      for (int i = 0; i < j; ++i)
	{
	  this->GetMatrixElement(i,j,Tmp);
	  LapackMatrix[TotalIndex].r = Tmp.Re;
	  LapackMatrix[TotalIndex].i = Tmp.Im;
	  ++TotalIndex;
	}
      LapackMatrix[TotalIndex].r = this->DiagonalElements[j];
      LapackMatrix[TotalIndex].i = 0.0;
      ++TotalIndex;      
    }
  FORTRAN_NAME(zhpev)(Jobz, UpperLower, &this->NbrRow, LapackMatrix, M.DiagonalElements, LapackEVMatrix, &this->NbrRow, LapackWorkingArea, LapackRealWorkingArea, &Information);  
  return M;
}

// Diagonalize a real symmetric matrix and evaluate transformation matrix using the LAPACK library (modifying current matrix)
//
// M = reference on real diagonal matrix where result has to be stored
// Q = matrix where transformation matrix has to be stored
// err = absolute error on matrix element
// maxIter = maximum number of iteration to fund an eigenvalue
// return value = reference on real tridiagonal symmetric matrix

RealDiagonalMatrix& BandDiagonalHermitianMatrix::LapackDiagonalize (RealDiagonalMatrix& M, ComplexMatrix& Q, double err, int maxIter)
{
  // Lazy butt method : full matrix diagonalization
  if (this->NbrRow>500)
    cout << "Attention, using slow diagonalization routine in BandDiagonalHermitianMatrix::LapackDiagonalize!"<<endl;
  if (M.GetNbrRow() != this->NbrRow)
    M.Resize(this->NbrRow, this->NbrColumn);
  if (Q.GetNbrRow() != this->NbrRow)
    Q.Resize(this->NbrRow, this->NbrColumn);
  Complex Tmp;
  if (this->LapackWorkAreaDimension<this->NbrRow)
    {
      if (this->LapackWorkAreaDimension>0)
	{
	  delete [] LapackMatrix;
	  delete [] LapackWorkingArea;
	  delete [] LapackRealWorkingArea;
	  if (LapackEVMatrix!=NULL)
	    delete [] LapackEVMatrix;	  
	}
      this->LapackMatrix = new doublecomplex [this->NbrRow * (this->NbrRow+1)/2];
      this->LapackEVMatrix = NULL;
      this->LapackWorkingArea = new doublecomplex [2*this->NbrRow-1];
      this->LapackRealWorkingArea = new double [3*this->NbrRow-2];
      this->LapackWorkAreaDimension=this->NbrRow;
    }
  if (LapackEVMatrix==NULL)
    LapackEVMatrix = new doublecomplex[this->NbrRow * this->NbrRow];
  int Information = 0;  
  char Jobz = 'V';
  char UpperLower = 'U';
  int TotalIndex = 0;
  for (int j = 0; j < this->NbrRow; ++j)
    {
      for (int i = 0; i < j; ++i)
	{
	  this->GetMatrixElement(i,j,Tmp);
	  LapackMatrix[TotalIndex].r = Tmp.Re;
	  LapackMatrix[TotalIndex].i = Tmp.Im;
	  ++TotalIndex;
	}
      LapackMatrix[TotalIndex].r = this->DiagonalElements[j];
      LapackMatrix[TotalIndex].i = 0.0;
      ++TotalIndex;      
    }
  FORTRAN_NAME(zhpev)(&Jobz, &UpperLower, &this->NbrRow, LapackMatrix, M.DiagonalElements, LapackEVMatrix, &this->NbrRow, LapackWorkingArea, LapackRealWorkingArea, &Information);
  
  TotalIndex=0;
  for (int i = 0; i < this->NbrRow; ++i)
    for (int j = 0; j < this->NbrRow; ++j)
      {
	Tmp.Re = LapackEVMatrix[TotalIndex].r;
	Tmp.Im = -LapackEVMatrix[TotalIndex].i;
	Q.SetMatrixElement(j, i, Tmp);
	++TotalIndex;
      }
  return M;
}

#endif

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const BandDiagonalHermitianMatrix& P)
{
  for (int i = 0; i < P.NbrRow; ++i)
    {
      int j = 0;
      for (; j < (i - P.NbrBands); ++j)
	Str << "0    ";
      for (; j < i; ++j)
	{
	  Str << P.RealUpperOffDiagonalElements[i - j - 1][j];
	  if (P.ImaginaryUpperOffDiagonalElements[i - j - 1][j] > 0.0)
	    Str << -P.ImaginaryUpperOffDiagonalElements[i - j - 1][j] << "i    ";
	  else
	    if (P.ImaginaryUpperOffDiagonalElements[i - j - 1][j] != 0.0)
	      Str <<  "+" << -P.ImaginaryUpperOffDiagonalElements[i - j - 1][j] << "i    ";
	    else
	      Str << "    ";
	}
      Str << P.DiagonalElements[i] << "    ";
      ++j;
      for (; ((j <= (P.NbrBands + i)) && (j < P.NbrColumn)); ++j)
	{
	  Str << P.RealUpperOffDiagonalElements[j - i - 1][i];
	  if (P.ImaginaryUpperOffDiagonalElements[j - i - 1][i] < 0.0)
	    Str << P.ImaginaryUpperOffDiagonalElements[j - i - 1][i] <<"i    ";
	  else
	    if (P.ImaginaryUpperOffDiagonalElements[j - i - 1][i] != 0.0)
	      Str << "+" << P.ImaginaryUpperOffDiagonalElements[j - i - 1][i] <<"i    ";
	    else
	      Str << "    ";
	}
      for (; j < P.NbrColumn; j++)
	Str << "0    ";
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

MathematicaOutput& operator << (MathematicaOutput& Str, const BandDiagonalHermitianMatrix& P)
{
  Str << "{";
  for (int i = 0; i < P.NbrRow; i++)
    {
      Str << "{";
      int j = 0;
      for (; j < (i - P.NbrBands); ++j)
	Str << "0,";
      for (; j < i; ++j)
	{
	  Str << P.RealUpperOffDiagonalElements[i - j - 1][j];
	  if (P.ImaginaryUpperOffDiagonalElements[i - j - 1][j] < 0.0)
	    Str << P.ImaginaryUpperOffDiagonalElements[i - j - 1][j] << "I, ";
	  else
	    if (P.ImaginaryUpperOffDiagonalElements[i - j - 1][j] != 0.0)
	      Str << "+" << P.ImaginaryUpperOffDiagonalElements[i - j - 1][j] << "I, ";
	    else
	      Str << ", ";

	}
      Str << P.DiagonalElements[i] << ",";
      j++;

      for (; ((j < (P.NbrBands + i)) && (j < P.NbrColumn)); ++j)
	{

	  Str << P.RealUpperOffDiagonalElements[j - i - 1][i];
	  if (P.ImaginaryUpperOffDiagonalElements[j - i - 1][i] < 0.0)
	    Str << P.ImaginaryUpperOffDiagonalElements[j - i - 1][i] << "I";
	  else
	    if (P.ImaginaryUpperOffDiagonalElements[j - i - 1][i] != 0.0)
	      Str << "+" << P.ImaginaryUpperOffDiagonalElements[j - i - 1][i] << "I";
	  if (j != (P.NbrColumn - 1))
	    Str << ",";
	}
      for (; j < (P.NbrColumn - 1); j++)
	Str << "0,";
      if (j == (P.NbrColumn - 1))
	Str << "0";
      Str << "}";
      if (i != (P.NbrRow - 1))
	Str << ",";
    }
  Str << "}";
  return Str;
}

#endif
