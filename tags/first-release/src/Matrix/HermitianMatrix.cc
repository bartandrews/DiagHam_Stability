////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of hermitian matrix                        //
//                                                                            //
//                        last modification : 23/01/2001                      //
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


#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/ComplexVector.h"

#include <stdlib.h>


using std::endl;


// default constructor
//

HermitianMatrix::HermitianMatrix() 
{
  this->DiagonalElements = 0;
  this->RealOffDiagonalElements = 0;
  this->ImaginaryOffDiagonalElements = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = 0;
  this->MatrixType = Matrix::ComplexElements | Matrix::Hermitian;
}

// constructor from matrix elements (without duplicating datas)
//
// diagonal = pointer to diagonal element array
// realOffDiagonal = pointer to real part of off-diagonal elements
// imaginaryOffDiagonal = pointer to imaginary part of off-diagonal elements
// dimension = matrix dimension

HermitianMatrix::HermitianMatrix(double* diagonal, double* realOffDiagonal, double* imaginaryOffDiagonal, int dimension) 
{
  this->DiagonalElements = diagonal;
  this->RealOffDiagonalElements = realOffDiagonal;
  this->ImaginaryOffDiagonalElements = imaginaryOffDiagonal;
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::ComplexElements | Matrix::Hermitian;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

HermitianMatrix::HermitianMatrix(const HermitianMatrix& M) 
{
  this->DiagonalElements = M.DiagonalElements;
  this->Flag = M.Flag;
  this->RealOffDiagonalElements = M.RealOffDiagonalElements;
  this->ImaginaryOffDiagonalElements = M.ImaginaryOffDiagonalElements;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::ComplexElements | Matrix::Hermitian;
}

// copy constructor from a real tridiagonal symmetric matrix (without duplicating diagonal elements)
//
// M = matrix to copy

HermitianMatrix::HermitianMatrix(const RealTriDiagonalSymmetricMatrix& M) 
{
}

// destructor
//

HermitianMatrix::~HermitianMatrix() 
{
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) 
      && (this->DiagonalElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->DiagonalElements;
	delete[] this->RealOffDiagonalElements;
	delete[] this->ImaginaryOffDiagonalElements;
      }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

HermitianMatrix& HermitianMatrix::operator = (const HermitianMatrix& M) 
{
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) 
      && (this->DiagonalElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->DiagonalElements;
	delete[] this->RealOffDiagonalElements;
	delete[] this->ImaginaryOffDiagonalElements;
      }
  this->Flag = M.Flag;
  this->DiagonalElements = M.DiagonalElements;
  this->RealOffDiagonalElements = M.RealOffDiagonalElements;
  this->ImaginaryOffDiagonalElements = M.ImaginaryOffDiagonalElements;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  return *this;
}

// assignement from  a real tridiagonal symmetric matrix (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix  

HermitianMatrix& HermitianMatrix::operator = (const RealTriDiagonalSymmetricMatrix& M) 
{
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* HermitianMatrix::Clone ()
{
  return ((Matrix*) new HermitianMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void HermitianMatrix::SetMatrixElement(int i, int j, double x)
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
      j -= i * (i - 2 * this->NbrRow - this->Increment + 3) / 2 + 1;
      this->RealOffDiagonalElements[j] = x;
      this->ImaginaryOffDiagonalElements[j] = 0.0;	        
    }
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void HermitianMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  if ((i == j) || (i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  else
    {
      if (i > j)
	{
	  int tmp = j;
	  j = i;
	  i = tmp;
	  j -= i * (i - 2 * this->NbrRow - this->Increment + 3) / 2 + 1;
	  this->RealOffDiagonalElements[j] = x.Re;
	  this->ImaginaryOffDiagonalElements[j] = -x.Im;	        
	}
      else
	{
	  j -= i * (i - 2 * this->NbrRow - this->Increment + 3) / 2 + 1;
	  this->RealOffDiagonalElements[j] = x.Re;
	  this->ImaginaryOffDiagonalElements[j] = x.Im;	        
	}
    }
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void HermitianMatrix::AddToMatrixElement(int i, int j, double x)
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
      j -= i * (i - 2 * this->NbrRow - this->Increment + 3) / 2 + 1;
      this->RealOffDiagonalElements[j] += x;
    }
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void HermitianMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  if ((i == j) || (i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  else
    {
      if (i > j)
	{
	  int tmp = j;
	  j = i;
	  i = tmp;
	  j -= i * (i - 2 * this->NbrRow - this->Increment + 3) / 2 + 1;
	  this->RealOffDiagonalElements[j] += x.Re;
	  this->ImaginaryOffDiagonalElements[j] -= x.Im;	        
	}
      else
	{
	  j -= i * (i - 2 * this->NbrRow - this->Increment + 3) / 2 + 1;
	  this->RealOffDiagonalElements[j] += x.Re;
	  this->ImaginaryOffDiagonalElements[j] += x.Im;	        
	}
    }
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void HermitianMatrix::Resize (int nbrRow, int nbrColumn)
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
  int Tot = (nbrRow * (nbrRow - 1)) / 2;
  double* TmpRealOffDiag = new double [Tot];
  double* TmpImaginaryOffDiag = new double [Tot];
  for (int i = 0; i < this->NbrRow; i++)
    TmpDiag [i] = this->DiagonalElements[i];
  for (int i = this->NbrRow; i < nbrRow; i++)
    TmpDiag [i]  = 0.0;
  int k = 0;
  int l = 0;
  for (int i = 0; i < (this->NbrRow - 1); ++i)
    {
      for (int j = 0; j < (this->NbrRow - 1 - i); ++j)
	{
	  TmpRealOffDiag[k] = this->RealOffDiagonalElements[l];
	  TmpImaginaryOffDiag[k] = this->ImaginaryOffDiagonalElements[l];
	  ++l;
	  ++k;
	}
      for (int j = this->NbrRow; j < nbrRow; j++)
	{
	  TmpRealOffDiag[k] = 0.0;
	  TmpImaginaryOffDiag[k] = 0.0;
	  ++k;
	}      
    }
  for (int i = (this->NbrRow * (this->NbrRow - 1)) / 2; i < Tot; ++i)
    {
      TmpRealOffDiag[i] = 0.0;
      TmpImaginaryOffDiag[i] = 0.0;
    }
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) 
      && (this->DiagonalElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->DiagonalElements;
	delete[] this->RealOffDiagonalElements;
	delete[] this->ImaginaryOffDiagonalElements;
      }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->Flag.Initialize();
  this->DiagonalElements = TmpDiag;
  this->RealOffDiagonalElements = TmpRealOffDiag;
  this->ImaginaryOffDiagonalElements = TmpImaginaryOffDiag;
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void HermitianMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  if (nbrRow <= this->TrueNbrRow)
    {
      if (this->NbrRow < nbrRow)
	{
	  int Tot = (nbrRow * (nbrRow - 1)) / 2;
	  for (int i = this->NbrRow; i < nbrRow; ++i)
	    this->DiagonalElements[i] = 0.0;
	  int k = 2 * (this->NbrRow - 1);
	  for (int i = 0; i < (this->NbrRow - 1); ++i)
	    {
	      for (int j = this->NbrRow; j < nbrRow; ++j)
		{
		  this->RealOffDiagonalElements[k] = 0.0;
		  this->ImaginaryOffDiagonalElements[k] = 0.0;
		  ++k;
		}
	      k += (this->NbrRow - 2 - i);
	    }
	  for (int i = (this->NbrRow * (this->NbrRow - 1)) / 2; i < Tot; i++)
	    {
	      this->RealOffDiagonalElements[i] = 0.0;
	      this->ImaginaryOffDiagonalElements[i] = 0.0;
	    }
	}
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      this->Increment = (this->TrueNbrRow - this->NbrRow);
      return;
    }
  double* TmpDiag = new double [nbrRow];
  int Tot = (nbrRow * (nbrRow - 1)) / 2;
  double* TmpRealOffDiag = new double [Tot];
  double* TmpImaginaryOffDiag = new double [Tot];
  for (int i = 0; i < this->NbrRow; i++)
    TmpDiag [i] = this->DiagonalElements[i];
  for (int i = this->NbrRow; i < nbrRow; i++)
    TmpDiag [i]  = 0.0;
  int k = 0;
  int l = 0;
  for (int i = 0; i < (this->NbrRow - 1); i++)
    {
      for (int j = 0; j < (this->NbrRow - 1 - i); j++)
	{
	  TmpRealOffDiag[k] = this->RealOffDiagonalElements[l];
	  TmpImaginaryOffDiag[k] = this->ImaginaryOffDiagonalElements[l];
	  ++l;
	  ++k;
	}
      for (int j = this->NbrRow; j < nbrRow; j++)
	{
	  TmpRealOffDiag[k] = 0.0;
	  TmpImaginaryOffDiag[k] = 0.0;
	  ++k;
	}      
    }
  for (int i = (this->NbrRow * (this->NbrRow - 1)) / 2; i < Tot; i++)
    {
      TmpRealOffDiag[i] = 0.0;
      TmpImaginaryOffDiag[i] = 0.0;
    }
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) 
      && (this->DiagonalElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->DiagonalElements;
	delete[] this->RealOffDiagonalElements;
	delete[] this->ImaginaryOffDiagonalElements;
      }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->Flag.Initialize();
  this->DiagonalElements = TmpDiag;
  this->RealOffDiagonalElements = TmpRealOffDiag;
  this->ImaginaryOffDiagonalElements = TmpImaginaryOffDiag;
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

HermitianMatrix operator + (const HermitianMatrix& M1, const HermitianMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return HermitianMatrix();
  double* Diagonal = new double [M1.NbrRow];
  int ReducedNbr = M1.NbrRow - 1;
  int max = (ReducedNbr * M1.NbrRow) / 2;
  double* RealOffDiagonal = new double [max];
  double* ImaginaryOffDiagonal = new double [max];
  for (int i = 0; i < M1.NbrRow; i++)
    {
      Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
    }
  for (int i = 0; i < max; ++i)
    {
      RealOffDiagonal[i] = M1.RealOffDiagonalElements[i] + M2.RealOffDiagonalElements[i];      
      ImaginaryOffDiagonal[i] = M1.ImaginaryOffDiagonalElements[i] + M2.ImaginaryOffDiagonalElements[i];      
    }
  return HermitianMatrix(Diagonal, RealOffDiagonal, ImaginaryOffDiagonal, M1.NbrRow);
}

// add two matrices where the left one is a real tridiagonal symmetric matrix
//
// M1 = left matrix
// M2 = right matrix
// return value = sum of the two matrices

HermitianMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const HermitianMatrix& M2)
{
  return HermitianMatrix();
}

// add two matrices where the right one is a real tridiagonal symmetric matrix
//
// M1 = left matrix
// M2 = right matrix
// return value = sum of the two matrices

HermitianMatrix operator + (const HermitianMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2)
{
  return HermitianMatrix();
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

HermitianMatrix operator - (const HermitianMatrix& M1, const HermitianMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return HermitianMatrix();
  double* Diagonal = new double [M1.NbrRow];
  int ReducedNbr = M1.NbrRow - 1;
  int max = (ReducedNbr * M1.NbrRow) / 2;
  double* RealOffDiagonal = new double [max];
  double* ImaginaryOffDiagonal = new double [max];
  for (int i = 0; i < M1.NbrRow; ++i)
    {
      Diagonal[i] = M1.DiagonalElements[i] - M2.DiagonalElements[i];
    }
  for (int i = 0; i < max; ++i)
    {
      RealOffDiagonal[i] = M1.RealOffDiagonalElements[i] - M2.RealOffDiagonalElements[i];      
      ImaginaryOffDiagonal[i] = M1.ImaginaryOffDiagonalElements[i] - M2.ImaginaryOffDiagonalElements[i];      
    }
  return HermitianMatrix(Diagonal, RealOffDiagonal, ImaginaryOffDiagonal, M1.NbrRow);
}

// substract two matrices where the left one is a real tridiagonal symmetric matrix
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

HermitianMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, const HermitianMatrix& M2)
{
  return HermitianMatrix();
}

// substract two matrices where the right one is a real tridiagonal symmetric matrix
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

HermitianMatrix operator - (const HermitianMatrix& M1,  const RealTriDiagonalSymmetricMatrix& M2)
{
  return HermitianMatrix();
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

HermitianMatrix operator * (const HermitianMatrix& M, double x) 
{
  double* Diagonal = new double [M.NbrRow];
  int ReducedNbr = M.NbrRow - 1;
  int max = (ReducedNbr * M.NbrRow) / 2;
  double* RealOffDiagonal = new double [max];
  double* ImaginaryOffDiagonal = new double [max];
  for (int i = 0; i < M.NbrRow; ++i)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
    }
  for (int i = 0; i < max; ++i)
    {
      RealOffDiagonal[i] = M.RealOffDiagonalElements[i] * x;      
      ImaginaryOffDiagonal[i] = M.ImaginaryOffDiagonalElements[i] * x;      
    }
  return HermitianMatrix(Diagonal, RealOffDiagonal, ImaginaryOffDiagonal, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

HermitianMatrix operator * (double x, const HermitianMatrix& M) 
{
  double* Diagonal = new double [M.NbrRow];
  int ReducedNbr = M.NbrRow - 1;
  int max = (ReducedNbr * M.NbrRow) / 2;
  double* RealOffDiagonal = new double [max];
  double* ImaginaryOffDiagonal = new double [max];
  for (int i = 0; i < M.NbrRow; ++i)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
    }
  for (int i = 0; i < max; ++i)
    {
      RealOffDiagonal[i] = M.RealOffDiagonalElements[i] * x;      
      ImaginaryOffDiagonal[i] = M.ImaginaryOffDiagonalElements[i] * x;      
    }
  return HermitianMatrix(Diagonal, RealOffDiagonal, ImaginaryOffDiagonal, M.NbrRow);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

HermitianMatrix operator / (const HermitianMatrix& M, double x) 
{
  double* Diagonal = new double [M.NbrRow];
  int ReducedNbr = M.NbrRow - 1;
  int max = (ReducedNbr * M.NbrRow) / 2;
  double* RealOffDiagonal = new double [max];
  double* ImaginaryOffDiagonal = new double [max];
  double inv = 1.0 / x;
  for (int i = 0; i < M.NbrRow; ++i)
    {
      Diagonal[i] = M.DiagonalElements[i] * inv;
    }
  for (int i = 0; i < max; ++i)
    {
      RealOffDiagonal[i] = M.RealOffDiagonalElements[i] * inv;      
      ImaginaryOffDiagonal[i] = M.ImaginaryOffDiagonalElements[i] * inv;      
    }
  return HermitianMatrix(Diagonal, RealOffDiagonal, ImaginaryOffDiagonal, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

HermitianMatrix& HermitianMatrix::operator += (const HermitianMatrix& M) 
{
  if (this->NbrRow == 0)
    return *this;
  for (int i = 0; i < this->NbrRow; ++i)
    {
      this->DiagonalElements[i] += M.DiagonalElements[i];
    }
  int max = (this->NbrRow - 1) * this->NbrRow;
  for (int i = 0; i < max; ++i)
    {
      this->RealOffDiagonalElements[i] += M.RealOffDiagonalElements[i];
      this->ImaginaryOffDiagonalElements[i] += M.ImaginaryOffDiagonalElements[i]; 
    }
  return *this;
}

// add two matrices where the right one is a real tridiagonal symmetric matrix
//
// M = matrix to add to current matrix
// return value = reference on current matrix

HermitianMatrix& HermitianMatrix::operator += (const RealTriDiagonalSymmetricMatrix& M) 
{
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

HermitianMatrix& HermitianMatrix::operator -= (const HermitianMatrix& M) 
{
  if (this->NbrRow == 0)
    return *this;
  for (int i = 0; i < this->NbrRow; ++i)
    {
      this->DiagonalElements[i] -= M.DiagonalElements[i];
    }
  int max = ((this->NbrRow - 1) * this->NbrRow) / 2;
  for (int i = 0; i < max; ++i)
    {
      this->RealOffDiagonalElements[i] -= M.RealOffDiagonalElements[i];
      this->ImaginaryOffDiagonalElements[i] -= M.ImaginaryOffDiagonalElements[i]; 
    }
  return *this;
}

// substract two matrices where the right one is a real tridiagonal symmetric matrix
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

HermitianMatrix& HermitianMatrix::operator -= (const RealTriDiagonalSymmetricMatrix& M) 
{
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

HermitianMatrix& HermitianMatrix::operator *= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  for (int i = 0; i < this->NbrRow; ++i)
    {
      this->DiagonalElements[i] *= x;
    }
  int max = ((this->NbrRow - 1) * this->NbrRow) / 2;
  for (int i = 0; i < max; ++i)
    {
      this->ImaginaryOffDiagonalElements[i] *= x; 
      this->RealOffDiagonalElements[i] *= x;      
    }
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

HermitianMatrix& HermitianMatrix::operator /= (double x)
{
  if (this->NbrRow == 0)
    return *this;
  x = 1.0 / x;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] *= x;
    }
  int max = ((this->NbrRow - 1) * this->NbrRow) / 2;
  for (int i = 0; i < max; ++i)
    {
      this->RealOffDiagonalElements[i] *= x; 
      this->ImaginaryOffDiagonalElements[i] *= x;      
    }
  return *this;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex HermitianMatrix::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  Complex Z(0.0, 0.0);
  if ((V1.Dimension != this->NbrRow) || (V2.Dimension != this->NbrColumn))
    return Z;
  for (int i = 0; i < this->NbrRow ; ++i)
    {
      Complex Z2(this->DiagonalElements[i] * V2.RealComponents[i], this->DiagonalElements[i] * V2.ImaginaryComponents[i]);
      int l = (i - 1);
      for (int k = 0; k < i; ++k)
	{
	  Z2.Re += (this->RealOffDiagonalElements[l] * V2.RealComponents[k] + 
		    this->ImaginaryOffDiagonalElements[l] * V2.ImaginaryComponents[k]);
	  Z2.Im += (this->RealOffDiagonalElements[l] * V2.ImaginaryComponents[k] - 
		    this->ImaginaryOffDiagonalElements[l] * V2.RealComponents[k]);
	  l += (this->NbrColumn - 2 - k) + this->Increment;
	}
      ++l;
      for (int k = i + 1; k < this->NbrColumn; ++k)
	{
	  Z2.Re += (this->RealOffDiagonalElements[l] * V2.RealComponents[k] - 
		   this->ImaginaryOffDiagonalElements[l] * V2.ImaginaryComponents[k]);
	  Z2.Im += (this->RealOffDiagonalElements[l] * V2.ImaginaryComponents[k] + 
		   this->ImaginaryOffDiagonalElements[l] * V2.RealComponents[k]);
	  ++l;
	}      
      Z += Z2 * Complex (V1.RealComponents[i], -V1.ImaginaryComponents[i]);
    }
  return Z;
}

// conjugate matrix with an unitary matrix (Ut M U)
//
// UnitaryM = unitary matrix to use
// return value = pointer to conjugated matrix

Matrix* HermitianMatrix::Conjugate(RealMatrix& UnitaryM)
{
  if (UnitaryM.NbrRow != this->NbrColumn)
    return 0;
  double* TmpDiag = new double [UnitaryM.NbrColumn];
  int NbrOffDiag = (UnitaryM.NbrColumn * (UnitaryM.NbrColumn - 1)) / 2;
  double* TmpRealOffDiag = new double [NbrOffDiag];
  double* TmpImaginaryOffDiag = new double [NbrOffDiag];
  for (int i = 0; i < UnitaryM.NbrColumn; i++)
    {
      TmpDiag[i] = 0.0;
      for (int j = 0; j < this->NbrColumn; ++j)
	{
	  double tmp = 0.0;
	  int k = 0;
	  int l = (j - 1);
	  for (; k < j; ++k)
	    {
	      tmp += this->RealOffDiagonalElements[l] * UnitaryM.Columns[i].Components[k];
	      l += (this->NbrColumn - 2 - k) + this->Increment;
	    }
	  ++l;
	  tmp += this->DiagonalElements[j] * UnitaryM.Columns[i].Components[j];
	  ++k;
	  for (; k < this->NbrColumn; k++)
	    {
	      tmp += this->RealOffDiagonalElements[l] * UnitaryM.Columns[i].Components[k];
	      ++l;
	    }
	  TmpDiag[i] += tmp * UnitaryM.Columns[i].Components[j];
	}
    }    
  int ReducedNbrColumn = UnitaryM.NbrColumn - 1;
  for (int i = 0; i < ReducedNbrColumn; i++)
    for (int m = i + 1; m < UnitaryM.NbrColumn; m++)
      {    
	TmpRealOffDiag[i] = 0.0;
	TmpImaginaryOffDiag[i] = 0.0;
	for (int j = 0; j < this->NbrColumn; j++)
	  {
	    double tmp1 = 0.0;
	    double tmp2 = 0.0;
	    int k = 0;
	    int l = (j - 1);
	    for (; k < j; k++)
	      {
		tmp1 += this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k];
		tmp2 -= this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k];
		l += (this->NbrColumn - 3 - k) + this->Increment + 1;
	      }
	    ++l;
	    tmp1 += this->DiagonalElements[j] * UnitaryM.Columns[m].Components[j];
	    ++k;
	    for (; k < this->NbrColumn; ++k)
	      {
		tmp1 += this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k];
		tmp2 += this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k];
		++l;
	      }
	    TmpRealOffDiag[i] += tmp1 * UnitaryM.Columns[i].Components[j];
	    TmpImaginaryOffDiag[i] += tmp2 * UnitaryM.Columns[i].Components[j];
	  }
      }    
  return new HermitianMatrix(TmpDiag, TmpRealOffDiag, TmpImaginaryOffDiag, UnitaryM.NbrColumn);
}

// conjugate an hermitian matrix with an unitary matrix (Ut M U)
//
// UnitaryM = unitary matrix to use
// return value = conjugated matrix

HermitianMatrix HermitianMatrix::Conjugate(ComplexMatrix& UnitaryM)
{
  if (UnitaryM.NbrRow != this->NbrColumn)
    return HermitianMatrix();
  double* TmpDiag = new double [UnitaryM.NbrColumn];
  int NbrOffDiag = (UnitaryM.NbrColumn * (UnitaryM.NbrColumn - 1)) / 2;
  double* TmpRealOffDiag = new double [NbrOffDiag];
  double* TmpImaginaryOffDiag = new double [NbrOffDiag];
  for (int i = 0; i < UnitaryM.NbrColumn; i++)
    {
      TmpDiag[i] = 0.0;
      for (int j = 0; j < this->NbrColumn; ++j)
	{
	  double tmp1 = 0.0;
	  double tmp2 = 0.0;
	  int k = 0;
	  int l = (j - 1);
	  for (; k < j; ++k)
	    {
	      tmp1 += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[i].RealComponents[k] + 
		      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[i].ImaginaryComponents[k]);
	      tmp2 += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[i].ImaginaryComponents[k] - 
		      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[i].RealComponents[k]);
	      ++k;
	      l += (this->NbrColumn - 2 - k) + this->Increment;
	    }
	  ++l;
	  tmp1 += this->DiagonalElements[j] * UnitaryM.Columns[i].RealComponents[k];
	  tmp2 += this->DiagonalElements[j] * UnitaryM.Columns[i].ImaginaryComponents[k];
	  ++k;
	  for (; k < this->NbrColumn; ++k)
	    {
	      tmp1 += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[i].RealComponents[k] - 
		      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[i].ImaginaryComponents[k]);
	      tmp2 += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[i].ImaginaryComponents[k] + 
		      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[i].RealComponents[k]);
	      ++l;
	      ++k;
	    }
	  TmpDiag[i] += tmp1 * UnitaryM.Columns[i].RealComponents[j] + tmp2 * UnitaryM.Columns[i].ImaginaryComponents[j];
	}
    }    
  int ReducedNbrColumn = UnitaryM.NbrColumn - 1;
  for (int i = 0; i < ReducedNbrColumn; i++)
    for (int m = i + 1; m < UnitaryM.NbrColumn; m++)
      {    
	TmpRealOffDiag[i] = 0.0;
	TmpImaginaryOffDiag[i] = 0.0;
	for (int j = 0; j < this->NbrColumn; j++)
	  {
	    double tmp1 = 0.0;
	    double tmp2 = 0.0;
	    int k = 0;
	    int l = (j - 1);
	    for (; k < j; k++)
	      {
		tmp1 += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].RealComponents[k] + 
			 this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].ImaginaryComponents[k]);
		tmp2 += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].ImaginaryComponents[k] - 
			 this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].RealComponents[k]);
		l += (this->NbrColumn - 2 - k) + this->Increment;
	      }
	    ++l;
	    tmp1 += this->DiagonalElements[j] * UnitaryM.Columns[m].RealComponents[k];
	    tmp2 += this->DiagonalElements[j] * UnitaryM.Columns[m].ImaginaryComponents[k];
	    ++k;
	    for (; k < this->NbrColumn; ++k)
	      {
		tmp1 += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].RealComponents[k] - 
			 this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].ImaginaryComponents[k]);
		tmp2 += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].ImaginaryComponents[k] + 
			 this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].RealComponents[k]);
		++l;
	      }
	    TmpRealOffDiag[i] += tmp1 * UnitaryM.Columns[i].RealComponents[j] + tmp2 * UnitaryM.Columns[i].ImaginaryComponents[j];
	    TmpImaginaryOffDiag[i] += tmp2 * UnitaryM.Columns[i].RealComponents[j] - tmp1 * UnitaryM.Columns[i].ImaginaryComponents[j] ;
	  }
      }    
  return HermitianMatrix(TmpDiag, TmpRealOffDiag, TmpImaginaryOffDiag, UnitaryM.NbrColumn);
}

// conjugate an hermitian matrix with an hermitian transposed unitary matrix (U M Ut)
//
// UnitaryM = unitary matrix to use
// return value = conjugated matrix

HermitianMatrix HermitianMatrix::InvConjugate(ComplexMatrix& UnitaryM)
{
  if (UnitaryM.NbrColumn != this->NbrColumn)
    return HermitianMatrix();
  double* TmpDiag = new double [UnitaryM.NbrRow];
  int NbrOffDiag = (UnitaryM.NbrRow * (UnitaryM.NbrRow - 1)) / 2;
  double* TmpRealOffDiag = new double [NbrOffDiag];
  double* TmpImaginaryOffDiag = new double [NbrOffDiag];
  for (int i = 0; i < UnitaryM.NbrRow; ++i)
    {
      TmpDiag[i] = 0.0;
      for (int j = 0; j < this->NbrColumn; ++j)
	{
	  double tmp1 = 0.0;
	  double tmp2 = 0.0;
	  int k = 0;
	  int l = (j - 1);
	  for (; k < j; k++)
	    {
	      tmp1 += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[k].RealComponents[i] - 
		      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[k].ImaginaryComponents[i]);
	      tmp2 -= (this->RealOffDiagonalElements[l] * UnitaryM.Columns[k].ImaginaryComponents[i] + 
		      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[k].RealComponents[i]);
	      l += (this->NbrColumn - 2 - k) + this->Increment;
	    }
	  ++l;
	  tmp1 += this->DiagonalElements[j] * UnitaryM.Columns[k].RealComponents[i];
	  tmp2 -= this->DiagonalElements[j] * UnitaryM.Columns[k].ImaginaryComponents[i];
	  ++k;
	  for (; k < this->NbrColumn; ++k)
	    {
	      tmp1 += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[k].RealComponents[i] + 
		      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[k].ImaginaryComponents[i]);
	      tmp2 += (this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[k].RealComponents[i] - 
		       this->RealOffDiagonalElements[l] * UnitaryM.Columns[k].ImaginaryComponents[i] );
	      ++l;
	    }
	  TmpDiag[i] += tmp1 * UnitaryM.Columns[i].RealComponents[j] - tmp2 * UnitaryM.Columns[i].ImaginaryComponents[j];
	}
    }    
  int ReducedNbrRow = UnitaryM.NbrRow - 1;
  for (int i = 0; i < ReducedNbrRow; i++)
    for (int m = i + 1; m < UnitaryM.NbrRow; m++)
      {    
	TmpRealOffDiag[i] = 0.0;
	TmpImaginaryOffDiag[i] = 0.0;
	for (int j = 0; j < this->NbrColumn; j++)
	  {
	    double tmp1 = 0.0;
	    double tmp2 = 0.0;
	    int k = 0;
	    int l = (j - 1);
	    for (; k < j; k++)
	      {
		tmp1 += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[k].RealComponents[m] - 
			 this->ImaginaryOffDiagonalElements[l + 1] * UnitaryM.Columns[k].ImaginaryComponents[m]);
		tmp2 -= (this->RealOffDiagonalElements[l] * UnitaryM.Columns[k].ImaginaryComponents[m] + 
			 this->ImaginaryOffDiagonalElements[l + 1] * UnitaryM.Columns[k].RealComponents[m]);
		l += (this->NbrColumn - 2 - k) + this->Increment;
	      }
	    ++l;
	    tmp1 += this->DiagonalElements[j] * UnitaryM.Columns[k].RealComponents[m];
	    tmp2 += this->DiagonalElements[j] * UnitaryM.Columns[k].ImaginaryComponents[m];
	    ++k;
	    for (; k < this->NbrColumn; ++k)
	      {
		tmp1 += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[k].RealComponents[m] + 
			 this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[k].ImaginaryComponents[m]);
		tmp2 += (this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[k].RealComponents[m] - 
			 this->RealOffDiagonalElements[l] * UnitaryM.Columns[k].ImaginaryComponents[m]);
		++l;
	      }
	    TmpRealOffDiag[i] += tmp1 * UnitaryM.Columns[j].RealComponents[i] - tmp2 * UnitaryM.Columns[j].ImaginaryComponents[i];
	    TmpImaginaryOffDiag[i] += tmp2 * UnitaryM.Columns[j].RealComponents[i] + tmp1 * UnitaryM.Columns[j].ImaginaryComponents[i] ;
	  }
      }    
  return HermitianMatrix(TmpDiag, TmpRealOffDiag, TmpImaginaryOffDiag, UnitaryM.NbrColumn);
}

// evaluate matrix trace
//
// return value = matrix trace 

double HermitianMatrix::Tr () 
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

double HermitianMatrix::Det () 
{
  return 1.0;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const HermitianMatrix& P)
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      int pos = 2 * (i - 1);
      for (int j = 0; j < i; ++j)
	{
	  Str << P.RealOffDiagonalElements[pos];
	  if (P.ImaginaryOffDiagonalElements[pos] > 0.0)
	    Str << -P.ImaginaryOffDiagonalElements[pos] << "i    ";
	  else
	    if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
	      Str << "+" << -P.ImaginaryOffDiagonalElements[pos] << "i    ";
	    else
	      Str << "    ";
	  pos += (P.NbrRow - j - 2) + P.Increment;
	}
      Str << P.DiagonalElements[i] << "    ";
      ++pos;
      for (int j = i + 1; j < P.NbrRow; ++j)
	{
	  Str << P.RealOffDiagonalElements[pos];
	  if (P.ImaginaryOffDiagonalElements[pos] < 0.0)
	    Str << P.ImaginaryOffDiagonalElements[pos] << "i    ";
	  else
	    if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
	      Str << "+" << P.ImaginaryOffDiagonalElements[pos] << "i    ";
	    else
	      Str << "    ";
	  pos += 2;
	}
      Str << endl;
    }
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// P = matrix to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, const HermitianMatrix& P)
{
  Str << "{";
  for (int i = 0; i < P.NbrRow; ++i)
    {
      Str << "{";
      int pos = i - 1;
      for (int j = 0; j < i; ++j)
	{
	  if ((P.RealOffDiagonalElements[pos] != 0) || (P.ImaginaryOffDiagonalElements[pos + 1] == 0))
	    Str << P.RealOffDiagonalElements[pos];
	  if (P.ImaginaryOffDiagonalElements[pos] > 0.0)
	    Str << -P.ImaginaryOffDiagonalElements[pos] << "I";
	  else
	    if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
	      Str << "+" << -P.ImaginaryOffDiagonalElements[pos] << "I";
	  Str << ",";
	  pos += (P.NbrRow - j - 2) + P.Increment;
	}
      Str << P.DiagonalElements[i];
      if (i != (P.NbrRow - 1))
	{
	  Str << ",";	  
	  ++pos;
	  for (int j = i + 1; j < (P.NbrRow - 1); ++j)
	    {
	      if ((P.RealOffDiagonalElements[pos] != 0) || (P.ImaginaryOffDiagonalElements[pos] == 0))
		Str << P.RealOffDiagonalElements[pos];
	      if (P.ImaginaryOffDiagonalElements[pos] < 0.0)
		Str << P.ImaginaryOffDiagonalElements[pos] << "I";
	      else
		if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
		  Str << "+" << P.ImaginaryOffDiagonalElements[pos] << "I";
	      Str << ",";
	      pos += 2;
	    }
	  Str << P.RealOffDiagonalElements[pos];
	  if (P.ImaginaryOffDiagonalElements[pos] < 0.0)
	    Str << P.ImaginaryOffDiagonalElements[pos] << "I";
	  else
	    if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
	      Str << "+" << P.ImaginaryOffDiagonalElements[pos] << "I";
	  Str << "},";
	}
      else
	Str << "}";
    }
  Str << "}";
  return Str;
}

// Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step
//
// dimension = maximum iteration number
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// V1 = reference on complex vector used as first vector (will contain last produced vector at the end)
// return value = reference on complex tridiagonal hermitian matrix

RealTriDiagonalSymmetricMatrix& HermitianMatrix::Lanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, ComplexVector& V1)
{  
  int Index = 0;
  V1 /= V1.Norm();
  ComplexVector V2(this->NbrRow);
  ComplexVector V3(this->NbrRow);
  V2.Multiply(*this, V1);
  M.DiagonalElements[Index] = (V1 * V2).Re;
  V2.AddLinearCombination(-M.DiagonalElements[Index], V1);
  V2 /= V2.Norm();
  dimension -= 2;
  for (int i = 0; i < dimension; i++)
    {
      V3.Multiply(*this, V2);
      M.UpperDiagonalElements[Index] = (V1 * V3).Re;
      M.DiagonalElements[Index + 1] = (V2 * V3).Re;
      V3.AddLinearCombination(-M.DiagonalElements[Index + 1], V2);
      V3.AddLinearCombination(-M.UpperDiagonalElements[Index], V1);
      V3 /= V3.Norm();
      Index++;
      ComplexVector TmpV = V1;
      V1 = V2;
      V2 = V3;
      V3 = TmpV;
    }  
  V3.Multiply(*this, V2);
  M.UpperDiagonalElements[Index] = (V1 * V3).Re;
  M.DiagonalElements[Index + 1] = (V2 * V3).Re;
  ComplexVector TmpV = V1;
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

RealTriDiagonalSymmetricMatrix& HermitianMatrix::Lanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, ComplexMatrix& Q)
{  
  int Index = 0;
  if ((Q.NbrColumn != dimension) || (Q.NbrRow != this->NbrRow))
    Q.Resize(this->NbrRow, dimension);
  if (M.NbrRow != dimension)
    M.Resize(dimension, dimension);
  Q.Columns[0] /= Q.Columns[0].Norm(); 
  Q.Columns[1].Multiply(*this, Q.Columns[0]);
  M.DiagonalElements[Index] = (Q.Columns[0] * Q.Columns[1]).Re;
  Q.Columns[1].AddLinearCombination(-M.DiagonalElements[Index], Q.Columns[0]);
  Q.Columns[1] /= Q.Columns[1].Norm(); 
  for (int i = 2; i < dimension; i++)
    {
      Q.Columns[i].Multiply(*this, Q.Columns[i - 1]);
      M.UpperDiagonalElements[Index] = (Q.Columns[i - 2] * Q.Columns[i]).Re;
      M.DiagonalElements[Index + 1] = (Q.Columns[i - 1] * Q.Columns[i]).Re;
      Q.Columns[i].AddLinearCombination(-M.DiagonalElements[Index + 1], Q.Columns[i - 1]);
      Q.Columns[i].AddLinearCombination(-M.UpperDiagonalElements[Index], Q.Columns[i - 2]);
      Q.Columns[i] /= Q.Columns[i].Norm();
      Index++;
    }  
  M.UpperDiagonalElements[Index] = this->MatrixElement(Q.Columns[dimension - 2], Q.Columns[dimension - 1]).Re;
  M.DiagonalElements[Index + 1] = this->MatrixElement(Q.Columns[dimension - 1], Q.Columns[dimension - 1]).Re;
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

RealTriDiagonalSymmetricMatrix& HermitianMatrix::OrthoLanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, ComplexMatrix& Q,
							       double err)
{  
/*  int Index = 0;
  if ((Q.NbrColumn != dimension) || (Q.NbrRow != this->NbrRow))
    Q.Resize(this->NbrRow, dimension);
  if (M.NbrRow != dimension)
    M.Resize(dimension, dimension);
  Q.Columns[0] /= Q.Columns[0].Norm(); 
  Q.Columns[1].Multiply(*this, Q.Columns[0]);
  M.DiagonalElements[Index] = (Q.Columns[0] * Q.Columns[1]).Re;
  Q.Columns[1].AddLinearCombination(-M.DiagonalElements[Index], Q.Columns[0]);
  Q.Columns[1] /= Q.Columns[1].Norm(); 
  for (int i = 2; i < dimension; i++)
    {
      Q.Columns[i].Multiply(*this, Q.Columns[i - 1]);
      M.UpperDiagonalElements[Index] = (Q.Columns[i - 2] * Q.Columns[i]).Re;
      M.DiagonalElements[Index + 1] = (Q.Columns[i - 1] * Q.Columns[i]).Re;
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
	      Q.Columns[i].AddLinearCombination(- (Q.Columns[j] * Q.Columns[i]).Re, Q.Columns[j]);
	    }
	  VectorNorm = Q.Columns[i].Norm();
	}
      Q.Columns[i] /= VectorNorm;
      Index++;
    }  
  M.UpperDiagonalElements[Index] = this->MatrixElement(Q.Columns[dimension - 2], Q.Columns[dimension - 1]).Re;
  M.DiagonalElements[Index + 1] = this->MatrixElement(Q.Columns[dimension - 1], Q.Columns[dimension - 1]).Re;*/
  return M;
}

// Tridiagonalize a hermitian matrix using Householder algorithm and evaluate transformation matrix
//
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// err = absolute error on matrix element
// Q = matrix where transformation matrix has to be stored
// return value = reference on real tridiagonal symmetric matrix

RealTriDiagonalSymmetricMatrix& HermitianMatrix::Householder (RealTriDiagonalSymmetricMatrix& M, double err, ComplexMatrix& Q)
{
  if (M.NbrRow != this->NbrRow)
    M.Resize(this->NbrRow, this->NbrColumn);
  if ((Q.NbrRow != M.NbrRow) || (Q.NbrColumn != M.NbrColumn))
    Q.Resize(this->NbrRow, this->NbrColumn);
   for (int i = 0; i < Q.NbrRow; i++)
     {
       for (int j = 0; j < i; j++)
	 {
	   Q.Columns[j].RealComponents[i] = 0.0;       
	   Q.Columns[j].ImaginaryComponents[i] = 0.0;       
	 }
       Q.Columns[i].RealComponents[i] = 1.0;
       Q.Columns[i].ImaginaryComponents[i] = 0.0;       
       for (int j = i + 1; j < Q.NbrColumn; j++)
	 {
	   Q.Columns[j].RealComponents[i] = 0.0;       
	   Q.Columns[j].ImaginaryComponents[i] = 0.0;       
	 }
     }
  int ReducedNbrRow = this->NbrRow -1;
  int ReducedNbrRow2 ;
  double* TmpVRe = new double [ReducedNbrRow];
  double* TmpVIm = new double [ReducedNbrRow];
  double* TmpCoefRe = new double [this->NbrRow];
  double* TmpCoefIm = new double [this->NbrRow];
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
  for (int i = 1; i < ReducedNbrRow; ++i)
    {
      ReducedNbrRow2 = this->NbrRow - i;
      TmpNorm = 0.0;
      TmpPos = Pos;
      // construct vector for Householder transformation
      for (int j = 0; j < ReducedNbrRow2; ++j)
	{
	  TmpNorm += this->RealOffDiagonalElements[TmpPos] * this->RealOffDiagonalElements[TmpPos];
	  TmpNorm += this->ImaginaryOffDiagonalElements[TmpPos] * this->ImaginaryOffDiagonalElements[TmpPos];
	  TmpVRe[j] = 0.0;
	  TmpVIm[j] = 0.0;
	  ++TmpPos;
	}
      TmpPos = Pos;
      Coef = sqrt(TmpNorm);
      M.UpperDiagonalElement(i- 1) = Coef;
/*      TmpNorm = sqrt(TmpNorm - Coef * this->OffDiagonalElements[Pos]);
      if  (TmpNorm > err)
	{
	  TmpNorm = 1.0 / TmpNorm;
 	  this->OffDiagonalElements[Pos] -= Coef;
 	  for (int j = 0; j < ReducedNbrRow2; ++j)
 	    {
 	      this->OffDiagonalElements[TmpPos] *= TmpNorm;
	      ++TmpPos;
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
	      TmpCoefRe[TmpPos5] = 0.0;
	      TmpCoefIm[TmpPos5] = 0.0;
	      for (int k = i; k < this->NbrRow; ++k)
		{
	 	  TmpCoefRe[TmpPos5] += Q.Columns[k].RealComponents[j] * this->RealOffDiagonalElements[TmpPos]
		    + Q.Columns[k].ImaginaryComponents[j] * this->ImaginaryOffDiagonalElements[TmpPos];
	 	  TmpCoefIm[TmpPos5] += Q.Columns[k].RealComponents[j] * this->ImaginaryOffDiagonalElements[TmpPos]
		    -  Q.Columns[k].ImaginaryComponents[j] * this->RealOffDiagonalElements[TmpPos];
		  ++TmpPos;
		}
	      ++TmpPos5;
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
	  
	  // conjugate Hermitian matrix
	  TmpPos2 = i * (this->NbrRow + this->Increment) - (i * (i + 1)) / 2;
	  TmpPos = Pos + 1;
	  TmpPos6 = Pos;
	  for (int j = 1; j < i; j++)
	    {
	      TmpPos5 = TmpPos6;
	      Coef = TmpCoef[j];
	      for (int k = i; k < this->NbrRow; ++k)
		{
		  Q.Columns[k].Components[j] -=  this->OffDiagonalElements[TmpPos5] * Coef;
		  ++TmpPos5;
		}	      
	    }
	  for (int j = i; j < this->NbrRow; ++j)
	    {	  
	      TmpPos4 = j - i;
	      TmpPos5 = TmpPos6;
	      Coef = TmpCoef[j];
	      for (int k = i; k < j; ++k)
		{
		  Q.Columns[k].Components[j] -=  this->OffDiagonalElements[TmpPos5++] * Coef;
		}	      
	      this->DiagonalElements[j] -=  this->OffDiagonalElements[Pos] * 2.0 * TmpV[TmpPos4];
	      Q.Columns[j].Components[j] -=  this->OffDiagonalElements[TmpPos5++] * Coef;
	      TmpPos3 = TmpPos;
	      for (int k = j + 1; k < this->NbrRow; ++k)
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
      Pos += this->Increment;*/
    }
//  M.UpperDiagonalElement(ReducedNbrRow - 1) = this->OffDiagonalElements[Pos - this->Increment];  
  M.DiagonalElement(ReducedNbrRow) = this->DiagonalElements[ReducedNbrRow];
  delete[] TmpVRe;
  delete[] TmpVIm;
  return M;
}

