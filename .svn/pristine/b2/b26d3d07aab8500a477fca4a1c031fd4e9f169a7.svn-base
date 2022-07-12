////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of complex lower triangular matrix                 //
//                                                                            //
//                        last modification : 20/08/2004                      //
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


#include "Matrix/ComplexLowerTriangularMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "GeneralTools/ListIterator.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "MathTools/Complex.h"

#include <stdlib.h>
#include <fstream>


using std::endl;


// default constructor
//

ComplexLowerTriangularMatrix::ComplexLowerTriangularMatrix() 
{
  this->RealDiagonalElements = 0;
  this->ImaginaryDiagonalElements = 0;
  this->RealOffDiagonalElements = 0;
  this->ImaginaryOffDiagonalElements = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = 0;
  this->MatrixType = Matrix::ComplexElements | Matrix::Triangular | Matrix::Lower;
  this->Dummy = 0.0;
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

ComplexLowerTriangularMatrix::ComplexLowerTriangularMatrix(int dimension, bool zero) 
{
  this->DiagonalFlag.Initialize();
  this->OffDiagonalFlag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::ComplexElements | Matrix::Triangular | Matrix::Lower;
  this->RealDiagonalElements = new double [this->NbrRow];
  this->ImaginaryDiagonalElements = new double [this->NbrRow];
  this->RealOffDiagonalElements = new double [(this->NbrRow * (this->NbrRow - 1)) / 2];
  this->ImaginaryOffDiagonalElements = new double [(this->NbrRow * (this->NbrRow - 1)) / 2];
  if (zero == true)
    {
      int pos = 0;
      for (int i = 0; i < this->NbrRow; i++)
	{
	  this->RealDiagonalElements[i] = 0.0;
	  this->ImaginaryDiagonalElements[i] = 0.0;
	  for (int j = i + 1; j < this->NbrRow; j++)
	    {
	      this->RealOffDiagonalElements[pos] = 0.0;
	      this->ImaginaryOffDiagonalElements[pos] = 0.0;
	      pos++;
	    }
	}
    }
  this->Dummy = 0.0;
}

// constructor from matrix elements (without duplicating datas)
//
// realDiagonal = pointer to real part of the diagonal elements
// imaginaryDiagonal = pointer to imaginary part of the diagonal elements
// realOffDiagonal = pointer to real part of the off-diagonal elements
// imaginaryOffDiagonal = pointer to imaginary part of the off-diagonal elements
// dimension = matrix dimension

ComplexLowerTriangularMatrix::ComplexLowerTriangularMatrix(double* realDiagonal, double* imaginaryDiagonal, 
							   double* realOffDiagonal, double* imaginaryOffDiagonal, int dimension) 
{
  this->RealDiagonalElements = realDiagonal;
  this->ImaginaryDiagonalElements = imaginaryDiagonal;
  this->RealOffDiagonalElements = realOffDiagonal;
  this->ImaginaryOffDiagonalElements = imaginaryOffDiagonal;
  this->DiagonalFlag.Initialize();
  this->OffDiagonalFlag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::ComplexElements | Matrix::Triangular | Matrix::Lower;
  this->Dummy = 0.0;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

ComplexLowerTriangularMatrix::ComplexLowerTriangularMatrix(const ComplexLowerTriangularMatrix& M) 
{
  this->RealDiagonalElements = M.RealDiagonalElements;
  this->ImaginaryDiagonalElements = M.ImaginaryDiagonalElements;
  this->DiagonalFlag = M.DiagonalFlag;
  this->RealOffDiagonalElements = M.RealOffDiagonalElements;
  this->ImaginaryOffDiagonalElements = M.ImaginaryOffDiagonalElements;
  this->OffDiagonalFlag = M.OffDiagonalFlag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::ComplexElements | Matrix::Triangular | Matrix::Lower;
  this->Dummy = 0.0;
}

// destructor
//

ComplexLowerTriangularMatrix::~ComplexLowerTriangularMatrix() 
{
  if ((this->RealDiagonalElements != 0) && (this->ImaginaryDiagonalElements != 0) && 
      (this->DiagonalFlag.Used() == true) && (this->DiagonalFlag.Shared() == false))
    {
      delete[] this->RealDiagonalElements;
      delete[] this->ImaginaryDiagonalElements;
    }
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) && 
      (this->OffDiagonalFlag.Used() == true) && (this->OffDiagonalFlag.Shared() == false))
    {
      delete[] this->RealOffDiagonalElements;
      delete[] this->ImaginaryOffDiagonalElements;
    }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

ComplexLowerTriangularMatrix& ComplexLowerTriangularMatrix::operator = (const ComplexLowerTriangularMatrix& M) 
{
  if ((this->RealDiagonalElements != 0) && (this->ImaginaryDiagonalElements != 0) && 
      (this->DiagonalFlag.Used() == true) && (this->DiagonalFlag.Shared() == false))
    {
      delete[] this->RealDiagonalElements;
      delete[] this->ImaginaryDiagonalElements;
    }
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) && 
      (this->OffDiagonalFlag.Used() == true) && (this->OffDiagonalFlag.Shared() == false))
    {
      delete[] this->RealOffDiagonalElements;
      delete[] this->ImaginaryOffDiagonalElements;
    }
  this->RealDiagonalElements = M.RealDiagonalElements;
  this->ImaginaryDiagonalElements = M.ImaginaryDiagonalElements;
  this->DiagonalFlag = M.DiagonalFlag;
  this->RealOffDiagonalElements = M.RealOffDiagonalElements;
  this->ImaginaryOffDiagonalElements = M.ImaginaryOffDiagonalElements;
  this->OffDiagonalFlag = M.OffDiagonalFlag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::ComplexElements | Matrix::Triangular | Matrix::Lower;
  this->Dummy = 0.0;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* ComplexLowerTriangularMatrix::Clone ()
{
  return ((Matrix*) new ComplexLowerTriangularMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexLowerTriangularMatrix::SetMatrixElement(int i, int j, double x)
{
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexLowerTriangularMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      this->RealDiagonalElements[i] = x.Re;
      this->ImaginaryDiagonalElements[i] = x.Im;
    }
  else
    if (i > j)
      {
	i += (j * (j - 1)) / 2;
	this->RealOffDiagonalElements[i] = x.Re;
	this->ImaginaryOffDiagonalElements[i] = x.Im;
      }
  return;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void ComplexLowerTriangularMatrix::AddToMatrixElement(int i, int j, double x)
{
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void ComplexLowerTriangularMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      this->RealDiagonalElements[i] += x.Re;
      this->ImaginaryDiagonalElements[i] += x.Im;
    }
  else
    if (i > j)
      {
	j += (i * (i - 1)) / 2;
	this->RealOffDiagonalElements[i] += x.Re;
	this->ImaginaryOffDiagonalElements[i] += x.Im;
      }
  return;
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexLowerTriangularMatrix::Resize (int nbrRow, int nbrColumn)
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
  double* TmpRealDiag = new double [nbrRow];
  double* TmpImaginaryDiag = new double [nbrRow];
  int Tot = (nbrRow * (nbrRow - 1)) / 2;
  double* TmpRealOffDiag = new double [Tot];
  double* TmpImaginaryOffDiag = new double [Tot];
  for (int i = 0; i < this->NbrRow; i++)
    {
      TmpRealDiag [i] = this->RealDiagonalElements[i];
      TmpImaginaryDiag [i] = this->ImaginaryDiagonalElements[i];
    }
  for (int i = this->NbrRow; i < nbrRow; i++)
    {
      TmpRealDiag [i]  = 0.0;
      TmpImaginaryDiag [i]  = 0.0;
    }
  int k = 0;
  int l = 0;
  for (int i = 0; i < (this->NbrRow - 1); i++)
    {
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  TmpRealOffDiag[k] = this->RealOffDiagonalElements[l];
	  TmpImaginaryOffDiag[k] = this->ImaginaryOffDiagonalElements[l];
	  ++k;
	  ++l;
	}
      l += this->Increment;
      for (int j = this->NbrRow; j < nbrRow; j++)
	{
	  TmpRealOffDiag[k] = 0.0;
	  TmpImaginaryOffDiag[k] = 0.0;
	  ++k;
	}      
    }
  for (int i = this->NbrRow * (this->NbrRow - 1); i < Tot; i++)
    {
      TmpRealOffDiag[i] = 0.0;
      TmpImaginaryOffDiag[i] = 0.0;
    }
  if ((this->RealDiagonalElements != 0) && (this->ImaginaryDiagonalElements != 0) && 
      (this->DiagonalFlag.Used() == true) && (this->DiagonalFlag.Shared() == false))
    {
      delete[] this->RealDiagonalElements;
      delete[] this->ImaginaryDiagonalElements;
    }
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) && 
      (this->OffDiagonalFlag.Used() == true) && (this->OffDiagonalFlag.Shared() == false))
    {
      delete[] this->RealOffDiagonalElements;
      delete[] this->ImaginaryOffDiagonalElements;
    }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = this->TrueNbrRow - this->NbrRow;
  this->RealDiagonalElements = TmpRealDiag;
  this->ImaginaryDiagonalElements = TmpImaginaryDiag;
  this->RealOffDiagonalElements = TmpRealOffDiag;
  this->ImaginaryOffDiagonalElements = TmpImaginaryOffDiag;
  this->DiagonalFlag.Initialize();
  this->OffDiagonalFlag.Initialize();
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexLowerTriangularMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  if (nbrRow <= this->TrueNbrRow)
    {
      if (this->NbrRow < nbrRow)
	{
	  int Tot = (nbrRow * (nbrRow - 1));
	  for (int i = this->NbrRow; i < nbrRow; i++)
	    {
	      this->RealDiagonalElements[i] = 0.0;
	      this->ImaginaryDiagonalElements[i] = 0.0;
	    }
	  int k = (this->NbrRow - 1);
	  for (int i = 0; i < (this->NbrRow - 1); i++)
	    {
	      for (int j = this->NbrRow; j < nbrRow; j++)
		{
		  this->RealOffDiagonalElements[k] = 0.0;
		  this->ImaginaryOffDiagonalElements[k] = 0.0;
		  ++k;
		}
	      k += (this->NbrRow - 2 - i);
	    }
	  for (int i = this->NbrRow * (this->NbrRow - 1); i < Tot; i++)
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
  double* TmpRealDiag = new double [nbrRow];
  double* TmpImaginaryDiag = new double [nbrRow];
  int Tot = (nbrRow * (nbrRow - 1)) / 2;
  double* TmpRealOffDiag = new double [Tot];
  double* TmpImaginaryOffDiag = new double [Tot];
  for (int i = 0; i < this->NbrRow; i++)
    {
      TmpRealDiag [i] = this->RealDiagonalElements[i];
      TmpImaginaryDiag [i] = this->ImaginaryDiagonalElements[i];
    }
  for (int i = this->NbrRow; i < nbrRow; i++)
    {
      TmpRealDiag [i]  = 0.0;
      TmpImaginaryDiag [i]  = 0.0;
    }
  int k = 0;
  int l = 0;
  for (int i = 0; i < (this->NbrRow - 1); i++)
    {
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  TmpRealOffDiag[k] = this->RealOffDiagonalElements[l];
	  TmpImaginaryOffDiag[k] = this->ImaginaryOffDiagonalElements[l];
	  ++k;
	  ++l;
	}
      l += this->Increment;
      for (int j = this->NbrRow; j < nbrRow; j++)
	{
	  TmpRealOffDiag[k] = 0.0;
	  TmpImaginaryOffDiag[k] = 0.0;
	  ++k;
	}      
    }
  for (int i = this->NbrRow * (this->NbrRow - 1); i < Tot; i++)
    {
      TmpRealOffDiag[i] = 0.0;
      TmpImaginaryOffDiag[i] = 0.0;
    }
  if ((this->RealDiagonalElements != 0) && (this->ImaginaryDiagonalElements != 0) && 
      (this->DiagonalFlag.Used() == true) && (this->DiagonalFlag.Shared() == false))
    {
      delete[] this->RealDiagonalElements;
      delete[] this->ImaginaryDiagonalElements;
    }
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) && 
      (this->OffDiagonalFlag.Used() == true) && (this->OffDiagonalFlag.Shared() == false))
    {
      delete[] this->RealOffDiagonalElements;
      delete[] this->ImaginaryOffDiagonalElements;
    }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = this->TrueNbrRow - this->NbrRow;
  this->RealDiagonalElements = TmpRealDiag;
  this->ImaginaryDiagonalElements = TmpImaginaryDiag;
  this->RealOffDiagonalElements = TmpRealOffDiag;
  this->ImaginaryOffDiagonalElements = TmpImaginaryOffDiag;
  this->DiagonalFlag.Initialize();
  this->OffDiagonalFlag.Initialize();
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

ComplexLowerTriangularMatrix operator + (const ComplexLowerTriangularMatrix& M1, const ComplexLowerTriangularMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexLowerTriangularMatrix();
  double* RealDiagonal = new double [M1.NbrRow];
  double* ImaginaryDiagonal = new double [M1.NbrRow];
  int ReducedNbr = M1.NbrRow - 1;
  double* RealOffDiagonal = new double [M1.NbrRow * ReducedNbr];
  double* ImaginaryOffDiagonal = new double [M1.NbrRow * ReducedNbr];
  for (int i = 0; i < M1.NbrRow; i++)
    {
      RealDiagonal[i] = M1.RealDiagonalElements[i] + M2.RealDiagonalElements[i];
      ImaginaryDiagonal[i] = M1.ImaginaryDiagonalElements[i] + M2.ImaginaryDiagonalElements[i];
    }
  int k = 0;
  int l1 = 0;
  int l2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = 0; j < i; j++)
	{
	  RealOffDiagonal[k] = M1.RealOffDiagonalElements[l1] + M2.RealOffDiagonalElements[l2];      
	  ImaginaryOffDiagonal[k] = M1.ImaginaryOffDiagonalElements[l1] + M2.ImaginaryOffDiagonalElements[l2];      
	  ++k;
	  ++l2;
	  ++l1;
	}
      l1 += M2.Increment;
      l2 += M2.Increment;
    }
  return ComplexLowerTriangularMatrix(RealDiagonal, ImaginaryDiagonal, RealOffDiagonal, ImaginaryOffDiagonal, M1.NbrRow);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexLowerTriangularMatrix operator - (const ComplexLowerTriangularMatrix& M1, const ComplexLowerTriangularMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexLowerTriangularMatrix();
  double* RealDiagonal = new double [M1.NbrRow];
  double* ImaginaryDiagonal = new double [M1.NbrRow];
  int ReducedNbr = M1.NbrRow - 1;
  double* RealOffDiagonal = new double [M1.NbrRow * ReducedNbr];
  double* ImaginaryOffDiagonal = new double [M1.NbrRow * ReducedNbr];
  for (int i = 0; i < M1.NbrRow; i++)
    {
      RealDiagonal[i] = M1.RealDiagonalElements[i] - M2.RealDiagonalElements[i];
      ImaginaryDiagonal[i] = M1.ImaginaryDiagonalElements[i] - M2.ImaginaryDiagonalElements[i];
    }
  int k = 0;
  int l1 = 0;
  int l2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = 0; j < i; j++)
	{
	  RealOffDiagonal[k] = M1.RealOffDiagonalElements[l1] - M2.RealOffDiagonalElements[l2];      
	  ImaginaryOffDiagonal[k] = M1.ImaginaryOffDiagonalElements[l1] - M2.ImaginaryOffDiagonalElements[l2];      
	  ++k;
	  ++l2;
	  ++l1;
	}
      l1 += M2.Increment;
      l2 += M2.Increment;
    }
  return ComplexLowerTriangularMatrix(RealDiagonal, ImaginaryDiagonal, RealOffDiagonal, ImaginaryOffDiagonal, M1.NbrRow);
}

// multiply a real matrix with a real lower triangular matrix
//
// m1 = real matrix
// m2 = real lower triangular matrix
// return value = product result

ComplexMatrix operator * (ComplexMatrix& m1, const ComplexLowerTriangularMatrix& m2)
{
  ComplexMatrix TmpM(m1.GetNbrRow(), m2.NbrRow);
  Complex Tmp = 0.0;
  int Pos = 0;
  for (int i = 0; i < m1.GetNbrRow(); ++i)
    {
      for (int j = 0; j < m2.NbrRow; ++j)
	{
	  Tmp.Re = (m1[j].Re(i) * m2.RealDiagonalElements[j]
		    - m1[j].Im(i) * m2.ImaginaryDiagonalElements[j]);
	  Tmp.Im = (m1[j].Im(i) * m2.RealDiagonalElements[j]
		    + m1[j].Re(i) * m2.ImaginaryDiagonalElements[j]);
	  Pos = ((j - 1) * j) / 2;
	  for (int k = 0; k < j; ++k)
	    {
	      Tmp.Re += (m1[k].Re(i) * m2.RealDiagonalElements[Pos]
			 - m1[k].Im(i) * m2.ImaginaryDiagonalElements[Pos]);
	      Tmp.Im += (m1[k].Im(i) * m2.RealDiagonalElements[Pos]
			 + m1[k].Re(i) * m2.ImaginaryDiagonalElements[Pos]);
	      ++Pos;
	    }
	  TmpM[j].Re(i) = Tmp.Re;	
	  TmpM[j].Im(i) = Tmp.Im;	
	}
    }
  return TmpM;
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexLowerTriangularMatrix operator * (const ComplexLowerTriangularMatrix& M, double x) 
{
  double* RealDiagonal = new double [M.NbrRow];
  double* ImaginaryDiagonal = new double [M.NbrRow];
  int ReducedNbr = M.NbrRow - 1;
  double* RealOffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  double* ImaginaryOffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  for (int i = 0; i < M.NbrRow; i++)
    {
      RealDiagonal[i] = M.RealDiagonalElements[i] * x;
      ImaginaryDiagonal[i] = M.ImaginaryDiagonalElements[i] * x;
    }
  int k = 0;
  int k2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  RealOffDiagonal[k] = M.RealOffDiagonalElements[k2] * x;
	  ImaginaryOffDiagonal[k] = M.ImaginaryOffDiagonalElements[k2] * x;
	  ++k;
	  ++k2;
	}
      k2 += M.Increment;
    }
  return ComplexLowerTriangularMatrix(RealDiagonal, ImaginaryDiagonal, RealOffDiagonal, ImaginaryOffDiagonal, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexLowerTriangularMatrix operator * (double x, const ComplexLowerTriangularMatrix& M) 
{
  return (M * x);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

ComplexLowerTriangularMatrix operator / (const ComplexLowerTriangularMatrix& M, double x) 
{
  x = 1.0 / x;
  return (M * x);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

ComplexLowerTriangularMatrix& ComplexLowerTriangularMatrix::operator += (const ComplexLowerTriangularMatrix& M) 
{
  if (this->NbrRow == 0)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  for (int i = 0; i < M.NbrRow; i++)
    {
      this->RealDiagonalElements[i] += M.RealDiagonalElements[i];
      this->ImaginaryDiagonalElements[i] += M.ImaginaryDiagonalElements[i];
    }
  int k = 0;
  int k2 = 0;  
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  this->RealOffDiagonalElements[k] += M.RealOffDiagonalElements[k2];
	  this->ImaginaryOffDiagonalElements[k] += M.ImaginaryOffDiagonalElements[k2];
	  ++k;
	  ++k2;
	}
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

ComplexLowerTriangularMatrix& ComplexLowerTriangularMatrix::operator -= (const ComplexLowerTriangularMatrix& M) 
{
  if (this->NbrRow == 0)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  for (int i = 0; i < M.NbrRow; i++)
    {
      this->RealDiagonalElements[i] -= M.RealDiagonalElements[i];
      this->ImaginaryDiagonalElements[i] -= M.ImaginaryDiagonalElements[i];
    }
  int k = 0;
  int k2 = 0;  
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  this->RealOffDiagonalElements[k] -= M.RealOffDiagonalElements[k2];
	  this->ImaginaryOffDiagonalElements[k] -= M.ImaginaryOffDiagonalElements[k2];
	  ++k;
	  ++k2;
	}
      k += this->Increment;
      k2 += M.Increment;
    }
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexLowerTriangularMatrix& ComplexLowerTriangularMatrix::operator *= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  int ReducedNbr = this->NbrRow - 1;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->RealDiagonalElements[i] *= x;
      this->ImaginaryDiagonalElements[i] *= x;
    }
  int k = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  this->RealOffDiagonalElements[k] *= x;
	  this->ImaginaryOffDiagonalElements[k] *= x;
	  ++k;
	}
      k += this->Increment;
    }
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexLowerTriangularMatrix& ComplexLowerTriangularMatrix::operator /= (double x)
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

Complex ComplexLowerTriangularMatrix::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  Complex x = 0.0;
  if ((V1.Dimension != this->NbrRow) || (V2.Dimension != this->NbrColumn))
    return x;
  for (int i = 0; i < this->NbrRow ; i++)
    {
      Complex x2;
      x2.Re = this->RealDiagonalElements[i] * V2.RealComponents[i] - this->ImaginaryDiagonalElements[i] * V2.ImaginaryComponents[i];
      x2.Im = this->ImaginaryDiagonalElements[i] * V2.RealComponents[i] + this->RealDiagonalElements[i] * V2.ImaginaryComponents[i];
      int l = (i - 1) ;
      for (int k = 0; k < i; k++)
	{
	  x2.Re += this->RealDiagonalElements[l] * V2.RealComponents[k] - this->ImaginaryDiagonalElements[l] * V2.ImaginaryComponents[k];
	  x2.Im += this->ImaginaryDiagonalElements[l] * V2.RealComponents[k] + this->RealDiagonalElements[l] * V2.ImaginaryComponents[k];
	  l += (this->NbrColumn - 2 - k) + this->Increment;
	}
      x.Re += (x2.Re * V1.RealComponents[i] + x2.Im * V1.ImaginaryComponents[i]);
      x.Im += (x2.Im * V1.RealComponents[i] - x2.Re * V1.ImaginaryComponents[i]);
    }
  return x;
}

// evaluate matrix trace
//
// return value = matrix trace 

double ComplexLowerTriangularMatrix::Tr () 
{
  return 0.0;
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double ComplexLowerTriangularMatrix::Det () 
{
  return 0.0;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const ComplexLowerTriangularMatrix& P)
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      int pos = ((i - 1) * i) / 2;
      for (int j = 0; j < i; j ++)
	{
	  Str << P.RealOffDiagonalElements[pos] << "    ";
	  if (P.ImaginaryOffDiagonalElements[pos] < 0.0)
	    Str << P.ImaginaryOffDiagonalElements[pos] << "i    ";
	  else
	    if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
	      Str << "+" << P.ImaginaryOffDiagonalElements[pos] << "i    ";
	    else
	      Str << "    ";
	}
      Str << P.RealDiagonalElements[i] << "    ";
      if (P.ImaginaryDiagonalElements[i] < 0.0)
	Str << P.ImaginaryDiagonalElements[i] << "i    ";
      else
	if (P.ImaginaryDiagonalElements[i] != 0.0)
	  Str << "+" << P.ImaginaryDiagonalElements[i] << "i    ";
	else
	  Str << "    ";
      for (int j = i + 1; j < P.NbrRow; j++)
	{
	  Str << 0.0 << "    ";
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

MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexLowerTriangularMatrix& P)
{
  Str << "{";
  for (int i = 0; i < P.NbrRow; i++)
    {
      Str << "{";
      int pos = (i - 1);
      for (int j = 0; j < i; j ++)
	{
	      if ((P.RealOffDiagonalElements[pos] != 0) || (P.ImaginaryOffDiagonalElements[pos] == 0))
		Str << P.RealOffDiagonalElements[pos];
	      if (P.ImaginaryOffDiagonalElements[pos] < 0.0)
		Str << P.ImaginaryOffDiagonalElements[pos] << "I";
	      else
		if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
		  Str << "+" << P.ImaginaryOffDiagonalElements[pos] << "I";
	      Str << ",";
	  pos += (P.NbrRow - j - 2) + P.Increment;
	}
      if ((P.RealDiagonalElements[i] != 0) || (P.ImaginaryDiagonalElements[i] == 0))
	Str << P.RealDiagonalElements[i];
      if (P.ImaginaryDiagonalElements[i] < 0.0)
	Str << P.ImaginaryDiagonalElements[i] << "I";
      else
	if (P.ImaginaryDiagonalElements[i] != 0.0)
	  Str << "+" << P.ImaginaryDiagonalElements[i] << "I";
      if (i != (P.NbrRow - 1))
	{
	  Str << ",";	  
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
