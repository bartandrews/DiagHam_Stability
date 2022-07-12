////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of complex skew symmetric matrix                  //
//                                                                            //
//                        last modification : 19/08/2004                      //
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


#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "GeneralTools/ListIterator.h"

#include <stdlib.h>


using std::cout;
using std::endl;


// default constructor
//

ComplexSkewSymmetricMatrix::ComplexSkewSymmetricMatrix() 
{
  this->Dummy = 0.0;
  this->RealOffDiagonalElements = 0;
  this->ImaginaryOffDiagonalElements =  0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = 0;
  this->MatrixType = Matrix::ComplexElements | Matrix::Antisymmetric;
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

ComplexSkewSymmetricMatrix::ComplexSkewSymmetricMatrix(int dimension, bool zero) 
{
  this->Dummy = 0.0;
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::ComplexElements | Matrix::Antisymmetric;
  if (this->NbrRow > 1)
    {
      this->RealOffDiagonalElements = new double [(this->NbrRow * (this->NbrRow - 1)) / 2];
      this->ImaginaryOffDiagonalElements = new double [(this->NbrRow * (this->NbrRow - 1)) / 2];
    }
  else
    {    
      this->RealOffDiagonalElements = new double [1];
      this->ImaginaryOffDiagonalElements = new double [1];
    }
  if (zero == true)
    {
      int pos = 0;
      for (int i = 0; i < this->NbrRow; i++)
	{
	  for (int j = i + 1; j < this->NbrRow; j++)
	    {
	      this->RealOffDiagonalElements[pos] = 0.0;
	      this->ImaginaryOffDiagonalElements[pos] = 0.0;
	      pos++;
	    }
	}
    }
}

// constructor from matrix elements (without duplicating datas)
//
// realUpperDiagonal = pointer to real part of the upper-diagonal elements
// imaginaryUpperDiagonal = pointer to imaginary part of the upper-diagonal elements
// dimension = matrix dimension

ComplexSkewSymmetricMatrix::ComplexSkewSymmetricMatrix(double* realUpperDiagonal, double* imaginaryUpperDiagonal, int dimension) 
{
  this->Dummy = 0.0;
  this->RealOffDiagonalElements = realUpperDiagonal;  
  this->ImaginaryOffDiagonalElements = imaginaryUpperDiagonal;  
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::ComplexElements | Matrix::Antisymmetric;
}

// copy constructor
//
// M = matrix to copy
// duplicateFlag = true if datas have to be duplicated

ComplexSkewSymmetricMatrix::ComplexSkewSymmetricMatrix(const ComplexSkewSymmetricMatrix& M, bool duplicateFlag) 
{
  this->Dummy = 0.0;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::ComplexElements | Matrix::Antisymmetric;
  if (duplicateFlag == false)
    {
      this->RealOffDiagonalElements = M.RealOffDiagonalElements;  
      this->ImaginaryOffDiagonalElements = M.ImaginaryOffDiagonalElements;  
      this->Flag = M.Flag;
    }
  else
    {
      this->Flag.Initialize();
      if (this->TrueNbrRow > 1)
	{
	  this->RealOffDiagonalElements = new double [(this->TrueNbrRow * (this->TrueNbrRow - 1)) / 2];
	  this->ImaginaryOffDiagonalElements = new double [(this->TrueNbrRow * (this->TrueNbrRow - 1)) / 2];
	}
      else
	{    
	  this->RealOffDiagonalElements = new double [1];
	  this->ImaginaryOffDiagonalElements = new double [1];
	}
      int pos = 0;
      for (int i = 0; i < this->TrueNbrRow; ++i)
	{
	  for (int j = i + 1; j < this->TrueNbrRow; ++j)
	    {
	      this->RealOffDiagonalElements[pos] = M.RealOffDiagonalElements[pos];
	      this->ImaginaryOffDiagonalElements[pos] = M.ImaginaryOffDiagonalElements[pos];
	      ++pos;
	    }
	}
    }
}

// destructor
//

ComplexSkewSymmetricMatrix::~ComplexSkewSymmetricMatrix() 
{
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) 
      && (this->Flag.Used() == true) && (this->Flag.Shared() == false))
      {
	delete[] this->RealOffDiagonalElements;
	delete[] this->ImaginaryOffDiagonalElements;
      }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

ComplexSkewSymmetricMatrix& ComplexSkewSymmetricMatrix::operator = (const ComplexSkewSymmetricMatrix& M) 
{
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) 
      && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealOffDiagonalElements;
	delete[] this->ImaginaryOffDiagonalElements;
      }
  this->RealOffDiagonalElements = M.RealOffDiagonalElements;
  this->ImaginaryOffDiagonalElements = M.ImaginaryOffDiagonalElements;
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->Increment = this->TrueNbrRow - this->NbrRow;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* ComplexSkewSymmetricMatrix::Clone ()
{
  return ((Matrix*) new ComplexSkewSymmetricMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexSkewSymmetricMatrix::SetMatrixElement(int i, int j, double x)
{
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element
void ComplexSkewSymmetricMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn) || (i == j))
    return;
  if (i > j)
    {
      i -= (j * (j + 1)) / 2 - j * (this->NbrRow + this->Increment - 1) + 1;
      this->RealOffDiagonalElements[i] = -x.Re;
      this->ImaginaryOffDiagonalElements[i] = -x.Im;
    }
  else
    {
      j -= (i * (i + 1)) / 2 - i * (this->NbrRow + this->Increment - 1) + 1;
      this->RealOffDiagonalElements[j] = x.Re;
      this->ImaginaryOffDiagonalElements[j] = x.Im;
    }
}

// get a matrix element (real part if complex)
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

void ComplexSkewSymmetricMatrix::GetMatrixElement(int i, int j, double& x) const
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn) || (i == j))
    {
      x = 0.0;
      return;
    }
  if (i > j)
    {
      i -= (j * (j + 1)) / 2 - j * (this->NbrRow + this->Increment - 1) + 1;
      x = -this->RealOffDiagonalElements[i];
    }
  else
    {
      j -= (i * (i + 1)) / 2 - i * (this->NbrRow + this->Increment - 1) + 1;
      x = this->RealOffDiagonalElements[j];
    }
}
    
// get a matrix element
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

void ComplexSkewSymmetricMatrix::GetMatrixElement(int i, int j, Complex& x) const
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn) || (i == j))
    {
      x.Re = 0.0;
      x.Im = 0.0;
      return;
    }
  if (i > j)
    {
      i -= (j * (j + 1)) / 2 - j * (this->NbrRow + this->Increment - 1) + 1;
      x.Re = -this->RealOffDiagonalElements[i];
      x.Im = -this->ImaginaryOffDiagonalElements[i];
    }
  else
    {
      j -= (i * (i + 1)) / 2 - i * (this->NbrRow + this->Increment - 1) + 1;
      x.Re = this->RealOffDiagonalElements[j];
      x.Im = this->ImaginaryOffDiagonalElements[j];
    }
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void ComplexSkewSymmetricMatrix::AddToMatrixElement(int i, int j, double x)
{
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void ComplexSkewSymmetricMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn) || (i == j))
    return;
  if (i > j)
    {
      i -= (j * (j + 1)) / 2 - j * (this->NbrRow + this->Increment - 1) + 1;
      this->RealOffDiagonalElements[i] -= x.Re;
      this->ImaginaryOffDiagonalElements[i] -= x.Im;
    }
  else
    {
      j -= (i * (i + 1)) / 2 - i * (this->NbrRow + this->Increment - 1) + 1;
      this->RealOffDiagonalElements[j] += x.Re;
      this->ImaginaryOffDiagonalElements[j] += x.Im;
    }
}

// get reference of a given matrix element supposing i < j
//
// i = line position
// j = column position
// return value = reference om matrix elememt

double& ComplexSkewSymmetricMatrix::operator () (int i, int j)
{
  return this->Dummy;
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexSkewSymmetricMatrix::Resize (int nbrRow, int nbrColumn)
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
  int Tot = (nbrRow * (nbrRow - 1)) / 2;
  double* TmpRealOffDiag = new double [Tot];
  double* TmpImaginaryOffDiag = new double [Tot];
  int k = 0;
  int l = 0;
  for (int i = 0; i < (this->NbrRow - 1); i++)
    {
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  TmpRealOffDiag[k] = this->RealOffDiagonalElements[l];
	  TmpImaginaryOffDiag[k] = this->ImaginaryOffDiagonalElements[l];
	  ++l;
	  ++k;
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
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) 
      && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealOffDiagonalElements;
	delete[] this->ImaginaryOffDiagonalElements;
      }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = this->TrueNbrRow - this->NbrRow;
  this->RealOffDiagonalElements = TmpRealOffDiag;
  this->ImaginaryOffDiagonalElements = TmpImaginaryOffDiag;
  this->Flag.Initialize();
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexSkewSymmetricMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  if (nbrRow <= this->TrueNbrRow)
    {
      if (this->NbrRow < nbrRow)
	{
	  int Tot = (nbrRow * (nbrRow - 1));
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
  int Tot = (nbrRow * (nbrRow - 1)) / 2;
  double* TmpRealOffDiag = new double [Tot];
  double* TmpImaginaryOffDiag = new double [Tot];
  int k = 0;
  int l = 0;
  for (int i = 0; i < (this->NbrRow - 1); i++)
    {
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  TmpRealOffDiag[k] = this->RealOffDiagonalElements[l];
	  TmpImaginaryOffDiag[k] = this->ImaginaryOffDiagonalElements[l];
	  ++l;
	  ++k;
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
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) 
      && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealOffDiagonalElements;
	delete[] this->ImaginaryOffDiagonalElements;
      }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = this->TrueNbrRow - this->NbrRow;
  this->RealOffDiagonalElements = TmpRealOffDiag;
  this->ImaginaryOffDiagonalElements = TmpImaginaryOffDiag;
  this->Flag.Initialize();
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

ComplexSkewSymmetricMatrix operator + (const ComplexSkewSymmetricMatrix& M1, const ComplexSkewSymmetricMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexSkewSymmetricMatrix();
  int ReducedNbr = M1.NbrRow - 1;
  double* RealOffDiagonal = new double [M1.NbrRow * ReducedNbr];
  double* ImaginaryOffDiagonal = new double [M1.NbrRow * ReducedNbr];
  int k = 0;
  int l1 = 0;
  int l2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = 0; j < i; j++)
	{
	  RealOffDiagonal[k] = M1.RealOffDiagonalElements[l1] + M2.RealOffDiagonalElements[l2];      
	  ImaginaryOffDiagonal[k++] = M1.ImaginaryOffDiagonalElements[l1++] + M2.ImaginaryOffDiagonalElements[l2++];      
	}
      l1 += M2.Increment;
      l2 += M2.Increment;
    }
  return ComplexSkewSymmetricMatrix(RealOffDiagonal, ImaginaryOffDiagonal, M1.NbrRow);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexSkewSymmetricMatrix operator - (const ComplexSkewSymmetricMatrix& M1, const ComplexSkewSymmetricMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexSkewSymmetricMatrix();
  int ReducedNbr = M1.NbrRow - 1;
  double* RealOffDiagonal = new double [M1.NbrRow * ReducedNbr];
  double* ImaginaryOffDiagonal = new double [M1.NbrRow * ReducedNbr];
  int k = 0;
  int l1 = 0;
  int l2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = 0; j < i; j++)
	{
	  RealOffDiagonal[k] = M1.RealOffDiagonalElements[l1] - M2.RealOffDiagonalElements[l2];      
	  ImaginaryOffDiagonal[k++] = M1.ImaginaryOffDiagonalElements[l1++] - M2.ImaginaryOffDiagonalElements[l2++];      
	}
      l1 += M2.Increment;
      l2 += M2.Increment;
    }
  return ComplexSkewSymmetricMatrix(RealOffDiagonal, ImaginaryOffDiagonal, M1.NbrRow);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexSkewSymmetricMatrix operator * (const ComplexSkewSymmetricMatrix& M, double x) 
{
  int ReducedNbr = M.NbrRow - 1;
  double* RealOffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  double* ImaginaryOffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  int k = 0;
  int k2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  RealOffDiagonal[k] = M.RealOffDiagonalElements[k2] * x;
	  ImaginaryOffDiagonal[k++] = M.ImaginaryOffDiagonalElements[k2++] * x;
	}
      k2 += M.Increment;
    }
  return ComplexSkewSymmetricMatrix(RealOffDiagonal, ImaginaryOffDiagonal, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexSkewSymmetricMatrix operator * (double x, const ComplexSkewSymmetricMatrix& M) 
{
  int ReducedNbr = M.NbrRow - 1;
  double* RealOffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  double* ImaginaryOffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  int k = 0;
  int k2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  RealOffDiagonal[k] = M.RealOffDiagonalElements[k2] * x;
	  ImaginaryOffDiagonal[k++] = M.ImaginaryOffDiagonalElements[k2++] * x;
	}
      k2 += M.Increment;
    }
  return ComplexSkewSymmetricMatrix(RealOffDiagonal, ImaginaryOffDiagonal, M.NbrRow);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

ComplexSkewSymmetricMatrix operator / (const ComplexSkewSymmetricMatrix& M, double x) 
{
  x = 1.0 / x;
  int ReducedNbr = M.NbrRow - 1;
  double* RealOffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  double* ImaginaryOffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  int k = 0;
  int k2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  RealOffDiagonal[k] = M.RealOffDiagonalElements[k2] * x;
	  ImaginaryOffDiagonal[k++] = M.ImaginaryOffDiagonalElements[k2++] * x;
	}
      k2 += M.Increment;
    }
  return ComplexSkewSymmetricMatrix(RealOffDiagonal, ImaginaryOffDiagonal, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

ComplexSkewSymmetricMatrix& ComplexSkewSymmetricMatrix::operator += (const ComplexSkewSymmetricMatrix& M) 
{
  if (this->NbrRow < M.NbrRow)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  int k = 0;
  int k2 = 0;  
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  this->RealOffDiagonalElements[k] += M.RealOffDiagonalElements[k2];
	  this->ImaginaryOffDiagonalElements[k++] += M.ImaginaryOffDiagonalElements[k2++];
	}
      k += this->Increment;
      k2 += M.Increment;
    }
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

ComplexSkewSymmetricMatrix& ComplexSkewSymmetricMatrix::operator -= (const ComplexSkewSymmetricMatrix& M) 
{
  if (this->NbrRow < M.NbrRow)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  int k = 0;
  int k2 = 0;  
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  this->RealOffDiagonalElements[k] -= M.RealOffDiagonalElements[k2];
	  this->ImaginaryOffDiagonalElements[k++] -= M.ImaginaryOffDiagonalElements[k2++];
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

ComplexSkewSymmetricMatrix& ComplexSkewSymmetricMatrix::operator *= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  int ReducedNbr = this->NbrRow - 1;
  int k = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  this->RealOffDiagonalElements[k] *= x;
	  this->ImaginaryOffDiagonalElements[k++] *= x;
	}
      k += this->Increment;
    }
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexSkewSymmetricMatrix& ComplexSkewSymmetricMatrix::operator /= (double x)
{
  if (this->NbrRow == 0)
    return *this;
  x = 1.0 / x;
  int ReducedNbr = this->NbrRow - 1;
  int k = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  this->RealOffDiagonalElements[k] *= x;
	  this->ImaginaryOffDiagonalElements[k++] *= x;
	}
      k += this->Increment;
    }
  return *this;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ComplexSkewSymmetricMatrix::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  Complex x;
  if ((V1.Dimension != this->NbrRow) || (V2.Dimension != this->NbrColumn))
    return x;
  for (int i = 0; i < this->NbrRow ; i++)
    {
      Complex x2;
      int l = (i - 1);
      for (int k = 0; k < i; k++)
	{
	  x2.Re -= (this->RealOffDiagonalElements[l] * V2.Components[k].Re - 
		    this->ImaginaryOffDiagonalElements[l] * V2.Components[k].Im);
	  x2.Im -= (this->ImaginaryOffDiagonalElements[l] * V2.Components[k].Re + 
		    this->RealOffDiagonalElements[l] * V2.Components[k].Im);
	  l += (this->NbrColumn - 2 - k) + this->Increment;
	}
      l++;
      for (int k = i + 1; k < this->NbrColumn; k++)
	{
	  x2.Re += (this->RealOffDiagonalElements[l] * V2.Components[k].Re -
		    this->ImaginaryOffDiagonalElements[l] * V2.Components[k].Im);
	  x2.Im += (this->ImaginaryOffDiagonalElements[l] * V2.Components[k].Re + 
		    this->RealOffDiagonalElements[l] * V2.Components[k].Im);
	  ++l;
	}      
      x.Re += (x2.Re * V1.Components[i].Re + x2.Im * V1.Components[i].Im);
      x.Re += (x2.Im * V1.Components[i].Re - x2.Re * V1.Components[i].Im);
    }
  return x;
}

// conjugate an hermitian matrix with an unitary matrix (Uh M U)
//
// UnitaryM = unitary matrix to use
// return value = conjugated matrix

Matrix* ComplexSkewSymmetricMatrix::Conjugate(ComplexMatrix& UnitaryM)
{
  if (UnitaryM.NbrRow != this->NbrColumn)
    return 0;
  int NbrOffDiag = (UnitaryM.NbrColumn * (UnitaryM.NbrColumn - 1)) / 2;
  double* TmpRealOffDiag = new double [NbrOffDiag];
  double* TmpImaginaryOffDiag = new double [NbrOffDiag];
  int i2 = 0;
  int ReducedNbrColumn = UnitaryM.NbrColumn - 1;
  int Inc = this->NbrColumn - 3 + this->Increment;
  Complex tmp1;
  int k;
  int l;
  for (int i = 0; i < ReducedNbrColumn; i++)
    {
      for (int m = i + 1; m < UnitaryM.NbrColumn; m++)
	{    
	  TmpRealOffDiag[i2] = 0.0;
	  TmpImaginaryOffDiag[i2] = 0.0;
	  for (int j = 0; j < this->NbrColumn; j++)
	    {
	      tmp1 = 0.0;
	      k = 0;
	      l = (j - 1);
	      for (; k < j; k++)
		{
		  tmp1.Re -= (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Re - 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Im);
		  tmp1.Im -= (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Im + 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Re);
		  l += Inc - k;
		}
	      l++;
	      k++;
	      for (; k < this->NbrColumn; k++)
		{
		  tmp1.Re += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Re - 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Im);
		  tmp1.Im += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Im + 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Re);
		  ++l;
		}
	      TmpRealOffDiag[i2] += tmp1.Re * UnitaryM.Columns[i].Components[j].Re + tmp1.Im * UnitaryM.Columns[i].Components[j].Im;
	      TmpImaginaryOffDiag[i2] += tmp1.Im * UnitaryM.Columns[i].Components[j].Re - tmp1.Re * UnitaryM.Columns[i].Components[j].Im ;
	    }
	  ++i2;
	}    
    }
  return new ComplexSkewSymmetricMatrix(TmpRealOffDiag, TmpImaginaryOffDiag, UnitaryM.NbrColumn);
}

// conjugate a block of the matrix with an unitary matrix (Ut M U)
//
// UnitaryM = unitary matrix to use
// sourcePosition = index of the row where the block to conjugate starts
// destinationPosition = index of the row where the conjugated block has to be stored
// matrix = matrix where result has to be stored

void ComplexSkewSymmetricMatrix::Conjugate(ComplexMatrix& UnitaryM, int sourcePosition, int destinationPosition,
					   ComplexSkewSymmetricMatrix& matrix)
{
  if (((UnitaryM.NbrRow + sourcePosition) > this->NbrColumn) || 
      ((UnitaryM.NbrColumn + destinationPosition) > matrix.NbrColumn))
    return;
  int i2 = (destinationPosition - (destinationPosition * (destinationPosition + 1)) / 2 +
	    destinationPosition * (matrix.NbrRow + matrix.Increment - 1));
  int ReducedNbrColumn = UnitaryM.NbrColumn - 1;
  for (int i = 0; i < ReducedNbrColumn; i++)
    {
      for (int m = i + 1; m < UnitaryM.NbrColumn; m++)
	{    
	  matrix.RealOffDiagonalElements[i2] = 0.0;
	  matrix.ImaginaryOffDiagonalElements[i2] = 0.0;
	  for (int j = 0; j < UnitaryM.NbrRow; j++)
	    {
	      Complex tmp1;
	      int k = 0;
	      int l = ((j + sourcePosition) - 1 - (sourcePosition * (sourcePosition + 1)) / 2 +
		       sourcePosition * (this->NbrRow + this->Increment - 1));
	      for (; k < j; k++)
		{
		  tmp1.Re -= (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Re + 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Im);
		  tmp1.Im -= (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Im - 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Re);
		  l += (this->NbrColumn - 2 - k - sourcePosition) + this->Increment;
		}
	      l++;
	      k++;
	      for (; k < UnitaryM.NbrRow; k++)
		{
		  tmp1.Re += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Re - 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Im);
		  tmp1.Im += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Im + 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].Components[k].Re);
		  ++l;
		}
	      matrix.RealOffDiagonalElements[i2] += tmp1.Re * UnitaryM.Columns[i].Components[j].Re + tmp1.Im * UnitaryM.Columns[i].Components[j].Im;
	      matrix.RealOffDiagonalElements[i2] += tmp1.Im * UnitaryM.Columns[i].Components[j].Re - tmp1.Re * UnitaryM.Columns[i].Components[j].Im;
	    }
	  i2++;
	}
      i2 += matrix.NbrColumn - destinationPosition - UnitaryM.NbrColumn + matrix.Increment;
    }
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

void ComplexSkewSymmetricMatrix::Conjugate(ComplexMatrix& UnitaryMl, ComplexMatrix& UnitaryMr, int sourceRowIndex, 
					   int sourceColumnIndex, int destinationRowIndex,
					   int destinationColumnIndex, ComplexSkewSymmetricMatrix& matrix)
{
  if (((UnitaryMr.NbrRow + sourceColumnIndex) > this->NbrColumn) || 
      ((UnitaryMl.NbrRow + sourceRowIndex) > this->NbrRow) || 
      ((UnitaryMr.NbrColumn + destinationColumnIndex) > matrix.NbrColumn) || 
      ((UnitaryMl.NbrColumn + destinationRowIndex) > matrix.NbrRow))
    return;
  for (int i = 0; i < UnitaryMl.NbrColumn; i++)
    {
      int i2 = (destinationColumnIndex - 1 - ((i + destinationRowIndex) * ((i + destinationRowIndex) + 1)) / 2 +
		(i + destinationRowIndex) * (matrix.NbrRow + matrix.Increment - 1));
      for (int m = 0; m < UnitaryMr.NbrColumn; m++)
	{    
	  matrix.RealOffDiagonalElements[i2] = 0.0;
	  matrix.ImaginaryOffDiagonalElements[i2] = 0.0;
	  for (int j = 0; j < UnitaryMl.NbrRow; j++)
	    {
	      Complex tmp1;
	      int l = (sourceColumnIndex - 1 - 
		       ((j + sourceRowIndex) * ((sourceRowIndex + j) + 1)) / 2 +
		       (sourceRowIndex + j) * (this->NbrRow + this->Increment - 1));
	      for (int k = 0; k < UnitaryMr.NbrRow; k++)
		{
		  tmp1.Re += (this->RealOffDiagonalElements[l] * UnitaryMr.Columns[m].Components[k].Re - 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryMr.Columns[m].Components[k].Im);
		  tmp1.Im += (this->RealOffDiagonalElements[l] * UnitaryMr.Columns[m].Components[k].Im + 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryMr.Columns[m].Components[k].Re);
		  ++l;
		}
	      matrix.RealOffDiagonalElements[i2] += tmp1.Re * UnitaryMl.Columns[i].Components[j].Re + tmp1.Im * UnitaryMl.Columns[i].Components[j].Im;
	      matrix.ImaginaryOffDiagonalElements[i2] += tmp1.Im * UnitaryMl.Columns[i].Components[j].Re - tmp1.Re * UnitaryMl.Columns[i].Components[j].Im ;
	    }
	  i2++;
	}
    }
  return;
}

// swap the i-th row/column with the j-th row/column (thus preserving the skew symmetric form)
//
// i = index of the first the row/column
// j = index of the second the row/column
// return value = reference on the current matrix

ComplexSkewSymmetricMatrix& ComplexSkewSymmetricMatrix::SwapRowColumn (int i, int j)
{
  if (i == j)  
    return *this;
  if (i > j)
    {
      int Tmp = j;
      j = i;
      i = Tmp;
    }
  double Tmp;
  int Pos1 = i - 1;
  int Pos2 = j - 1;
  int k = 0;
  int TmpInc = this->Increment + this->NbrColumn - 2;
  for (; k < i; ++k)
    {
      Tmp = this->RealOffDiagonalElements[Pos1]; 
      this->RealOffDiagonalElements[Pos1] = this->RealOffDiagonalElements[Pos2]; 
      this->RealOffDiagonalElements[Pos2] = Tmp;
      Tmp = this->ImaginaryOffDiagonalElements[Pos1]; 
      this->ImaginaryOffDiagonalElements[Pos1] = this->ImaginaryOffDiagonalElements[Pos2]; 
      this->ImaginaryOffDiagonalElements[Pos2] = Tmp;
      Pos1 += TmpInc - k; 
      Pos2 += TmpInc - k; 
    }
  this->RealOffDiagonalElements[Pos2] *= -1.0;
  this->ImaginaryOffDiagonalElements[Pos2] *= -1.0;
  ++Pos1;
  Pos2 +=  TmpInc - k;
  ++k;
  for (; k < j; ++k)
    {
      Tmp = this->RealOffDiagonalElements[Pos1]; 
      this->RealOffDiagonalElements[Pos1] = -this->RealOffDiagonalElements[Pos2]; 
      this->RealOffDiagonalElements[Pos2] = -Tmp;
      Tmp = this->ImaginaryOffDiagonalElements[Pos1]; 
      this->ImaginaryOffDiagonalElements[Pos1] = -this->ImaginaryOffDiagonalElements[Pos2]; 
      this->ImaginaryOffDiagonalElements[Pos2] = -Tmp;
      ++Pos1; 
      Pos2 += TmpInc - k;
    }
  ++k;
  ++Pos1;
  ++Pos2;
  for (; k < this->NbrColumn; ++k)
    {
      Tmp = this->RealOffDiagonalElements[Pos1]; 
      this->RealOffDiagonalElements[Pos1] = this->RealOffDiagonalElements[Pos2]; 
      this->RealOffDiagonalElements[Pos2] = Tmp;
      Tmp = this->ImaginaryOffDiagonalElements[Pos1]; 
      this->ImaginaryOffDiagonalElements[Pos1] = this->ImaginaryOffDiagonalElements[Pos2]; 
      this->ImaginaryOffDiagonalElements[Pos2] = Tmp;
      ++Pos1; 
      ++Pos2;
    }
  return *this;
}


// evaluate matrix trace
//
// return value = matrix trace 

double ComplexSkewSymmetricMatrix::Tr () 
{
  return 0.0;
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double ComplexSkewSymmetricMatrix::Det () 
{
  return 1.0;
}

// evaluate matrix pfaffian
//
// return value = matrix pfaffian 

Complex ComplexSkewSymmetricMatrix::Pfaffian()
{
  Complex Pfaffian;
  if (this->NbrColumn == 2)
    {
      Pfaffian.Re = this->RealOffDiagonalElements[0];
      Pfaffian.Im = this->ImaginaryOffDiagonalElements[0];
      return Pfaffian;
    }
  if (this->NbrColumn == 4)
    {
      Pfaffian.Re = (this->RealOffDiagonalElements[0] * this->RealOffDiagonalElements[5 + (2 * this->Increment)]
		 - this->ImaginaryOffDiagonalElements[0] * this->ImaginaryOffDiagonalElements[5 + (2 * this->Increment)]);
      Pfaffian.Im = (this->ImaginaryOffDiagonalElements[0] * this->RealOffDiagonalElements[5 + (2 * this->Increment)]
		+ this->RealOffDiagonalElements[0] * this->ImaginaryOffDiagonalElements[5 + (2 * this->Increment)]);
      Pfaffian.Re -= (this->RealOffDiagonalElements[1] * this->RealOffDiagonalElements[4 + this->Increment]
		 - this->ImaginaryOffDiagonalElements[1] * this->ImaginaryOffDiagonalElements[4 + this->Increment]);
      Pfaffian.Im -= (this->ImaginaryOffDiagonalElements[1] * this->RealOffDiagonalElements[4 + this->Increment]
		 + this->RealOffDiagonalElements[1] * this->ImaginaryOffDiagonalElements[4 + this->Increment]);
      Pfaffian.Re += (this->RealOffDiagonalElements[2] * this->RealOffDiagonalElements[3 + this->Increment]
		 - this->ImaginaryOffDiagonalElements[2] * this->ImaginaryOffDiagonalElements[3 + this->Increment]);
      Pfaffian.Im += (this->ImaginaryOffDiagonalElements[2] * this->RealOffDiagonalElements[3 + this->Increment]
		 + this->RealOffDiagonalElements[2] * this->ImaginaryOffDiagonalElements[3 + this->Increment]);
//      cout << Pfaffian << " " << (Pfaffian * Pfaffian) << endl;
//      return Pfaffian;
    }

  Pfaffian = 1.0;

  ComplexSkewSymmetricMatrix TmpMatrix(*this, true);
  int Pos;
  Complex Pivot;
  double PivotNorm;
  double TmpNorm;
  int PivotColumn;
  int PivotRow;
  int PivotPos;
  int ReducedK;
  int TmpInc;
  for (int k = this->NbrColumn; k >= 4; k -= 2)
    {      
      // find pivot (greater real value)
      Pos = 0;
      PivotNorm = ((TmpMatrix.RealOffDiagonalElements[0] * TmpMatrix.RealOffDiagonalElements[0]) + 
		   (TmpMatrix.ImaginaryOffDiagonalElements[0] * TmpMatrix.ImaginaryOffDiagonalElements[0]));
      PivotPos = 0;
      PivotColumn = 1;
      PivotRow = 0;
      for (int i = 0; i < k; ++i)
	{
	  for (int j = i + 1; j < k; ++j) 
	    {
	      TmpNorm = ((TmpMatrix.RealOffDiagonalElements[Pos] * TmpMatrix.RealOffDiagonalElements[Pos]) + 
			 (TmpMatrix.ImaginaryOffDiagonalElements[Pos] * TmpMatrix.ImaginaryOffDiagonalElements[Pos]));
	      if (TmpNorm > PivotNorm)
		{
		  PivotNorm = TmpNorm;
		  PivotColumn = j;
		  PivotRow = i;
		  PivotPos = Pos;
		}
  	      Pos++; 
	    }
	  Pos += TmpMatrix.Increment;
	}
      if (PivotNorm == 0.0)
	return Complex();
      Pivot.Re = TmpMatrix.RealOffDiagonalElements[PivotPos];
      Pivot.Im = TmpMatrix.ImaginaryOffDiagonalElements[PivotPos];
      // add contribution to the pfaffian
      Pfaffian *= Pivot;
      Pivot = 1.0 / Pivot;
      // apply permutation to put pivot in position
      if (PivotColumn == (TmpMatrix.NbrRow - 2))
	Pos = PivotRow;
      else
	Pos = PivotColumn;
      if (PivotRow != (TmpMatrix.NbrRow - 2))
	{
	  TmpMatrix.SwapRowColumn(PivotRow, TmpMatrix.NbrRow - 2);
	  Pfaffian *= -1.0;
	}
      if (Pos != (TmpMatrix.NbrRow - 1))
	{
	  TmpMatrix.SwapRowColumn(Pos, TmpMatrix.NbrRow - 1);
	  Pfaffian *= -1.0;
	}
      // write Schur complement instead of the rest of the skew symmetric matrix
      Pos = 0;
      Complex Factor1;
      Complex Factor2;
      ReducedK = k - 2;
      int Pos1 = ReducedK - 1;
      int Pos2;
      TmpInc = TmpMatrix.NbrColumn + TmpMatrix.Increment - 2;
      for (int i = 0; i < ReducedK; ++i)
	{
	  Pos2 = Pos1 + TmpMatrix.NbrColumn + TmpMatrix.Increment - i - 2;
	  Factor1.Re = (Pivot.Re * TmpMatrix.RealOffDiagonalElements[Pos1]) - (Pivot.Im * TmpMatrix.ImaginaryOffDiagonalElements[Pos1]);
	  Factor1.Im = (Pivot.Im * TmpMatrix.RealOffDiagonalElements[Pos1]) + (Pivot.Re * TmpMatrix.ImaginaryOffDiagonalElements[Pos1]);
	  Factor2.Re = (Pivot.Re * TmpMatrix.RealOffDiagonalElements[Pos1 + 1]) - (Pivot.Im * TmpMatrix.ImaginaryOffDiagonalElements[Pos1 + 1]);
	  Factor2.Im = (Pivot.Im * TmpMatrix.RealOffDiagonalElements[Pos1 + 1]) + (Pivot.Re * TmpMatrix.ImaginaryOffDiagonalElements[Pos1 + 1]);
	  for (int j = i + 1; j < ReducedK; ++j)
	    {
	      TmpMatrix.RealOffDiagonalElements[Pos] += ((Factor2.Re * TmpMatrix.RealOffDiagonalElements[Pos2]) 
							 + (Factor1.Im * TmpMatrix.ImaginaryOffDiagonalElements[Pos2 + 1])
							 - (Factor2.Im * TmpMatrix.ImaginaryOffDiagonalElements[Pos2]) 
							 - (Factor1.Re * TmpMatrix.RealOffDiagonalElements[Pos2 + 1]));
	      TmpMatrix.ImaginaryOffDiagonalElements[Pos] += ((Factor2.Im * TmpMatrix.RealOffDiagonalElements[Pos2]) 
							      + (Factor2.Re * TmpMatrix.ImaginaryOffDiagonalElements[Pos2]) 
							      - (Factor1.Re * TmpMatrix.ImaginaryOffDiagonalElements[Pos2 + 1])
							      - (Factor1.Im * TmpMatrix.RealOffDiagonalElements[Pos2 + 1]));
	      ++Pos;
	      Pos2 += TmpInc - j;
	    }
	  Pos += 2 + TmpMatrix.Increment;
	  Pos1 += TmpInc - i;
	}
      TmpMatrix.Resize(ReducedK, ReducedK);
    }
  TmpMatrix.Resize(this->NbrColumn, this->NbrColumn);
  Pivot.Re = TmpMatrix.RealOffDiagonalElements[0];
  Pivot.Im = TmpMatrix.ImaginaryOffDiagonalElements[0];
  Pfaffian *= Pivot;
  return Pfaffian;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const ComplexSkewSymmetricMatrix& P)
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      int pos = (i - 1);
      for (int j = 0; j < i; j ++)
	{
	  Str << -(P.RealOffDiagonalElements[pos]) << "    ";
	  if (P.ImaginaryOffDiagonalElements[pos] > 0.0)
	    Str << -P.ImaginaryOffDiagonalElements[pos] << "i    ";
	  else
	    if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
	      Str << "+" << -P.ImaginaryOffDiagonalElements[pos] << "i    ";
	    else
	      Str << "    ";
	  pos += (P.NbrRow - j - 2) + P.Increment;
	}
      Str << "0    ";
      pos++;
      for (int j = i + 1; j < P.NbrRow; j++)
	{
	  Str << P.RealOffDiagonalElements[pos] << "    ";
	  if (P.ImaginaryOffDiagonalElements[pos] < 0.0)
	    Str << P.ImaginaryOffDiagonalElements[pos] << "i    ";
	  else
	    if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
	      Str << "+" << P.ImaginaryOffDiagonalElements[pos] << "i    ";
	    else
	      Str << "    ";
	  ++pos;
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

MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexSkewSymmetricMatrix& P)
{
  Str << "{";
  for (int i = 0; i < P.NbrRow; i++)
    {
      Str << "{";
      int pos = (i - 1);
      for (int j = 0; j < i; j ++)
	{
	  if ((P.RealOffDiagonalElements[pos] != 0) || (P.ImaginaryOffDiagonalElements[pos] == 0))
	    Str << -P.RealOffDiagonalElements[pos];
	  if (P.ImaginaryOffDiagonalElements[pos] > 0.0)
	    Str << -P.ImaginaryOffDiagonalElements[pos] << "I";
	  else
	    if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
	      Str << "+" << -P.ImaginaryOffDiagonalElements[pos] << "I";
	  Str << ",";
	  pos += (P.NbrRow - j - 2) + P.Increment;
	}
      Str << "0";
      if (i != (P.NbrRow - 1))
	{
	  Str << ",";	  
	  pos++;
	  for (int j = i + 1; j < (P.NbrRow - 1); j++)
	    {
	      if ((P.RealOffDiagonalElements[pos] != 0) || (P.ImaginaryOffDiagonalElements[pos] == 0))
		Str << P.RealOffDiagonalElements[pos];
	      if (P.ImaginaryOffDiagonalElements[pos] < 0.0)
		Str << P.ImaginaryOffDiagonalElements[pos] << "I";
	      else
		if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
		  Str << "+" << P.ImaginaryOffDiagonalElements[pos] << "I";
	      Str << ",";
	      ++pos;
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

#endif
