////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of complex tridiagonal hermitian matrix              //
//                                                                            //
//                        last modification : 18/01/2001                      //
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


#include "Matrix/ComplexTriDiagonalHermitianMatrix.h"
#include "Vector/ComplexVector.h"


using std::endl;


// default constructor
//

ComplexTriDiagonalHermitianMatrix::ComplexTriDiagonalHermitianMatrix() 
{
  this->DiagonalElements = 0;
  this->RealUpperDiagonalElements = 0;
  this->ImaginaryUpperDiagonalElements = 0;
  this->Flag.Initialize();
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->MatrixType = Matrix::ComplexElements | Matrix::TriDiagonal | Matrix::Hermitian;
}

// constructor from matrix elements (without duplicating datas)
//
// diagonal = pointer to diagonal element array
// realUpperDiagonal = pointer to real part of upper diagonal element
// imaginaryUpperDiagonal = pointer to imaginary part of upper diagonal element
// dimension = matrix dimension

ComplexTriDiagonalHermitianMatrix::ComplexTriDiagonalHermitianMatrix(double* diagonal, double* realUpperDiagonal, 
								     double* imaginaryUpperDiagonal, int dimension) 
{
  this->DiagonalElements = diagonal;
  this->RealUpperDiagonalElements = realUpperDiagonal;
  this->ImaginaryUpperDiagonalElements = imaginaryUpperDiagonal;
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->MatrixType = Matrix::ComplexElements | Matrix::TriDiagonal | Matrix::Hermitian;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

ComplexTriDiagonalHermitianMatrix::ComplexTriDiagonalHermitianMatrix(const ComplexTriDiagonalHermitianMatrix& M) 
{
  this->DiagonalElements = M.DiagonalElements;
  this->Flag = M.Flag;
  this->RealUpperDiagonalElements = M.RealUpperDiagonalElements;
  this->ImaginaryUpperDiagonalElements = M.ImaginaryUpperDiagonalElements;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::TriDiagonal | Matrix::Hermitian;
}

// copy constructor from a real tridiagonal symmetric matrix
//
// M = matrix to copy

ComplexTriDiagonalHermitianMatrix::ComplexTriDiagonalHermitianMatrix(const RealTriDiagonalSymmetricMatrix& M) 
{
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::TriDiagonal | Matrix::Hermitian;
  this->Flag.Initialize();
  if (this->NbrRow == 0)
    {
      this->DiagonalElements = 0;
      this->RealUpperDiagonalElements = 0;
      this->ImaginaryUpperDiagonalElements = 0;
    }
  else
    {
      int ReducedNbr = M.NbrRow - 1;
      this->RealUpperDiagonalElements = new double [ReducedNbr];
      this->ImaginaryUpperDiagonalElements = new double [ReducedNbr];
      for (int i = 0; i < ReducedNbr; i++)
	{
	  this->DiagonalElements[i] = M.DiagonalElements[i];      
	  this->RealUpperDiagonalElements[i] = M.UpperDiagonalElements[i];
	  this->ImaginaryUpperDiagonalElements[i] = 0.0;
	}
      this->DiagonalElements[ReducedNbr] = M.DiagonalElements[ReducedNbr];
    }
}

// destructor
//

ComplexTriDiagonalHermitianMatrix::~ComplexTriDiagonalHermitianMatrix() 
{
  if ((this->RealUpperDiagonalElements != 0) && (this->ImaginaryUpperDiagonalElements != 0) 
      && (this->DiagonalElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
	{
	  delete[] this->DiagonalElements;
	  delete[] this->RealUpperDiagonalElements;
	  delete[] this->ImaginaryUpperDiagonalElements;
	}
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

ComplexTriDiagonalHermitianMatrix& ComplexTriDiagonalHermitianMatrix::operator = (const ComplexTriDiagonalHermitianMatrix& M) 
{
  if ((this->RealUpperDiagonalElements != 0) && (this->ImaginaryUpperDiagonalElements != 0) 
      && (this->DiagonalElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
	{
	  delete[] this->DiagonalElements;
	  delete[] this->RealUpperDiagonalElements;
	  delete[] this->ImaginaryUpperDiagonalElements;
	}
  this->DiagonalElements = M.DiagonalElements;
  this->Flag = M.Flag;
  this->RealUpperDiagonalElements = M.RealUpperDiagonalElements;
  this->ImaginaryUpperDiagonalElements = M.ImaginaryUpperDiagonalElements;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  return *this;
}

// assignement from  a real tridiagonal symmetric matrix (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

ComplexTriDiagonalHermitianMatrix& ComplexTriDiagonalHermitianMatrix::operator = (const RealTriDiagonalSymmetricMatrix& M) 
{
  if ((this->RealUpperDiagonalElements != 0) && (this->ImaginaryUpperDiagonalElements != 0) 
      && (this->DiagonalElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
	{
	  delete[] this->DiagonalElements;
	  delete[] this->RealUpperDiagonalElements;
	  delete[] this->ImaginaryUpperDiagonalElements;
	}
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::TriDiagonal | Matrix::Hermitian;
  this->Flag.Initialize();
  if (this->NbrRow == 0)
    {
      this->DiagonalElements = 0;
      this->RealUpperDiagonalElements = 0;
      this->ImaginaryUpperDiagonalElements = 0;
    }
  else
    {
      int ReducedNbr = M.NbrRow - 1;
      this->RealUpperDiagonalElements = new double [ReducedNbr];
      this->ImaginaryUpperDiagonalElements = new double [ReducedNbr];
      for (int i = 0; i < ReducedNbr; i++)
	{
	  this->DiagonalElements[i] = M.DiagonalElements[i];      
	  this->RealUpperDiagonalElements[i] = M.UpperDiagonalElements[i];
	  this->ImaginaryUpperDiagonalElements[i] = 0.0;
	}
      this->DiagonalElements[ReducedNbr] = M.DiagonalElements[ReducedNbr];
    }
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* ComplexTriDiagonalHermitianMatrix::Clone ()
{
  return ((Matrix*) new ComplexTriDiagonalHermitianMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexTriDiagonalHermitianMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i == j) && (i < this->NbrRow))
    {
      this->DiagonalElements[i] = x;
    }
  else
    {
      j -= i;
      if ((j == 1) && (i < (this->NbrRow - 1)))
	{
	  this->RealUpperDiagonalElements[i] = x;
	  this->ImaginaryUpperDiagonalElements[i] = 0.0;	  
	}
      else
	if ((j == -1) && (i < this->NbrRow))
	  {
	    this->RealUpperDiagonalElements[i - 1] = x;
	    this->ImaginaryUpperDiagonalElements[i - 1] = 0.0;	  
	  }	
    }    
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexTriDiagonalHermitianMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  j -= i;
  if ((j == 1) && (i < (this->NbrRow - 1)))
    {
      this->RealUpperDiagonalElements[i] = x.Re;
      this->ImaginaryUpperDiagonalElements[i] = x.Im;	  
    }
  else
    if ((j == -1) && (i < this->NbrRow))
      {
	this->RealUpperDiagonalElements[i - 1] = x.Re;
	this->ImaginaryUpperDiagonalElements[i - 1] = -x.Im;	  
      }	
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void ComplexTriDiagonalHermitianMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i == j) && (i < this->NbrRow))
    {
      this->DiagonalElements[i] += x;
    }
  else
    {
      j -= i;
      if ((j == 1) && (i < (this->NbrRow - 1)))
	{
	  this->RealUpperDiagonalElements[i] += x;
	}
      else
	if ((j == -1) && (i < this->NbrRow))
	  {
	    this->RealUpperDiagonalElements[i - 1] += x;
	  }	
    }    
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void ComplexTriDiagonalHermitianMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  j -= i;
  if ((j == 1) && (i < (this->NbrRow - 1)))
    {
      this->RealUpperDiagonalElements[i] += x.Re;
      this->ImaginaryUpperDiagonalElements[i + 1] += x.Im;	  
    }
  else
    if ((j == -1) && (i < this->NbrRow))
      {
	this->RealUpperDiagonalElements[i - 1] += x.Re;
	this->ImaginaryUpperDiagonalElements[i - 1] -= x.Im;	  
      }	
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

ComplexTriDiagonalHermitianMatrix operator + (const ComplexTriDiagonalHermitianMatrix& M1, const ComplexTriDiagonalHermitianMatrix& M2) 
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexTriDiagonalHermitianMatrix();
  int ReducedNbr = M1.NbrRow - 1;
  double* Diagonal = new double [M1.NbrRow];
  double* RealUpperDiagonal = new double [ReducedNbr];
  double* ImaginaryUpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; ++i)
    {
      Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
      RealUpperDiagonal[i] = M1.RealUpperDiagonalElements[i] + M2.RealUpperDiagonalElements[i];      
      ImaginaryUpperDiagonal[i] = M1.ImaginaryUpperDiagonalElements[i] + M2.ImaginaryUpperDiagonalElements[i];      
    }
  Diagonal[ReducedNbr] = M1.DiagonalElements[ReducedNbr] + M2.DiagonalElements[ReducedNbr];
  return ComplexTriDiagonalHermitianMatrix(Diagonal, RealUpperDiagonal, ImaginaryUpperDiagonal, M1.NbrRow);
}

// add two matrices where the left one is a real tridiagonal symmetric matrix
//
// M1 = left matrix
// M2 = right matrix
// return value = sum of the two matrices

ComplexTriDiagonalHermitianMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const ComplexTriDiagonalHermitianMatrix& M2) 
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexTriDiagonalHermitianMatrix();
  int ReducedNbr = M1.NbrRow - 1;
  double* Diagonal = new double [M1.NbrRow];
  double* RealUpperDiagonal = new double [ReducedNbr];
  double* ImaginaryUpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; ++i)
    {
      Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
      RealUpperDiagonal[i] = M1.UpperDiagonalElements[i] + M2.RealUpperDiagonalElements[i];      
      ImaginaryUpperDiagonal[i] = M2.ImaginaryUpperDiagonalElements[i];      
    }
  Diagonal[ReducedNbr] = M1.DiagonalElements[ReducedNbr] + M2.DiagonalElements[ReducedNbr];
  return ComplexTriDiagonalHermitianMatrix(Diagonal, RealUpperDiagonal, ImaginaryUpperDiagonal, M1.NbrRow);
}

// add two matrices where the right one is a real tridiagonal symmetric matrix
//
// M1 = left matrix
// M2 = right matrix
// return value = sum of the two matrices

ComplexTriDiagonalHermitianMatrix operator + (const ComplexTriDiagonalHermitianMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2) 
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexTriDiagonalHermitianMatrix();
  int ReducedNbr = M1.NbrRow - 1;
  double* Diagonal = new double [M1.NbrRow];
  double* RealUpperDiagonal = new double [ReducedNbr];
  double* ImaginaryUpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; ++i)
    {
      Diagonal[i] = M2.DiagonalElements[i] + M1.DiagonalElements[i];
      RealUpperDiagonal[i] = M2.UpperDiagonalElements[i] + M1.RealUpperDiagonalElements[i];      
      ImaginaryUpperDiagonal[i] = M1.ImaginaryUpperDiagonalElements[i];      
    }
  Diagonal[ReducedNbr] = M1.DiagonalElements[ReducedNbr] + M2.DiagonalElements[ReducedNbr];
  return ComplexTriDiagonalHermitianMatrix(Diagonal, RealUpperDiagonal, ImaginaryUpperDiagonal, M1.NbrRow);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexTriDiagonalHermitianMatrix operator - (const ComplexTriDiagonalHermitianMatrix& M1, const ComplexTriDiagonalHermitianMatrix& M2) 
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexTriDiagonalHermitianMatrix();
  int ReducedNbr = M1.NbrRow - 1;
  double* Diagonal = new double [M1.NbrRow];
  double* RealUpperDiagonal = new double [ReducedNbr];
  double* ImaginaryUpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; ++i)
    {
      Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
      RealUpperDiagonal[i] = M1.RealUpperDiagonalElements[i] - M2.RealUpperDiagonalElements[i];      
      ImaginaryUpperDiagonal[i] = M1.ImaginaryUpperDiagonalElements[i] - M2.ImaginaryUpperDiagonalElements[i];      
    }
  Diagonal[ReducedNbr] = M1.DiagonalElements[ReducedNbr] - M2.DiagonalElements[ReducedNbr];
  return ComplexTriDiagonalHermitianMatrix(Diagonal, RealUpperDiagonal, ImaginaryUpperDiagonal, M1.NbrRow);
}

// substract two matrices where the left one is a real tridiagonal symmetric matrix
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexTriDiagonalHermitianMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, const ComplexTriDiagonalHermitianMatrix& M2) 
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexTriDiagonalHermitianMatrix();
  int ReducedNbr = M1.NbrRow - 1;
  double* Diagonal = new double [M1.NbrRow];
  double* RealUpperDiagonal = new double [ReducedNbr];
  double* ImaginaryUpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; ++i)
    {
      Diagonal[i] = M1.DiagonalElements[i] - M2.DiagonalElements[i];
      RealUpperDiagonal[i] = M1.UpperDiagonalElements[i] - M2.RealUpperDiagonalElements[i];      
      ImaginaryUpperDiagonal[i] = -M2.ImaginaryUpperDiagonalElements[i];      
    }
  Diagonal[ReducedNbr] = M1.DiagonalElements[ReducedNbr] - M2.DiagonalElements[ReducedNbr];
  return ComplexTriDiagonalHermitianMatrix(Diagonal, RealUpperDiagonal, ImaginaryUpperDiagonal, M1.NbrRow);
}

// substract two matrices where the right one is a real tridiagonal symmetric matrix
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexTriDiagonalHermitianMatrix operator - (const ComplexTriDiagonalHermitianMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2) 
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexTriDiagonalHermitianMatrix();
  int ReducedNbr = M1.NbrRow - 1;
  double* Diagonal = new double [M1.NbrRow];
  double* RealUpperDiagonal = new double [ReducedNbr];
  double* ImaginaryUpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; ++i)
    {
      Diagonal[i] = M1.DiagonalElements[i] - M2.DiagonalElements[i];
      RealUpperDiagonal[i] = M1.RealUpperDiagonalElements[i] - M2.UpperDiagonalElements[i];      
      ImaginaryUpperDiagonal[i] = M1.ImaginaryUpperDiagonalElements[i];      
    }
  Diagonal[ReducedNbr] = M1.DiagonalElements[ReducedNbr] - M2.DiagonalElements[ReducedNbr];
  return ComplexTriDiagonalHermitianMatrix(Diagonal, RealUpperDiagonal, ImaginaryUpperDiagonal, M1.NbrRow);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexTriDiagonalHermitianMatrix operator * (const ComplexTriDiagonalHermitianMatrix& M, double x) 
{
  int ReducedNbr = M.NbrRow - 1;
  double* Diagonal = new double [M.NbrRow];
  double* RealUpperDiagonal = new double [ReducedNbr];
  double* ImaginaryUpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; ++i)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
      RealUpperDiagonal[i] = M.RealUpperDiagonalElements[i] * x;      
      ImaginaryUpperDiagonal[i] = M.ImaginaryUpperDiagonalElements[i] * x;      
    }
  Diagonal[ReducedNbr] = M.DiagonalElements[ReducedNbr] * x;
  return ComplexTriDiagonalHermitianMatrix(Diagonal, RealUpperDiagonal, ImaginaryUpperDiagonal, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexTriDiagonalHermitianMatrix operator * (double x, const ComplexTriDiagonalHermitianMatrix& M) 
{
  int ReducedNbr = M.NbrRow - 1;
  double* Diagonal = new double [M.NbrRow];
  double* RealUpperDiagonal = new double [ReducedNbr];
  double* ImaginaryUpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; ++i)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
      RealUpperDiagonal[i] = M.RealUpperDiagonalElements[i] * x;      
      ImaginaryUpperDiagonal[i] = M.ImaginaryUpperDiagonalElements[i] * x;      
    }
  Diagonal[ReducedNbr] = M.DiagonalElements[ReducedNbr] * x;
  return ComplexTriDiagonalHermitianMatrix(Diagonal, RealUpperDiagonal, ImaginaryUpperDiagonal, M.NbrRow);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

ComplexTriDiagonalHermitianMatrix operator / (const ComplexTriDiagonalHermitianMatrix& M, double x) 
{
  x = 1.0 / x;
  int ReducedNbr = M.NbrRow - 1;
  double* Diagonal = new double [M.NbrRow];
  double* RealUpperDiagonal = new double [ReducedNbr];
  double* ImaginaryUpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; ++i)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
      RealUpperDiagonal[i] = M.RealUpperDiagonalElements[i] * x;      
      ImaginaryUpperDiagonal[i] = M.ImaginaryUpperDiagonalElements[i] * x;      
    }
  Diagonal[ReducedNbr] = M.DiagonalElements[ReducedNbr] * x;
  return ComplexTriDiagonalHermitianMatrix(Diagonal, RealUpperDiagonal, ImaginaryUpperDiagonal, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

ComplexTriDiagonalHermitianMatrix& ComplexTriDiagonalHermitianMatrix::operator += (const ComplexTriDiagonalHermitianMatrix& M) 
{
  if (this->NbrRow != M.NbrRow)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  for (int i = 0; i < ReducedNbr; ++i)
    {
      this->DiagonalElements[i] += M.DiagonalElements[i];
      this->RealUpperDiagonalElements[i] += M.RealUpperDiagonalElements[i];      
      this->ImaginaryUpperDiagonalElements[i] += M.ImaginaryUpperDiagonalElements[i];      
    }
  this->DiagonalElements[ReducedNbr] += DiagonalElements[ReducedNbr];
  return *this;
}

// add two matrices where the right one is a real tridiagonal symmetric matrix
//
// M = matrix to add to current matrix
// return value = reference on current matrix

ComplexTriDiagonalHermitianMatrix& ComplexTriDiagonalHermitianMatrix::operator += (const RealTriDiagonalSymmetricMatrix& M) 
{
  if (this->NbrRow != M.NbrRow)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  for (int i = 0; i < ReducedNbr; ++i)
    {
      this->DiagonalElements[i] += M.DiagonalElements[i];
      this->RealUpperDiagonalElements[i] += M.UpperDiagonalElements[i];      
    }
  this->DiagonalElements[ReducedNbr] += DiagonalElements[ReducedNbr];
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

ComplexTriDiagonalHermitianMatrix& ComplexTriDiagonalHermitianMatrix::operator -= (const ComplexTriDiagonalHermitianMatrix& M) 
{
  if (this->NbrRow != M.NbrRow)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  for (int i = 0; i < ReducedNbr; ++i)
    {
      this->DiagonalElements[i] -= M.DiagonalElements[i];
      this->RealUpperDiagonalElements[i] -= M.RealUpperDiagonalElements[i];      
      this->ImaginaryUpperDiagonalElements[i] -= M.ImaginaryUpperDiagonalElements[i];      
    }
  this->DiagonalElements[ReducedNbr] -= DiagonalElements[ReducedNbr];
  return *this;
}

// substract two matrices where the right one is a real tridiagonal symmetric matrix
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

ComplexTriDiagonalHermitianMatrix& ComplexTriDiagonalHermitianMatrix::operator -= (const RealTriDiagonalSymmetricMatrix& M) 
{
  if (this->NbrRow != M.NbrRow)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  for (int i = 0; i < ReducedNbr; ++i)
    {
      this->DiagonalElements[i] -= M.DiagonalElements[i];
      this->RealUpperDiagonalElements[i] -= M.UpperDiagonalElements[i];      
    }
  this->DiagonalElements[ReducedNbr] -= DiagonalElements[ReducedNbr];
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexTriDiagonalHermitianMatrix& ComplexTriDiagonalHermitianMatrix::operator *= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  int ReducedNbr = this->NbrRow - 1;
  for (int i = 0; i < ReducedNbr; ++i)
    {
      this->DiagonalElements[i] *= x;
      this->RealUpperDiagonalElements[i] *= x;      
      this->ImaginaryUpperDiagonalElements[i] *= x;      
    }
  this->DiagonalElements[ReducedNbr] *= x;
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexTriDiagonalHermitianMatrix& ComplexTriDiagonalHermitianMatrix::operator /= (double x)
{
  if (this->NbrRow == 0)
    return *this;
  x = 1.0 / x;
  int ReducedNbr = this->NbrRow - 1;
  for (int i = 0; i < ReducedNbr; ++i)
    {
      this->DiagonalElements[i] *= x;
      this->RealUpperDiagonalElements[i] *= x;      
      this->ImaginaryUpperDiagonalElements[i] *= x;      
    }
  this->DiagonalElements[ReducedNbr] *= x;
  return *this;
}

// evaluate matrix trace
//
// return value = matrix trace 

double ComplexTriDiagonalHermitianMatrix::Tr () 
{
  if (this->NbrRow == 0)
    return 0.0;
  double x = this->DiagonalElements[0];
  for (int i = 1; i < this->NbrRow; ++i)
    {
      x += this->DiagonalElements[i];
    }
  return x;
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double ComplexTriDiagonalHermitianMatrix::Det () 
{
  if (this->NbrRow == 0)
    return 0.0;
  double d0 = this->DiagonalElements[0];
  if (this->NbrRow == 1)
    return d0;
  double d1 = d0 * this->DiagonalElements[0] - (this->RealUpperDiagonalElements[0] * this->RealUpperDiagonalElements[0] + 
						this->ImaginaryUpperDiagonalElements[0] * this->ImaginaryUpperDiagonalElements[0]);
  if (this->NbrRow == 2)
    return d0;
  double d = 0.0;
  for (int i = 2; i < this->NbrRow; i++)
    {
      d = this->DiagonalElements[i] * d1 - (this->RealUpperDiagonalElements[i] * this->RealUpperDiagonalElements[i] +
					    this->ImaginaryUpperDiagonalElements[i] * this->ImaginaryUpperDiagonalElements[i]) * d0;
      d0 = d1;
      d1 = d;
    }
  return d1;
}

#ifdef USE_POLYNOMIAL

// return matrix characteritic equation
//
// return value =  reference one polynomial corresponding to matrix characteritic equation

Polynomial& ComplexTriDiagonalHermitianMatrix::CharacteristicEquation()
{
  double* P0 = new double [this->NbrRow + 1];
  double* P1 = new double [this->NbrRow + 1];
  double* TmpP;
  P0[1] = -1;
  P0[0] = this->DiagonalElements[0];
  P1[2] = 1;
  P1[1] = - this->DiagonalElements[0] - this->DiagonalElements[1];
  P1[0] = - this->RealUpperDiagonalElements[0] * this->RealUpperDiagonalElements[0] 
    - this->ImaginaryUpperDiagonalElements[0] * this->ImaginaryUpperDiagonalElements[0] +  
    this->DiagonalElements[0] * this->DiagonalElements[1];

  for (int i = 2; i < this->NbrRow; ++i)
    {
      int j = (i  - 1);
      double fac = - this->RealUpperDiagonalElements[j] * this->RealUpperDiagonalElements[j] -
	this->ImaginaryUpperDiagonalElements[j] * this->ImaginaryUpperDiagonalElements[j];
      for (int k = 0; k < i; ++k)
	P0[k] *= fac;
      P0[0] += this->DiagonalElements[i] * P1[0];
      for (int k = 1; k < i; ++k)
	{
	  P0[k] += -P1[k - 1];
	  P0[k] += this->DiagonalElements[i] * P1[k];
	} 
      P0[i] = -P1[i - 1];
      P0[i] += this->DiagonalElements[i] * P1[i];
      P0[i + 1] = -P1[i];
      TmpP = P0;
      P0 = P1;
      P1 = TmpP;
    }
  Polynomial* P = new Polynomial (this->NbrRow, P1, true);
  delete[] P0;
  return *P;
}

#endif

// evaluate a normalized eigenvector for a given eigenvalue (supposing the eigenvalue is non-degenerate)
//
// eigenvalue = eigenvalue to use
// eigenvector = vector where the eigenvector has to be stored
// return value = reference on eigenvector

ComplexVector& ComplexTriDiagonalHermitianMatrix::Eigenvector(double eigenvalue, ComplexVector& eigenvector)
{
  double Norm = 1.0;
  eigenvector.Components[0].Re = 1.0;
  eigenvector.Components[0].Im = 0.0;
  if (this->RealUpperDiagonalElements[0] == 0.0)
    {
      eigenvector.Components[1].Re = 0.0;
      eigenvector.Components[1].Im = 0.0;
    }
  else
    {
      double fac = 1.0 / (this->RealUpperDiagonalElements[0] * this->RealUpperDiagonalElements[0]
			  + this->ImaginaryUpperDiagonalElements[0] * this->ImaginaryUpperDiagonalElements[0]);
      eigenvector.Components[1].Re = (eigenvalue - this->DiagonalElements[0]) * this->RealUpperDiagonalElements[0] * fac;
      eigenvector.Components[1].Im = (this->DiagonalElements[0] - eigenvalue) * this->ImaginaryUpperDiagonalElements[0] * fac;
    }
  for (int i = 2; i < this->NbrRow; i++)
    {
      if ((this->RealUpperDiagonalElements[i - 1] == 0.0) && (this->ImaginaryUpperDiagonalElements[i - 1] == 0.0))
	{
	  eigenvector.Components[i].Re = 0.0;
	  eigenvector.Components[i].Im = 0.0;
	}
      else
	{
	  double fac = 1.0 / (this->RealUpperDiagonalElements[i - 1] * this->RealUpperDiagonalElements[i - 1]
			      + this->ImaginaryUpperDiagonalElements[i - 1] * this->ImaginaryUpperDiagonalElements[i - 1]);
	  eigenvector.Components[i].Re = (((eigenvalue - this->DiagonalElements[i - 1]) * eigenvector.Components[i - 1].Re
					   - this->ImaginaryUpperDiagonalElements[i - 2] * eigenvector.Components[i - 1].Im 
					    - this->RealUpperDiagonalElements[i - 2] * eigenvector.Components[i - 2].Re) * fac);
	  eigenvector.Components[i].Im = (((eigenvalue - this->DiagonalElements[i - 1]) * eigenvector.Components[i - 1].Im
						 + this->ImaginaryUpperDiagonalElements[i - 2] * eigenvector.Components[i - 2].Re 
						 - this->RealUpperDiagonalElements[i - 2] * eigenvector.Components[i - 2].Im) * fac);
	  Norm += (eigenvector.Components[i].Re * eigenvector.Components[i].Re 
		   + eigenvector.Components[i].Im * eigenvector.Components[i].Im);
	}
    }      
  Norm = 1.0 / sqrt (Norm);
  for (int i = 0; i < eigenvector.Dimension; ++i)
    {
      eigenvector.Components[i].Re *= Norm;
      eigenvector.Components[i].Im *= Norm;
    }
  return eigenvector;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const ComplexTriDiagonalHermitianMatrix& P)
{
  for (int i = 0; i < P.NbrRow; ++i)
    {
      int j = 0;
      for (; j < (i-1); ++j)
	Str << "0    ";
      if (i > 0)
	{
	  Str << P.RealUpperDiagonalElements[i - 1];
	  if (P.ImaginaryUpperDiagonalElements[i - 1] > 0)
	    Str << -P.ImaginaryUpperDiagonalElements[i - 1] << "i    ";
	  else
	    Str << "+" << -P.ImaginaryUpperDiagonalElements[i - 1] << "i    ";
	  j++;
	}
      Str << P.DiagonalElements[i] << "    ";
      j++;
      if (i < (P.NbrRow -1))
	{
	  Str << P.RealUpperDiagonalElements[i];
	  if (P.ImaginaryUpperDiagonalElements[i] > 0)
	    Str << "+" << P.ImaginaryUpperDiagonalElements[i] << "i    ";
	  else
	    Str << P.RealUpperDiagonalElements[i * 2 + 1] << "i    ";
	  j++;
	}
      for (; j < P.NbrColumn; j++)
	Str << "0    ";
      Str << endl;
    }
  return Str;
}
