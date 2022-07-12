////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of complex symmetric matrix                   //
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


#include "Matrix/ComplexDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#ifdef USE_HILBERT_SPACE
#include "HilbertSpace/SubspaceSpaceConverter.h"
#endif

#include <stdlib.h>


using std::endl;


// default constructor
//

ComplexDiagonalMatrix::ComplexDiagonalMatrix() 
{
  this->DiagonalElements = 0;
  this->DiagonalGarbageFlag =  0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Symmetric | Matrix::Diagonal;
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

ComplexDiagonalMatrix::ComplexDiagonalMatrix(int dimension, bool zero) 
{
  this->DiagonalGarbageFlag =  new int;
  *(this->DiagonalGarbageFlag) = 1;
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Symmetric | Matrix::Diagonal;
  this->DiagonalElements = new Complex[this->NbrRow];
  if (zero == true)
    {
      for (int i = 0; i < this->NbrRow; i++)
	this->DiagonalElements[i] = 0.0;
    }
}

// constructor from matrix elements (without duplicating datas)
//
// diagonal = pointer to diagonal element array
// dimension = matrix dimension

ComplexDiagonalMatrix::ComplexDiagonalMatrix(Complex* diagonal, int dimension) 
{
  this->DiagonalElements = diagonal;
  this->DiagonalGarbageFlag =  new int;
  *(this->DiagonalGarbageFlag) = 1;
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Symmetric | Matrix::Diagonal;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

ComplexDiagonalMatrix::ComplexDiagonalMatrix(const ComplexDiagonalMatrix& M) 
{
  this->DiagonalElements = M.DiagonalElements;
  this->DiagonalGarbageFlag = M.DiagonalGarbageFlag;
  if (this->DiagonalGarbageFlag != 0)
    (*(this->DiagonalGarbageFlag))++;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Symmetric | Matrix::Diagonal;
}

// constructor from a full complex matrix
// M = matrix to copy
// isDiagonal = returns whether off-diagonal elements were non-zero
// tolerance = maximal value of off-diagonal elements before isDiagonal is set false
//
ComplexDiagonalMatrix::ComplexDiagonalMatrix(const ComplexMatrix& M, bool &isDiagonal, double tolerance) 
{
  this->DiagonalGarbageFlag =  new int;
  *(this->DiagonalGarbageFlag) = 1;  
  this->NbrRow = M.GetNbrRow();
  this->NbrColumn = M.GetNbrColumn();
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Symmetric | Matrix::Diagonal;
  int dim = (NbrRow<NbrColumn?NbrRow:NbrColumn);
  this->DiagonalElements = new Complex[dim];
  for (int i = 0; i < dim; i++)
    M.GetMatrixElement(i,i,this->DiagonalElements[i]);
  isDiagonal=true;
  Complex Tmp;
  for (int i = 0; (i < M.GetNbrColumn())&(isDiagonal); ++i)
    for (int j = 0; (j < M.GetNbrRow())&(isDiagonal); ++j)
      {
	M.GetMatrixElement(j,i,Tmp);
	if ((i!=j)&&(Norm(Tmp)>tolerance)) isDiagonal=false;
      }
}

// destructor
//

ComplexDiagonalMatrix::~ComplexDiagonalMatrix() 
{
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

ComplexDiagonalMatrix& ComplexDiagonalMatrix::operator = (const ComplexDiagonalMatrix& M) 
{
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
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* ComplexDiagonalMatrix::Clone ()
{
  return ((Matrix*) new ComplexDiagonalMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexDiagonalMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i != j) || (i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  this->DiagonalElements[i].Re = x;
  this->DiagonalElements[i].Im = 0.0;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element
void ComplexDiagonalMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  if ((i != j) || (i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  this->DiagonalElements[i] = x;
}

// set a matrix element
//
// i = line position
// j = column position
// real = real part of new value for matrix element
// imag = imaginary part of new value for matrix element
void ComplexDiagonalMatrix::SetMatrixElement(int i, int j, double real, double imag)
{
  if ((i != j) || (i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  this->DiagonalElements[i].Re = real;
  this->DiagonalElements[i].Im = imag;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void ComplexDiagonalMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i != j) || (i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  this->DiagonalElements[i].Re += x;
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void ComplexDiagonalMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  if ((i != j) || (i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  this->DiagonalElements[i] += x;  
}

// get reference of a given matrix element (supposing i == j)
//
// i = line position
// j = column position
// return value = reference om matrix elememt

double& ComplexDiagonalMatrix::operator () (int i, int j)
{
  return this->DiagonalElements[i].Re;
}

// get reference of a given matrix diagonal element
//
// i = line position
// return value = reference om matrix elememt

Complex& ComplexDiagonalMatrix::operator [] (int i)
{
  return this->DiagonalElements[i];
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexDiagonalMatrix::Resize (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  if (nbrRow <= this->TrueNbrRow)
    {
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      return;
    }
  Complex* TmpDiag = new Complex [nbrRow];
  for (int i = 0; i < this->NbrRow; i++)
    TmpDiag [i] = this->DiagonalElements[i];
  for (int i = this->NbrRow; i < nbrRow; i++)
    TmpDiag [i]  = 0.0;
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
  this->DiagonalElements = TmpDiag;
  this->DiagonalGarbageFlag =  new int;
  *(this->DiagonalGarbageFlag) = 1;
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexDiagonalMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  if (nbrRow <= this->TrueNbrRow)
    {
      if (this->NbrRow < nbrRow)
	{
	  for (int i = this->NbrRow; i < nbrRow; i++)
	    this->DiagonalElements[i] = 0.0;
	}
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      return;
    }
  Complex* TmpDiag = new Complex [nbrRow];
  for (int i = 0; i < this->NbrRow; i++)
    TmpDiag [i] = this->DiagonalElements[i];
  for (int i = this->NbrRow; i < nbrRow; i++)
    TmpDiag [i]  = 0.0;
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
  this->DiagonalElements = TmpDiag;
  this->DiagonalGarbageFlag =  new int;
  *(this->DiagonalGarbageFlag) = 1;
}

#ifdef USE_HILBERT_SPACE
// project matrix into a given subspace
//
// subspace = reference on subspace structure
// return value = pointer to projected matrix

Matrix* ComplexDiagonalMatrix::Project (SubspaceSpaceConverter& subspace)
{
  ComplexDiagonalMatrix* TmpM = new ComplexDiagonalMatrix (subspace.SubspaceDimension);
  for (int i = 0; i < subspace.SubspaceDimension; i++)
    TmpM->DiagonalElements[i] = this->DiagonalElements[subspace.SubspaceSpaceConverterArray[i]];
  return TmpM;
}
#endif

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

ComplexDiagonalMatrix operator + (const ComplexDiagonalMatrix& M1, const ComplexDiagonalMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexDiagonalMatrix();
  Complex* Diagonal = new Complex [M1.NbrRow];
  for (int i = 0; i < M1.NbrRow; i++)
    {
      Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
    }
  return ComplexDiagonalMatrix(Diagonal, M1.NbrRow);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexDiagonalMatrix operator - (const ComplexDiagonalMatrix& M1, const ComplexDiagonalMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexDiagonalMatrix();
  Complex* Diagonal = new Complex [M1.NbrRow];
  for (int i = 0; i < M1.NbrRow; i++)
    {
      Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
    }
  return ComplexDiagonalMatrix(Diagonal, M1.NbrRow);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexDiagonalMatrix operator * (const ComplexDiagonalMatrix& M, double x) 
{
  Complex* Diagonal = new Complex [M.NbrRow];
  for (int i = 0; i < M.NbrRow; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
    }
  return ComplexDiagonalMatrix(Diagonal, M.NbrRow);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = complex number to use
// return value = product result

ComplexDiagonalMatrix operator * (const ComplexDiagonalMatrix& M, const Complex& x) 
{
  Complex* Diagonal = new Complex [M.NbrRow];
  for (int i = 0; i < M.NbrRow; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
    }
  return ComplexDiagonalMatrix(Diagonal, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexDiagonalMatrix operator * (double x, const ComplexDiagonalMatrix& M) 
{
  Complex* Diagonal = new Complex [M.NbrRow];
  for (int i = 0; i < M.NbrRow; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
    }
  return ComplexDiagonalMatrix(Diagonal, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = complexnumber to use
// return value = product result

ComplexDiagonalMatrix operator * (const Complex &x, const ComplexDiagonalMatrix& M) 
{
  Complex* Diagonal = new Complex [M.NbrRow];
  for (int i = 0; i < M.NbrRow; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
    }
  return ComplexDiagonalMatrix(Diagonal, M.NbrRow);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

ComplexDiagonalMatrix operator / (const ComplexDiagonalMatrix& M, double x) 
{
  x = 1.0 / x;
  Complex* Diagonal = new Complex [M.NbrRow];
  for (int i = 0; i < M.NbrRow; i++)
    {
      Diagonal[i] = (Complex)M.DiagonalElements[i] * x;
    }
  return ComplexDiagonalMatrix(Diagonal, M.NbrRow);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = complex number to use
// return value = division result

ComplexDiagonalMatrix operator / (const ComplexDiagonalMatrix& M, const Complex &x) 
{
  Complex y = (Complex)1.0 / x;
  Complex* Diagonal = new Complex [M.NbrRow];
  for (int i = 0; i < M.NbrRow; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * y;
    }
  return ComplexDiagonalMatrix(Diagonal, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

ComplexDiagonalMatrix& ComplexDiagonalMatrix::operator += (const RealDiagonalMatrix& M) 
{
  if (this->NbrRow < M.NbrRow)
    return *this;
  for (int i = 0; i < M.NbrRow; i++)
    {
      this->DiagonalElements[i] += M.DiagonalElements[i];
    }
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

ComplexDiagonalMatrix& ComplexDiagonalMatrix::operator -= (const RealDiagonalMatrix& M) 
{
  if (this->NbrRow < M.NbrRow)
    return *this;
  for (int i = 0; i < M.NbrRow; i++)
    {
      this->DiagonalElements[i] -= M.DiagonalElements[i];
    }
  return *this;
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

ComplexDiagonalMatrix& ComplexDiagonalMatrix::operator += (const ComplexDiagonalMatrix& M) 
{
  if (this->NbrRow < M.NbrRow)
    return *this;
  for (int i = 0; i < M.NbrRow; i++)
    {
      this->DiagonalElements[i] += M.DiagonalElements[i];
    }
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

ComplexDiagonalMatrix& ComplexDiagonalMatrix::operator -= (const ComplexDiagonalMatrix& M) 
{
  if (this->NbrRow < M.NbrRow)
    return *this;
  for (int i = 0; i < M.NbrRow; i++)
    {
      this->DiagonalElements[i] -= M.DiagonalElements[i];
    }
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexDiagonalMatrix& ComplexDiagonalMatrix::operator *= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] *= x;
    }
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexDiagonalMatrix& ComplexDiagonalMatrix::operator /= (double x)
{
  if (this->NbrRow == 0)
    return *this;
  x = 1.0 / x;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] *= x;
    }
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexDiagonalMatrix& ComplexDiagonalMatrix::operator *= (Complex &x) 
{
  if (this->NbrRow == 0)
    return *this;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] *= x;
    }
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexDiagonalMatrix& ComplexDiagonalMatrix::operator /= (Complex &x)
{
  if (this->NbrRow == 0)
    return *this;
  Complex y = (Complex)1.0/x;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] *= y;
    }
  return *this;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ComplexDiagonalMatrix::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  Complex x = 0.0;
  double TmpRe, TmpIm;
  if ((V1.Dimension != this->NbrRow) || (V2.Dimension != this->NbrColumn))
    return x;
  for (int i = 0; i < this->NbrRow ; i++)
    {
      TmpRe=this->DiagonalElements[i].Re * V2.Components[i].Re - this->DiagonalElements[i].Im * V2.Components[i].Im;
      TmpIm=this->DiagonalElements[i].Im * V2.Components[i].Re + this->DiagonalElements[i].Re * V2.Components[i].Im;
      x.Re += V1.Components[i].Re * TmpRe + V1.Components[i].Im * TmpIm;
      x.Im += V1.Components[i].Re * TmpIm - V1.Components[i].Im * TmpRe;
    }
  return x;
}

// evaluate matrix trace
//
// return value = matrix trace 

Complex ComplexDiagonalMatrix::Trace () 
{
  if (this->NbrRow == 0)
    return 0.0;
  Complex x = this->DiagonalElements[0];
  for (int i = 1; i < this->NbrRow; i++)
    {
      x += this->DiagonalElements[i];
    }
  return x;
}

// evaluate matrix determinant
//
// return value = matrix determinant 

Complex ComplexDiagonalMatrix::Determinant () 
{
  if (this->NbrRow == 0)
    return 0.0;
  Complex x = this->DiagonalElements[0];
  for (int i = 1; i < this->NbrRow; i++)
    {
      x *= this->DiagonalElements[i];
    }
  return x;
}

// Sort Matrix such that diagnonal elements are sort in decreasing order
//
// return value = reference on current Matrix

ComplexDiagonalMatrix& ComplexDiagonalMatrix::SortMatrixDownOrder()
{
  int ReducedDim = this->NbrColumn - 2;
  Complex tmp;
  int MinPos;
  Complex MinValue;
  for (int i = 0; i <= ReducedDim; i++)
    {
      MinPos = this->NbrColumn - 1;
      MinValue = this->DiagonalElements[MinPos];
      for (int j = ReducedDim; j >= i; j--)
	if (this->DiagonalElements[j].Re > MinValue.Re)
	{
	  MinValue = this->DiagonalElements[j];
	  MinPos = j;
	}
      tmp = this->DiagonalElements[i];
      this->DiagonalElements[i] = MinValue;
      this->DiagonalElements[MinPos] = tmp;
    }
  return *this;
}
 
// Sort Matrix such that diagnonal elements are sort in decreasing order
// and apply corresponding transformation to column of a given complex matrix 
//
// matrix = matrix on which transformation has to be applied
// return value = reference on current Matrix

ComplexDiagonalMatrix& ComplexDiagonalMatrix::SortMatrixDownOrder(ComplexMatrix& matrix)
{
  int ReducedDim = this->NbrColumn - 2;
  ComplexVector TmpV;
  Complex tmp;
  int MinPos;
  Complex MinValue;
  for (int i = 0; i <= ReducedDim; i++)
    {
      MinPos = this->NbrColumn - 1;
      MinValue = this->DiagonalElements[MinPos];
      for (int j = ReducedDim; j >= i; j--)
	if (this->DiagonalElements[j].Re > MinValue.Re)
	  {
	    MinValue = this->DiagonalElements[j];
	    MinPos = j;
	  }
      tmp = this->DiagonalElements[i];
      this->DiagonalElements[i] = MinValue;
      this->DiagonalElements[MinPos] = tmp;
      TmpV = matrix[i];
      matrix[i] = matrix[MinPos];
      matrix[MinPos] = TmpV;
    }
  return *this;
}

// Sort Matrix such that diagnonal elements are sort in increasing order
//
// return value = reference on current Matrix

ComplexDiagonalMatrix& ComplexDiagonalMatrix::SortMatrixUpOrder()
{
  int ReducedDim = this->NbrColumn - 2;
  Complex tmp;
  int MinPos;
  Complex MinValue;
  for (int i = 0; i <= ReducedDim; i++)
    {
      MinPos = this->NbrColumn - 1;
      MinValue = this->DiagonalElements[MinPos];
      for (int j = ReducedDim; j >= i; j--)
	if (this->DiagonalElements[j].Re < MinValue.Re)
	{
	  MinValue = this->DiagonalElements[j];
	  MinPos = j;
	}
      tmp = this->DiagonalElements[i];
      this->DiagonalElements[i] = MinValue;
      this->DiagonalElements[MinPos] = tmp;
    }
  return *this;
}


// Sort Matrix such that diagnonal elements are sort in increasing order
// and apply corresponding transformation to column of a given complex matrix 
//
// matrix = matrix on which transformation has to be applied
// return value = reference on current Matrix

ComplexDiagonalMatrix& ComplexDiagonalMatrix::SortMatrixUpOrder(ComplexMatrix& matrix)
{
  int ReducedDim = this->NbrColumn - 2;
  ComplexVector TmpV;
  Complex tmp;
  int MinPos;
  Complex MinValue;
  for (int i = 0; i <= ReducedDim; i++)
    {
      MinPos = this->NbrColumn - 1;
      MinValue = this->DiagonalElements[MinPos];
      for (int j = ReducedDim; j >= i; j--)
	if (this->DiagonalElements[j].Re < MinValue.Re)
	{
	  MinValue = this->DiagonalElements[j];
	  MinPos = j;
	}
      tmp = this->DiagonalElements[i];
      this->DiagonalElements[i] = MinValue;
      this->DiagonalElements[MinPos] = tmp;
      TmpV = matrix[i];
      matrix[i] = matrix[MinPos];
      matrix[MinPos] = TmpV;
    }
  return *this;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const ComplexDiagonalMatrix& P)
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      for (int j = 0; j < i; j ++)
	{
	  Str << "0    ";
	}
      Str << P.DiagonalElements[i].Re<<"+"<<P.DiagonalElements[i].Im << "i    ";
      for (int j = i + 1; j < P.NbrRow; j++)
	{
	  Str << "0    ";
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

MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexDiagonalMatrix& P)
{
  Str << "{";
  for (int i = 0; i < P.NbrRow; i++)
    {
      Str << "{";
      for (int j = 0; j < i; j ++)
	{
	  Str << "0,";
	}
      Str << P.DiagonalElements[i];
      if (i != (P.NbrRow - 1))
	{
	  Str << ",";	  
	  for (int j = i + 1; j < (P.NbrRow - 1); j++)
	    {
	      Str << "0,";
	    }
	  Str << "0},";
	}
      else
	Str << "}";
    }
  Str << "}";
  return Str;
}


#endif
