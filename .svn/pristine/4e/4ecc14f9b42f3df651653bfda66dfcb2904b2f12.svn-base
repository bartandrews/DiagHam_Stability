////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of real symmetric matrix                      //
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


#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#ifdef USE_HILBERT_SPACE
#include "HilbertSpace/SubspaceSpaceConverter.h"
#endif

#include <stdlib.h>


using std::endl;


// default constructor
//

RealDiagonalMatrix::RealDiagonalMatrix() 
{
  this->DiagonalElements = 0;
  this->DiagonalGarbageFlag =  0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal;
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

RealDiagonalMatrix::RealDiagonalMatrix(int dimension, bool zero) 
{
  this->DiagonalGarbageFlag =  new int;
  *(this->DiagonalGarbageFlag) = 1;
   this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal;
  this->DiagonalElements = new double [this->NbrRow];
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

RealDiagonalMatrix::RealDiagonalMatrix(double* diagonal, int dimension) 
{
  this->DiagonalElements = diagonal;
  this->DiagonalGarbageFlag =  new int;
  *(this->DiagonalGarbageFlag) = 1;
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

RealDiagonalMatrix::RealDiagonalMatrix(const RealDiagonalMatrix& M) 
{
  this->DiagonalElements = M.DiagonalElements;
  this->DiagonalGarbageFlag = M.DiagonalGarbageFlag;
  if (this->DiagonalGarbageFlag != 0)
    (*(this->DiagonalGarbageFlag))++;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal;
}

// destructor
//

RealDiagonalMatrix::~RealDiagonalMatrix() 
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

RealDiagonalMatrix& RealDiagonalMatrix::operator = (const RealDiagonalMatrix& M) 
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

Matrix* RealDiagonalMatrix::Clone ()
{
  return ((Matrix*) new RealDiagonalMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void RealDiagonalMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i != j) || (i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  this->DiagonalElements[i] = x;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element
void RealDiagonalMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  return;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void RealDiagonalMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i != j) || (i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  this->DiagonalElements[i] += x;
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void RealDiagonalMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  return;
}

// get reference of a given matrix element (supposing i == j)
//
// i = line position
// j = column position
// return value = reference om matrix elememt

double& RealDiagonalMatrix::operator () (int i, int j)
{
  return this->DiagonalElements[i];
}

// get reference of a given matrix diagonal element
//
// i = line position
// return value = reference om matrix elememt

double& RealDiagonalMatrix::operator [] (int i)
{
  return this->DiagonalElements[i];
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void RealDiagonalMatrix::Resize (int nbrRow, int nbrColumn)
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

void RealDiagonalMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
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
  double* TmpDiag = new double [nbrRow];
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

Matrix* RealDiagonalMatrix::Project (SubspaceSpaceConverter& subspace)
{
  RealDiagonalMatrix* TmpM = new RealDiagonalMatrix (subspace.SubspaceDimension);
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

RealDiagonalMatrix operator + (const RealDiagonalMatrix& M1, const RealDiagonalMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return RealDiagonalMatrix();
  double* Diagonal = new double [M1.NbrRow];
  for (int i = 0; i < M1.NbrRow; i++)
    {
      Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
    }
  return RealDiagonalMatrix(Diagonal, M1.NbrRow);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

RealDiagonalMatrix operator - (const RealDiagonalMatrix& M1, const RealDiagonalMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return RealDiagonalMatrix();
  double* Diagonal = new double [M1.NbrRow];
  for (int i = 0; i < M1.NbrRow; i++)
    {
      Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
    }
  return RealDiagonalMatrix(Diagonal, M1.NbrRow);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealDiagonalMatrix operator * (const RealDiagonalMatrix& M, double x) 
{
  double* Diagonal = new double [M.NbrRow];
  for (int i = 0; i < M.NbrRow; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
    }
  return RealDiagonalMatrix(Diagonal, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealDiagonalMatrix operator * (double x, const RealDiagonalMatrix& M) 
{
  double* Diagonal = new double [M.NbrRow];
  for (int i = 0; i < M.NbrRow; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
    }
  return RealDiagonalMatrix(Diagonal, M.NbrRow);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

RealDiagonalMatrix operator / (const RealDiagonalMatrix& M, double x) 
{
  x = 1.0 / x;
  double* Diagonal = new double [M.NbrRow];
  for (int i = 0; i < M.NbrRow; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
    }
  return RealDiagonalMatrix(Diagonal, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

RealDiagonalMatrix& RealDiagonalMatrix::operator += (const RealDiagonalMatrix& M) 
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

RealDiagonalMatrix& RealDiagonalMatrix::operator -= (const RealDiagonalMatrix& M) 
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

RealDiagonalMatrix& RealDiagonalMatrix::operator *= (double x) 
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

RealDiagonalMatrix& RealDiagonalMatrix::operator /= (double x)
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

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

double RealDiagonalMatrix::MatrixElement (RealVector& V1, RealVector& V2)
{
  double x = 0.0;
  if ((V1.Dimension != this->NbrRow) || (V2.Dimension != this->NbrColumn))
    return x;
  for (int i = 0; i < this->NbrRow ; i++)
    {
      x += V1.Components[i] * this->DiagonalElements[i] * V2.Components[i];
    }
  return x;
}

// evaluate matrix trace
//
// return value = matrix trace 

double RealDiagonalMatrix::Tr () 
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

double RealDiagonalMatrix::Det () 
{
  if (this->NbrRow == 0)
    return 0.0;
  double x = this->DiagonalElements[0];
  for (int i = 1; i < this->NbrRow; i++)
    {
      x *= this->DiagonalElements[i];
    }
  return x;
}

// Sort Matrix such that diagnonal elements are sort in decreasing order
//
// return value = reference on current Matrix

RealDiagonalMatrix& RealDiagonalMatrix::SortMatrixDownOrder()
{
  int ReducedDim = this->NbrColumn - 2;
  double tmp;
  int MinPos;
  double MinValue;
  for (int i = 0; i <= ReducedDim; i++)
    {
      MinPos = this->NbrColumn - 1;
      MinValue = this->DiagonalElements[MinPos];
      for (int j = ReducedDim; j >= i; j--)
	if (this->DiagonalElements[j] > MinValue)
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
// and apply corresponding transformation to column of a given real matrix 
//
// matrix = matrix on which transformation has to be applied
// return value = reference on current Matrix

RealDiagonalMatrix& RealDiagonalMatrix::SortMatrixDownOrder(RealMatrix& matrix)
{
  int ReducedDim = this->NbrColumn - 2;
  RealVector TmpV;
  double tmp;
  int MinPos;
  double MinValue;
  for (int i = 0; i <= ReducedDim; i++)
    {
      MinPos = this->NbrColumn - 1;
      MinValue = this->DiagonalElements[MinPos];
      for (int j = ReducedDim; j >= i; j--)
	if (this->DiagonalElements[j] > MinValue)
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

// Sort Matrix such that diagnonal elements are sort in decreasing order
// and apply corresponding transformation to column of a given complex matrix 
//
// matrix = matrix on which transformation has to be applied
// return value = reference on current Matrix

RealDiagonalMatrix& RealDiagonalMatrix::SortMatrixDownOrder(ComplexMatrix& matrix)
{
  int ReducedDim = this->NbrColumn - 2;
  RealVector TmpV;
  double tmp;
  int MinPos;
  double MinValue;
  for (int i = 0; i <= ReducedDim; i++)
    {
      MinPos = this->NbrColumn - 1;
      MinValue = this->DiagonalElements[MinPos];
      for (int j = ReducedDim; j >= i; j--)
	if (this->DiagonalElements[j] > MinValue)
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

RealDiagonalMatrix& RealDiagonalMatrix::SortMatrixUpOrder()
{
  int ReducedDim = this->NbrColumn - 2;
  double tmp;
  int MinPos;
  double MinValue;
  for (int i = 0; i <= ReducedDim; i++)
    {
      MinPos = this->NbrColumn - 1;
      MinValue = this->DiagonalElements[MinPos];
      for (int j = ReducedDim; j >= i; j--)
	if (this->DiagonalElements[j] < MinValue)
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
// and apply corresponding transformation to column of a given real matrix 
//
// matrix = matrix on which transformation has to be applied
// return value = reference on current Matrix

RealDiagonalMatrix& RealDiagonalMatrix::SortMatrixUpOrder(RealMatrix& matrix)
{
  int ReducedDim = this->NbrColumn - 2;
  RealVector TmpV;
  double tmp;
  int MinPos;
  double MinValue;
  for (int i = 0; i <= ReducedDim; i++)
    {
      MinPos = this->NbrColumn - 1;
      MinValue = this->DiagonalElements[MinPos];
      for (int j = ReducedDim; j >= i; j--)
	if (this->DiagonalElements[j] < MinValue)
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
// and apply corresponding transformation to column of a given complex matrix 
//
// matrix = matrix on which transformation has to be applied
// return value = reference on current Matrix

RealDiagonalMatrix& RealDiagonalMatrix::SortMatrixUpOrder(ComplexMatrix& matrix)
{
  int ReducedDim = this->NbrColumn - 2;
  ComplexVector TmpV;
  double tmp;
  int MinPos;
  double MinValue;
  for (int i = 0; i <= ReducedDim; i++)
    {
      MinPos = this->NbrColumn - 1;
      MinValue = this->DiagonalElements[MinPos];
      for (int j = ReducedDim; j >= i; j--)
	if (this->DiagonalElements[j] < MinValue)
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

ostream& operator << (ostream& Str, const RealDiagonalMatrix& P)
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      for (int j = 0; j < i; j ++)
	{
	  Str << "0    ";
	}
      Str << P.DiagonalElements[i] << "    ";
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

MathematicaOutput& operator << (MathematicaOutput& Str, const RealDiagonalMatrix& P)
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
