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


// default constructor
//

ComplexMatrix::ComplexMatrix() 
{
  this->Columns = 0;
  this->ColumnGarbageFlag = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;  
  this->MatrixType = Matrix::ComplexElements;
}

// constructor from matrix elements (without duplicating datas)
//
// columns = pointer an array of vector
// nbrColumn = number of columns

ComplexMatrix::ComplexMatrix(ComplexVector* columns, int nbrColumn) 
{
  this->Columns = columns;
  this->ColumnGarbageFlag = new int;
  *(this->ColumnGarbageFlag) = 1;
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
  this->ColumnGarbageFlag = M.ColumnGarbageFlag;
  (*(this->ColumnGarbageFlag))++;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;  
  this->MatrixType = Matrix::ComplexElements;
}

// destructor
//

ComplexMatrix::~ComplexMatrix() 
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

ComplexMatrix& ComplexMatrix::operator = (const ComplexMatrix& M) 
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
  this->Columns = M.Columns;
  this->ColumnGarbageFlag = M.ColumnGarbageFlag;
  (*(this->ColumnGarbageFlag))++;
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
  this->Columns[j].RealComponents[i] = x;
  this->Columns[j].ImaginaryComponents[i] = 0.0;
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
  this->Columns[j].RealComponents[i] = x.Re;
  this->Columns[j].ImaginaryComponents[i] = x.Im;
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
  this->Columns[j].RealComponents[i] += x;
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
  this->Columns[j].RealComponents[i] += x.Re;
  this->Columns[j].ImaginaryComponents[i] += x.Im;
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
      delete[] this->Columns;
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
      delete[] this->Columns;
      this->Columns = Tmp;
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
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] + M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j] + M2.Columns[i].ImaginaryComponents[j];
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
	  TmpColumns[i].RealComponents[j] = M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = M2.Columns[i].ImaginaryComponents[j];
	}
      if (i > 0)
	{
	  TmpColumns[i].RealComponents[j] = M1.UpperDiagonalElements[i - 1] + M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = M2.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      TmpColumns[i].RealComponents[j] = M1.DiagonalElements[i] + M2.Columns[i].RealComponents[j];
      TmpColumns[i].ImaginaryComponents[j] = M2.Columns[i].ImaginaryComponents[j];
      ++j;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].RealComponents[j] = M1.UpperDiagonalElements[i + 1] + M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = M2.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      ++j;
      for (; j < M1.NbrColumn; ++j)
	{
	  TmpColumns[i].RealComponents[j] = M2.Columns[i].RealComponents[j];	
	  TmpColumns[i].ImaginaryComponents[j] = M2.Columns[i].ImaginaryComponents[j];	
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
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
	}
      if (i > 0)
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] + M2.UpperDiagonalElements[i - 1];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] + M2.DiagonalElements[i];
      TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
      ++j;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] + M2.UpperDiagonalElements[i + 1];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      ++j;
      for (; j < M1.NbrColumn; ++j)
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j];	
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];	
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
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] - M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j] - M2.Columns[i].ImaginaryComponents[j];	  
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
	  TmpColumns[i].RealComponents[j] = -M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = -M2.Columns[i].ImaginaryComponents[j];
	}
      if (i > 0)
	{
	  TmpColumns[i].RealComponents[j] = M1.UpperDiagonalElements[i - 1] - M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = -M2.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      TmpColumns[i].RealComponents[j] = M1.DiagonalElements[i] - M2.Columns[i].RealComponents[j];
      TmpColumns[i].ImaginaryComponents[j] = -M2.Columns[i].ImaginaryComponents[j];
      ++j;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].RealComponents[j] = M1.UpperDiagonalElements[i + 1] - M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = -M2.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      ++j;
      for (; j < M1.NbrColumn; ++j)
	{
	  TmpColumns[i].RealComponents[j] = -M2.Columns[i].RealComponents[j];	
	  TmpColumns[i].ImaginaryComponents[j] = -M2.Columns[i].ImaginaryComponents[j];	
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
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
	}
      if (i > 0)
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] - M2.UpperDiagonalElements[i - 1];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] - M2.DiagonalElements[i];
      TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
      ++j;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] - M2.UpperDiagonalElements[i + 1];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      ++j;
      for (; j < M1.NbrColumn; ++j)
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j];	
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];	
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
	  TmpColumns[i].RealComponents[j] = 0.0;
	  TmpColumns[i].ImaginaryComponents[j + 1] = 0.0;
	  for (int k = 0; k < M2.NbrRow; ++k)	
	    {
	      TmpColumns[i].RealComponents[j] += (M1.Columns[k].RealComponents[j] * M2.Columns[i].RealComponents[k] - 
						  M1.Columns[k].ImaginaryComponents[j] * M2.Columns[i].ImaginaryComponents[k]);
	      TmpColumns[i].ImaginaryComponents[j] += (M1.Columns[k].RealComponents[j] * M2.Columns[i].ImaginaryComponents[k] + 
						       M1.Columns[k].ImaginaryComponents[j] * M2.Columns[i].RealComponents[k]);
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
	  TmpColumns[i].RealComponents[j] = M.Columns[i].RealComponents[j] * x;
	  TmpColumns[i].ImaginaryComponents[j] = M.Columns[i].ImaginaryComponents[j] * x;
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
	  TmpColumns[i].RealComponents[j] = M.Columns[i].RealComponents[j] * x;
	  TmpColumns[i].ImaginaryComponents[j] = M.Columns[i].ImaginaryComponents[j] * x;
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
	  TmpColumns[i].RealComponents[j] = M.Columns[i].RealComponents[j] * x;
	  TmpColumns[i].ImaginaryComponents[j] = M.Columns[i].ImaginaryComponents[j] * x;
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
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow) || (this->ColumnGarbageFlag == 0))
    return *this;  
  this->Columns[0].RealComponents[0] += M.DiagonalElements[0];
  for (int i = 1; i < this->NbrColumn; i++)
    {
      this->Columns[i].RealComponents[i] += M.DiagonalElements[i];
      this->Columns[i].RealComponents[i - 1] += M.UpperDiagonalElements[i - 1];
      this->Columns[i - 1].RealComponents[i] += M.UpperDiagonalElements[i - 1];
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
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow) || (this->ColumnGarbageFlag == 0))
    return *this;  
  this->Columns[0].RealComponents[0] -= M.DiagonalElements[0];
  for (int i = 1; i < this->NbrColumn; i++)
    {
      this->Columns[i].RealComponents[i] -= M.DiagonalElements[i];
      this->Columns[i].RealComponents[i - 1] -= M.UpperDiagonalElements[i - 1];
      this->Columns[i - 1].RealComponents[i] -= M.UpperDiagonalElements[i - 1];
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
	  Str << P.Columns[j].RealComponents[i];      
	  if (P.Columns[j].ImaginaryComponents[i] < 0.0)
	    Str << P.Columns[j].ImaginaryComponents[i] << "i    ";
	  else
	    if (P.Columns[j].ImaginaryComponents[i] != 0.0)
	      Str << "+" << P.Columns[j].ImaginaryComponents[i] << "i    ";
	    else
	      Str << "    ";
	}
      Str << P.Columns[P.NbrColumn - 1].RealComponents[i];      
      if (P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] < 0.0)
	Str << P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] << "i";
      else
	if (P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] != 0.0)
	  Str << "+" << P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] << "i";
      Str << endl;
    }
  return Str;
}

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
	  Str << P.Columns[j].RealComponents[i];      
	  if (P.Columns[j].ImaginaryComponents[i] < 0.0)
	    Str << P.Columns[j].ImaginaryComponents[i] << "I,";
	  else
	    if (P.Columns[j].ImaginaryComponents[i] != 0.0)
	      Str << "+" << P.Columns[j].ImaginaryComponents[i] << "I,";
	    else
	      Str << ",";
	}
      Str << P.Columns[P.NbrColumn - 1].RealComponents[i];      
      if (P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] < 0.0)
	Str << P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] << "I";
      else
	if (P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] != 0.0)
	  Str << "+" << P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] << "I";
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (P.NbrColumn - 1); ++j)
    {
      Str << P.Columns[j].RealComponents[P.NbrRow - 1];      
      if (P.Columns[j].ImaginaryComponents[P.NbrRow - 1] < 0.0)
	Str << P.Columns[j].ImaginaryComponents[P.NbrRow - 1] << "I,";
      else
	if (P.Columns[j].ImaginaryComponents[P.NbrRow - 1] != 0.0)
	  Str << "+" << P.Columns[j].ImaginaryComponents[P.NbrRow - 1] << "I,";
	else
	  Str << ",";
    }
  Str << P.Columns[P.NbrColumn - 1].RealComponents[P.NbrRow - 1];      
  if (P.Columns[P.NbrColumn - 1].ImaginaryComponents[P.NbrRow - 1] < 0.0)
    Str << P.Columns[P.NbrColumn - 1].ImaginaryComponents[P.NbrRow - 1] << "I";
  else
    if (P.Columns[P.NbrColumn - 1].ImaginaryComponents[P.NbrRow - 1] != 0.0)
      Str << "+" << P.Columns[P.NbrColumn - 1].ImaginaryComponents[P.NbrRow - 1] << "I";
  Str << "}}";
  return Str;
}
