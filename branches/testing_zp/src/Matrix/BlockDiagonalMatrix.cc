////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of block diagonal matrix                      //
//                                                                            //
//                        last modification : 16/05/2001                      //
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


#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "GeneralTools/ListIterator.h"

#include <stdlib.h>


using std::endl;


// default constructor
//

BlockDiagonalMatrix::BlockDiagonalMatrix() 
{
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::BlockDiagonal;
  this->BlockRowPosition = 0;
  this->BlockColumnPosition = 0;
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

BlockDiagonalMatrix::BlockDiagonalMatrix(int dimension, bool zero) 
{
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::BlockDiagonal;
  this->Blocks += new RealMatrix(this->NbrRow, this->NbrColumn, zero);
  this->BlockRowPosition = new int [1];
  this->BlockColumnPosition = new int [1];
  this->BlockRowPosition[0] = 0;
  this->BlockColumnPosition[0] = 0;
  this->Flag.Initialize();
}

// constructor from matrix elements (without duplicating datas)
//
// blocks = reference on list of pointers to matirx to use as block

BlockDiagonalMatrix::BlockDiagonalMatrix(List<Matrix*>& blocks) 
{
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->MatrixType = 0xffffffff;
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(blocks);
  this->BlockRowPosition = new int [blocks.GetNbrElement()];
  this->BlockColumnPosition = new int [blocks.GetNbrElement()];
  int BlockIndex = 0;
  while ((TmpM = IterMatrix()))
    {
      this->BlockRowPosition[BlockIndex] = this->NbrRow;
      this->BlockColumnPosition[BlockIndex++] = this->NbrColumn;      
      this->NbrRow += (*TmpM)->GetNbrRow();
      this->NbrColumn += (*TmpM)->GetNbrColumn();
      this->MatrixType &= (*TmpM)->GetMatrixType();
      this->Blocks += (*TmpM)->Clone();
    }
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType |= Matrix::BlockDiagonal;
  this->Flag.Initialize();
}

// constructor from matrix elements and block positions (without duplicating datas)
//
// blocks =  reference on list of pointers to matirx to use as block
// rowPosition = array containing block row positions
// columnPosition = array containing block column positions

BlockDiagonalMatrix::BlockDiagonalMatrix(List<Matrix*>& blocks, int* rowPosition, int* columnPosition) 
{
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->MatrixType = 0xffffffff;
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(blocks);
  this->BlockRowPosition = rowPosition;
  this->BlockColumnPosition = columnPosition;
  while ((TmpM = IterMatrix()))
    {
      this->MatrixType &= (*TmpM)->GetMatrixType();
      this->Blocks += (*TmpM)->Clone();
    }
  this->NbrRow += this->BlockRowPosition[this->Blocks.GetNbrElement() - 1] + 
    this->Blocks[this->Blocks.GetNbrElement() - 1]->GetNbrRow();
  this->NbrColumn += this->BlockColumnPosition[this->Blocks.GetNbrElement() - 1] + 
    this->Blocks[this->Blocks.GetNbrElement() - 1]->GetNbrColumn();
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType |= Matrix::BlockDiagonal;
  this->Flag.Initialize();
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

BlockDiagonalMatrix::BlockDiagonalMatrix(const BlockDiagonalMatrix& M) 
{
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = M.MatrixType;
  this->BlockRowPosition = M.BlockRowPosition;
  this->BlockColumnPosition = M.BlockColumnPosition;
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(M.Blocks);
  while ((TmpM = IterMatrix()))
    {
      this->Blocks += (*TmpM)->Clone();
    }
  this->Flag = M.Flag;
}

// destructor
//

BlockDiagonalMatrix::~BlockDiagonalMatrix() 
{
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(this->Blocks);
  while ((TmpM = IterMatrix()))
    {
      delete *TmpM;
    }
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->BlockRowPosition;
      delete[] this->BlockColumnPosition;
    }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

BlockDiagonalMatrix& BlockDiagonalMatrix::operator = (const BlockDiagonalMatrix& M) 
{
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(this->Blocks);
  while ((TmpM = IterMatrix()))
    {
      delete *TmpM;
    }
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->BlockRowPosition;
      delete[] this->BlockColumnPosition;
    }
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = M.MatrixType;
  this->BlockRowPosition = M.BlockRowPosition;
  this->BlockColumnPosition = M.BlockColumnPosition;
  this->Flag = M.Flag;
  IterMatrix.DefineList(M.Blocks);
  while ((TmpM = IterMatrix()))
    {
      this->Blocks += (*TmpM)->Clone();
    }
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* BlockDiagonalMatrix::Clone ()
{
  return ((Matrix*) new BlockDiagonalMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void BlockDiagonalMatrix::SetMatrixElement(int i, int j, double x)
{
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(this->Blocks);
  while ((TmpM = IterMatrix()))
    {
      if ((i < (*TmpM)->GetNbrRow()) && (j < (*TmpM)->GetNbrColumn()))
	{
	  (*TmpM)->SetMatrixElement(i, j, x);
	  return;
	}
      i -= (*TmpM)->GetNbrRow();
      j -= (*TmpM)->GetNbrColumn();
      if ((i < 0) || (j < 0))
	return;
    }
  return;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element
void BlockDiagonalMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(this->Blocks);
  while ((TmpM = IterMatrix()))
    {
      if ((i < (*TmpM)->GetNbrRow()) && (j < (*TmpM)->GetNbrColumn()))
	{
	  (*TmpM)->SetMatrixElement(i, j, x);
	  return;
	}
      i -= (*TmpM)->GetNbrRow();
      j -= (*TmpM)->GetNbrColumn();
      if ((i < 0) || (j < 0))
	return;
    }
  return;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void BlockDiagonalMatrix::AddToMatrixElement(int i, int j, double x)
{
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(this->Blocks);
  while ((TmpM = IterMatrix()))
    {
      if ((i < (*TmpM)->GetNbrRow()) && (j < (*TmpM)->GetNbrColumn()))
	{
	  (*TmpM)->AddToMatrixElement(i, j, x);
	  return;
	}
      i -= (*TmpM)->GetNbrRow();
      j -= (*TmpM)->GetNbrColumn();
      if ((i < 0) || (j < 0))
	return;
    }
  return;
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void BlockDiagonalMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(this->Blocks);
  while ((TmpM = IterMatrix()))
    {
      if ((i < (*TmpM)->GetNbrRow()) && (j < (*TmpM)->GetNbrColumn()))
	{
	  (*TmpM)->AddToMatrixElement(i, j, x);
	  return;
	}
      i -= (*TmpM)->GetNbrRow();
      j -= (*TmpM)->GetNbrColumn();
      if ((i < 0) || (j < 0))
	return;
    }
  return;
}

// get reference of a given matrix element (supposing i == j)
//
// i = line position
// j = column position
// return value = reference om matrix elememt

double& BlockDiagonalMatrix::operator () (int i, int j)
{
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(this->Blocks);
  while ((TmpM = IterMatrix()))
    {
      if ((i < (*TmpM)->GetNbrRow()) && (j < (*TmpM)->GetNbrColumn()))
	{
	  return (**TmpM)(i, j);
	}
      i -= (*TmpM)->GetNbrRow();
      j -= (*TmpM)->GetNbrColumn();
      if ((i < 0) || (j < 0))
	return this->Dummy;
    }
  return this->Dummy;
}

// get pointer to a given block
//
// i = block position
// return value = pointer to the block

Matrix* BlockDiagonalMatrix::operator [] (int i)
{
  return this->Blocks[i];
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void BlockDiagonalMatrix::Resize (int nbrRow, int nbrColumn)
{
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void BlockDiagonalMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
}

// add a block
//
// M = pointer to the matrix to add as a block
// return value = reference on current matrix

BlockDiagonalMatrix& BlockDiagonalMatrix::operator += (Matrix* M)
{
  if (this->Blocks.GetNbrElement() == 0)
    {  
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->MatrixType = 0xffffffff;
    }
  this->Blocks += M->Clone();
  this->NbrRow += M->GetNbrRow();
  this->NbrColumn += M->GetNbrColumn();
  this->MatrixType &= M->GetMatrixType();  
  this->MatrixType |= Matrix::BlockDiagonal; 
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  return *this;
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

BlockDiagonalMatrix operator + (const BlockDiagonalMatrix& M1, const BlockDiagonalMatrix& M2)
{
  return BlockDiagonalMatrix();
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

BlockDiagonalMatrix operator - (const BlockDiagonalMatrix& M1, const BlockDiagonalMatrix& M2)
{
  return BlockDiagonalMatrix();
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

BlockDiagonalMatrix operator * (const BlockDiagonalMatrix& M, double x) 
{
  return BlockDiagonalMatrix();
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

BlockDiagonalMatrix operator * (double x, const BlockDiagonalMatrix& M) 
{
  return BlockDiagonalMatrix();
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

BlockDiagonalMatrix operator / (const BlockDiagonalMatrix& M, double x) 
{
  x = 1.0 / x;
  return BlockDiagonalMatrix();
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

BlockDiagonalMatrix& BlockDiagonalMatrix::operator += (const BlockDiagonalMatrix& M) 
{
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

BlockDiagonalMatrix& BlockDiagonalMatrix::operator -= (const BlockDiagonalMatrix& M) 
{
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

BlockDiagonalMatrix& BlockDiagonalMatrix::operator *= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(this->Blocks);
  while ((TmpM = IterMatrix()))
    {
      //      (**TmpM) *= x;
    }
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

BlockDiagonalMatrix& BlockDiagonalMatrix::operator /= (double x)
{
  if (this->NbrRow == 0)
    return *this;
  x = 1.0 / x;
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(this->Blocks);
  while ((TmpM = IterMatrix()))
    {
      //  (**TmpM) *= x;
    }
  return *this;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

double BlockDiagonalMatrix::MatrixElement (RealVector& V1, RealVector& V2)
{
  double x = 0.0;
  /*  if ((V1.Dimension != this->NbrRow) || (V2.Dimension != this->NbrColumn))
    return x;
  for (int i = 0; i < this->NbrRow ; i++)
    {
      x += V1.Components[i] * this->DiagonalElements[i] * V2.Components[i];
      }*/
  return x;
}

// evaluate matrix trace
//
// return value = matrix trace 

double BlockDiagonalMatrix::Tr () 
{
  if (this->Blocks.GetNbrElement() == 0)
    return 0.0;
  double x = 0.0;
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(this->Blocks);
  while ((TmpM = IterMatrix()))
    {
      x += (*TmpM)->Tr();
    }
  return x;
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double BlockDiagonalMatrix::Det () 
{
  if (this->Blocks.GetNbrElement() == 0)
    return 1.0;
  double x = 1.0;
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(this->Blocks);
  while ((TmpM = IterMatrix()))
    {
      x *= (*TmpM)->Det();
    }
  return x;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const BlockDiagonalMatrix& P)
{
  int BlockID = 0 ;
  Matrix** TmpM;
  ListIterator<Matrix*> IterMatrix(P.Blocks);
  while ((TmpM = IterMatrix()))
    {
      Str << "block " << BlockID++ << ":" << endl;
      Str << **TmpM;
    }
  return Str;
}

#ifdef USE_OUTPUT

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// P = matrix to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, const BlockDiagonalMatrix& P)
{
  /*  Str << "{";
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
    Str << "}";*/
  return Str;
}

#endif
