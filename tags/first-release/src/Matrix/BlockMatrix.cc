////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                             class of block matrix                          //
//                                                                            //
//                        last modification : 28/05/2001                      //
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


#include "Matrix/BlockMatrix.h"
#include "Matrix/RealMatrix.h"
#include "GeneralTools/ListIterator.h"

#include <stdlib.h>


// default constructor
//

BlockMatrix::BlockMatrix() 
{
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->NbrBlockRow = 0;
  this->NbrBlockColumn = 0;
  this->GarbageFlag = 0;
  this->Blocks = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::Block;
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

BlockMatrix::BlockMatrix(int nbrRow, int nbrColumn, bool zero) 
{
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->NbrBlockRow = 1;
  this->NbrBlockColumn = 1;
  this->GarbageFlag = new int;
  (*(this->GarbageFlag)) = 1;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::Block;
  this->Blocks = new Matrix** [1];
  this->Blocks[0] = new Matrix* [1];
  this->Blocks[0][0] = new RealMatrix(this->NbrRow, this->NbrColumn, zero);
}

// constructor from matrix elements (without duplicating datas)
//
// blocks = list of pointers to matirx to use as block

BlockMatrix::BlockMatrix(Matrix* matrix)//, const SpaceDecomposition& rowDecomposition) 
{
  this->NbrRow = matrix->GetNbrRow();
  this->NbrColumn = matrix->GetNbrColumn();
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->GarbageFlag = new int;
  (*(this->GarbageFlag))++;
  this->MatrixType = matrix->GetMatrixType() | Matrix::Block;

}

// constructor from real symmetric matrix
//
// blocks = list of pointers to matirx to use as block
/*
BlockMatrix::BlockMatrix(RealSymmetricMatrix& matrix, const SpaceDecomposition& decomposition) 
{
  this->NbrRow = matrix->GetNbrRow();
  this->NbrColumn = matrix->GetNbrColumn();
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->GarbageFlag = new int;
  (*(this->GarbageFlag))++;
  this->MatrixType = matrix.GetMatrixType() | Matrix::Block;
  this->NbrBlockRow = decomposition.GetNbrSubspace();
  
}
*/

// copy constructor (without duplicating datas)
//
// M = matrix to copy

BlockMatrix::BlockMatrix(const BlockMatrix& M) 
{
  if (M.GarbageFlag == 0)
    {
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->NbrBlockRow = 0;
      this->NbrBlockColumn = 0;
      this->GarbageFlag = 0;
      this->Blocks = 0;
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;
      this->MatrixType = Matrix::Block;
    }
  else
    {
      this->NbrRow = M.NbrRow;
      this->NbrColumn = M.NbrColumn;
      this->NbrBlockRow = M.NbrBlockRow;
      this->NbrBlockColumn = M.NbrBlockColumn;
      this->TrueNbrRow = M.TrueNbrRow;
      this->TrueNbrColumn = M.TrueNbrColumn;
      this->MatrixType = M.MatrixType;
      this->Blocks = M.Blocks;
    }
}

// destructor
//

BlockMatrix::~BlockMatrix() 
{
  if (this->GarbageFlag != 0)
    {
      if ((*(this->GarbageFlag)) == 1)
	{
	  for (int i = 0; i < this->NbrBlockRow; i++)
	    {
	      for (int j = 0; j < this->NbrBlockColumn; j++)
		delete this->Blocks[i][j];
	      delete[] this->Blocks[i];
	    }
	  delete[] this->Blocks;
	  delete this->GarbageFlag;
	}
      else
	(*(this->GarbageFlag))--;
    }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

BlockMatrix& BlockMatrix::operator = (const BlockMatrix& M) 
{
  if (this->GarbageFlag != 0)
    {
      if ((*(this->GarbageFlag)) == 1)
	{
	  for (int i = 0; i < this->NbrBlockRow; i++)
	    {
	      for (int j = 0; j < this->NbrBlockColumn; j++)
		delete this->Blocks[i][j];
	      delete[] this->Blocks[i];
	    }
	  delete[] this->Blocks;
	  delete this->GarbageFlag;
	}
      else
	(*(this->GarbageFlag))--;
    }
  if (M.GarbageFlag == 0)
    {
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->NbrBlockRow = 0;
      this->NbrBlockColumn = 0;
      this->GarbageFlag = 0;
      this->Blocks = 0;
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;
      this->MatrixType = Matrix::Block;
    }
  else
    {
      this->NbrRow = M.NbrRow;
      this->NbrColumn = M.NbrColumn;
      this->NbrBlockRow = M.NbrBlockRow;
      this->NbrBlockColumn = M.NbrBlockColumn;
      this->TrueNbrRow = M.TrueNbrRow;
      this->TrueNbrColumn = M.TrueNbrColumn;
      this->MatrixType = M.MatrixType;
      this->Blocks = M.Blocks;
    }
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* BlockMatrix::Clone ()
{
  return ((Matrix*) new BlockMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void BlockMatrix::SetMatrixElement(int i, int j, double x)
{
  return;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element
void BlockMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  return;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void BlockMatrix::AddToMatrixElement(int i, int j, double x)
{
  return;
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void BlockMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  return;
}

// get reference of a given matrix element (supposing i == j)
//
// i = line position
// j = column position
// return value = reference om matrix elememt

double& BlockMatrix::operator () (int i, int j)
{
  
  return this->Dummy;
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void BlockMatrix::Resize (int nbrRow, int nbrColumn)
{
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void BlockMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

BlockMatrix operator + (const BlockMatrix& M1, const BlockMatrix& M2)
{
  return BlockMatrix();
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

BlockMatrix operator - (const BlockMatrix& M1, const BlockMatrix& M2)
{
  return BlockMatrix();
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

BlockMatrix operator * (const BlockMatrix& M, double x) 
{
  return BlockMatrix();
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

BlockMatrix operator * (double x, const BlockMatrix& M) 
{
  return BlockMatrix();
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

BlockMatrix operator / (const BlockMatrix& M, double x) 
{
  x = 1.0 / x;
  return BlockMatrix();
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

BlockMatrix& BlockMatrix::operator += (const BlockMatrix& M) 
{
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

BlockMatrix& BlockMatrix::operator -= (const BlockMatrix& M) 
{
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

BlockMatrix& BlockMatrix::operator *= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

BlockMatrix& BlockMatrix::operator /= (double x)
{
  if (this->GarbageFlag == 0)
    return *this;
  return *this;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

double BlockMatrix::MatrixElement (RealVector& V1, RealVector& V2)
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

double BlockMatrix::Tr () 
{
  if (this->GarbageFlag == 0) 
    return 0.0;
  double x = 0.0;
  return x;
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double BlockMatrix::Det () 
{
  if (this->GarbageFlag == 0)
    return 1.0;
  return 1.0;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const BlockMatrix& P)
{
  /*  for (int i = 0; i < P.NbrRow; i++)
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
      }*/
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// P = matrix to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, const BlockMatrix& P)
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

