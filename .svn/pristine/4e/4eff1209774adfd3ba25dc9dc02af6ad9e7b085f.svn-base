////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of real antisymmetric matrix                    //
//                                                                            //
//                        last modification : 03/04/2001                      //
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


#include "Matrix/RealAntisymmetricMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "GeneralTools/ListIterator.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include <stdlib.h>


using std::endl;


// default constructor
//

RealAntisymmetricMatrix::RealAntisymmetricMatrix() 
{
  this->Dummy = 0.0;
  this->OffDiagonalElements = 0;
  this->OffDiagonalGarbageFlag =  0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = 0;
  this->MatrixType = Matrix::RealElements | Matrix::Antisymmetric;
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

RealAntisymmetricMatrix::RealAntisymmetricMatrix(int dimension, bool zero) 
{
  this->Dummy = 0.0;
  this->OffDiagonalGarbageFlag =  new int;
  *(this->OffDiagonalGarbageFlag) = 1;
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::RealElements | Matrix::Antisymmetric;
  if (this->NbrRow > 1)
    this->OffDiagonalElements = new double [(this->NbrRow * (this->NbrRow - 1)) / 2];
  else    
    this->OffDiagonalElements = new double [1];
  if (zero == true)
    {
      int pos = 0;
      for (int i = 0; i < this->NbrRow; i++)
	{
	  for (int j = i + 1; j < this->NbrRow; j++)
	    {
	      this->OffDiagonalElements[pos] = 0.0;
	      pos++;
	    }
	}
    }
}

// constructor from matrix elements (without duplicating datas)
//
// upperDiagonal = pointer to upper-diagonal element array (with real part in even position and imaginary part in odd position)
// dimension = matrix dimension

RealAntisymmetricMatrix::RealAntisymmetricMatrix(double* upperDiagonal, int dimension) 
{
  this->Dummy = 0.0;
  this->OffDiagonalElements = upperDiagonal;
  this->OffDiagonalGarbageFlag =  new int;
  *(this->OffDiagonalGarbageFlag) = 1;
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::RealElements | Matrix::Antisymmetric;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

RealAntisymmetricMatrix::RealAntisymmetricMatrix(const RealAntisymmetricMatrix& M) 
{
  this->Dummy = 0.0;
  this->OffDiagonalElements = M.OffDiagonalElements;  
  this->OffDiagonalGarbageFlag = M.OffDiagonalGarbageFlag;
  if (this->OffDiagonalGarbageFlag != 0)
    (*(this->OffDiagonalGarbageFlag))++;  
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::RealElements | Matrix::Antisymmetric;
}

// destructor
//

RealAntisymmetricMatrix::~RealAntisymmetricMatrix() 
{
  if (this->OffDiagonalGarbageFlag != 0)
    {
      if ((*(this->OffDiagonalGarbageFlag)) == 1)
	{
	  delete this->OffDiagonalGarbageFlag;
	  delete[] this->OffDiagonalElements;
	}
      else
	(*(this->OffDiagonalGarbageFlag))--;
    }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

RealAntisymmetricMatrix& RealAntisymmetricMatrix::operator = (const RealAntisymmetricMatrix& M) 
{
  if (this->OffDiagonalGarbageFlag != 0)
    {
      if ((*(this->OffDiagonalGarbageFlag)) == 1)
	{
	  delete this->OffDiagonalGarbageFlag;
	  delete[] this->OffDiagonalElements;
	}
      else
	(*(this->OffDiagonalGarbageFlag))--;
    }
  this->OffDiagonalElements = M.OffDiagonalElements;
  this->OffDiagonalGarbageFlag = M.OffDiagonalGarbageFlag;
  if (this->OffDiagonalGarbageFlag != 0)
    (*(this->OffDiagonalGarbageFlag))++;  
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

Matrix* RealAntisymmetricMatrix::Clone ()
{
  return ((Matrix*) new RealAntisymmetricMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void RealAntisymmetricMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn) || (i == j))
    return;
  if (i > j)
    {
      i -= (j * (j + 1)) / 2 - j * (this->NbrRow + this->Increment - 1) + 1;
      this->OffDiagonalElements[i] = -x;
    }
  else
    {
      j -= (i * (i + 1)) / 2 - i * (this->NbrRow + this->Increment - 1) + 1;
      this->OffDiagonalElements[j] = x;
    }
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element
void RealAntisymmetricMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  return;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void RealAntisymmetricMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn) || (i == j))
    return;
  if (i > j)
    {
      i -= (j * (j + 1)) / 2 - j * (this->NbrRow + this->Increment - 1) + 1;
      this->OffDiagonalElements[i] -= x;
    }
  else
    {
      j -= (i * (i + 1)) / 2 - i * (this->NbrRow + this->Increment - 1) + 1;
      this->OffDiagonalElements[j] += x;
    }
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void RealAntisymmetricMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  return;
}

// get reference of a given matrix element supposing i < j
//
// i = line position
// j = column position
// return value = reference om matrix elememt

double& RealAntisymmetricMatrix::operator () (int i, int j)
{
  if (i >= j)
    {
      return this->Dummy;
    }
  j -= (i * (i + 1)) / 2 - i * (this->NbrRow + this->Increment - 1) + 1;
  return this->OffDiagonalElements[j];
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void RealAntisymmetricMatrix::Resize (int nbrRow, int nbrColumn)
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
  double* TmpOffDiag = new double [Tot];
  int k = 0;
  int l = 0;
  for (int i = 0; i < (this->NbrRow - 1); i++)
    {
      for (int j = i + 1; j < this->NbrRow; j++)
	TmpOffDiag[k++] = this->OffDiagonalElements[l++];
      l += this->Increment;
      for (int j = this->NbrRow; j < nbrRow; j++)
	{
	  TmpOffDiag[k++] = 0.0;
	}      
    }
  for (int i = this->NbrRow * (this->NbrRow - 1); i < Tot; i++)
    TmpOffDiag[i] = 0.0;
  if (this->OffDiagonalGarbageFlag != 0)
    {
      if ((*(this->OffDiagonalGarbageFlag)) == 1)
	{
	  delete this->OffDiagonalGarbageFlag;
	  delete[] this->OffDiagonalElements;
	}
      else
	(*(this->OffDiagonalGarbageFlag))--;
    }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = this->TrueNbrRow - this->NbrRow;
  this->OffDiagonalElements = TmpOffDiag;
  this->OffDiagonalGarbageFlag =  new int;
  *(this->OffDiagonalGarbageFlag) = 1;
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void RealAntisymmetricMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
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
		  this->OffDiagonalElements[k++] = 0.0;
		}
	      k += (this->NbrRow - 2 - i);
	    }
	  for (int i = this->NbrRow * (this->NbrRow - 1); i < Tot; i++)
	    this->OffDiagonalElements[i] = 0.0;
	}
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      this->Increment = (this->TrueNbrRow - this->NbrRow);
      return;
    }
  int Tot = (nbrRow * (nbrRow - 1)) / 2;
  double* TmpOffDiag = new double [Tot];
  int k = 0;
  int l = 0;
  for (int i = 0; i < (this->NbrRow - 1); i++)
    {
      for (int j = i + 1; j < this->NbrRow; j++)
	TmpOffDiag[k++] = this->OffDiagonalElements[l++];
      l += this->Increment;
      for (int j = this->NbrRow; j < nbrRow; j++)
	{
	  TmpOffDiag[k++] = 0.0;
	}      
    }
  for (int i = (this->NbrRow * (this->NbrRow - 1)) / 2; i < Tot; i++)
    TmpOffDiag[i] = 0.0;
  if (this->OffDiagonalGarbageFlag != 0)
    {
      if ((*(this->OffDiagonalGarbageFlag)) == 1)
	{
	  delete this->OffDiagonalGarbageFlag;
	  delete[] this->OffDiagonalElements;
	}
      else
	(*(this->OffDiagonalGarbageFlag))--;
    }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = 2 * (this->TrueNbrRow - this->NbrRow);
  this->OffDiagonalElements = TmpOffDiag;
  this->OffDiagonalGarbageFlag =  new int;
  *(this->OffDiagonalGarbageFlag) = 1;
}

// project matrix into a given subspace
//
// subspace = reference on subspace structure
// return value = pointer to projected matrix

Matrix* RealAntisymmetricMatrix::Project (SubspaceSpaceConverter& subspace)
{
  RealAntisymmetricMatrix* TmpM = new RealAntisymmetricMatrix (subspace.SubspaceDimension);
  int Pos;
  int Pos2;
  int Pos3 = 0;
  for (int i = 0; i < subspace.SubspaceDimension; i++)
    {
      Pos = subspace.SubspaceSpaceConverterArray[i];
      for (int j = i + 1; j < subspace.SubspaceDimension; j++)
	{
	  Pos2 = subspace.SubspaceSpaceConverterArray[j];
	  if (Pos2 > Pos)
	    TmpM->OffDiagonalElements[Pos3++] = this->OffDiagonalElements[Pos2 - 1 - (Pos * (Pos + 3)) / 2 
									 + Pos * this->NbrRow];	
	  else
	    TmpM->OffDiagonalElements[Pos3++] = -this->OffDiagonalElements[Pos - 1 - (Pos2 * (Pos2 + 3)) / 2
									  + Pos2 * this->NbrRow];	
	}
    }
  return TmpM;
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

RealAntisymmetricMatrix operator + (const RealAntisymmetricMatrix& M1, const RealAntisymmetricMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return RealAntisymmetricMatrix();
  int ReducedNbr = M1.NbrRow - 1;
  double* OffDiagonal = new double [M1.NbrRow * ReducedNbr];
  int k = 0;
  int l1 = 0;
  int l2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = 0; j < i; j++)
	OffDiagonal[k++] = M1.OffDiagonalElements[l1++] + M2.OffDiagonalElements[l2++];      
      l1 += M2.Increment;
      l2 += M2.Increment;
    }
  return RealAntisymmetricMatrix(OffDiagonal, M1.NbrRow);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

RealAntisymmetricMatrix operator - (const RealAntisymmetricMatrix& M1, const RealAntisymmetricMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return RealAntisymmetricMatrix();
  int ReducedNbr = M1.NbrRow - 1;
  double* OffDiagonal = new double [M1.NbrRow * ReducedNbr];
  int k = 0;
  int l1 = 0;
  int l2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = 0; j < i; j++)
	OffDiagonal[k++] = M1.OffDiagonalElements[l1++] - M2.OffDiagonalElements[l2++];      
      l1 += M2.Increment;
      l2 += M2.Increment;
    }
  return RealAntisymmetricMatrix(OffDiagonal, M1.NbrRow);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealAntisymmetricMatrix operator * (const RealAntisymmetricMatrix& M, double x) 
{
  int ReducedNbr = M.NbrRow - 1;
  double* OffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  int k = 0;
  int k2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	OffDiagonal[k++] = M.OffDiagonalElements[k2++] * x;
      k2 += M.Increment;
    }
  return RealAntisymmetricMatrix(OffDiagonal, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealAntisymmetricMatrix operator * (double x, const RealAntisymmetricMatrix& M) 
{
  int ReducedNbr = M.NbrRow - 1;
  double* OffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  int k = 0;
  int k2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	OffDiagonal[k++] = M.OffDiagonalElements[k2++] * x;
      k2 += M.Increment;
    }
  return RealAntisymmetricMatrix(OffDiagonal, M.NbrRow);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

RealAntisymmetricMatrix operator / (const RealAntisymmetricMatrix& M, double x) 
{
  x = 1.0 / x;
  int ReducedNbr = M.NbrRow - 1;
  double* OffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  int k = 0;
  int k2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	OffDiagonal[k++] = M.OffDiagonalElements[k2++] * x;
      k2 += M.Increment;
    }
  return RealAntisymmetricMatrix(OffDiagonal, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

RealAntisymmetricMatrix& RealAntisymmetricMatrix::operator += (const RealAntisymmetricMatrix& M) 
{
  if (this->NbrRow < M.NbrRow)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  int k = 0;
  int k2 = 0;  
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	this->OffDiagonalElements[k++] += M.OffDiagonalElements[k2++];
      k += this->Increment;
      k2 += M.Increment;
    }
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

RealAntisymmetricMatrix& RealAntisymmetricMatrix::operator -= (const RealAntisymmetricMatrix& M) 
{
  if (this->NbrRow < M.NbrRow)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  int k = 0;
  int k2 = 0;  
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	this->OffDiagonalElements[k++] -= M.OffDiagonalElements[k2++];
      k += this->Increment;
      k2 += M.Increment;
    }
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

RealAntisymmetricMatrix& RealAntisymmetricMatrix::operator *= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  int ReducedNbr = this->NbrRow - 1;
  int k = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	this->OffDiagonalElements[k++] *= x;
      k += this->Increment;
    }
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

RealAntisymmetricMatrix& RealAntisymmetricMatrix::operator /= (double x)
{
  if (this->NbrRow == 0)
    return *this;
  x = 1.0 / x;
  int ReducedNbr = this->NbrRow - 1;
  int k = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	this->OffDiagonalElements[k++] *= x;
      k += this->Increment;
    }
  return *this;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

double RealAntisymmetricMatrix::MatrixElement (RealVector& V1, RealVector& V2)
{
  double x = 0.0;
  if ((V1.Dimension != this->NbrRow) || (V2.Dimension != this->NbrColumn))
    return x;
  for (int i = 0; i < this->NbrRow ; i++)
    {
      double x2 = 0.0;
      int l = (i - 1);
      for (int k = 0; k < i; k++)
	{
	  x2 -= this->OffDiagonalElements[l] * V2.Components[k];
	  l += (this->NbrColumn - 2 - k) + this->Increment;
	}
      l++;
      for (int k = i + 1; k < this->NbrColumn; k++)
	{
	  x2 += this->OffDiagonalElements[l++] * V2.Components[k];
	}      
      x += x2 * V1.Components[i];
    }
  return x;
}

// conjugate an hermitian matrix with an unitary matrix (Ut M U)
//
// UnitaryM = unitary matrix to use
// return value = conjugated matrix

Matrix* RealAntisymmetricMatrix::Conjugate(RealMatrix& UnitaryM)
{
  if (UnitaryM.NbrRow != this->NbrColumn)
    return 0;
  int NbrOffDiag = (UnitaryM.NbrColumn * (UnitaryM.NbrColumn - 1)) / 2;
  double* TmpOffDiag = new double [NbrOffDiag];
  int i2 = 0;
  int ReducedNbrColumn = UnitaryM.NbrColumn - 1;
  int Inc = this->NbrColumn - 3 + this->Increment;
  double tmp1;
  int k;
  int l;
  for (int i = 0; i < ReducedNbrColumn; i++)
    {
      for (int m = i + 1; m < UnitaryM.NbrColumn; m++)
	{    
	  TmpOffDiag[i2] = 0.0;
	  for (int j = 0; j < this->NbrColumn; j++)
	    {
	      tmp1 = 0.0;
	      k = 0;
	      l = (j - 1);
	      for (; k < j; k++)
		{
		  tmp1 -= this->OffDiagonalElements[l++] * UnitaryM.Columns[m].Components[k];
		  l += Inc - k;
		}
	      l++;
	      k++;
	      for (; k < this->NbrColumn; k++)
		{
		  tmp1 += this->OffDiagonalElements[l++] * UnitaryM.Columns[m].Components[k];
		}
	      TmpOffDiag[i2] += tmp1 * UnitaryM.Columns[i].Components[j];
	    }
	  i2++;
	}    
    }
  return new RealAntisymmetricMatrix(TmpOffDiag, UnitaryM.NbrColumn);
}

// conjugate a matrix with an unitary block-diagonal matrix (Ut M U)
//
// UnitaryM = unitary matrix to use
// return value = pointer to conjugated matrix

Matrix* RealAntisymmetricMatrix::Conjugate(BlockDiagonalMatrix& UnitaryM)
{
  if (UnitaryM.NbrRow != this->NbrColumn)
    return 0;
  RealAntisymmetricMatrix* Result = new RealAntisymmetricMatrix(UnitaryM.NbrColumn);
  int NbrBlock = UnitaryM.Blocks.GetNbrElement();
  Matrix** Blocks = new Matrix* [NbrBlock];
  Matrix** TmpM;
  ListIterator<Matrix*> IterBlocks (UnitaryM.Blocks);
  int i = 0;
  while ((TmpM = IterBlocks()))
    {
      Blocks[i++] = *TmpM;
    }
  int SrcRowIndex = 0;
  int DestRowIndex = 0;
  int SrcColumnIndex = 0;
  int DestColumnIndex = 0;
  for (int i = 0; i < NbrBlock; i++)
    {
//      SrcColumnIndex = SrcRowIndex;
//      DestColumnIndex = DestRowIndex;
      SrcRowIndex = UnitaryM.BlockRowPosition[i];
      DestRowIndex = UnitaryM.BlockColumnPosition[i];
      SrcColumnIndex =  UnitaryM.BlockRowPosition[i];
      DestColumnIndex = UnitaryM.BlockColumnPosition[i];
      if (Blocks[i]->GetMatrixType() == Matrix::RealElements)
	{
	  this->Conjugate(*((RealMatrix*) Blocks[i]), SrcRowIndex, DestRowIndex, *Result);
	}
//      SrcColumnIndex += Blocks[i]->GetNbrRow();
//      DestColumnIndex += Blocks[i]->GetNbrColumn();
      for (int j = i + 1; j < NbrBlock; j++)
	{
	  SrcColumnIndex = UnitaryM.BlockRowPosition[j];
	  DestColumnIndex = UnitaryM.BlockColumnPosition[j];
	  if ((Blocks[i]->GetMatrixType() == Matrix::RealElements) && 
	      (Blocks[j]->GetMatrixType() == Matrix::RealElements))
	    {
	      this->Conjugate(*((RealMatrix*) Blocks[i]), *((RealMatrix*) Blocks[j]), 
			      SrcRowIndex, SrcColumnIndex, DestRowIndex, DestColumnIndex, 
			      *Result);
	    }
//	  SrcColumnIndex += Blocks[j]->GetNbrRow();
//	  DestColumnIndex += Blocks[j]->GetNbrColumn();
	}
      SrcRowIndex += Blocks[i]->GetNbrRow();
      DestRowIndex += Blocks[i]->GetNbrColumn();
    }
  delete[] Blocks;
  return Result;
}

// conjugate a block of the matrix with an unitary matrix (Ut M U)
//
// UnitaryM = unitary matrix to use
// sourcePosition = index of the row where the block to conjugate starts
// destinationPosition = index of the row where the conjugated block has to be stored
// matrix = matrix where result has to be stored

void RealAntisymmetricMatrix::Conjugate(RealMatrix& UnitaryM, int sourcePosition, int destinationPosition,
					RealAntisymmetricMatrix& matrix)
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
	  matrix.OffDiagonalElements[i2] = 0.0;
	  for (int j = 0; j < UnitaryM.NbrRow; j++)
	    {
	      double tmp1 = 0.0;
	      int k = 0;
	      int l = ((j + sourcePosition) - 1 - (sourcePosition * (sourcePosition + 1)) / 2 +
		       sourcePosition * (this->NbrRow + this->Increment - 1));
	      for (; k < j; k++)
		{
		  tmp1 -= this->OffDiagonalElements[l++] * UnitaryM.Columns[m].Components[k];
		  l += (this->NbrColumn - 3 - k - sourcePosition) + this->Increment;
		}
	      l++;
	      k++;
	      for (; k < UnitaryM.NbrRow; k++)
		{
		  tmp1 += this->OffDiagonalElements[l++] * UnitaryM.Columns[m].Components[k];
		}
	      matrix.OffDiagonalElements[i2] += tmp1 * UnitaryM.Columns[i].Components[j];
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

void RealAntisymmetricMatrix::Conjugate(RealMatrix& UnitaryMl, RealMatrix& UnitaryMr, int sourceRowIndex, 
					int sourceColumnIndex, int destinationRowIndex,
					int destinationColumnIndex, RealAntisymmetricMatrix& matrix)
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
	  matrix.OffDiagonalElements[i2] = 0.0;
	  for (int j = 0; j < UnitaryMl.NbrRow; j++)
	    {
	      double tmp1 = 0.0;
	      int l = (sourceColumnIndex - 1 - 
		       ((j + sourceRowIndex) * ((sourceRowIndex + j) + 1)) / 2 +
		       (sourceRowIndex + j) * (this->NbrRow + this->Increment - 1));
	      for (int k = 0; k < UnitaryMr.NbrRow; k++)
		{
		  tmp1 += this->OffDiagonalElements[l++] * UnitaryMr.Columns[m].Components[k];
		}
	      matrix.OffDiagonalElements[i2] += tmp1 * UnitaryMl.Columns[i].Components[j];
	    }
	  i2++;
	}
    }
  return;
}

// evaluate matrix trace
//
// return value = matrix trace 

double RealAntisymmetricMatrix::Tr () 
{
  return 0.0;
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double RealAntisymmetricMatrix::Det () 
{
  return 1.0;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const RealAntisymmetricMatrix& P)
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      int pos = (i - 1);
      for (int j = 0; j < i; j ++)
	{
	  Str << -(P.OffDiagonalElements[pos]) << "    ";
	  pos += (P.NbrRow - j - 2) + P.Increment;
	}
      Str << "0    ";
      pos++;
      for (int j = i + 1; j < P.NbrRow; j++)
	{
	  Str << P.OffDiagonalElements[pos++] << "    ";
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

MathematicaOutput& operator << (MathematicaOutput& Str, const RealAntisymmetricMatrix& P)
{
  Str << "{";
  for (int i = 0; i < P.NbrRow; i++)
    {
      Str << "{";
      int pos = (i - 1);
      for (int j = 0; j < i; j ++)
	{
	  Str << -P.OffDiagonalElements[pos] << ",";
	  pos += (P.NbrRow - j - 2) + P.Increment;
	}
      Str << "0";
      if (i != (P.NbrRow - 1))
	{
	  Str << ",";	  
	  pos++;
	  for (int j = i + 1; j < (P.NbrRow - 1); j++)
	    {
	      Str << P.OffDiagonalElements[pos++] << ",";
	    }
	  Str << P.OffDiagonalElements[pos];
	  Str << "},";
	}
      else
	Str << "}";
    }
  Str << "}";
  return Str;
}

