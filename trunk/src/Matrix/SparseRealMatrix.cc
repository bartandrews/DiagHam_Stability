////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of real matrix with sparse storage                  //
//                                                                            //
//                        last modification : 17/10/2012                      //
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


#include "Matrix/SparseRealMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Vector/RealVector.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"
#include "Architecture/ArchitectureOperation/SparseMatrixMatrixMultiplyOperation.h"

#include <iostream>
#include <cstdlib>


using std::endl;
using std::cout;




// default constructor
//

SparseRealMatrix::SparseRealMatrix() 
{
  this->MatrixElements = 0;
  this->ColumnIndices = 0;
  this->RowPointers = 0;
  this->RowLastPointers = 0;
  this->NbrMatrixElements = 0l;
  this->MaximumNbrMatrixElements = 0l;
  this->NbrMatrixElementPacketSize = 0l;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;  
  this->MatrixType = Matrix::RealElements | Matrix::Sparse;
}

// constructor for a sparse matrix without any specific struture
//
// nbrRow = number of rows
// nbrColumn = number of columns

SparseRealMatrix::SparseRealMatrix(int nbrRow, int nbrColumn)
{
  this->Flag.Initialize();
  this->MatrixType = Matrix::RealElements | Matrix::Sparse;
  this->NbrColumn = nbrColumn;
  this->NbrRow = nbrRow;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->NbrMatrixElements = 0;
  this->RowPointers = new long[this->NbrRow];
  this->RowLastPointers = new long[this->NbrRow];
  this->NbrMatrixElementPacketSize = 1024l;
  this->MaximumNbrMatrixElements = this->NbrMatrixElementPacketSize;
  this->MatrixElements = new double[this->MaximumNbrMatrixElements];
  this->ColumnIndices = new int[this->MaximumNbrMatrixElements];
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->RowPointers[i] = -1l;
      this->RowLastPointers[i] = -1l;
    }
}

// constructor for a sparse matrix without any specific struture but a given number of non-zero matrix elements
//
// nbrRow = number of rows
// nbrColumn = number of columns
// nbrMatrixElements = number of non-zero matrix elements
// zero = true if matrix elements have to be set to zero

SparseRealMatrix::SparseRealMatrix(int nbrRow, int nbrColumn, long nbrMatrixElements, bool zero)
{
  this->Flag.Initialize();
  this->MatrixType = Matrix::RealElements | Matrix::Sparse;
  this->NbrColumn = nbrColumn;
  this->NbrRow = nbrRow;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->NbrMatrixElements = nbrMatrixElements;
  this->MaximumNbrMatrixElements = nbrMatrixElements;
  this->NbrMatrixElementPacketSize = 0l;
  this->RowPointers = new long[this->NbrRow];
  this->RowLastPointers = new long[this->NbrRow];
  if (this->NbrMatrixElements > 0)
    {
      this->MatrixElements = new double[this->NbrMatrixElements];
      this->ColumnIndices = new int[this->NbrMatrixElements];
      for (int i = 0; i < this->NbrRow; i++)
	{
	  this->RowPointers[i] = -1l;
	  this->RowLastPointers[i] = -1l;
	}
      if (zero == true)
	{
	  for (long i = 0l; i < this->NbrMatrixElements; ++i)
	    {
	      this->MatrixElements[i] = 0.0;
	      this->ColumnIndices[i] = -1;
	    }
	}
   }
  else
    {
      this->MatrixElements = 0;
      this->ColumnIndices = 0;
    }
}

// constructor for a sparse matrix knowing how many non-zero elements per row will be required
//
// nbrRow = number of rows
// nbrColumn = number of columns
// nbrElementPerRow = number of non-zero matrix elements per row

SparseRealMatrix::SparseRealMatrix(int nbrRow, int nbrColumn, int* nbrElementPerRow)
{
  this->Flag.Initialize();
  this->MatrixType = Matrix::RealElements | Matrix::Sparse;
  this->NbrColumn = nbrColumn;
  this->NbrRow = nbrRow;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->NbrMatrixElements = 0l;
  for (int i = 0; i < this->NbrRow; ++i)
    this->NbrMatrixElements += (long) nbrElementPerRow[i];
  this->MaximumNbrMatrixElements = this->NbrMatrixElements;
  this->NbrMatrixElementPacketSize = -1l;
  this->RowPointers = new long[this->NbrRow];
  this->RowLastPointers = new long[this->NbrRow];
  if (this->NbrMatrixElements > 0)
    {
      this->MatrixElements = new double[this->NbrMatrixElements];
      this->ColumnIndices = new int[this->NbrMatrixElements];
      long Index = 0l;
      for (int i = 0; i < this->NbrRow; i++)
	{
	  if (nbrElementPerRow[i] != 0)
	    {
	      this->RowPointers[i] = Index;
	      this->RowLastPointers[i] = this->RowPointers[i] - 1l;
	      Index += (long) nbrElementPerRow[i];
	    }
	  else
	    {
	      this->RowPointers[i] = -1l;
	      this->RowLastPointers[i] = -1l;
	    }
	}
   }
  else
    {
      this->MatrixElements = 0;
      this->ColumnIndices = 0;
    }
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

SparseRealMatrix::SparseRealMatrix(const SparseRealMatrix& M) 
{
  this->Flag = M.Flag;
  this->MatrixElements = M.MatrixElements;
  this->ColumnIndices = M.ColumnIndices;
  this->RowPointers = M.RowPointers;
  this->RowLastPointers = M.RowLastPointers;
  this->NbrMatrixElements = M.NbrMatrixElements;
  this->MaximumNbrMatrixElements = M.MaximumNbrMatrixElements;
  this->NbrMatrixElementPacketSize = M.NbrMatrixElementPacketSize;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;  
  this->MatrixType = Matrix::RealElements | Matrix::Sparse;
}

// copy constructor (duplicating all datas)
//
// M = matrix to copy
// accuracy = value below which a matrix element is considered to be zero

SparseRealMatrix::SparseRealMatrix(Matrix& M, double accuracy)
{
  if ((M.GetNbrRow() == 0) || (M.GetNbrColumn() == 0))
    {
      this->MatrixElements = 0;
      this->ColumnIndices = 0;
      this->RowPointers = 0;
      this->RowLastPointers = 0;
      this->NbrMatrixElements = 0l;
      this->MaximumNbrMatrixElements = 0l;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::RealElements | Matrix::Sparse;
    }
  else
    {
      this->Flag.Initialize();
      this->NbrColumn = M.GetNbrColumn();
      this->NbrRow = M.GetNbrRow();
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;
      this->NbrMatrixElements = 0l;
      double Tmp;
      this->RowPointers = new long[this->NbrRow];
      this->RowLastPointers = new long[this->NbrRow];
      for (int i = 0; i < this->NbrRow; i++)
	{
	  long PreviousNbrMatrixElements = this->NbrMatrixElements;
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      M.GetMatrixElement(i, j, Tmp);
	      if (fabs(Tmp) > accuracy)
		++this->NbrMatrixElements;
	    }
	  if (PreviousNbrMatrixElements == this->NbrMatrixElements)
	    {
	      this->RowPointers[i] = -1l;
	      this->RowLastPointers[i] = -1l;
	    }
	  else
	    {
	      this->RowPointers[i] = PreviousNbrMatrixElements;	      
	      this->RowLastPointers[i] = this->NbrMatrixElements - 1l;
	    }
	}
      this->MatrixElements = new double[this->NbrMatrixElements];
      this->ColumnIndices = new int[this->NbrMatrixElements];
      this->NbrMatrixElements = 0l;
      for (int i = 0; i < this->NbrRow; i++)
	{
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      M.GetMatrixElement(i, j, Tmp);
	      if (fabs(Tmp) > accuracy)
		{
		  this->MatrixElements[this->NbrMatrixElements] = Tmp;
		  this->ColumnIndices[this->NbrMatrixElements] = j;
		  ++this->NbrMatrixElements;
		}
	    }
	}
      this->MaximumNbrMatrixElements = this->NbrMatrixElements;
      this->NbrMatrixElementPacketSize = 0l;
      this->MatrixType = Matrix::RealElements | Matrix::Sparse;
    }
}

// destructor
//

SparseRealMatrix::~SparseRealMatrix() 
{
  if ((this->MatrixElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->MatrixElements;
	delete[] this->ColumnIndices;
	delete[] this->RowPointers;
	delete[] this->RowLastPointers;
      }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

SparseRealMatrix& SparseRealMatrix::operator = (const SparseRealMatrix& M) 
{
  if ((this->MatrixElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->MatrixElements;
	delete[] this->ColumnIndices;
	delete[] this->RowPointers;
	delete[] this->RowLastPointers;
      }
  this->Flag = M.Flag;
  this->MatrixElements = M.MatrixElements;
  this->ColumnIndices = M.ColumnIndices;
  this->RowPointers = M.RowPointers;
  this->RowLastPointers = M.RowLastPointers;
  this->NbrMatrixElements = M.NbrMatrixElements;
  this->MaximumNbrMatrixElements = M.MaximumNbrMatrixElements;
  this->NbrMatrixElementPacketSize = M.NbrMatrixElementPacketSize;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;  
  this->MatrixType = Matrix::RealElements | Matrix::Sparse;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* SparseRealMatrix::Clone ()
{
  return ((Matrix*) new SparseRealMatrix (*this));
}

// copy a matrix into another (duplicating data)
//
// matrix = matrix to copy
// return value = reference on current matrix

SparseRealMatrix& SparseRealMatrix::Copy (SparseRealMatrix& matrix)
{
  if ((this->MatrixElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->MatrixElements;
	delete[] this->ColumnIndices;
	delete[] this->RowPointers;
	delete[] this->RowLastPointers;
      }
  this->Flag.Initialize();
  this->NbrMatrixElements = matrix.NbrMatrixElements;
  this->NbrRow = matrix.NbrRow;
  this->NbrColumn = matrix.NbrColumn;
  this->TrueNbrRow = matrix.TrueNbrRow;
  this->TrueNbrColumn = matrix.TrueNbrColumn;  
  this->MaximumNbrMatrixElements = matrix.MaximumNbrMatrixElements;
  this->NbrMatrixElementPacketSize = matrix.NbrMatrixElementPacketSize;
  this->RowPointers = new long[this->NbrRow];
  this->RowLastPointers = new long[this->NbrRow];
  this->MatrixElements = new double[this->NbrMatrixElements];
  this->ColumnIndices = new int[this->NbrMatrixElements];
  for (int j = 0; j < this->NbrRow; ++j)
    {
      this->RowPointers[j] = matrix.RowPointers[j];
      this->RowLastPointers[j] = matrix.RowLastPointers[j];
    }
  for (long j = 0l; j < this->NbrMatrixElements; ++j)
    {
      this->MatrixElements[j] = matrix.MatrixElements[j];
      this->ColumnIndices[j] = matrix.ColumnIndices[j];
    }
  return *this;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void SparseRealMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  if (this->NbrMatrixElementPacketSize == 0l)
    {
      if (this->RowPointers[i] >= 0l)
	{
	  long TmpIndex = this->FindColumnIndexPosition(j, this->RowPointers[i], this->RowLastPointers[i]);
	  if (TmpIndex >= 0l)
	    {
	      this->MatrixElements[TmpIndex] = x;
	      return;
	    }
	}     
      return;
    }
  if (this->NbrMatrixElementPacketSize < 0l)
    {
      long& TmpIndex = this->RowLastPointers[i];
      ++TmpIndex;
      this->ColumnIndices[TmpIndex] = j;
      this->MatrixElements[TmpIndex] = x;
      return;
    }
  if (this->RowPointers[i] >= 0l)
    {
      long TmpIndex = this->FindColumnIndexPosition(j, this->RowPointers[i], this->RowLastPointers[i]);
      if (TmpIndex >= 0l)
	{
	  this->MatrixElements[TmpIndex] = x;
	  return;
	}
      this->IncreaseNbrMatrixElements();
      for (int k = i + 1; k < this->NbrRow; ++k)
	{
	  if (this->RowPointers[k] >= 0l)
	    {
	      ++this->RowPointers[k];
	      ++this->RowLastPointers[k];
	    }
	}
      TmpIndex = this->RowPointers[i];
      while ((TmpIndex <= this->RowLastPointers[i]) && (this->ColumnIndices[TmpIndex] < j))
	++TmpIndex;
      for (long k = this->NbrMatrixElements - 1; k > TmpIndex;  --k)
	{
	  this->MatrixElements[k] = this->MatrixElements[k - 1];
	  this->ColumnIndices[k] = this->ColumnIndices[k - 1];	  
	}
      this->MatrixElements[TmpIndex] = x;
      this->ColumnIndices[TmpIndex] = j;
      ++this->RowLastPointers[i];
      return;
    }
  int TmpIndex1 = i;
  while ((TmpIndex1 >= 0) && (this->RowPointers[TmpIndex1] == -1l))
    --TmpIndex1;
  long TmpIndex = 0l;
  if (TmpIndex1 >= 0)
    TmpIndex = this->RowLastPointers[TmpIndex1] + 1l;

  this->IncreaseNbrMatrixElements();
  for (int k = i + 1; k < this->NbrRow; ++k)
    {
      if (this->RowPointers[k] >= 0l)
	{
	  ++this->RowPointers[k];
	  ++this->RowLastPointers[k];
	}
    }
  for (long k = this->NbrMatrixElements - 1; k > TmpIndex;  --k)
    {
      this->MatrixElements[k] = this->MatrixElements[k - 1];
      this->ColumnIndices[k] = this->ColumnIndices[k - 1];	  
    }
  this->MatrixElements[TmpIndex] = x;
  this->ColumnIndices[TmpIndex] = j;
  this->RowPointers[i] = TmpIndex;
  this->RowLastPointers[i] = TmpIndex;
  return;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void SparseRealMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn) || (this->RowPointers[i] == -1l))
    return;
  long TmpIndex = this->FindColumnIndexPosition(j, this->RowPointers[i], this->RowLastPointers[i]);
  if (TmpIndex == -1l)
    {
      return;
    }
  this->MatrixElements[TmpIndex] += x;
}



// increase the number of matrix elements
//
// nbrElements = number of elements to add

void SparseRealMatrix::IncreaseNbrMatrixElements(long nbrElements)
{
  long TmpNbrMatrixElements = this->NbrMatrixElements + nbrElements;
  if (TmpNbrMatrixElements <= this->MaximumNbrMatrixElements)
    {
      this->NbrMatrixElements = TmpNbrMatrixElements;
      return;
    }
  this->MaximumNbrMatrixElements += this->NbrMatrixElementPacketSize;
  double* TmpMatrixElements = new double[this->MaximumNbrMatrixElements];
  int* TmpColumnIndices = new int[this->MaximumNbrMatrixElements];
  for (long i = 0l; i < this->NbrMatrixElements; ++i)
    {
      TmpMatrixElements[i] = this->MatrixElements[i];
      TmpColumnIndices[i] = this->ColumnIndices[i];
    }
  this->NbrMatrixElements = TmpNbrMatrixElements;
  if (this->Flag.Shared() == false)
    {
      delete[] this->MatrixElements;
      delete[] this->ColumnIndices;
    }
  else
    {
      this->Flag.Initialize();
    }
  this->MatrixElements = TmpMatrixElements;
  this->ColumnIndices = TmpColumnIndices;
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void SparseRealMatrix::Resize (int nbrRow, int nbrColumn)
{
//   if (this->NbrRow != nbrRow)
//     {
//       for (int i = 0; i < this->NbrColumn; i++)
// 	this->Columns[i].Resize(nbrRow);
//       if (this->TrueNbrRow >= nbrRow)
// 	{
// 	  this->NbrRow = nbrRow;
// 	}
//       else
// 	{
// 	  this->NbrRow = nbrRow;
// 	  this->TrueNbrRow = nbrRow;
// 	}
//     }
//   if (this->TrueNbrColumn >= nbrColumn)
//     {
//       for (int i = this->NbrColumn; i < nbrColumn; i++)
// 	this->Columns[i].Resize(nbrRow);
//       this->NbrColumn = nbrColumn;
//     }
//   else
//     {
//       RealVector* Tmp = new RealVector[nbrColumn];
//       for (int i = 0; i < this->NbrColumn; i++)
// 	Tmp[i] = this->Columns[i];      
//       for (int i = this->NbrColumn; i < nbrColumn; i++)
// 	Tmp[i] = RealVector(nbrRow);
//       if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
// 	{
// 	  delete[] this->Columns;
// 	}
//       this->Columns = Tmp;
//       this->TrueNbrColumn = nbrColumn;
//       this->NbrColumn = nbrColumn;
//       this->Flag = GarbageFlag();
//       this->Flag.Initialize();
//     }
//   return;
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void SparseRealMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{

  if ((this->MatrixElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->MatrixElements;
	delete[] this->ColumnIndices;
	delete[] this->RowPointers;
	delete[] this->RowLastPointers;
      }

  this->MatrixElements = 0;
  this->ColumnIndices = 0;
  this->RowPointers = 0;
  this->RowLastPointers = 0;
  this->NbrMatrixElements = 0l;
  this->NbrRow = 0;
  this->NbrColumn = 0;

//   if (this->NbrRow != nbrRow)
//     {
//       for (int i = 0; i < this->NbrColumn; i++)
// 	this->Columns[i].ResizeAndClean(nbrRow);
//       if (this->TrueNbrRow >= nbrRow)
// 	{
// 	  this->NbrRow = nbrRow;
// 	}
//       else
// 	{
// 	  this->NbrRow = nbrRow;
// 	  this->TrueNbrRow = nbrRow;
// 	}
//     }
//   if (this->TrueNbrColumn >= nbrColumn)
//     {
//       for (int i = this->NbrColumn; i < nbrColumn; i++)
// 	this->Columns[i].ResizeAndClean(nbrRow);
//       this->TrueNbrColumn = nbrColumn;
//     }
//   else
//     {
//       RealVector* Tmp = new RealVector[nbrColumn];
//       for (int i = 0; i < this->NbrColumn; i++)
// 	Tmp[i] = this->Columns[i];      
//       for (int i = this->NbrColumn; i < nbrColumn; i++)
// 	Tmp[i] = RealVector(nbrRow, true);
//       if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
// 	{
// 	  delete[] this->Columns;
// 	}
//       this->Columns = Tmp;
//       this->Flag = GarbageFlag();
//       this->Flag.Initialize();
//       this->TrueNbrColumn = nbrColumn;
//       this->NbrColumn = nbrColumn;
//     }
//   return;
}

// Set all entries in matrix to zero
//

void SparseRealMatrix::ClearMatrix ()
{
  for (long i = 0l; i < this->NbrMatrixElements; ++i)
    this->MatrixElements[i] = 0.0;
  return;
}

// set matrix to identity 
//

void SparseRealMatrix::SetToIdentity()
{
  if ((this->MatrixElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->MatrixElements;
	delete[] this->ColumnIndices;
	delete[] this->RowPointers;
	delete[] this->RowLastPointers;
      }
  
  this->NbrMatrixElements = (long) this->NbrRow;
  if (this->NbrRow > this->NbrColumn)
    {
      this->NbrMatrixElements = (long) this->NbrColumn;
    }
  this->MaximumNbrMatrixElements = this->NbrMatrixElements;
  this->NbrMatrixElementPacketSize = 0l;
  this->MatrixElements = new double [this->NbrMatrixElements];
  this->ColumnIndices = new int [this->NbrMatrixElements];
  this->RowPointers = new long[this->NbrRow];
  this->RowLastPointers = new long[this->NbrRow];
  for (long i = 0l; i < this->NbrMatrixElements; ++i)
    {
      this->MatrixElements[i] = 1.0;
      this->ColumnIndices[i] = (int) i;
      this->RowPointers[i] = i;
      this->RowLastPointers[i] = i;
    }
  for (int i = (int) this->NbrMatrixElements; i < this->NbrRow; ++i)
    {
      this->RowPointers[i] = -1l;
      this->RowLastPointers[i] = -1l;
    }
  return;  
}

// check if a sparse matrix is correctly formed
//
// return value = true if the sparse matrix is correctly formed

bool SparseRealMatrix::CheckSanity()
{
  long TmpNbrNonZeroMatrixElement = 0l;
  for (int i = 0; i < this->NbrRow; ++i)
    {
      if (this->RowPointers[i] >= 0l)
	{
	  long MinPos = this->RowPointers[i];
	  long MaxPos = this->RowLastPointers[i];
	  TmpNbrNonZeroMatrixElement += MaxPos - MinPos + 1;
	  for (; MinPos < MaxPos; ++MinPos)
	    if (this->ColumnIndices[MinPos] > this->ColumnIndices[MinPos + 1])
	      return false;
	}
    }
  if (TmpNbrNonZeroMatrixElement == this->NbrMatrixElements)
    return true;
  else
    return false;
}

// add two matrices
//
// matrix1 = first matrix
// matrix2 = second matrix
// return value = sum of the two matrices

SparseRealMatrix operator + (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2)
{
  return SparseRealMatrixLinearCombination(1.0, matrix1, 1.0, matrix2);
}

// difference of two matrices
//
// matrix1 = first matrix
// matrix2 = second matrix
// return value = difference of the two matrices

SparseRealMatrix operator - (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2)
{
  return SparseRealMatrixLinearCombination(1.0, matrix1, -1.0, matrix2);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result
SparseRealMatrix operator * (SparseRealMatrix& M, double x)
{ 
  SparseRealMatrix TmpMatrix;
  TmpMatrix.Copy(M);
  TmpMatrix *= x;
  return  TmpMatrix;
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result
SparseRealMatrix operator * (double x,SparseRealMatrix& M)
{
  return operator * (M,x);
}



// create the linear combination of two matrices
//
// x1 = prefactor of the first matrix
// matrix1 = first matrix
// x2 = prefactor of the second matrix
// matrix2 = second matrix
// return value = linear combination

SparseRealMatrix SparseRealMatrixLinearCombination(const double& x1, const SparseRealMatrix& matrix1, const double& x2, const SparseRealMatrix& matrix2)
{
  if ((matrix1.NbrRow != matrix2.NbrRow) || (matrix1.NbrColumn != matrix2.NbrColumn))
    return SparseRealMatrix();
  long TmpNbrMatrixElements = 0l;
  for (int i = 0; i < matrix1.NbrRow; ++i)
    {
      if (matrix1.RowPointers[i] < 0l)
	{
	  if (matrix2.RowPointers[i] >= 0l)
	    TmpNbrMatrixElements += matrix2.RowLastPointers[i] - matrix2.RowPointers[i] + 1l;	 	  
	}
      else
	{
	  if (matrix2.RowPointers[i] < 0l)
	    TmpNbrMatrixElements += matrix1.RowLastPointers[i] - matrix1.RowPointers[i] + 1l;	 	  
	  else
	    {
	      long MinPos1 = matrix1.RowPointers[i];
	      long MinPos2 = matrix2.RowPointers[i];
	      long MaxPos1 = matrix1.RowLastPointers[i];
	      long MaxPos2 = matrix2.RowLastPointers[i];
	      while (MinPos1 <= MaxPos1)
		{
		  while ((MinPos1 <= MaxPos1) && (matrix1.ColumnIndices[MinPos1] < matrix2.ColumnIndices[MinPos2]))
		    {
		      ++TmpNbrMatrixElements;
		      ++MinPos1;
		    }
		  if (MinPos1 <= MaxPos1)
		    {
		      if (matrix1.ColumnIndices[MinPos1] == matrix2.ColumnIndices[MinPos2])
			{
			  ++TmpNbrMatrixElements;
			  ++MinPos1;
			  ++MinPos2;
			}
		      while ((MinPos2 <= MaxPos2) && (matrix2.ColumnIndices[MinPos2] < matrix1.ColumnIndices[MinPos1]))
			{
			  ++TmpNbrMatrixElements;
			  ++MinPos2;
			}
		      if ((MinPos2 > MaxPos2) && (MinPos1 <= MaxPos1))
			{
			  TmpNbrMatrixElements += MaxPos1 - MinPos1 + 1l;
			  MinPos1 = MaxPos1 + 1l;
			}
		    }
		}
	      if (MinPos2 <= MaxPos2)
		{
		  TmpNbrMatrixElements += MaxPos2 - MinPos2 + 1l;
		}		  
	    }
	}
    }  

  SparseRealMatrix TmpMatrix (matrix1.NbrRow, matrix1.NbrColumn, TmpNbrMatrixElements);
  TmpNbrMatrixElements = 0l;
  for (int i = 0; i < matrix1.NbrRow; ++i)
    {
      if (matrix1.RowPointers[i] < 0l)
	{
	  if (matrix2.RowPointers[i] >= 0l)
	    {
	      TmpMatrix.RowPointers[i] = TmpNbrMatrixElements;
	      long MinPos2 = matrix2.RowPointers[i];
	      long MaxPos2 = matrix2.RowLastPointers[i];
	      for (; MinPos2 <= MaxPos2; ++MinPos2)
		{
		  TmpMatrix.MatrixElements[TmpNbrMatrixElements] = x2 * matrix2.MatrixElements[MinPos2];
		  TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = matrix2.ColumnIndices[MinPos2];
		  ++TmpNbrMatrixElements;
		}
	      TmpMatrix.RowLastPointers[i] = TmpNbrMatrixElements - 1l;
	    }
	  else
	    {
	      TmpMatrix.RowLastPointers[i] = -1l;
	      TmpMatrix.RowPointers[i] = -1l;
	    }
	}
      else
	{
	  if (matrix2.RowPointers[i] < 0l)
	    {
	      TmpMatrix.RowPointers[i] = TmpNbrMatrixElements;
	      long MinPos1 = matrix1.RowPointers[i];
	      long MaxPos1 = matrix1.RowLastPointers[i];
	      for (; MinPos1 <= MaxPos1; ++MinPos1)
		{
		  TmpMatrix.MatrixElements[TmpNbrMatrixElements] = x1 * matrix1.MatrixElements[MinPos1];
		  TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = matrix1.ColumnIndices[MinPos1];
		  ++TmpNbrMatrixElements;
		}
	      TmpMatrix.RowLastPointers[i] = TmpNbrMatrixElements - 1l;
	    }
	  else
	    {
	      long MinPos1 = matrix1.RowPointers[i];
	      long MinPos2 = matrix2.RowPointers[i];
	      long MaxPos1 = matrix1.RowLastPointers[i];
	      long MaxPos2 = matrix2.RowLastPointers[i];
	      TmpMatrix.RowPointers[i] = TmpNbrMatrixElements;
	      while (MinPos1 <= MaxPos1)
		{
		  while ((MinPos1 <= MaxPos1) && (matrix1.ColumnIndices[MinPos1] < matrix2.ColumnIndices[MinPos2]))
		    {
		      TmpMatrix.MatrixElements[TmpNbrMatrixElements] = x1 * matrix1.MatrixElements[MinPos1];
		      TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = matrix1.ColumnIndices[MinPos1];
		      ++TmpNbrMatrixElements;
		      ++MinPos1;
		    }
		  if (MinPos1 <= MaxPos1)
		    {
		      if (matrix1.ColumnIndices[MinPos1] == matrix2.ColumnIndices[MinPos2])
			{
			  TmpMatrix.MatrixElements[TmpNbrMatrixElements] = x1 * matrix1.MatrixElements[MinPos1] + x2 * matrix2.MatrixElements[MinPos2];
			  TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = matrix1.ColumnIndices[MinPos1];
			  ++TmpNbrMatrixElements;
			  ++MinPos1;
			  ++MinPos2;
			}
		      while ((MinPos2 <= MaxPos2) && (matrix2.ColumnIndices[MinPos2] < matrix1.ColumnIndices[MinPos1]))
			{
			  TmpMatrix.MatrixElements[TmpNbrMatrixElements] = x2 * matrix2.MatrixElements[MinPos2];
			  TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = matrix2.ColumnIndices[MinPos2];
			  ++TmpNbrMatrixElements;
			  ++MinPos2;
			}		  
		      if (MinPos2 > MaxPos2)		  
			{
			  while (MinPos1 <= MaxPos1)
			    {
			      TmpMatrix.MatrixElements[TmpNbrMatrixElements] = x1 * matrix1.MatrixElements[MinPos1];
			      TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = matrix1.ColumnIndices[MinPos1];
			      ++TmpNbrMatrixElements;
			      ++MinPos1;
			    }
			}
		    }
		}
	      while (MinPos2 <= MaxPos2)
		{
		  TmpMatrix.MatrixElements[TmpNbrMatrixElements] = x2 * matrix2.MatrixElements[MinPos2];
		  TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = matrix2.ColumnIndices[MinPos2];
		  ++TmpNbrMatrixElements;
		  ++MinPos2;
		}		  
	      TmpMatrix.RowLastPointers[i] = TmpNbrMatrixElements - 1l;
	    }
	}
    }      
  return TmpMatrix;
}

// create the linear combination of several matrices
//
// nbrMatrices = number of matrices that should be added
// prefactors = array of prefactors
// matrices = array of matrices
// return value = linear combination

SparseRealMatrix SparseRealMatrixLinearCombination(int nbrMatrices, double* prefactors, SparseRealMatrix* matrices)
{
  if (nbrMatrices <= 0)
    {
      return SparseRealMatrix();
    }
  if (nbrMatrices == 1)
    {
      SparseRealMatrix TmpMatrix;
      TmpMatrix.Copy(matrices[0]);
      TmpMatrix *= prefactors[0];
      return TmpMatrix;
    }
  int TmpNbrRow = matrices[0].GetNbrRow();
  if (TmpNbrRow == 0)
    {
      return SparseRealMatrix();
    }
  for (int i = 1 ; i < nbrMatrices; ++i)
    {
      if (TmpNbrRow != matrices[i].GetNbrRow())
	{
	  return SparseRealMatrix();
	}
    }
  if (nbrMatrices == 2)
    {
      return SparseRealMatrixLinearCombination(prefactors[0], matrices[0], prefactors[1], matrices[1]);
    }

  int TmpNbrMatrices = nbrMatrices / 2;
  if ((nbrMatrices & 1) != 0)
    ++TmpNbrMatrices;
  SparseRealMatrix* TmpMatrices = new SparseRealMatrix[TmpNbrMatrices];
  double* TmpPrefactors = new double[TmpNbrMatrices]; 
  if ((nbrMatrices & 1) != 0)
    --TmpNbrMatrices;
  for (int i = 0 ; i < TmpNbrMatrices; ++i)
    {
      TmpMatrices[i] = SparseRealMatrixLinearCombination(prefactors[2 * i], matrices[2 * i], prefactors[2 * i + 1], matrices[2 * i + 1]);
      TmpPrefactors[i] = 1.0;
    }
  if ((nbrMatrices & 1) != 0)
    {
      ++TmpNbrMatrices;
      TmpMatrices[TmpNbrMatrices - 1] = matrices[nbrMatrices - 1];
      TmpPrefactors[TmpNbrMatrices - 1] = prefactors[nbrMatrices - 1];
    }
  SparseRealMatrix TmpMatrix = SparseRealMatrixLinearCombination(TmpNbrMatrices, TmpPrefactors, TmpMatrices);
  delete[] TmpMatrices;
  delete[] TmpPrefactors;
  return TmpMatrix;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

SparseRealMatrix& SparseRealMatrix::operator *= (double x) 
{
  for (long i = 0; i < this->NbrMatrixElements; ++i)
    this->MatrixElements[i] *= x;
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

SparseRealMatrix& SparseRealMatrix::operator /= (double x)
{
  x = 1.0 / x;;
  for (long i = 0; i < this->NbrMatrixElements; ++i)
    this->MatrixElements[i] *= x;
  return *this;
}

// multiply a matrix to the right by another matrix
//
// matrix = matrix used as multiplicator
// return value = reference on current matrix

SparseRealMatrix& SparseRealMatrix::Multiply (const SparseRealMatrix& matrix)
{
  if (matrix.NbrRow != this->NbrColumn)
    {
      cout << "error, cannot multiply the two matrices" << endl;
      return *this; 
    }
  long TmpNbrMatrixElements = 0l;
  long PreviousTmpNbrMatrixElements = 0l;
  double* TmpElements = new double [matrix.NbrColumn];
  for (int i = 0; i < matrix.NbrColumn; ++i)
    {
      TmpElements[i] = 0.0;
    }

  for (int i = 0; i < this->NbrRow; ++i)
    {
      long MinPos =  this->RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = this->RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      int TmpIndex = this->ColumnIndices[MinPos];
	      long MinPos2 = matrix.RowPointers[TmpIndex];
	      if (MinPos2 >= 0)
		{
		  double Tmp = this->MatrixElements[MinPos];
		  long MaxPos2 = matrix.RowLastPointers[TmpIndex];
		  for (; MinPos2 <= MaxPos2; ++MinPos2)
		    {
		      ++TmpElements[matrix.ColumnIndices[MinPos2]];
		    }      
		}
	    }	    
	  for (int j = 0; j < matrix.NbrColumn; ++j)
	    if (TmpElements[j] != 0)
	      {
		TmpElements[j] = 0.0;
		++TmpNbrMatrixElements;
	      }	  
	}
    }
  
  double* TmpMatrixElements = new double[TmpNbrMatrixElements];
  int* TmpColumnIndices = new int[TmpNbrMatrixElements];
  TmpNbrMatrixElements = 0l;

  for (int i = 0; i < this->NbrRow; ++i)
    {
      long MinPos =  this->RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = this->RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      int TmpIndex = this->ColumnIndices[MinPos];
	      long MinPos2 = matrix.RowPointers[TmpIndex];
	      if (MinPos2 >= 0)
		{
		  double Tmp = this->MatrixElements[MinPos];
		  long MaxPos2 = matrix.RowLastPointers[TmpIndex];
		  for (; MinPos2 <= MaxPos2; ++MinPos2)
		    {
		      TmpElements[matrix.ColumnIndices[MinPos2]] += Tmp * matrix.MatrixElements[MinPos2];
		    }      
		}
	    }	 
   
	  PreviousTmpNbrMatrixElements = TmpNbrMatrixElements;
	  for (int j = 0; j < matrix.NbrColumn; ++j)
	    if (TmpElements[j] != 0)
	      {
		TmpMatrixElements[TmpNbrMatrixElements] = TmpElements[j];
		TmpColumnIndices[TmpNbrMatrixElements] = j;
		TmpElements[j] = 0.0;
		++TmpNbrMatrixElements;
	      }	  
	  this->RowPointers[i] = PreviousTmpNbrMatrixElements;
	  this->RowLastPointers[i] = TmpNbrMatrixElements - 1;
	}
    }
  delete[] TmpElements;
  delete[] this->MatrixElements;
  delete[] this->ColumnIndices;
  this->NbrMatrixElements = TmpNbrMatrixElements;
  this->MatrixElements = TmpMatrixElements;
  this->ColumnIndices = TmpColumnIndices;
//   this->MatrixElements = new double[this->NbrMatrixElements];
//   this->ColumnIndices = new int[this->NbrMatrixElements];
//   for (long i = 0l; i < this->NbrMatrixElements; ++i)
//     {
//       this->MatrixElements[i] = TmpMatrixElements[i];
//       this->ColumnIndices[i] = TmpColumnIndices[i];
//     }
//   delete[] TmpMatrixElements;
//   delete[] TmpColumnIndices;
  return *this;
}

// multiply two matrices
//
// matrix1 = left matrix
// matrix2 = right matrix
// return value = reference on current matrix

SparseRealMatrix Multiply (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2)
{
  if (matrix2.NbrRow != matrix1.NbrColumn)
    {
      cout << "error, cannot multiply the two matrices" << endl;
      return SparseRealMatrix();      
    }
  double* TmpMatrixElements = new double[matrix1.NbrMatrixElements * matrix2.NbrMatrixElements];
  int* TmpColumnIndices = new int[matrix1.NbrMatrixElements * matrix2.NbrMatrixElements];
  double* TmpElements = new double [matrix2.NbrColumn];
  SparseRealMatrix TmpM = Multiply(matrix1, matrix2, TmpMatrixElements, TmpColumnIndices, TmpElements);
  delete[] TmpMatrixElements;
  delete[] TmpColumnIndices;
  delete[] TmpElements;
  return TmpM;
}

// multiply two matrices, minimizing the amount of temporary storage
//
// matrix1 = left matrix
// matrix2 = right matrix
// return value = reference on current matrix

SparseRealMatrix MemoryEfficientMultiply (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2)
{
  if (matrix2.NbrRow != matrix1.NbrColumn)
    {
      cout << "error, cannot multiply the two matrices" << endl;
      return SparseRealMatrix(); 
    }
  long TmpNbrMatrixElements = 0l;
  long PreviousTmpNbrMatrixElements = 0l;
  double* TmpElements = new double [matrix2.NbrColumn];
  for (int i = 0; i < matrix2.NbrColumn; ++i)
    {
      TmpElements[i] = 0.0;
    }
  SparseRealMatrix TmpMatrix(matrix1.NbrRow, matrix2.NbrColumn, 0l);
  for (int i = 0; i < matrix1.NbrRow; ++i)
    {
      long MinPos =  matrix1.RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = matrix1.RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      int TmpIndex = matrix1.ColumnIndices[MinPos];
	      long MinPos2 = matrix2.RowPointers[TmpIndex];
	      if (MinPos2 >= 0)
		{
		  double Tmp = matrix1.MatrixElements[MinPos];
		  long MaxPos2 = matrix2.RowLastPointers[TmpIndex];
		  for (; MinPos2 <= MaxPos2; ++MinPos2)
		    {
		      TmpElements[matrix2.ColumnIndices[MinPos2]] += Tmp * matrix2.MatrixElements[MinPos2];
		    }      
		}
	    }	 
   
	  PreviousTmpNbrMatrixElements = TmpNbrMatrixElements;
	  for (int j = 0; j < matrix2.NbrColumn; ++j)
	    if (TmpElements[j] != 0.0)
	      {
		TmpElements[j] = 0.0;
		++TmpNbrMatrixElements;
	      }	  
	  if (TmpNbrMatrixElements == PreviousTmpNbrMatrixElements)
	    {
	      TmpMatrix.RowPointers[i] = -1l;
	      TmpMatrix.RowLastPointers[i] = 1l;
	    }
	  else
	    {
	      TmpMatrix.RowPointers[i] = PreviousTmpNbrMatrixElements;
	      TmpMatrix.RowLastPointers[i] = TmpNbrMatrixElements - 1;
	    }
	}
      else
	{
	  TmpMatrix.RowPointers[i] = -1l;
	  TmpMatrix.RowLastPointers[i] = -1;
	}
    }
  TmpMatrix.NbrMatrixElements = TmpNbrMatrixElements;
  TmpMatrix.MatrixElements = new double[TmpNbrMatrixElements];
  TmpMatrix.ColumnIndices = new int[TmpNbrMatrixElements];
  TmpNbrMatrixElements = 0l;
  for (int i = 0; i < matrix1.NbrRow; ++i)
    {
      if (TmpMatrix.RowPointers[i] != -1l)
	{
	  long MinPos =  matrix1.RowPointers[i];
	  long MaxPos = matrix1.RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      int TmpIndex = matrix1.ColumnIndices[MinPos];
	      long MinPos2 = matrix2.RowPointers[TmpIndex];
	      if (MinPos2 >= 0)
		{
		  double Tmp = matrix1.MatrixElements[MinPos];
		  long MaxPos2 = matrix2.RowLastPointers[TmpIndex];
		  for (; MinPos2 <= MaxPos2; ++MinPos2)
		    {
		      TmpElements[matrix2.ColumnIndices[MinPos2]] += Tmp * matrix2.MatrixElements[MinPos2];
		    }      
		}
	    }	 
	  
	  PreviousTmpNbrMatrixElements = TmpNbrMatrixElements;
	  for (int j = 0; j < matrix2.NbrColumn; ++j)
	    if (TmpElements[j] != 0.0)
	      {
		TmpMatrix.MatrixElements[TmpNbrMatrixElements] = TmpElements[j];
		TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = j;
		TmpElements[j] = 0.0;
		++TmpNbrMatrixElements;
	      }
	}
    }
  delete[] TmpElements;
  return TmpMatrix;
}

// multiply two matrices, providing all the required temporary arrays
//
// matrix1 = left matrix
// matrix2 = right matrix
// tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
// tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
// tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
// return value = reference on current matrix

SparseRealMatrix Multiply (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2, 
			   double* tmpMatrixElements, int* tmpColumnIndices, double* tmpElements)
{
  if (matrix2.NbrRow != matrix1.NbrColumn)
    {
      cout << "error, cannot multiply the two matrices" << endl;
      return SparseRealMatrix(); 
    }
  long TmpNbrMatrixElements = 0l;
  long PreviousTmpNbrMatrixElements = 0l;
  for (int i = 0; i < matrix2.NbrColumn; ++i)
    {
      tmpElements[i] = 0.0;
    }
  SparseRealMatrix TmpMatrix(matrix1.NbrRow, matrix2.NbrColumn, 0l);
  for (int i = 0; i < matrix1.NbrRow; ++i)
    {
      long MinPos =  matrix1.RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = matrix1.RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      int TmpIndex = matrix1.ColumnIndices[MinPos];
	      long MinPos2 = matrix2.RowPointers[TmpIndex];
	      if (MinPos2 >= 0)
		{
		  double Tmp = matrix1.MatrixElements[MinPos];
		  long MaxPos2 = matrix2.RowLastPointers[TmpIndex];
		  for (; MinPos2 <= MaxPos2; ++MinPos2)
		    {
		      tmpElements[matrix2.ColumnIndices[MinPos2]] += Tmp * matrix2.MatrixElements[MinPos2];
		    }      
		}
	    }	 
   
	  PreviousTmpNbrMatrixElements = TmpNbrMatrixElements;
	  for (int j = 0; j < matrix2.NbrColumn; ++j)
	    if (tmpElements[j] != 0.0)
	      {
		tmpMatrixElements[TmpNbrMatrixElements] = tmpElements[j];
		tmpColumnIndices[TmpNbrMatrixElements] = j;
		tmpElements[j] = 0.0;
		++TmpNbrMatrixElements;
	      }	  
	  if (TmpNbrMatrixElements == PreviousTmpNbrMatrixElements)
	    {
	      TmpMatrix.RowPointers[i] = -1l;
	      TmpMatrix.RowLastPointers[i] = 1l;
	    }
	  else
	    {
	      TmpMatrix.RowPointers[i] = PreviousTmpNbrMatrixElements;
	      TmpMatrix.RowLastPointers[i] = TmpNbrMatrixElements - 1;
	    }
	}
      else
	{
	  TmpMatrix.RowPointers[i] = -1l;
	  TmpMatrix.RowLastPointers[i] = -1;
	}
    }
  TmpMatrix.NbrMatrixElements = TmpNbrMatrixElements;
  TmpMatrix.MatrixElements = new double[TmpNbrMatrixElements];
  TmpMatrix.ColumnIndices = new int[TmpNbrMatrixElements];
  for (long i = 0l; i < TmpNbrMatrixElements; ++i)
    {
      TmpMatrix.MatrixElements[i] = tmpMatrixElements[i];
      TmpMatrix.ColumnIndices[i] = tmpColumnIndices[i];
    }
  return TmpMatrix;
}

// multiply two matrices, providing all the required temporary arrays and using architecture optimisation
//
// matrix1 = pointer to the left matrix
// matrix2 = pointer to the right matrix
// tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
// tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
// nbrTmpMatrixElements = maximum number of elements available in tmpMatrixElements
// architecture = pointer to the architecture
// return value = reference on current matrix

SparseRealMatrix Multiply (SparseRealMatrix* matrix1, SparseRealMatrix* matrix2, 
			   double* tmpMatrixElements, int* tmpColumnIndices, 
			   long nbrTmpMatrixElements, AbstractArchitecture* architecture)
{
  if (matrix2->NbrRow != matrix1->NbrColumn)
    {
      cout << "error, cannot multiply the two matrices" << endl;
      return SparseRealMatrix(); 
    }
  SparseMatrixMatrixMultiplyOperation Operation (matrix1, matrix2, tmpMatrixElements, tmpColumnIndices, nbrTmpMatrixElements);
  Operation.ApplyOperation(architecture);  
  SparseRealMatrix TmpMatrix(Operation.GetDestinationMatrix());
  return TmpMatrix;
}

// multiply a matrix to the right by another matrix, providing all the required temporary arrays
//
// matrix = matrix used as multiplicator
// tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
// tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
// tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
// return value = reference on current matrix

SparseRealMatrix& SparseRealMatrix::Multiply (const SparseRealMatrix& matrix, 
					      double* tmpMatrixElements, int* tmpColumnIndices, double* tmpElements)
{
  if (matrix.NbrRow != this->NbrColumn)
    {
      cout << "error, cannot multiply the two matrices" << endl;
      return *this; 
    }
  long TmpNbrMatrixElements = 0l;
  long PreviousTmpNbrMatrixElements = 0l;
  for (int i = 0; i < matrix.NbrColumn; ++i)
    {
      tmpElements[i] = 0.0;
    }
  for (int i = 0; i < this->NbrRow; ++i)
    {
      long MinPos =  this->RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = this->RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      int TmpIndex = this->ColumnIndices[MinPos];
	      long MinPos2 = matrix.RowPointers[TmpIndex];
	      if (MinPos2 >= 0)
		{
		  double Tmp = this->MatrixElements[MinPos];
		  long MaxPos2 = matrix.RowLastPointers[TmpIndex];
		  for (; MinPos2 <= MaxPos2; ++MinPos2)
		    {
		      tmpElements[matrix.ColumnIndices[MinPos2]] += Tmp * matrix.MatrixElements[MinPos2];
		    }      
		}
	    }	 
   
	  PreviousTmpNbrMatrixElements = TmpNbrMatrixElements;
	  for (int j = 0; j < matrix.NbrColumn; ++j)
	    if (tmpElements[j] != 0.0)
	      {
		tmpMatrixElements[TmpNbrMatrixElements] = tmpElements[j];
		tmpColumnIndices[TmpNbrMatrixElements] = j;
		tmpElements[j] = 0.0;
		++TmpNbrMatrixElements;
	      }	  
	  this->RowPointers[i] = PreviousTmpNbrMatrixElements;
	  this->RowLastPointers[i] = TmpNbrMatrixElements - 1;
	}
    }
  delete[] this->MatrixElements;
  delete[] this->ColumnIndices;
  this->NbrMatrixElements = TmpNbrMatrixElements;
  this->MatrixElements = new double[this->NbrMatrixElements];
  this->ColumnIndices = new int[this->NbrMatrixElements];
  for (long i = 0l; i < this->NbrMatrixElements; ++i)
    {
      this->MatrixElements[i] = tmpMatrixElements[i];
      this->ColumnIndices[i] = tmpColumnIndices[i];
    }
  return *this;
}


// multiply a matrix to the right by another matrix, providing all the required temporary arrays, extend their capacity if needed
//
// matrix = matrix used as multiplicator
// tmpMatrixElements = reference on the temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
// tmpColumnIndices = reference on the temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
// nbrElements = reference ont the number of elements in tmpMatrixElements and tmpColumnIndices
// tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
// return value = reference on current matrix

SparseRealMatrix& SparseRealMatrix::Multiply (const SparseRealMatrix& matrix, 
					      double*& tmpMatrixElements, int*& tmpColumnIndices, 
					      long& nbrElements, double* tmpElements)
{
  if (matrix.NbrRow != this->NbrColumn)
    {
      cout << "error, cannot multiply the two matrices" << endl;
      return *this; 
    }
  long TmpNbrMatrixElements = 0l;
  long PreviousTmpNbrMatrixElements = 0l;
  for (int i = 0; i < matrix.NbrColumn; ++i)
    {
      tmpElements[i] = 0.0;
    }
  for (int i = 0; i < this->NbrRow; ++i)
    {
      long MinPos =  this->RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = this->RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      int TmpIndex = this->ColumnIndices[MinPos];
	      long MinPos2 = matrix.RowPointers[TmpIndex];
	      if (MinPos2 >= 0)
		{
		  double Tmp = this->MatrixElements[MinPos];
		  long MaxPos2 = matrix.RowLastPointers[TmpIndex];
		  for (; MinPos2 <= MaxPos2; ++MinPos2)
		    {
		      tmpElements[matrix.ColumnIndices[MinPos2]] += Tmp * matrix.MatrixElements[MinPos2];
		    }      
		}
	    }	 
   
	  PreviousTmpNbrMatrixElements = TmpNbrMatrixElements;
	  for (int j = 0; j < matrix.NbrColumn; ++j)
	    if (tmpElements[j] != 0)
	      {
		tmpMatrixElements[TmpNbrMatrixElements] = tmpElements[j];
		tmpColumnIndices[TmpNbrMatrixElements] = j;
		tmpElements[j] = 0.0;
		++TmpNbrMatrixElements;
	      }	  
	  this->RowPointers[i] = PreviousTmpNbrMatrixElements;
	  this->RowLastPointers[i] = TmpNbrMatrixElements - 1;
	}
    }
  delete[] this->MatrixElements;
  delete[] this->ColumnIndices;
  this->NbrMatrixElements = TmpNbrMatrixElements;
  this->MatrixElements = new double[this->NbrMatrixElements];
  this->ColumnIndices = new int[this->NbrMatrixElements];
  for (long i = 0l; i < this->NbrMatrixElements; ++i)
    {
      this->MatrixElements[i] = tmpMatrixElements[i];
      this->ColumnIndices[i] = tmpColumnIndices[i];
    }
  return *this;
}

// compute the tensor product of two sparse matrices (matrix1 x matrix2), and store the result in a sparse matrix
//
// matrix1 = reference on the left matrix
// matrix2 = reference on the right matrix
// return value = tensor product

SparseRealMatrix TensorProduct (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2)
{
  SparseRealMatrix TmpMatrix (matrix1.NbrRow * matrix2.NbrRow, 
			      matrix1.NbrColumn * matrix2.NbrColumn, 
			      matrix1.NbrMatrixElements * matrix2.NbrMatrixElements);
  for (int i = 0; i < TmpMatrix.NbrRow; ++i)
    {
      TmpMatrix.RowPointers[i] = -1;
      TmpMatrix.RowLastPointers[i] = -1;
    }
  long TmpPosition = 0l;
  for (int i = 0; i < matrix1.NbrRow; ++i)
    {
      if (matrix1.RowPointers[i] >= 0l)
	{
	  long MinM1 = matrix1.RowPointers[i];
	  long MaxM1 = matrix1.RowLastPointers[i];
	  for (int j = 0; j < matrix2.NbrRow; ++j)
	    {
	      if (matrix2.RowPointers[j] >= 0l)
		{
		  int RowIndex = i *  matrix2.NbrRow + j;
		  TmpMatrix.RowPointers[RowIndex] = TmpPosition;
		  long MinM2 = matrix2.RowPointers[j];
		  long MaxM2 = matrix2.RowLastPointers[j];
		  for (long k1 = MinM1; k1 <= MaxM1; ++k1)
		    {
		      int Shift =  matrix1.ColumnIndices[k1] * matrix2.NbrColumn;
		      double Tmp =  matrix1.MatrixElements[k1];;
		      for (long k2 = MinM2; k2 <= MaxM2; ++k2)
			{
			  TmpMatrix.MatrixElements[TmpPosition] =  Tmp * matrix2.MatrixElements[k2];
			  TmpMatrix.ColumnIndices[TmpPosition] =  matrix2.ColumnIndices[k2] + Shift;			
			  ++TmpPosition;
			}
		    }
		  TmpMatrix.RowLastPointers[RowIndex] = TmpPosition - 1l;
		}
	    }
	}
    }
  return TmpMatrix;  
}

// conjugate the current sparse matrix (M1^+ A M2), assuming A is symmetric
//
// matrix1 = left matrix used for the conjugation
// matrix2 = left matrix used for the conjugation
// return value = conjugated symmetric matrix

RealSymmetricMatrix SparseRealMatrix::Conjugate (RealMatrix& matrix1, RealMatrix& matrix2)
{
  if ((matrix1.GetNbrRow() != this->NbrRow) || (matrix2.GetNbrRow() != this->NbrColumn) || 
      (matrix1.GetNbrColumn() != matrix2.GetNbrColumn()))
    {
      cout << "error, cannot conjugate the matrices" << endl;
      return RealSymmetricMatrix(); 
    }
  RealSymmetricMatrix TmpMatrix(matrix1.GetNbrColumn(), true);
  for (int i = 0; i < matrix1.GetNbrColumn(); ++i)
    {
      for (int k = 0; k < this->NbrRow; ++k)
	{
	  long MinPos =  this->RowPointers[k];
	  if (MinPos >= 0l)
	    {
	      long MaxPos = this->RowLastPointers[k];
	      for (; MinPos <= MaxPos; ++MinPos)
		{
		  double Tmp = matrix1[i][k] * this->MatrixElements[MinPos];
		  int TmpIndex = this->ColumnIndices[MinPos];
		  for (int j = i; j < matrix2.GetNbrColumn(); ++j)
		    {
		      double Tmp2 = Tmp * matrix2[j][TmpIndex];
		      TmpMatrix.AddToMatrixElement(i, j, Tmp2);
		    }
		}
	    }
	}
    }
  return TmpMatrix;
}

// conjugate a matrix
//
// matrix1 = left matrix
// matrix2 = matrix to conjugate
// matrix3 = right matrix
// return value = reference on conjugated matrix

SparseRealMatrix Conjugate (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2, 			    
			    const SparseRealMatrix& matrix3)
{
  double* TmpMatrixElements = new double[matrix1.NbrMatrixElements * matrix2.NbrMatrixElements * matrix3.NbrMatrixElements];
  int* TmpColumnIndices = new int[matrix1.NbrMatrixElements * matrix2.NbrMatrixElements * matrix3.NbrMatrixElements];
  double* TmpElements = new double [matrix3.NbrColumn];
  SparseRealMatrix TmpMatrix = Conjugate(matrix1, matrix2, matrix3, TmpMatrixElements, TmpColumnIndices, TmpElements);
  delete[] TmpMatrixElements;
  delete[] TmpColumnIndices;
  delete[] TmpElements;
  return TmpMatrix;
}

// multiply three matrices, providing all the required temporary arrays
//
// matrix1 = left matrix
// matrix2 = matrix to conjugate
// matrix3 = right matrix
// tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
// tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
// tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
// return value = reference on current matrix

SparseRealMatrix Conjugate (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2, const SparseRealMatrix& matrix3, 
			    double* tmpMatrixElements, int* tmpColumnIndices, double* tmpElements)
{
  if ((matrix2.NbrRow != matrix1.NbrColumn) || (matrix3.NbrRow != matrix2.NbrColumn))
    {
      cout << "error, cannot conjugate the matrices" << endl;
      return SparseRealMatrix(); 
    }
  long TmpNbrMatrixElements = 0l;
  long PreviousTmpNbrMatrixElements = 0l;
  for (int i = 0; i < matrix3.NbrColumn; ++i)
    {
      tmpElements[i] = 0.0;
    }
  SparseRealMatrix TmpMatrix(matrix1.NbrRow, matrix3.NbrColumn, 0l);
  for (int i = 0; i < matrix1.NbrRow; ++i)
    {
      long MinPos =  matrix1.RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = matrix1.RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      int TmpIndex = matrix1.ColumnIndices[MinPos];
	      long MinPos2 = matrix2.RowPointers[TmpIndex];
	      if (MinPos2 >= 0)
		{
		  double Tmp = matrix1.MatrixElements[MinPos];
		  long MaxPos2 = matrix2.RowLastPointers[TmpIndex];
		  for (; MinPos2 <= MaxPos2; ++MinPos2)
		    {
		      int TmpIndex2 = matrix2.ColumnIndices[MinPos2];
		      long MinPos3 = matrix3.RowPointers[TmpIndex2];
		      if (MinPos3 >= 0)
			{
			  double Tmp2 = Tmp * matrix2.MatrixElements[MinPos2];
			  long MaxPos3 = matrix3.RowLastPointers[TmpIndex2];
			  for (; MinPos3 <= MaxPos3; ++MinPos3)
			    {
			      tmpElements[matrix3.ColumnIndices[MinPos3]] += Tmp2 * matrix3.MatrixElements[MinPos3];
			    }
			}
		    }      
		}
	    }	 
   
	  PreviousTmpNbrMatrixElements = TmpNbrMatrixElements;
	  for (int j = 0; j < matrix3.NbrColumn; ++j)
	    if (tmpElements[j] != 0.0)
	      {
		tmpMatrixElements[TmpNbrMatrixElements] = tmpElements[j];
		tmpColumnIndices[TmpNbrMatrixElements] = j;
		tmpElements[j] = 0.0;
		++TmpNbrMatrixElements;
	      }	  
	  if (TmpNbrMatrixElements == PreviousTmpNbrMatrixElements)
	    {
	      TmpMatrix.RowPointers[i] = -1l;
	      TmpMatrix.RowLastPointers[i] = 1l;
	    }
	  else
	    {
	      TmpMatrix.RowPointers[i] = PreviousTmpNbrMatrixElements;
	      TmpMatrix.RowLastPointers[i] = TmpNbrMatrixElements - 1;
	    }
	}
      else
	{
	  TmpMatrix.RowPointers[i] = -1l;
	  TmpMatrix.RowLastPointers[i] = -1;
	}
    }
  TmpMatrix.NbrMatrixElements = TmpNbrMatrixElements;
  TmpMatrix.MatrixElements = new double[TmpNbrMatrixElements];
  TmpMatrix.ColumnIndices = new int[TmpNbrMatrixElements];
  for (long i = 0l; i < TmpNbrMatrixElements; ++i)
    {
      TmpMatrix.MatrixElements[i] = tmpMatrixElements[i];
      TmpMatrix.ColumnIndices[i] = tmpColumnIndices[i];
    }
  return TmpMatrix;
}

// multiply three matrices, providing all the required temporary arrays
//
// matrix1 = pointer to the left matrix
// matrix2 = pointer to the matrix to conjugate
// matrix3 = pointer to the right matrix
// tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
// tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
// nbrTmpMatrixElements = maximum number of elements available in tmpMatrixElements
// architecture = pointer to the architecture
// return value = reference on current matrix

SparseRealMatrix Conjugate (SparseRealMatrix* matrix1, SparseRealMatrix* matrix2, SparseRealMatrix* matrix3, 
			    double* tmpMatrixElements, int* tmpColumnIndices, 
			    long nbrTmpMatrixElements, AbstractArchitecture* architecture)
{
  if ((matrix2->NbrRow != matrix1->NbrColumn) || (matrix3->NbrRow != matrix2->NbrColumn))
    {
      cout << "error, cannot multiply the two matrices" << endl;
      return SparseRealMatrix(); 
    }
  SparseMatrixMatrixMultiplyOperation Operation (matrix1, matrix2, matrix3,
						 tmpMatrixElements, tmpColumnIndices, nbrTmpMatrixElements);
  Operation.ApplyOperation(architecture);  
  SparseRealMatrix TmpMatrix(Operation.GetDestinationMatrix());
  return TmpMatrix;
}

// matrix-vector multiplication action to the right (i.e. v^t M)
//
// inputVector = vector that will be multiplied
// outputVector = vector where the result will be stored

void SparseRealMatrix::RightMultiply (RealVector& inputVector, RealVector& outputVector)
{
  if ((this->NbrRow != inputVector.GetVectorDimension()) || (outputVector.GetVectorDimension() != this->NbrColumn))
    {
      cout << "error, cannot conjugate the matrices" << endl;
      return; 
    }
  for (int i = 0; i < this->NbrRow; ++i)
    {
      double& Tmp = outputVector[i];
      Tmp = 0.0;
      long MinPos =  this->RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = this->RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      Tmp += this->MatrixElements[MinPos] * inputVector[this->ColumnIndices[MinPos]];
	    }	    
	}
    }
}

// matrix-vector multiplication action to the right including a global scaling factor (i.e. alpha v^t M)
//
// coefficient = global multiplicative coefficient 
// inputVector = vector that will be multiplied
// outputVector = vector where the result will be stored

void SparseRealMatrix::RightMultiply (double coefficient, RealVector& inputVector, RealVector& outputVector)
{
  if ((this->NbrRow != inputVector.GetVectorDimension()) || (outputVector.GetVectorDimension() != this->NbrColumn))
    {
      cout << "error, cannot conjugate the matrices" << endl;
      return; 
    }
  for (int i = 0; i < this->NbrRow; ++i)
    {
      double& Tmp = outputVector[i];
      Tmp = 0.0;
      long MinPos =  this->RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = this->RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      Tmp += this->MatrixElements[MinPos] * inputVector[this->ColumnIndices[MinPos]];
	    }	    
	}
      Tmp *= coefficient;
    }
}

// matrix-vector multiplication action to the right (i.e. v^t M), adding the result to another vector
//
// inputVector = vector that will be multiplied
// outputVector = vector where the result will be added

void SparseRealMatrix::RightAddMultiply (RealVector& inputVector, RealVector& outputVector)
{
  if ((this->NbrRow != inputVector.GetVectorDimension()) || (outputVector.GetVectorDimension() != this->NbrColumn))
    {
      cout << "error, cannot conjugate the matrices" << endl;
      return; 
    }
  for (int i = 0; i < this->NbrRow; ++i)
    {
      double& Tmp = outputVector[i];
      long MinPos =  this->RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = this->RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      Tmp += this->MatrixElements[MinPos] * inputVector[this->ColumnIndices[MinPos]];
	    }	    
	}
    }
}

// matrix-vector multiplication action to the right including a global scaling factor (i.e. alpha v^t M), adding the result to another vector
//
// coefficient = global multiplicative coefficient 
// inputVector = vector that will be multiplied
// outputVector = vector where the result will be added

void SparseRealMatrix::RightAddMultiply (double coefficient, RealVector& inputVector, RealVector& outputVector)
{
  if ((this->NbrRow != inputVector.GetVectorDimension()) || (outputVector.GetVectorDimension() != this->NbrColumn))
    {
      cout << "error, cannot conjugate the matrices" << endl;
      return; 
    }
  for (int i = 0; i < this->NbrRow; ++i)
    {
      double& Tmp = outputVector[i];
      long MinPos =  this->RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = this->RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      Tmp += coefficient * this->MatrixElements[MinPos] * inputVector[this->ColumnIndices[MinPos]];
	    }	    
	}
    }
}


// compute the transpose of the current matrix
//
// return value = hermitian transposed matrix

SparseRealMatrix SparseRealMatrix::Transpose ()
{
  if (this->NbrMatrixElements == 0)
    return SparseRealMatrix(this->NbrColumn, this->NbrRow);
  SparseRealMatrix TmpMatrix (this->NbrColumn, this->NbrRow, this->NbrMatrixElements);
  for (int i = 0; i < TmpMatrix.NbrRow; ++i)
    {
      TmpMatrix.RowPointers[i] = 0l;
    }
  for (long i = 0l; i < this->NbrMatrixElements; ++i)
    ++TmpMatrix.RowPointers[this->ColumnIndices[i]];
  int PreviousIndex = 0;
  while ((PreviousIndex < TmpMatrix.NbrRow) && (TmpMatrix.RowPointers[PreviousIndex] == 0l))
    {
      TmpMatrix.RowPointers[PreviousIndex] = -1l;
      TmpMatrix.RowLastPointers[PreviousIndex] = -1l;
      ++PreviousIndex;        
    }
  int Index = PreviousIndex + 1;
  TmpMatrix.RowLastPointers[PreviousIndex] = TmpMatrix.RowPointers[PreviousIndex] - 1l;
  TmpMatrix.RowPointers[PreviousIndex] = 0l;
  while (Index < TmpMatrix.NbrRow)
    {
      while ((Index < TmpMatrix.NbrRow) && (TmpMatrix.RowPointers[Index] == 0l))
	{
	  TmpMatrix.RowPointers[Index] = -1l;
	  TmpMatrix.RowLastPointers[Index] = -1l;
	  ++Index;
	}
      if (Index < TmpMatrix.NbrRow)
	{
	  TmpMatrix.RowLastPointers[Index] = TmpMatrix.RowPointers[Index] + TmpMatrix.RowLastPointers[PreviousIndex];;
	  TmpMatrix.RowPointers[Index] = TmpMatrix.RowLastPointers[PreviousIndex] + 1l;
	  PreviousIndex = Index;	  
	  ++Index;
	}
    }
  for (int i = 0; i < this->NbrRow; ++i)
    {
      long MinPos =  this->RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = this->RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      long& TmpPointer = TmpMatrix.RowPointers[this->ColumnIndices[MinPos]];
	      TmpMatrix.MatrixElements[TmpPointer] = this->MatrixElements[MinPos];
	      TmpMatrix.ColumnIndices[TmpPointer] = i;
	      ++TmpPointer;
	    }
	}
    }
  Index = 0;
  while ((Index < TmpMatrix.NbrRow) && (TmpMatrix.RowPointers[Index] < 0l))
    {
      ++Index;        
    }
  long PreviousPointer =  TmpMatrix.RowPointers[Index];
  TmpMatrix.RowPointers[Index] = 0l;
  ++Index;
   while (Index < TmpMatrix.NbrRow)
    {
      while ((Index < TmpMatrix.NbrRow) && (TmpMatrix.RowPointers[Index] < 0l))
	{
	  ++Index;
	}
      if (Index < TmpMatrix.NbrRow)
	{
	  long Tmp = TmpMatrix.RowPointers[Index];
	  TmpMatrix.RowPointers[Index] = PreviousPointer;
	  PreviousPointer = Tmp;
	  ++Index;
	}
    }
	
  return TmpMatrix;
}

// create a block diagonal matrix from two matrices 
//
// matrix1 = first matrix (i.e. the one at starting from the first row, first column)
// matrix2 = second matrix
// coefficient1 = optional multiplicative coefficient in front of matrix1
// coefficient2 = optional multiplicative coefficient in front of matrix2
// return value = sparse block diagonal matrix

SparseRealMatrix CreateBlockDiagonalMatrix(const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2, double coefficient1, double coefficient2)
{
  SparseRealMatrix TmpMatrix(matrix1.NbrRow + matrix2.NbrRow, matrix1.NbrColumn + matrix2.NbrColumn, matrix1.NbrMatrixElements + matrix2.NbrMatrixElements);
  for (long i = 0; i < matrix1.NbrMatrixElements; ++i)
    TmpMatrix.MatrixElements[i] = coefficient1 * matrix1.MatrixElements[i];
  for (long i = 0; i < matrix2.NbrMatrixElements; ++i)
    TmpMatrix.MatrixElements[i + matrix1.NbrMatrixElements] = coefficient2 * matrix2.MatrixElements[i];
  for (long i = 0; i < matrix1.NbrMatrixElements; ++i)
    TmpMatrix.ColumnIndices[i] = matrix1.ColumnIndices[i];
  for (long i = 0; i < matrix2.NbrMatrixElements; ++i)
    TmpMatrix.ColumnIndices[i + matrix1.NbrMatrixElements] = matrix1.NbrColumn + matrix2.ColumnIndices[i];
  for (int i = 0; i < matrix1.NbrRow; ++i)
    {
      TmpMatrix.RowPointers[i] = matrix1.RowPointers[i];
      TmpMatrix.RowLastPointers[i] = matrix1.RowLastPointers[i];
    }
  for (int i = 0; i < matrix2.NbrRow; ++i)
    {
      if (matrix2.RowPointers[i] >= 0l)
	{
	  TmpMatrix.RowPointers[matrix1.NbrRow + i] = matrix1.NbrMatrixElements + matrix2.RowPointers[i];
	  TmpMatrix.RowLastPointers[matrix1.NbrRow + i] = matrix1.NbrMatrixElements + matrix2.RowLastPointers[i];
	}
      else
	{
	  TmpMatrix.RowPointers[matrix1.NbrRow + i] = -1l;
	  TmpMatrix.RowLastPointers[matrix1.NbrRow + i] = -1l;
	}
    }
  return  TmpMatrix;
}
  
// create an  block off-diagonal matrix from two matrices 
//
// matrix1 = first matrix (i.e. the one at starting from the first row)
// matrix2 = second matrix
// coefficient1 = optional multiplicative coefficient in front of matrix1
// coefficient2 = optional multiplicative coefficient in front of matrix2
// return value = sparse block diagonal matrix

SparseRealMatrix CreateBlockOffDiagonalMatrix(const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2, double coefficient1, double coefficient2)
{
  SparseRealMatrix TmpMatrix(matrix1.NbrRow + matrix2.NbrRow, matrix1.NbrColumn + matrix2.NbrColumn, matrix1.NbrMatrixElements + matrix2.NbrMatrixElements);
  for (long i = 0; i < matrix1.NbrMatrixElements; ++i)
    TmpMatrix.MatrixElements[i] = coefficient1 * matrix1.MatrixElements[i];
  for (long i = 0; i < matrix2.NbrMatrixElements; ++i)
    TmpMatrix.MatrixElements[i + matrix1.NbrMatrixElements] = coefficient2 * matrix2.MatrixElements[i];
  for (long i = 0; i < matrix1.NbrMatrixElements; ++i)
    TmpMatrix.ColumnIndices[i] = matrix2.NbrColumn + matrix1.ColumnIndices[i];
  for (long i = 0; i < matrix2.NbrMatrixElements; ++i)
    TmpMatrix.ColumnIndices[i + matrix1.NbrMatrixElements] = matrix2.ColumnIndices[i];
  for (int i = 0; i < matrix1.NbrRow; ++i)
    {
      TmpMatrix.RowPointers[i] = matrix1.RowPointers[i];
      TmpMatrix.RowLastPointers[i] = matrix1.RowLastPointers[i];
    }
  for (int i = 0; i < matrix2.NbrRow; ++i)
    {
      if (matrix2.RowPointers[i] >= 0l)
	{
	  TmpMatrix.RowPointers[matrix1.NbrRow + i] = matrix1.NbrMatrixElements + matrix2.RowPointers[i];
	  TmpMatrix.RowLastPointers[matrix1.NbrRow + i] = matrix1.NbrMatrixElements + matrix2.RowLastPointers[i];
	}
      else
	{
	  TmpMatrix.RowPointers[matrix1.NbrRow + i] = -1l;
	  TmpMatrix.RowLastPointers[matrix1.NbrRow + i] = -1l;
	}
    }
  return  TmpMatrix;
}
  
// extract a submatrix 
//
// nbrRow = number of rows for the submatrix
// nbrColumn = number of columns for the submatrix
// rowFlags = array that indicated if a row index is part of the submtrix
// columnFlags = array that indicated if a column index is part of the submtrix
// return value = extracted matrix

SparseRealMatrix SparseRealMatrix::ExtractMatrix(int nbrRow, int nbrColumn, bool* rowFlags, bool* columnFlags)
{
  int* TmpNbrElementPerRow = new int[nbrRow];
  int TmpIndex = 0;
  int* TmpColumnIndices = new int [this->NbrColumn];
  for (int i = 0; i < this->NbrColumn; ++i)
    {
      if (columnFlags[i] == true)
	{
	  TmpColumnIndices[i] = TmpIndex;
	  ++TmpIndex;
	}
    }
  TmpIndex = 0;
  long TmpNbrMatrixElements = 0l;
  for (int i = 0; i < this->NbrRow; ++i)
    {
      if (rowFlags[i] == true)
	{
	  int& TmpIndex2 = TmpNbrElementPerRow[TmpIndex];
	  TmpIndex2 = 0;
	  long MinPos = this->RowPointers[i];
	  if (MinPos >=  0l)
	    {
	      long MaxPos = this->RowLastPointers[i];
	      for (; MinPos <= MaxPos; ++MinPos)
		if (columnFlags[this->ColumnIndices[MinPos]] == true)
		  ++TmpIndex2;
	      TmpNbrMatrixElements += TmpIndex2;
	    }
	  ++TmpIndex;
	}
    }
  if (TmpNbrMatrixElements == 0l)
    {
      SparseRealMatrix TmpMatrix;
      return TmpMatrix;
    }
  SparseRealMatrix TmpMatrix(nbrRow, nbrColumn, TmpNbrElementPerRow);
  TmpIndex = 0;
  for (int i = 0; i < this->NbrRow; ++i)
    {
      if (rowFlags[i] == true)
	{
	  long MinPos = this->RowPointers[i];
	  if (MinPos >=  0l)
	    {
	      long MaxPos = this->RowLastPointers[i];
	      for (; MinPos <= MaxPos; ++MinPos)
		if (columnFlags[this->ColumnIndices[MinPos]] == true)
		  {
		    TmpMatrix.SetMatrixElement(TmpIndex, TmpColumnIndices[this->ColumnIndices[MinPos]], this->MatrixElements[MinPos]);
		  }
	    }
	  ++TmpIndex;
	}
    }
  return TmpMatrix;
}
  
  
// extract a submatrix 
//
// nbrRow = number of rows for the submatrix
// nbrColumn = number of columns for the submatrix and onto which index it has to be mapped (negative if it should not be kept)
// rowKeptIndices = array that lists the row indices that have to be kept
// columnKeptIndices = array that lists the column indices that have to be kept
// return value = extracted matrix

SparseRealMatrix SparseRealMatrix::ExtractMatrix(int nbrRow, int nbrColumn, int* rowKeptIndices, int* columnKeptIndices)
{
  int* TmpNbrElementPerRow = new int[nbrRow];
  int* SortedColumnKeptIndices = new int [nbrColumn];
  int* TranslationColumnIndices =  new int [nbrColumn];
  for (int i = 0; i < nbrColumn; ++i)
    {
      TranslationColumnIndices[i] = i;
      SortedColumnKeptIndices[i] = columnKeptIndices[i];
    }
  SortArrayUpOrdering<int>(SortedColumnKeptIndices, TranslationColumnIndices, nbrColumn);
  long TmpNbrMatrixElements = 0l;
  for (int i = 0; i < nbrRow; ++i)
    {
      long MinPos = this->RowPointers[rowKeptIndices[i]];
      int& TmpIndex2 = TmpNbrElementPerRow[i];
      TmpIndex2 = 0;
      if (MinPos >=  0l)
	{
	  long MaxPos = this->RowLastPointers[rowKeptIndices[i]];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      if (SearchInArray<int>(this->ColumnIndices[MinPos], SortedColumnKeptIndices, nbrColumn) >= 0)
		++TmpIndex2;
	    }
	  TmpNbrMatrixElements += TmpIndex2;
	}
    }
  if (TmpNbrMatrixElements == 0l)
    {
      delete[] SortedColumnKeptIndices;
      delete[] TranslationColumnIndices;
      SparseRealMatrix TmpMatrix;
      return TmpMatrix;
    }
  SparseRealMatrix TmpMatrix(nbrRow, nbrColumn, TmpNbrElementPerRow);
  for (int i = 0; i < nbrRow; ++i)
    {
      long MinPos = this->RowPointers[rowKeptIndices[i]];
      if (MinPos >=  0l)
	{
	  long MaxPos = this->RowLastPointers[rowKeptIndices[i]];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      int TmpPos = SearchInArray<int>(this->ColumnIndices[MinPos], SortedColumnKeptIndices, nbrColumn);
	      if (TmpPos >= 0)
		TmpMatrix.SetMatrixElement(i, TranslationColumnIndices[TmpPos], this->MatrixElements[MinPos]);
	    }
	}
    }
  delete[] SortedColumnKeptIndices;
  delete[] TranslationColumnIndices;
  return TmpMatrix;
}
  
// compute the number of non-zero matrix elements (zero having strictly zero square norm)
//
// return value = number of non-zero matrix elements

long SparseRealMatrix::ComputeNbrNonZeroMatrixElements()
{
  long NbrNonZero = 0l;
  for (long i = 0; i < this->NbrMatrixElements; ++i)
    if (this->MatrixElements[i] != 0.0)
      ++NbrNonZero;
  return NbrNonZero;
}

// compute the total amount of memory needed to store the sparse matrix
//
// return value = amount of memory (in bytes)

unsigned long SparseRealMatrix::GetAllocatedMemory()
{
  return (((sizeof(double)+ sizeof(int)) * this->NbrMatrixElements)
	  + ((2ul * sizeof(int)) * this->NbrRow));
}

// evaluate the real part of the matrix trace
//
// return value = real part of the matrix trace 

double SparseRealMatrix::Tr ()
{
  double Trace = 0.0;
  long TmpIndex;
  for (long i = 0; i < this->NbrRow; ++i)
    {
      if (this->RowPointers[i] >= 0l)
	{
	  TmpIndex = this->FindColumnIndexPosition(i, this->RowPointers[i], this->RowLastPointers[i]);
	  if (TmpIndex >= 0l)
	    Trace += this->MatrixElements[TmpIndex];
	}
    }
  return Trace;
}

// evaluate the real part of the matrix partial trace 
//
// indices = array of indices that describes the partial trace 
// nbrIndices = number of indices
// return value = real part of the matrix partial trace 

double SparseRealMatrix::PartialTr(int* indices, int nbrIndices)
{
  double Trace = 0.0;
  long TmpIndex;
  int TmpIndex2;
  for (long i = 0; i < nbrIndices; ++i)
    {
      TmpIndex2 = indices[i];
      if (this->RowPointers[TmpIndex2] >= 0l)
	{
	  TmpIndex = this->FindColumnIndexPosition(TmpIndex2, this->RowPointers[TmpIndex2], this->RowLastPointers[TmpIndex2]);
	  if (TmpIndex >= 0l)
	    Trace += this->MatrixElements[TmpIndex];
	}
    }
  return Trace;
}

// write matrix in a file 
//
// file = reference on the output file stream
// return value = true if no error occurs

bool SparseRealMatrix::WriteMatrix (ofstream& file)
{
  WriteLittleEndian(file, this->MatrixType);
  WriteLittleEndian(file, this->NbrRow);
  WriteLittleEndian(file, this->NbrColumn);
  WriteLittleEndian(file, this->NbrMatrixElements);
  for (int j = 0; j < this->NbrRow; ++j)
    {
      WriteLittleEndian(file, this->RowPointers[j]);
    }
  for (int j = 0; j < this->NbrRow; ++j)
    {
      WriteLittleEndian(file, this->RowLastPointers[j]);
    }
  for (long j = 0l; j < this->NbrMatrixElements; ++j)
    {
      WriteLittleEndian(file, this->MatrixElements[j]);
    }
  for (long j = 0l; j < this->NbrMatrixElements; ++j)
    {
      WriteLittleEndian(file, this->ColumnIndices[j]);
    }
  return true;
}

// read matrix from a file 
//
// file = reference  on the input file stream
// return value = true if no error occurs

bool SparseRealMatrix::ReadMatrix (ifstream& file)
{
  int TmpType = Matrix::RealElements;
  ReadLittleEndian(file, TmpType);
  if ((this->MatrixType & TmpType & Matrix::RealElements) == 0)
    {
      file.close();
      return false;
    }
  if ((this->MatrixElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->MatrixElements;
	delete[] this->ColumnIndices;
	delete[] this->RowPointers;
	delete[] this->RowLastPointers;
      }
  this->Flag.Initialize();

  int TmpNbrRow;
  int TmpNbrColumn;
  ReadLittleEndian(file, TmpNbrRow);
  ReadLittleEndian(file, TmpNbrColumn);
  this->NbrRow = TmpNbrRow;
  this->NbrColumn = TmpNbrColumn;
  this->TrueNbrRow = TmpNbrRow;
  this->TrueNbrColumn = TmpNbrColumn;  
  long Tmp;
  ReadLittleEndian(file, Tmp);
  this->NbrMatrixElements = Tmp;
  this->MaximumNbrMatrixElements = Tmp;
  this->NbrMatrixElementPacketSize = 0;

  this->RowPointers = new long[this->NbrRow];
  this->RowLastPointers = new long[this->NbrRow];
  this->MatrixElements = new double[this->NbrMatrixElements];
  this->ColumnIndices = new int[this->NbrMatrixElements];
  for (int j = 0; j < this->NbrRow; ++j)
    {
      ReadLittleEndian(file, Tmp);
      this->RowPointers[j] = Tmp;
    }

  for (int j = 0; j < this->NbrRow; ++j)
    {
      ReadLittleEndian(file, Tmp);
      this->RowLastPointers[j] = Tmp;
    }
  double Tmp2;
  for (long j = 0l; j < this->NbrMatrixElements; ++j)
    {
      ReadLittleEndian(file, Tmp2);
      this->MatrixElements[j] = Tmp2;
    }
  for (long j = 0l; j < this->NbrMatrixElements; ++j)
    {
      ReadLittleEndian(file, TmpNbrColumn);
      this->ColumnIndices[j] =TmpNbrColumn;
    }
  return true;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const SparseRealMatrix& P) 
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      if (P.RowPointers[i] == -1l)
	{
	  for (int j = 0; j < P.NbrColumn; ++j)
	    Str << "0    ";
	}
      else
	{
	  long MinPos = P.RowPointers[i];
	  long MaxPos = P.RowLastPointers[i];
	  for (int j = 0; j < P.NbrColumn; ++j)
	    {
	      if ((MinPos <= MaxPos) && (P.ColumnIndices[MinPos] == j))
		{
		  Str << P.MatrixElements[MinPos] << "    ";
		  ++MinPos;
		}
	      else
		Str << "0    ";
	    }
	}
      Str << endl;
    }
  return Str;
}

// output the matrix in a sparse display (column formatted output)
//
// str = reference on output stream
// error = numerical accuracy below which a matrix element is considered to be equal to zero (discarded fro sparse matrices)
// return value = reference on output stream

ostream& SparseRealMatrix::PrintNonZero (ostream& str, double error) 
{
  for (long i = 0; i < this->NbrRow; ++i)
    if (this->RowPointers[i] >= 0l)
      {
	long MinPos = this->RowPointers[i];
	long MaxPos = this->RowLastPointers[i];
	while (MinPos <= MaxPos)
	  {
	    str << i << " " << this->ColumnIndices[MinPos] << " " << this->MatrixElements[MinPos] << endl;
	    ++MinPos;
	  }
      }
  return str;
}

// output the matrix in a sparse display (column formatted output), using labels for the row and column indices
//
// str = reference on output stream
// rowLabels = array of labels for the row indices
// columnLabels = array of labels for the column indices
// error = numerical accuracy below which a matrix element is considered to be equal to zero
// return value = reference on output stream  

ostream& SparseRealMatrix::PrintNonZero (ostream& str, char** rowLabels, char** columnLabels, double error) 
{
  for (long i = 0; i < this->NbrRow; ++i)
    if (this->RowPointers[i] >= 0l)
      {
	long MinPos = this->RowPointers[i];
	long MaxPos = this->RowLastPointers[i];
	while (MinPos <= MaxPos)
	  {
	    str << rowLabels[i] << " " << columnLabels[this->ColumnIndices[MinPos]] << " " << this->MatrixElements[MinPos] << endl;
	    ++MinPos;
	  }
      }
  return str;
}

#ifdef USE_OUTPUT

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// P = matrix to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, const SparseRealMatrix& P) 
{
  Str << "{";
  for (int i = 0; i < P.NbrRow; ++i)
    {
      Str << "{";
      if (P.RowPointers[0] == -1l)
	{
	  Str << "0,";
	}
      else
	{
	  long MinPos = P.RowPointers[i];
	  long MaxPos = P.RowLastPointers[i];
	  for (int j = 0; j < (P.NbrColumn - 1); ++j)
	    {
	      if ((MinPos <= MaxPos) && (P.ColumnIndices[MinPos] == j))
		{
		  Str << P.MatrixElements[MinPos] << ",";
		  ++MinPos;
		}
	      else
		Str << "0,";
	    }
	  if ((MinPos <= MaxPos) && (P.ColumnIndices[MinPos] == (P.NbrColumn - 1)))
	    {
	      Str << P.MatrixElements[MinPos];
	      ++MinPos;
	    }
	  else
	    Str << "0";	  
	  Str << "}";
	}
      Str << "}";
      if (i != (P.NbrRow - 1))
	Str << ",";
    }
  Str << "}";
  return Str;
}

#endif

//returns the array with indices of rows

void SparseRealMatrix::GetRowIndices(int* RowIndices)
{
  int counter = 0;
    for (long i = 0; i < this->NbrRow; ++i)
      for (long j = 0; j < this->NbrColumn; ++j)
        {
           double Tmp;
           this->GetMatrixElement(i,j,Tmp);
           if ((fabs(Tmp) != 0) && (this->MatrixElements[counter] == Tmp) && (counter < this->NbrMatrixElements))
             {
               RowIndices[counter] = i;
               counter++; 
             }
        }
}
