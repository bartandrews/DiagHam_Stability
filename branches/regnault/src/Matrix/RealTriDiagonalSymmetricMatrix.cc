////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of real tridiagoonal symmetric matrix                 //
//                                                                            //
//                        last modification : 25/01/2001                      //
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


#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealUpperTriangularMatrix.h"
#include "Matrix/RealLowerTriangularMatrix.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealMatrix.h"

#include <math.h>


using std::endl;


// default constructor
//

RealTriDiagonalSymmetricMatrix::RealTriDiagonalSymmetricMatrix() 
{
  this->DiagonalElements = 0;
  this->UpperDiagonalElements = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = 0;
  this->TrueNbrColumn = 0;
  this->MatrixType = Matrix::RealElements | Matrix::TriDiagonal | Matrix::Symmetric;
  this->Dummy = 0.0;
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

RealTriDiagonalSymmetricMatrix::RealTriDiagonalSymmetricMatrix(int dimension, bool zero)
{
  this->DiagonalElements = new double [dimension];
  this->UpperDiagonalElements = new double [dimension];
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = dimension;
  this->TrueNbrColumn = dimension;
  this->MatrixType = Matrix::RealElements | Matrix::TriDiagonal | Matrix::Symmetric;
  if (zero == true)
    {
      for (int i = 0; i < this->NbrRow; i++)
	{
	  this->DiagonalElements[i] = 0.0;
	  this->UpperDiagonalElements[i] = 0.0;
	}
    }
  this->Dummy = 0.0;
}

// constructor from matrix elements (without duplicating datas)
//
// diagonal = pointer to diagonal element array
// upperDiagonal = pointer to upper diagonal element arra
// dimension = matrix dimension

RealTriDiagonalSymmetricMatrix::RealTriDiagonalSymmetricMatrix(double* diagonal, double* upperDiagonal, int dimension) 
{
  this->DiagonalElements = diagonal;
  this->UpperDiagonalElements = upperDiagonal;
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = dimension;
  this->TrueNbrColumn = dimension;
  this->MatrixType = Matrix::RealElements | Matrix::TriDiagonal | Matrix::Symmetric;
  this->Dummy = 0.0;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

RealTriDiagonalSymmetricMatrix::RealTriDiagonalSymmetricMatrix(const RealTriDiagonalSymmetricMatrix& M) 
{  
  this->DiagonalElements = M.DiagonalElements;
  this->UpperDiagonalElements = M.UpperDiagonalElements;
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::TriDiagonal | Matrix::Symmetric;
  this->Dummy = 0.0;
}

// destructor
//

RealTriDiagonalSymmetricMatrix::~RealTriDiagonalSymmetricMatrix() 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->UpperDiagonalElements;
      delete[] this->DiagonalElements;
    }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

RealTriDiagonalSymmetricMatrix& RealTriDiagonalSymmetricMatrix::operator = (const RealTriDiagonalSymmetricMatrix& M) 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->UpperDiagonalElements;
      delete[] this->DiagonalElements;
    }
  this->DiagonalElements = M.DiagonalElements;
  this->UpperDiagonalElements = M.UpperDiagonalElements;
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::TriDiagonal | Matrix::Symmetric;
  this->Dummy = 0.0;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* RealTriDiagonalSymmetricMatrix::Clone ()
{
  return ((Matrix*) new RealTriDiagonalSymmetricMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void RealTriDiagonalSymmetricMatrix::SetMatrixElement(int i, int j, double x)
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
	  this->UpperDiagonalElements[i] = x;
	}
      else
	if ((j == -1) && (i < this->NbrRow))
	  {
	    this->UpperDiagonalElements[i - 1] = x;
	  }	
    }    
}

// return refernce on real part of a given matrix element
//
// i = line position
// j = column position
// return value = reference on real part

double& RealTriDiagonalSymmetricMatrix::operator () (int i, int j)
{
  if ((i == j) && (i < this->NbrRow))
    {
      return this->DiagonalElements[i];
    }
  else
    {
      j -= i;
      if ((j == 1) && (i < (this->NbrRow - 1)))
	{
	  return this->UpperDiagonalElements[i];
	}
      else
	if ((j == -1) && (i < this->NbrRow))
	  {
	    return this->UpperDiagonalElements[i - 1];
	  }	
    }    
  return this->Dummy;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void RealTriDiagonalSymmetricMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
}

// access to i-th diagonal element
// 
// i = position 
// return value = reference on i-th diagonal element

double& RealTriDiagonalSymmetricMatrix::DiagonalElement(int i)
{
  return this->DiagonalElements[i];
}

// access to i-th upper diagonal element
// 
// i = position 
// return value = reference on i-th upper diagonal element 

double& RealTriDiagonalSymmetricMatrix::UpperDiagonalElement(int i)
{
  return this->UpperDiagonalElements[i];
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void RealTriDiagonalSymmetricMatrix::AddToMatrixElement(int i, int j, double x)
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
	  this->UpperDiagonalElements[i] += x;
	}
      else
	if ((j == -1) && (i < this->NbrRow))
	  {
	    this->UpperDiagonalElements[i - 1] += x;
	  }	
    }    
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void RealTriDiagonalSymmetricMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void RealTriDiagonalSymmetricMatrix::Resize (int nbrRow, int nbrColumn)
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
  double* TmpUpperDiag = new double [nbrRow];
  if (this->Flag.Used() == true)
    {
      for (int i = 0; i < this->NbrRow; i++)
	{
	  TmpDiag[i] = this->DiagonalElements[i];
	  TmpUpperDiag[i] = this->UpperDiagonalElements[i]; 
	}
    }
   if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->UpperDiagonalElements;
      delete[] this->DiagonalElements;
    }
  this->DiagonalElements = TmpDiag;
  this->UpperDiagonalElements = TmpUpperDiag;
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = nbrRow;
  this->TrueNbrColumn = nbrColumn;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
  return;
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void RealTriDiagonalSymmetricMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  if (nbrRow <= this->TrueNbrRow)
    {
      for (int i = this->NbrRow; i < nbrRow; i++)
	{
	  this->DiagonalElements[i] = 0.0;
	  this->UpperDiagonalElements[i] = 0.0;
	}
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      return;
    }
  double* TmpDiag = new double [nbrRow];
  double* TmpUpperDiag = new double [nbrRow];
  if (this->Flag.Used() == true)
    {
      int i = 0;
      for (; i < this->NbrRow; i++)
	{
	  TmpDiag[i] = this->DiagonalElements[i];
	  TmpUpperDiag[i] = this->UpperDiagonalElements[i]; 
	}
      for (; i < nbrRow; i++)
	{
	  TmpDiag[i] = 0.0;
	  TmpUpperDiag[i] = 0.0; 
	}
    }
   if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->UpperDiagonalElements;
      delete[] this->DiagonalElements;
    }
  this->DiagonalElements = TmpDiag;
  this->UpperDiagonalElements = TmpUpperDiag;
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = nbrRow;
  this->TrueNbrColumn = nbrColumn;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
  return;
}

// copy matrix
//
// M = matrix to copy
// return value = refence on current matrix

RealTriDiagonalSymmetricMatrix& RealTriDiagonalSymmetricMatrix::Copy (RealTriDiagonalSymmetricMatrix& M)
{
  if (this->NbrRow != M.NbrRow)
    this->Resize(M.NbrRow, M.NbrColumn);
  for (int i = 0; i < M.NbrColumn; i++)
    {
      this->DiagonalElements[i] = M.DiagonalElements[i];
      this->UpperDiagonalElements[i] = M.UpperDiagonalElements[i];
    }
  return *this;
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

RealTriDiagonalSymmetricMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2) 
{
  if (M1.NbrRow != M2.NbrRow)
    return RealTriDiagonalSymmetricMatrix();
  int ReducedNbr = M1.NbrRow - 1;
  double* Diagonal = new double [M1.NbrRow];
  double* UpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; i++)
    {
      Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
      UpperDiagonal[i] = M1.UpperDiagonalElements[i] + M2.UpperDiagonalElements[i];      
    }
  Diagonal[ReducedNbr] = M1.DiagonalElements[ReducedNbr] + M2.DiagonalElements[ReducedNbr];
  return RealTriDiagonalSymmetricMatrix(Diagonal, UpperDiagonal, M1.NbrRow);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

RealTriDiagonalSymmetricMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2) 
{
  if (M1.NbrRow != M2.NbrRow)
    return RealTriDiagonalSymmetricMatrix();
  int ReducedNbr = M1.NbrRow - 1;
  double* Diagonal = new double [M1.NbrRow];
  double* UpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; i++)
    {
      Diagonal[i] = M1.DiagonalElements[i] - M2.DiagonalElements[i];
      UpperDiagonal[i] = M1.UpperDiagonalElements[i] - M2.UpperDiagonalElements[i];      
    }
  Diagonal[ReducedNbr] = M1.DiagonalElements[ReducedNbr] - M2.DiagonalElements[ReducedNbr];
  return RealTriDiagonalSymmetricMatrix(Diagonal, UpperDiagonal, M1.NbrRow);
}

// multiply a matrix with a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealTriDiagonalSymmetricMatrix operator * (const RealTriDiagonalSymmetricMatrix& M, double x) 
{
  int ReducedNbr = M.NbrRow - 1;
  double* Diagonal = new double [M.NbrRow];
  double* UpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
      UpperDiagonal[i] = M.UpperDiagonalElements[i] * x;      
    }
  Diagonal[ReducedNbr] = M.DiagonalElements[ReducedNbr] * x;
  return RealTriDiagonalSymmetricMatrix(Diagonal, UpperDiagonal, M.NbrRow);
}

// multiply a matrix with a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealTriDiagonalSymmetricMatrix operator * (double x, const RealTriDiagonalSymmetricMatrix& M) 
{
  int ReducedNbr = M.NbrRow - 1;
  double* Diagonal = new double [M.NbrRow];
  double* UpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
      UpperDiagonal[i] = M.UpperDiagonalElements[i] * x;      
    }
  Diagonal[ReducedNbr] = M.DiagonalElements[ReducedNbr] * x;
  return RealTriDiagonalSymmetricMatrix(Diagonal, UpperDiagonal, M.NbrRow);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

RealTriDiagonalSymmetricMatrix operator / (const RealTriDiagonalSymmetricMatrix& M, double x) 
{
  int ReducedNbr = M.NbrRow - 1;
  double* Diagonal = new double [M.NbrRow];
  double* UpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] / x;
      UpperDiagonal[i] = M.UpperDiagonalElements[i] / x;      
    }
  Diagonal[ReducedNbr] = M.DiagonalElements[ReducedNbr] / x;
  return RealTriDiagonalSymmetricMatrix(Diagonal, UpperDiagonal, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

RealTriDiagonalSymmetricMatrix& RealTriDiagonalSymmetricMatrix::operator += (const RealTriDiagonalSymmetricMatrix& M) 
{
  if (this->NbrRow != M.NbrRow)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  for (int i = 0; i < ReducedNbr; i++)
    {
      this->DiagonalElements[i] += M.DiagonalElements[i];
      this->UpperDiagonalElements[i] += M.UpperDiagonalElements[i];      
    }
  this->DiagonalElements[ReducedNbr] += DiagonalElements[ReducedNbr];
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

RealTriDiagonalSymmetricMatrix& RealTriDiagonalSymmetricMatrix::operator -= (const RealTriDiagonalSymmetricMatrix& M) 
{
  if (this->NbrRow != M.NbrRow)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  for (int i = 0; i < ReducedNbr; i++)
    {
      this->DiagonalElements[i] -= M.DiagonalElements[i];
      this->UpperDiagonalElements[i] -= M.UpperDiagonalElements[i];      
    }
  this->DiagonalElements[ReducedNbr] -= DiagonalElements[ReducedNbr];
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

RealTriDiagonalSymmetricMatrix& RealTriDiagonalSymmetricMatrix::operator *= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  int ReducedNbr = this->NbrRow - 1;
  for (int i = 0; i < ReducedNbr; i++)
    {
      this->DiagonalElements[i] *= x;
      this->UpperDiagonalElements[i] *= x;      
    }
  this->DiagonalElements[ReducedNbr] *= x;
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

RealTriDiagonalSymmetricMatrix& RealTriDiagonalSymmetricMatrix::operator /= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  int ReducedNbr = this->NbrRow - 1;
  for (int i = 0; i < ReducedNbr; i++)
    {
      this->DiagonalElements[i] /= x;
      this->UpperDiagonalElements[i] /= x;      
    }
  this->DiagonalElements[ReducedNbr] /= x;
  return *this;
}

// get a matrix element 
// 
// i = Row number
// j = Column number
// return value = matrix element M_(i,j)

double RealTriDiagonalSymmetricMatrix::GetElement(int i, int j)
{
  if (i == j)
    return this->DiagonalElements[i];
  if (j < i)
    {
      int tmp = i;
      i = j;
      j = tmp;
    }
  if ((j - i) == 1)
    return this->UpperDiagonalElements[i];
  return 0.0;
}

// evaluate matrix trace
//
// return value = matrix trace 

double RealTriDiagonalSymmetricMatrix::Tr () 
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

double RealTriDiagonalSymmetricMatrix::Det () 
{
  if (this->NbrRow == 0)
    return 0.0;
  double d0 = this->DiagonalElements[0];
  if (this->NbrRow == 1)
    return d0;
  double d1 = d0 * this->DiagonalElements[0] - this->UpperDiagonalElements[0] * this->UpperDiagonalElements[0];
  if (this->NbrRow == 2)
    return d0;
  double d;
  for (int i = 2; i < this->NbrRow; i++)
    {
      d = this->DiagonalElements[i] * d1 - this->UpperDiagonalElements[i - 1] * this->UpperDiagonalElements[i - 1] * d0;
      d0 = d1;
      d1 = d;
    }
  return d1;
}

// return matrix characteritic equation
//
// return value =  reference one polynomial corresponding to matrix characteritic equation

Polynomial& RealTriDiagonalSymmetricMatrix::CharacteristicEquation()
{
  double* P0 = new double [this->NbrRow + 1];
  double* P1 = new double [this->NbrRow + 1];
  double* TmpP;
  P0[1] = -1;
  P0[0] = this->DiagonalElements[0];
  P1[2] = 1;
  P1[1] = - this->DiagonalElements[0] - this->DiagonalElements[1];
  P1[0] = - this->UpperDiagonalElements[0] * this->UpperDiagonalElements[0] + this->DiagonalElements[0] * this->DiagonalElements[1];

  for (int i = 2; i < this->NbrRow; i++)
    {
      int j = (i  - 1);
      double fac = - this->UpperDiagonalElements[j] * this->UpperDiagonalElements[j];
      for (int k = 0; k < i; k++)
	P0[k] *= fac;
      P0[0] += this->DiagonalElements[i] * P1[0];
      for (int k = 1; k < i; k++)
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

// Diagonalize RealTridiagonal Symmetric Matrix using QL algorithm with implicit shift
// current matrix is replaced by its corresponding diagonalized matrix
//
// maxIter = maximum number of iteration to find an eigenvalue
// return value = reference on current Matrix

RealTriDiagonalSymmetricMatrix& RealTriDiagonalSymmetricMatrix::Diagonalize(int maxIter)
{
  int ReducedDimension = this->NbrRow - 1;
  double Cos;
  double Sin;
  double Theta;
  double T, R, P, F, B;

  for (int i = 0; i < ReducedDimension; i++)
    {
      int iter = 0;
      while (iter < maxIter)
	{
	  // find bloc matrices so that QL algorithm will be applied on submatrix from i to j
	  int j = i;
	  bool Flag = false;
	  while ((j < ReducedDimension) && (Flag == false))
	    {
	      double d2 = fabs(this->DiagonalElements[j]) + fabs(this->DiagonalElements[j + 1]);
	      if ((d2 + fabs(this->UpperDiagonalElements[j])) == d2)
		Flag = true;
	      else
		j++;
	    }
	  // if i != j, i-th eigenvalue has not been obtained yet, apply diagonalization on submatrix      
	  if (j != i)
	    {
	      iter++;
	      // evaluate shift
	      Theta = (this->DiagonalElements[i + 1] - this->DiagonalElements[i]) 
		/ (2.0 * this->UpperDiagonalElements[i]);
	      R = sqrt (1.0 + Theta * Theta);
	      T = this->DiagonalElements[j] - this->DiagonalElements[i];	      	      
	      if (Theta >= 0)
		T += this->UpperDiagonalElements[i] / (Theta + R);
	      else
		T += this->UpperDiagonalElements[i] / (Theta - R);
	      Cos = 1.0;
	      Sin = 1.0;
	      P = 0.0;
	      // apply shift and conjugation with Jacobi and Givens rotations
	      for (int k = j - 1; k >= i; k--)
		{
		  F = Sin * this->UpperDiagonalElements[k];
		  B = Cos * this->UpperDiagonalElements[k];
		  R = sqrt (F * F + T * T);
		  this->UpperDiagonalElements[k + 1] = R;
		  if (R == 0.0)
		    {
		      this->DiagonalElements[k + 1] -= P;
		      this->UpperDiagonalElements[j] = 0.0;
		      k = i - 1;
		    }
		  else
		    {
		      Sin = 1.0 / R;
		      Cos = Sin * T;
		      Sin *= F;
		      T = this->DiagonalElements[k + 1] - P;
		      R = (this->DiagonalElements[k] - T) * Sin + 2.0 * Cos * B;
		      P = Sin * R;
		      this->DiagonalElements[k + 1] = T + P;
		      T = Cos * R - B;
		    }
		}
	      this->DiagonalElements[i] -= P;
	      this->UpperDiagonalElements[i] = T;
	      this->UpperDiagonalElements[j] = 0.0;
	    }
	  else
	    iter = maxIter;
	}
    }
  return *this;
}

// Diagonalize RealTridiagonal Symmetric Matrix using QL algorithm with implicit shift, evaluating eigenvectors in a given base
// current matrix is replaced by its corresponding diagonalized matrix
//
// Q = matrix initialized with corresponding base in which eigenvectors have to be calculated
// maxIter = maximum number of iteration to fund an eigenvalue
// return value = reference on current Matrix

RealTriDiagonalSymmetricMatrix& RealTriDiagonalSymmetricMatrix::Diagonalize(ComplexMatrix& Q, int maxIter)
{
  int ReducedDimension = this->NbrRow - 1;
  double Cos;
  double Sin;
  double Theta;
  double T, R, P, F, B;

  for (int i = 0; i < ReducedDimension; i++)
    {
      int iter = 0;
      while (iter < maxIter)
	{
	  // find bloc matrices so that QL algorithm will be applied on submatrix from i to j
	  int j = i;
	  bool Flag = false;
	  while ((j < ReducedDimension) && (Flag == false))
	    {
	      double d2 = fabs(this->DiagonalElements[j]) + fabs(this->DiagonalElements[j + 1]);
	      if ((d2 + fabs(this->UpperDiagonalElements[j])) == d2)
		Flag = true;
	      else
		j++;
	    }
	  // if i != j, i-th eigenvalue has not been obtained yet, apply diagonalization on submatrix      
	  if (j != i)
	    {
	      iter++;
	      // evaluate shift
	      Theta = (this->DiagonalElements[i + 1] - this->DiagonalElements[i]) 
		/ (2.0 * this->UpperDiagonalElements[i]);
	      R = sqrt (1.0 + Theta * Theta);
	      T = this->DiagonalElements[j] - this->DiagonalElements[i];	      	      
	      if (Theta >= 0)
		T += this->UpperDiagonalElements[i] / (Theta + R);
	      else
		T += this->UpperDiagonalElements[i] / (Theta - R);
	      Cos = 1.0;
	      Sin = 1.0;
	      P = 0.0;
	      // apply shift and conjugation with Jacobi and Givens rotations
	      for (int k = j - 1; k >= i; k--)
		{
		  F = Sin * this->UpperDiagonalElements[k];
		  B = Cos * this->UpperDiagonalElements[k];
		  R = sqrt (F * F + T * T);
		  this->UpperDiagonalElements[k + 1] = R;
		  if (R == 0.0)
		    {
		      this->DiagonalElements[k + 1] -= P;
		      this->UpperDiagonalElements[j] = 0.0;
		      k = i - 1;
		    }
		  else
		    {
		      Sin = 1.0 / R;
		      Cos = Sin * T;
		      Sin *= F;
		      T = this->DiagonalElements[k + 1] - P;
		      R = (this->DiagonalElements[k] - T) * Sin + 2.0 * Cos * B;
		      P = Sin * R;
		      this->DiagonalElements[k + 1] = T + P;
		      T = Cos * R - B;
		      // apply transformation to vectors
		      double tmp;
		      for (int n = 0; n < Q.NbrRow; n++)
			{
			  tmp = Q.Columns[k + 1].RealComponents[n];
			  Q.Columns[k + 1].RealComponents[n] *= Cos;
			  Q.Columns[k + 1].RealComponents[n] += Sin * Q.Columns[k].RealComponents[n];
			  Q.Columns[k].RealComponents[n] *= Cos;
			  Q.Columns[k].RealComponents[n] -= Sin * tmp;
			  tmp = Q.Columns[k + 1].ImaginaryComponents[n];
			  Q.Columns[k + 1].ImaginaryComponents[n] *= Cos;
			  Q.Columns[k + 1].ImaginaryComponents[n] += Sin * Q.Columns[k].ImaginaryComponents[n];
			  Q.Columns[k].ImaginaryComponents[n] *= Cos;
			  Q.Columns[k].ImaginaryComponents[n] -= Sin * tmp;
			}
		    }
		}
	      this->DiagonalElements[i] -= P;
	      this->UpperDiagonalElements[i] = T;
	      this->UpperDiagonalElements[j] = 0.0;
	    }
	  else
	    iter = maxIter;
	}
    }
  return *this;
}

// Diagonalize RealTridiagonal Symmetric Matrix using QL algorithm with implicit shift, evaluating eigenvectors in a given base
// current matrix is replaced by its corresponding diagonalized matrix
//
// Q = matrix initialized with corresponding base in which eigenvectors have to be calculated
// maxIter = maximum number of iteration to fund an eigenvalue
// return value = reference on current Matrix

RealTriDiagonalSymmetricMatrix& RealTriDiagonalSymmetricMatrix::Diagonalize(RealMatrix& Q, int maxIter)
{
  int ReducedDimension = this->NbrRow - 1;
  double Cos;
  double Sin;
  double Theta;
  double T, R, P, F, B;

  for (int i = 0; i < ReducedDimension; i++)
    {
      int iter = 0;
      while (iter < maxIter)
	{
	  // find bloc matrices so that QL algorithm will be applied on submatrix from i to j
	  int j = i;
	  bool Flag = false;
	  while ((j < ReducedDimension) && (Flag == false))
	    {
	      double d2 = fabs(this->DiagonalElements[j]) + fabs(this->DiagonalElements[j + 1]);
	      if ((d2 + fabs(this->UpperDiagonalElements[j])) == d2)
		Flag = true;
	      else
		j++;
	    }
	  // if i != j, i-th eigenvalue has not been obtained yet, apply diagonalization on submatrix      
	  if (j != i)
	    {
	      iter++;
	      // evaluate shift
	      Theta = (this->DiagonalElements[i + 1] - this->DiagonalElements[i]) 
		/ (2.0 * this->UpperDiagonalElements[i]);
	      R = sqrt (1.0 + Theta * Theta);
	      T = this->DiagonalElements[j] - this->DiagonalElements[i];	      	      
	      if (Theta >= 0)
		T += this->UpperDiagonalElements[i] / (Theta + R);
	      else
		T += this->UpperDiagonalElements[i] / (Theta - R);
	      Cos = 1.0;
	      Sin = 1.0;
	      P = 0.0;
	      // apply shift and conjugation with Jacobi and Givens rotations
	      for (int k = j - 1; k >= i; k--)
		{
		  F = Sin * this->UpperDiagonalElements[k];
		  B = Cos * this->UpperDiagonalElements[k];
		  R = sqrt (F * F + T * T);
		  this->UpperDiagonalElements[k + 1] = R;
		  if (R == 0.0)
		    {
		      this->DiagonalElements[k + 1] -= P;
		      this->UpperDiagonalElements[j] = 0.0;
		      k = i - 1;
		    }
		  else
		    {
		      Sin = 1.0 / R;
		      Cos = Sin * T;
		      Sin *= F;
		      T = this->DiagonalElements[k + 1] - P;
		      R = (this->DiagonalElements[k] - T) * Sin + 2.0 * Cos * B;
		      P = Sin * R;
		      this->DiagonalElements[k + 1] = T + P;
		      T = Cos * R - B;
		      // apply transformation to vectors
		      double tmp;
		      for (int n = 0; n < Q.NbrRow; n++)
			{
			  tmp = Q.Columns[k + 1].Components[n];
			  Q.Columns[k + 1].Components[n] *= Cos;
			  Q.Columns[k + 1].Components[n] += Sin * Q.Columns[k].Components[n];
			  Q.Columns[k].Components[n] *= Cos;
			  Q.Columns[k].Components[n] -= Sin * tmp;
			}
		    }
		}
	      this->DiagonalElements[i] -= P;
	      this->UpperDiagonalElements[i] = T;
	      this->UpperDiagonalElements[j] = 0.0;
	    }
	  else
	    iter = maxIter;
	}
    }
  return *this;
}

// find QR factorization
//
// Q = matrix initialized with corresponding base in which unitary matrix of QR factorization has to be calculated
// return value = upper triangular matrix corresponding to the QR factorization of the matrix

RealUpperTriangularMatrix RealTriDiagonalSymmetricMatrix::QRFactorization(RealMatrix& Q)
{
  RealUpperTriangularMatrix TmpR(this->NbrRow, true);
  double Tmp = 0.0;
  double CosTheta = 0.0;
  double SinTheta = 0.0;
  double PreviousDiagonal = this->DiagonalElements[0];
  double OffPreviousDiagonal = this->UpperDiagonalElements[0];  
  int ReducedSize = this->NbrRow - 2;  
  int Inci = 1;
  int i = 0;
  for (; i < ReducedSize; ++i)
    {
      // find Householder transformation and apply it to lower diagonal element
      Tmp =  sqrt ((PreviousDiagonal * PreviousDiagonal) + (this->UpperDiagonalElements[i] *  this->UpperDiagonalElements[i]));
      TmpR[i] = Tmp;
      Tmp = 1.0 / Tmp;
      CosTheta = PreviousDiagonal * Tmp;
      SinTheta = this->UpperDiagonalElements[i] * Tmp;
      // apply Householder transformation to the other elements
      TmpR(i, Inci) = OffPreviousDiagonal * CosTheta + this->DiagonalElements[Inci] * SinTheta;
      TmpR(i, i + 2) = this->UpperDiagonalElements[Inci] * SinTheta;
      PreviousDiagonal = OffPreviousDiagonal * SinTheta - this->DiagonalElements[Inci] * CosTheta;
      OffPreviousDiagonal = -this->UpperDiagonalElements[Inci] * CosTheta;
      // apply transformation on the Q matrix
      for (int j = 0; j < this->NbrRow; ++j)
	{
	  Tmp = Q(j, i);
	  Q(j, i) *= CosTheta;
	  Q(j, i) += Q(j, Inci) * SinTheta;
	  Q(j, Inci) *= -CosTheta;
	  Q(j, Inci) += Tmp * SinTheta;
	}
      ++Inci;
    }
  // find Householder transformation and apply it to lower diagonal element
  Tmp =  sqrt ((PreviousDiagonal * PreviousDiagonal) + (this->UpperDiagonalElements[i] *  this->UpperDiagonalElements[i]));
  TmpR[i] = Tmp;
  Tmp = 1.0 / Tmp;
  CosTheta = PreviousDiagonal * Tmp;
  SinTheta = this->UpperDiagonalElements[i] * Tmp;
  // apply Householder transformation to the other elements
  TmpR(i, Inci) = OffPreviousDiagonal * CosTheta + this->DiagonalElements[Inci] * SinTheta;
  TmpR[Inci] = OffPreviousDiagonal * SinTheta - this->DiagonalElements[Inci] * CosTheta;
  // apply transformation on the Q matrix
  for (int j = 0; j < this->NbrRow; ++j)
    {
      Tmp = Q(j, i);
      Q(j, i) *= CosTheta;
      Q(j, i) += Q(j, Inci) * SinTheta;
      Q(j, Inci) *= -CosTheta;
      Q(j, Inci) += Tmp * SinTheta;
    }

  return TmpR;
}

// find QL factorization
//
// Q = matrix initialized with corresponding base in which unitary matrix of QR factorization has to be calculated
// return value = lower triangular matrix corresponding to the QL factorization of the matrix

RealLowerTriangularMatrix RealTriDiagonalSymmetricMatrix::QLFactorization(RealMatrix& Q)
{
  RealLowerTriangularMatrix TmpL(this->NbrRow, true);
  double Tmp = 0.0;
  double CosTheta = 0.0;
  double SinTheta = 0.0;
  double PreviousDiagonal = this->DiagonalElements[this->NbrRow - 1];
  double OffPreviousDiagonal = this->UpperDiagonalElements[this->NbrRow - 2];  
  for (int i = this->NbrRow - 2; i > 0; --i)
    {
      // find Householder transformation and apply it to lower diagonal element
      Tmp =  sqrt ((PreviousDiagonal * PreviousDiagonal) + (this->UpperDiagonalElements[i] *  this->UpperDiagonalElements[i]));
      TmpL[i + 1] = Tmp;
      Tmp = 1.0 / Tmp;
      CosTheta = -PreviousDiagonal * Tmp;
      SinTheta = this->UpperDiagonalElements[i] * Tmp;
      // apply Householder transformation to the other elements
      TmpL(i + 1, i) = this->DiagonalElements[i] * SinTheta - CosTheta * OffPreviousDiagonal;
      TmpL(i + 1, i - 1) = this->UpperDiagonalElements[i - 1] * SinTheta;
      PreviousDiagonal = this->DiagonalElements[i] * CosTheta + SinTheta * OffPreviousDiagonal;
      OffPreviousDiagonal = this->UpperDiagonalElements[i - 1] * CosTheta;
      // apply transformation on the Q matrix
      for (int j = 0; j < this->NbrRow; ++j)
	{
	  Tmp = Q(j, i);
	  Q(j, i) *= CosTheta;
	  Q(j, i) += Q(j, i + 1) * SinTheta;
	  Q(j, i + 1) *= -CosTheta;
	  Q(j, i + 1) += Tmp * SinTheta;
	}
    }
  // evaluate last transformation
  Tmp =  sqrt ((PreviousDiagonal * PreviousDiagonal) + (this->UpperDiagonalElements[0] *  this->UpperDiagonalElements[0]));
  TmpL[1] = Tmp;
  Tmp = 1.0 / Tmp;
  CosTheta = -PreviousDiagonal * Tmp;
  SinTheta = this->UpperDiagonalElements[0] * Tmp;
  TmpL(1, 0) = this->DiagonalElements[0] * SinTheta - CosTheta * OffPreviousDiagonal;
  TmpL[0] = this->DiagonalElements[0] * CosTheta + SinTheta * OffPreviousDiagonal;
  for (int j = 0; j < this->NbrRow; ++j)
    {
      Tmp = Q(j, 0);
      Q(j, 0) *= CosTheta;
      Q(j, 0) += Q(j, 1) * SinTheta;
      Q(j, 1) *= -CosTheta;
      Q(j, 1) += Tmp * SinTheta;
    }
  return TmpL;
}

// find QL factorization with shift (aka M - x 1) 
//
// Q = matrix initialized with corresponding base in which unitary matrix of QR factorization has to be calculated
// shift = shift value
// return value = lower triangular matrix corresponding to the QL factorization of the matrix

RealLowerTriangularMatrix RealTriDiagonalSymmetricMatrix::QLFactorization(RealMatrix& Q, double shift)
{
  RealLowerTriangularMatrix TmpL(this->NbrRow, true);
  double Tmp = 0.0;
  double CosTheta = 0.0;
  double SinTheta = 0.0;
  double PreviousDiagonal = this->DiagonalElements[this->NbrRow - 1] - shift;
  double OffPreviousDiagonal = this->UpperDiagonalElements[this->NbrRow - 2];  
  for (int i = this->NbrRow - 2; i > 0; --i)
    {
      // find Householder transformation and apply it to lower diagonal element
      Tmp =  sqrt ((PreviousDiagonal * PreviousDiagonal) + (this->UpperDiagonalElements[i] *  this->UpperDiagonalElements[i]));
      TmpL[i + 1] = Tmp;
      Tmp = 1.0 / Tmp;
      CosTheta = -PreviousDiagonal * Tmp;
      SinTheta = this->UpperDiagonalElements[i] * Tmp;
      // apply Householder transformation to the other elements
      TmpL(i + 1, i) = (this->DiagonalElements[i]  - shift) * SinTheta - CosTheta * OffPreviousDiagonal;
      TmpL(i + 1, i - 1) = this->UpperDiagonalElements[i - 1] * SinTheta;
      PreviousDiagonal = (this->DiagonalElements[i]  - shift) * CosTheta + SinTheta * OffPreviousDiagonal;
      OffPreviousDiagonal = this->UpperDiagonalElements[i - 1] * CosTheta;
      // apply transformation on the Q matrix
      for (int j = 0; j < this->NbrRow; ++j)
	{
	  Tmp = Q(j, i);
	  Q(j, i) *= CosTheta;
	  Q(j, i) += Q(j, i + 1) * SinTheta;
	  Q(j, i + 1) *= -CosTheta;
	  Q(j, i + 1) += Tmp * SinTheta;
	}
    }
  // evaluate last transformation
  Tmp =  sqrt ((PreviousDiagonal * PreviousDiagonal) + (this->UpperDiagonalElements[0] *  this->UpperDiagonalElements[0]));
  TmpL[1] = Tmp;
  Tmp = 1.0 / Tmp;
  CosTheta = -PreviousDiagonal * Tmp;
  SinTheta = this->UpperDiagonalElements[0] * Tmp;
  TmpL(1, 0) = (this->DiagonalElements[0]  - shift) * SinTheta - CosTheta * OffPreviousDiagonal;
  TmpL[0] = (this->DiagonalElements[0]  - shift) * CosTheta + SinTheta * OffPreviousDiagonal;
  for (int j = 0; j < this->NbrRow; ++j)
    {
      Tmp = Q(j, 0);
      Q(j, 0) *= CosTheta;
      Q(j, 0) += Q(j, 1) * SinTheta;
      Q(j, 1) *= -CosTheta;
      Q(j, 1) += Tmp * SinTheta;
    }
  return TmpL;
}

// find QL factorization and evaluate LQ (ie Qt H Q)
//
// Q = matrix initialized with corresponding base in which unitary matrix of QR factorization has to be calculated
// return value = Qt H Q

RealTriDiagonalSymmetricMatrix RealTriDiagonalSymmetricMatrix::QLConjugaison(RealMatrix& Q)
{
  RealTriDiagonalSymmetricMatrix TmpL(this->NbrRow, true);
  double Tmp = 0.0;
  double CosTheta = 0.0;
  double SinTheta = 0.0;
  double PreviousSinTheta = 0.0;
  int i = this->NbrRow - 2;
  int Inci = i + 1;
  // perform transformation on the latest elements of the tridiagonal matrix
  // find Householder transformation and apply it to lower diagonal element
  Tmp =  1.0 / sqrt ((this->DiagonalElements[Inci] * this->DiagonalElements[Inci]) + (this->UpperDiagonalElements[i] *  this->UpperDiagonalElements[i]));
  CosTheta = - this->DiagonalElements[Inci] * Tmp;  
  SinTheta = this->UpperDiagonalElements[i] * Tmp;
  TmpL.UpperDiagonalElements[i] = this->DiagonalElements[i] * CosTheta + SinTheta * this->UpperDiagonalElements[i];
  TmpL.DiagonalElements[i]  =  TmpL.UpperDiagonalElements[i] * CosTheta;
  TmpL.DiagonalElements[Inci] = this->DiagonalElements[Inci] + ((this->DiagonalElements[i] + this->DiagonalElements[Inci]) * SinTheta * SinTheta);
  TmpL.UpperDiagonalElements[i - 1] = this->UpperDiagonalElements[i - 1] * CosTheta;
  PreviousSinTheta = SinTheta;
  // apply transformation on the Q matrix
  for (int j = 0; j < this->NbrRow; ++j)
    {
      Tmp = Q(j, i);
      Q(j, i) *= CosTheta;
      Q(j, i) += Q(j, Inci) * SinTheta;
      Q(j, Inci) *= -CosTheta;
      Q(j, Inci) += Tmp * SinTheta;
    }
  --i;
  --Inci;
  for (; i > 0; --i)
    {
      // find Householder transformation
      Tmp = 1.0 / sqrt ((TmpL.UpperDiagonalElements[Inci] * TmpL.UpperDiagonalElements[Inci]) + (this->UpperDiagonalElements[i] *  this->UpperDiagonalElements[i]));
      CosTheta = -TmpL.UpperDiagonalElements[Inci] * Tmp;
      SinTheta = this->UpperDiagonalElements[i] * Tmp;
      // apply Householder transformation to the other elements

      TmpL.DiagonalElements[Inci] = TmpL.DiagonalElements[Inci] + ((this->DiagonalElements[i] + TmpL.DiagonalElements[Inci]) * SinTheta * SinTheta);
      TmpL.UpperDiagonalElements[i] *= SinTheta;
      TmpL.UpperDiagonalElements[i] += this->DiagonalElements[i] * CosTheta;
      TmpL.DiagonalElements[i] = TmpL.UpperDiagonalElements[i] * CosTheta;
      TmpL.UpperDiagonalElements[i - 1] = this->UpperDiagonalElements[i - 1] * CosTheta;
      TmpL.UpperDiagonalElements[Inci] *= -CosTheta;
      TmpL.UpperDiagonalElements[Inci] += this->UpperDiagonalElements[i] * SinTheta;
      TmpL.UpperDiagonalElements[Inci] *= PreviousSinTheta;
      PreviousSinTheta = SinTheta;
      // apply transformation on the Q matrix
      for (int j = 0; j < this->NbrRow; ++j)
	{
	  Tmp = Q(j, i);
	  Q(j, i) *= CosTheta;
	  Q(j, i) += Q(j, Inci) * SinTheta;
	  Q(j, Inci) *= -CosTheta;
	  Q(j, Inci) += Tmp * SinTheta;
	}
      --Inci;
    }
  // evaluate last transformation
  Tmp = 1.0 / sqrt ((TmpL.UpperDiagonalElements[1] * TmpL.UpperDiagonalElements[1]) + (this->UpperDiagonalElements[0] *  this->UpperDiagonalElements[0]));
  CosTheta = -TmpL.UpperDiagonalElements[1] * Tmp;
  SinTheta = this->UpperDiagonalElements[0] * Tmp;
  TmpL.DiagonalElements[1] = TmpL.DiagonalElements[1] + ((this->DiagonalElements[0] + TmpL.DiagonalElements[1]) * SinTheta * SinTheta);
  TmpL.UpperDiagonalElements[0] *= SinTheta;
  TmpL.UpperDiagonalElements[0] += this->DiagonalElements[0] * CosTheta;
  TmpL.DiagonalElements[0] = TmpL.UpperDiagonalElements[0] * CosTheta;
  TmpL.UpperDiagonalElements[0] *= SinTheta;
  TmpL.UpperDiagonalElements[1] *= -CosTheta;
  TmpL.UpperDiagonalElements[1] += this->UpperDiagonalElements[0] * SinTheta;
  TmpL.UpperDiagonalElements[1] *= PreviousSinTheta;
  for (int j = 0; j < this->NbrRow; ++j)
    {
      Tmp = Q(j, 0);
      Q(j, 0) *= CosTheta;
      Q(j, 0) += Q(j, 1) * SinTheta;
      Q(j, 1) *= -CosTheta;
      Q(j, 1) += Tmp * SinTheta;
    }
  return TmpL;
}

// find QL factorization and evaluate LQ (ie Qt H Q), shifting initial matrix diagonal elements
//
// Q = matrix initialized with corresponding base in which unitary matrix of QR factorization has to be calculated
// shift = shift value
// return value = Qt H Q

RealTriDiagonalSymmetricMatrix RealTriDiagonalSymmetricMatrix::QLConjugaison(RealMatrix& Q, double shift)
{
  RealTriDiagonalSymmetricMatrix TmpL(this->NbrRow, true);
  double Tmp = 0.0;
  double CosTheta = 0.0;
  double SinTheta = 0.0;
  double PreviousSinTheta = 0.0;
  int i = this->NbrRow - 2;
  int Inci = i + 1;
  // perform transformation on the latest elements of the tridiagonal matrix
  // find Householder transformation and apply it to lower diagonal element
  Tmp =  1.0 / sqrt (((this->DiagonalElements[Inci] - shift) * (this->DiagonalElements[Inci] - shift)) + (this->UpperDiagonalElements[i] *  this->UpperDiagonalElements[i]));
  CosTheta = - (this->DiagonalElements[Inci]  - shift) * Tmp;  
  SinTheta = this->UpperDiagonalElements[i] * Tmp;
  TmpL.UpperDiagonalElements[i] = (this->DiagonalElements[i] - shift) * CosTheta + SinTheta * this->UpperDiagonalElements[i];
  TmpL.DiagonalElements[i]  =  TmpL.UpperDiagonalElements[i] * CosTheta;
  TmpL.DiagonalElements[Inci] = (this->DiagonalElements[Inci] - shift) + ((this->DiagonalElements[i] + this->DiagonalElements[Inci] - 2.0 * shift) * SinTheta * SinTheta);
  TmpL.UpperDiagonalElements[i - 1] = this->UpperDiagonalElements[i - 1] * CosTheta;
  PreviousSinTheta = SinTheta;
  // apply transformation on the Q matrix
  for (int j = 0; j < this->NbrRow; ++j)
    {
      Tmp = Q(j, i);
      Q(j, i) *= CosTheta;
      Q(j, i) += Q(j, Inci) * SinTheta;
      Q(j, Inci) *= -CosTheta;
      Q(j, Inci) += Tmp * SinTheta;
    }
  --i;
  --Inci;
  for (; i > 0; --i)
    {
      // find Householder transformation
      Tmp = 1.0 / sqrt ((TmpL.UpperDiagonalElements[Inci] * TmpL.UpperDiagonalElements[Inci]) + (this->UpperDiagonalElements[i] *  this->UpperDiagonalElements[i]));
      CosTheta = -TmpL.UpperDiagonalElements[Inci] * Tmp;
      SinTheta = this->UpperDiagonalElements[i] * Tmp;
      // apply Householder transformation to the other elements

      TmpL.DiagonalElements[Inci] = TmpL.DiagonalElements[Inci] + (((this->DiagonalElements[i] - shift) + TmpL.DiagonalElements[Inci]) * SinTheta * SinTheta);
      TmpL.UpperDiagonalElements[i] *= SinTheta;
      TmpL.UpperDiagonalElements[i] += (this->DiagonalElements[i] - shift) * CosTheta;
      TmpL.DiagonalElements[i] = TmpL.UpperDiagonalElements[i] * CosTheta;
      TmpL.UpperDiagonalElements[i - 1] = this->UpperDiagonalElements[i - 1] * CosTheta;
      TmpL.UpperDiagonalElements[Inci] *= -CosTheta;
      TmpL.UpperDiagonalElements[Inci] += this->UpperDiagonalElements[i] * SinTheta;
      TmpL.UpperDiagonalElements[Inci] *= PreviousSinTheta;
      PreviousSinTheta = SinTheta;
      // apply transformation on the Q matrix
      for (int j = 0; j < this->NbrRow; ++j)
	{
	  Tmp = Q(j, i);
	  Q(j, i) *= CosTheta;
	  Q(j, i) += Q(j, Inci) * SinTheta;
	  Q(j, Inci) *= -CosTheta;
	  Q(j, Inci) += Tmp * SinTheta;
	}
      --Inci;
    }
  // evaluate last transformation
  Tmp = 1.0 / sqrt ((TmpL.UpperDiagonalElements[1] * TmpL.UpperDiagonalElements[1]) + (this->UpperDiagonalElements[0] *  this->UpperDiagonalElements[0]));
  CosTheta = -TmpL.UpperDiagonalElements[1] * Tmp;
  SinTheta = this->UpperDiagonalElements[0] * Tmp;
  TmpL.DiagonalElements[1] = TmpL.DiagonalElements[1] + (((this->DiagonalElements[0] - shift) + TmpL.DiagonalElements[1]) * SinTheta * SinTheta);
  TmpL.UpperDiagonalElements[0] *= SinTheta;
  TmpL.UpperDiagonalElements[0] += (this->DiagonalElements[0] - shift) * CosTheta;
  TmpL.DiagonalElements[0] = TmpL.UpperDiagonalElements[0] * CosTheta;
  TmpL.UpperDiagonalElements[0] *= SinTheta;
  TmpL.UpperDiagonalElements[1] *= -CosTheta;
  TmpL.UpperDiagonalElements[1] += this->UpperDiagonalElements[0] * SinTheta;
  TmpL.UpperDiagonalElements[1] *= PreviousSinTheta;
  for (int j = 0; j < this->NbrRow; ++j)
    {
      Tmp = Q(j, 0);
      Q(j, 0) *= CosTheta;
      Q(j, 0) += Q(j, 1) * SinTheta;
      Q(j, 1) *= -CosTheta;
      Q(j, 1) += Tmp * SinTheta;
      // shift diagonal values of Qt H Q
      TmpL.DiagonalElements[j] += shift;
    }
  return TmpL;
}

// find QR factorization and evaluate RQ (ie Qt H Q), shifting initial matrix diagonal elements and shifting back after evaluating RQ
//
// Q = matrix initialized with corresponding base in which unitary matrix of QR factorization has to be calculated
// shift = shift value
// return value = Qt H Q

RealTriDiagonalSymmetricMatrix RealTriDiagonalSymmetricMatrix::ConjugateQR(RealMatrix& Q, double shift)
{
  RealTriDiagonalSymmetricMatrix TmpL(this->NbrRow, true);
  double Tmp = 0.0;
  double CosTheta = 0.0;
  double SinTheta = 0.0;
  double PreviousSinTheta = 0.0;
  // perform transformation on the latest elements of the tridiagonal matrix
  // find Householder transformation and apply it to lower diagonal element
  Tmp =  1.0 / sqrt (((this->DiagonalElements[0] - shift) * (this->DiagonalElements[0] - shift)) + (this->UpperDiagonalElements[0] *  this->UpperDiagonalElements[0]));
  CosTheta = (this->DiagonalElements[0]  - shift) * Tmp;  
  SinTheta = this->UpperDiagonalElements[0] * Tmp;
  TmpL.UpperDiagonalElements[0] = SinTheta * this->UpperDiagonalElements[0] - (this->DiagonalElements[1] - shift) * CosTheta;
  TmpL.DiagonalElements[1]  =  -TmpL.UpperDiagonalElements[0] * CosTheta;
  TmpL.DiagonalElements[0] = (this->DiagonalElements[0] - shift) + ((this->DiagonalElements[1] + this->DiagonalElements[0] - 2.0 * shift) * SinTheta * SinTheta);
  TmpL.UpperDiagonalElements[1] = -this->UpperDiagonalElements[1] * CosTheta;
  PreviousSinTheta = SinTheta;
  // apply transformation on the Q matrix
  for (int j = 0; j < this->NbrRow; ++j)
    {
      Tmp = Q(j, 0);
      Q(j, 0) *= CosTheta;
      Q(j, 0) += Q(j, 1) * SinTheta;
      Q(j, 1) *= -CosTheta;
      Q(j, 1) += Tmp * SinTheta;
    }
  int i = 1;
  int Inci = 2;
  int Deci = 0;
  int ReducedSize = this->NbrRow - 2;
  for (; i < ReducedSize; ++i)
    {
      // find Householder transformation
      Tmp = 1.0 / sqrt ((this->UpperDiagonalElements[i] * this->UpperDiagonalElements[i]) + (TmpL.UpperDiagonalElements[Deci] *  TmpL.UpperDiagonalElements[Deci]));
      CosTheta = TmpL.UpperDiagonalElements[Deci] * Tmp;
      SinTheta =  this->UpperDiagonalElements[i] * Tmp;
      // apply Householder transformation to the other elements

      TmpL.UpperDiagonalElements[Deci] *= CosTheta;
      TmpL.UpperDiagonalElements[Deci] += this->UpperDiagonalElements[i] * SinTheta;
      TmpL.UpperDiagonalElements[Deci] *= PreviousSinTheta;
      TmpL.DiagonalElements[i] += (((this->DiagonalElements[Inci] - shift) + TmpL.DiagonalElements[i]) * SinTheta * SinTheta);
      TmpL.UpperDiagonalElements[i] *= SinTheta;
      TmpL.UpperDiagonalElements[i] -= (this->DiagonalElements[Inci] - shift) * CosTheta;
      TmpL.DiagonalElements[Inci] = -TmpL.UpperDiagonalElements[i] * CosTheta;
      TmpL.UpperDiagonalElements[Inci] = -this->UpperDiagonalElements[Inci] * CosTheta;
      PreviousSinTheta = SinTheta;
      // apply transformation on the Q matrix
      for (int j = 0; j < this->NbrRow; ++j)
	{
	  Tmp = Q(j, i);
	  Q(j, i) *= CosTheta;
	  Q(j, i) += Q(j, Inci) * SinTheta;
	  Q(j, Inci) *= -CosTheta;
	  Q(j, Inci) += Tmp * SinTheta;
	}
      ++Deci;
      ++Inci;
    }
  // evaluate last transformation
  // find Householder transformation
  Tmp = 1.0 / sqrt ((this->UpperDiagonalElements[i] * this->UpperDiagonalElements[i]) + (TmpL.UpperDiagonalElements[Deci] *  TmpL.UpperDiagonalElements[Deci]));
  CosTheta = TmpL.UpperDiagonalElements[Deci] * Tmp;
  SinTheta =  this->UpperDiagonalElements[i] * Tmp;
  // apply Householder transformation to the other elements
  
  TmpL.UpperDiagonalElements[Deci] *= CosTheta;
  TmpL.UpperDiagonalElements[Deci] += this->UpperDiagonalElements[i] * SinTheta;
  TmpL.UpperDiagonalElements[Deci] *= PreviousSinTheta;
  TmpL.DiagonalElements[i] += (((this->DiagonalElements[Inci] - shift) + TmpL.DiagonalElements[i]) * SinTheta * SinTheta);
  TmpL.UpperDiagonalElements[i] *= SinTheta;
  TmpL.UpperDiagonalElements[i] -= (this->DiagonalElements[Inci] - shift) * CosTheta;
  TmpL.DiagonalElements[Inci] = -TmpL.UpperDiagonalElements[i] * CosTheta;
  TmpL.UpperDiagonalElements[i] *= SinTheta;
  // apply transformation on the Q matrix
  for (int j = 0; j < this->NbrRow; ++j)
    {
      Tmp = Q(j, i);
      Q(j, i) *= CosTheta;
      Q(j, i) += Q(j, Inci) * SinTheta;
      Q(j, Inci) *= -CosTheta;
      Q(j, Inci) += Tmp * SinTheta;
      // shift diagonal values of Qt H Q
      TmpL.DiagonalElements[j] += shift;
    }
  return TmpL;
}

// apply polynomial filter assumuing shifts are exacte shift (i.e. eigenvalues of the initial matrix)
//
// Q = unitary matrix encoding the transformation
// shift = array of shift values
// nbrShift = number of shifts
// return value = filtered matrix stored as [[H 0], [0 D]] where D is a diagonal matrix with shifts as element and H a real tridiagonal matrix

RealTriDiagonalSymmetricMatrix RealTriDiagonalSymmetricMatrix::PolynomialFilterWithExactShift(RealMatrix& Q, double* shift, int nbrShift)
{
  RealTriDiagonalSymmetricMatrix TmpInitialMatrix(this->NbrRow);
  TmpInitialMatrix.Copy(*this);
  RealTriDiagonalSymmetricMatrix TmpR(this->NbrRow);
  double Tmp = 0.0;
  double CosTheta = 0.0;
  double SinTheta = 0.0;
  double PreviousSinTheta = 0.0;
  int EffectiveSize = this->NbrRow;
//      for (int i = 0; i < (this->NbrRow - 1); ++i)
//	cout << i << " " << TmpInitialMatrix.DiagonalElement(i) << " " << TmpInitialMatrix(i, i + 1) << endl;
//      cout << "---------------" << endl;
  for (int k = 0; k < nbrShift; ++k)
    {      
      int i = 0;
      int Inci = i + 1;
      // perform transformation on the latest elements of the tridiagonal matrix
      // find Householder transformation and apply it to lower diagonal element
      Tmp =  1.0 / sqrt (((TmpInitialMatrix.DiagonalElements[i] - shift[k]) * (TmpInitialMatrix.DiagonalElements[i] - shift[k])) + 
			 (TmpInitialMatrix.UpperDiagonalElements[i] *  TmpInitialMatrix.UpperDiagonalElements[i]));
      CosTheta = (TmpInitialMatrix.DiagonalElements[i]  - shift[k]) * Tmp;  
      SinTheta = TmpInitialMatrix.UpperDiagonalElements[i] * Tmp;
      TmpR.UpperDiagonalElements[i] = SinTheta * TmpInitialMatrix.UpperDiagonalElements[i] - (TmpInitialMatrix.DiagonalElements[Inci] - shift[k]) * CosTheta;
      TmpR.DiagonalElements[Inci]  =  -TmpR.UpperDiagonalElements[0] * CosTheta;
      TmpR.DiagonalElements[i] = ((TmpInitialMatrix.DiagonalElements[i] - shift[k]) + 
				  ((TmpInitialMatrix.DiagonalElements[Inci] + TmpInitialMatrix.DiagonalElements[i] - 2.0 * shift[k]) * SinTheta * SinTheta));
      TmpR.UpperDiagonalElements[Inci] = -TmpInitialMatrix.UpperDiagonalElements[Inci] * CosTheta;
      PreviousSinTheta = SinTheta;
      // apply transformation on the Q matrix
      for (int j = 0; j < this->NbrRow; ++j)
	{
	  Tmp = Q(j, i);
	  Q(j, i) *= CosTheta;
	  Q(j, i) += Q(j, Inci) * SinTheta;
	  Q(j, Inci) *= -CosTheta;
	  Q(j, Inci) += Tmp * SinTheta;
	}
      int ReducedSize = EffectiveSize - 2;
      ++i;
      ++Inci;
      int Deci = i- 1;
      for (; i < ReducedSize; ++i)
	{
	  // find Householder transformation
	  Tmp = 1.0 / sqrt ((TmpInitialMatrix.UpperDiagonalElements[i] * TmpInitialMatrix.UpperDiagonalElements[i]) + 
			    (TmpR.UpperDiagonalElements[Deci] *  TmpR.UpperDiagonalElements[Deci]));
	  CosTheta = TmpR.UpperDiagonalElements[Deci] * Tmp;
	  SinTheta =  TmpInitialMatrix.UpperDiagonalElements[i] * Tmp;
	  if (SinTheta >= 1.0)
	    {
	      CosTheta = 0.0;
	      if (SinTheta > 0.0)
		SinTheta = 1.0;
	      else
		SinTheta = -1.0;
	      TmpR.UpperDiagonalElements[Deci] = 0.0;
	    }
//	  cout << "cos = " << CosTheta << " sin = " << SinTheta << " tmp = " << Tmp << " " << ((CosTheta * CosTheta) + (SinTheta * SinTheta))<< endl;
	  // apply Householder transformation to the other elements

	  TmpR.UpperDiagonalElements[Deci] *= CosTheta;
	  TmpR.UpperDiagonalElements[Deci] += TmpInitialMatrix.UpperDiagonalElements[i] * SinTheta;
	  TmpR.UpperDiagonalElements[Deci] *= PreviousSinTheta;
	  TmpR.DiagonalElements[i] += (((TmpInitialMatrix.DiagonalElements[Inci] - shift[k]) + TmpR.DiagonalElements[i]) * SinTheta * SinTheta);
	  TmpR.UpperDiagonalElements[i] *= SinTheta;
	  TmpR.UpperDiagonalElements[i] -= (TmpInitialMatrix.DiagonalElements[Inci] - shift[k]) * CosTheta;
	  TmpR.DiagonalElements[Inci] = -TmpR.UpperDiagonalElements[i] * CosTheta;
	  TmpR.UpperDiagonalElements[Inci] = -TmpInitialMatrix.UpperDiagonalElements[Inci] * CosTheta;
	  PreviousSinTheta = SinTheta;
	  // apply transformation on the Q matrix
	  for (int j = 0; j < this->NbrRow; ++j)
	    {
	      Tmp = Q(j, i);
	      Q(j, i) *= CosTheta;
	      Q(j, i) += Q(j, Inci) * SinTheta;
	      Q(j, Inci) *= -CosTheta;
	      Q(j, Inci) += Tmp * SinTheta;
	    }
	  ++Deci;
	  ++Inci;
//      for (int i = 0; i < (this->NbrRow - 1); ++i)
//	cout << i << " " << TmpR.DiagonalElement(i) << " " << TmpR(i, i + 1) << endl;
//      cout << "---------------" << endl;
	}
      // evaluate last transformation
      // find Householder transformation
      Tmp = 1.0 / sqrt ((TmpInitialMatrix.UpperDiagonalElements[i] * TmpInitialMatrix.UpperDiagonalElements[i]) + 
			(TmpR.UpperDiagonalElements[Deci] *  TmpR.UpperDiagonalElements[Deci]));
      CosTheta = TmpR.UpperDiagonalElements[Deci] * Tmp;
      SinTheta =  TmpInitialMatrix.UpperDiagonalElements[i] * Tmp;
      // apply Householder transformation to the other elements
  
      TmpR.UpperDiagonalElements[Deci] *= CosTheta;
      TmpR.UpperDiagonalElements[Deci] += TmpInitialMatrix.UpperDiagonalElements[i] * SinTheta;
      TmpR.UpperDiagonalElements[Deci] *= PreviousSinTheta;
      TmpR.DiagonalElements[i] += (((TmpInitialMatrix.DiagonalElements[Inci] - shift[k]) + TmpR.DiagonalElements[i]) * SinTheta * SinTheta);
      TmpR.UpperDiagonalElements[i] *= SinTheta;
      TmpR.UpperDiagonalElements[i] -= (TmpInitialMatrix.DiagonalElements[Inci] - shift[k]) * CosTheta;
      TmpR.DiagonalElements[Inci] = -TmpR.UpperDiagonalElements[i] * CosTheta;
      TmpR.UpperDiagonalElements[i] *= SinTheta;
//      for (int i = 0; i < (this->NbrRow - 1); ++i)
//	cout << i << " " << TmpR.DiagonalElement(i) << " " << TmpR(i, i + 1) << endl;
//      cout << "---------------" << endl;
      // apply transformation on the Q matrix
      for (int j = 0; j < this->NbrRow; ++j)
	{
	  Tmp = Q(j, i);
	  Q(j, i) *= CosTheta;
	  Q(j, i) += Q(j, Inci) * SinTheta;
	  Q(j, Inci) *= -CosTheta;
	  Q(j, Inci) += Tmp * SinTheta;
	}
      // shift diagonal values of Qt H Q
      for (int j = 0; j < EffectiveSize; ++j)
	{
	  TmpR.DiagonalElements[j] += shift[k];
	}
//      cout << TmpR << endl;
//      for (int i = 0; i < (this->NbrRow - 1); ++i)
//	cout << i << " " << TmpR.DiagonalElement(i) << " " << TmpR(i, i + 1) << endl;
//      cout << "---------------" << endl;
      RealTriDiagonalSymmetricMatrix TmpM (TmpR);
      TmpR = TmpInitialMatrix;
      TmpInitialMatrix = TmpM;
      --EffectiveSize;
      TmpR.DiagonalElements[EffectiveSize] = TmpInitialMatrix.DiagonalElements[EffectiveSize];
      TmpR.UpperDiagonalElements[EffectiveSize - 1] = TmpInitialMatrix.UpperDiagonalElements[EffectiveSize - 1];
      return TmpInitialMatrix;
    }
  return TmpInitialMatrix;
}

// evaluate a normalized eigenvector for a given eigenvalue (supposing the eigenvalue is non-degenerate)
//
// eigenvalue = eigenvalue to use
// eigenvector = vector where the eigenvector has to be stored
// return value = reference on eigenvector

RealVector& RealTriDiagonalSymmetricMatrix::Eigenvector(double eigenvalue, RealVector& eigenvector)
{
  double Norm = 1.0;
  eigenvector.Components[0] = 1.0;
  if (this->UpperDiagonalElements[0] == 0.0)
    eigenvector.Components[1] = 0.0;
  else
    eigenvector.Components[1] = (eigenvalue - this->DiagonalElements[0]) / this->UpperDiagonalElements[0];
  int ReducedNbrRow = this->NbrRow - 1;
  for (int i = 1; i < ReducedNbrRow; i++)
    {
      if (this->UpperDiagonalElements[i] == 0.0)
	eigenvector.Components[i + 1] = 0.0;
      else
	{
	  eigenvector.Components[i + 1] = (((eigenvalue - this->DiagonalElements[i]) * eigenvector.Components[i]
					    - this->UpperDiagonalElements[i - 1] * eigenvector.Components[i - 1])
					   / this->UpperDiagonalElements[i]);
	  Norm += eigenvector.Components[i + 1] * eigenvector.Components[i + 1];
	}
    }      
/*  eigenvector.Components[ReducedNbrRow + 1] = ((eigenvector.Components[ReducedNbrRow] * this->UpperDiagonalElements[ReducedNbrRow]) /
					       (eigenvalue - this->DiagonalElements[ReducedNbrRow + 1]));
  Norm += eigenvector.Components[ReducedNbrRow + 1] * eigenvector.Components[ReducedNbrRow + 1];*/
  Norm = 1.0 / sqrt (Norm);
  for (int i = 0; i < this->NbrRow; i++)
    eigenvector.Components[i] *= Norm;
  return eigenvector;
}

// evaluate a normalized eigenvector for a given eigenvalue (supposing the eigenvalue is non-degenerate)
//
// eigenvalue = eigenvalue to use
// eigenvector = vector where the eigenvector has to be stored
// return value = reference on eigenvector

ComplexVector& RealTriDiagonalSymmetricMatrix::Eigenvector(double eigenvalue, ComplexVector& eigenvector)
{
  double Norm = 1.0;
  eigenvector.RealComponents[0] = 1.0;
  eigenvector.ImaginaryComponents[0] = 0.0;
  eigenvector.ImaginaryComponents[1] = 0.0;
  if (this->UpperDiagonalElements[0] == 0.0)
    eigenvector.RealComponents[1] = 0.0;
  else
    eigenvector.RealComponents[1] = (eigenvalue - this->DiagonalElements[0])/ this->UpperDiagonalElements[0];
  for (int i = 2; i < this->NbrRow; i++)
    {
      if (this->UpperDiagonalElements[i - 1] == 0.0)
	eigenvector.RealComponents[i] = 0.0;
      else
	{
	  eigenvector.RealComponents[i] = - (((this->DiagonalElements[i - 1] - eigenvalue) * eigenvector.RealComponents[i - 1]
					      + this->UpperDiagonalElements[i - 2] * eigenvector.RealComponents[i - 2])
					     / this->UpperDiagonalElements[i - 1]);
	  Norm += eigenvector.RealComponents[i] * eigenvector.RealComponents[i];
	}
      eigenvector.ImaginaryComponents[i] = 0.0;
    }      
  Norm = 1.0 / sqrt (Norm);
  for (int i = 0; i < eigenvector.Dimension; i++)
    eigenvector.RealComponents[i] *= Norm;
  return eigenvector;
}

// Sort Matrix such that diagnonal elements are sort in increasing order (offdiagonal elements left unchanged)
//
// return value = reference on current Matrix

RealTriDiagonalSymmetricMatrix& RealTriDiagonalSymmetricMatrix::SortMatrixUpOrder()
{
  int ReducedDim = this->NbrColumn - 1;
  double tmp;
  for (int i = 0; i < ReducedDim; i++)
    for (int j = 0; j < (ReducedDim - i); j++)
      if (this->DiagonalElements[j] > this->DiagonalElements[j + 1])
	{
	  tmp = this->DiagonalElements[j];
	  this->DiagonalElements[j] = this->DiagonalElements[j + 1];
	  this->DiagonalElements[j + 1] = tmp;
	}
  return *this;
}

// Sort Matrix such that diagnonal elements are sort in increasing order (offdiagonal elements left unchanged) 
// and apply corresponding transformation to column of a given real matrix 
//
// matrix = matrix on which transformation has to be applied
// return value = reference on current Matrix

RealTriDiagonalSymmetricMatrix& RealTriDiagonalSymmetricMatrix::SortMatrixUpOrder(RealMatrix& matrix)
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
      TmpV = matrix.Columns[i];
      matrix.Columns[i] = matrix.Columns[MinPos];
      matrix.Columns[MinPos] = TmpV;
    }
  return *this;
}

// Sort Matrix such that diagnonal elements are sort in decreasing order (offdiagonal elements left unchanged)
//
// return value = reference on current Matrix

RealTriDiagonalSymmetricMatrix& RealTriDiagonalSymmetricMatrix::SortMatrixDownOrder()
{
  int ReducedDim = this->NbrColumn - 1;
  double tmp;
  for (int i = 0; i < ReducedDim; i++)
    for (int j = 0; j < (ReducedDim - i); j++)
      if (this->DiagonalElements[j] < this->DiagonalElements[j + 1])
	{
	  tmp = this->DiagonalElements[j];
	  this->DiagonalElements[j] = this->DiagonalElements[j + 1];
	  this->DiagonalElements[j + 1] = tmp;
	}
  return *this;
}

// Sort Matrix such that diagnonal elements are sort in decreasing order (offdiagonal elements left unchanged) 
// and apply corresponding transformation to column of a given real matrix 
//
// matrix = matrix on which transformation has to be applied
// return value = reference on current Matrix

RealTriDiagonalSymmetricMatrix& RealTriDiagonalSymmetricMatrix::SortMatrixDownOrder(RealMatrix& matrix)
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
      TmpV = matrix.Columns[i];
      matrix.Columns[i] = matrix.Columns[MinPos];
      matrix.Columns[MinPos] = TmpV;
    }
  return *this;
}

// Sort Matrix such that diagnonal elements are sort in increasing order and apply corresponding transformation to column
// of a given complex matrix (offdiagonal elements left unchanged)
//
// Q = matrix on which transformation has to be applied

// return value = reference on current Matrix

RealTriDiagonalSymmetricMatrix& RealTriDiagonalSymmetricMatrix::SortMatrix(ComplexMatrix& Q)
{
  int ReducedDim = this->NbrColumn - 1;
  ComplexVector TmpV;
  double tmp;
  for (int i = 0; i < ReducedDim; i++)
    for (int j = 0; j < (ReducedDim - i); j++)
      if (this->DiagonalElements[j] > this->DiagonalElements[j + 1])
	{
	  tmp = this->DiagonalElements[j];
	  this->DiagonalElements[j] = this->DiagonalElements[j + 1];
	  this->DiagonalElements[j + 1] = tmp;
	  TmpV = Q.Columns[j];
	  Q.Columns[j] = Q.Columns[j + 1];
	  Q.Columns[j + 1] = TmpV;
	}
  return *this;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const RealTriDiagonalSymmetricMatrix& P)
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      int j = 0;
      for (; j < (i-1); j++)
	Str << "0    ";
      if (i > 0)
	{
	  Str << P.UpperDiagonalElements[(i - 1)] << "    ";
	  j++;
	}
      Str << P.DiagonalElements[i] << "    ";
      j++;
      if (i < (P.NbrRow -1))
	{
	  Str << P.UpperDiagonalElements[i] << "    ";
	  j++;
	}
      for (; j < P.NbrColumn; j++)
	Str << "0    ";
      Str << endl;
    }
/*  for (int i = 0; i < (P.NbrColumn - 1); i++)
    if (i > 0)
      {
	Str << P.DiagonalElements[i] << "    " << P.UpperDiagonalElements[i] << endl;
      }
  Str << P.DiagonalElements[P.NbrColumn - 1] << endl;*/
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// P = matrix to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, const RealTriDiagonalSymmetricMatrix& P)
{
  Str << "{";
  for (int i = 0; i < (P.NbrRow - 1); i++)
    {
      Str << "{";
      int j = 0;
      for (; j < (i-1); j++)
	Str << "0,";
      if (i > 0)
	{
	  Str << P.UpperDiagonalElements[(i - 1)] << ",";
	  j++;
	}
      Str << P.DiagonalElements[i] << ",";
      j++;
      if (i < (P.NbrRow -1))
	{
	  Str << P.UpperDiagonalElements[i];
	  if (j != (P.NbrColumn - 1))
	    Str << ",";
	  j++;
	}
      for (; j < (P.NbrColumn - 1); j++)
	Str << "0,";
      if (j == (P.NbrColumn - 1))
	Str << "0";
      Str << "},";
    }
  Str << "{";
  int j = 0;
  for (; j < (P.NbrRow-2); j++)
    Str << "0,";
  if (P.NbrRow > 0)
    {
      Str << P.UpperDiagonalElements[P.NbrRow - 2] << ",";
    }
  Str << P.DiagonalElements[P.NbrRow - 1];
  Str << "}}";
  return Str;
}

