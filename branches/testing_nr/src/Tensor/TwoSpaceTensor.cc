////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of tensor acting on two consecutive spaces            //
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


#include "Tensor/TwoSpaceTensor.h"
#include "Tensor/OneSpaceTensor.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "HilbertSpace/SpaceDecomposition.h"


using std::endl;


// default constructor
//

TwoSpaceTensor::TwoSpaceTensor() 
{
  this->Structure = 0;
  this->ElementaryMatrix = 0;
  this->FirstTargetSpace = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = 0;
  this->TrueNbrColumn = 0;
  this->MatrixType = 0;
}

// constructor from standard datas
//
// struture = reference on tensor product structure
// elementaryMatrix = matrix to use
// firstTargetSpace = first space on which matrix acts

TwoSpaceTensor::TwoSpaceTensor(AbstractTensorProductStructure* structure, Matrix& elementaryMatrix, int firstTargetSpace) 
{
  this->Structure = structure;
  this->FirstTargetSpace = firstTargetSpace;
  this->ElementaryMatrix = elementaryMatrix.Clone();
  this->NbrRow = this->ElementaryMatrix->GetNbrRow();
  this->NbrColumn = this->NbrRow;
  this->TrueNbrRow = this->Structure->GetTotalDimension();
  this->TrueNbrColumn = this->TrueNbrRow;
  this->MatrixType = this->ElementaryMatrix->GetMatrixType();
}

// constructor from the tensor product of two matrices
//
// struture = reference on tensor product structure
// m1 = pointer to the matrix to use for first target space
// m2 = pointer to the matrix to use for second target space
// firstTargetSpace = first space on which matrix acts
// coef = global multiplicative factor

TwoSpaceTensor::TwoSpaceTensor(AbstractTensorProductStructure* structure, Matrix* m1, Matrix* m2, 
			       int firstTargetSpace, double coef)
{
  this->Structure = structure;
  int Dim1 = this->Structure->GetDimension(firstTargetSpace);
  int Dim2 = this->Structure->GetDimension(firstTargetSpace + 1);
  if ((m1->GetMatrixType() != m2->GetMatrixType()) || (m1->GetNbrRow() != Dim1) || 
      (m2->GetNbrRow() != Dim2)) 
    {
      this->Structure = 0;
      this->ElementaryMatrix = 0;
      this->FirstTargetSpace = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = 0;
    }
  else
    {
      this->FirstTargetSpace = firstTargetSpace;
      this->MatrixType = m1->GetMatrixType();
      this->NbrRow = Dim1 * Dim2;
      this->NbrColumn = Dim1 * Dim2;
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;
      switch (m1->GetMatrixType())
	{
	case (Matrix::RealElements | Matrix::Symmetric):
	  {    
	    RealSymmetricMatrix* TmpM = new RealSymmetricMatrix (Dim1 * Dim2);
	    for (int i2 = 0; i2 < Dim2; i2++)
	      for (int j2 = 0; j2 <= i2; j2++)
		{
		  double Coef = coef * (*((RealSymmetricMatrix*) m2))(i2, j2);
		  for (int i1 = 0; i1 < Dim1; i1++)
		    {
		      int Pos =i1 + Dim1 * i2;
		      int Pos2 = Dim1 * j2;
		      int Lim;
		      if (j2 == i2)
			Lim = i1;
		      else
			Lim = Dim1 - 1;
		      for (int j1 = 0; j1 <= Lim; j1++)
			{
			  (*TmpM)(Pos, Pos2++) = (*((RealSymmetricMatrix*) m1))(i1, j1) * Coef;
			}
		    }
		  this->ElementaryMatrix = (Matrix*) TmpM;
		}
	  }
	  break;
	case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
	  {
	    RealDiagonalMatrix* TmpM = new RealDiagonalMatrix (Dim1 * Dim2);
	    int Pos = 0;
	    for (int i2 = 0; i2 < Dim2; i2++)
	      {
		double Coef = coef * (*((RealDiagonalMatrix*) m2))[i2];
		for (int i1 = 0; i1 < Dim1; i1++)
		  {
		    (*TmpM) [Pos++] = (Coef * (*((RealDiagonalMatrix*) m1))[i1]);
		  }
	      }
	    this->ElementaryMatrix = (Matrix*) TmpM;
	  }
	  break;
	case (Matrix::RealElements | Matrix::Antisymmetric):
	  {
	    this->MatrixType = Matrix::RealElements | Matrix::Symmetric;
	    RealSymmetricMatrix* TmpM = new RealSymmetricMatrix (Dim1 * Dim2);
	    for (int i2 = 0; i2 < Dim2; i2++)
	      for (int j2 = 0; j2 < i2; j2++)
		{
		  double Coef = (*((RealAntisymmetricMatrix*) m2))(j2, i2) * coef;
		  for (int i1 = 0; i1 < Dim1; i1++)
		    {
		      int Pos2 = Dim1 * j2;
		      int Pos =i1 + Dim1 * i2;
		      int j1 = 0;
		      for (; j1 < i1; j1++)
			(*TmpM)(Pos, Pos2++) = ((*((RealAntisymmetricMatrix*) m1))(j1, i1)
						* Coef);
		      j1++;
		      Pos2++;
		      for (; j1 < Dim1; j1++)
			(*TmpM)(Pos, Pos2++) = -((*((RealAntisymmetricMatrix*) m1))(i1, j1)
						 * Coef);
		    }
		}
	    this->ElementaryMatrix = (Matrix*) TmpM;
	  }
	  break;
	default:
	  this->ElementaryMatrix = new Matrix ();
	  break;
	}
    }
}

// constructor from the tensor product of two matrices with a given space decomposition
//
// struture = reference on tensor product structure
// m1 = reference on matrix to use for first target space
// m2 = reference on matrix to use for second target space
// decomposition1 = first space decomposition
// decomposition2 = second space decomposition
// firstTargetSpace = first space on which matrix acts
// coef = global multiplicative factor

TwoSpaceTensor::TwoSpaceTensor(AbstractTensorProductStructure* structure, Matrix& m1, 
			       Matrix& m2, SpaceDecomposition& decomposition1,
			       SpaceDecomposition& decomposition2, int firstTargetSpace, double coef)
{
  if (m2.GetMatrixType() != m1.GetMatrixType())
    *this = TwoSpaceTensor();
  else
    {
      switch (m1.GetMatrixType())
	{
	case (Matrix::RealElements | Matrix::Symmetric):
	  {
	    *this = TwoSpaceTensor(structure, (RealSymmetricMatrix&) m1, 
				   (RealSymmetricMatrix&) m2, decomposition1, decomposition2,
				   firstTargetSpace, coef);
	  }
	  break;
	default:
	  *this = TwoSpaceTensor();
	}
    }
}

// constructor from the tensor product of two matrices with a given space decomposition
//
// struture = reference on tensor product structure
// m1 = reference on matrix to use for first target space
// m2 = reference on matrix to use for second target space
// decomposition1 = first space decomposition
// decomposition2 = second space decomposition
// firstTargetSpace = first space on which matrix acts
// coef = global multiplicative factor

TwoSpaceTensor::TwoSpaceTensor(AbstractTensorProductStructure* structure, RealSymmetricMatrix& m1, 
			       RealSymmetricMatrix& m2, SpaceDecomposition& decomposition1,
			       SpaceDecomposition& decomposition2, int firstTargetSpace, double coef)
{
  this->NbrRow = (decomposition1.GetSubspaceDimension(decomposition1.GetNbrSubspace() - 1) * 
		  decomposition2.GetSubspaceDimension(decomposition1.GetNbrSubspace() - 1));
  int* SubspacePosition = new int [decomposition1.GetNbrSubspace()];
  SubspacePosition[0] = 0;
  for (int i = 1; i < decomposition1.GetNbrSubspace(); i++)
    {
      SubspacePosition[i] = (decomposition1.GetSubspaceDimension(i - 1) * 
			     decomposition2.GetSubspaceDimension(i - 1));
      this->NbrRow += SubspacePosition[i];
      SubspacePosition[i] += SubspacePosition[i - 1];
    }
  this->NbrColumn = this->NbrRow;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::Symmetric;
  this->Structure = structure;
  this->FirstTargetSpace = firstTargetSpace;
  RealSymmetricMatrix* TmpM = new RealSymmetricMatrix (this->NbrColumn);
  double Coef;
  int I = 0;
  int J = 0;
  int Lim = 0;
  for (int b1 = 0; b1 < decomposition1.GetNbrSubspace(); b1++)
    {
      for (int i2 = 0; i2 < decomposition2.GetSubspaceDimension(b1); i2 ++)
	{
	  I = SubspacePosition[b1] + i2 * decomposition1.GetSubspaceDimension(b1);
	  for (int j2 = i2; j2 < decomposition2.GetSubspaceDimension(b1); j2 ++)
	    {
	      Coef = m2 (i2 + decomposition2.GetSubspacePosition(b1), 
			 j2 + decomposition2.GetSubspacePosition(b1)) * coef;
	      J = SubspacePosition[b1] + j2 * decomposition1.GetSubspaceDimension(b1);
	      for (int i1 = 0; i1 < decomposition1.GetSubspaceDimension(b1); i1 ++)
		{
		  if (I == J)
		    Lim = i1;
		  else
		    Lim = 0;
		  for (int j1 = Lim; j1 < decomposition1.GetSubspaceDimension(b1); j1 ++)
		    {
		      (*TmpM)(I + i1, J + j1) = m1 (i1 + decomposition1.GetSubspacePosition(b1), 
						    j1 + decomposition1.GetSubspacePosition(b1)) * Coef;
		    }
		}
	    }
	}
      for (int b2 = b1 + 1; b2 < decomposition1.GetNbrSubspace(); b2++)
	for (int i2 = 0; i2 < decomposition2.GetSubspaceDimension(b1); i2 ++)
	  {
	    I = SubspacePosition[b1] + i2 * decomposition1.GetSubspaceDimension(b1);
	    for (int j2 = 0; j2 < decomposition2.GetSubspaceDimension(b2); j2 ++)
	      {
		Coef = m2 (i2 + decomposition2.GetSubspacePosition(b1), 
			   j2 + decomposition2.GetSubspacePosition(b2)) * coef;
		J = SubspacePosition[b2] + j2 * decomposition1.GetSubspaceDimension(b2);
		for (int i1 = 0; i1 < decomposition1.GetSubspaceDimension(b1); i1 ++)
		  for (int j1 = 0; j1 < decomposition1.GetSubspaceDimension(b2); j1 ++)
		    {
		      (*TmpM)(I + i1, J + j1) = m1 (i1 + decomposition1.GetSubspacePosition(b1), 
						    j1 + decomposition1.GetSubspacePosition(b2)) * Coef;
		    }
	      }
	  }
    }
  this->ElementaryMatrix = (Matrix*) TmpM;
}

// constructor from a given matrix assuming the other is identity (with a space index permutation)
//
// struture = reference on tensor product structure
// m = reference on matrix to use
// permutation = two dimensional array describing permutation
// firstTargetSpace = first space on which matrix acts
// coef = global multiplicative factor

TwoSpaceTensor::TwoSpaceTensor(AbstractTensorProductStructure* structure, Matrix& m, 
			       int** permutation, int firstTargetSpace, double coef)
{
      this->Structure = structure;
      int Dim1 = this->Structure->GetDimension(firstTargetSpace);
      int Dim2 = this->Structure->GetDimension(firstTargetSpace + 1);
      this->NbrRow = Dim1 * Dim2;
      this->NbrColumn = this->NbrRow;
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;
      this->FirstTargetSpace = firstTargetSpace;
      this->MatrixType = m.GetMatrixType();
      switch (this->MatrixType)
	{
	case (Matrix::RealElements | Matrix::Symmetric):
	  {
	    RealSymmetricMatrix* TmpM = new RealSymmetricMatrix (this->NbrColumn, true);
	    for (int j1 = 0; j1 < this->Structure->GetDimension(firstTargetSpace); j1++)
	      for (int j2 = 0; j2 < this->Structure->GetDimension(firstTargetSpace + 1); j2++)
		for (int i2 = j2; i2 < this->Structure->GetDimension(firstTargetSpace + 1); i2++)
		  (*TmpM)(permutation[j1][i2], permutation[j1][j2]) = ((RealSymmetricMatrix&) m)(i2, j2) * coef;
	    this->ElementaryMatrix = (Matrix*) TmpM;
	  }
	  break;
	case (Matrix::RealElements | Matrix::Antisymmetric):
	  {
	    int PosI;
	    int PosJ;
	    RealAntisymmetricMatrix* TmpM = new RealAntisymmetricMatrix (this->NbrColumn, true);
	    for (int j1 = 0; j1 < this->Structure->GetDimension(firstTargetSpace); j1++)
	      for (int i2 = 0 ; i2 < this->Structure->GetDimension(firstTargetSpace + 1); i2++)
		{
		  PosI = permutation[j1][i2];		  
		  for (int j2 = i2 + 1; j2 < this->Structure->GetDimension(firstTargetSpace + 1); j2++)
		    {
		      PosJ = permutation[j1][j2];
		      if (PosI < PosJ)
			(*TmpM)(PosI, PosJ) = ((RealAntisymmetricMatrix&) m)(i2, j2) * coef;
		      else
			(*TmpM)(PosI, PosJ) = -(((RealAntisymmetricMatrix&) m)(i2, j2) * coef);	
		    }
		}		
	    this->ElementaryMatrix = (Matrix*) TmpM;
	  }
	  break;
	}
}

// constructor from a given matrix asociated to a given space 
// and array containing positions where elements are to be put (assuming identity on other spaces)
//
// struture = reference on tensor product structure
// m = reference on matrix to use
// totalSpaceDimension = dimension of the total space
// arraySize = size of used array
// subspaceSize = indicate size of each subspace associated to a fixed coordinate for spaces where identity acts
// spaceIndex = array containing index in space where matrix acts for each fixed coordinate
// globalSpaceIndex = array containing global index in space where matrix acts for each fixed coordinate
// targetSpace = target space index
// coef = global multiplicative factor

TwoSpaceTensor::TwoSpaceTensor(AbstractTensorProductStructure* structure, Matrix& m, 
			       int totalSpaceDimension, int arraySize, int* subspaceSize,
			       int** spaceIndex, int** globalSpaceIndex, int targetSpace, double coef)
{
  this->Structure = structure;
  this->NbrRow =totalSpaceDimension;
  this->NbrColumn = totalSpaceDimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->FirstTargetSpace = targetSpace;
  this->MatrixType = m.GetMatrixType();
  switch (this->MatrixType)
    {
    case (Matrix::RealElements | Matrix::Symmetric):
      {
	RealSymmetricMatrix* TmpM = new RealSymmetricMatrix (this->NbrColumn, true);
	int TmpDim;
	int TmpIndex;
	int TmpGlobalIndex;
	int* TmpSpaceIndex;
	int* TmpGlobalSpaceIndex;
	for (int i = 0; i < arraySize; i++)
	  {
	    TmpDim = subspaceSize[i];		
	    if (TmpDim > 0)
	      {
		TmpSpaceIndex = spaceIndex[i];
		TmpGlobalSpaceIndex = globalSpaceIndex[i];
		for (int j2 = 0; j2 < TmpDim; j2++)
		  {
		    TmpIndex = TmpSpaceIndex[j2];
		    TmpGlobalIndex = TmpGlobalSpaceIndex[j2];
		    for (int i2 = j2; i2 < TmpDim; i2++)
		      (*TmpM)(TmpGlobalSpaceIndex[i2], TmpGlobalIndex) = 
			((RealSymmetricMatrix&) m)(TmpSpaceIndex[i2], TmpIndex) * coef;
		  }
	      }
	  }
	this->ElementaryMatrix = (Matrix*) TmpM;
      }
    break;
    case (Matrix::RealElements | Matrix::Antisymmetric):
      {
	RealAntisymmetricMatrix* TmpM = new RealAntisymmetricMatrix (this->NbrColumn, true);
	int TmpDim;
	int TmpIndex;
	int TmpGlobalIndex;
	int* TmpSpaceIndex;
	int* TmpGlobalSpaceIndex;
	int Pos;
	for (int i = 0; i < arraySize; i++)
	  {
	    TmpDim = subspaceSize[i];		
	    if (TmpDim > 0)
	      {
		TmpSpaceIndex = spaceIndex[i];
		TmpGlobalSpaceIndex = globalSpaceIndex[i];
		for (int j2 = 0; j2 < TmpDim; j2++)
		  {
		    TmpIndex = TmpSpaceIndex[j2];
		    TmpGlobalIndex = TmpGlobalSpaceIndex[j2];
		    for (int i2 = j2 + 1; i2 < TmpDim; i2++)
		      {
			Pos = TmpGlobalSpaceIndex[i2];
			if (Pos > TmpGlobalIndex)
			  (*TmpM)(TmpGlobalIndex, TmpGlobalSpaceIndex[i2]) = 
			    ((RealAntisymmetricMatrix&) m)(TmpIndex, TmpSpaceIndex[i2]) * coef;
			else
			  (*TmpM)(TmpGlobalSpaceIndex[i2], TmpGlobalIndex) = 
			    -((RealAntisymmetricMatrix&) m)(TmpIndex, TmpSpaceIndex[i2]) * coef;
		      }
		  }
	      }
	  }
	this->ElementaryMatrix = (Matrix*) TmpM;
      }
    break;
    }
}

// constructor from the product of two one space tensors
//
// t1 = first one space tensor
// t2 = second one space tensor
// coef = global multiplicative factor

TwoSpaceTensor::TwoSpaceTensor(const OneSpaceTensor& t1, const OneSpaceTensor& t2, double coef)
{
  if ((t1.Structure != t2.Structure) || ((t1.TargetSpace != t2.TargetSpace + 1) && (t1.TargetSpace != t2.TargetSpace - 1)) 
      || (t1.MatrixType != t2.MatrixType))
    {
      this->Structure = 0;
      this->ElementaryMatrix = 0;
      this->FirstTargetSpace = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = 0;
    }
  else
    {
      this->Structure = t1.Structure;
      this->NbrRow = t1.NbrRow;
      this->NbrColumn = t1.NbrColumn;
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;
      if (t1.TargetSpace < t2.TargetSpace)
	this->FirstTargetSpace = t1.TargetSpace;
      else
	this->FirstTargetSpace = t2.TargetSpace;
      this->MatrixType = t1.MatrixType;
      switch (this->MatrixType)
	{
	case (Matrix::RealElements | Matrix::Symmetric):
	  {
	    int Dim1 = this->Structure->GetDimension(this->FirstTargetSpace);
	    int Dim2 = this->Structure->GetDimension(this->FirstTargetSpace + 1);
	    RealSymmetricMatrix* TmpM = new RealSymmetricMatrix (Dim1 * Dim2);
	    for (int i2 = 0; i2 < Dim2; i2++)
	      for (int j2 = 0; j2 <= i2; j2++)
		{
		  double Coef = coef * (*((RealSymmetricMatrix*) t2.ElementaryMatrix))(i2, j2);
		  for (int i1 = 0; i1 < Dim1; i1++)
		    {
		      int Pos =i1 + Dim1 * i2;
		      int Pos2 = Dim1 * j2;
		      int Lim;
		      if (j2 == i2)
			Lim = i1;
		      else
			Lim = Dim1 - 1;
		      for (int j1 = 0; j1 <= Lim; j1++)
			(*TmpM)(Pos, Pos2++) += (*((RealSymmetricMatrix*) t1.ElementaryMatrix))(i1, j1) * Coef;
		    }
		  this->ElementaryMatrix = (Matrix*) TmpM;
		}
	  }
	  break;
	case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
	  {
	    int Dim1 = this->Structure->GetDimension(this->FirstTargetSpace);
	    int Dim2 = this->Structure->GetDimension(this->FirstTargetSpace + 1);
	    RealDiagonalMatrix* TmpM = new RealDiagonalMatrix (Dim1 * Dim2);
	    int Pos = 0;
	    for (int i2 = 0; i2 < Dim2; i2++)
	      {
		double Coef = coef * (*((RealDiagonalMatrix*) t2.ElementaryMatrix))[i2];
		for (int i1 = 0; i1 < Dim1; i1++)
		  {
		    (*TmpM) [Pos++] = (Coef * (*((RealDiagonalMatrix*) t1.ElementaryMatrix))[i1]);
		  }
	      }
	  }
	  break;
	case (Matrix::RealElements | Matrix::Antisymmetric):
	  {
	    this->MatrixType = Matrix::RealElements | Matrix::Symmetric;
	    int Dim1 = this->Structure->GetDimension(this->FirstTargetSpace);
	    int Dim2 = this->Structure->GetDimension(this->FirstTargetSpace + 1);
	    RealSymmetricMatrix* TmpM = new RealSymmetricMatrix (Dim1 * Dim2);
	    for (int i2 = 0; i2 < Dim2; i2++)
	      for (int j2 = 0; j2 < i2; j2++)
		{
		  double Coef = (*((RealAntisymmetricMatrix*) t2.ElementaryMatrix))(j2, i2) * coef;
		  for (int i1 = 0; i1 < Dim1; i1++)
		    {
		      int Pos2 = Dim1 * j2;
		      int Pos =i1 + Dim1 * i2;
		      int j1 = 0;
		      for (; j1 < i1; j1++)
			(*TmpM)(Pos, Pos2++) = ((*((RealAntisymmetricMatrix*) t1.ElementaryMatrix))(j1, i1)
						* Coef);
		      j1++;
		      Pos2++;
		      for (; j1 < Dim1; j1++)
			(*TmpM)(Pos, Pos2++) = -((*((RealAntisymmetricMatrix*) t1.ElementaryMatrix))(i1, j1)
						 * Coef);
		    }
		}
	    this->ElementaryMatrix = (Matrix*) TmpM;
	  }
	  break;
	default:
	  this->ElementaryMatrix = new Matrix ();
	  break;
	}
    }
}
  
// copy constructor (without duplicating datas)
//
// M = matrix to copy

TwoSpaceTensor::TwoSpaceTensor(const TwoSpaceTensor& T) 
{
  this->Structure = T.Structure;
  this->FirstTargetSpace = T.FirstTargetSpace;
  this->ElementaryMatrix = T.ElementaryMatrix->Clone();
  this->NbrRow = T.NbrRow;
  this->NbrColumn = T.NbrColumn;
  this->TrueNbrRow = T.TrueNbrRow;
  this->TrueNbrColumn = T.TrueNbrColumn;
  this->MatrixType = T.MatrixType;
}

// destructor
//

TwoSpaceTensor::~TwoSpaceTensor() 
{
  if (this->ElementaryMatrix != 0)
    delete this->ElementaryMatrix;
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

TwoSpaceTensor& TwoSpaceTensor::operator = (const TwoSpaceTensor& T) 
{
  this->Structure = T.Structure;
  this->FirstTargetSpace = T.FirstTargetSpace;
  this->ElementaryMatrix = T.ElementaryMatrix->Clone();
  this->NbrRow = T.NbrRow;
  this->NbrColumn = T.NbrColumn;
  this->TrueNbrRow = T.TrueNbrRow;
  this->TrueNbrColumn = T.TrueNbrColumn;
  this->MatrixType = T.MatrixType;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* TwoSpaceTensor::Clone ()
{
  return (Matrix*) new TwoSpaceTensor (*this);
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void TwoSpaceTensor::SetMatrixElement(int i, int j, double x) 
{
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void TwoSpaceTensor::SetMatrixElement(int i, int j, const Complex& x) 
{
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void TwoSpaceTensor::AddToMatrixElement(int i, int j, double x) 
{
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void TwoSpaceTensor::AddToMatrixElement(int i, int j, const Complex& x) 
{
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void TwoSpaceTensor::Resize (int nbrRow, int nbrColumn) 
{
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void TwoSpaceTensor::ResizeAndClean (int nbrRow, int nbrColumn) 
{
}

// Resize tensor
//
// structure = new product tensor structure

void TwoSpaceTensor::Resize (const TensorProductStructure& structure) 
{
}

// Resize tensor and set to zero all components that have been added
//
// structure = new product tensor structure

void TwoSpaceTensor::ResizeAndClean (const TensorProductStructure& structure) 
{
}

// Add product of two one space tensors to the current tensor
//
// t1 = first one space tensor
// t2 = second one space tensor
// return value = reference on current tensor

TwoSpaceTensor& TwoSpaceTensor::AddProduct (const OneSpaceTensor& t1, const OneSpaceTensor& t2)
{
  if (((*(t1.Structure)) != (*(t2.Structure))) || ((*(this->Structure)) != (*(t1.Structure))) ||
      ((t1.TargetSpace != t2.TargetSpace + 1) && (t1.TargetSpace != t2.TargetSpace - 1)) ||
      ((t1.TargetSpace != t1.TargetSpace) && (t1.TargetSpace != t2.TargetSpace))|| (t1.MatrixType != t2.MatrixType))
    return *this;
  return *this;
}

// add tensor product of two matrices to the cuurent tensor (assuming identical space tensor pattern)
//
// m1 = pointer to the matrix associated to first target space
// m2 = pointer to the matrix associated to first target space
// coef = global multiplicative factor
// return value = reference to current tensor

TwoSpaceTensor& TwoSpaceTensor::AddTensorProductMatrices (Matrix* m1, Matrix* m2, double coef)
{
  int Dim1 = this->Structure->GetDimension(this->FirstTargetSpace);
  int Dim2 = this->Structure->GetDimension(this->FirstTargetSpace + 1);
  if ((m1->GetMatrixType() != m2->GetMatrixType()) || (m1->GetNbrRow() != Dim1) || (m2->GetNbrRow() != Dim2)) 
    return *this;
  switch (m1->GetMatrixType())
    {
    case (Matrix::RealElements | Matrix::Symmetric):
      {
	  for (int i2 = 0; i2 < Dim2; i2++)
	    for (int j2 = 0; j2 <= i2; j2++)
	      {
		double Coef = coef * (*((RealSymmetricMatrix*) m2))(i2, j2);
		for (int i1 = 0; i1 < Dim1; i1++)
		  {
		    int Pos = i1 + Dim1 * i2;
		    int Pos2 = Dim1 * j2;
		    int Lim;
		    if (j2 == i2)
		      Lim = i1;
		    else
		      Lim = Dim1 - 1;
		    for (int j1 = 0; j1 <= Lim; j1++)
		      (*(this->ElementaryMatrix))(Pos, Pos2++) += (*((RealSymmetricMatrix*) m1))(i1, j1) * Coef;
									   }
	      }
      }
      break;
    case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
      {
	int Pos = 0;
	for (int i1 = 0; i1 < Dim1; i1++)
	  {
	    double Coef = coef * (*((RealDiagonalMatrix*) m1))[i1];
	    for (int i2 = 0; i2 < Dim2; i2++)
	      {
		(*(this->ElementaryMatrix))(Pos, Pos) += (Coef * (*((RealDiagonalMatrix*) m2))[i2]);
		Pos++;
	      }
	  }
      }
      break;
    case (Matrix::RealElements | Matrix::Antisymmetric):
      {
	if ((Dim2 > 1) && (Dim1 > 1))
	  {
	    for (int i2 = 0; i2 < Dim2; i2++)
	      for (int j2 = 0; j2 < i2; j2++)
		{
		  double Coef = (*((RealAntisymmetricMatrix*) m2))(j2, i2) * coef;
		  for (int i1 = 0; i1 < Dim1; i1++)
		    {
		      int Pos =i1 + Dim1 * i2;
		      int Pos2 = Dim1 * j2;
		      int j1 = 0;
		      for (; j1 < i1; j1++)
			(*(this->ElementaryMatrix))(Pos, Pos2++) += ((*((RealAntisymmetricMatrix*) m1))(j1, i1)
								     * Coef);
		      j1++;
		      Pos2++;
		      for (; j1 < Dim1; j1++)
			(*(this->ElementaryMatrix))(Pos, Pos2++) -= ((*((RealAntisymmetricMatrix*) m1))(i1, j1)
								     * Coef);
		    }
		}
	  }
      }
      break;
    }
  return *this;
}

// add tensor product of two matrices to the current tensor with a given space decomposition
//
// m1 = reference on the matrix associated to first target space
// m2 = reference on  the matrix associated to second  target space
// decomposition1 = first space decomposition
// decomposition2 = second space decomposition
// coef = global multiplicative factor
// return value = reference to current tensor

TwoSpaceTensor& TwoSpaceTensor::AddTensorProductMatrices(Matrix& m1, Matrix& m2, 
							 SpaceDecomposition& decomposition1, 
							 SpaceDecomposition& decomposition2, double coef)
{
  if (m2.GetMatrixType() != m1.GetMatrixType())
    return *this;
  else
    {
      switch (m1.GetMatrixType())
	{
	case (Matrix::RealElements | Matrix::Symmetric):
	  {
	    return this->AddTensorProductMatrices((RealSymmetricMatrix&) m1, 
						  (RealSymmetricMatrix&) m2, decomposition1, 
						  decomposition2, coef);
	  }
	  break;
	case (Matrix::RealElements | Matrix::Antisymmetric):
	  {
	    return this->AddTensorProductMatrices((RealAntisymmetricMatrix&) m1, 
						  (RealAntisymmetricMatrix&) m2, decomposition1, 
						  decomposition2, coef);
	  }
	  break;
	}
    }
  return *this;
}

// add tensor product of two matrices to the current tensor with a given space decomposition
//
// m1 = reference on the matrix associated to first target space
// m2 = reference on  the matrix associated to second  target space
// decomposition1 = first space decomposition
// decomposition2 = second space decomposition
// coef = global multiplicative factor
// return value = reference to current tensor

TwoSpaceTensor& TwoSpaceTensor::AddTensorProductMatrices(RealSymmetricMatrix& m1, RealSymmetricMatrix& m2, 
							 SpaceDecomposition& decomposition1, 
							 SpaceDecomposition& decomposition2, double coef)
{
  int TmpNbrRow = (decomposition1.GetSubspaceDimension(decomposition1.GetNbrSubspace() - 1) * 
		   decomposition2.GetSubspaceDimension(decomposition1.GetNbrSubspace() - 1));
  int* SubspacePosition = new int [decomposition1.GetNbrSubspace()];
  SubspacePosition[0] = 0;
  for (int i = 1; i < decomposition1.GetNbrSubspace(); i++)
    {
      SubspacePosition[i] = (decomposition1.GetSubspaceDimension(i - 1) * 
			     decomposition2.GetSubspaceDimension(i - 1));
      TmpNbrRow += SubspacePosition[i];
      SubspacePosition[i] += SubspacePosition[i - 1];
    }
  if ((TmpNbrRow != this->NbrRow) || (TmpNbrRow != this->NbrColumn))
    return *this;
  double Coef;
  int I = 0;
  int J = 0;
  int Lim = 0;
  for (int b1 = 0; b1 < decomposition1.GetNbrSubspace(); b1++)
    {
      for (int i2 = 0; i2 < decomposition2.GetSubspaceDimension(b1); i2 ++)
	{
	  I = SubspacePosition[b1] + i2 * decomposition1.GetSubspaceDimension(b1);
	  for (int j2 = i2; j2 < decomposition2.GetSubspaceDimension(b1); j2 ++)
	    {
	      Coef = m2 (i2 + decomposition2.GetSubspacePosition(b1), 
			 j2 + decomposition2.GetSubspacePosition(b1)) * coef;
	      J = SubspacePosition[b1] + j2 * decomposition1.GetSubspaceDimension(b1);
	      for (int i1 = 0; i1 < decomposition1.GetSubspaceDimension(b1); i1 ++)
		{
		  if (I == J)
		    Lim = i1;
		  else
		    Lim = 0;
		  for (int j1 = Lim; j1 < decomposition1.GetSubspaceDimension(b1); j1 ++)
		    {
		      (*((RealSymmetricMatrix*) (this->ElementaryMatrix)))(I + i1, J + j1) += 
			m1 (i1 + decomposition1.GetSubspacePosition(b1), 
			    j1 + decomposition1.GetSubspacePosition(b1)) * Coef;
		    }
		}
	    }
	}
      for (int b2 = b1 + 1; b2 < decomposition1.GetNbrSubspace(); b2++)
	for (int i2 = 0; i2 < decomposition2.GetSubspaceDimension(b1); i2 ++)
	  {
	    I = SubspacePosition[b1] + i2 * decomposition1.GetSubspaceDimension(b1);
	    for (int j2 = 0; j2 < decomposition2.GetSubspaceDimension(b2); j2 ++)
	      {
		Coef = m2 (i2 + decomposition2.GetSubspacePosition(b1), 
			   j2 + decomposition2.GetSubspacePosition(b2)) * coef;
		J = SubspacePosition[b2] + j2 * decomposition1.GetSubspaceDimension(b2);
		for (int i1 = 0; i1 < decomposition1.GetSubspaceDimension(b1); i1 ++)
		  for (int j1 = 0; j1 < decomposition1.GetSubspaceDimension(b2); j1 ++)
		    {
		      (*((RealSymmetricMatrix*) (this->ElementaryMatrix)))(I + i1, J + j1) += 
			m1 (i1 + decomposition1.GetSubspacePosition(b1), 
			    j1 + decomposition1.GetSubspacePosition(b2)) * Coef;
		    }
	      }
	  }
    }
  return *this;
}

// add tensor product of two matrices to the current tensor with a given space decomposition
//
// m1 = reference on the matrix associated to first target space
// m2 = reference on  the matrix associated to second  target space
// decomposition1 = first space decomposition
// decomposition2 = second space decomposition
// coef = global multiplicative factor
// return value = reference to current tensor

TwoSpaceTensor& TwoSpaceTensor::AddTensorProductMatrices(RealAntisymmetricMatrix& m1, 
							 RealAntisymmetricMatrix& m2, 
							 SpaceDecomposition& decomposition1, 
							 SpaceDecomposition& decomposition2, double coef)
{
  int TmpNbrRow = (decomposition1.GetSubspaceDimension(decomposition1.GetNbrSubspace() - 1) * 
		   decomposition2.GetSubspaceDimension(decomposition1.GetNbrSubspace() - 1));
  int* SubspacePosition = new int [decomposition1.GetNbrSubspace()];
  SubspacePosition[0] = 0;
  for (int i = 1; i < decomposition1.GetNbrSubspace(); i++)
    {
      SubspacePosition[i] = (decomposition1.GetSubspaceDimension(i - 1) * 
			     decomposition2.GetSubspaceDimension(i - 1));
      TmpNbrRow += SubspacePosition[i];
      SubspacePosition[i] += SubspacePosition[i - 1];
    }
  if ((TmpNbrRow != this->NbrRow) || (TmpNbrRow != this->NbrColumn))
    return *this;
  double Coef;
  int I = 0;
  int J = 0;
  for (int b1 = 0; b1 < decomposition1.GetNbrSubspace(); b1++)
    {
      for (int i2 = 0; i2 < decomposition2.GetSubspaceDimension(b1); i2 ++)
	{
	  I = SubspacePosition[b1] + i2 * decomposition1.GetSubspaceDimension(b1);
	  for (int j2 = i2 + 1; j2 < decomposition2.GetSubspaceDimension(b1); j2 ++)
	    {
	      Coef = m2 (i2 + decomposition2.GetSubspacePosition(b1), 
			 j2 + decomposition2.GetSubspacePosition(b1)) * coef;
	      J = SubspacePosition[b1] + j2 * decomposition1.GetSubspaceDimension(b1);
	      for (int i1 = 0; i1 < decomposition1.GetSubspaceDimension(b1); i1 ++)
		{
		  if (I == J)
		    for (int j1 = i1 + 1; j1 < decomposition1.GetSubspaceDimension(b1); j1 ++)
		      {
			(*((RealSymmetricMatrix*) (this->ElementaryMatrix)))(I + i1, J + j1) += 
			  m1 (i1 + decomposition1.GetSubspacePosition(b1), 
			      j1 + decomposition1.GetSubspacePosition(b1)) * Coef;
		      }
		  else
		    {
		      for (int j1 = 0; j1 < i1; j1 ++)
			(*((RealSymmetricMatrix*) (this->ElementaryMatrix)))(I + i1, J + j1) -= 
			  m1 (j1 + decomposition1.GetSubspacePosition(b1), 
			      i1 + decomposition1.GetSubspacePosition(b1)) * Coef;
		      for (int j1 = i1 + 1; j1 < decomposition1.GetSubspaceDimension(b1); j1 ++)
			(*((RealSymmetricMatrix*) (this->ElementaryMatrix)))(I + i1, J + j1) += 
			  m1 (i1 + decomposition1.GetSubspacePosition(b1), 
			      j1 + decomposition1.GetSubspacePosition(b1)) * Coef;
		    }
		}
	    }
	}
      int TmpI = 0;
      int TmpJ = 0;
      for (int b2 = b1 + 1; b2 < decomposition1.GetNbrSubspace(); b2++)
	for (int i2 = 0; i2 < decomposition2.GetSubspaceDimension(b1); i2 ++)
	  {
	    I = SubspacePosition[b1] + i2 * decomposition1.GetSubspaceDimension(b1);
	    for (int j2 = 0; j2 < decomposition2.GetSubspaceDimension(b2); j2 ++)
	      {
		TmpI = i2 + decomposition2.GetSubspacePosition(b1);
		TmpJ = j2 + decomposition2.GetSubspacePosition(b2);
		/*		if (TmpI < TmpJ)
		  Coef = m2 (TmpI, TmpJ) * coef;
		else
		if (TmpI > TmpJ)*/
		    Coef = m2 (TmpI, TmpJ) * coef;
		J = SubspacePosition[b2] + j2 * decomposition1.GetSubspaceDimension(b2);
		for (int i1 = 0; i1 < decomposition1.GetSubspaceDimension(b1); i1 ++)
		  for (int j1 = 0; j1 < decomposition1.GetSubspaceDimension(b2); j1 ++)
		    {
		      (*((RealSymmetricMatrix*) (this->ElementaryMatrix)))(I + i1, J + j1) += 
			m1 (i1 + decomposition1.GetSubspacePosition(b1), 
			    j1 + decomposition1.GetSubspacePosition(b2)) * Coef;
		    }
	      }
	  }
    }
  return *this;
}

// substract tensor product of two matrices to the cuurent tensor (assuming identical space tensor pattern)
//
// m1 = pointer to the matrix associated to first target space
// m2 = pointer to the matrix associated to first target space
// return value = reference to current tensor

TwoSpaceTensor& TwoSpaceTensor::SubTensorProductMatrices (Matrix* m1, Matrix* m2)
{
  int Dim1 = this->Structure->GetDimension(this->FirstTargetSpace);
  int Dim2 = this->Structure->GetDimension(this->FirstTargetSpace + 1);
  if ((m1->GetMatrixType() != m2->GetMatrixType()) || (m1->GetNbrRow() != Dim1) || (m2->GetNbrRow() != Dim2)) 
    return *this;
  switch (m1->GetMatrixType())
    {
    case (Matrix::RealElements | Matrix::Symmetric):
      {
	for (int i1 = 0; i1 < Dim1; i1++)
	  for (int i2 = 0; i2 < Dim2; i2++)
	    {
	      int Pos =i1 + Dim1 * i2;
	      for (int j2 = 0; j2 <= i2; j2++)
		{
		  int Lim;
		  if (j2 == i2)
		    Lim = i1;
		  else
		    Lim = Dim1 - 1;
		  for (int j1 = 0; j1 <= Lim; j1++)
		    (*(this->ElementaryMatrix))(Pos, j1 + Dim1 * j2) -= ((*((RealSymmetricMatrix*) m1))(i1, j1)
									 * (*((RealSymmetricMatrix*) m2))(i2, j2));
		    }
		}
	  }
	  break;
	case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
	  {
	    int Dim1 = this->Structure->GetDimension(this->FirstTargetSpace);
	    int Dim2 = this->Structure->GetDimension(this->FirstTargetSpace + 1);
	    for (int i1 = 0; i1 < Dim1; i1++)
	      for (int i2 = 0; i2 < Dim2; i2++)
		{
		  int Pos = i1 + Dim1 * i2;
		  (*(this->ElementaryMatrix))(Pos, Pos) -= ((*((RealDiagonalMatrix*) m1))[i1]
							    * (*((RealDiagonalMatrix*) m2))[i2]);
		}
	  }
	  break;
	case (Matrix::RealElements | Matrix::Antisymmetric):
	  {
	    if ((Dim2 > 1) && (Dim1 > 1))
	      {
		for (int i2 = 0; i2 < Dim2; i2++)
		  for (int j2 = 0; j2 < i2; j2++)
		    {
		      double Coef = (*((RealAntisymmetricMatrix*) m2))(j2, i2);
		      for (int i1 = 0; i1 < Dim1; i1++)
			{
			  int Pos =i1 + Dim1 * i2;
			  int j1 = 0;
			  for (; j1 < i1; j1++)
			    (*(this->ElementaryMatrix))(Pos, j1 + Dim1 * j2) -= ((*((RealAntisymmetricMatrix*) m1))(j1, i1)
										 * Coef);
			  j1++;
			  for (; j1 < Dim1; j1++)
			    (*(this->ElementaryMatrix))(Pos, j1 + Dim1 * j2) += ((*((RealAntisymmetricMatrix*) m1))(i1, j1)
										 * Coef);
			}
		    }
	      }
	  }
	  break;
    }
  return *this;
}

// add a one space tensor to the current two space tensor
//
// tensor = one space tensor to add 
// return value = reference to current tensor

TwoSpaceTensor& TwoSpaceTensor::operator += (const OneSpaceTensor& tensor)
{
  if (((*(tensor.Structure)) != (*(this->Structure))) || ((tensor.TargetSpace != this->FirstTargetSpace) && 
							  (tensor.TargetSpace != (this->FirstTargetSpace + 1))))
    return *this;
  if (tensor.TargetSpace == this->FirstTargetSpace)
    switch (tensor.ElementaryMatrix->GetMatrixType())
      {
      case Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal :
	{
	  this->AddToFirstSpace((RealDiagonalMatrix&) (*(tensor.ElementaryMatrix)));
	  return *this;
	}
      break;
      case Matrix::RealElements | Matrix::TriDiagonal | Matrix::Symmetric :
	{
	  this->AddToFirstSpace((RealTriDiagonalSymmetricMatrix&)  (*(tensor.ElementaryMatrix)));
	  return *this;
	}
      break;
      default:
	return *this;
      }
  else
    switch (tensor.ElementaryMatrix->GetMatrixType())
      {
      case Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal :
	{
	  this->AddToSecondSpace((RealDiagonalMatrix&) (*(tensor.ElementaryMatrix)));
	  return *this;
	}
      break;
      case Matrix::RealElements | Matrix::TriDiagonal | Matrix::Symmetric :
	{
	  this->AddToSecondSpace((RealTriDiagonalSymmetricMatrix&)  (*(tensor.ElementaryMatrix)));
	  return *this;
	}
      break;
      default:
	return *this;
      }
}

// add a real diagonal matrix to matrix of first target space
//
// matrix = matrix to add

void TwoSpaceTensor::AddToFirstSpace(RealDiagonalMatrix& matrix)
{
  int Dim1 = this->Structure->GetDimension(this->FirstTargetSpace);
  int Dim2 = this->Structure->GetDimension(this->FirstTargetSpace + 1);
  int Pos;
  for (int i2 = 0; i2 < Dim2; i2++)
    {
      Pos = i2 * Dim1;
      for (int i1 = 0; i1 < Dim1; i1++)
	{
	  (*(this->ElementaryMatrix)) (Pos, Pos) += matrix [i1];
	  Pos++;
	}
    }
}

void TwoSpaceTensor::AddToFirstSpace(RealDiagonalMatrix& matrix, SpaceDecomposition& decomposition1, 
				     SpaceDecomposition& decomposition2)
{
  int I = 0;
  int TmpNbrRow = (decomposition1.GetSubspaceDimension(decomposition1.GetNbrSubspace() - 1) * 
		   decomposition2.GetSubspaceDimension(decomposition1.GetNbrSubspace() - 1));
  int* SubspacePosition = new int [decomposition1.GetNbrSubspace()];
  SubspacePosition[0] = 0;
  for (int i = 1; i < decomposition1.GetNbrSubspace(); i++)
    {
      SubspacePosition[i] = (decomposition1.GetSubspaceDimension(i - 1) * 
			     decomposition2.GetSubspaceDimension(i - 1));
      TmpNbrRow += SubspacePosition[i];
      SubspacePosition[i] += SubspacePosition[i - 1];
    }
  for (int b1 = 0; b1 < decomposition1.GetNbrSubspace(); b1++)
    {
      for (int i2 = 0; i2 < decomposition2.GetSubspaceDimension(b1); i2 ++)
	{
	  I = SubspacePosition[b1] + i2 * decomposition1.GetSubspaceDimension(b1);
	  for (int i1 = 0; i1 < decomposition1.GetSubspaceDimension(b1); i1 ++)
	    {
	      (*(this->ElementaryMatrix))(I + i1, I + i1) += 
		matrix [i1 + decomposition1.GetSubspacePosition(b1)]; 
	    }
	}
    }
}

// add a matrix to matrix of first target space
//
// matrix = matrix to add

void TwoSpaceTensor::AddToFirstSpace(Matrix& matrix)
{
  switch (matrix.GetMatrixType())
    {
    case (Matrix::RealElements | Matrix::Symmetric):
      {
	this->AddToFirstSpace((RealSymmetricMatrix&) matrix);
      }
      break;
    case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
      {
	this->AddToFirstSpace((RealDiagonalMatrix&) matrix);
      }
      break;
    }
}
  
// add a real symmetric matrix to matrix of first target space
//
// matrix = matrix to add

void TwoSpaceTensor::AddToFirstSpace(RealSymmetricMatrix& matrix)
{
  int Dim1 = this->Structure->GetDimension(this->FirstTargetSpace);
  if (matrix.GetNbrRow() != Dim1) 
    return;
  int Dim2 = this->Structure->GetDimension(this->FirstTargetSpace + 1);
  double Coef;
  int Pos;
  for (int i1 = 0; i1 < Dim1; i1++)
    for (int j1 = 0; j1 <= i1; j1++)
      {
	Coef = matrix(i1, j1);
	Pos = 0;
	for (int i2 = 0; i2 < Dim2; i2++)
	  {
	    (*(this->ElementaryMatrix))(i1 + Pos, Pos) += Coef;
	    Pos += Dim1; 
	  }
      }
}

void TwoSpaceTensor::AddToSecondSpace(RealDiagonalMatrix& matrix, SpaceDecomposition& decomposition1, 
				      SpaceDecomposition& decomposition2)
{
  int TmpNbrRow = (decomposition1.GetSubspaceDimension(decomposition1.GetNbrSubspace() - 1) * 
		   decomposition2.GetSubspaceDimension(decomposition1.GetNbrSubspace() - 1));
  int* SubspacePosition = new int [decomposition1.GetNbrSubspace()];
  SubspacePosition[0] = 0;
  for (int i = 1; i < decomposition1.GetNbrSubspace(); i++)
    {
      SubspacePosition[i] = (decomposition1.GetSubspaceDimension(i - 1) * 
			     decomposition2.GetSubspaceDimension(i - 1));
      TmpNbrRow += SubspacePosition[i];
      SubspacePosition[i] += SubspacePosition[i - 1];
    }
  int I = 0;
  double Coef = 0.0;
  for (int b1 = 0; b1 < decomposition1.GetNbrSubspace(); b1++)
    {
      for (int i2 = 0; i2 < decomposition2.GetSubspaceDimension(b1); i2 ++)
	{
	  I = SubspacePosition[b1] + i2 * decomposition1.GetSubspaceDimension(b1);
	  Coef = matrix [i2 + decomposition2.GetSubspacePosition(b1)];
	  for (int i1 = 0; i1 < decomposition1.GetSubspaceDimension(b1); i1 ++)
	      (*(this->ElementaryMatrix))(I + i1, I + i1) += Coef;
	}
    }
}

// add a real diagonal matrix to matrix of second target space
//
// matrix = matrix to add

void TwoSpaceTensor::AddToSecondSpace(RealDiagonalMatrix& matrix)
{
  int Dim1 = this->Structure->GetDimension(this->FirstTargetSpace);
  int Dim2 = this->Structure->GetDimension(this->FirstTargetSpace + 1);
  int Pos;
  for (int i1 = 0; i1 < Dim1; i1++)
    {
      for (int i2 = 0; i2 < Dim2; i2++)
	{
	  Pos = i1 + i2 * Dim1;
	  (*(this->ElementaryMatrix)) (Pos, Pos) += matrix [i2];
	}
    }
}

// add a matrix to matrix of second target space
//
// matrix = matrix to add

void TwoSpaceTensor::AddToSecondSpace(Matrix& matrix)
{
  switch (matrix.GetMatrixType())
    {
    case (Matrix::RealElements | Matrix::Symmetric):
      {
	this->AddToSecondSpace((RealSymmetricMatrix&) matrix);
      }
      break;
    case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
      {
	this->AddToSecondSpace((RealDiagonalMatrix&) matrix);
      }
      break;
    }
}
  
// add a real symmetric matrix to matrix of second target space
//
// matrix = matrix to add

void TwoSpaceTensor::AddToSecondSpace(RealSymmetricMatrix& matrix)
{
  int Dim2 = this->Structure->GetDimension(this->FirstTargetSpace + 1);
  if (matrix.GetNbrRow() != Dim2) 
    return;
  int Dim1 = this->Structure->GetDimension(this->FirstTargetSpace);
  int Pos;
  for (int i2 = 0; i2 < Dim2; i2++)
    for (int j2 = 0; j2 <= i2; j2++)
      {
	double Coef = matrix(i2, j2);
	Pos = Dim1 * i2;
	for (int i1 = 0; i1 < Dim1; i1++)
	  {
	    (*(this->ElementaryMatrix))(Pos, Pos) += Coef;
	    Pos++;
	  }
      }
}

// add a real tridiagonal symmetric matrix to matrix of first target space
//
// matrix = matrix to add

void TwoSpaceTensor::AddToFirstSpace(RealTriDiagonalSymmetricMatrix& matrix)
{
}

// add a real tridiagonal symmetric matrix to matrix of second target space
//
// matrix = matrix to add

void TwoSpaceTensor::AddToSecondSpace(RealTriDiagonalSymmetricMatrix& matrix)
{
}

// evaluate matrix trace
//
// return value = matrix trace 

double TwoSpaceTensor::Tr () 
{
  return ((this->ElementaryMatrix->Tr () * (double) this->Structure->GetTotalDimension()) / 
	  (double) (this->Structure->GetDimension(this->FirstTargetSpace) * this->Structure->GetDimension(this->FirstTargetSpace + 1)));
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double TwoSpaceTensor::Det () 
{
  return 1.0;
}

// Output Stream overload
//
// str = reference on output stream
// tensor = tensor to print
// return value = reference on output stream

ostream& operator << (ostream& str, const TwoSpaceTensor& tensor)
{
  str << " tensor acting on spaces " << tensor.FirstTargetSpace << " and " << (tensor.FirstTargetSpace + 1) << endl;
  str << *(tensor.ElementaryMatrix) << endl;
  return str;
}

