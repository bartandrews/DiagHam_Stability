////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                     class for tensor product real vector                   //
//                                                                            //
//                        last modification : 23/03/2001                      //
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


#include "TensorProduct/TensorProductRealVector.h"
#include "TensorProduct/TensorProductStructure.h"
#include "TensorProduct/CompositeTensorProductStructure.h"
#include "Tensor/OneSpaceTensor.h"
#include "Tensor/TwoSpaceTensor.h"
#include "GeneralTools/ListIterator.h"
#include "Matrix/BlockDiagonalMatrix.h"

#include <iostream>




// default constructor
//

TensorProductRealVector::TensorProductRealVector() 
{
}

// constructor for an empty tensor product real vector 
//
// struture = reference on tensor product structure
// zeroFlag = true if all coordinates have to be set to zero

TensorProductRealVector::TensorProductRealVector(AbstractTensorProductStructure* structure, bool zeroFlag) 
{
  this->Structure = structure;
  this->Dimension = this->Structure->GetTotalDimension();
  this->TrueDimension = this->Structure->GetTotalDimension();
  this->GlobalVector = RealVector(this->Dimension, zeroFlag);
}

// constructor from an array of doubles
//
// struture = reference on tensor product structure
// array = array of doublesn

TensorProductRealVector::TensorProductRealVector(AbstractTensorProductStructure* structure, double* array) 
{
  this->Structure = structure;
  this->Dimension = this->Structure->GetTotalDimension();
  this->TrueDimension = this->Structure->GetTotalDimension();
  this->GlobalVector = RealVector(array, this->Dimension);  
}

// constructor from a real vector
//
// struture = reference on tensor product structure
// v = reference on vector containing datas

TensorProductRealVector::TensorProductRealVector(AbstractTensorProductStructure* structure, RealVector& v) 
{
  this->Structure = structure;
  this->Dimension = v.GetVectorDimension();
  this->TrueDimension = this->Dimension;
  this->GlobalVector = v;
}

// copy constructor
//
// vector = vector to copy
// duplicateFlag = true if datas have to be duplicated

TensorProductRealVector::TensorProductRealVector(const TensorProductRealVector& vector, bool duplicateFlag) 
{
  this->Structure = vector.Structure;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  this->GlobalVector = RealVector(vector.GlobalVector, duplicateFlag);   
}

// destructor
//

TensorProductRealVector::~TensorProductRealVector () 
{  
}

// assignement
//
// vector = vector to assign
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::operator = (const TensorProductRealVector& vector) 
{
  this->Structure = vector.Structure;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  this->GlobalVector = RealVector(vector.GlobalVector); 
  return *this;
}

// Resize vector
//
// dimension = new dimension

void TensorProductRealVector::Resize (int dimension) 
{
}

// Resize vector and set to zero all components that have been added
//
// dimension = new dimension

void TensorProductRealVector::ResizeAndClean (int dimension) 
{
}

// Resize vector
//
// structure = new product tensor structure

void TensorProductRealVector::Resize (AbstractTensorProductStructure* structure) 
{
}

// Resize vector and set to zero all components that have been added
//
// structure = new product tensor structure

void TensorProductRealVector::ResizeAndClean (AbstractTensorProductStructure* structure) 
{
}

// change sign of a vector
//
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::operator - () 
{
  this->GlobalVector *= -1.0;
  return *this;
}

// return a new vector with opposite sign form a given source vector
//
// V1 = source vector
// return value = new vector

TensorProductRealVector operator - (const TensorProductRealVector& V1) 
{
  return TensorProductRealVector();
}

// scalar product between two vectors
//
// V1 = first vector
// V2 = second vector
// return value = result of scalar product

double operator * (TensorProductRealVector& V1, TensorProductRealVector& V2) 
{
  return (V1.GlobalVector * V2.GlobalVector);
}

// sum two vectors
//
// V1 = vector to add
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::operator += (TensorProductRealVector& V1) 
{
  if (this->Structure != V1.Structure)
    return *this;
  this->GlobalVector += V1.GlobalVector;
  return *this;  
}

// substract two vectors
//
// V1 = first vector
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::operator -= (TensorProductRealVector& V1) 
{
  if (this->Structure != V1.Structure)
    return *this;
  this->GlobalVector -= V1.GlobalVector;
  return *this;  
}

// sum two vectors
//
// V1 = first vector
// V2 = second vector
// return value = resulting vector

TensorProductRealVector operator + (TensorProductRealVector& V1, TensorProductRealVector& V2) 
{
  return TensorProductRealVector();
}

// substract two vectors
//
// V1 = first vector
// V2 = second vector
// return value = resulting vector

TensorProductRealVector operator - (TensorProductRealVector& V1, TensorProductRealVector& V2) 
{
  return TensorProductRealVector();
}

// add a linear combination to a given vector
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::AddLinearCombination (double x, TensorProductRealVector& V) 
{
  if (this->Structure != V.Structure)
    return *this;
  this->GlobalVector.AddLinearCombination(x, V.GlobalVector);
  return *this;  
}

// add a linear combination of two vectors to a given vector
//
// x1 = multiplicative coefficient of first vector
// v1 = first vector to add
// x2 = multiplicative coefficient of first vector
// v2 = first vector to add
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::AddLinearCombination (double x1, TensorProductRealVector& v1, double x2, 
									TensorProductRealVector& v2)
{
  if ((this->Structure != v1.Structure) || (this->Structure != v2.Structure))
    return *this;
  this->GlobalVector.AddLinearCombination(x1, v1.GlobalVector, x2, v2.GlobalVector);
  return *this;  
}

// multiply a vector with a real number on the right hand side
//
// V1 = vector to multiply
// d = real to use
// return value = resulting vector

TensorProductRealVector operator * (TensorProductRealVector& V1, double d) 
{
  RealVector Tmp = V1.GlobalVector * d;
  return TensorProductRealVector(V1.Structure, Tmp);
}

// multiply a vector with a real number on the left hand side
//
// V1 = vector to multiply
// d = real to use
// return value = resulting vector

TensorProductRealVector operator * (double d, TensorProductRealVector& V1) 
{
  RealVector Tmp = V1.GlobalVector * d;
  return TensorProductRealVector(V1.Structure, Tmp);
}

// multiply a vector with a real number on the right hand side
//
// d = real to use
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::operator *= (double d) 
{
  this->GlobalVector *= d;
  return *this;  
}

// divide a vector by a real number on the right hand side
//
// d = real to use
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::operator /= (double d) 
{
  this->GlobalVector /= d;
  return *this;  
}

// left multiply a vector with a real symmetric matrix (without using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::operator *= (const RealSymmetricMatrix&  M) 
{
  this->GlobalVector *= M;
  return *this;
}

// left multiply a vector with a real tridiagonal symmetric matrix (without using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::operator *= (const RealTriDiagonalSymmetricMatrix&  M) 
{
  this->GlobalVector *= M;
  return *this;
}

// left multiply a vector with a symmetric matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply  
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::Multiply (const RealSymmetricMatrix&  M, TensorProductRealVector& V) 
{
  if (this->Structure != V.Structure)
    return *this;
  this->GlobalVector.Multiply(M, V.GlobalVector);
  return *this;
}

// left multiply a vector with a real matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::Multiply (const RealMatrix&  M, TensorProductRealVector& V) 
{
  if (this->Structure != V.Structure)
    return *this;
  this->GlobalVector.Multiply(M, V.GlobalVector);
  return *this;
}

// left multiply a vector with a one space tensor and use to store result
//  in current vector (without creating temporary vector)
//
// T = tensor to use
// V = vector to multiply  
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::Multiply (const OneSpaceTensor& T, 
							    TensorProductRealVector& V) 
{
  if (((*(this->Structure)) != (*(V.Structure))) || ((*(T.Structure)) != (*(V.Structure))))
    return *this;
  switch (this->Structure->GetTensorProductStructureType())
    {
    case AbstractTensorProductStructure::Simple:
      {
	for (int i = 0; i < this->Structure->GetTotalDimension(); i += 
	       ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace + 1))
	  for (int j = 0; j < ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace); j++)
	    this->GlobalVector.Multiply(*(T.ElementaryMatrix), V.GlobalVector, i + j, 
					((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace), 
					i + j, ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace));
      }
      break;
    case AbstractTensorProductStructure::Composite:
      {
	if (T.MatrixType & Matrix::BlockDiagonal)
	  {
	    Matrix** TmpMatrix;
	    TensorProductStructure* TmpStructure;
	    int Pos = 0;
	    int Inc = ((CompositeTensorProductStructure*) this->Structure)->GetSubspaceIncrement(Pos++);
	    ListIterator<Matrix*> IterMatrix(((BlockDiagonalMatrix*) T.ElementaryMatrix)->Blocks);
	    ListIterator<TensorProductStructure> IterStructure(((CompositeTensorProductStructure*) 
								this->Structure)->
							       GetTensorProductStructures());
	    if ((TmpMatrix = IterMatrix()))
	      {
		TmpStructure = IterStructure();
		for (int i = 0; i < TmpStructure->GetTotalDimension(); 
		     i += TmpStructure->GetIncrement(T.TargetSpace + 1))
		  for (int j = 0; j < TmpStructure->GetIncrement(T.TargetSpace); j++)
		    this->GlobalVector.Multiply(**TmpMatrix, V.GlobalVector, i + j + Inc, 
						TmpStructure->GetIncrement(T.TargetSpace), 
						i + j + Inc, TmpStructure->GetIncrement(T.TargetSpace));
	      }
	    while ((TmpMatrix = IterMatrix()))
	      {
		Inc = ((CompositeTensorProductStructure*) this->Structure)->GetSubspaceIncrement(Pos++);
		TmpStructure = IterStructure();
		for (int i = 0; i < TmpStructure->GetTotalDimension(); 
		     i += TmpStructure->GetIncrement(T.TargetSpace + 1))
		  for (int j = 0; j < TmpStructure->GetIncrement(T.TargetSpace); j++)
		    this->GlobalVector.Multiply(**TmpMatrix, V.GlobalVector, i + j + Inc, 
						TmpStructure->GetIncrement(T.TargetSpace), 
						i + j + Inc, 
						TmpStructure->GetIncrement(T.TargetSpace));
	      }
	  }
      }
      break;
    }
  return *this;
}

// left multiply a vector with a one space tensor for a given range of indices 
// and use to store result in current vector (without creating temporary vector)
//
// T = tensor to use
// V = vector to multiply  
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::Multiply (const OneSpaceTensor& T, 
							    TensorProductRealVector& V,
							    int firstComponent, int nbrComponent)
{
  if (((*(this->Structure)) != (*(V.Structure))) || ((*(T.Structure)) != (*(V.Structure))))
    return *this;
  switch (this->Structure->GetTensorProductStructureType())
    {
    case AbstractTensorProductStructure::Simple:
      {
	for (int i = 0; i < this->Structure->GetTotalDimension(); i += 
	       ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace + 1))
	  for (int j = 0; j < ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace); j++)
	    this->GlobalVector.Multiply(*(T.ElementaryMatrix), V.GlobalVector, i + j, 
					((TensorProductStructure*) this->Structure)->
					GetIncrement(T.TargetSpace), 
					i + j, ((TensorProductStructure*) this->Structure)->
					GetIncrement(T.TargetSpace));
      }
      break;
    case AbstractTensorProductStructure::Composite:
      {
	if (T.MatrixType & Matrix::BlockDiagonal)
	  {
	    Matrix** TmpMatrix;
	    TensorProductStructure* TmpStructure;
	    int Pos = 0;
	    int Inc = ((CompositeTensorProductStructure*) this->Structure)->GetSubspaceIncrement(Pos++);
	    ListIterator<Matrix*> IterMatrix(((BlockDiagonalMatrix*) T.ElementaryMatrix)->Blocks);
	    ListIterator<TensorProductStructure> IterStructure(((CompositeTensorProductStructure*) 
								this->Structure)->
							       GetTensorProductStructures());
	    if ((TmpMatrix = IterMatrix()))
	      {
		TmpStructure = IterStructure();
		for (int i = 0; i < TmpStructure->GetTotalDimension(); 
		     i += TmpStructure->GetIncrement(T.TargetSpace + 1))
		  for (int j = 0; j < TmpStructure->GetIncrement(T.TargetSpace); j++)
		    this->GlobalVector.Multiply(**TmpMatrix, V.GlobalVector, i + j + Inc, 
						TmpStructure->GetIncrement(T.TargetSpace), 
						i + j + Inc, TmpStructure->GetIncrement(T.TargetSpace));
	      }
	    while ((TmpMatrix = IterMatrix()))
	      {
		Inc = ((CompositeTensorProductStructure*) this->Structure)->GetSubspaceIncrement(Pos++);
		TmpStructure = IterStructure();
		for (int i = 0; i < TmpStructure->GetTotalDimension(); 
		     i += TmpStructure->GetIncrement(T.TargetSpace + 1))
		  for (int j = 0; j < TmpStructure->GetIncrement(T.TargetSpace); j++)
		    this->GlobalVector.Multiply(**TmpMatrix, V.GlobalVector, i + j + Inc, 
						TmpStructure->GetIncrement(T.TargetSpace), 
						i + j + Inc, 
						TmpStructure->GetIncrement(T.TargetSpace));
	      }
	  }
      }
      break;
    }
  return *this;
}

// left multiply a vector with a one space tensor and add result to current vector (without creating temporary vector)
//
// T = tensor to use
// V = vector to multiply  
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::AddMultiply (const OneSpaceTensor& T, TensorProductRealVector& V) 
{
  if (((*(this->Structure)) != (*(V.Structure))) || ((*(T.Structure)) != (*(V.Structure))))
    return *this;
  switch (this->Structure->GetTensorProductStructureType())
    {
    case AbstractTensorProductStructure::Simple:
      {
	for (int i = 0; i < this->Structure->GetTotalDimension(); i += 
	       ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace + 1))
	  for (int j = 0; j < ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace); j++)
	    this->GlobalVector.AddMultiply(*(T.ElementaryMatrix), V.GlobalVector, i + j, 
					   ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace), 
					   i + j, ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace));
      }
      break;
    case AbstractTensorProductStructure::Composite:
      {
	if (T.MatrixType & Matrix::BlockDiagonal)
	  {
	    Matrix** TmpMatrix;
	    TensorProductStructure* TmpStructure;
	    int Pos = 0;
	    int Inc = ((CompositeTensorProductStructure*) this->Structure)->GetSubspaceIncrement(Pos++);
	    ListIterator<Matrix*> IterMatrix(((BlockDiagonalMatrix*) T.ElementaryMatrix)->Blocks);
	    ListIterator<TensorProductStructure> IterStructure(((CompositeTensorProductStructure*) 
								this->Structure)->
							       GetTensorProductStructures());
	    if ((TmpMatrix = IterMatrix()))
	      {
		TmpStructure = IterStructure();
		for (int i = 0; i < TmpStructure->GetTotalDimension(); 
		     i += TmpStructure->GetIncrement(T.TargetSpace + 1))
		  for (int j = 0; j < TmpStructure->GetIncrement(T.TargetSpace); j++)
		    this->GlobalVector.AddMultiply(**TmpMatrix, V.GlobalVector, i + j + Inc, 
						   TmpStructure->GetIncrement(T.TargetSpace), 
						   i + j + Inc, TmpStructure->GetIncrement(T.TargetSpace));
	      }
	    while ((TmpMatrix = IterMatrix()))
	      {
		Inc = ((CompositeTensorProductStructure*) this->Structure)->GetSubspaceIncrement(Pos++);
		TmpStructure = IterStructure();
		for (int i = 0; i < TmpStructure->GetTotalDimension(); 
		     i += TmpStructure->GetIncrement(T.TargetSpace + 1))
		  for (int j = 0; j < TmpStructure->GetIncrement(T.TargetSpace); j++)
		    this->GlobalVector.AddMultiply(**TmpMatrix, V.GlobalVector, i + j + Inc, 
						   TmpStructure->GetIncrement(T.TargetSpace), 
						   i + j + Inc, 
						   TmpStructure->GetIncrement(T.TargetSpace));
	      }
	  }
      }
      break;
    }
  return *this;
}

// left multiply a vector with a one space tensor for a given range of indices 
// and add result to current vector (without creating temporary vector)
//
// T = tensor to use
// V = vector to multiply  
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::AddMultiply (const OneSpaceTensor& T, 
							       TensorProductRealVector& V,
							       int firstComponent, int nbrComponent)

{
  if (((*(this->Structure)) != (*(V.Structure))) || ((*(T.Structure)) != (*(V.Structure))))
    return *this;
  switch (this->Structure->GetTensorProductStructureType())
    {
    case AbstractTensorProductStructure::Simple:
      {
	int StartI = firstComponent / ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace + 1);
	int i = StartI * ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace + 1);
	int StartJ = firstComponent - i;
	int EndI = (nbrComponent - firstComponent) / ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace + 1);
	int EndJ = nbrComponent - firstComponent - EndI * ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace + 1);
	if (StartI == EndI)
	  {
	    for (int j = StartJ; j < EndJ; ++j)
	      this->GlobalVector.AddMultiply(*(T.ElementaryMatrix), V.GlobalVector, i + j, 
					     ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace), 
					     i + j, ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace));
	  }
	else
	  {
	    for (int j = StartJ; j < ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace); ++j)
	      this->GlobalVector.AddMultiply(*(T.ElementaryMatrix), V.GlobalVector, i + j, 
					     ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace), 
					     i + j, ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace));
	    i += ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace + 1);
	    --EndI;
	    EndI *= ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace + 1);
	    for (; i < EndI; i += ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace + 1))
	      for (int j = 0; j < ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace); j++)
		this->GlobalVector.AddMultiply(*(T.ElementaryMatrix), V.GlobalVector, i + j, 
					       ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace), 
					       i + j, ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace));
	    for (int j = 0; j < EndJ; ++j)
	      this->GlobalVector.AddMultiply(*(T.ElementaryMatrix), V.GlobalVector, i + j, 
					     ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace), 
					     i + j, ((TensorProductStructure*) this->Structure)->GetIncrement(T.TargetSpace));
	  }
      }
      break;
    case AbstractTensorProductStructure::Composite:
      {
	if (T.MatrixType & Matrix::BlockDiagonal)
	  {
	    Matrix** TmpMatrix;
	    TensorProductStructure* TmpStructure;
	    int Pos = 0;
	    int Inc = ((CompositeTensorProductStructure*) this->Structure)->GetSubspaceIncrement(Pos++);
	    ListIterator<Matrix*> IterMatrix(((BlockDiagonalMatrix*) T.ElementaryMatrix)->Blocks);
	    ListIterator<TensorProductStructure> IterStructure(((CompositeTensorProductStructure*) 
								this->Structure)->
							       GetTensorProductStructures());
	    if ((TmpMatrix = IterMatrix()))
	      {
		TmpStructure = IterStructure();
		for (int i = 0; i < TmpStructure->GetTotalDimension(); 
		     i += TmpStructure->GetIncrement(T.TargetSpace + 1))
		  for (int j = 0; j < TmpStructure->GetIncrement(T.TargetSpace); j++)
		    this->GlobalVector.AddMultiply(**TmpMatrix, V.GlobalVector, i + j + Inc, 
						   TmpStructure->GetIncrement(T.TargetSpace), 
						   i + j + Inc, TmpStructure->GetIncrement(T.TargetSpace));
	      }
	    while ((TmpMatrix = IterMatrix()))
	      {
		Inc = ((CompositeTensorProductStructure*) this->Structure)->GetSubspaceIncrement(Pos++);
		TmpStructure = IterStructure();
		for (int i = 0; i < TmpStructure->GetTotalDimension(); 
		     i += TmpStructure->GetIncrement(T.TargetSpace + 1))
		  for (int j = 0; j < TmpStructure->GetIncrement(T.TargetSpace); j++)
		    this->GlobalVector.AddMultiply(**TmpMatrix, V.GlobalVector, i + j + Inc, 
						   TmpStructure->GetIncrement(T.TargetSpace), 
						   i + j + Inc, 
						   TmpStructure->GetIncrement(T.TargetSpace));
	      }
	  }
      }
      break;
    }
  return *this;
}

// left multiply a vector with a two space tensor and use to store result in current vector (without creating temporary vector)
//
// T = tensor to use
// V = vector to multiply  
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::Multiply (const TwoSpaceTensor& T, TensorProductRealVector& V)
{
  if (((*(this->Structure)) != (*(V.Structure))) || ((*(T.Structure)) != (*(V.Structure))))
    return *this;
  switch (this->Structure->GetTensorProductStructureType())
    {
    case AbstractTensorProductStructure::Simple:
      {
	for (int i = 0; i < this->Structure->GetTotalDimension(); i += 
	       ((TensorProductStructure*) this->Structure)->GetIncrement(T.FirstTargetSpace + 2))
	  for (int j = 0; j < ((TensorProductStructure*) this->Structure)->GetIncrement(T.FirstTargetSpace); j++)
	    this->GlobalVector.Multiply(*(T.ElementaryMatrix), V.GlobalVector, i + j, 
					((TensorProductStructure*) this->Structure)->GetIncrement(T.FirstTargetSpace), 
					i + j, ((TensorProductStructure*) this->Structure)->GetIncrement(T.FirstTargetSpace));
      }
      break;
    case AbstractTensorProductStructure::Composite:
      {
	if (T.MatrixType & Matrix::BlockDiagonal)
	  {
	    Matrix** TmpMatrix;
	    TensorProductStructure* TmpStructure;
	    int Pos = 0;
	    int Inc = ((CompositeTensorProductStructure*) this->Structure)->GetSubspaceIncrement(Pos++);
	    ListIterator<Matrix*> IterMatrix(((BlockDiagonalMatrix*) T.ElementaryMatrix)->Blocks);
	    ListIterator<TensorProductStructure> IterStructure(((CompositeTensorProductStructure*) this->Structure)->
							       GetTensorProductStructures());
	    if ((TmpMatrix = IterMatrix()))
	      {
		TmpStructure = IterStructure();
		for (int i = 0; i < TmpStructure->GetTotalDimension(); i += TmpStructure->GetIncrement(T.FirstTargetSpace + 2))
		  for (int j = 0; j < TmpStructure->GetIncrement(T.FirstTargetSpace); j++)
		    this->GlobalVector.Multiply(**TmpMatrix, V.GlobalVector, i + j + Inc, TmpStructure->GetIncrement(T.FirstTargetSpace), 
						i + j + Inc, TmpStructure->GetIncrement(T.FirstTargetSpace));
		Inc = ((CompositeTensorProductStructure*) this->Structure)->GetSubspaceIncrement(Pos++);
	      }
	    while ((TmpMatrix = IterMatrix()))
	      {
		TmpStructure = IterStructure();
		for (int i = 0; i < TmpStructure->GetTotalDimension(); i += TmpStructure->GetIncrement(T.FirstTargetSpace + 2))
		  for (int j = 0; j < TmpStructure->GetIncrement(T.FirstTargetSpace); j++)
		    this->GlobalVector.AddMultiply(**TmpMatrix, V.GlobalVector, i + j + Inc, 
						   TmpStructure->GetIncrement(T.FirstTargetSpace), 
						   i + j + Inc, TmpStructure->GetIncrement(T.FirstTargetSpace));
		Inc = ((CompositeTensorProductStructure*) this->Structure)->GetSubspaceIncrement(Pos++);
	      }
	  }
      }
      break;
    }
  return *this;
}

// return vector i-th coordinate (without testing if position is valid)
//
// i = coordinate position

double& TensorProductRealVector::operator [] (int i) 
{
  return this->GlobalVector[i];
}

// return vector coordinate corresponding to a tensor product index
//
// index = tensor product index
// return value = reference on indexed coordinate

double& TensorProductRealVector::operator [] (TensorProductIndex& index)
{
  return this->GlobalVector[index.GetGlobalIndex()];
}

// get vector norm
//
// return value = vector norm

double TensorProductRealVector::Norm() 
{
  return this->GlobalVector.Norm();
}
  
// get square of vector norm
//
// return value = square of vector norm

double TensorProductRealVector::SqrNorm () 
{
  return this->GlobalVector.SqrNorm();
}
  
// normalize vector
//
// return value = reference on current vector

TensorProductRealVector& TensorProductRealVector::Normalize() 
{
  this->GlobalVector.Normalize();
  return *this;
}

