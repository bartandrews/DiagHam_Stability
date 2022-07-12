////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of hamiltonian defined as a linear combination            //
//                           of triple tensor products                        //
//                                                                            //
//                        last modification : 31/07/2016                      //
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


#include "Hamiltonian/TripleTensorProductSparseMatrixHamiltonian.h"
#include "MathTools/Complex.h" 
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"
#include "Architecture/ArchitectureOperation/VectorTensorMultiplicationCoreOperation.h"
#include "Architecture/ArchitectureOperation/VectorSparseTensorMultiplyOperation.h"


#include <iostream>


using std::cout;
using std::endl;


// default contructor 
//

TripleTensorProductSparseMatrixHamiltonian::TripleTensorProductSparseMatrixHamiltonian()
{
}

// contructor 
//
// nbrTensorProducts = number of tensor products whose linear combination defined the Hamiltonian 
// leftMatrices = left matrices of each tensor product
  // middleMatrices = middle matrices of each tensor product
// rightMatrices = right matrices of each tensor product
// coefficients = coefficients of the ensor product linear combination
// architecture = architecture to use for precalculation

TripleTensorProductSparseMatrixHamiltonian::TripleTensorProductSparseMatrixHamiltonian(int nbrTensorProducts, SparseRealMatrix* leftMatrices, SparseRealMatrix* middleMatrices, 
										       SparseRealMatrix* rightMatrices, double* coefficients, AbstractArchitecture* architecture)
{
  this->NbrTensorProducts = nbrTensorProducts;
  this->LeftMatrices = new SparseRealMatrix[this->NbrTensorProducts];
  this->MiddleMatrices = new SparseRealMatrix[this->NbrTensorProducts];
  this->RightMatrices = new SparseRealMatrix[this->NbrTensorProducts];
  this->Coefficients = new double[this->NbrTensorProducts];
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      this->LeftMatrices[i] = leftMatrices[i];
      this->MiddleMatrices[i] = middleMatrices[i];
      this->RightMatrices[i] = rightMatrices[i];
      this->Coefficients[i] = coefficients[i];
    }
  this->RightMatrixNbrRow = this->RightMatrices[0].GetNbrRow();
  this->MiddleMatrixNbrRow = this->MiddleMatrices[0].GetNbrRow();
  this->LeftMatrixNbrRow  = this->LeftMatrices[0].GetNbrRow();
  this->HamiltonianShift = 0.0;
  this->Architecture = architecture;
  long HamiltonianDimension = this->RightMatrixNbrRow * this->MiddleMatrixNbrRow * this->LeftMatrixNbrRow;
  this->HilbertSpace = new UndescribedHilbertSpace(HamiltonianDimension);
  this->LeftHamiltonianVectorMultiplicationFlag = true;
  this->InitializeTemporaryArrays();
}

// destructor
//

TripleTensorProductSparseMatrixHamiltonian::~TripleTensorProductSparseMatrixHamiltonian() 
{
  if (this->LeftMatrices != 0)
    {
      delete[] this->MiddleMatrices;
    }
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& TripleTensorProductSparseMatrixHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									    int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  RealVector TmpVector (vSource.GetVectorDimension());
  RealVector TmpVector2 (vSource.GetVectorDimension());
  int TmpLeftRowIndex;
  int TmpMiddleRowIndex;
  int TmpRightRowIndex;
  long TmpRowPointer;
  long TmpRowLastPointer;
  for (int k = 0; k <  this->NbrTensorProducts; ++k)
    {
      SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[k];
      SparseRealMatrix& TmpMiddleMatrix = this->MiddleMatrices[k];
      SparseRealMatrix& TmpRightMatrix = this->RightMatrices[k];
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  TmpVector[i] = 0.0;
	  this->GetIndicesFromLinearizedIndex(i, TmpLeftRowIndex, TmpMiddleRowIndex, TmpRightRowIndex);
	  TmpRowPointer = TmpRightMatrix.RowPointers[TmpRightRowIndex];
	  if (TmpRowPointer >= 0l)
	    {
	      TmpRowLastPointer = TmpRightMatrix.RowLastPointers[TmpRightRowIndex];
	      for (; TmpRowPointer <= TmpRowLastPointer; ++TmpRowPointer)
		{
		  TmpVector[i] += TmpRightMatrix.MatrixElements[TmpRowPointer] * vSource[this->GetLinearizedIndex(TmpLeftRowIndex, TmpMiddleRowIndex, 
														  TmpRightMatrix.ColumnIndices[TmpRowPointer])];
// 		  cout << "right " << TmpRightMatrix.MatrixElements[TmpRowPointer]  << " " << vSource[this->GetLinearizedIndex(TmpLeftRowIndex, TmpMiddleRowIndex, TmpRightMatrix.ColumnIndices[TmpRowPointer])] 
// 		       << " " << this->GetLinearizedIndex(TmpLeftRowIndex, TmpMiddleRowIndex, TmpRightMatrix.ColumnIndices[TmpRowPointer]) << endl;
		}
	      TmpVector[i] *= this->Coefficients[k];
// 	      cout << "check right " << k << " " << TmpVector[i] << " " << this->Coefficients[k] << endl;
	    }
	}
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  TmpVector2[i] = 0.0;
	  this->GetIndicesFromLinearizedIndex(i, TmpLeftRowIndex, TmpMiddleRowIndex, TmpRightRowIndex);
	  TmpRowPointer = TmpMiddleMatrix.RowPointers[TmpMiddleRowIndex];
// 	  cout << "TmpMiddleRowIndex=" << TmpMiddleRowIndex << " " << TmpMiddleMatrix.RowPointers[TmpMiddleRowIndex] << endl;
	  if (TmpRowPointer >= 0l)
	    {
	      TmpRowLastPointer = TmpMiddleMatrix.RowLastPointers[TmpMiddleRowIndex];
	      for (; TmpRowPointer <= TmpRowLastPointer; ++TmpRowPointer)
		{
		  TmpVector2[i] += TmpMiddleMatrix.MatrixElements[TmpRowPointer] * TmpVector[this->GetLinearizedIndex(TmpLeftRowIndex, 
														      TmpMiddleMatrix.ColumnIndices[TmpRowPointer], TmpRightRowIndex)];
// 		  cout << "middle " << TmpMiddleMatrix.MatrixElements[TmpRowPointer] << " " <<  TmpVector[this->GetLinearizedIndex(TmpLeftRowIndex, 
// 																   TmpMiddleMatrix.ColumnIndices[TmpRowPointer], TmpRightRowIndex)]
// 		       << " " << this->GetLinearizedIndex(TmpLeftRowIndex, TmpMiddleMatrix.ColumnIndices[TmpRowPointer], TmpRightRowIndex) << endl;
		}
	    }
	}
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  this->GetIndicesFromLinearizedIndex(i, TmpLeftRowIndex, TmpMiddleRowIndex, TmpRightRowIndex);
	  TmpRowPointer = TmpLeftMatrix.RowPointers[TmpLeftRowIndex];
	  if (TmpRowPointer >= 0l)
	    {
	      TmpRowLastPointer = TmpLeftMatrix.RowLastPointers[TmpLeftRowIndex];
	      for (; TmpRowPointer <= TmpRowLastPointer; ++TmpRowPointer)
		{
		  vDestination[i] += TmpLeftMatrix.MatrixElements[TmpRowPointer] * TmpVector2[this->GetLinearizedIndex(TmpLeftMatrix.ColumnIndices[TmpRowPointer], TmpMiddleRowIndex,
														       TmpRightRowIndex)];
// 		  cout << "left " << TmpLeftMatrix.MatrixElements[TmpRowPointer] << " " << TmpVector2[this->GetLinearizedIndex(TmpLeftMatrix.ColumnIndices[TmpRowPointer], TmpMiddleRowIndex,
// 																TmpRightRowIndex)] 
// 		       << " " << this->GetLinearizedIndex(TmpLeftMatrix.ColumnIndices[TmpRowPointer], TmpMiddleRowIndex,
// 							  TmpRightRowIndex) << endl;	
		}
	    }
	}
    }
  
  if (this->HamiltonianShift != 0.0)
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	vDestination[i] += this->HamiltonianShift * vSource[i];
    }
  return vDestination;
//   int RightMatrixDimension = this->RightMatrices[0].GetNbrRow();
//   int LeftMatrixDimension = this->LeftMatrices[0].GetNbrRow();
//   int IndexStep = this->RightMatrices[0].GetNbrColumn();
//   int LastComponent = firstComponent + nbrComponent - 1;
//   int LeftMatrixLastIndex = LastComponent / this->RightMatrices[0].GetNbrRow();
//   int RightMatrixLastIndex = LastComponent % this->RightMatrices[0].GetNbrRow();
//   long TmpARowPointer;
//   long TmpARowLastPointer;
//   long TmpBRowPointer;
//   long TmpBRowLastPointer;

//   double** LocalTemporaryMatrix = this->TemporaryArray;

//   for (int i = 0; i < this->NbrTensorProducts; ++i)
//     {
//       SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
//       SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];

//       VectorTensorMultiplicationCoreOperation Operation(this, i, vSource);
//       Operation.ApplyOperation(this->Architecture);

//       VectorSparseTensorMultiplyOperation Operation2(this, i, &vDestination);
//       Operation2.ApplyOperation(this->Architecture);
//     }
  
//   if (this->HamiltonianShift != 0.0)
//     {
//       for (int i = firstComponent; i < LastComponent; ++i)
// 	vDestination[i] += this->HamiltonianShift * vSource[i];
//     }
//   return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* TripleTensorProductSparseMatrixHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  int IndexStep = this->RightMatrices[0].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent - 1;
  int AMatrixLastIndex = LastComponent / this->RightMatrices[0].GetNbrRow();
  int BMatrixLastIndex = LastComponent % this->RightMatrices[0].GetNbrRow();
  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  double* Tmp = new double[nbrVectors];
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
      SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
      int AMatrixStartingIndex = firstComponent / this->RightMatrices[0].GetNbrRow();
      int BMatrixStartingIndex = firstComponent % this->RightMatrices[0].GetNbrRow();
      int TotalIndex = firstComponent;
      for (; AMatrixStartingIndex <=  AMatrixLastIndex; ++AMatrixStartingIndex)
	{
	  TmpARowPointer = TmpLeftMatrix.RowPointers[AMatrixStartingIndex];
	  if (TmpARowPointer >= 0l)
	    {
	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[AMatrixStartingIndex];
	      int TmpBMatrixLastIndex = TmpRightMatrix.GetNbrRow() - 1;
	      if (AMatrixStartingIndex == (AMatrixLastIndex - 1))
		TmpBMatrixLastIndex = BMatrixLastIndex;
	      for (; BMatrixStartingIndex <=  TmpBMatrixLastIndex; ++BMatrixStartingIndex)
		{
		  TmpBRowPointer = TmpRightMatrix.RowPointers[BMatrixStartingIndex];
		  if (TmpBRowPointer >= 0l)
		    {
		      TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[BMatrixStartingIndex];
		      for (int l = 0; l < nbrVectors; ++l)
			Tmp[l] = 0.0;
		      for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
			{
			  double Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
			  int TmpIndex = TmpLeftMatrix.ColumnIndices[k] * IndexStep;
			  for (long j = TmpBRowPointer; j <= TmpBRowLastPointer; ++j)
			    {
			      int InputIndex = TmpIndex + TmpRightMatrix.ColumnIndices[j];
			      for (int l = 0; l < nbrVectors; ++l)			      
				Tmp[l] += Tmp2 * TmpRightMatrix.MatrixElements[j] * vSources[l][InputIndex];
			    }
			}
		      int OutputIndex = AMatrixStartingIndex * IndexStep + BMatrixStartingIndex;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][OutputIndex] += Tmp[l];
		    }
		}
	    }
	  BMatrixStartingIndex = 0;
	}
    }
  delete[] Tmp;
  if (this->HamiltonianShift != 0.0)
    {
      for (int k= 0; k < nbrVectors; ++k)
	{
	  RealVector& TmpDestination = vDestinations[k];
	  RealVector& TmpSource = vSources[k];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    TmpDestination[i] += this->HamiltonianShift * TmpSource[i];
	}
    }
  return vDestinations;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& TripleTensorProductSparseMatrixHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  ComplexVector TmpVector (vSource.GetVectorDimension());
  ComplexVector TmpVector2 (vSource.GetVectorDimension());
  int TmpLeftRowIndex;
  int TmpMiddleRowIndex;
  int TmpRightRowIndex;
  long TmpRowPointer;
  long TmpRowLastPointer;
  for (int k = 0; k <  this->NbrTensorProducts; ++k)
    {
      SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[k];
      SparseRealMatrix& TmpMiddleMatrix = this->MiddleMatrices[k];
      SparseRealMatrix& TmpRightMatrix = this->RightMatrices[k];
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  TmpVector[i] = 0.0;
	  this->GetIndicesFromLinearizedIndex(i, TmpLeftRowIndex, TmpMiddleRowIndex, TmpRightRowIndex);
	  TmpRowPointer = TmpRightMatrix.RowPointers[TmpRightRowIndex];
	  if (TmpRowPointer >= 0l)
	    {
	      TmpRowLastPointer = TmpRightMatrix.RowLastPointers[TmpRightRowIndex];
	      for (; TmpRowPointer <= TmpRowLastPointer; ++TmpRowPointer)
		{
		  TmpVector[i] += TmpRightMatrix.MatrixElements[TmpRowPointer] * vSource[this->GetLinearizedIndex(TmpLeftRowIndex, TmpMiddleRowIndex, 
														  TmpRightMatrix.ColumnIndices[TmpRowPointer])];
		}
	      TmpVector[i] *= this->Coefficients[k];
	    }
	}
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  TmpVector2[i] = 0.0;
	  this->GetIndicesFromLinearizedIndex(i, TmpLeftRowIndex, TmpMiddleRowIndex, TmpRightRowIndex);
	  TmpRowPointer = TmpMiddleMatrix.RowPointers[TmpMiddleRowIndex];
	  if (TmpRowPointer >= 0l)
	    {
	      TmpRowLastPointer = TmpMiddleMatrix.RowLastPointers[TmpMiddleRowIndex];
	      for (; TmpRowPointer <= TmpRowLastPointer; ++TmpRowPointer)
		{
		  TmpVector2[i] += TmpMiddleMatrix.MatrixElements[TmpRowPointer] * TmpVector[this->GetLinearizedIndex(TmpLeftRowIndex, 
														      TmpMiddleMatrix.ColumnIndices[TmpRowPointer], TmpRightRowIndex)];
		}
	    }
	}
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  this->GetIndicesFromLinearizedIndex(i, TmpLeftRowIndex, TmpMiddleRowIndex, TmpRightRowIndex);
	  TmpRowPointer = TmpLeftMatrix.RowPointers[TmpLeftRowIndex];
	  if (TmpRowPointer >= 0l)
	    {
	      TmpRowLastPointer = TmpLeftMatrix.RowLastPointers[TmpLeftRowIndex];
	      for (; TmpRowPointer <= TmpRowLastPointer; ++TmpRowPointer)
		{
		  vDestination[i] += TmpLeftMatrix.MatrixElements[TmpRowPointer] * TmpVector2[this->GetLinearizedIndex(TmpLeftMatrix.ColumnIndices[TmpRowPointer], TmpMiddleRowIndex,
														       TmpRightRowIndex)];
		}
	    }
	}
    }
    
//   int RightMatrixDimension = this->RightMatrices[0].GetNbrRow();
//   int LeftMatrixDimension = this->LeftMatrices[0].GetNbrRow();
//   int IndexStep = this->RightMatrices[0].GetNbrColumn();
//   int LastComponent = firstComponent + nbrComponent - 1;
//   int LeftMatrixLastIndex = LastComponent / this->RightMatrices[0].GetNbrRow();
//   int RightMatrixLastIndex = LastComponent % this->RightMatrices[0].GetNbrRow();
//   long TmpARowPointer;
//   long TmpARowLastPointer;
//   long TmpBRowPointer;
//   long TmpBRowLastPointer;

//   double** LocalTemporaryMatrix = this->TemporaryArray;

//   for (int i = 0; i < this->NbrTensorProducts; ++i)
//     {
//       SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
//       SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];

//       VectorTensorMultiplicationCoreOperation Operation(this, i, vSource);
//       Operation.ApplyOperation(this->Architecture);

//       VectorSparseTensorMultiplyOperation Operation2(this, i, &vDestination);
//       Operation2.ApplyOperation(this->Architecture);
//     }
  
  if (this->HamiltonianShift != 0.0)
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	vDestination[i] += this->HamiltonianShift * vSource[i];
    }
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* TripleTensorProductSparseMatrixHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, 
										       int nbrVectors, int firstComponent, int nbrComponent)
{
  int IndexStep = this->RightMatrices[0].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent - 1;
  int AMatrixLastIndex = LastComponent / this->RightMatrices[0].GetNbrRow();
  int BMatrixLastIndex = LastComponent % this->RightMatrices[0].GetNbrRow();
  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  Complex* Tmp = new Complex[nbrVectors];
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
      SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
      int AMatrixStartingIndex = firstComponent / this->RightMatrices[0].GetNbrRow();
      int BMatrixStartingIndex = firstComponent % this->RightMatrices[0].GetNbrRow();
      int TotalIndex = firstComponent;
      for (; AMatrixStartingIndex <=  AMatrixLastIndex; ++AMatrixStartingIndex)
	{
	  TmpARowPointer = TmpLeftMatrix.RowPointers[AMatrixStartingIndex];
	  if (TmpARowPointer >= 0l)
	    {
	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[AMatrixStartingIndex];
	      int TmpBMatrixLastIndex = TmpRightMatrix.GetNbrRow() - 1;
	      if (AMatrixStartingIndex == (AMatrixLastIndex - 1))
		TmpBMatrixLastIndex = BMatrixLastIndex;
	      for (; BMatrixStartingIndex <=  TmpBMatrixLastIndex; ++BMatrixStartingIndex)
		{
		  TmpBRowPointer = TmpRightMatrix.RowPointers[BMatrixStartingIndex];
		  if (TmpBRowPointer >= 0l)
		    {
		      TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[BMatrixStartingIndex];
		      for (int l = 0; l < nbrVectors; ++l)
			Tmp[l] = 0.0;
		      for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
			{
			  Complex Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
			  int TmpIndex = TmpLeftMatrix.ColumnIndices[k] * IndexStep;
			  for (long j = TmpBRowPointer; j <= TmpBRowLastPointer; ++j)
			    {
			      int InputIndex = TmpIndex + TmpRightMatrix.ColumnIndices[j];
			      for (int l = 0; l < nbrVectors; ++l)			      
				Tmp[l] += Tmp2 * TmpRightMatrix.MatrixElements[j] * vSources[l][InputIndex];
			    }
			}
		      int OutputIndex = AMatrixStartingIndex * IndexStep + BMatrixStartingIndex;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][OutputIndex] += Tmp[l];
		    }
		}
	    }
	  BMatrixStartingIndex = 0;
	}
    }
  delete[] Tmp;
  if (this->HamiltonianShift != 0.0)
    {
      for (int k= 0; k < nbrVectors; ++k)
	{
	  ComplexVector& TmpDestination = vDestinations[k];
	  ComplexVector& TmpSource = vSources[k];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    TmpDestination[i] += this->HamiltonianShift * TmpSource[i];
	}
    }
  return vDestinations;
}

// initialize the temporary arrays
//

void TripleTensorProductSparseMatrixHamiltonian::InitializeTemporaryArrays()
{
  int RightMatrixDimension = this->RightMatrices[0].GetNbrRow();
  int LeftMatrixDimension = this->LeftMatrices[0].GetNbrRow();
  this->TemporaryArray = new double*[RightMatrixDimension];
  for (int i = 0; i < RightMatrixDimension; ++i)
    this->TemporaryArray[i] = new double[LeftMatrixDimension];
  this->ComplexTemporaryArray = new Complex*[RightMatrixDimension];
  for (int i = 0; i < RightMatrixDimension; ++i)
    this->ComplexTemporaryArray[i] = new Complex[LeftMatrixDimension];
}

// core part of the tensor-multiplication
//
// tensorIndex = index of tensore to consider
// localTemporaryArray = temporary array used to store the partial multiplication
// vSource = vector to be multiplied
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate

void TripleTensorProductSparseMatrixHamiltonian::LowLevelAddMultiplyTensorCore(int tensorIndex, double** localTemporaryArray, RealVector& vSource, 
									       int firstComponent, int nbrComponent)
{
  int IndexStep = this->RightMatrices[tensorIndex].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent - 1;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  int RightMatrixDimension = this->RightMatrices[tensorIndex].GetNbrRow();
  int LeftMatrixDimension = this->LeftMatrices[tensorIndex].GetNbrRow();
  SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[tensorIndex];
  SparseRealMatrix& TmpRightMatrix = this->RightMatrices[tensorIndex];
  for (int j = firstComponent; j <= LastComponent; ++j)
    {
      TmpBRowPointer = TmpRightMatrix.RowPointers[j];
      if (TmpBRowPointer >= 0l)
	{
	  double* Tmp2 = localTemporaryArray[j];
	  for (int k = 0; k < LeftMatrixDimension; ++k)
	    Tmp2[k] = 0.0;
	  TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[j];
	  for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
	    {
	      double Tmp = TmpRightMatrix.MatrixElements[l] * this->Coefficients[tensorIndex];
	      int TmpIndex = TmpRightMatrix.ColumnIndices[l];
	      for (int k = 0; k < LeftMatrixDimension; ++k)
		{
		  Tmp2[k] += Tmp * vSource[k * IndexStep + TmpIndex];
		}
	    }
	}
    }
}

// core part of the tensor-multiplication
//
// tensorIndex = index of tensore to consider
// localTemporaryArray = temporary array used to store the partial multiplication
// vSource = vector to be multiplied
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate

void TripleTensorProductSparseMatrixHamiltonian::LowLevelAddMultiplyTensorCore(int tensorIndex, Complex** localTemporaryArray, ComplexVector& vSource, 
									       int firstComponent, int nbrComponent)
{
  int IndexStep = this->RightMatrices[tensorIndex].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent - 1;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  int RightMatrixDimension = this->RightMatrices[tensorIndex].GetNbrRow();
  int LeftMatrixDimension = this->LeftMatrices[tensorIndex].GetNbrRow();
  SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[tensorIndex];
  SparseRealMatrix& TmpRightMatrix = this->RightMatrices[tensorIndex];
  for (int j = firstComponent; j <= LastComponent; ++j)
    {
      TmpBRowPointer = TmpRightMatrix.RowPointers[j];
      if (TmpBRowPointer >= 0l)
	{
	  Complex* Tmp2 = localTemporaryArray[j];
	  for (int k = 0; k < LeftMatrixDimension; ++k)
	    Tmp2[k] = 0.0;
	  TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[j];
	  for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
	    {
	      double Tmp = TmpRightMatrix.MatrixElements[l] * this->Coefficients[tensorIndex];
	      int TmpIndex = TmpRightMatrix.ColumnIndices[l];
	      for (int k = 0; k < LeftMatrixDimension; ++k)
		{
		  Tmp2[k] += Tmp * vSource[k * IndexStep + TmpIndex];
		}
	    }
	}
    }
}

// core part of the tensor-multiplication (second part computing the final result for one tensor product)
//
// tensorIndex = index of tensore to consider
// localTemporaryArray = temporary array used to store the partial multiplication
// vDestination = vector where the result will be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate

void TripleTensorProductSparseMatrixHamiltonian::LowLevelAddMultiplyTensorCoreDestination(int tensorIndex, double** localTemporaryArray, 
											  RealVector& vDestination, 
											  int firstComponent, int nbrComponent)
{
  int IndexStep = this->RightMatrices[0].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent - 1;
  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  int RightMatrixDimension = this->RightMatrices[tensorIndex].GetNbrRow();
  int LeftMatrixDimension = this->LeftMatrices[tensorIndex].GetNbrRow();
  SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[tensorIndex];
  SparseRealMatrix& TmpRightMatrix = this->RightMatrices[tensorIndex];
  for (int LeftMatrixStartingIndex = 0; LeftMatrixStartingIndex < LeftMatrixDimension; ++LeftMatrixStartingIndex)
    {
      TmpARowPointer = TmpLeftMatrix.RowPointers[LeftMatrixStartingIndex];
      if (TmpARowPointer >= 0l)
	{
	  TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[LeftMatrixStartingIndex];
	  for (int RightMatrixStartingIndex = firstComponent; RightMatrixStartingIndex <= LastComponent; ++RightMatrixStartingIndex)
	    {
	      if (TmpRightMatrix.RowPointers[RightMatrixStartingIndex] >= 0)
		{
		  double Tmp = 0.0;
		  double* Tmp2 = localTemporaryArray[RightMatrixStartingIndex];
		  for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
		    {
		      Tmp += TmpLeftMatrix.MatrixElements[k] * Tmp2[TmpLeftMatrix.ColumnIndices[k]];
		    }
		  vDestination[LeftMatrixStartingIndex * IndexStep + RightMatrixStartingIndex] += Tmp;
		}
	    }
	}
    }
}

// core part of the tensor-multiplication (second part computing the final result for one tensor product)
//
// tensorIndex = index of tensore to consider
// localTemporaryArray = temporary array used to store the partial multiplication
// vDestination = vector where the result will be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate

void TripleTensorProductSparseMatrixHamiltonian::LowLevelAddMultiplyTensorCoreDestination(int tensorIndex, Complex** localTemporaryArray, 
											  ComplexVector& vDestination, 
											  int firstComponent, int nbrComponent)
{
  int IndexStep = this->RightMatrices[0].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent - 1;
  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  int RightMatrixDimension = this->RightMatrices[tensorIndex].GetNbrRow();
  int LeftMatrixDimension = this->LeftMatrices[tensorIndex].GetNbrRow();
  SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[tensorIndex];
  SparseRealMatrix& TmpRightMatrix = this->RightMatrices[tensorIndex];
  for (int LeftMatrixStartingIndex = 0; LeftMatrixStartingIndex < LeftMatrixDimension; ++LeftMatrixStartingIndex)
    {
      TmpARowPointer = TmpLeftMatrix.RowPointers[LeftMatrixStartingIndex];
      if (TmpARowPointer >= 0l)
	{
	  TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[LeftMatrixStartingIndex];
	  for (int RightMatrixStartingIndex = firstComponent; RightMatrixStartingIndex <= LastComponent; ++RightMatrixStartingIndex)
	    {
	      if (TmpRightMatrix.RowPointers[RightMatrixStartingIndex] >= 0)
		{
		  Complex Tmp = 0.0;
		  Complex* Tmp2 = localTemporaryArray[RightMatrixStartingIndex];
		  for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
		    {
		      Tmp += TmpLeftMatrix.MatrixElements[k] * Tmp2[TmpLeftMatrix.ColumnIndices[k]];
		    }
		  vDestination[LeftMatrixStartingIndex * IndexStep + RightMatrixStartingIndex] += Tmp;
		}
	    }
	}
    }
}
