////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//  class of hamiltonian defined as a linear combination of tensor products   //
//                                                                            //
//                        last modification : 08/11/2012                      //
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


#include "Hamiltonian/MultipleTensorProductSparseMatrixHamiltonian.h"
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

MultipleTensorProductSparseMatrixHamiltonian::MultipleTensorProductSparseMatrixHamiltonian()
{
  this->Matrices = 0;
}

// contructor 
//
// nbrTensorProducts = number of multiple tensor products whose linear combination defined the Hamiltonian 
// nbrMatrixPerTensorProducts = number of matrices within each multiple tensor products
// matrices = matrices of each multiple tensor product (first entry) and whithin each multiple tensor product (second entry)
// coefficients = coefficients of the ensor product linear combination
// architecture = architecture to use for precalculation

MultipleTensorProductSparseMatrixHamiltonian::MultipleTensorProductSparseMatrixHamiltonian(int nbrTensorProducts, int nbrMatrixPerTensorProducts, SparseRealMatrix** matrices,
											   double* coefficients, AbstractArchitecture* architecture)
{
  this->NbrTensorProducts = nbrTensorProducts;
  this->NbrMatrixPerTensorProducts = nbrMatrixPerTensorProducts;
  this->Matrices = new SparseRealMatrix*[this->NbrTensorProducts];
  this->Coefficients = new double[this->NbrTensorProducts];
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      this->Coefficients[i] = coefficients[i];
      this->Matrices[i] = new SparseRealMatrix[this->NbrMatrixPerTensorProducts];
      for (int j = 0; j < this->NbrMatrixPerTensorProducts; ++j)
	{
	  this->Matrices[i][j] = matrices[i][j];
	}
    }
  this->HamiltonianShift = 0.0;
  this->Architecture = architecture;
  long HamiltonianDimension = this->Matrices[0][0].GetNbrRow();
  for (int i = 1; i < this->NbrMatrixPerTensorProducts; ++i)
    {
      HamiltonianDimension *=  this->Matrices[0][i].GetNbrRow();
    }
  this->HilbertSpace = new UndescribedHilbertSpace(HamiltonianDimension);
  this->LeftHamiltonianVectorMultiplicationFlag = true;
  this->InitializeTemporaryArrays();
}

// destructor
//

MultipleTensorProductSparseMatrixHamiltonian::~MultipleTensorProductSparseMatrixHamiltonian() 
{
  if (this->Matrices != 0)
    {
      int TmpMatrixDimension = this->Matrices[0][0].GetNbrRow();
      for (int j = 0; j < TmpMatrixDimension; ++j)
	delete[] this->TemporaryArray[j];
      delete[] this->TemporaryArray;
      for (int j = 0; j < TmpMatrixDimension; ++j)
	delete[] this->ComplexTemporaryArray[j];
      delete[] this->ComplexTemporaryArray;
      for (int i = 0; i < this->NbrTensorProducts; ++i)
	{
	  delete[] this->Matrices[i];
	}
      delete[] this->Matrices;
      delete[] this->Coefficients;

    }
  delete this->HilbertSpace;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void MultipleTensorProductSparseMatrixHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace) 
{
  this->HilbertSpace = hilbertSpace;
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* MultipleTensorProductSparseMatrixHamiltonian::GetHilbertSpace ()
{
  return this->HilbertSpace;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int MultipleTensorProductSparseMatrixHamiltonian::GetHilbertSpaceDimension () 
{
  return this->HilbertSpace->GetHilbertSpaceDimension();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value

void MultipleTensorProductSparseMatrixHamiltonian::ShiftHamiltonian (double shift) 
{
  this->HamiltonianShift = shift;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& MultipleTensorProductSparseMatrixHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									      int firstComponent, int nbrComponent)
{
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

RealVector* MultipleTensorProductSparseMatrixHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
//   int IndexStep = this->RightMatrices[0].GetNbrColumn();
//   int LastComponent = firstComponent + nbrComponent - 1;
//   int AMatrixLastIndex = LastComponent / this->RightMatrices[0].GetNbrRow();
//   int BMatrixLastIndex = LastComponent % this->RightMatrices[0].GetNbrRow();
//   long TmpARowPointer;
//   long TmpARowLastPointer;
//   long TmpBRowPointer;
//   long TmpBRowLastPointer;
//   double* Tmp = new double[nbrVectors];
//   for (int i = 0; i < this->NbrTensorProducts; ++i)
//     {
//       SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
//       SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
//       int AMatrixStartingIndex = firstComponent / this->RightMatrices[0].GetNbrRow();
//       int BMatrixStartingIndex = firstComponent % this->RightMatrices[0].GetNbrRow();
//       int TotalIndex = firstComponent;
//       for (; AMatrixStartingIndex <=  AMatrixLastIndex; ++AMatrixStartingIndex)
// 	{
// 	  TmpARowPointer = TmpLeftMatrix.RowPointers[AMatrixStartingIndex];
// 	  if (TmpARowPointer >= 0l)
// 	    {
// 	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[AMatrixStartingIndex];
// 	      int TmpBMatrixLastIndex = TmpRightMatrix.GetNbrRow() - 1;
// 	      if (AMatrixStartingIndex == (AMatrixLastIndex - 1))
// 		TmpBMatrixLastIndex = BMatrixLastIndex;
// 	      for (; BMatrixStartingIndex <=  TmpBMatrixLastIndex; ++BMatrixStartingIndex)
// 		{
// 		  TmpBRowPointer = TmpRightMatrix.RowPointers[BMatrixStartingIndex];
// 		  if (TmpBRowPointer >= 0l)
// 		    {
// 		      TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[BMatrixStartingIndex];
// 		      for (int l = 0; l < nbrVectors; ++l)
// 			Tmp[l] = 0.0;
// 		      for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
// 			{
// 			  double Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
// 			  int TmpIndex = TmpLeftMatrix.ColumnIndices[k] * IndexStep;
// 			  for (long j = TmpBRowPointer; j <= TmpBRowLastPointer; ++j)
// 			    {
// 			      int InputIndex = TmpIndex + TmpRightMatrix.ColumnIndices[j];
// 			      for (int l = 0; l < nbrVectors; ++l)			      
// 				Tmp[l] += Tmp2 * TmpRightMatrix.MatrixElements[j] * vSources[l][InputIndex];
// 			    }
// 			}
// 		      int OutputIndex = AMatrixStartingIndex * IndexStep + BMatrixStartingIndex;
// 		      for (int l = 0; l < nbrVectors; ++l)
// 			vDestinations[l][OutputIndex] += Tmp[l];
// 		    }
// 		}
// 	    }
// 	  BMatrixStartingIndex = 0;
// 	}
//     }
//   delete[] Tmp;
//   if (this->HamiltonianShift != 0.0)
//     {
//       for (int k= 0; k < nbrVectors; ++k)
// 	{
// 	  RealVector& TmpDestination = vDestinations[k];
// 	  RealVector& TmpSource = vSources[k];
// 	  for (int i = firstComponent; i < LastComponent; ++i)
// 	    TmpDestination[i] += this->HamiltonianShift * TmpSource[i];
// 	}
//     }
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

ComplexVector& MultipleTensorProductSparseMatrixHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									 int firstComponent, int nbrComponent)
{
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

ComplexVector* MultipleTensorProductSparseMatrixHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
//   int IndexStep = this->RightMatrices[0].GetNbrColumn();
//   int LastComponent = firstComponent + nbrComponent - 1;
//   int AMatrixLastIndex = LastComponent / this->RightMatrices[0].GetNbrRow();
//   int BMatrixLastIndex = LastComponent % this->RightMatrices[0].GetNbrRow();
//   long TmpARowPointer;
//   long TmpARowLastPointer;
//   long TmpBRowPointer;
//   long TmpBRowLastPointer;
//   Complex* Tmp = new Complex[nbrVectors];
//   for (int i = 0; i < this->NbrTensorProducts; ++i)
//     {
//       SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
//       SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
//       int AMatrixStartingIndex = firstComponent / this->RightMatrices[0].GetNbrRow();
//       int BMatrixStartingIndex = firstComponent % this->RightMatrices[0].GetNbrRow();
//       int TotalIndex = firstComponent;
//       for (; AMatrixStartingIndex <=  AMatrixLastIndex; ++AMatrixStartingIndex)
// 	{
// 	  TmpARowPointer = TmpLeftMatrix.RowPointers[AMatrixStartingIndex];
// 	  if (TmpARowPointer >= 0l)
// 	    {
// 	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[AMatrixStartingIndex];
// 	      int TmpBMatrixLastIndex = TmpRightMatrix.GetNbrRow() - 1;
// 	      if (AMatrixStartingIndex == (AMatrixLastIndex - 1))
// 		TmpBMatrixLastIndex = BMatrixLastIndex;
// 	      for (; BMatrixStartingIndex <=  TmpBMatrixLastIndex; ++BMatrixStartingIndex)
// 		{
// 		  TmpBRowPointer = TmpRightMatrix.RowPointers[BMatrixStartingIndex];
// 		  if (TmpBRowPointer >= 0l)
// 		    {
// 		      TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[BMatrixStartingIndex];
// 		      for (int l = 0; l < nbrVectors; ++l)
// 			Tmp[l] = 0.0;
// 		      for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
// 			{
// 			  Complex Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
// 			  int TmpIndex = TmpLeftMatrix.ColumnIndices[k] * IndexStep;
// 			  for (long j = TmpBRowPointer; j <= TmpBRowLastPointer; ++j)
// 			    {
// 			      int InputIndex = TmpIndex + TmpRightMatrix.ColumnIndices[j];
// 			      for (int l = 0; l < nbrVectors; ++l)			      
// 				Tmp[l] += Tmp2 * TmpRightMatrix.MatrixElements[j] * vSources[l][InputIndex];
// 			    }
// 			}
// 		      int OutputIndex = AMatrixStartingIndex * IndexStep + BMatrixStartingIndex;
// 		      for (int l = 0; l < nbrVectors; ++l)
// 			vDestinations[l][OutputIndex] += Tmp[l];
// 		    }
// 		}
// 	    }
// 	  BMatrixStartingIndex = 0;
// 	}
//     }
//   delete[] Tmp;
//   if (this->HamiltonianShift != 0.0)
//     {
//       for (int k= 0; k < nbrVectors; ++k)
// 	{
// 	  ComplexVector& TmpDestination = vDestinations[k];
// 	  ComplexVector& TmpSource = vSources[k];
// 	  for (int i = firstComponent; i < LastComponent; ++i)
// 	    TmpDestination[i] += this->HamiltonianShift * TmpSource[i];
// 	}
//     }
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

RealVector& MultipleTensorProductSparseMatrixHamiltonian::ConjugateLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
										       int firstComponent, int nbrComponent)
{
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* MultipleTensorProductSparseMatrixHamiltonian::ConjugateLowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											    int firstComponent, int nbrComponent)
{
  return this->LowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& MultipleTensorProductSparseMatrixHamiltonian::ConjugateLowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										       int firstComponent, int nbrComponent)
{
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* MultipleTensorProductSparseMatrixHamiltonian::ConjugateLowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
											       int firstComponent, int nbrComponent)
{
  return this->LowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
}

// initialize the temporary arrays
//

void MultipleTensorProductSparseMatrixHamiltonian::InitializeTemporaryArrays()
{
//   int RightMatrixDimension = this->RightMatrices[0].GetNbrRow();
//   int LeftMatrixDimension = this->LeftMatrices[0].GetNbrRow();
//   this->TemporaryArray = new double*[RightMatrixDimension];
//   for (int i = 0; i < RightMatrixDimension; ++i)
//     this->TemporaryArray[i] = new double[LeftMatrixDimension];
//   this->ComplexTemporaryArray = new Complex*[RightMatrixDimension];
//   for (int i = 0; i < RightMatrixDimension; ++i)
//     this->ComplexTemporaryArray[i] = new Complex[LeftMatrixDimension];
}

// core part of the tensor-multiplication
//
// tensorIndex = index of tensore to consider
// localTemporaryArray = temporary array used to store the partial multiplication
// vSource = vector to be multiplied
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate

void MultipleTensorProductSparseMatrixHamiltonian::LowLevelAddMultiplyTensorCore(int tensorIndex, double** localTemporaryArray, RealVector& vSource, 
										 int firstComponent, int nbrComponent)
{
//   int IndexStep = this->RightMatrices[tensorIndex].GetNbrColumn();
//   int LastComponent = firstComponent + nbrComponent - 1;
//   long TmpBRowPointer;
//   long TmpBRowLastPointer;
//   int RightMatrixDimension = this->RightMatrices[tensorIndex].GetNbrRow();
//   int LeftMatrixDimension = this->LeftMatrices[tensorIndex].GetNbrRow();
//   SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[tensorIndex];
//   SparseRealMatrix& TmpRightMatrix = this->RightMatrices[tensorIndex];
//   for (int j = firstComponent; j <= LastComponent; ++j)
//     {
//       TmpBRowPointer = TmpRightMatrix.RowPointers[j];
//       if (TmpBRowPointer >= 0l)
// 	{
// 	  double* Tmp2 = localTemporaryArray[j];
// 	  for (int k = 0; k < LeftMatrixDimension; ++k)
// 	    Tmp2[k] = 0.0;
// 	  TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[j];
// 	  for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
// 	    {
// 	      double Tmp = TmpRightMatrix.MatrixElements[l] * this->Coefficients[tensorIndex];
// 	      int TmpIndex = TmpRightMatrix.ColumnIndices[l];
// 	      for (int k = 0; k < LeftMatrixDimension; ++k)
// 		{
// 		  Tmp2[k] += Tmp * vSource[k * IndexStep + TmpIndex];
// 		}
// 	    }
// 	}
//     }
}

// core part of the tensor-multiplication
//
// tensorIndex = index of tensore to consider
// localTemporaryArray = temporary array used to store the partial multiplication
// vSource = vector to be multiplied
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate

void MultipleTensorProductSparseMatrixHamiltonian::LowLevelAddMultiplyTensorCore(int tensorIndex, Complex** localTemporaryArray, ComplexVector& vSource, 
										 int firstComponent, int nbrComponent)
{
//   int IndexStep = this->RightMatrices[tensorIndex].GetNbrColumn();
//   int LastComponent = firstComponent + nbrComponent - 1;
//   long TmpBRowPointer;
//   long TmpBRowLastPointer;
//   int RightMatrixDimension = this->RightMatrices[tensorIndex].GetNbrRow();
//   int LeftMatrixDimension = this->LeftMatrices[tensorIndex].GetNbrRow();
//   SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[tensorIndex];
//   SparseRealMatrix& TmpRightMatrix = this->RightMatrices[tensorIndex];
//   for (int j = firstComponent; j <= LastComponent; ++j)
//     {
//       TmpBRowPointer = TmpRightMatrix.RowPointers[j];
//       if (TmpBRowPointer >= 0l)
// 	{
// 	  Complex* Tmp2 = localTemporaryArray[j];
// 	  for (int k = 0; k < LeftMatrixDimension; ++k)
// 	    Tmp2[k] = 0.0;
// 	  TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[j];
// 	  for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
// 	    {
// 	      double Tmp = TmpRightMatrix.MatrixElements[l] * this->Coefficients[tensorIndex];
// 	      int TmpIndex = TmpRightMatrix.ColumnIndices[l];
// 	      for (int k = 0; k < LeftMatrixDimension; ++k)
// 		{
// 		  Tmp2[k] += Tmp * vSource[k * IndexStep + TmpIndex];
// 		}
// 	    }
// 	}
//     }
}

// core part of the tensor-multiplication (second part computing the final result for one tensor product)
//
// tensorIndex = index of tensore to consider
// localTemporaryArray = temporary array used to store the partial multiplication
// vDestination = vector where the result will be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate

void MultipleTensorProductSparseMatrixHamiltonian::LowLevelAddMultiplyTensorCoreDestination(int tensorIndex, double** localTemporaryArray, 
											    RealVector& vDestination, 
											    int firstComponent, int nbrComponent)
{
//   int IndexStep = this->RightMatrices[0].GetNbrColumn();
//   int LastComponent = firstComponent + nbrComponent - 1;
//   long TmpARowPointer;
//   long TmpARowLastPointer;
//   long TmpBRowPointer;
//   long TmpBRowLastPointer;
//   int RightMatrixDimension = this->RightMatrices[tensorIndex].GetNbrRow();
//   int LeftMatrixDimension = this->LeftMatrices[tensorIndex].GetNbrRow();
//   SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[tensorIndex];
//   SparseRealMatrix& TmpRightMatrix = this->RightMatrices[tensorIndex];
//   for (int LeftMatrixStartingIndex = 0; LeftMatrixStartingIndex < LeftMatrixDimension; ++LeftMatrixStartingIndex)
//     {
//       TmpARowPointer = TmpLeftMatrix.RowPointers[LeftMatrixStartingIndex];
//       if (TmpARowPointer >= 0l)
// 	{
// 	  TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[LeftMatrixStartingIndex];
// 	  for (int RightMatrixStartingIndex = firstComponent; RightMatrixStartingIndex <= LastComponent; ++RightMatrixStartingIndex)
// 	    {
// 	      if (TmpRightMatrix.RowPointers[RightMatrixStartingIndex] >= 0)
// 		{
// 		  double Tmp = 0.0;
// 		  double* Tmp2 = localTemporaryArray[RightMatrixStartingIndex];
// 		  for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
// 		    {
// 		      Tmp += TmpLeftMatrix.MatrixElements[k] * Tmp2[TmpLeftMatrix.ColumnIndices[k]];
// 		    }
// 		  vDestination[LeftMatrixStartingIndex * IndexStep + RightMatrixStartingIndex] += Tmp;
// 		}
// 	    }
// 	}
//     }
}

// core part of the tensor-multiplication (second part computing the final result for one tensor product)
//
// tensorIndex = index of tensore to consider
// localTemporaryArray = temporary array used to store the partial multiplication
// vDestination = vector where the result will be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate

void MultipleTensorProductSparseMatrixHamiltonian::LowLevelAddMultiplyTensorCoreDestination(int tensorIndex, Complex** localTemporaryArray, 
											    ComplexVector& vDestination, 
											    int firstComponent, int nbrComponent)
{
//   int IndexStep = this->RightMatrices[0].GetNbrColumn();
//   int LastComponent = firstComponent + nbrComponent - 1;
//   long TmpARowPointer;
//   long TmpARowLastPointer;
//   long TmpBRowPointer;
//   long TmpBRowLastPointer;
//   int RightMatrixDimension = this->RightMatrices[tensorIndex].GetNbrRow();
//   int LeftMatrixDimension = this->LeftMatrices[tensorIndex].GetNbrRow();
//   SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[tensorIndex];
//   SparseRealMatrix& TmpRightMatrix = this->RightMatrices[tensorIndex];
//   for (int LeftMatrixStartingIndex = 0; LeftMatrixStartingIndex < LeftMatrixDimension; ++LeftMatrixStartingIndex)
//     {
//       TmpARowPointer = TmpLeftMatrix.RowPointers[LeftMatrixStartingIndex];
//       if (TmpARowPointer >= 0l)
// 	{
// 	  TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[LeftMatrixStartingIndex];
// 	  for (int RightMatrixStartingIndex = firstComponent; RightMatrixStartingIndex <= LastComponent; ++RightMatrixStartingIndex)
// 	    {
// 	      if (TmpRightMatrix.RowPointers[RightMatrixStartingIndex] >= 0)
// 		{
// 		  Complex Tmp = 0.0;
// 		  Complex* Tmp2 = localTemporaryArray[RightMatrixStartingIndex];
// 		  for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
// 		    {
// 		      Tmp += TmpLeftMatrix.MatrixElements[k] * Tmp2[TmpLeftMatrix.ColumnIndices[k]];
// 		    }
// 		  vDestination[LeftMatrixStartingIndex * IndexStep + RightMatrixStartingIndex] += Tmp;
// 		}
// 	    }
// 	}
//     }
}
