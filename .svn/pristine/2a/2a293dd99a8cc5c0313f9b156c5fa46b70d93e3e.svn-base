////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//  class of hamiltonian defined as a linear combination of tensor products   //
//               focusing on a single block of  tensor product                //
//                                                                            //
//                        last modification : 07/01/2013                      //
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


#include "Hamiltonian/TensorProductSparseMatrixSelectedBlockHamiltonian.h"
#include "MathTools/Complex.h" 
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/GenericHamiltonianPrecalculationOperation.h"
#include "Architecture/ArchitectureOperation/VectorTensorMultiplicationCoreOperation.h"
#include "Architecture/ArchitectureOperation/VectorSparseTensorMultiplyOperation.h"


#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;


// contructor 
//
// nbrTensorProducts = number of tensor products whose linear combination defined the Hamiltonian 
// leftMatrices = left matrices of each tensor product
// rightMatrices = right matrices of each tensor product
// coefficients = coefficients of the ensor product linear combination
// blockSize = number of indices in the selected block
// blockIndices = pairs of indices (for resp. the left and right matrix) that define the selected block 
// architecture = architecture to use for precalculation
// memory = amount of memory that can be used for precalculations (in bytes)

TensorProductSparseMatrixSelectedBlockHamiltonian::TensorProductSparseMatrixSelectedBlockHamiltonian(int nbrTensorProducts, SparseRealMatrix* leftMatrices,  
												     SparseRealMatrix* rightMatrices, double* coefficients,
												     int blockSize, long* blockIndices, AbstractArchitecture* architecture, long memory)
{
  this->NbrTensorProducts = nbrTensorProducts;
  this->LeftMatrices = new SparseRealMatrix[this->NbrTensorProducts];
  this->RightMatrices = new SparseRealMatrix[this->NbrTensorProducts];
  this->Coefficients = new double[this->NbrTensorProducts];
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      this->LeftMatrices[i] = leftMatrices[i];
      this->RightMatrices[i] = rightMatrices[i];
      this->Coefficients[i] = coefficients[i];
    }
  this->HamiltonianShift = 0.0;
  this->Architecture = architecture;
  this->HilbertSpace = new UndescribedHilbertSpace(blockSize);
  this->LeftHamiltonianVectorMultiplicationFlag = true;
  this->BlockIndices = blockIndices;
  long TmpMatrixSize = this->LeftMatrices[0].GetNbrRow();
  this->InitializeTemporaryArrays();
  this->ExternalBlockIndexProductTable = false;
  this->BlockIndexProductTable = new long* [TmpMatrixSize];
  this->BlockIndexProductTableNbrElements = new int [TmpMatrixSize];
  this->BlockIndexProductTableShift = new int [TmpMatrixSize];
  long* TmpBlockIndices = new long [TmpMatrixSize];
  for (long i = 0; i < TmpMatrixSize; ++i)
    {
      this->BlockIndexProductTableNbrElements[i] = 0;
      this->BlockIndexProductTableShift[i] = -1;
      for (int j = 0; j < blockSize; ++j)
	{
	  if ((this->BlockIndices[j] / TmpMatrixSize) == i)
	    {
	      if (this->BlockIndexProductTableShift[i] < 0)
		this->BlockIndexProductTableShift[i] = j;
	      TmpBlockIndices[this->BlockIndexProductTableNbrElements[i]] = this->BlockIndices[j];
	      ++this->BlockIndexProductTableNbrElements[i];
	    }
	}
      this->BlockIndexProductTable[i] = new long[this->BlockIndexProductTableNbrElements[i]];
      for (int j = 0; j < this->BlockIndexProductTableNbrElements[i]; ++j)
	{
	  this->BlockIndexProductTable[i][j] = TmpBlockIndices[j];
	}
    }
  delete[] TmpBlockIndices;
  this->BlockSize = blockSize;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  if (memory > 0l)
    {
      this->EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
      this->TemporaryRowPointers = new long[this->EffectiveHilbertSpaceDimension];
      timeval TotalStartingTime;
      timeval TotalEndingTime;
      gettimeofday (&(TotalStartingTime), 0);
      long NbrNonZeroMatrixElements = this->FastMultiplicationMemory() >> 3;
      gettimeofday (&(TotalEndingTime), 0);
      double DTime = ((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
		      ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));
      cout << "nbr non-zero matrix elements = " << NbrNonZeroMatrixElements << " (done in " << DTime << "s)" <<  endl;
      if (memory > (NbrNonZeroMatrixElements << 3))
	{
	  this->TemporaryRowLastPointers = new long[this->EffectiveHilbertSpaceDimension];
	  this->FastMultiplicationFlag = true;
	  long TmpPointer = 0;
	  for (int i = 0; i < this->EffectiveHilbertSpaceDimension; ++i)
	    {
	      long Tmp = this->TemporaryRowPointers[i]; 
	      if (this->TemporaryRowPointers[i] > 0l)
		{     
		  this->TemporaryRowPointers[i] = TmpPointer;
		  this->TemporaryRowLastPointers[i] = TmpPointer + Tmp - 1l;
		}
	      else
		{
		  this->TemporaryRowPointers[i] = -1l;
		  this->TemporaryRowLastPointers[i] = -1l;
		}
	      TmpPointer += Tmp;
	    }
	  this->TemporaryMatrixElements = new double[NbrNonZeroMatrixElements];
	  this->TemporaryMatrixColumnIndices = new int[NbrNonZeroMatrixElements];
	  gettimeofday (&(TotalStartingTime), 0);
	  this->EnableFastMultiplication(); 
	  gettimeofday (&(TotalEndingTime), 0);
	  DTime = ((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
		   ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));
	  cout << "using ";
	  PrintMemorySize(cout, (NbrNonZeroMatrixElements << 3)) << " for hamiltonian precalculations (done in " << DTime << "s)" << endl;
	}
      else
	{
	  cout << "cannot use hamiltonian precalculations, ";
	  PrintMemorySize(cout, (NbrNonZeroMatrixElements << 3)) << " are required" << endl;
	  this->FastMultiplicationFlag = false;
	  delete[] this->TemporaryRowPointers;
	}
    }
  else
    {
      this->FastMultiplicationFlag = false;
    }
}

// contructor providing an efficient block index scheme
//
// nbrTensorProducts = number of tensor products whose linear combination defined the Hamiltonian 
// leftMatrices = left matrices of each tensor product
// rightMatrices = right matrices of each tensor product
// coefficients = coefficients of the ensor product linear combination
// blockSize = number of indices in the selected block
// blockIndices = pairs of indices (for resp. the left and right matrix) that define the selected block 
// blockIndexProductTable = table that contains all the linearized indices attached to one left matrix index
// blockIndexProductTableNbrElements = number of  linearized indices attached to one left matrix index
// blockIndexProductTableShift = first index in the hilbert space where a given left matrix index occurs
// architecture = architecture to use for precalculation
// memory = amount of memory that can be used for precalculations (in bytes)

TensorProductSparseMatrixSelectedBlockHamiltonian::TensorProductSparseMatrixSelectedBlockHamiltonian(int nbrTensorProducts, SparseRealMatrix* leftMatrices,  
												     SparseRealMatrix* rightMatrices, double* coefficients,
												     int blockSize, long* blockIndices, 
												     long** blockIndexProductTable, int* blockIndexProductTableNbrElements,
												     int* blockIndexProductTableShift,
												     AbstractArchitecture* architecture, long memory)
{
  this->NbrTensorProducts = nbrTensorProducts;
  this->LeftMatrices = new SparseRealMatrix[this->NbrTensorProducts];
  this->RightMatrices = new SparseRealMatrix[this->NbrTensorProducts];
  this->Coefficients = new double[this->NbrTensorProducts];
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      this->LeftMatrices[i] = leftMatrices[i];
      this->RightMatrices[i] = rightMatrices[i];
      this->Coefficients[i] = coefficients[i];
      cout << "Sparsity of the left matrix " << i << " = " << ((leftMatrices[i].ComputeNbrNonZeroMatrixElements() * 100.0) / (((double) leftMatrices[i].GetNbrRow()) * ((double) leftMatrices[i].GetNbrRow()))) << "%" << endl;
      
    }
  this->HamiltonianShift = 0.0;
  this->Architecture = architecture;
  this->HilbertSpace = new UndescribedHilbertSpace(blockSize);
  this->LeftHamiltonianVectorMultiplicationFlag = true;
  this->BlockIndices = blockIndices;
  this->BlockSize = blockSize;
  long TmpMatrixSize = this->LeftMatrices[0].GetNbrRow();
  this->InitializeTemporaryArrays();
  long FullEMatrixSize = ((long) TmpMatrixSize) *  ((long) TmpMatrixSize);
  this->InvertBlockIndices = new int [FullEMatrixSize];
  for (long i = 0l; i < FullEMatrixSize; ++i)
    this->InvertBlockIndices[i] = -1;
  for (int i = 0; i < this->BlockSize; ++i)
    this->InvertBlockIndices[this->BlockIndices[i]] = i;
  this->ExternalBlockIndexProductTable = true;
  this->BlockIndexProductTable = blockIndexProductTable;
  this->BlockIndexProductTableNbrElements = blockIndexProductTableNbrElements;
  this->BlockIndexProductTableShift = blockIndexProductTableShift;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  if (memory > 0l)
    {
      this->EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
      this->TemporaryRowPointers = new long[this->EffectiveHilbertSpaceDimension];
      timeval TotalStartingTime;
      timeval TotalEndingTime;
      gettimeofday (&(TotalStartingTime), 0);
      long NbrNonZeroMatrixElements = this->FastMultiplicationMemory() >> 3;
      gettimeofday (&(TotalEndingTime), 0);
      double DTime = ((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
		      ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));
      cout << "nbr non-zero matrix elements = " << NbrNonZeroMatrixElements << " (done in " << DTime << "s)" <<  endl;
      if (memory > (NbrNonZeroMatrixElements << 3))
	{
	  this->TemporaryRowLastPointers = new long[this->EffectiveHilbertSpaceDimension];
	  this->FastMultiplicationFlag = true;
	  long TmpPointer = 0;
	  for (int i = 0; i < this->EffectiveHilbertSpaceDimension; ++i)
	    {
	      long Tmp = this->TemporaryRowPointers[i]; 
	      if (this->TemporaryRowPointers[i] > 0l)
		{     
		  this->TemporaryRowPointers[i] = TmpPointer;
		  this->TemporaryRowLastPointers[i] = TmpPointer + Tmp - 1l;
		}
	      else
		{
		  this->TemporaryRowPointers[i] = -1l;
		  this->TemporaryRowLastPointers[i] = -1l;
		}
	      TmpPointer += Tmp;
	    }
	  this->TemporaryMatrixElements = new double[NbrNonZeroMatrixElements];
	  this->TemporaryMatrixColumnIndices = new int[NbrNonZeroMatrixElements];
	  gettimeofday (&(TotalStartingTime), 0);
	  this->EnableFastMultiplication(); 
	  gettimeofday (&(TotalEndingTime), 0);
	  DTime = ((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
		   ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));
	  cout << "using ";
	  PrintMemorySize(cout, (NbrNonZeroMatrixElements << 3)) << " for hamiltonian precalculations (done in " << DTime << "s)" << endl;
	}
      else
	{
	  cout << "cannot use hamiltonian precalculations, ";
	  PrintMemorySize(cout, (NbrNonZeroMatrixElements << 3)) << " are required" << endl;
	  this->FastMultiplicationFlag = false;
	  delete[] this->TemporaryRowPointers;
	}
    }
  else
    {
      this->FastMultiplicationFlag = false;
    }
}

// destructor
//

TensorProductSparseMatrixSelectedBlockHamiltonian::~TensorProductSparseMatrixSelectedBlockHamiltonian() 
{
  long TmpMatrixSize = this->LeftMatrices[0].GetNbrRow();
  if (this->ExternalBlockIndexProductTable == false)
    {
      for (long i = 0; i < TmpMatrixSize; ++i)
	if (this->BlockIndexProductTableNbrElements[i] > 0)
	  delete[] this->BlockIndexProductTable[i];
      delete[] this->BlockIndexProductTable;
      delete[] this->BlockIndexProductTableNbrElements;
      delete[] this->BlockIndexProductTableShift;
    }
  if (this->FastMultiplicationFlag == true)
    {
      delete[] this->TemporaryMatrixElements;
      delete[] this->TemporaryMatrixColumnIndices;
      delete[] this->TemporaryRowPointers;
      delete[] this->TemporaryRowLastPointers;
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

RealVector& TensorProductSparseMatrixSelectedBlockHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
										   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == true)
    {
      for (int j = firstComponent; j < LastComponent; ++j)
	{
	  long TmpRowPointer = this->TemporaryRowPointers[j - this->PrecalculationShift];
	  if (TmpRowPointer >= 0l)
	    {
	      long TmpRowLastPointer = this->TemporaryRowLastPointers[j - this->PrecalculationShift];
	      double Tmp = 0.0;
	      for (; TmpRowPointer <= TmpRowLastPointer; ++TmpRowPointer)
		{
		  Tmp += (this->TemporaryMatrixElements[TmpRowPointer] 
			  * vSource[this->TemporaryMatrixColumnIndices[TmpRowPointer]]);
		}		
	      vDestination[j] += Tmp;
	    }
	}
      return vDestination;
    }

  int IndexStep = this->RightMatrices[0].GetNbrColumn();
  int LeftMatrixLastIndex = LastComponent / this->RightMatrices[0].GetNbrRow();
  int RightMatrixLastIndex = LastComponent % this->RightMatrices[0].GetNbrRow();

  int RightMatrixDimension = this->RightMatrices[0].GetNbrRow();
  int LeftMatrixDimension = this->LeftMatrices[0].GetNbrRow();


  double** LocalTemporaryMatrix = this->TemporaryArray;

  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
      SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];

      VectorTensorMultiplicationCoreOperation Operation(this, i, vSource);
      Operation.ApplyOperation(this->Architecture);

      VectorSparseTensorMultiplyOperation Operation2(this, i, &vDestination, true);
      Operation2.ApplyOperation(this->Architecture);

    }

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

RealVector* TensorProductSparseMatrixSelectedBlockHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  int IndexStep = this->LeftMatrices[0].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent;
  int AMatrixLastIndex = LastComponent / this->LeftMatrices[0].GetNbrRow();
  int BMatrixLastIndex = LastComponent % this->LeftMatrices[0].GetNbrRow();
  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  double* Tmp = new double[nbrVectors];
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
      SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
      int AMatrixStartingIndex = firstComponent / this->LeftMatrices[0].GetNbrRow();
      int BMatrixStartingIndex = firstComponent % this->LeftMatrices[0].GetNbrRow();
      int TotalIndex = firstComponent;
      for (; AMatrixStartingIndex <  AMatrixLastIndex; ++AMatrixStartingIndex)
	{
	  TmpARowPointer = TmpLeftMatrix.RowPointers[AMatrixStartingIndex];
	  if (TmpARowPointer >= 0l)
	    {
	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[AMatrixStartingIndex];
	      int TmpBMatrixLastIndex = TmpRightMatrix.GetNbrRow();
	      if (AMatrixStartingIndex == (AMatrixLastIndex - 1))
		TmpBMatrixLastIndex = BMatrixLastIndex;
	      for (; BMatrixStartingIndex <  TmpBMatrixLastIndex; ++BMatrixStartingIndex)
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

ComplexVector& TensorProductSparseMatrixSelectedBlockHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										      int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == true)
    {
      for (int j = firstComponent; j < LastComponent; ++j)
	{
	  long TmpRowPointer = this->TemporaryRowPointers[j];
	  if (TmpRowPointer >= 0l)
	    {
	      long TmpRowLastPointer = this->TemporaryRowLastPointers[j];
	      Complex Tmp = 0.0;
	      for (; TmpRowPointer <= TmpRowLastPointer; ++TmpRowPointer)
		{
		  Tmp += (this->TemporaryMatrixElements[TmpRowPointer] 
			  * vSource[this->TemporaryMatrixColumnIndices[TmpRowPointer]]);
		}		
	      vDestination[j] += Tmp;
	    }
	}
      return vDestination;
    }


  int IndexStep = this->RightMatrices[0].GetNbrColumn();
  int LeftMatrixLastIndex = LastComponent / this->RightMatrices[0].GetNbrRow();
  int RightMatrixLastIndex = LastComponent % this->RightMatrices[0].GetNbrRow();

  int RightMatrixDimension = this->RightMatrices[0].GetNbrRow();
  int LeftMatrixDimension = this->LeftMatrices[0].GetNbrRow();


  Complex** LocalTemporaryMatrix = this->ComplexTemporaryArray;

  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
      SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];

      VectorTensorMultiplicationCoreOperation Operation(this, i, vSource);
      Operation.ApplyOperation(this->Architecture);

      for (int j = firstComponent; j < LastComponent; ++j)
	{
	  TmpARowPointer = TmpLeftMatrix.RowPointers[this->BlockIndices[j] / IndexStep];
	  if (TmpARowPointer >= 0l)
	    {
	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[this->BlockIndices[j] / IndexStep];
	      int TmpRightMatrixIndex = this->BlockIndices[j] % IndexStep;
	      if (TmpRightMatrix.RowPointers[TmpRightMatrixIndex] >= 0l)
		{
		  Complex Tmp = 0.0;
		  Complex* Tmp2 = LocalTemporaryMatrix[TmpRightMatrixIndex];
		  for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
		    {
		      Tmp += TmpLeftMatrix.MatrixElements[k] * Tmp2[TmpLeftMatrix.ColumnIndices[k]];
		    }
		  vDestination[j] += Tmp;
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

ComplexVector* TensorProductSparseMatrixSelectedBlockHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, 
											      int nbrVectors, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int IndexStep = this->LeftMatrices[0].GetNbrColumn();
  int AMatrixLastIndex = LastComponent / this->LeftMatrices[0].GetNbrRow();
  int BMatrixLastIndex = LastComponent % this->LeftMatrices[0].GetNbrRow();
  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  Complex* Tmp = new Complex[nbrVectors];
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
      SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
      int AMatrixStartingIndex = firstComponent / this->LeftMatrices[0].GetNbrRow();
      int BMatrixStartingIndex = firstComponent % this->LeftMatrices[0].GetNbrRow();
      int TotalIndex = firstComponent;
      for (; AMatrixStartingIndex <  AMatrixLastIndex; ++AMatrixStartingIndex)
	{
	  TmpARowPointer = TmpLeftMatrix.RowPointers[AMatrixStartingIndex];
	  if (TmpARowPointer >= 0l)
	    {
	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[AMatrixStartingIndex];
	      int TmpBMatrixLastIndex = TmpRightMatrix.GetNbrRow();
	      if (AMatrixStartingIndex == (AMatrixLastIndex - 1))
		TmpBMatrixLastIndex = BMatrixLastIndex;
	      for (; BMatrixStartingIndex <  TmpBMatrixLastIndex; ++BMatrixStartingIndex)
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

// test the amount of memory needed for fast multiplication algorithm
//
// return value = amount of memory needed

long TensorProductSparseMatrixSelectedBlockHamiltonian::FastMultiplicationMemory()
{
  GenericHamiltonianPrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);
  long NbrNonZeroMatrixElements = 0l;
  for (int i = 0; i < this->EffectiveHilbertSpaceDimension; ++i)
    {
      NbrNonZeroMatrixElements += this->TemporaryRowPointers[i];
    }
  return (NbrNonZeroMatrixElements << 3);
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that has to be precalcualted
// return value = number of non-zero matrix elements that have to be stored

long TensorProductSparseMatrixSelectedBlockHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int nbrComponent)
{
//   int LastComponent = firstComponent + nbrComponent;
//   int IndexStep = this->LeftMatrices[0].GetNbrColumn();
//   int AMatrixLastIndex = LastComponent / this->LeftMatrices[0].GetNbrRow();
//   int BMatrixLastIndex = LastComponent % this->LeftMatrices[0].GetNbrRow();
//   long TmpARowPointer;
//   long TmpARowLastPointer;
//   long TmpBRowPointer;
//   long TmpBRowLastPointer;
//   long Count = 0l;
//   for (int i = 0; i < this->NbrTensorProducts; ++i)
//     {
//       SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
//       SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
//       for (int j = firstComponent; j < LastComponent; ++j)
// 	{
// 	  TmpARowPointer = TmpLeftMatrix.RowPointers[this->BlockIndices[j] / IndexStep];
// 	  if (TmpARowPointer >= 0l)
// 	    {
// 	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[this->BlockIndices[j] / IndexStep];
// 	      TmpBRowPointer = TmpRightMatrix.RowPointers[this->BlockIndices[j] % IndexStep];
// 	      if (TmpBRowPointer >= 0l)
// 		{
// 		  TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[this->BlockIndices[j] % IndexStep];
// 		  double Tmp= 0.0;
// 		  for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
// 		    {
// 		      int TmpLeftMatrixColumnIndex = TmpLeftMatrix.ColumnIndices[k];
// 		      double Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
// 		      int* TmpInvertBlockIndices = this->InvertBlockIndices + (((long) TmpLeftMatrixColumnIndex) * IndexStep);
// 		      for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
// 			{
// 			  int TmpIndex3 = TmpInvertBlockIndices[TmpRightMatrix.ColumnIndices[l]];
// 			  if (TmpIndex3 >= 0)
// 			    {
// 			      ++Count;
// 			    }
// 			}
// 		    }
// 		}
// 	    }
// 	}
//     }

//   if (this->HamiltonianShift != 0.0)
//     {
//       for (int i = firstComponent; i < LastComponent; ++i)
// 	++Count;
//     }
//   cout << Count << endl;
//   return Count;


  int IndexStep = this->LeftMatrices[0].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent;
  int AMatrixLastIndex = LastComponent / this->LeftMatrices[0].GetNbrRow();
  int BMatrixLastIndex = LastComponent % this->LeftMatrices[0].GetNbrRow();
  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  long NbrNonZeroMatrixElements = 0l;
  int TmpDimension = this->HilbertSpace->GetHilbertSpaceDimension();
  int* TmpNonZeroMatrixIndices = new int [TmpDimension * this->NbrTensorProducts];
  for (int j = firstComponent; j < LastComponent; ++j)
    {
      int TmpNbrNonZeroMatrixElements = 0;
      for (int i = 0; i < this->NbrTensorProducts; ++i)
	{
	  SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
	  SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
	  TmpARowPointer = TmpLeftMatrix.RowPointers[this->BlockIndices[j] / IndexStep];
	  if (TmpARowPointer >= 0l)
	    {
	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[this->BlockIndices[j] / IndexStep];
	      TmpBRowPointer = TmpRightMatrix.RowPointers[this->BlockIndices[j] % IndexStep];
	      if (TmpBRowPointer >= 0l)
		{
		  TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[this->BlockIndices[j] % IndexStep];
		  double Tmp= 0.0;
		  for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
		    {
		      int TmpLeftMatrixColumnIndex = TmpLeftMatrix.ColumnIndices[k];
		      int LocalBlockSize = this->BlockIndexProductTableNbrElements[TmpLeftMatrixColumnIndex];
		      if (LocalBlockSize > 0)
			{
			  double Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
			  if (Tmp2 != 0.0)
			    {
			      int* TmpInvertBlockIndices = this->InvertBlockIndices + (((long) TmpLeftMatrixColumnIndex) * IndexStep);
			      for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
				{
				  int TmpIndex3 = TmpInvertBlockIndices[TmpRightMatrix.ColumnIndices[l]];
				  if ((TmpIndex3 >= 0) && (TmpRightMatrix.MatrixElements[l] != 0.0))
				    {
				      TmpNonZeroMatrixIndices[TmpNbrNonZeroMatrixElements] = TmpIndex3;
				      ++TmpNbrNonZeroMatrixElements;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
      if (this->HamiltonianShift != 0.0)
	{
	  TmpNonZeroMatrixIndices[TmpNbrNonZeroMatrixElements] = j;
	  ++TmpNbrNonZeroMatrixElements;
	}
      if (TmpNbrNonZeroMatrixElements > 0)
	{
	  SortArrayUpOrdering<int>(TmpNonZeroMatrixIndices, TmpNbrNonZeroMatrixElements);
	  long Tmp = 1;
	  int i = 1;
	  while (i < TmpNbrNonZeroMatrixElements)
	    {
	      ++Tmp;
	      while ((i < TmpNbrNonZeroMatrixElements) && (TmpNonZeroMatrixIndices[i] == TmpNonZeroMatrixIndices[i - 1]))
		{
		  ++i;
		}
	      ++i;
	    }
	  this->TemporaryRowPointers[j - this->PrecalculationShift] = Tmp;
	  NbrNonZeroMatrixElements += Tmp;
	}
      else
	{
	  this->TemporaryRowPointers[j - this->PrecalculationShift] = 0;
	}
    }
  delete[] TmpNonZeroMatrixIndices;
  return NbrNonZeroMatrixElements;
}

// enable fast multiplication algorithm
//
  
void TensorProductSparseMatrixSelectedBlockHamiltonian::EnableFastMultiplication()
{
  GenericHamiltonianPrecalculationOperation Operation(this, false);
  Operation.ApplyOperation(this->Architecture);
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that has to be precalcualted

void TensorProductSparseMatrixSelectedBlockHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{
  int IndexStep = this->LeftMatrices[0].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent;
  int AMatrixLastIndex = LastComponent / this->LeftMatrices[0].GetNbrRow();
  int BMatrixLastIndex = LastComponent % this->LeftMatrices[0].GetNbrRow();
  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  long NbrNonZeroMatrixElements = 0l;
  int TmpDimension = this->HilbertSpace->GetHilbertSpaceDimension();
  double* TmpNonZeroMatrixElements = new double [TmpDimension * this->NbrTensorProducts];
  int* TmpNonZeroMatrixIndices = new int [TmpDimension * this->NbrTensorProducts];
  for (int j = firstComponent; j < LastComponent; ++j)
    {
      int TmpNbrNonZeroMatrixElements = 0;
      for (int i = 0; i < this->NbrTensorProducts; ++i)
	{
	  SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
	  SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
	  TmpARowPointer = TmpLeftMatrix.RowPointers[this->BlockIndices[j] / IndexStep];
	  if (TmpARowPointer >= 0l)
	    {
	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[this->BlockIndices[j] / IndexStep];
	      TmpBRowPointer = TmpRightMatrix.RowPointers[this->BlockIndices[j] % IndexStep];
	      if (TmpBRowPointer >= 0l)
		{
		  TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[this->BlockIndices[j] % IndexStep];
		  double Tmp= 0.0;
		  for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
		    {
		      int TmpLeftMatrixColumnIndex = TmpLeftMatrix.ColumnIndices[k];
		      int LocalBlockSize = this->BlockIndexProductTableNbrElements[TmpLeftMatrixColumnIndex];
		      if (LocalBlockSize > 0)
			{
			  double Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
			  if (Tmp2 != 0.0)
			    {
			      long TmpIndex = ((long) TmpLeftMatrixColumnIndex) * IndexStep;
			      long* LocalBlockIndices = this->BlockIndexProductTable[TmpLeftMatrixColumnIndex];
			      int LocalShift = this->BlockIndexProductTableShift[TmpLeftMatrixColumnIndex];
			      for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
				{
				  int TmpIndex2 = TmpIndex + TmpRightMatrix.ColumnIndices[l];
				  int TmpIndex3 = SearchInUnsortedArray<long>(TmpIndex + TmpRightMatrix.ColumnIndices[l], 
								      LocalBlockIndices, LocalBlockSize);
				  if ((TmpIndex3 >= 0) && (TmpRightMatrix.MatrixElements[l] != 0.0))
				    {
				      TmpNonZeroMatrixElements[TmpNbrNonZeroMatrixElements] = Tmp2 * TmpRightMatrix.MatrixElements[l];
				      TmpNonZeroMatrixIndices[TmpNbrNonZeroMatrixElements] = LocalShift + TmpIndex3;
				      ++TmpNbrNonZeroMatrixElements;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
      if (this->HamiltonianShift != 0.0)
	{
	  TmpNonZeroMatrixIndices[TmpNbrNonZeroMatrixElements] = j;
	  ++TmpNbrNonZeroMatrixElements;
	}
      if (TmpNbrNonZeroMatrixElements > 0)
	{
	  SortArrayUpOrdering<double>(TmpNonZeroMatrixIndices, TmpNonZeroMatrixElements, TmpNbrNonZeroMatrixElements);
	  long Shift = this->TemporaryRowPointers[j - this->PrecalculationShift];
	  int i = 0;
	  while (i < TmpNbrNonZeroMatrixElements)
	    {
	      this->TemporaryMatrixElements[Shift] = TmpNonZeroMatrixElements[i];
	      this->TemporaryMatrixColumnIndices[Shift] = TmpNonZeroMatrixIndices[i];
	      ++i;
	      while ((i < TmpNbrNonZeroMatrixElements) && (TmpNonZeroMatrixIndices[i] == TmpNonZeroMatrixIndices[i - 1]))
		{
		  this->TemporaryMatrixElements[Shift] += TmpNonZeroMatrixElements[i];
		  ++i;
		}
	      ++Shift;
	    }
	  NbrNonZeroMatrixElements += Shift - this->TemporaryRowPointers[j - this->PrecalculationShift];
	}
      else
	{
	  this->TemporaryRowPointers[j - this->PrecalculationShift] = 0;
	}
    }
  delete[] TmpNonZeroMatrixElements;
  delete[] TmpNonZeroMatrixIndices;
}
  
// core part of the tensor-multiplication
//
// tensorIndex = index of tensor to consider
// localTemporaryArray = temporary array used to store the partial multiplication
// vSource = vector to be multiplied
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate

void TensorProductSparseMatrixSelectedBlockHamiltonian::LowLevelAddMultiplyTensorCore(int tensorIndex, double** localTemporaryArray, 
										      RealVector& vSource, 
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
      double* Tmp2 = localTemporaryArray[j];
      TmpBRowPointer = TmpRightMatrix.RowPointers[j];
      if (TmpBRowPointer >= 0l)
	{
	  for (int k = 0; k < LeftMatrixDimension; ++k)
	    Tmp2[k] = 0.0;
	  TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[j];
	  for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
	    {
	      double Tmp = TmpRightMatrix.MatrixElements[l] * this->Coefficients[tensorIndex];
	      int TmpIndex = TmpRightMatrix.ColumnIndices[l];
	      for (int k = 0; k < LeftMatrixDimension; ++k)
		{
		  int TmpIndex2 = this->InvertBlockIndices[k * IndexStep + TmpIndex];
		  if (TmpIndex2 >= 0)
		    Tmp2[k] += Tmp * vSource[TmpIndex2];
		}
	    }
	}
    }
}

// core part of the tensor-multiplication
//
// tensorIndex = index of tensor to consider
// localTemporaryArray = temporary array used to store the partial multiplication
// vSource = vector to be multiplied
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate

void TensorProductSparseMatrixSelectedBlockHamiltonian::LowLevelAddMultiplyTensorCore(int tensorIndex, Complex** localTemporaryArray, 
										      ComplexVector& vSource, 
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
      Complex* Tmp2 = localTemporaryArray[j];
      TmpBRowPointer = TmpRightMatrix.RowPointers[j];
      if (TmpBRowPointer >= 0l)
	{
	  for (int k = 0; k < LeftMatrixDimension; ++k)
	    Tmp2[k] = 0.0;
	  TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[j];
	  for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
	    {
	      double Tmp = TmpRightMatrix.MatrixElements[l] * this->Coefficients[tensorIndex];
	      int TmpIndex = TmpRightMatrix.ColumnIndices[l];
	      for (int k = 0; k < LeftMatrixDimension; ++k)
		{
		  int TmpIndex2 = this->InvertBlockIndices[k * IndexStep + TmpIndex];
		  if (TmpIndex2 >= 0)
		    Tmp2[k] += Tmp * vSource[TmpIndex2];
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

void TensorProductSparseMatrixSelectedBlockHamiltonian::LowLevelAddMultiplyTensorCoreDestination(int tensorIndex, double** localTemporaryArray, 
												 RealVector& vDestination, 
												 int firstComponent, int nbrComponent)
{
  int IndexStep = this->RightMatrices[tensorIndex].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent - 1;
  long TmpARowPointer;
  long TmpARowLastPointer;
  SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[tensorIndex];
  SparseRealMatrix& TmpRightMatrix = this->RightMatrices[tensorIndex];
  for (int j = firstComponent; j <= LastComponent; ++j)
    {
      TmpARowPointer = TmpLeftMatrix.RowPointers[this->BlockIndices[j] / IndexStep];
      if (TmpARowPointer >= 0l)
	{
	  TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[this->BlockIndices[j] / IndexStep];
	  int TmpRightMatrixIndex = this->BlockIndices[j] % IndexStep;
	  if (TmpRightMatrix.RowPointers[TmpRightMatrixIndex] >= 0l)
	    {
	      double Tmp = 0.0;
	      double* Tmp2 = localTemporaryArray[TmpRightMatrixIndex];
	      for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
		{
		  Tmp += TmpLeftMatrix.MatrixElements[k] * Tmp2[TmpLeftMatrix.ColumnIndices[k]];
		}
	      vDestination[j] += Tmp;
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

void TensorProductSparseMatrixSelectedBlockHamiltonian::LowLevelAddMultiplyTensorCoreDestination(int tensorIndex, Complex** localTemporaryArray, 
												 ComplexVector& vDestination, 
												 int firstComponent, int nbrComponent)
{
  int IndexStep = this->RightMatrices[tensorIndex].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent - 1;
  long TmpARowPointer;
  long TmpARowLastPointer;
  SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[tensorIndex];
  SparseRealMatrix& TmpRightMatrix = this->RightMatrices[tensorIndex];
  for (int j = firstComponent; j <= LastComponent; ++j)
    {
      TmpARowPointer = TmpLeftMatrix.RowPointers[this->BlockIndices[j] / IndexStep];
      if (TmpARowPointer >= 0l)
	{
	  TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[this->BlockIndices[j] / IndexStep];
	  int TmpRightMatrixIndex = this->BlockIndices[j] % IndexStep;
	  if (TmpRightMatrix.RowPointers[TmpRightMatrixIndex] >= 0l)
	    {
	      Complex Tmp = 0.0;
	      Complex* Tmp2 = localTemporaryArray[TmpRightMatrixIndex];
	      for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
		{
		  Tmp += TmpLeftMatrix.MatrixElements[k] * Tmp2[TmpLeftMatrix.ColumnIndices[k]];
		}
	      vDestination[j] += Tmp;
	    }
	}
    }
}
