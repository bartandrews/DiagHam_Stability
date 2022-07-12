////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of hamiltonian vector multiplication operation           //
//                                                                            //
//                        last modification : 23/10/2002                      //
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


#include "config.h"
#include "Architecture/ArchitectureOperation/VectorSparseTensorMultiplyOperation.h"
#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "Hamiltonian/TensorProductSparseMatrixSelectedBlockHamiltonian.h"

#include <iostream>
#include <sys/time.h>
#include <stdlib.h>


using std::cout;
using std::endl;


// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// tensorIndex = index of tensor to consider
// destinationVector = vector where the result has to be stored
// fullHilbertSpace = split the workload with respect to the full Hilbert space dimension instead of the auxillary space dimension

VectorSparseTensorMultiplyOperation::VectorSparseTensorMultiplyOperation (TensorProductSparseMatrixHamiltonian* hamiltonian,
									  int tensorIndex, Vector* destinationVector, bool fullHilbertSpace)
{
  this->TensorIndex = tensorIndex;
  this->FirstComponent = 0;
  this->TensorHamiltonian = hamiltonian;
  if (fullHilbertSpace == false)
    this->NbrComponent = this->TensorHamiltonian->RightMatrices[this->TensorIndex].GetNbrRow();
  else
    this->NbrComponent = ((TensorProductSparseMatrixSelectedBlockHamiltonian*) this->TensorHamiltonian)->BlockSize;
  this->DestinationVector = destinationVector;
  this->OperationType = AbstractArchitectureOperation::VectorSparseTensorMultiply;
  this->ExecutionTime=0.0;
}

// copy constructor 
//
// operation = reference on operation to copy

VectorSparseTensorMultiplyOperation::VectorSparseTensorMultiplyOperation(const VectorSparseTensorMultiplyOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->TensorIndex = operation.TensorIndex;
  this->TensorHamiltonian = operation.TensorHamiltonian;
  this->DestinationVector = operation.DestinationVector;  
  this->OperationType = AbstractArchitectureOperation::VectorSparseTensorMultiply;
  this->ExecutionTime=operation.ExecutionTime;
}
  
// destructor
//

VectorSparseTensorMultiplyOperation::~VectorSparseTensorMultiplyOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void VectorSparseTensorMultiplyOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void VectorSparseTensorMultiplyOperation::SetDestinationVector (Vector* DestinationVector)
{
  this->DestinationVector = DestinationVector;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* VectorSparseTensorMultiplyOperation::Clone()
{
  return new VectorSparseTensorMultiplyOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool VectorSparseTensorMultiplyOperation::RawApplyOperation()
{
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  gettimeofday (&(TotalStartingTime2), 0);
  if ((this->DestinationVector->GetVectorType() & Vector::RealDatas) != 0)
    {
      this->TensorHamiltonian->LowLevelAddMultiplyTensorCoreDestination(this->TensorIndex, this->TensorHamiltonian->TemporaryArray, 
									*((RealVector*) this->DestinationVector), 
									this->FirstComponent, this->NbrComponent);
    }
  else
    {
      this->TensorHamiltonian->LowLevelAddMultiplyTensorCoreDestination(this->TensorIndex, this->TensorHamiltonian->ComplexTemporaryArray, 
									*((ComplexVector*) this->DestinationVector), 
									this->FirstComponent, this->NbrComponent);
    }

  gettimeofday (&(TotalEndingTime2), 0);
  this->ExecutionTime = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool VectorSparseTensorMultiplyOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  long *SegmentIndices=0;
  int TmpNbrThreads = architecture->GetNbrThreads();
  bool CleanUp = false;
  long TmpMinimumIndex = 0;
  long TmpMaximumIndex = 0;
  SegmentIndices = new long[TmpNbrThreads + 1];
  CleanUp = true;
  int Step = this->NbrComponent / TmpNbrThreads;
  SegmentIndices[0] = this->FirstComponent;
  for (int i = 0; i < TmpNbrThreads; ++i)
    SegmentIndices[i] = this->FirstComponent + i * Step;
  SegmentIndices[TmpNbrThreads] = this->FirstComponent + this->NbrComponent;

  VectorSparseTensorMultiplyOperation** TmpOperations = new VectorSparseTensorMultiplyOperation* [architecture->GetNbrThreads()];
  for (int i = 0; i < TmpNbrThreads; ++i)
    {
      TmpOperations[i] = (VectorSparseTensorMultiplyOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpOperations[i]->SetIndicesRange(SegmentIndices[i], SegmentIndices[i + 1] - SegmentIndices[i]);
    }

  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      TmpOperations[i]->SetDestinationVector(this->DestinationVector->EmptyClone(true));
    }

  architecture->SendJobs();
  if (architecture->VerboseMode() == true)
    {
      char TmpString[512];
      for (int i = 0; i < architecture->GetNbrThreads(); ++i)
	{
	  sprintf (TmpString, "VectorSparseTensorMultiply core operation on SMP id %d done in %.3f seconds", i, TmpOperations[i]->ExecutionTime);
	  architecture->AddToLog(TmpString);
	}
    }
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      (*(this->DestinationVector)) += (*(TmpOperations[i]->DestinationVector));
      delete TmpOperations[i]->DestinationVector;
    }
  delete TmpOperations[0];
  delete[] TmpOperations;
  if (CleanUp)
    delete [] SegmentIndices;
  return true;  
}
