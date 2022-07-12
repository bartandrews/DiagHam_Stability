////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of hamiltonian-sparse tensor multiplication operation      //
//                                                                            //
//                        last modification : 21/07/2013                      //
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


#ifndef MPOAPPLYONTENSOROPERATION_H
#define MPOAPPLYONTENSOROPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"

#include "Tensor/Tensor3.h"
#include "MPSObjects/AbstractMPOperatorOBC.h"
#include <sys/time.h>

class AbstractSparseTensor;
class Vector;


template <typename T>
class MPOApplyOnTensorOperation : public AbstractArchitectureOperation
{

 protected:

  Tensor3<T> * SourceTensor;
  Tensor3<T> * DestinationTensor;  
  AbstractMPOperatorOBC * MPOperator;
  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;
  bool RightFlag;
  // execution time measured in RawApply
  double ExecutionTime;
  
 public:
  
  // constructor 
  //
  // hamiltonian = pointer to the hamiltonian to use
  // tensorIndex = index of tensor to consider
  // destinationVector = vector where the result has to be stored
  // fullHilbertSpace = split the workload with respect to the full Hilbert space dimension instead of the auxillary space dimension
  MPOApplyOnTensorOperation(Tensor3<T> * sourceTensor, Tensor3<T> * destinationTensor,  AbstractMPOperatorOBC * mPOperator, bool rightFlag);

  // copy constructor 
  //
  // operation = reference on operation to copy
  MPOApplyOnTensorOperation(const MPOApplyOnTensorOperation & operation);

  // destructor
  //
  ~MPOApplyOnTensorOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();

 protected:

  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
};


// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// tensorIndex = index of tensor to consider
// destinationVector = vector where the result has to be stored
// fullHilbertSpace = split the workload with respect to the full Hilbert space dimension instead of the auxillary space dimension

template <typename T>
MPOApplyOnTensorOperation<T>::MPOApplyOnTensorOperation( Tensor3<T> * sourceTensor, Tensor3<T> * destinationTensor,  AbstractMPOperatorOBC * mPOperator, bool rightFlag)
{
  this->SourceTensor = sourceTensor;
  this->DestinationTensor =  destinationTensor;
  this->MPOperator = mPOperator;
  this->RightFlag = rightFlag;
//  this->OperationType = AbstractArchitectureOperation::TensorMatrixContractionOperation;
  this->ExecutionTime=0.0;
}

// copy constructor 
//
// operation = reference on operation to copy

template <typename T>
MPOApplyOnTensorOperation<T>::MPOApplyOnTensorOperation(const MPOApplyOnTensorOperation & operation)
{
  this->SourceTensor = operation.SourceTensor;
  this->DestinationTensor = operation.DestinationTensor;
  this->MPOperator = operation.MPOperator;
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->RightFlag = operation.RightFlag;
//  this->OperationType = AbstractArchitectureOperation::TensorMatrixContractionOperation;
   this->ExecutionTime = operation.ExecutionTime;
}
  
// destructor
//
template <typename T>
MPOApplyOnTensorOperation<T>::~MPOApplyOnTensorOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component
template <typename T>
void MPOApplyOnTensorOperation<T>::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent; 
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation
template <typename T>
AbstractArchitectureOperation* MPOApplyOnTensorOperation<T>::Clone()
{
  return new MPOApplyOnTensorOperation<T> (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs
template <typename T>
bool MPOApplyOnTensorOperation<T>::RawApplyOperation()
{
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  gettimeofday (&(TotalStartingTime2), 0);

if(this->RightFlag)
  this->MPOperator->MPOApplyOnTensorOnTheRightCore(this->DestinationTensor,this->SourceTensor, this->FirstComponent, this->NbrComponent);
else
  this->MPOperator->MPOApplyOnTensorOnTheLeftCore(this->DestinationTensor,this->SourceTensor, this->FirstComponent, this->NbrComponent);

  gettimeofday (&(TotalEndingTime2), 0);
  this->ExecutionTime = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs
template <typename T>
bool MPOApplyOnTensorOperation<T>::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int TmpNbrThreads = architecture->GetNbrThreads();
  long TmpMinimumIndex = 0;
  long TmpMaximumIndex = 0;
  long * SegmentIndices = new long[TmpNbrThreads + 1];
  int Step = this->NbrComponent / TmpNbrThreads;
  SegmentIndices[0] = this->FirstComponent;
  for (int i = 0; i < TmpNbrThreads; ++i)
    SegmentIndices[i] = this->FirstComponent + i * Step;
  SegmentIndices[TmpNbrThreads] = this->FirstComponent + this->NbrComponent;

  MPOApplyOnTensorOperation** TmpOperations = new MPOApplyOnTensorOperation * [architecture->GetNbrThreads()];
  for (int i = 0; i < TmpNbrThreads; ++i)
    {
      TmpOperations[i] = (MPOApplyOnTensorOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpOperations[i]->SetIndicesRange(SegmentIndices[i], SegmentIndices[i + 1] - SegmentIndices[i]);
    }

  architecture->SendJobs();

/*  if (architecture->VerboseMode() == true)
    {
      char TmpString[512];
      for (int i = 0; i < architecture->GetNbrThreads(); ++i)
	{
	  sprintf (TmpString, "VectorSparseTensorMultiply core operation on SMP id %d done in %.3f seconds", i, TmpOperations[i]->ExecutionTime);
	  architecture->AddToLog(TmpString);
	}
    }
 */
  delete TmpOperations[0];
  delete[] TmpOperations;
  delete[] SegmentIndices;
  return true;  
}


#endif
