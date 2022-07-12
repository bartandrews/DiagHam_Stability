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


#ifndef TENSORVECTORCONTRACTIONPERATION_H
#define TENSORVECTORCONTRACTIONPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"

#include "Vector/Vector.h"
#include "Tensor/Tensor3.h"
#include "MPSObjects/AbstractMPOperatorOBC.h"
#include <sys/time.h>


template <typename T>
class TensorVectorContractionOperation : public AbstractArchitectureOperation
{ };

template <>
class TensorVectorContractionOperation<double>: public AbstractArchitectureOperation
{ 
 protected:

  Tensor3<double> * SourceTensor;
  Tensor3<double> * DestinationTensor;  
  AbstractMPOperatorOBC * MPOperator;
  RealVector * InputVector;
  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;
  bool RightFlag;
  bool TwoSitesFlag;
  // execution time measured in RawApply
  double ExecutionTime;
  
 public:
  
  // constructor 
  //
  // hamiltonian = pointer to the hamiltonian to use
  // tensorIndex = index of tensor to consider
  // destinationVector = vector where the result has to be stored
  // fullHilbertSpace = split the workload with respect to the full Hilbert space dimension instead of the auxillary space dimension
  inline TensorVectorContractionOperation(Tensor3<double> * sourceTensor, Tensor3<double> * destinationTensor,  RealVector * vector,  AbstractMPOperatorOBC * mPOperator,bool rightFlag, bool twoSites = false);

  // copy constructor 
  //
  // operation = reference on operation to copy
  inline TensorVectorContractionOperation(const TensorVectorContractionOperation & operation);

  // destructor
  //
  inline ~TensorVectorContractionOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  inline void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // clone operation
  //
  // return value = pointer to cloned operation
  inline AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  inline bool RawApplyOperation();

 protected:

  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  inline bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
};


// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// tensorIndex = index of tensor to consider
// destinationVector = vector where the result has to be stored
// fullHilbertSpace = split the workload with respect to the full Hilbert space dimension instead of the auxillary space dimension


TensorVectorContractionOperation<double>::TensorVectorContractionOperation( Tensor3<double> * sourceTensor, Tensor3<double> * destinationTensor, RealVector * vector,  AbstractMPOperatorOBC * mPOperator , bool rightFlag, bool twoSites)
{
  this->SourceTensor = sourceTensor;
  this->DestinationTensor =  destinationTensor;
  this->InputVector  = vector;
  this->MPOperator = mPOperator;
  this->RightFlag = rightFlag;
  this->TwoSitesFlag = twoSites;
//  this->OperationType = AbstractArchitectureOperation::TensorMatrixContractionOperation;
  this->ExecutionTime=0.0;
}

// copy constructor 
//
// operation = reference on operation to copy

TensorVectorContractionOperation<double>::TensorVectorContractionOperation(const TensorVectorContractionOperation & operation)
{
  this->SourceTensor = operation.SourceTensor;
  this->DestinationTensor = operation.DestinationTensor;
  this->InputVector = operation.InputVector;
  this->MPOperator = operation.MPOperator;
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->RightFlag = operation.RightFlag;
  this->TwoSitesFlag = operation.TwoSitesFlag;
//  this->OperationType = AbstractArchitectureOperation::TensorMatrixContractionOperation;
   this->ExecutionTime = operation.ExecutionTime;
}
  
// destructor
//
TensorVectorContractionOperation<double>::~TensorVectorContractionOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component
void TensorVectorContractionOperation<double>::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent; 
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation
AbstractArchitectureOperation* TensorVectorContractionOperation<double>::Clone()
{
  return new TensorVectorContractionOperation<double> (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs
bool TensorVectorContractionOperation<double>::RawApplyOperation()
{
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  gettimeofday (&(TotalStartingTime2), 0);
if (this->TwoSitesFlag == true)
{
if(this->RightFlag)
  this->MPOperator->LowLevelMultiplyCoreTwoSitesFirst(this->DestinationTensor,this->SourceTensor, *this->InputVector, this->FirstComponent, this->NbrComponent);
else
  this->MPOperator->LowLevelMultiplyCoreTwoSitesSecond(this->SourceTensor,this->DestinationTensor,*this->InputVector,  this->FirstComponent, this->NbrComponent);
}
else
{
if(this->RightFlag)
  this->MPOperator->LowLevelMultiplyCoreFirst(this->DestinationTensor,this->SourceTensor, *this->InputVector, this->FirstComponent, this->NbrComponent);
else
  this->MPOperator->LowLevelMultiplyCoreSecond(this->SourceTensor,this->DestinationTensor,*this->InputVector,  this->FirstComponent, this->NbrComponent);
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
bool TensorVectorContractionOperation<double>::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
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

  TensorVectorContractionOperation ** TmpOperations = new TensorVectorContractionOperation * [architecture->GetNbrThreads()];
  for (int i = 0; i < TmpNbrThreads; ++i)
    {
      TmpOperations[i] = (TensorVectorContractionOperation*) this->Clone();
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


template < >
class TensorVectorContractionOperation <Complex>: public AbstractArchitectureOperation
{ 
 protected:

  Tensor3<Complex> * SourceTensor;
  Tensor3<Complex> * DestinationTensor;  
  AbstractMPOperatorOBC * MPOperator;
  ComplexVector * InputVector;
  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;
  bool RightFlag;
  bool TwoSitesFlag;
  // execution time measured in RawApply
  double ExecutionTime;
  
 public:
  
  // constructor 
  //
  // hamiltonian = pointer to the hamiltonian to use
  // tensorIndex = index of tensor to consider
  // destinationVector = vector where the result has to be stored
  // fullHilbertSpace = split the workload with respect to the full Hilbert space dimension instead of the auxillary space dimension
  inline  TensorVectorContractionOperation(Tensor3<Complex> * sourceTensor, Tensor3<Complex> * destinationTensor, ComplexVector * vector,  AbstractMPOperatorOBC * mPOperator, bool rightFlag, bool twoSites = false);

  // copy constructor 
  //
  // operation = reference on operation to copy
  inline TensorVectorContractionOperation(const TensorVectorContractionOperation & operation);

  // destructor
  //
  inline ~TensorVectorContractionOperation();
   
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  inline void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // clone operation
  //
  // return value = pointer to cloned operation
  inline AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  inline bool RawApplyOperation();

 protected:

  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  inline bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
};


// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// tensorIndex = index of tensor to consider
// destinationVector = vector where the result has to be stored
// fullHilbertSpace = split the workload with respect to the full Hilbert space dimension instead of the auxillary space dimension

TensorVectorContractionOperation<Complex>::TensorVectorContractionOperation( Tensor3<Complex> * sourceTensor, Tensor3<Complex> * destinationTensor, ComplexVector * vector,  AbstractMPOperatorOBC * mPOperator , bool rightFlag, bool twoSites)
{
  this->SourceTensor = sourceTensor;
  this->DestinationTensor =  destinationTensor;
  this->InputVector  = vector;
  this->MPOperator = mPOperator;
  this->RightFlag = rightFlag;
  this->TwoSitesFlag = twoSites;
//  this->OperationType = AbstractArchitectureOperation::TensorMatrixContractionOperation;
  this->ExecutionTime=0.0;
}

// copy constructor 
//
// operation = reference on operation to copy

TensorVectorContractionOperation<Complex>::TensorVectorContractionOperation(const TensorVectorContractionOperation<Complex> & operation)
{
  this->SourceTensor = operation.SourceTensor;
  this->DestinationTensor = operation.DestinationTensor;
  this->InputVector = operation.InputVector;
  this->MPOperator = operation.MPOperator;
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->RightFlag = operation.RightFlag;
  this->TwoSitesFlag = operation.TwoSitesFlag;
//  this->OperationType = AbstractArchitectureOperation::TensorMatrixContractionOperation;
   this->ExecutionTime = operation.ExecutionTime;
}
  
// destructor
//
TensorVectorContractionOperation<Complex>::~TensorVectorContractionOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void TensorVectorContractionOperation<Complex>::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent; 
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* TensorVectorContractionOperation<Complex>::Clone()
{
  return new TensorVectorContractionOperation<Complex> (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool TensorVectorContractionOperation<Complex>::RawApplyOperation()
{
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  gettimeofday (&(TotalStartingTime2), 0);

if (this->TwoSitesFlag == true)
{
if(this->RightFlag)
  this->MPOperator->LowLevelMultiplyCoreTwoSitesFirst(this->DestinationTensor,this->SourceTensor, *this->InputVector, this->FirstComponent, this->NbrComponent);
else
  this->MPOperator->LowLevelMultiplyCoreTwoSitesSecond(this->SourceTensor,this->DestinationTensor,*this->InputVector,  this->FirstComponent, this->NbrComponent);
}
else
{
if(this->RightFlag)
  this->MPOperator->LowLevelMultiplyCoreFirst(this->DestinationTensor,this->SourceTensor, *this->InputVector, this->FirstComponent, this->NbrComponent);
else
  this->MPOperator->LowLevelMultiplyCoreSecond(this->SourceTensor,this->DestinationTensor,*this->InputVector,  this->FirstComponent, this->NbrComponent);
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

bool TensorVectorContractionOperation<Complex>::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
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

  TensorVectorContractionOperation ** TmpOperations = new TensorVectorContractionOperation * [architecture->GetNbrThreads()];
  for (int i = 0; i < TmpNbrThreads; ++i)
    {
      TmpOperations[i] = (TensorVectorContractionOperation*) this->Clone();
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
