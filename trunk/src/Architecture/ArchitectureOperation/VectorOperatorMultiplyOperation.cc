////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of operator vector multiplication operation             //
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
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"
#include "Operator/AbstractOperator.h"
#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"

#include <iostream>
#include <sys/time.h>
#include <stdlib.h>


using std::cout;
using std::endl;


// constructor 
//
// operator = pointer to the operator to use
// sourceVector = vector to be multiplied by the operator
// destinationVector = vector where the result has to be stored

VectorOperatorMultiplyOperation::VectorOperatorMultiplyOperation (AbstractOperator* theOperator, Vector* sourceVector, Vector* destinationVector)
{
  this->FirstComponent = 0;
  this->NbrComponent = sourceVector->GetVectorDimension();
  this->Operator = theOperator;
  this->SourceVector = sourceVector;
  this->DestinationVector = destinationVector;  
  this->OperationType = AbstractArchitectureOperation::VectorOperatorMultiply;
}

// copy constructor 
//
// operation = reference on operation to copy

VectorOperatorMultiplyOperation::VectorOperatorMultiplyOperation(const VectorOperatorMultiplyOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->Operator = operation.Operator->Clone ();
  this->SourceVector = operation.SourceVector;
  this->DestinationVector = operation.DestinationVector;  
  this->OperationType = AbstractArchitectureOperation::VectorOperatorMultiply;
}
  
// constructor from a master node information
//
// operator = pointer to the operator to use
// architecture = pointer to the distributed architecture to use for communications

VectorOperatorMultiplyOperation::VectorOperatorMultiplyOperation(AbstractOperator* theOperator, SimpleMPIArchitecture* architecture)
{
  this->Operator = theOperator;
  this->OperationType = AbstractArchitectureOperation::VectorOperatorMultiply;
  long TmpMinimumIndex = 0;
  long TmpMaximumIndex = 0;
  architecture->GetTypicalRange(TmpMinimumIndex, TmpMaximumIndex);
  this->FirstComponent = (int) TmpMinimumIndex;  
  this->NbrComponent = (int) (TmpMaximumIndex - TmpMinimumIndex + 1l);
  this->SourceVector = architecture->ScatterVector();
  this->DestinationVector = architecture->BroadcastVectorType();  
}
  
// destructor
//

VectorOperatorMultiplyOperation::~VectorOperatorMultiplyOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void VectorOperatorMultiplyOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void VectorOperatorMultiplyOperation::SetDestinationVector (Vector* DestinationVector)
{
  this->DestinationVector = DestinationVector;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* VectorOperatorMultiplyOperation::Clone()
{
  return new VectorOperatorMultiplyOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool VectorOperatorMultiplyOperation::RawApplyOperation()
{
  this->Operator->Multiply((*(this->SourceVector)), (*(this->DestinationVector)), this->FirstComponent, 
			      this->NbrComponent);
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool VectorOperatorMultiplyOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrComponent / architecture->GetNbrThreads();
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  VectorOperatorMultiplyOperation** TmpOperations = new VectorOperatorMultiplyOperation* [architecture->GetNbrThreads()];
  this->DestinationVector->ClearVector();
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      TmpOperations[i] = (VectorOperatorMultiplyOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      TmpFirstComponent += Step;
    }
  TmpOperations[ReducedNbrThreads] = (VectorOperatorMultiplyOperation*) this->Clone();
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, 
						    this->FirstComponent + this->NbrComponent - TmpFirstComponent);  
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      TmpOperations[i]->SetDestinationVector(this->DestinationVector->EmptyClone(true));
    }
  architecture->SendJobs();
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      (*(this->DestinationVector)) += (*(TmpOperations[i]->DestinationVector));
      delete TmpOperations[i]->DestinationVector;
      delete TmpOperations[i];
    }
  delete TmpOperations[0];
  delete[] TmpOperations;
  return true;  
}
  
// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool VectorOperatorMultiplyOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __MPI__
   if (architecture->IsMasterNode())
     {
       if (architecture->RequestOperation(this->OperationType) == false)
 	{
 	  return false;
 	}
       architecture->ScatterVector(this->SourceVector);
       architecture->BroadcastVectorType(this->DestinationVector);  
     }
   long TmpMinimumIndex = 0;
   long TmpMaximumIndex = 0;
   architecture->GetTypicalRange(TmpMinimumIndex, TmpMaximumIndex);
   this->FirstComponent = (int) TmpMinimumIndex;  
   this->NbrComponent = (int) (TmpMaximumIndex - TmpMinimumIndex + 1l);
   timeval TotalStartingTime;
   if (architecture->VerboseMode())
     gettimeofday (&TotalStartingTime, 0);
   if (architecture->GetLocalArchitecture()->GetArchitectureID() == AbstractArchitecture::SMP)
     this->ArchitectureDependentApplyOperation((SMPArchitecture*) architecture->GetLocalArchitecture());
   else
     this->RawApplyOperation();
   if (architecture->VerboseMode())
     {
       timeval TotalEndingTime;
       gettimeofday (&TotalEndingTime, 0);
       double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
		     (((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
       char TmpString[256];
       sprintf (TmpString, "VectorOperatorMultiply core operation done in %.3f seconds", Dt);
       architecture->AddToLog(TmpString);
     }
   if ((architecture->IsMasterNode()) && (architecture->VerboseMode()))
     gettimeofday (&TotalStartingTime, 0);
   architecture->SumVector(*(this->DestinationVector));
   if ((architecture->IsMasterNode()) && (architecture->VerboseMode()))
     {
       timeval TotalEndingTime;
       gettimeofday (&TotalEndingTime, 0);
       double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
		     (((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
       char TmpString[256];
       sprintf (TmpString, "VectorOperatorMultiply sum operation done in %.3f seconds", Dt);
       architecture->AddToLog(TmpString, true);
     }
   if (architecture->IsMasterNode() == false)
     {
       delete this->DestinationVector;
       delete this->SourceVector;
     }

  return true;
#else
  return this->RawApplyOperation();
#endif
}

  
