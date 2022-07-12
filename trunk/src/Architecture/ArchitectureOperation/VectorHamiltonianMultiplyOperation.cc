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
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"

#include <iostream>
#include <sys/time.h>
#include <stdlib.h>


using std::cout;
using std::endl;


// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// sourceVector = vector to be multiplied by the hamiltonian
// destinationVector = vector where the result has to be stored
// forceNormalMultiplication = force to use direct hamiltonian multiplication (override hamiltonian option)
// forceConjugateMultiplication = force to use hamiltonian for multiplication (override hamiltonian option)
// forceHermitianMultiplication = force to use hermitian property when multiplying  (override hamiltonian option)

VectorHamiltonianMultiplyOperation::VectorHamiltonianMultiplyOperation (AbstractHamiltonian* hamiltonian, Vector* sourceVector, Vector* destinationVector, bool forceNormalMultiplication, bool forceConjugateMultiplication, bool forceHermitianeMultiplication)
{
  this->FirstComponent = 0;
  this->NbrComponent = sourceVector->GetVectorDimension();
  this->Hamiltonian = hamiltonian;
  this->SourceVector = sourceVector;
  this->DestinationVector = destinationVector;
  this->OperationType = AbstractArchitectureOperation::VectorHamiltonianMultiply;
  this->ExecutionTime=0.0;
  if ((forceNormalMultiplication == false) && (forceConjugateMultiplication == false) && (forceHermitianeMultiplication == false))
    {
      this->UseConjugateFlag = false;
      this->UseHermitianFlag = false;
      if (this->Hamiltonian->IsHermitian() == true)
	{
	  this->UseHermitianFlag = true;
	}
      else
	{
	  if (this->Hamiltonian->IsConjugate() == true)
	    {
	      this->UseConjugateFlag = true;
	    }
	}
    }
  else
    {
      if (forceNormalMultiplication == true)
	{
	  this->UseConjugateFlag = false;
	  this->UseHermitianFlag = false;
	}
      else
	{
	  if (forceConjugateMultiplication == true)
	    {
	      this->UseConjugateFlag = true;
	      this->UseHermitianFlag = false;
	    }
	  else
	    {
	      this->UseConjugateFlag = false;
	      this->UseHermitianFlag = true;
	    }
	}
    }
}

// copy constructor 
//
// operation = reference on operation to copy

VectorHamiltonianMultiplyOperation::VectorHamiltonianMultiplyOperation(const VectorHamiltonianMultiplyOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->Hamiltonian = operation.Hamiltonian;
  this->SourceVector = operation.SourceVector;
  this->DestinationVector = operation.DestinationVector;  
  this->OperationType = AbstractArchitectureOperation::VectorHamiltonianMultiply;
  this->UseConjugateFlag = operation.UseConjugateFlag;
  this->UseHermitianFlag = operation.UseHermitianFlag;
  this->ExecutionTime=operation.ExecutionTime;
}
  
// constructor from a master node information
//
// hamiltonian = pointer to the hamiltonian to use
// architecture = pointer to the distributed architecture to use for communications

VectorHamiltonianMultiplyOperation::VectorHamiltonianMultiplyOperation(AbstractHamiltonian* hamiltonian, SimpleMPIArchitecture* architecture)
{
  this->Hamiltonian = hamiltonian;
  this->OperationType = AbstractArchitectureOperation::VectorHamiltonianMultiply;
  long TmpMinimumIndex = 0;
  long TmpMaximumIndex = 0;
  architecture->GetTypicalRange(TmpMinimumIndex, TmpMaximumIndex);
  this->FirstComponent = (int) TmpMinimumIndex;  
  this->NbrComponent = (int) (TmpMaximumIndex - TmpMinimumIndex + 1l);
  int TmpFlag = 0;
  architecture->BroadcastToSlaves(TmpFlag);
  if (TmpFlag == 1)
    this->UseConjugateFlag = true;
  else
    this->UseConjugateFlag = false;
  TmpFlag = 0;
  architecture->BroadcastToSlaves(TmpFlag);
  if (TmpFlag == 1)
    this->UseHermitianFlag = true;
  else
    this->UseHermitianFlag = false;
  if (this->UseConjugateFlag == true)
    {
      this->SourceVector = architecture->BroadcastVector();
      this->DestinationVector = architecture->ScatterVector();
    }
  else
    {
      if (this->UseHermitianFlag == true)
	{
	  this->SourceVector = architecture->BroadcastVector();
	  this->DestinationVector = architecture->BroadcastVectorType();  	  
	}
      else
	{
	  this->SourceVector = architecture->ScatterVector();
	  this->DestinationVector = architecture->BroadcastVectorType();  
	}
    }
}
  
// destructor
//

VectorHamiltonianMultiplyOperation::~VectorHamiltonianMultiplyOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void VectorHamiltonianMultiplyOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void VectorHamiltonianMultiplyOperation::SetDestinationVector (Vector* DestinationVector)
{
  this->DestinationVector = DestinationVector;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* VectorHamiltonianMultiplyOperation::Clone()
{
  return new VectorHamiltonianMultiplyOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool VectorHamiltonianMultiplyOperation::RawApplyOperation()
{
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  gettimeofday (&(TotalStartingTime2), 0);

  if (this->UseConjugateFlag == true)
    {
      this->Hamiltonian->ConjugateMultiply((*(this->SourceVector)), (*(this->DestinationVector)), this->FirstComponent, 
					   this->NbrComponent);
    }
  else if (this->UseHermitianFlag == true)
    {
      this->Hamiltonian->HermitianMultiply((*(this->SourceVector)), (*(this->DestinationVector)), this->FirstComponent, 
					   this->NbrComponent);
    }
  else
    { 
      this->Hamiltonian->Multiply((*(this->SourceVector)), (*(this->DestinationVector)), this->FirstComponent, 
				  this->NbrComponent);
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

bool VectorHamiltonianMultiplyOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  if (this->Hamiltonian->IsHamiltonianVectorOperationCompatible() == false)
    {
      this->RawApplyOperation();
      return true;
    }
  long* SegmentIndices = 0;
  int TmpNbrThreads = architecture->GetNbrThreads();
  bool CleanUp = false;
  long TmpMinimumIndex = 0l;
  long TmpMaximumIndex = 0l;
  architecture->GetTypicalRange(TmpMinimumIndex, TmpMaximumIndex);
  if (TmpMaximumIndex >= 0l)
    {
      this->FirstComponent = (int) TmpMinimumIndex;  
      this->NbrComponent = (int) (TmpMaximumIndex - TmpMinimumIndex + 1l);
    }
  if (Hamiltonian->GetLoadBalancing(TmpNbrThreads, SegmentIndices) == false)
    {
      SegmentIndices = new long[TmpNbrThreads+1];
      CleanUp = true;
      int Step = this->NbrComponent / TmpNbrThreads;
      SegmentIndices[0] = this->FirstComponent;
      for (int i = 0; i < TmpNbrThreads; ++i)
	SegmentIndices[i] = this->FirstComponent + i * Step;
      SegmentIndices[TmpNbrThreads] = this->FirstComponent + this->NbrComponent;
    }
//   else
//     {
//       cout << "Hamiltonian supplied load balancing: ["<<SegmentIndices[0];
//       for (int i=1; i<=TmpNbrThreads; ++i)
// 	cout <<" "<<SegmentIndices[i];
//       cout <<"]"<<endl;
//     }
  VectorHamiltonianMultiplyOperation** TmpOperations = new VectorHamiltonianMultiplyOperation* [architecture->GetNbrThreads()];
  //this->DestinationVector->ClearVector();
  for (int i = 0; i < TmpNbrThreads; ++i)
    {
      TmpOperations[i] = (VectorHamiltonianMultiplyOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpOperations[i]->SetIndicesRange(SegmentIndices[i], SegmentIndices[i+1]-SegmentIndices[i]);
    }
  if (this->UseConjugateFlag == false)
    {
      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
	{
	  TmpOperations[i]->SetDestinationVector(this->DestinationVector->EmptyClone(true));
	}
    }
  architecture->SendJobs();
  if (architecture->VerboseMode() == true)
    {
      char TmpString[512];
      for (int i = 0; i < architecture->GetNbrThreads(); ++i)
	{
	  sprintf (TmpString, "VectorHamiltonianMultiply core operation on SMP id %d done in %.3f seconds", i, TmpOperations[i]->ExecutionTime);
	  architecture->AddToLog(TmpString);
	}
    }
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      if (this->UseConjugateFlag == false)
	{
	  (*(this->DestinationVector)) += (*(TmpOperations[i]->DestinationVector));
	  delete TmpOperations[i]->DestinationVector;
	}
      delete TmpOperations[i];
    }
  delete TmpOperations[0];
  delete[] TmpOperations;
  if (CleanUp)
    delete [] SegmentIndices;
  return true;  
}
  
// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool VectorHamiltonianMultiplyOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __MPI__
   if (architecture->IsMasterNode())
     {
       if (architecture->RequestOperation(this->OperationType) == false)
 	{
 	  return false;
 	}
       int TmpFlag = 0;
       if (this->UseConjugateFlag == true)
	 TmpFlag = 1;
       architecture->BroadcastToSlaves(TmpFlag);
       TmpFlag = 0;
       if (this->UseHermitianFlag == true)
         TmpFlag = 1;
       architecture->BroadcastToSlaves(TmpFlag);
       this->DestinationVector->ClearVector();
       if (this->UseConjugateFlag == true)
	 {
	   architecture->BroadcastVector(this->SourceVector);
	   architecture->ScatterVector(this->DestinationVector);  
	 }
       else
	 {
	   if (this->UseHermitianFlag == true)
	     {
	       architecture->BroadcastVector(this->SourceVector);
	       architecture->BroadcastVectorType(this->DestinationVector);  
	     }
	   else
	     {
	       architecture->ScatterVector(this->SourceVector);
	       architecture->BroadcastVectorType(this->DestinationVector);  
	     }
	 }
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
     {
       this->ArchitectureDependentApplyOperation((SMPArchitecture*) architecture->GetLocalArchitecture());
     }
   else
     this->RawApplyOperation();
   if (architecture->VerboseMode())
     {
       timeval TotalEndingTime;
       gettimeofday (&TotalEndingTime, 0);
       double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
		     (((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
       char TmpString[256];
       sprintf (TmpString, "VectorHamiltonianMultiply core operation done in %.3f seconds", Dt);
       architecture->AddToLog(TmpString);
     }
   if ((architecture->IsMasterNode()) && (architecture->VerboseMode()))
     gettimeofday (&TotalStartingTime, 0);
  if (this->UseConjugateFlag == true)
    {
      architecture->ReassembleVector(*(this->DestinationVector));
    }
  else
    {
      architecture->SumVector(*(this->DestinationVector));
    }
   if ((architecture->IsMasterNode()) && (architecture->VerboseMode()))
     {
       timeval TotalEndingTime;
       gettimeofday (&TotalEndingTime, 0);
       double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
		     (((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
       char TmpString[256];
       sprintf (TmpString, "VectorHamiltonianMultiply sum operation done in %.3f seconds", Dt);
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

  
