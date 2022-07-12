////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of multiple hamiltonian-vector multiplication operation        //
//                                                                            //
//                        last modification : 15/03/2005                      //
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
#include "Architecture/ArchitectureOperation/MultipleVectorHamiltonianMultiplyOperation.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Architecture/SMPArchitecture.h"


#include <sys/time.h>
#include <stdlib.h>


// constructor for real vectors
//
// hamiltonian = pointer to the hamiltonian to use
// sourceVectors = array of vectors to be multiplied by the hamiltonian
// destinationVectors = array of vectors where the result has to be stored
// nbrVectors = number of vectors that have to be evaluated together

MultipleVectorHamiltonianMultiplyOperation::MultipleVectorHamiltonianMultiplyOperation (AbstractHamiltonian* hamiltonian, RealVector* sourceVectors, 
											RealVector* destinationVectors, int nbrVectors)
{
  this->FirstComponent = 0;
  this->NbrComponent = sourceVectors[0].GetVectorDimension();
  this->Hamiltonian = hamiltonian;
  this->RealSourceVectors = sourceVectors;
  this->RealDestinationVectors = destinationVectors;
  this->ComplexSourceVectors = 0;
  this->ComplexDestinationVectors = 0; 
  this->NbrVectors = nbrVectors;
  this->OperationType = AbstractArchitectureOperation::MultipleVectorHamiltonianMultiply;
}

// constructor for complex vectors
//
// hamiltonian = pointer to the hamiltonian to use
// sourceVectors = array of vectors to be multiplied by the hamiltonian
// destinationVectors = array of vectors where the result has to be stored
// nbrVectors = number of vectors that have to be evaluated together

MultipleVectorHamiltonianMultiplyOperation::MultipleVectorHamiltonianMultiplyOperation (AbstractHamiltonian* hamiltonian, ComplexVector* sourceVectors, 
											ComplexVector* destinationVectors, int nbrVectors)
{
  this->FirstComponent = 0;
  this->NbrComponent = sourceVectors[0].GetVectorDimension();
  this->Hamiltonian = hamiltonian;
  this->RealSourceVectors = 0;
  this->RealDestinationVectors = 0;
  this->ComplexSourceVectors = sourceVectors;
  this->ComplexDestinationVectors = destinationVectors; 
  this->NbrVectors = nbrVectors;
  this->OperationType = AbstractArchitectureOperation::MultipleVectorHamiltonianMultiply;
}

// copy constructor 
//
// operation = reference on operation to copy

MultipleVectorHamiltonianMultiplyOperation::MultipleVectorHamiltonianMultiplyOperation(const MultipleVectorHamiltonianMultiplyOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->Hamiltonian = operation.Hamiltonian;
  this->RealSourceVectors = operation.RealSourceVectors;
  this->RealDestinationVectors = operation.RealDestinationVectors;  
  this->ComplexSourceVectors = operation.ComplexSourceVectors;
  this->ComplexDestinationVectors = operation.ComplexDestinationVectors;  
  this->NbrVectors = operation.NbrVectors;
  this->OperationType = AbstractArchitectureOperation::MultipleVectorHamiltonianMultiply;
}
  
// constructor from a master node information
//
// hamiltonian = pointer to the hamiltonian to use
// architecture = pointer to the distributed architecture to use for communications

MultipleVectorHamiltonianMultiplyOperation::MultipleVectorHamiltonianMultiplyOperation(AbstractHamiltonian* hamiltonian, SimpleMPIArchitecture* architecture)
{
  this->OperationType = AbstractArchitectureOperation::MultipleVectorHamiltonianMultiply;
  this->Hamiltonian = hamiltonian;
  long TmpMinimumIndex = 0;
  long TmpMaximumIndex = 0;
  architecture->GetTypicalRange(TmpMinimumIndex, TmpMaximumIndex);
  this->FirstComponent = (int) TmpMinimumIndex;  
  this->NbrComponent = (int) (TmpMaximumIndex - TmpMinimumIndex + 1l);
  Vector** TmpSourceVectors = architecture->BroadcastVectorArray(this->NbrVectors);
  Vector** TmpDestinationVectors = architecture->BroadcastVectorTypeArray(this->NbrVectors);
  if (TmpSourceVectors[0]->GetVectorType() == Vector::RealDatas)
    {
      this->RealSourceVectors = new RealVector[this->NbrVectors];
      this->RealDestinationVectors = new RealVector[this->NbrVectors];
      for (int i = 0; i < this->NbrVectors; ++i)
	{
	  this->RealSourceVectors[i] = *((RealVector*) TmpSourceVectors[i]);
	  this->RealDestinationVectors[i] = *((RealVector*) TmpDestinationVectors[i]);
	  delete TmpSourceVectors[i];
	  delete TmpDestinationVectors[i];
	}
      delete[] TmpSourceVectors;
      delete[] TmpDestinationVectors;
      this->ComplexSourceVectors = 0;
      this->ComplexDestinationVectors = 0;       
    }
  else
    {
      this->RealSourceVectors = 0;
      this->RealDestinationVectors = 0;
      this->ComplexSourceVectors = (ComplexVector*) TmpSourceVectors;
      this->ComplexDestinationVectors = (ComplexVector*) TmpDestinationVectors;       
    }
}
  
// destructor
//

MultipleVectorHamiltonianMultiplyOperation::~MultipleVectorHamiltonianMultiplyOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void MultipleVectorHamiltonianMultiplyOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set destination vectors 
// 
// destinationVectors = array of vector where the result has to be stored

void MultipleVectorHamiltonianMultiplyOperation::SetDestinationVectors (RealVector* destinationVectors)
{
  this->RealDestinationVectors = (RealVector*) destinationVectors;
  this->ComplexDestinationVectors = 0;
}

// set destination vectors 
// 
// destinationVectors = array of vector where the result has to be stored

void MultipleVectorHamiltonianMultiplyOperation::SetDestinationVectors (ComplexVector* destinationVectors)
{
  this->ComplexDestinationVectors = (ComplexVector*) destinationVectors;
  this->RealDestinationVectors = 0;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* MultipleVectorHamiltonianMultiplyOperation::Clone()
{
  return new MultipleVectorHamiltonianMultiplyOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool MultipleVectorHamiltonianMultiplyOperation::RawApplyOperation()
{
  if (this->RealSourceVectors != 0)
    {
      this->Hamiltonian->LowLevelMultipleMultiply(this->RealSourceVectors, this->RealDestinationVectors, this->NbrVectors, this->FirstComponent, 
						  this->NbrComponent);
    }
  else
    {
      this->Hamiltonian->LowLevelMultipleMultiply(this->ComplexSourceVectors, this->ComplexDestinationVectors, this->NbrVectors, this->FirstComponent, 
						  this->NbrComponent);
    }
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool MultipleVectorHamiltonianMultiplyOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  bool RealFlag = false;
  if (this->ComplexDestinationVectors == 0)
    RealFlag = true;
  int Step = this->NbrComponent / architecture->GetNbrThreads();
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  MultipleVectorHamiltonianMultiplyOperation** TmpOperations = new MultipleVectorHamiltonianMultiplyOperation* [architecture->GetNbrThreads()];
  if (RealFlag == true)
    {
      for (int i = 0; i < this->NbrVectors; ++i)
	this->RealDestinationVectors[i].ClearVector();
    }
  else
    {
      for (int i = 0; i < this->NbrVectors; ++i)
	this->ComplexDestinationVectors[i].ClearVector();
    }
   for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      TmpOperations[i] = (MultipleVectorHamiltonianMultiplyOperation*) this->Clone();
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpFirstComponent += Step;
    }
  TmpOperations[ReducedNbrThreads] = (MultipleVectorHamiltonianMultiplyOperation*) this->Clone();
  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->NbrComponent + this->FirstComponent - TmpFirstComponent);  
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      if (RealFlag == true)
	TmpOperations[i]->SetDestinationVectors((RealVector*) this->RealDestinationVectors[0].EmptyCloneArray(this->NbrVectors, true));
      else
	TmpOperations[i]->SetDestinationVectors((ComplexVector*) this->ComplexDestinationVectors[0].EmptyCloneArray(this->NbrVectors, true));

    }
  architecture->SendJobs();
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      if (RealFlag == true)
	{
	  for (int j = 0; j < this->NbrVectors; ++j)
	    this->RealDestinationVectors[j] += TmpOperations[i]->RealDestinationVectors[j];
	  delete[] TmpOperations[i]->RealDestinationVectors;
	}
      else
	{
	  for (int j = 0; j < this->NbrVectors; ++j)
	    this->ComplexDestinationVectors[j] += TmpOperations[i]->ComplexDestinationVectors[j];
	  delete[] TmpOperations[i]->ComplexDestinationVectors;
	}
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

bool MultipleVectorHamiltonianMultiplyOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __MPI__
  if (architecture->IsMasterNode())
    {
      if (architecture->RequestOperation(this->OperationType) == false)
	{
	  return false;
	}
      if (this->RealDestinationVectors != 0)
	{
	  architecture->BroadcastVectorArray(this->NbrVectors, this->RealSourceVectors);
	  architecture->BroadcastVectorTypeArray(this->NbrVectors, this->RealDestinationVectors);  
	}
      else
	{
	  architecture->BroadcastVectorArray(this->NbrVectors, this->ComplexSourceVectors);
	  architecture->BroadcastVectorTypeArray(this->NbrVectors, this->ComplexDestinationVectors);  
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
      sprintf (TmpString, "MultipleVectorHamiltonianMultiply core operation done in %.3f seconds", Dt);
      architecture->AddToLog(TmpString);
    }
  if (this->RealDestinationVectors != 0)
    {
      for (int i = 0; i < this->NbrVectors; ++i)
	architecture->SumVector(this->RealDestinationVectors[i]);
    }
  else
    {
      for (int i = 0; i < this->NbrVectors; ++i)
	architecture->SumVector(this->ComplexDestinationVectors[i]);
    }
  if (architecture->IsMasterNode() == false)
    {
      if (this->RealDestinationVectors != 0)
	{
	  delete[] this->RealDestinationVectors;
	  delete[] this->RealSourceVectors;
	}
      else
	{
	  delete[] this->ComplexDestinationVectors;
	  delete[] this->ComplexSourceVectors;
	}
    }
  return true;
#else
  return this->RawApplyOperation();
#endif
}
  

