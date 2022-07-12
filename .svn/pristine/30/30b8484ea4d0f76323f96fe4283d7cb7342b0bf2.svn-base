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
#include "Vector/PartialRealVector.h"
#include "Vector/PartialComplexVector.h"
#include "Architecture/SMPArchitecture.h"


#include <sys/time.h>
#include <stdlib.h>


// constructor for real vectors
//
// hamiltonian = pointer to the hamiltonian to use
// sourceVectors = array of vectors to be multiplied by the hamiltonian
// destinationVectors = array of vectors where the result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// forceNormalMultiplication = force to use direct hamiltonian multiplication (override hamiltonian option)
// forceConjugateMultiplication = force to use hamiltonian for multiplication (override hamiltonian option)
// forceHermitianMultiplication = force to use hermitian property when multiplying  (override hamiltonian option)

MultipleVectorHamiltonianMultiplyOperation::MultipleVectorHamiltonianMultiplyOperation (AbstractHamiltonian* hamiltonian, RealVector* sourceVectors, 
											RealVector* destinationVectors, int nbrVectors, bool forceNormalMultiplication, bool forceConjugateMultiplication, bool forceHermitianeMultiplication)
{
  this->FirstComponent = 0;
  this->NbrComponent = sourceVectors[0].GetVectorDimension();
  this->Hamiltonian = hamiltonian;
  this->RealSourceVectors = sourceVectors;
  this->RealDestinationVectors = destinationVectors;
  this->ComplexSourceVectors = 0;
  this->ComplexDestinationVectors = 0; 
  this->RealSourcePartialVectors = 0;
  this->RealDestinationPartialVectors = 0;
  this->ComplexSourcePartialVectors = 0;
  this->ComplexDestinationPartialVectors = 0; 
  this->NbrVectors = nbrVectors;
  this->OperationType = AbstractArchitectureOperation::MultipleVectorHamiltonianMultiply;
  this->ExecutionTime = 0.0;
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

// constructor for complex vectors
//
// hamiltonian = pointer to the hamiltonian to use
// sourceVectors = array of vectors to be multiplied by the hamiltonian
// destinationVectors = array of vectors where the result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// forceNormalMultiplication = force to use direct hamiltonian multiplication (override hamiltonian option)
// forceConjugateMultiplication = force to use hamiltonian for multiplication (override hamiltonian option)
// forceHermitianMultiplication = force to use hermitian property when multiplying  (override hamiltonian option)

MultipleVectorHamiltonianMultiplyOperation::MultipleVectorHamiltonianMultiplyOperation (AbstractHamiltonian* hamiltonian, ComplexVector* sourceVectors, 
											ComplexVector* destinationVectors, int nbrVectors, bool forceNormalMultiplication, bool forceConjugateMultiplication, bool forceHermitianeMultiplication)
{

  this->FirstComponent = 0;
  this->NbrComponent = sourceVectors[0].GetVectorDimension();
  this->Hamiltonian = hamiltonian;
  this->RealSourceVectors = 0;
  this->RealDestinationVectors = 0;
  this->ComplexSourceVectors = sourceVectors;
  this->ComplexDestinationVectors = destinationVectors; 
  this->RealSourcePartialVectors = 0;
  this->RealDestinationPartialVectors = 0;
  this->ComplexSourcePartialVectors = 0;
  this->ComplexDestinationPartialVectors = 0; 
  this->NbrVectors = nbrVectors;
  this->OperationType = AbstractArchitectureOperation::MultipleVectorHamiltonianMultiply;
  this->ExecutionTime=0.0;
  // what is the following line for?
  // PartialComplexVector* DummyComplexDestinationPartialVectors = new PartialComplexVector[this->NbrVectors];
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

MultipleVectorHamiltonianMultiplyOperation::MultipleVectorHamiltonianMultiplyOperation(const MultipleVectorHamiltonianMultiplyOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->Hamiltonian = operation.Hamiltonian;
  this->RealSourceVectors = operation.RealSourceVectors;
  this->RealDestinationVectors = operation.RealDestinationVectors;  
  this->ComplexSourceVectors = operation.ComplexSourceVectors;
  this->ComplexDestinationVectors = operation.ComplexDestinationVectors;  
  this->RealSourcePartialVectors = operation.RealSourcePartialVectors;
  this->RealDestinationPartialVectors = operation.RealDestinationPartialVectors;
  this->ComplexSourcePartialVectors = operation.ComplexSourcePartialVectors;
  this->ComplexDestinationPartialVectors = operation.ComplexDestinationPartialVectors; 
  this->NbrVectors = operation.NbrVectors;
  this->UseConjugateFlag = operation.UseConjugateFlag;
  this->UseHermitianFlag = operation.UseHermitianFlag;
  this->ExecutionTime = operation.ExecutionTime;
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
  // copied from VectorHamiltonianMultiplyOperation, please check!
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

  Vector** TmpSourceVectors;
  Vector** TmpDestinationVectors;


  if (this->UseConjugateFlag == true)
    {
      TmpSourceVectors = architecture->BroadcastVectorArray(this->NbrVectors);  
      TmpDestinationVectors = architecture->ScatterVectorArray(this->NbrVectors);
    }
  else
    {
      if (this->UseHermitianFlag == true)
	{
	  TmpSourceVectors = architecture->BroadcastVectorArray(this->NbrVectors);  
	  TmpDestinationVectors = architecture->BroadcastVectorTypeArray(this->NbrVectors);
	}
      else
	{
	  // think one needs to scatter source vectors, rather than broadcast, here
	  // TmpSourceVectors = architecture->BroadcastVectorArray(this->NbrVectors);
	  TmpSourceVectors = architecture->ScatterVectorArray(this->NbrVectors);  
	  TmpDestinationVectors = architecture->BroadcastVectorTypeArray(this->NbrVectors);
	}
    }

  if (TmpSourceVectors[0]->GetVectorType() & Vector::RealDatas)
    {
      this->ComplexSourceVectors = 0;
      this->ComplexDestinationVectors = 0;       
      this->ComplexSourcePartialVectors = 0;
      this->ComplexDestinationPartialVectors = 0; 
      if ((TmpSourceVectors[0]->GetVectorType() & Vector::PartialData) != 0)
	{
	  this->RealSourceVectors = 0;
	  this->RealSourcePartialVectors = new PartialRealVector[this->NbrVectors];
	  for (int i = 0; i < this->NbrVectors; ++i)
	    {
	      this->RealSourcePartialVectors[i] = *((PartialRealVector*) TmpSourceVectors[i]);
	      delete TmpSourceVectors[i];
	    }
	}
      else
	{
	  this->RealSourcePartialVectors = 0;
	  this->RealSourceVectors = new RealVector[this->NbrVectors];
	  for (int i = 0; i < this->NbrVectors; ++i)
	    {
	      this->RealSourceVectors[i] = *((RealVector*) TmpSourceVectors[i]);
	      delete TmpSourceVectors[i];
	    }
	}
      if ((TmpDestinationVectors[0]->GetVectorType() & Vector::PartialData) != 0)
	{
	  this->RealDestinationVectors = 0;
	  this->RealDestinationPartialVectors = new PartialRealVector[this->NbrVectors];
	  for (int i = 0; i < this->NbrVectors; ++i)
	    {
	      this->RealDestinationPartialVectors[i] = *((PartialRealVector*) TmpDestinationVectors[i]);
	      delete TmpDestinationVectors[i];
	    }
	}
      else
	{
	  this->RealDestinationPartialVectors = 0;
	  this->RealDestinationVectors = new RealVector[this->NbrVectors];
	  for (int i = 0; i < this->NbrVectors; ++i)
	    {
	      this->RealDestinationVectors[i] = *((RealVector*) TmpDestinationVectors[i]);
	      delete TmpDestinationVectors[i];
	    }
	}
      delete[] TmpSourceVectors;
      delete[] TmpDestinationVectors;
    }
  else
    {
      this->RealSourceVectors = 0;
      this->RealDestinationVectors = 0;
      this->RealSourcePartialVectors = 0;
      this->RealDestinationPartialVectors = 0;
      if ((TmpSourceVectors[0]->GetVectorType() & Vector::PartialData) != 0)
	{
	  this->ComplexSourceVectors = 0;
	  this->ComplexSourcePartialVectors = new PartialComplexVector[this->NbrVectors];
	  for (int i = 0; i < this->NbrVectors; ++i)
	    {
	      this->ComplexSourcePartialVectors[i] = *((PartialComplexVector*) TmpSourceVectors[i]);
	      delete TmpSourceVectors[i];
	    }
	}
      else
	{
	  this->ComplexSourcePartialVectors = 0;
	  this->ComplexSourceVectors = new ComplexVector[this->NbrVectors];
	  for (int i = 0; i < this->NbrVectors; ++i)
	    {
	      this->ComplexSourceVectors[i] = *((ComplexVector*) TmpSourceVectors[i]);
	      delete TmpSourceVectors[i];
	    }
	}
      if ((TmpDestinationVectors[0]->GetVectorType() & Vector::PartialData) != 0)
	{
	  this->ComplexDestinationVectors = 0;
	  this->ComplexDestinationPartialVectors = new PartialComplexVector[this->NbrVectors];
	  for (int i = 0; i < this->NbrVectors; ++i)
	    {
	      this->ComplexDestinationPartialVectors[i] = *((PartialComplexVector*) (TmpDestinationVectors[i]));
	      delete TmpDestinationVectors[i];
	    }
	}
      else
	{
	  this->ComplexDestinationPartialVectors = 0;
	  this->ComplexDestinationVectors = new ComplexVector[this->NbrVectors];
	  for (int i = 0; i < this->NbrVectors; ++i)
	    {
	      this->ComplexDestinationVectors[i] = *((ComplexVector*) TmpDestinationVectors[i]);
	      delete TmpDestinationVectors[i];
	    }
	}
       delete[] TmpSourceVectors;
       delete[] TmpDestinationVectors;
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
  if ((destinationVectors[0].GetVectorType() & Vector::PartialData) != 0)
    {
      this->RealDestinationVectors = 0;
      this->RealDestinationPartialVectors = (PartialRealVector*) destinationVectors;
    }
  else
    {
      this->RealDestinationPartialVectors = 0;
      this->RealDestinationVectors = (RealVector*) destinationVectors;
    }
  this->ComplexDestinationVectors = 0;
}

// set destination vectors 
// 
// destinationVectors = array of vector where the result has to be stored

void MultipleVectorHamiltonianMultiplyOperation::SetDestinationVectors (ComplexVector* destinationVectors)
{
  if ((destinationVectors[0].GetVectorType() & Vector::PartialData) != 0)
    {
      this->ComplexDestinationVectors = 0;
      this->ComplexDestinationPartialVectors = (PartialComplexVector*) destinationVectors;
    }
  else
    {
      this->ComplexDestinationPartialVectors = 0;
      this->ComplexDestinationVectors = (ComplexVector*) destinationVectors;
    }
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
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  gettimeofday (&(TotalStartingTime2), 0);
  if (this->RealSourceVectors != 0)
    {
      if (this->UseConjugateFlag == true)
	{
	  this->Hamiltonian->ConjugateLowLevelMultipleMultiply(this->RealSourceVectors, this->RealDestinationVectors, this->NbrVectors,
							       this->FirstComponent, this->NbrComponent);
	}
      else if (this->UseHermitianFlag == true)
	{
	  this->Hamiltonian->HermitianLowLevelMultipleMultiply(this->RealSourceVectors, this->RealDestinationVectors,
							       this->NbrVectors, this->FirstComponent, this->NbrComponent);
	}
      else 
	{
	  this->Hamiltonian->LowLevelMultipleMultiply(this->RealSourceVectors, this->RealDestinationVectors, this->NbrVectors,
						      this->FirstComponent, this->NbrComponent);
	}
    }
  else
    {
      if (this->UseConjugateFlag == true)
	{
	  if (this->ComplexSourceVectors != 0)
	    {
	      if (this->ComplexDestinationVectors != 0)
		{
		  this->Hamiltonian->ConjugateLowLevelMultipleMultiply(this->ComplexSourceVectors, this->ComplexDestinationVectors, this->NbrVectors,
								       this->FirstComponent, this->NbrComponent);
		}
	      else
		{
		  this->Hamiltonian->ConjugateLowLevelMultipleMultiply(this->ComplexSourceVectors, this->ComplexDestinationPartialVectors, this->NbrVectors,
								       this->FirstComponent, this->NbrComponent);
		}
	    }
	  else
	    {
	      if (this->ComplexDestinationVectors != 0)
		{
		  this->Hamiltonian->ConjugateLowLevelMultipleMultiply(this->ComplexSourcePartialVectors, this->ComplexDestinationVectors, this->NbrVectors,
								       this->FirstComponent, this->NbrComponent);
		}
	      else
		{
		  this->Hamiltonian->ConjugateLowLevelMultipleMultiply(this->ComplexSourcePartialVectors, this->ComplexDestinationPartialVectors, this->NbrVectors,
								       this->FirstComponent, this->NbrComponent);
		}
	    }
	}
      else 
	{
	  if (this->UseHermitianFlag == true)
	    {
	      if (this->ComplexSourceVectors != 0)
		{
		  if (this->ComplexDestinationVectors != 0)
		    {
		      this->Hamiltonian->HermitianLowLevelMultipleMultiply(this->ComplexSourceVectors, this->ComplexDestinationVectors,
									   this->NbrVectors, this->FirstComponent, this->NbrComponent);
		    }
		  else
		    {
		      this->Hamiltonian->HermitianLowLevelMultipleMultiply(this->ComplexSourceVectors, this->ComplexDestinationPartialVectors,
									   this->NbrVectors, this->FirstComponent, this->NbrComponent);
		    }
		}
	      else
		{
		  if (this->ComplexDestinationVectors != 0)
		    {
		      this->Hamiltonian->HermitianLowLevelMultipleMultiply(this->ComplexSourcePartialVectors, this->ComplexDestinationVectors,
									   this->NbrVectors, this->FirstComponent, this->NbrComponent);
		    }
		  else
		    {
		      this->Hamiltonian->HermitianLowLevelMultipleMultiply(this->ComplexSourcePartialVectors, this->ComplexDestinationPartialVectors,
									   this->NbrVectors, this->FirstComponent, this->NbrComponent);
		    }
		}
	    }
	  else 
	    {
	      if (this->ComplexSourceVectors != 0)
		{
		  if (this->ComplexDestinationVectors != 0)
		    {
		      this->Hamiltonian->LowLevelMultipleMultiply(this->ComplexSourceVectors, this->ComplexDestinationVectors, this->NbrVectors, this->FirstComponent, 
								  this->NbrComponent);
		    }
		  else
		    {
		      this->Hamiltonian->LowLevelMultipleMultiply(this->ComplexSourceVectors, this->ComplexDestinationPartialVectors, this->NbrVectors, this->FirstComponent, 
								  this->NbrComponent);
		    }
		}
	      else
		{
		  if (this->ComplexDestinationVectors != 0)
		    {
		      this->Hamiltonian->LowLevelMultipleMultiply(this->ComplexSourcePartialVectors, this->ComplexDestinationVectors, this->NbrVectors, this->FirstComponent, 
								  this->NbrComponent);
		    }
		  else
		    {
		      this->Hamiltonian->LowLevelMultipleMultiply(this->ComplexSourcePartialVectors, this->ComplexDestinationPartialVectors, this->NbrVectors, this->FirstComponent, 
								  this->NbrComponent);
		    }
		}
	    }
	}
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

bool MultipleVectorHamiltonianMultiplyOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  bool RealFlag = false;
  if (this->ComplexDestinationVectors == 0)
    RealFlag = true;
  long *SegmentIndices=0;
  bool CleanUp = false;
  int TmpNbrThreads = architecture->GetNbrThreads();
  if (Hamiltonian->GetLoadBalancing(TmpNbrThreads, SegmentIndices) == false)
    {
      SegmentIndices = new long[TmpNbrThreads+1];
      CleanUp=true;
      int Step = this->NbrComponent / TmpNbrThreads;
      SegmentIndices[0] = this->FirstComponent;
      for (int i = 0; i < TmpNbrThreads; ++i)
	SegmentIndices[i] = this->FirstComponent+i*Step;
      SegmentIndices[TmpNbrThreads] = this->FirstComponent + this->NbrComponent;
    }

  MultipleVectorHamiltonianMultiplyOperation** TmpOperations = new MultipleVectorHamiltonianMultiplyOperation* [architecture->GetNbrThreads()];
  for (int i = 0; i < TmpNbrThreads; ++i)
    {
      TmpOperations[i] = (MultipleVectorHamiltonianMultiplyOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpOperations[i]->SetIndicesRange(SegmentIndices[i], SegmentIndices[i+1]-SegmentIndices[i]);
    }

  if (this->UseConjugateFlag == false)
    {
      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
	{
	  if (RealFlag == true)
	    TmpOperations[i]->SetDestinationVectors((RealVector*) this->RealDestinationVectors[0].EmptyCloneArray(this->NbrVectors, true));
	  else
	    TmpOperations[i]->SetDestinationVectors((ComplexVector*) this->ComplexDestinationVectors[0].EmptyCloneArray(this->NbrVectors, true));
	}
    }

  architecture->SendJobs();
  if (architecture->VerboseMode() == true)
    {
      char TmpString[512];
      double MinTime = 1e14;
      double MaxTime = 0.0;
      for (int i = 0; i < architecture->GetNbrThreads(); ++i)
	{
	  sprintf (TmpString, "MultipleVectorHamiltonianMultiply core operation on SMP id %d done in %.3f seconds", i, TmpOperations[i]->ExecutionTime);
	  if (MinTime > TmpOperations[i]->ExecutionTime)
	    {
	      MinTime = TmpOperations[i]->ExecutionTime;
	    }
	  if (MaxTime < TmpOperations[i]->ExecutionTime)
	    {
	      MaxTime = TmpOperations[i]->ExecutionTime;
	    }
	  architecture->AddToLog(TmpString);
	}
      sprintf (TmpString, "MultipleVectorHamiltonianMultiply core operation min time=%.3f sec, max time=%.3f sec", MinTime, MaxTime);
      architecture->AddToLog(TmpString);
   }
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      if (this->UseConjugateFlag == false)
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

bool MultipleVectorHamiltonianMultiplyOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
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
       if (this->RealDestinationVectors != 0)
	 {
	   for (int i = 0; i < this->NbrVectors; ++i)
	     this->RealDestinationVectors[i].ClearVector();
	 }
       if (this->ComplexDestinationVectors != 0)
	 {
	   for (int i = 0; i < this->NbrVectors; ++i)
	     this->ComplexDestinationVectors[i].ClearVector();
	 }
       if (this->RealDestinationPartialVectors != 0)
	 {
	   for (int i = 0; i < this->NbrVectors; ++i)
	     this->RealDestinationPartialVectors[i].ClearVector();
	 }
       if (this->ComplexDestinationPartialVectors != 0)
	 {
	   for (int i = 0; i < this->NbrVectors; ++i)
	     this->ComplexDestinationPartialVectors[i].ClearVector();
	 }

       if (this->UseConjugateFlag == true)
	 {
	   if (this->RealDestinationVectors != 0)
	     {
	       architecture->BroadcastVectorArray(this->NbrVectors, this->RealSourceVectors);
	       architecture->ScatterVectorArray(this->NbrVectors, this->RealDestinationVectors);  
	     }
	   else
	     {
	       architecture->BroadcastVectorArray(this->NbrVectors, this->ComplexSourceVectors);
	       architecture->ScatterVectorArray(this->NbrVectors, this->ComplexDestinationVectors);  
	     }
	 }
       else
	 {
	   if (this->UseHermitianFlag == true)
	     {
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
	   else
	     {
	       if (this->RealDestinationVectors != 0)
		 {
		   architecture->ScatterVectorArray(this->NbrVectors, this->RealSourceVectors);
		   architecture->BroadcastVectorTypeArray(this->NbrVectors, this->RealDestinationVectors);  
		 }
	       else
		 {
		   architecture->ScatterVectorArray(this->NbrVectors, this->ComplexSourceVectors);
		   architecture->BroadcastVectorTypeArray(this->NbrVectors, this->ComplexDestinationVectors);  
		 }
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
  if ((architecture->IsMasterNode()) && (architecture->VerboseMode()))
    gettimeofday (&TotalStartingTime, 0);
  if (this->UseConjugateFlag == true)
    {
      if (this->RealDestinationVectors != 0)
	{
	  for (int i = 0; i < this->NbrVectors; ++i)
	    architecture->ReassembleVector(this->RealDestinationVectors[i]);
	}
      else
	{
	  for (int i = 0; i < this->NbrVectors; ++i)
	    architecture->ReassembleVector(this->ComplexDestinationVectors[i]);
	}
    }
  else
    {
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
    }
  if ((architecture->IsMasterNode()) && (architecture->VerboseMode()))
    {
      timeval TotalEndingTime;
      gettimeofday (&TotalEndingTime, 0);
      double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
		    (((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
      char TmpString[256];
      sprintf (TmpString, "MultipleVectorHamiltonianMultiply sum operation done in %.3f seconds", Dt);
      architecture->AddToLog(TmpString, true);
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
  

