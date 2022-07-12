////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Antoine Sterdyniak                //
//                                                                            //
//                                                                            //
//                   class of multiplication by Jastrow Operation	      //
//                                                                            //
//                        last modification : 03/03/2010                      //
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
#include "Architecture/ArchitectureOperation/FQHESphereJastrowMultiplicationOperation.h"
#include "Vector/RealVector.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "Architecture/MonoProcessorArchitecture.h"

#include <iostream>
#include <sys/time.h>
#include <stdlib.h>

using std::cout;
using std::endl;


// constructor 
//
// space = pointer to the HilbertSpace
// sourceVector = vector to be multiplied by the operator
// destinationVector = vector where the result has to be stored

FQHESphereJastrowMultiplicationOperation::FQHESphereJastrowMultiplicationOperation (BosonOnSphereShort* space, RealVector * sourceVector, RealVector* destinationVector,int nbrState)
{
  this->FirstComponent = 0;
  this->Space =(BosonOnSphereShort*) space->Clone();
  this->SourceVector = sourceVector;
  this->DestinationVector = destinationVector;
  this->JackVector=new RealVector(this->Space->GetHilbertSpaceDimension(),true);
  this->OperationType = AbstractArchitectureOperation::FQHESphereJastrowMultiplication;
  this->NbrState = nbrState;
}


// copy constructor 
//
// operation = reference on operation to copy

FQHESphereJastrowMultiplicationOperation::FQHESphereJastrowMultiplicationOperation(const FQHESphereJastrowMultiplicationOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->Space = (BosonOnSphereShort*) operation.Space->Clone();
  this->SourceVector = operation.SourceVector;
  this->DestinationVector = operation.DestinationVector;
  this->JackVector = new RealVector(this->Space->GetHilbertSpaceDimension(),true);
  this->OperationType = AbstractArchitectureOperation::FQHESphereJastrowMultiplication;
  this->NbrState = operation.NbrState;
}
  
// constructor from a master node information
//
// operator = pointer to the operator to use
// architecture = pointer to the distributed architecture to use for communications

FQHESphereJastrowMultiplicationOperation::FQHESphereJastrowMultiplicationOperation(BosonOnSphereShort* space, SimpleMPIArchitecture* architecture,int nbrState)
{
  this->Space =(BosonOnSphereShort*) space->Clone();
  this->OperationType = AbstractArchitectureOperation::FQHESphereJastrowMultiplication;
  long TmpMinimumIndex = 0;
  long TmpMaximumIndex=0;
  architecture->GetTypicalRange(TmpMinimumIndex, TmpMaximumIndex);
  this->FirstComponent = (int) TmpMinimumIndex;
  this->SourceVector = (RealVector*) architecture->ScatterVector();
  this->DestinationVector = (RealVector*) architecture->BroadcastVectorType();
  this->NbrState = nbrState;
}



// destructor
//

FQHESphereJastrowMultiplicationOperation::~FQHESphereJastrowMultiplicationOperation()
{
  delete this->Space;
  delete this->JackVector;
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHESphereJastrowMultiplicationOperation::SetFirstComponent (const long& firstComponent)
{
  this->FirstComponent = firstComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void FQHESphereJastrowMultiplicationOperation::SetDestinationVector (RealVector* destinationVector)
{
  this->DestinationVector = destinationVector;
}

// set Jack Vector
//
// vectors where the Jack Polynomial will be stored

void FQHESphereJastrowMultiplicationOperation::SetJackVector (RealVector* jackVector)
{
  this->JackVector = jackVector;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHESphereJastrowMultiplicationOperation::Clone()
{
  return new FQHESphereJastrowMultiplicationOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESphereJastrowMultiplicationOperation::RawApplyOperation()
{
  //timeval TotalStartingTime;
  //gettimeofday (&TotalStartingTime, 0);
	
  this->JackVector->ClearVector();
  this->Space->KostkaForOneSchur(this->FirstComponent,*(this->JackVector));
  
  /*timeval TotalEndingTime;
    gettimeofday (&TotalEndingTime, 0);
    double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
    (((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
    cout << this->FirstComponent << " " <<  this->NbrComponent << " : " << Dt << "s" << endl;*/
  return true;
}

// apply operation for MonoProcessor architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereJastrowMultiplicationOperation::ArchitectureDependentApplyOperation(MonoProcessorArchitecture* architecture)
{
  for (long i=0;i<this->Space->GetHilbertSpaceDimension();i++)
    {
      this->SetFirstComponent(i);
      this->JackVector->ClearVector();
      this->Space->KostkaForOneSchur(this->FirstComponent, *(this->JackVector));
      for(int k = 0; k< this->NbrState; k++)
	{
	  this->DestinationVector[k][i]=this->SourceVector[k][i];
	  this->SourceVector->AddLinearCombination(-this->DestinationVector[k][i], (*this->JackVector), i + 1,
						   this->Space->GetHilbertSpaceDimension() - i - 1);
	}
    }
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereJastrowMultiplicationOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  this->DestinationVector->ClearVector();
  long Step = 0;
  
  long Limit = this->Space->GetHilbertSpaceDimension() - architecture->GetNbrThreads();
  
  FQHESphereJastrowMultiplicationOperation** TmpOperations = new FQHESphereJastrowMultiplicationOperation* [architecture->GetNbrThreads()];
  for(int i = 0; i < architecture->GetNbrThreads();i++)
    {
      TmpOperations[i] = (FQHESphereJastrowMultiplicationOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  
  while(Step < Limit)
    {
      for(int i = 0; i <architecture->GetNbrThreads(); i++)
	{
	  TmpOperations[i]->SetFirstComponent(Step);
	  Step++;
	}		
      architecture->SendJobs();		
      for(long i = Step - architecture->GetNbrThreads(); i < Step; i++)
	{
	  int Tmp= i + architecture->GetNbrThreads() - Step;
	  for(int k = 0; k < this->NbrState; k++)
	    {
	      this->DestinationVector[k][i]=this->SourceVector[k][i];
	      this->SourceVector[k].AddLinearCombination(-this->DestinationVector[k][i], (*TmpOperations[Tmp]->JackVector), i + 1, 
							 this->Space->GetHilbertSpaceDimension() - i - 1);
	    }
	}
      
    }
	
  int Remainder = this->Space->GetHilbertSpaceDimension() - Step;  
  if(architecture->GetNbrThreads() < Remainder)
    cout << "Warning, too few tasks called in FQHESphereJastrowMultiplicationOperation" << endl;
  
  for(int i = 0; i < Remainder; i++)
    {
      TmpOperations[i]->SetFirstComponent(Step);
      Step++;
    }

  architecture->SendJobs();
  for(int i = Step - Remainder; i < Step; i++)
    {
      for(int k=0 ; k < this->NbrState; k++)
	this->DestinationVector[k][i] = this->SourceVector[k][i];
      if(i != (this->Space->GetHilbertSpaceDimension() - 1))
	{
	  int Tmp= i + Remainder - Step;
	  for(int k = 0; k < this->NbrState; k++)
	    this->SourceVector[k].AddLinearCombination(-this->DestinationVector[k][i], (*TmpOperations[Tmp]->JackVector), i + 1,
						       this->Space->GetHilbertSpaceDimension() - i - 1);
	}
    }
  
  for (int i = 0; i < architecture->GetNbrThreads(); i++)
    {
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  return true;
}




