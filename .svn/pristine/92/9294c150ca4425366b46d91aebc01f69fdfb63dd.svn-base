////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Antoine Sterdyniak                //
//                                                                            //
//                                                                            //
//                  class of division by a Jastrow operation                  //
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
#include "Architecture/ArchitectureOperation/FQHESphereJastrowDivisionOperation.h"
#include "Vector/RealVector.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "HilbertSpace/BosonOnSphereShort.h"

#include <iostream>
#include <sys/time.h>
#include <stdlib.h>

using std::cout;
using std::endl;


// constructor 
//
// space = pointer to the HilbertSpace to use
// sourceVector = array of vectors describing the fermionic states
// destinationVector = array of vectors where the resulting bosonic states have to be stored
// nbrStates = number of states to handle

FQHESphereJastrowDivisionOperation::FQHESphereJastrowDivisionOperation(BosonOnSphereShort* space, RealVector* sourceVector, RealVector* destinationVector, int nbrStates)
{
  this->FirstComponent = 0;
  this->NbrComponent = space->GetHilbertSpaceDimension();
  this->Space =(BosonOnSphereShort*) space->Clone();
  this->SourceVector = sourceVector;
  this->DestinationVector = destinationVector;
  this->OperationType = AbstractArchitectureOperation::FQHESphereJastrowDivision;
  this->NbrStates = nbrStates;
}


// copy constructor 
//
// operation = reference on operation to copy

FQHESphereJastrowDivisionOperation::FQHESphereJastrowDivisionOperation(const FQHESphereJastrowDivisionOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->Space = (BosonOnSphereShort*) operation.Space->Clone();
  this->SourceVector = operation.SourceVector;
  this->DestinationVector = operation.DestinationVector;
  this->OperationType = AbstractArchitectureOperation::FQHESphereJastrowDivision;
  this->NbrStates = operation.NbrStates;
}

// destructor
//

FQHESphereJastrowDivisionOperation::~FQHESphereJastrowDivisionOperation()
{
  delete this->Space;
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHESphereJastrowDivisionOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void FQHESphereJastrowDivisionOperation::SetDestinationVector (RealVector* destinationVector)
{
  this->DestinationVector = destinationVector;
}


// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHESphereJastrowDivisionOperation::Clone()
{
  return new FQHESphereJastrowDivisionOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESphereJastrowDivisionOperation::RawApplyOperation()
{
//   timeval TotalStartingTime;
//   gettimeofday (&TotalStartingTime, 0);
  this->Space->DivideByJastrow(this->SourceVector,this->DestinationVector,this->FirstComponent, this->NbrComponent,this->NbrStates);
//   timeval TotalEndingTime;
//   gettimeofday (&TotalEndingTime, 0);
//   double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
// 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
//   cout << this->FirstComponent << " " <<  this->NbrComponent << " : " << Dt << "s" << endl;
  return true;
}



// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereJastrowDivisionOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step;
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHESphereJastrowDivisionOperation** TmpOperations = new FQHESphereJastrowDivisionOperation* [architecture->GetNbrThreads()];
  unsigned long Sum;
  unsigned long Sigma = 0x0ul;
  for (int i= this->FirstComponent; i < this->NbrComponent; i++)
    {
      Sigma += this->TimeEvaluationFunction(i);
    }
  Sigma = Sigma / architecture->GetNbrThreads();
	
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      Step = 0;
      Sum = 0;
      while(Sum < Sigma)
	{
	  Sum += this->TimeEvaluationFunction(TmpFirstComponent + Step);
	  Step++;
	}
      TmpOperations[i] = (FQHESphereJastrowDivisionOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      TmpFirstComponent += Step;
    }
  TmpOperations[ReducedNbrThreads] = (FQHESphereJastrowDivisionOperation*) this->Clone();
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->FirstComponent + this->NbrComponent - TmpFirstComponent);
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      RealVector * TmpVector=new RealVector[this->NbrStates];
      for(int k = 0; k < this->NbrStates; k++)
	TmpVector[k] = RealVector(this->Space->GetHilbertSpaceDimension(), true);
      TmpOperations[i]->SetDestinationVector(TmpVector);
    }
  architecture->SendJobs();
  for (int i = 1; i < architecture->GetNbrThreads(); i++)
    {
      for(int k = 0; k < this->NbrStates; k++)
	this->DestinationVector[k] += TmpOperations[i]->DestinationVector[k];
      delete[] TmpOperations[i]->DestinationVector;
      delete TmpOperations[i];
    }

  delete TmpOperations[0];
  delete[] TmpOperations;
  return true;
}


