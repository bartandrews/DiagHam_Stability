////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of operations that apply a one body               //
//                        transformation to a many body state                 //
//                                                                            //
//                        last modification : 24/11/2015                      //
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
#include "Architecture/ArchitectureOperation/FQHESphereWithSpinApplyOneBodyTransformationOperation.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"



#include <iostream>
#include <sys/time.h>
#include <stdlib.h>

using std::cout;
using std::endl;


// constructor 
//
// inputState = vector where the initial state is stored
// outputState = vector where the rotated state is stored
// rotationMatrices =  matrices describing the one-body tranformation per orbital
// inputSpace = pointer to the Hilbert space

FQHESphereWithSpinApplyOneBodyTransformationOperation::FQHESphereWithSpinApplyOneBodyTransformationOperation(RealVector* inputState, RealVector* outputState, 
													     RealMatrix* rotationMatrices, 
													     ParticleOnSphereWithSpin* inputSpace)
{
  this->InputSpace =  (ParticleOnSphereWithSpin*) inputSpace->Clone();
  this->RealInputState = inputState;
  this->RealOutputState = outputState;
  this->RealRotationMatrices = rotationMatrices;
  this->InputState = 0;
  this->OutputState = 0;
  this->RotationMatrices = 0;
  this->FirstComponent = 0l;
  this->NbrComponent = this->InputSpace->GetLargeHilbertSpaceDimension();
  this->OperationType = AbstractArchitectureOperation::FQHESphereWithSpinApplyOneBodyTransformationOperation;
}

// constructor 
//
// inputState = vector where the initial state is stored
// outputState = vector where the rotated state is stored
// rotationMatrices =  matrices describing the one-body tranformation per orbital
// inputSpace = pointer to the Hilbert space

FQHESphereWithSpinApplyOneBodyTransformationOperation::FQHESphereWithSpinApplyOneBodyTransformationOperation(ComplexVector* inputState, ComplexVector* outputState, 
													     ComplexMatrix* rotationMatrices, 
													     ParticleOnSphereWithSpin* inputSpace)
{
  this->InputSpace =  (ParticleOnSphereWithSpin*) inputSpace->Clone();
  this->InputState = inputState;
  this->OutputState = outputState;
  this->RotationMatrices = rotationMatrices;
  this->RealInputState = 0;
  this->RealOutputState = 0;
  this->RealRotationMatrices = 0;
  this->FirstComponent = 0l;
  this->NbrComponent = this->InputSpace->GetLargeHilbertSpaceDimension();
  this->OperationType = AbstractArchitectureOperation::FQHESphereWithSpinApplyOneBodyTransformationOperation;
}

// copy constructor 
//
// operation = reference on operation to copy

FQHESphereWithSpinApplyOneBodyTransformationOperation::FQHESphereWithSpinApplyOneBodyTransformationOperation(const FQHESphereWithSpinApplyOneBodyTransformationOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->InputSpace =  (ParticleOnSphereWithSpin*) operation.InputSpace->Clone();
  this->InputState = operation.InputState;
  this->OutputState = operation.OutputState;
  this->RotationMatrices = operation.RotationMatrices;
  this->RealInputState = operation.RealInputState;
  this->RealOutputState = operation.RealOutputState;
  this->RealRotationMatrices = operation.RealRotationMatrices;
  this->OperationType = AbstractArchitectureOperation::FQHESphereWithSpinApplyOneBodyTransformationOperation;	
}

// destructor
//

FQHESphereWithSpinApplyOneBodyTransformationOperation::~FQHESphereWithSpinApplyOneBodyTransformationOperation()
{
  delete this->InputSpace;
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHESphereWithSpinApplyOneBodyTransformationOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void FQHESphereWithSpinApplyOneBodyTransformationOperation::SetDestinationVector (ComplexVector* destinationVector)
{
  this->OutputState = destinationVector;
}

// set destination vector 
// 
// vector where the result has to be stored

void FQHESphereWithSpinApplyOneBodyTransformationOperation::SetDestinationVector (RealVector* destinationVector)
{
  this->RealOutputState = destinationVector;
}


// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHESphereWithSpinApplyOneBodyTransformationOperation::Clone()
{
  return new FQHESphereWithSpinApplyOneBodyTransformationOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESphereWithSpinApplyOneBodyTransformationOperation::RawApplyOperation()
{
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);
  if (this->RealOutputState == 0)
    {
      this->InputSpace->TransformOneBodyBasis(*(this->InputState), *(this->OutputState), this->RotationMatrices, this->FirstComponent, this->NbrComponent);
    }
  else
    {
      this->InputSpace->TransformOneBodyBasis(*(this->RealInputState), *(this->RealOutputState), this->RealRotationMatrices, this->FirstComponent, this->NbrComponent);
    }
  timeval TotalEndingTime;
  gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout << "contribution from " << this->FirstComponent << " to " <<  (this->NbrComponent + this->FirstComponent - 1) << " computed in " << Dt << "s" << endl;
  return true;
}



// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereWithSpinApplyOneBodyTransformationOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  long Step = this->NbrComponent / ((long) architecture->GetNbrThreads());
  long TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHESphereWithSpinApplyOneBodyTransformationOperation** TmpOperations = new FQHESphereWithSpinApplyOneBodyTransformationOperation * [architecture->GetNbrThreads()];  
   
  for( int i = 0; i <  architecture->GetNbrThreads() ; i++)
    {
      TmpOperations[i] = (FQHESphereWithSpinApplyOneBodyTransformationOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }

  if (this->RealOutputState == 0)
    {
      for( int i = 1; i <  architecture->GetNbrThreads() ; i++)
	TmpOperations[i]->SetDestinationVector((ComplexVector*) this->OutputState->EmptyClone(true));
    }
  else
    {
      for( int i = 1; i <  architecture->GetNbrThreads() ; i++)
	TmpOperations[i]->SetDestinationVector((RealVector*) this->RealOutputState->EmptyClone(true));
    }
  for( int i = 0; i <  ReducedNbrThreads ; i++)
    {
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      TmpFirstComponent += Step;
    }

  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->FirstComponent + this->NbrComponent - TmpFirstComponent);
 

  architecture->SendJobs();

  for (int i = 1; i < architecture->GetNbrThreads(); i++)
    {
      if (this->RealOutputState == 0)
	{
	  (*this->OutputState) += (*TmpOperations[i]->OutputState);
	}
      else
	{
	  (*this->RealOutputState) += (*TmpOperations[i]->RealOutputState);
	}	
      delete  TmpOperations[i]->OutputState;
      delete TmpOperations[i];
    }

  delete TmpOperations[0];
  delete[] TmpOperations;
  return true;
}
