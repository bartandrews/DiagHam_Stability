////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of operations that apply a Cn rotation             //
//                                                                            //
//                        last modification : 18/05/2014                      //
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
#include "Architecture/ArchitectureOperation/FQHETorusApplyCNRotationOperation.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"



#include <iostream>
#include <sys/time.h>
#include <stdlib.h>

using std::cout;
using std::endl;


// constructor 
//
// nValue = N Value of the rotation (negative if clockwise)
// inputState = pointer to the Hilbert space of the initial state
// outputState = pointer to the Hilbert space of the rotated state
// inputSpace = vector where the initial state is stored
// outputSpace = vector where the rotated state is stored

FQHETorusApplyCNRotationOperation::FQHETorusApplyCNRotationOperation(int nValue, ComplexVector* inputState, ComplexVector* outputState, ParticleOnTorus* inputSpace, ParticleOnTorus* outputSpace)
{
  this->NValue = nValue;
  this->InputSpace =  (ParticleOnTorus*) inputSpace->Clone();
  this->OutputSpace = (ParticleOnTorus*) outputSpace->Clone();
  this->InputSpaceWithMagneticTranslations = 0;
  this->OutputSpaceWithMagneticTranslations = 0;
  this->InputState = inputState;
  this->OutputState = outputState;
  this->FirstComponent = 0;
  this->NbrComponent = this->OutputSpace->GetHilbertSpaceDimension();
  this->OperationType = AbstractArchitectureOperation::FQHETorusApplyCNRotationOperation;
}

// constructor 
//
// nValue = N Value of the rotation (negative if clockwise)
// inputState = pointer to the Hilbert space of the initial state
// outputState = pointer to the Hilbert space of the rotated state
// inputSpace = vector where the initial state is stored
// outputSpace = vector where the rotated state is stored

FQHETorusApplyCNRotationOperation::FQHETorusApplyCNRotationOperation(int nValue, ComplexVector* inputState, ComplexVector* outputState, ParticleOnTorusWithMagneticTranslations* inputSpace, ParticleOnTorusWithMagneticTranslations* outputSpace)
{
  this->NValue = nValue;
  this->InputSpaceWithMagneticTranslations =  (ParticleOnTorusWithMagneticTranslations*) inputSpace->Clone();
  this->OutputSpaceWithMagneticTranslations = (ParticleOnTorusWithMagneticTranslations*) outputSpace->Clone();
  this->InputSpace = 0;
  this->OutputSpace = 0;
  this->InputState = inputState;
  this->OutputState = outputState;
  this->FirstComponent = 0;
  this->NbrComponent = this->OutputSpaceWithMagneticTranslations->GetHilbertSpaceDimension();
  this->OperationType = AbstractArchitectureOperation::FQHETorusApplyCNRotationOperation;
}
    

// copy constructor 
//
// operation = reference on operation to copy

FQHETorusApplyCNRotationOperation::FQHETorusApplyCNRotationOperation(const FQHETorusApplyCNRotationOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  if (operation.InputSpaceWithMagneticTranslations == 0)
    {
      this->InputSpace =  (ParticleOnTorus*) operation.InputSpace->Clone();
      this->OutputSpace = (ParticleOnTorus*) operation.OutputSpace->Clone();
      this->InputSpaceWithMagneticTranslations = 0;
      this->OutputSpaceWithMagneticTranslations = 0;
    }
  else
    {
      this->InputSpaceWithMagneticTranslations =  (ParticleOnTorusWithMagneticTranslations*) operation.InputSpaceWithMagneticTranslations->Clone();
      this->OutputSpaceWithMagneticTranslations = (ParticleOnTorusWithMagneticTranslations*) operation.OutputSpaceWithMagneticTranslations->Clone();
      this->InputSpace = 0;
      this->OutputSpace = 0;
    }
  this->InputState = operation.InputState;
  this->OutputState = operation.OutputState;
  this->NValue = operation.NValue;
  this->OperationType = AbstractArchitectureOperation::FQHETorusApplyCNRotationOperation;	
}

// destructor
//

FQHETorusApplyCNRotationOperation::~FQHETorusApplyCNRotationOperation()
{
  if (this->InputSpaceWithMagneticTranslations == 0)
    {
      delete this->InputSpace;
      delete this->OutputSpace;
    }
  else
    {
      delete this->InputSpaceWithMagneticTranslations;
      delete this->OutputSpaceWithMagneticTranslations;
    }
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHETorusApplyCNRotationOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->FirstComponent = (int) firstComponent;
  this->NbrComponent = (int) nbrComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void FQHETorusApplyCNRotationOperation::SetDestinationVector (ComplexVector* destinationVector)
{
  this->OutputState = destinationVector;
}


// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHETorusApplyCNRotationOperation::Clone()
{
  return new FQHETorusApplyCNRotationOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHETorusApplyCNRotationOperation::RawApplyOperation()
{
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);  
  if (this->InputSpaceWithMagneticTranslations == 0)
    {
      if (this->NValue > 0)
	this->OutputSpace->CoreC4Rotation(*(this->InputState), this->InputSpace, *(this->OutputState), this->FirstComponent, this->NbrComponent, false);
      else
	this->OutputSpace->CoreC4Rotation(*(this->InputState), this->InputSpace, *(this->OutputState), this->FirstComponent, this->NbrComponent, true);
    }
  else
    {
      if (this->NValue > 0)
	this->OutputSpaceWithMagneticTranslations->CoreC4Rotation(*(this->InputState), this->InputSpaceWithMagneticTranslations, *(this->OutputState), this->FirstComponent, this->NbrComponent, false);
      else
	this->OutputSpaceWithMagneticTranslations->CoreC4Rotation(*(this->InputState), this->InputSpaceWithMagneticTranslations, *(this->OutputState), this->FirstComponent, this->NbrComponent, true);
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

bool FQHETorusApplyCNRotationOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrComponent / architecture->GetNbrThreads();
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHETorusApplyCNRotationOperation** TmpOperations = new FQHETorusApplyCNRotationOperation * [architecture->GetNbrThreads()];  
   
  for( int i = 0; i <  architecture->GetNbrThreads() ; i++)
    {
      TmpOperations[i] = (FQHETorusApplyCNRotationOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  
  for( int i = 1; i <  architecture->GetNbrThreads() ; i++)
    TmpOperations[i]->SetDestinationVector((ComplexVector*) this->OutputState->EmptyClone(true));

  for( int i = 0; i <  ReducedNbrThreads ; i++)
    {
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      TmpFirstComponent += Step;
    }

  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->FirstComponent + this->NbrComponent - TmpFirstComponent);
 

  architecture->SendJobs();

  for (int i = 1; i < architecture->GetNbrThreads(); i++)
    {
      (*this->OutputState) += (*TmpOperations[i]->OutputState);
      delete  TmpOperations[i]->OutputState;
      delete TmpOperations[i];
    }

  delete TmpOperations[0];
  delete[] TmpOperations;
  return true;
}
