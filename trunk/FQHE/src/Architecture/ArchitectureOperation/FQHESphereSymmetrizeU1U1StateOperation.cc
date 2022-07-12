////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Antoine Sterdyniak                //
//                                                                            //
//                                                                            //
//                   class of U1U1 states symmetrization Operation	      //
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
#include "Architecture/ArchitectureOperation/FQHESphereSymmetrizeU1U1StateOperation.h"
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
// finalSpace = pointer to the Hilbert space of the target space
// leftSpace = pointer to the Hilbert space of the first state
// rightSpace = pointer to the Hilbert space of the second state
// destinationVector = vector where the result has to be stored
// leftVector = vector that contains the first state
// rightVector = vector that contains the second state
// unnormalizedBasisFlag = true if the states are expressed in the unnormalized basis

FQHESphereSymmetrizeU1U1StateOperation::FQHESphereSymmetrizeU1U1StateOperation(BosonOnSphereShort* finalSpace,  BosonOnSphereShort* leftSpace, BosonOnSphereShort* rightSpace, 
									       RealVector* destinationVector , RealVector* leftVector, RealVector* rightVector, 
									       bool unnormalizedBasisFlag)
{
  this->FirstComponent = 0;
  this->NbrComponent = leftSpace->GetHilbertSpaceDimension();
  this->FinalSpace = (BosonOnSphereShort*) finalSpace->Clone();
  this->LeftSpace = (BosonOnSphereShort*) leftSpace->Clone();
  this->RightSpace = (BosonOnSphereShort*) rightSpace->Clone();
  this->FinalSpaceWithSpin = 0;
  this->LeftSpaceWithSpin = 0;
  this->RightSpaceWithSpin = 0;
  this->LeftVector = leftVector;
  this->RightVector = rightVector;
  this->DestinationVector = destinationVector;
  this->RationalLeftVector = 0;
  this->RationalRightVector = 0;
  this->RationalDestinationVector = 0;
  this->UnnormalizedBasisFlag = unnormalizedBasisFlag;
  this->OperationType = AbstractArchitectureOperation::FQHESphereSymmetrizeU1U1StateOperation;
}

// constructor for long rational vector input
//
// finalSpace = pointer to the Hilbert space of the target space
// leftSpace = pointer to the Hilbert space of the first state
// rightSpace = pointer to the Hilbert space of the second state
// destinationVector = vector where the result has to be stored
// leftVector = vector that contains the first state
// rightVector = vector that contains the second state

FQHESphereSymmetrizeU1U1StateOperation::FQHESphereSymmetrizeU1U1StateOperation(BosonOnSphereShort* finalSpace, BosonOnSphereShort* leftSpace , BosonOnSphereShort* rightSpace, 
									       LongRationalVector* destinationVector, LongRationalVector* leftVector, LongRationalVector* rightVector)
{
  this->FirstComponent = 0;
  this->NbrComponent = leftSpace->GetHilbertSpaceDimension();
  this->FinalSpace = (BosonOnSphereShort*) finalSpace->Clone();
  this->LeftSpace = (BosonOnSphereShort*) leftSpace->Clone();
  this->RightSpace = (BosonOnSphereShort*) rightSpace->Clone();
  this->FinalSpaceWithSpin = 0;
  this->LeftSpaceWithSpin = 0;
  this->RightSpaceWithSpin = 0;
  this->LeftVector = 0;
  this->RightVector = 0;
  this->DestinationVector = 0;
  this->RationalLeftVector = leftVector;
  this->RationalRightVector = rightVector;
  this->RationalDestinationVector = destinationVector;
  this->UnnormalizedBasisFlag = true;
  this->OperationType = AbstractArchitectureOperation::FQHESphereSymmetrizeU1U1StateOperation;
}

// constructor for spinful states
//
// finalSpace = pointer to the Hilbert space of the target space
// leftSpace = pointer to the Hilbert space of the first state
// rightSpace = pointer to the Hilbert space of the second state
// destinationVector = vector where the result has to be stored
// leftVector = vector that contains the first state
// rightVector = vector that contains the second state
// unnormalizedBasisFlag = true if the states are expressed in the unnormalized basis

FQHESphereSymmetrizeU1U1StateOperation::FQHESphereSymmetrizeU1U1StateOperation(ParticleOnSphereWithSpin* finalSpace, ParticleOnSphereWithSpin* leftSpace, ParticleOnSphereWithSpin* rightSpace, 
									       RealVector* destinationVector , RealVector* leftVector, RealVector* rightVector, 
									       bool unnormalizedBasisFlag)
{
  this->FirstComponent = 0;
  this->NbrComponent = leftSpace->GetHilbertSpaceDimension();
  this->FinalSpaceWithSpin = (ParticleOnSphereWithSpin*) finalSpace->Clone();
  this->LeftSpaceWithSpin = (ParticleOnSphereWithSpin*) leftSpace->Clone();
  this->RightSpaceWithSpin = (ParticleOnSphereWithSpin*) rightSpace->Clone();
  this->FinalSpace = 0;
  this->LeftSpace = 0;
  this->RightSpace = 0;
  this->LeftVector = leftVector;
  this->RightVector = rightVector;
  this->DestinationVector = destinationVector;
  this->RationalLeftVector = 0;
  this->RationalRightVector = 0;
  this->RationalDestinationVector = 0;
  this->UnnormalizedBasisFlag = unnormalizedBasisFlag;
  this->OperationType = AbstractArchitectureOperation::FQHESphereSymmetrizeU1U1StateOperation;
}

// copy constructor 
//
// operation = reference on operation to copy

FQHESphereSymmetrizeU1U1StateOperation::FQHESphereSymmetrizeU1U1StateOperation(const FQHESphereSymmetrizeU1U1StateOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  
  if (operation.FinalSpace != 0)
    {
      this->FinalSpace = (BosonOnSphereShort*) operation.FinalSpace->Clone();
      this->LeftSpace =  (BosonOnSphereShort*) operation.LeftSpace->Clone();
      this->RightSpace = (BosonOnSphereShort*) operation.RightSpace->Clone();
      this->FinalSpaceWithSpin = 0;
      this->LeftSpaceWithSpin = 0;
      this->RightSpaceWithSpin = 0;      
    }
  else
    {
      this->FinalSpaceWithSpin = (ParticleOnSphereWithSpin*) operation.FinalSpaceWithSpin->Clone();
      this->LeftSpaceWithSpin = (ParticleOnSphereWithSpin*) operation.LeftSpaceWithSpin->Clone();
      this->RightSpaceWithSpin = (ParticleOnSphereWithSpin*) operation.RightSpaceWithSpin->Clone();
      this->FinalSpace = 0;
      this->LeftSpace = 0;
      this->RightSpace = 0;
    }

  this->LeftVector = operation.LeftVector;
  this->RightVector = operation.RightVector;
  this->DestinationVector = operation.DestinationVector;
  this->RationalLeftVector = operation.RationalLeftVector;
  this->RationalRightVector = operation.RationalRightVector;
  this->RationalDestinationVector = operation.RationalDestinationVector;
  this->UnnormalizedBasisFlag = operation.UnnormalizedBasisFlag;
  this->OperationType = AbstractArchitectureOperation::FQHESphereSymmetrizeU1U1StateOperation;	
}

// destructor
//

FQHESphereSymmetrizeU1U1StateOperation::~FQHESphereSymmetrizeU1U1StateOperation()
{
  if (this->FinalSpace != 0)
    {
      delete this->FinalSpace;
      delete this->LeftSpace;
      delete this->RightSpace;
    }
  else
    {
      delete this->FinalSpaceWithSpin;
      delete this->LeftSpaceWithSpin;
      delete this->RightSpaceWithSpin;
    }
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHESphereSymmetrizeU1U1StateOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void FQHESphereSymmetrizeU1U1StateOperation::SetDestinationVector (RealVector* destinationVector)
{
  this->DestinationVector = destinationVector;
}


// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHESphereSymmetrizeU1U1StateOperation::Clone()
{
  return new FQHESphereSymmetrizeU1U1StateOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESphereSymmetrizeU1U1StateOperation::RawApplyOperation()
{
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);
  
  if (this->FinalSpace != 0)
    {
      if (this->DestinationVector != 0)
	{
	  this->FinalSpace->SymmetrizeU1U1StateCore (*this->DestinationVector, (*this->LeftVector), (*this->RightVector), 
						     this->LeftSpace,  this->RightSpace, this->UnnormalizedBasisFlag, this->FirstComponent, this->NbrComponent);
	}
      else
	{
	  this->FinalSpace->SymmetrizeU1U1StateCore (*this->RationalDestinationVector ,(*this->RationalLeftVector) , (*this->RationalRightVector),  
						     this->LeftSpace,  this->RightSpace, this->FirstComponent, this->NbrComponent);
	}
    }
  else
    {
      if (this->DestinationVector != 0)
	{
	  this->FinalSpaceWithSpin->SymmetrizeSU2SU2StateCore (*this->DestinationVector, (*this->LeftVector), (*this->RightVector), 
							       this->LeftSpaceWithSpin,  this->RightSpaceWithSpin, this->UnnormalizedBasisFlag, this->FirstComponent, this->NbrComponent);
	}      
    }
  
  
  timeval TotalEndingTime;
  gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout << this->FirstComponent << " " <<  this->NbrComponent << " : " << Dt << "s" << endl;
  return true;
}



// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereSymmetrizeU1U1StateOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrComponent/architecture->GetNbrThreads();
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHESphereSymmetrizeU1U1StateOperation** TmpOperations = new FQHESphereSymmetrizeU1U1StateOperation* [architecture->GetNbrThreads()];
  
   
  for( int i = 0; i <  architecture->GetNbrThreads() ; i++)
    {
      TmpOperations[i] = (FQHESphereSymmetrizeU1U1StateOperation *) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  
  for( int i = 1; i <  architecture->GetNbrThreads() ; i++)
    TmpOperations[i]->SetDestinationVector((RealVector*)this->DestinationVector->EmptyClone(true));

  for( int i = 0; i <  ReducedNbrThreads ; i++)
    {
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      TmpFirstComponent += Step;
    }

  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->FirstComponent + this->NbrComponent - TmpFirstComponent);
 

  architecture->SendJobs();

  for (int i = 1; i < architecture->GetNbrThreads(); i++)
    {
      (*this->DestinationVector) += (*TmpOperations[i]->DestinationVector);
      delete  TmpOperations[i]->DestinationVector;
      delete TmpOperations[i];
    }

  delete TmpOperations[0];
  delete[] TmpOperations;
  return true;
}
