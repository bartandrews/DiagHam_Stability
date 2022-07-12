////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2002 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//     class of U(1) states symmetrization operation for the torus geometry   //
//                                                                            //
//                        last modification : 11/11/2015                      //
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
#include "Architecture/ArchitectureOperation/FQHETorusSymmetrizeU1U1StateOperation.h"
#include "Vector/RealVector.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "HilbertSpace/ParticleOnTorus.h"

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

FQHETorusSymmetrizeU1U1StateOperation::FQHETorusSymmetrizeU1U1StateOperation(ParticleOnTorus* finalSpace,  ParticleOnTorus* leftSpace, ParticleOnTorus* rightSpace, 
									     RealVector* destinationVector , RealVector* leftVector, RealVector* rightVector)
{
  this->FirstComponent = 0;
  this->NbrComponent = leftSpace->GetHilbertSpaceDimension();
  this->FinalSpace = (ParticleOnTorus*) finalSpace->Clone();
  this->LeftSpace = (ParticleOnTorus*) leftSpace->Clone();
  this->RightSpace = (ParticleOnTorus*) rightSpace->Clone();
  this->LeftVector = leftVector;
  this->RightVector = rightVector;
  this->DestinationVector = destinationVector;
  this->ComplexLeftVector = 0;
  this->ComplexRightVector = 0;
  this->ComplexDestinationVector = 0;
  this->OperationType = AbstractArchitectureOperation::FQHETorusSymmetrizeU1U1StateOperation;
}

// constructor for long complex vector input
//
// finalSpace = pointer to the Hilbert space of the target space
// leftSpace = pointer to the Hilbert space of the first state
// rightSpace = pointer to the Hilbert space of the second state
// destinationVector = vector where the result has to be stored
// leftVector = vector that contains the first state
// rightVector = vector that contains the second state

FQHETorusSymmetrizeU1U1StateOperation::FQHETorusSymmetrizeU1U1StateOperation(ParticleOnTorus* finalSpace, ParticleOnTorus* leftSpace , ParticleOnTorus* rightSpace, 
									     ComplexVector* destinationVector, ComplexVector* leftVector, ComplexVector* rightVector)
{
  this->FirstComponent = 0;
  this->NbrComponent = leftSpace->GetHilbertSpaceDimension();
  this->FinalSpace = (ParticleOnTorus*) finalSpace->Clone();
  this->LeftSpace = (ParticleOnTorus*) leftSpace->Clone();
  this->RightSpace = (ParticleOnTorus*) rightSpace->Clone();
  this->LeftVector = 0;
  this->RightVector = 0;
  this->DestinationVector = 0;
  this->ComplexLeftVector = leftVector;
  this->ComplexRightVector = rightVector;
  this->ComplexDestinationVector = destinationVector;
  this->OperationType = AbstractArchitectureOperation::FQHETorusSymmetrizeU1U1StateOperation;
}

// copy constructor 
//
// operation = reference on operation to copy

FQHETorusSymmetrizeU1U1StateOperation::FQHETorusSymmetrizeU1U1StateOperation(const FQHETorusSymmetrizeU1U1StateOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;

  this->FinalSpace = (ParticleOnTorus*) operation.FinalSpace->Clone();
  this->LeftSpace =  (ParticleOnTorus*) operation.LeftSpace->Clone();
  this->RightSpace = (ParticleOnTorus*) operation.RightSpace->Clone();
  this->LeftVector = operation.LeftVector;
  this->RightVector = operation.RightVector;
  this->DestinationVector = operation.DestinationVector;
  this->ComplexLeftVector = operation.ComplexLeftVector;
  this->ComplexRightVector = operation.ComplexRightVector;
  this->ComplexDestinationVector = operation.ComplexDestinationVector;
  this->OperationType = AbstractArchitectureOperation::FQHETorusSymmetrizeU1U1StateOperation;	
}

// destructor
//

FQHETorusSymmetrizeU1U1StateOperation::~FQHETorusSymmetrizeU1U1StateOperation()
{
  delete this->FinalSpace;
  delete this->LeftSpace;
  delete this->RightSpace;
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHETorusSymmetrizeU1U1StateOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void FQHETorusSymmetrizeU1U1StateOperation::SetDestinationVector (RealVector* destinationVector)
{
  this->DestinationVector = destinationVector;
}


// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHETorusSymmetrizeU1U1StateOperation::Clone()
{
  return new FQHETorusSymmetrizeU1U1StateOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHETorusSymmetrizeU1U1StateOperation::RawApplyOperation()
{
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);
  
  if (this->DestinationVector != 0)
    {
      this->FinalSpace->SymmetrizeU1U1StateCore (*this->DestinationVector, (*this->LeftVector), (*this->RightVector), 
						 LeftSpace,  RightSpace, false, this->FirstComponent, this->NbrComponent);
    }
  else
    {
       this->FinalSpace->SymmetrizeU1U1StateCore (*this->ComplexDestinationVector ,(*this->ComplexLeftVector) , (*this->ComplexRightVector),  
						  LeftSpace,  RightSpace, false, this->FirstComponent, this->NbrComponent);
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

bool FQHETorusSymmetrizeU1U1StateOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrComponent/architecture->GetNbrThreads();
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHETorusSymmetrizeU1U1StateOperation** TmpOperations = new FQHETorusSymmetrizeU1U1StateOperation* [architecture->GetNbrThreads()];
  
   
  for( int i = 0; i <  architecture->GetNbrThreads() ; i++)
    {
      TmpOperations[i] = (FQHETorusSymmetrizeU1U1StateOperation*) this->Clone();
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
