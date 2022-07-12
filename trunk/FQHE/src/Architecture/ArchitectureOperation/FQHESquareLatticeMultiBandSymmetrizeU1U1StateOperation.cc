////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Antoine Sterdyniak                //
//                                                                            //
//                                                                            //
//                   class of U1U1 states symmetrization Operation            //
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
#include "Architecture/ArchitectureOperation/FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation.h"
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
// space = pointer to the HilbertSpace to use
// sourceVector = array of vectors describing the fermionic states
// destinationVector = array of vectors where the resulting bosonic states have to be stored
// nbrStates = number of states to handle

FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation::FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation( BosonOnSquareLatticeMultiBandMomentumSpace * finalSpace,  BosonOnSquareLatticeMultiBandMomentumSpace * leftSpace, BosonOnSquareLatticeMultiBandMomentumSpace * rightSpace, ComplexVector * destinationVector , ComplexVector * leftVector, ComplexVector * rightVector ,bool unnormalizedBasisFlag)
{
  this->FirstComponent = 0;
  this->NbrComponent = leftSpace->GetHilbertSpaceDimension();
  this->FinalSpace = finalSpace;
  this->LeftSpace = leftSpace;
  this->RightSpace = rightSpace;
  this->FinalSpaceLong = 0;
  this->LeftSpaceLong = 0;
  this->RightSpaceLong = 0;
  this->LeftVector = leftVector;
  this->RightVector = rightVector;
  this->DestinationVector = destinationVector;
  this->UnnormalizedBasisFlag = unnormalizedBasisFlag;
  this->OperationType = AbstractArchitectureOperation::FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation;
}

// constructor using long Hilbert spaces
//
// space = pointer to the HilbertSpace to use
// sourceVector = array of vectors describing the fermionic states
// destinationVector = array of vectors where the resulting bosonic states have to be stored
// nbrStates = number of states to handle

FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation::FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation(BosonOnSquareLatticeMomentumSpaceLong* finalSpace,  BosonOnSquareLatticeMomentumSpaceLong* leftSpace, BosonOnSquareLatticeMomentumSpaceLong* rightSpace, ComplexVector * destinationVector , ComplexVector * leftVector, ComplexVector * rightVector ,bool unnormalizedBasisFlag)
{
  this->FirstComponent = 0;
  this->NbrComponent = leftSpace->GetHilbertSpaceDimension();
  this->FinalSpaceLong = finalSpace;
  this->LeftSpaceLong = leftSpace;
  this->RightSpaceLong = rightSpace;
  this->FinalSpaceLong = 0;
  this->LeftSpaceLong = 0;
  this->RightSpaceLong = 0;
  this->LeftVector = leftVector;
  this->RightVector = rightVector;
  this->DestinationVector = destinationVector;
  this->UnnormalizedBasisFlag = unnormalizedBasisFlag;
  this->OperationType = AbstractArchitectureOperation::FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation;
}

// copy constructor 
//
// operation = reference on operation to copy

FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation::FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation(const FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;

  if (this->FinalSpace != 0)
    {
      this->FinalSpace = (BosonOnSquareLatticeMultiBandMomentumSpace*) operation.FinalSpace->Clone();
      this->LeftSpace =  (BosonOnSquareLatticeMultiBandMomentumSpace*) operation.LeftSpace->Clone();
      this->RightSpace = (BosonOnSquareLatticeMultiBandMomentumSpace*) operation.RightSpace->Clone();
      this->FinalSpaceLong = 0;
      this->LeftSpaceLong = 0;
      this->RightSpaceLong = 0;
    }
  else
    {
      this->FinalSpaceLong = (BosonOnSquareLatticeMomentumSpaceLong*) operation.FinalSpaceLong->Clone();
      this->LeftSpaceLong =  (BosonOnSquareLatticeMomentumSpaceLong*) operation.LeftSpaceLong->Clone();
      this->RightSpaceLong = (BosonOnSquareLatticeMomentumSpaceLong*) operation.RightSpaceLong->Clone();
      this->FinalSpace = 0;
      this->LeftSpace = 0;
      this->RightSpace = 0;
    }
  this->LeftVector = operation.LeftVector;
  this->RightVector = operation.RightVector;
  this->DestinationVector = operation.DestinationVector;
  this->UnnormalizedBasisFlag = operation.UnnormalizedBasisFlag;
  this->OperationType = AbstractArchitectureOperation::FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation;	
}

// destructor
//

FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation::~FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation()
{
  //delete this->FinalSpace;
  //delete this->LeftSpace;
  //delete this->RightSpace;
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation::SetDestinationVector (ComplexVector* destinationVector)
{
  this->DestinationVector = destinationVector;
}


// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation::Clone()
{
  return new FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation::RawApplyOperation()
{
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);

  if (this->FinalSpace != 0)
    {
      this->FinalSpace->SymmetrizeU1U1StateCore ( *this->DestinationVector ,(*this->LeftVector) , (*this->RightVector) ,  LeftSpace,  RightSpace , this->UnnormalizedBasisFlag, this->FirstComponent, this->NbrComponent);
    }
  else
    {
      this->FinalSpaceLong->SymmetrizeU1U1StateCore ( *this->DestinationVector ,(*this->LeftVector) , (*this->RightVector) ,  LeftSpaceLong,  RightSpaceLong, this->UnnormalizedBasisFlag, this->FirstComponent, this->NbrComponent);
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

bool FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrComponent/architecture->GetNbrThreads();
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation** TmpOperations = new FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation * [architecture->GetNbrThreads()];
  
   
  for( int i = 0; i <  architecture->GetNbrThreads() ; i++)
    {
      TmpOperations[i] = (FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation *) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  
  for( int i = 1; i <  architecture->GetNbrThreads() ; i++)
    TmpOperations[i]->SetDestinationVector((ComplexVector*)this->DestinationVector->EmptyClone(true));

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
