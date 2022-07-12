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
#include "Architecture/ArchitectureOperation/FQHELatticeFourierTransformOperation.h"
#include "Vector/ComplexVector.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/BosonOnLatticeKy.h"
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

FQHELatticeFourierTransformOperation::FQHELatticeFourierTransformOperation (BosonOnLattice * fullSpace, BosonOnLatticeKy * kySpace, ComplexVector* sourceVector, int nbrComponent,ComplexVector* destinationVector, int nbrVectors)
{
  this->FirstComponent = 0;
  this->NbrComponent = nbrComponent;
  this->FullSpace =(BosonOnLattice*) fullSpace->Clone();
  this->KySpace = (BosonOnLatticeKy*) kySpace->Clone();
  this->SourceVector = sourceVector;
  this->DestinationVector = destinationVector;
  this->OperationType = AbstractArchitectureOperation::FQHELatticeFourierTransformOperation;
  this->NbrVectors = nbrVectors;
}


// copy constructor 
//
// operation = reference on operation to copy

FQHELatticeFourierTransformOperation::FQHELatticeFourierTransformOperation(const FQHELatticeFourierTransformOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->FullSpace = (BosonOnLattice*) operation.FullSpace->Clone();
  this->KySpace = (BosonOnLatticeKy*) operation.KySpace->Clone();
  this->SourceVector = operation.SourceVector;
  this->DestinationVector = operation.DestinationVector;
  this->OperationType = AbstractArchitectureOperation::FQHELatticeFourierTransformOperation;
  this->NbrVectors = operation.NbrVectors;
}

// destructor
//

FQHELatticeFourierTransformOperation::~FQHELatticeFourierTransformOperation()
{
  delete this->FullSpace;
  delete this->KySpace;
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHELatticeFourierTransformOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void FQHELatticeFourierTransformOperation::SetDestinationVector (ComplexVector* destinationVector)
{
  this->DestinationVector = destinationVector;
}


// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHELatticeFourierTransformOperation::Clone()
{
  return new FQHELatticeFourierTransformOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHELatticeFourierTransformOperation::RawApplyOperation()
{
   timeval TotalStartingTime;
      gettimeofday (&TotalStartingTime, 0);
  this->KySpace->ConvertFromNbodyBasis( this->SourceVector,this->DestinationVector, *this->FullSpace,this->NbrVectors, this->FirstComponent, this->NbrComponent);

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

bool FQHELatticeFourierTransformOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrComponent/architecture->GetNbrThreads();
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHELatticeFourierTransformOperation** TmpOperations = new FQHELatticeFourierTransformOperation* [architecture->GetNbrThreads()];
	
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      TmpOperations[i] = (FQHELatticeFourierTransformOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      TmpFirstComponent += Step;
    }
  TmpOperations[ReducedNbrThreads] = (FQHELatticeFourierTransformOperation*) this->Clone();
  TmpOperations[ReducedNbrThreads]->SetDestinationVector((ComplexVector *) this->DestinationVector->EmptyClone(true));
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->FirstComponent + this->NbrComponent - TmpFirstComponent);
  
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      ComplexVector * TmpVector=new ComplexVector[this->NbrVectors];
      for(int k = 0; k < this->NbrVectors; k++)
	TmpVector[k] = ComplexVector(this->KySpace->GetHilbertSpaceDimension(), true);
      TmpOperations[i]->SetDestinationVector(TmpVector);
    }
  architecture->SendJobs();
  
  for (int i = 1; i < architecture->GetNbrThreads(); i++)
    {
      for(int k = 0; k < this->NbrVectors; k++)
	this->DestinationVector[k] += TmpOperations[i]->DestinationVector[k];
      delete[] TmpOperations[i]->DestinationVector;
      delete TmpOperations[i];
    }
    
  delete TmpOperations[0];
  delete[] TmpOperations;
   return true;
}

