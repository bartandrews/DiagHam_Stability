////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    class of FQHE particle entanglement matrix parallelization operation    //
//                                                                            //
//                        last modification : 23/12/2016                      //
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
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementMatrixOperation.h"

#include <sys/time.h>


// constructor 
//
// fullSpace = pointer to the full Hilbert space to use
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// groundState = reference on the total system ground state
// entanglementMatrix = reference on the entanglement matrix where result has to stored
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut

FQHESphereParticleEntanglementMatrixOperation::FQHESphereParticleEntanglementMatrixOperation(ParticleOnSphere* fullSpace, ParticleOnSphere* destinationSpace, 
											     ParticleOnSphere* complementarySpace, RealVector& groundState, 
											     RealMatrix& entanglementMatrix, bool removeBinomialCoefficient)
{
  this->FullSpace  = (ParticleOnSphere*) fullSpace->Clone();
  this->DestinationHilbertSpace = (ParticleOnSphere*) destinationSpace->Clone();
  this->ComplementaryHilbertSpace = (ParticleOnSphere*) complementarySpace->Clone();
  this->GroundState = groundState;
  this->EntanglementMatrix = entanglementMatrix;
  this->RemoveBinomialCoefficientFlag = removeBinomialCoefficient;
  this->NbrNonZeroElements = 0l;
  this->FirstComponent = 0;
  this->NbrComponent = 0;
  this->LargeFirstComponent = 0;
  this->LargeNbrComponent = this->ComplementaryHilbertSpace->GetLargeHilbertSpaceDimension();
  this->OperationType = AbstractArchitectureOperation::FQHESphereParticleEntanglementMatrixOperation;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}

// copy constructor 
//
// operation = reference on operation to copy

FQHESphereParticleEntanglementMatrixOperation::FQHESphereParticleEntanglementMatrixOperation(const FQHESphereParticleEntanglementMatrixOperation& operation)
{
  this->FullSpace  = (ParticleOnSphere*) operation.FullSpace->Clone();
  this->DestinationHilbertSpace = (ParticleOnSphere*) operation.DestinationHilbertSpace->Clone();
  this->ComplementaryHilbertSpace = (ParticleOnSphere*) operation.ComplementaryHilbertSpace->Clone();
  this->GroundState = operation.GroundState;
  this->RemoveBinomialCoefficientFlag = operation.RemoveBinomialCoefficientFlag;
  this->EntanglementMatrix = operation.EntanglementMatrix;
  this->NbrNonZeroElements = operation.NbrNonZeroElements;
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->LargeFirstComponent = operation.LargeFirstComponent;
  this->LargeNbrComponent = operation.LargeNbrComponent;
  this->OperationType = AbstractArchitectureOperation::FQHESphereParticleEntanglementMatrixOperation;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}
  
// destructor
//

FQHESphereParticleEntanglementMatrixOperation::~FQHESphereParticleEntanglementMatrixOperation()
{
  if (this->LocalOperations != 0)
    {
      for (int i = 0; i < this->NbrLocalOperations; ++i)
	{
	  delete this->LocalOperations[i];
	}
      delete[] this->LocalOperations;
    }
  delete this->FullSpace;
  delete this->DestinationHilbertSpace;
  delete this->ComplementaryHilbertSpace;
}
  
// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHESphereParticleEntanglementMatrixOperation::Clone()
{
  return new FQHESphereParticleEntanglementMatrixOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESphereParticleEntanglementMatrixOperation::RawApplyOperation()
{
  this->NbrNonZeroElements = this->FullSpace->EvaluatePartialEntanglementMatrixParticlePartitionCore(this->LargeFirstComponent, this->LargeNbrComponent, this->ComplementaryHilbertSpace, 
												     this->DestinationHilbertSpace, this->GroundState, &this->EntanglementMatrix,
												     this->RemoveBinomialCoefficientFlag);
  return true;
}


// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereParticleEntanglementMatrixOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  long Step = this->LargeNbrComponent / ((long) architecture->GetNbrThreads());
  long TmpFirstComponent = this->LargeFirstComponent;
  if (this->LocalOperations == 0)
    {
      this->NbrLocalOperations = architecture->GetNbrThreads();
      this->LocalOperations = new FQHESphereParticleEntanglementMatrixOperation* [this->NbrLocalOperations];
      for (int i = 0; i < this->NbrLocalOperations; ++i)
	this->LocalOperations[i] = (FQHESphereParticleEntanglementMatrixOperation*) this->Clone();
    }
  int ReducedNbrThreads = this->NbrLocalOperations - 1;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->LocalOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      architecture->SetThreadOperation(this->LocalOperations[i], i);
      TmpFirstComponent += Step;
    }
  this->LocalOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->LargeNbrComponent + this->LargeFirstComponent - TmpFirstComponent);  
  architecture->SetThreadOperation(this->LocalOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs();
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->LocalOperations[ReducedNbrThreads]->NbrNonZeroElements += this->LocalOperations[i]->NbrNonZeroElements;
    }
  this->NbrNonZeroElements = this->LocalOperations[ReducedNbrThreads]->NbrNonZeroElements;
  return true;
}
  
