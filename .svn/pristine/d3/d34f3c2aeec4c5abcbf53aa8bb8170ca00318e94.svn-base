////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of FQHE on disk quasihole propagator operation          //
//                                                                            //
//                        last modification : 08/03/2009                      //
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
#include "Architecture/ArchitectureOperation/FQHESphereJackGeneratorOperation.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"


// constructor 
//
// space = pointer to the Hilbert space to use
// invAlpha = inverse of the Jack polynomial alpha coefficient
// rootPartition = root partition (in fermionic binary representation)
// indexArray = array where state indices are stored
// stateArray = array use to store computed state description
// componentArray = array where computed component numerical factors are stored
// rhoArray = rho factor associated to each state
// nbrComputedComponentArray = number of connected components associated to each state through the Jack generator
// fermionicFlag = true if we are dealing with fermions
// symmetricFlag = true if the state is Lz<->-Lz symmetric

FQHESphereJackGeneratorOperation::FQHESphereJackGeneratorOperation (ParticleOnSphere* space, double invAlpha, unsigned long rootPartition, long** indexArray, unsigned long** stateArray, double** componentArray, double* rhoArray, int* nbrComputedComponentArray, bool fermionicFlag, bool symmetricFlag)
{
  this->InvAlpha = invAlpha;
  this->RootPartition = rootPartition;
  this->LocalShift = 0l;
  this->IndexArray = indexArray;
  this->StateArray = stateArray; 
  this->ComponentArray = componentArray;
  this->RhoArray = rhoArray;
  this->NbrComputedComponentArray = nbrComputedComponentArray;
  this->FirstComponent = 0;
  this->NbrComponent = 0;
  this->LargeFirstComponent = 0;
  this->LargeNbrComponent = space->GetLargeHilbertSpaceDimension();
  this->HilbertSpace = (ParticleOnSphere*) space->Clone();
  this->FermionicFlag = fermionicFlag;
  this->SymmetricFlag = symmetricFlag;
  this->OperationType = AbstractArchitectureOperation::FQHESphereJackGenerator;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}

// copy constructor 
//
// operation = reference on operation to copy

FQHESphereJackGeneratorOperation::FQHESphereJackGeneratorOperation(const FQHESphereJackGeneratorOperation& operation)
{
  this->InvAlpha = operation.InvAlpha;
  this->RootPartition = operation.RootPartition;
  this->LocalShift = 0l;
  this->IndexArray = operation.IndexArray;
  this->StateArray = operation.StateArray; 
  this->ComponentArray = operation.ComponentArray;
  this->RhoArray = operation.RhoArray;
  this->NbrComputedComponentArray = operation.NbrComputedComponentArray;
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->LargeFirstComponent = operation.LargeFirstComponent;
  this->LargeNbrComponent = operation.LargeNbrComponent;
  this->HilbertSpace = (ParticleOnSphere*) operation.HilbertSpace->Clone();
  this->FermionicFlag = operation.FermionicFlag;
  this->SymmetricFlag = operation.SymmetricFlag;
  this->OperationType = AbstractArchitectureOperation::FQHESphereJackGenerator;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}
  
// destructor
//

FQHESphereJackGeneratorOperation::~FQHESphereJackGeneratorOperation()
{
  if (this->LocalOperations != 0)
    {
      for (int i = 0; i < this->NbrLocalOperations; ++i)
	{
	  delete this->LocalOperations[i];
	}
      delete[] this->LocalOperations;
    }
  delete this->HilbertSpace;
}
  
// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHESphereJackGeneratorOperation::Clone()
{
  return new FQHESphereJackGeneratorOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESphereJackGeneratorOperation::RawApplyOperation()
{
  if (this->FermionicFlag == false)
    {
      if (this->SymmetricFlag == true)
	{
	  ((BosonOnSphereHaldaneHugeBasisShort*) this->HilbertSpace)->GenerateSymmetrizedJackPolynomialFactorizedCore(this->InvAlpha, this->RootPartition, this->LargeFirstComponent, this->LargeFirstComponent + this->LargeNbrComponent - 1l, this->StateArray + this->LocalShift, this->ComponentArray + this->LocalShift, this->IndexArray + this->LocalShift, this->NbrComputedComponentArray + this->LocalShift, this->RhoArray + this->LocalShift);
	}
      else
	{
	  ((BosonOnSphereHaldaneHugeBasisShort*) this->HilbertSpace)->GenerateJackPolynomialFactorizedCore(this->InvAlpha, this->RootPartition, this->LargeFirstComponent, this->LargeFirstComponent + this->LargeNbrComponent - 1l, this->StateArray + this->LocalShift, this->ComponentArray + this->LocalShift, this->IndexArray + this->LocalShift, this->NbrComputedComponentArray + this->LocalShift, this->RhoArray + this->LocalShift);
	}
    }
  else
    {
      if (this->SymmetricFlag == true)
	{
	  ((FermionOnSphereHaldaneHugeBasis*) this->HilbertSpace)->GenerateSymmetrizedJackPolynomialFactorizedCore(this->InvAlpha, this->RootPartition, this->LargeFirstComponent, this->LargeFirstComponent + this->LargeNbrComponent - 1l, this->StateArray + this->LocalShift, this->ComponentArray + this->LocalShift, this->IndexArray + this->LocalShift, this->NbrComputedComponentArray + this->LocalShift, this->RhoArray + this->LocalShift);
	}
      else
	{
	  ((FermionOnSphereHaldaneHugeBasis*) this->HilbertSpace)->GenerateJackPolynomialFactorizedCore(this->InvAlpha, this->RootPartition, this->LargeFirstComponent, this->LargeFirstComponent + this->LargeNbrComponent - 1l, this->StateArray + this->LocalShift, this->ComponentArray + this->LocalShift, this->IndexArray + this->LocalShift, this->NbrComputedComponentArray + this->LocalShift, this->RhoArray + this->LocalShift);
	}
    }
  return true;
}


// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereJackGeneratorOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  long Step = this->LargeNbrComponent / ((long) architecture->GetNbrThreads());
  long TmpFirstComponent = this->LargeFirstComponent;
  if (this->LocalOperations == 0)
    {
      this->NbrLocalOperations = architecture->GetNbrThreads();
      this->LocalOperations = new FQHESphereJackGeneratorOperation* [this->NbrLocalOperations];
      for (int i = 0; i < this->NbrLocalOperations; ++i)
	this->LocalOperations[i] = (FQHESphereJackGeneratorOperation*) this->Clone();
    }
  int ReducedNbrThreads = this->NbrLocalOperations - 1;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->LocalOperations[i]->LocalShift = TmpFirstComponent - this->LargeFirstComponent;
      this->LocalOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      architecture->SetThreadOperation(this->LocalOperations[i], i);
      TmpFirstComponent += Step;
    }
  this->LocalOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->LargeNbrComponent + this->LargeFirstComponent - TmpFirstComponent);  
  this->LocalOperations[ReducedNbrThreads]->LocalShift = TmpFirstComponent - this->LargeFirstComponent;
  architecture->SetThreadOperation(this->LocalOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs();
  return true;
}
  
