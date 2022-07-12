////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                     spin and generic 3-body interaction                    //
//                                                                            //
//                        last modification : 26/08/2008                      //
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
#include "Hamiltonian/AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Hamiltonian/ParticleOnSphereWithSpinL2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"

  
#include <stdio.h>
#include <iostream>
#include <stdlib.h>


using std::cout;
using std::endl;


// destructor
//

AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing::~AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing()
{
}


// need to overload multiplication routines to add new terms:

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored
RealVector& AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
												int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
      if (this->FullTwoBodyFlag == true)
	for (int i = firstComponent; i < LastComponent; ++i)
	  this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
      if (this->NbrPairingSectorSums != 0)
	for (int i = firstComponent; i < LastComponent; ++i)
	  this->EvaluatePairingAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
      for (int k = 3; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true)
	 {
	   for (int i = firstComponent; i < LastComponent; ++i)
	     this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination, k) ;
	 }	
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSource, vDestination);
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      Coefficient = vSource[k];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      vDestination[k++] += this->HamiltonianShift * Coefficient;
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
	      int* TmpIndexArray;
	      double* TmpCoefficientArray; 
	      int j;
	      int TmpNbrInteraction;
	      firstComponent -= this->PrecalculationShift;
	      LastComponent -= this->PrecalculationShift;
	      int Pos = firstComponent / this->FastMultiplicationStep; 
	      int PosMod = firstComponent % this->FastMultiplicationStep;
	      if (PosMod != 0)
		{
		  ++Pos;
		  PosMod = this->FastMultiplicationStep - PosMod;
		}
	      int l =  PosMod + firstComponent + this->PrecalculationShift;
	      for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
		{
		  TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
		  TmpIndexArray = this->InteractionPerComponentIndex[Pos];
		  TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
		  Coefficient = vSource[l];
		  for (j = 0; j < TmpNbrInteraction; ++j)
		    vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
		  vDestination[l] += this->HamiltonianShift * Coefficient;
		  l += this->FastMultiplicationStep;
		  ++Pos;
		}
	      firstComponent += this->PrecalculationShift;
	      LastComponent += this->PrecalculationShift;
	      for (l = 0; l < this->FastMultiplicationStep; ++l)
		if (PosMod != l)
		  {	
		    
		    if (this->FullTwoBodyFlag == true)
		      for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
		    for (int k = 2; k <= this->MaxNBody; ++k)
		      if (this->NBodyFlags[k] == true)
			{
			  for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			    this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination, k) ;
			}
		    this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSource, vDestination);
		  }
	      delete TmpParticles;
	    }
	  else
	    this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);	  
	}
    }
  if (this->L2Hamiltonian != 0)
    this->L2Hamiltonian->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
  if (this->S2Hamiltonian != 0)
    this->S2Hamiltonian->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
  return vDestination;

}




// multiply a et of vectors by the current hamiltonian for a given range of indices
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      double* Coefficient2 = new double [nbrVectors];
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
      if (this->FullTwoBodyFlag == true)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2) ;
	}
      if (this->NbrPairingSectorSums != 0)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    this->EvaluatePairingAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2) ;
	}
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true)
	  {
	    for (int i = firstComponent; i < LastComponent; ++i)
	      this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, k, Coefficient2);
	  }	    
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSources, vDestinations, nbrVectors);
      delete TmpParticles;
      delete[] Coefficient2;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  double* Coefficient2 = new double [nbrVectors];
	  int j;
	  int Pos;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  Coefficient2[l] = vSources[l][k];
		  vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
		}
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		    }
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
	      int* TmpIndexArray;
	      double* TmpCoefficientArray; 
	      double* Coefficient2 = new double [nbrVectors];
	      int j;
	      int Pos2;
	      int TmpNbrInteraction;
	      firstComponent -= this->PrecalculationShift;
	      LastComponent -= this->PrecalculationShift;
	      int Pos = firstComponent / this->FastMultiplicationStep; 
	      int PosMod = firstComponent % this->FastMultiplicationStep;
	      if (PosMod != 0)
		{
		  ++Pos;
		  PosMod = this->FastMultiplicationStep - PosMod;
		}
	      int l =  PosMod + firstComponent + this->PrecalculationShift;
	      for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
		{
		  TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
		  TmpIndexArray = this->InteractionPerComponentIndex[Pos];
		  TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      Coefficient2[k] = vSources[k][l];
		    }
		  for (j = 0; j < TmpNbrInteraction; ++j)
		    {
		      Pos2 = TmpIndexArray[j];
		      Coefficient = TmpCoefficientArray[j];
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][Pos2] += Coefficient  * Coefficient2[k];
			}
		    }
		  l += this->FastMultiplicationStep;
		  ++Pos;
		}
	      firstComponent += this->PrecalculationShift;
	      LastComponent += this->PrecalculationShift;
	      for (l = 0; l < this->FastMultiplicationStep; ++l)
		if (PosMod != l)
		  {	
		    if (this->FullTwoBodyFlag == true)
		      {
			for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			  this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2) ;
 		      }
		    for (int k = 2; k <= this->MaxNBody; ++k)
		      if (this->NBodyFlags[k] == true)
			{
			  for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			    this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, k, Coefficient2) ;
			}
		    this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSources, vDestinations, nbrVectors);
		  }
	      delete[] Coefficient2;
	      delete TmpParticles;
	    }
	  else
	    this->LowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	}
    }
  if (this->L2Hamiltonian != 0)
    this->L2Hamiltonian->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
  if (this->S2Hamiltonian != 0)
    this->S2Hamiltonian->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
  return vDestinations;
}
    


// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element
long AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  long Memory = 0;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int Dim = TmpParticles->GetHilbertSpaceDimension();
  int LastComponent = lastComponent + firstComponent;
  int Index;
  double Coefficient;
  
  this->EvaluateMNOneBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);
  
  if (this->FullTwoBodyFlag == true)
    this->EvaluateMNTwoBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);
  if (this->NbrPairingSectorSums != 0)
    this->EvaluatePairingFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);

  for (int k = 3; k <= this->MaxNBody; ++k)
    if (this->NBodyFlags[k] == true)
      EvaluateMNNBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, k, Memory);
  delete TmpParticles;
  
  return Memory;

}
  
// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
void AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{
  int LastComponent = nbrComponent + firstComponent;

  int* TmpIndexArray;
  double* TmpCoefficientArray;
  long ColumnIndex;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int Dim = TmpParticles->GetHilbertSpaceDimension();

  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  long Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
      ColumnIndex = 0l;
      if (this->FullTwoBodyFlag == true)
	this->EvaluateMNTwoBodyFastMultiplicationComponent(TmpParticles, i, TmpIndexArray, TmpCoefficientArray, ColumnIndex);
      if (this->NbrPairingSectorSums != 0)
	this->EvaluatePairingFastMultiplicationComponent(TmpParticles, i, TmpIndexArray, TmpCoefficientArray, ColumnIndex);
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true)
	  {
	    this->EvaluateMNNBodyFastMultiplicationComponent(TmpParticles, i, k, TmpIndexArray, TmpCoefficientArray, ColumnIndex);
	  }
      this->EvaluateMNOneBodyFastMultiplicationComponent(TmpParticles, i, TmpIndexArray, TmpCoefficientArray, ColumnIndex);
      ++Pos;
    }
  delete TmpParticles;

}
