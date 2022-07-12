////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//                                                                            //
//                       class of FQHE spinless hamiltonian                   //
//                     with real coefficients and no symmetry                 //
//                                                                            //
//                                                                            //
//                        last modification : 23/10/2012                      //
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
#include "Hamiltonian/AbstractQHEOnSphereFullHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;



// default constructor
//

AbstractQHEOnSphereFullHamiltonian::AbstractQHEOnSphereFullHamiltonian()
{
  this->HermitianSymmetryFlag = true;
}

// destructor
//

AbstractQHEOnSphereFullHamiltonian::~AbstractQHEOnSphereFullHamiltonian()
{
  for (int i = 0; i < this->NbrSectorSums; ++i)
    {
      if (this->NbrSectorIndicesPerSum[i] > 0)
	{
	  delete[] this->SectorIndicesPerSum[i];
	  delete[] this->InteractionFactors[i];
	}
    }
  if (this->NbrSectorSums != 0)
   {
     delete[] this->InteractionFactors;
     delete[] this->NbrSectorIndicesPerSum;
     delete[] this->SectorIndicesPerSum;
   }
  if (this->OneBodyTermFlag == true)
    delete[] this->OneBodyInteractionFactors;
  if (this->FastMultiplicationFlag == true)
    {
      long MinIndex;
      long MaxIndex;
      this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
      int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
      int ReducedDim = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      if ((ReducedDim * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
	++ReducedDim;
      for (int i = 0; i < ReducedDim; ++i)
	{
	  if (this->NbrInteractionPerComponent[i] > 0)
	    {
	      delete[] this->InteractionPerComponentIndex[i];
	      delete[] this->InteractionPerComponentCoefficient[i];
	    }
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
    }
}
  
// ask if Hamiltonian implements hermitian symmetry operations
//

bool AbstractQHEOnSphereFullHamiltonian::IsHermitian()
{
  return this->HermitianSymmetryFlag;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void AbstractQHEOnSphereFullHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (ParticleOnSphere*) hilbertSpace;
  this->EvaluateInteractionFactors();
}
 
// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* AbstractQHEOnSphereFullHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int AbstractQHEOnSphereFullHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value

void AbstractQHEOnSphereFullHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereFullHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								    int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
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
	  double Coefficient;
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
	  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  double Coefficient;
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
		for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
		  this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
		this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSource, vDestination);
	      }
	  delete TmpParticles;
	}
    }
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

RealVector* AbstractQHEOnSphereFullHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									    int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      double* Coefficient2 = new double [nbrVectors];
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSources, vDestinations, nbrVectors);
      delete[] Coefficient2;
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double Coefficient;
	  double* Coefficient2 = new double [nbrVectors];
	  double* TmpCoefficientArray; 
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
		    vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	}
      else
	{
	  this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	}
    }
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereFullHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
  int j;
  int TmpNbrInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int Pos2;
  double Coefficient;
  int PosMod = firstComponent % this->FastMultiplicationStep;
  double* Coefficient2 = new double [nbrVectors];
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
	  vDestinations[k][l] += this->HamiltonianShift * Coefficient2[k];
	}
      for (j = 0; j < TmpNbrInteraction; ++j)
	{
	  Pos2 = TmpIndexArray[j];
	  Coefficient = TmpCoefficientArray[j];
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient2[k];
	}
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (l = 0; l < this->FastMultiplicationStep; ++l)
    if (PosMod != l)
      {	
	for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
	  this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
	this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSources, vDestinations, nbrVectors);
      }
  delete[] Coefficient2;
  delete TmpParticles;
  return vDestinations;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereFullHamiltonian::HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									     int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
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
	  double Coefficient;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      Coefficient = vSource[k];
	      double TmpSum = 0.0;
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  TmpSum += (TmpCoefficientArray[j]) * vSource[TmpIndexArray[j]];
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
		}
	      TmpSum += this->HamiltonianShift * Coefficient;
	      vDestination[k++] += TmpSum;
	    }
	}
      else
	{
	  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  double Coefficient;
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
	      double TmpSum = 0.0;
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
		  TmpSum += (TmpCoefficientArray[j]) * vSource[TmpIndexArray[j]];
		}
	      TmpSum += this->HamiltonianShift * Coefficient;	      
	      vDestination[l] += TmpSum;
	      l += this->FastMultiplicationStep;
	      ++Pos;
	    }
	  firstComponent += this->PrecalculationShift;
	  LastComponent += this->PrecalculationShift;
	  for (l = 0; l < this->FastMultiplicationStep; ++l)
	    if (PosMod != l)
	      {	
		for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
		  this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
		this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSource, vDestination);
	      }
	  delete TmpParticles;
	}
    }
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

RealVector* AbstractQHEOnSphereFullHamiltonian::HermitianLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										     int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      double* Coefficient2 = new double [nbrVectors];
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSources, vDestinations, nbrVectors);
      delete[] Coefficient2;
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double Coefficient;
	  double* Coefficient2 = new double [nbrVectors];
	  double* TmpCoefficientArray; 
	  int j;
	  int Pos;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  double* TmpSum = new double [nbrVectors];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  TmpSum[l] = 0.0;
		  Coefficient2[l] = vSources[l][k];
		}
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		      TmpSum[l] += (Coefficient) * vSources[l][Pos];
		    }
		}
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  TmpSum[l] += this->HamiltonianShift * Coefficient2[l];
		  vDestinations[l][k] += TmpSum[l];
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	  delete[] TmpSum;
	}
      else
	{
	  this->HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	}
    }
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereFullHamiltonian::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
													int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
  int j;
  int TmpNbrInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int Pos2;
  double Coefficient;
  int PosMod = firstComponent % this->FastMultiplicationStep;
  double* Coefficient2 = new double [nbrVectors];
  double* TmpSum = new double [nbrVectors];
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
	  TmpSum[k] = 0.0;
	  Coefficient2[k] = vSources[k][l];
	}
      for (j = 0; j < TmpNbrInteraction; ++j)
	{
	  Pos2 = TmpIndexArray[j];
	  Coefficient = TmpCoefficientArray[j];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient2[k];
	      TmpSum[k] += (Coefficient) * vSources[k][Pos2];
	    }
	}
     for (int k = 0; k < nbrVectors; ++k)
       {
	 TmpSum[k] += this->HamiltonianShift * Coefficient2[k];
	 vDestinations[k][l] += TmpSum[k];
       }
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (l = 0; l < this->FastMultiplicationStep; ++l)
    if (PosMod != l)
      {	
	for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
	  this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
	this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSources, vDestinations, nbrVectors);
      }
  delete[] TmpSum;
  delete[] Coefficient2;
  delete TmpParticles;
  return vDestinations;
}

// evaluate all interaction factors
//   

void AbstractQHEOnSphereFullHamiltonian::EvaluateInteractionFactors()
{
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long AbstractQHEOnSphereFullHamiltonian::FastMultiplicationMemory(long allowedMemory)
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  
  this->NbrInteractionPerComponent = new int [EffectiveHilbertSpaceDimension];
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    this->NbrInteractionPerComponent[i] = 0;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;

  QHEParticlePrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);

  long Memory = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    Memory += this->NbrInteractionPerComponent[i];

  cout << "nbr interaction = " << Memory << endl;
  long TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(double*)) * EffectiveHilbertSpaceDimension;
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(double)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(double)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if (this->Particles->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	    Memory += this->NbrInteractionPerComponent[i];
	}
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension) + (Memory * (sizeof (int) + sizeof(double)));
      long ResidualMemory = allowedMemory - Memory;
      if (ResidualMemory > 0)
	{
	  int TotalReducedSpaceDimension = ReducedSpaceDimension;
	  int* TmpNbrInteractionPerComponent = new int [TotalReducedSpaceDimension];
	  int i = 0;
	  int Pos = 0;
	  for (; i < ReducedSpaceDimension; ++i)
	    {
	      TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[Pos];
	      Pos += this->FastMultiplicationStep;
	    }
	  delete[] this->NbrInteractionPerComponent;
	      this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
	}
      else
	{
	  int* TmpNbrInteractionPerComponent = new int [ReducedSpaceDimension];
	  for (int i = 0; i < ReducedSpaceDimension; ++i)
	    TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[i * this->FastMultiplicationStep];
	  delete[] this->NbrInteractionPerComponent;
	  this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
	}
    }
  else
    {
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(double*)) * EffectiveHilbertSpaceDimension) + (Memory * (sizeof (int) + sizeof(double)));
      this->FastMultiplicationStep = 1;
    }

  cout << "reduction factor=" << this->FastMultiplicationStep << endl;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  return Memory;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long AbstractQHEOnSphereFullHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  long Memory = 0;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int LastComponent = lastComponent + firstComponent;
  this->EvaluateMNTwoBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);
  this->EvaluateMNOneBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);
  
  delete TmpParticles;

  return Memory;
}

// enable fast multiplication algorithm
//

void AbstractQHEOnSphereFullHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  long Pos;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficient = new double* [ReducedSpaceDimension];

  for (int i = 0; i < ReducedSpaceDimension; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [this->NbrInteractionPerComponent[i]];
      this->InteractionPerComponentCoefficient[i] = new double [this->NbrInteractionPerComponent[i]];
    }

  QHEParticlePrecalculationOperation Operation(this, false);
  Operation.ApplyOperation(this->Architecture);

  this->FastMultiplicationFlag = true;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = index of the last component that has to be precalcualted

void AbstractQHEOnSphereFullHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{  
  int LastComponent = nbrComponent + firstComponent;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();

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
      long TotalPos = 0;
      this->EvaluateMNTwoBodyFastMultiplicationComponent(TmpParticles, i, this->InteractionPerComponentIndex[Pos], 
							 this->InteractionPerComponentCoefficient[Pos], TotalPos);
      this->EvaluateMNOneBodyFastMultiplicationComponent(TmpParticles, i, this->InteractionPerComponentIndex[Pos], 
							 this->InteractionPerComponentCoefficient[Pos], TotalPos);
      ++Pos;
    }
  delete TmpParticles;
}

