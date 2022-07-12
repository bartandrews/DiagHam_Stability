////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                to particles on a torus with magnetic translations          //
//                                                                            //
//                        last modification : 18/02/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeWithSpinChernInsulatorHamiltonian.h"
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

ParticleOnLatticeWithSpinChernInsulatorHamiltonian::ParticleOnLatticeWithSpinChernInsulatorHamiltonian()
{
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->NbrIntraSectorSums = 0;
  this->NbrInterSectorSums = 0;
  this->HermitianSymmetryFlag = false;
  this->S2Hamiltonian = 0;
  this->LoadBalancingArray = 0;
  this->NbrInteractionPerComponent = 0;
  this->FastMultiplicationStep = 0;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kineticFactor = multiplicative factor in front of the kinetic term
// uPotential = Hubbard potential strength
// bandParameter = band parameter
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinChernInsulatorHamiltonian::ParticleOnLatticeWithSpinChernInsulatorHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, 
												       int nbrSiteY, double kineticFactor, double uPotential, double bandParameter, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->UPotential = 2.0 * uPotential / ((double) this->NbrParticles);
  this->HamiltonianShift = 0.0;//4.0 * uPotential;
  this->KineticFactor = kineticFactor;
  this->BandParameter = bandParameter;
  this->FlatBand = flatBandFlag;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->HermitianSymmetryFlag = false;
  this->EvaluateInteractionFactors();
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      if (TmpMemory < 1024)
	cout  << "fast = " <<  TmpMemory << "b ";
      else
	if (TmpMemory < (1 << 20))
	  cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	else
	  if (TmpMemory < (1 << 30))
	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	  else
	    {
	      cout  << "fast = " << (TmpMemory >> 30) << ".";
	      TmpMemory -= ((TmpMemory >> 30) << 30);
	      TmpMemory *= 100l;
	      TmpMemory >>= 30;
	      if (TmpMemory < 10l)
		cout << "0";
	      cout  << TmpMemory << " Gb ";
	    }
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeWithSpinChernInsulatorHamiltonian::~ParticleOnLatticeWithSpinChernInsulatorHamiltonian()
{
  if (this->InteractionFactorsupup != 0)
    {
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  delete[] this->InteractionFactorsupup[i];
	}
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  delete[] this->InteractionFactorsupdown[i];
	}
      delete[] this->InteractionFactorsupup;
      delete[] this->InteractionFactorsupdown;
    }
  if (this->InteractionFactorsdowndown != 0)
    {
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  delete[] this->InteractionFactorsdowndown[i];
	}
      delete[] this->InteractionFactorsdowndown;
    }
  if (this->OneBodyInteractionFactorsupup != 0)
    {
      delete[] this->OneBodyInteractionFactorsupup;
    }
  if (this->OneBodyInteractionFactorsdowndown != 0)
    {
      delete[] this->OneBodyInteractionFactorsdowndown;
    }
  if (this->OneBodyInteractionFactorsupdown != 0)
    {
      delete[] this->OneBodyInteractionFactorsupdown;
    }
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    {
      if (this->NbrIntraSectorIndicesPerSum[i] > 0)
	delete[] this->IntraSectorIndicesPerSum[i];
    }
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      if (this->NbrInterSectorIndicesPerSum[i] > 0)
	delete[] this->InterSectorIndicesPerSum[i];
    }
  if (this->NbrIntraSectorSums > 0)
    {
      delete[] this->NbrIntraSectorIndicesPerSum;
      delete[] this->IntraSectorIndicesPerSum;
    }
  if (this->NbrInterSectorSums > 0)
    {
      delete[] this->NbrInterSectorIndicesPerSum;
      delete[] this->InterSectorIndicesPerSum;
    }
  if (this->S2Hamiltonian != 0)
    {
      delete this->S2Hamiltonian;
    }
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
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
      this->FastMultiplicationFlag = false;
    }
}
  
// ask if Hamiltonian implements hermitian symmetry operations
//

bool ParticleOnLatticeWithSpinChernInsulatorHamiltonian::IsHermitian()
{
  return this->HermitianSymmetryFlag;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeWithSpinChernInsulatorHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (ParticleOnSphereWithSpin*) hilbertSpace;
  this->EvaluateInteractionFactors();
}
 
// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnLatticeWithSpinChernInsulatorHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeWithSpinChernInsulatorHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnLatticeWithSpinChernInsulatorHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}

// add an additional S^2 term to the Hamiltonian
//
// factor = factor in front of the S^2
// fixedSz = flag indicating whether Sz needs to be evaluated
// memory = amount of memory that can be used for S^2  precalculations

void ParticleOnLatticeWithSpinChernInsulatorHamiltonian::AddS2(double factor, bool fixedSz, long memory)
{
  cout << "warning  ParticleOnLatticeWithSpinChernInsulatorHamiltonian::AddS2 is not implemented" << endl;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnLatticeWithSpinChernInsulatorHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
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
	  Complex* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  Complex Coefficient;
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
	  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  Complex Coefficient;
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
  if (this->S2Hamiltonian != 0)
    {
      this->S2Hamiltonian->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
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

ComplexVector* ParticleOnLatticeWithSpinChernInsulatorHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
											       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      Complex* Coefficient2 = new Complex [nbrVectors];
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
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
	  Complex Coefficient;
	  Complex* Coefficient2 = new Complex [nbrVectors];
	  Complex* TmpCoefficientArray; 
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
	  firstComponent += this->PrecalculationShift;
	  LastComponent += this->PrecalculationShift;
	}
      else
	{
	  this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	}
    }
  if (this->S2Hamiltonian != 0)
    {
      this->S2Hamiltonian->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
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

ComplexVector* ParticleOnLatticeWithSpinChernInsulatorHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
												   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int* TmpIndexArray;
  Complex* TmpCoefficientArray; 
  int j;
  int TmpNbrInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int Pos2;
  Complex Coefficient;
  int PosMod = firstComponent % this->FastMultiplicationStep;
  Complex* Coefficient2 = new Complex [nbrVectors];
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
  if (this->S2Hamiltonian != 0)
    {
      this->S2Hamiltonian->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources,  vDestinations, nbrVectors,firstComponent, nbrComponent);
    }
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

ComplexVector& ParticleOnLatticeWithSpinChernInsulatorHamiltonian::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
												int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
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
	  Complex* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  Complex Coefficient;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      Coefficient = vSource[k];
	      Complex TmpSum = 0.0;
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  TmpSum += Conj(TmpCoefficientArray[j]) * vSource[TmpIndexArray[j]];
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
		}
	      TmpSum += this->HamiltonianShift * Coefficient;
	      vDestination[k++] += TmpSum;
	    }
	  firstComponent += this->PrecalculationShift;
	  LastComponent += this->PrecalculationShift;
	}
      else
	{
	  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  Complex Coefficient;
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
	      Complex TmpSum = 0.0;
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
		  TmpSum += Conj(TmpCoefficientArray[j]) * vSource[TmpIndexArray[j]];
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
  if (this->S2Hamiltonian != 0)
    {
      this->S2Hamiltonian->HermitianLowLevelAddMultiply(vSource,  vDestination, firstComponent, nbrComponent);
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

ComplexVector* ParticleOnLatticeWithSpinChernInsulatorHamiltonian::HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
													int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      Complex* Coefficient2 = new Complex [nbrVectors];
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
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
	  Complex Coefficient;
	  Complex* Coefficient2 = new Complex [nbrVectors];
	  Complex* TmpCoefficientArray; 
	  int j;
	  int Pos;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  Complex* TmpSum = new Complex [nbrVectors];
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
		      TmpSum[l] += Conj(Coefficient) * vSources[l][Pos];
		    }
		}
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  TmpSum[l] += this->HamiltonianShift * Coefficient2[l];
		  vDestinations[l][k] += TmpSum[l];
		}
	      ++k;
	    }
	  firstComponent += this->PrecalculationShift;
	  LastComponent += this->PrecalculationShift;
	  delete[] Coefficient2;
	  delete[] TmpSum;
	}
      else
	{
	  this->HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	}
    }
  if (this->S2Hamiltonian != 0)
    {
      this->S2Hamiltonian->HermitianLowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
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

ComplexVector* ParticleOnLatticeWithSpinChernInsulatorHamiltonian::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
																     int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int* TmpIndexArray;
  Complex* TmpCoefficientArray; 
  int j;
  int TmpNbrInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int Pos2;
  Complex Coefficient;
  int PosMod = firstComponent % this->FastMultiplicationStep;
  Complex* Coefficient2 = new Complex [nbrVectors];
  Complex* TmpSum = new Complex [nbrVectors];
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
	      TmpSum[k] += Conj(Coefficient) * vSources[k][Pos2];
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

void ParticleOnLatticeWithSpinChernInsulatorHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  this->NbrInterSectorSums = this->NbrSiteX * this->NbrSiteY;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)      
	  ++this->NbrInterSectorIndicesPerSum[(((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)];    
  this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      if (this->NbrInterSectorIndicesPerSum[i] > 0)
	{
	  this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
	  this->NbrInterSectorIndicesPerSum[i] = 0;
	}
    }
  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)    
	  {
	    int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = (kx1 * this->NbrSiteY) + ky1;
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = (kx2 * this->NbrSiteY) + ky2;
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	  }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->NbrIntraSectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = (kx1 * this->NbrSiteY) + ky1;
		int Index2 = (kx2 * this->NbrSiteY) + ky2;
		if (Index1 < Index2)
		  ++this->NbrIntraSectorIndicesPerSum[(((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)];    
	      }
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
	    {
	      this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	      this->NbrIntraSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
	      {
		int Index1 = (kx1 * this->NbrSiteY) + ky1;
		int Index2 = (kx2 * this->NbrSiteY) + ky2;
		if (Index1 < Index2)
		  {
		    int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int kx1 = this->IntraSectorIndicesPerSum[i][j1 << 1] / this->NbrSiteY;
	      int ky1 = this->IntraSectorIndicesPerSum[i][j1 << 1] % this->NbrSiteY;
	      int kx2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1] / this->NbrSiteY;
	      int ky2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1] % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int kx3 = this->IntraSectorIndicesPerSum[i][j2 << 1] / this->NbrSiteY;
		  int ky3 = this->IntraSectorIndicesPerSum[i][j2 << 1] % this->NbrSiteY;
		  int kx4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1] / this->NbrSiteY;
		  int ky4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1] % this->NbrSiteY;
		  this->InteractionFactorsupup[i][Index] = - this->UPotential * ((cos (2.0 * M_PI * ((double) kx2 - kx4) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky2 - ky4) / ((double) this->NbrSiteY)))
										 - (cos (2.0 * M_PI * ((double) kx1 - kx4) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky1 - ky4) / ((double) this->NbrSiteY)))
										 - (cos (2.0 * M_PI * ((double) kx2 - kx3) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky2 - ky3) / ((double) this->NbrSiteY)))
										 + (cos (2.0 * M_PI * ((double) kx1 - kx3) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky1 - ky3) / ((double) this->NbrSiteY))));
		  this->InteractionFactorsdowndown[i][Index] = this->InteractionFactorsupup[i][Index];
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}
      this->InteractionFactorsupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int kx1 = this->InterSectorIndicesPerSum[i][j1 << 1] / this->NbrSiteY;
	      int ky1 = this->InterSectorIndicesPerSum[i][j1 << 1] % this->NbrSiteY;
	      int kx2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1] / this->NbrSiteY;
	      int ky2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1] % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int kx3 = this->InterSectorIndicesPerSum[i][j2 << 1] / this->NbrSiteY;
		  int ky3 = this->InterSectorIndicesPerSum[i][j2 << 1] % this->NbrSiteY;
		  int kx4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1] / this->NbrSiteY;
		  int ky4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1] % this->NbrSiteY;
		  this->InteractionFactorsupdown[i][Index] = - this->UPotential * (0.5 + (cos (2.0 * M_PI * ((double) kx2 - kx4) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky2 - ky4) / ((double) this->NbrSiteY)))
										   - (cos (2.0 * M_PI * ((double) kx1 - kx4) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky1 - ky4) / ((double) this->NbrSiteY)))
										   - (cos (2.0 * M_PI * ((double) kx2 - kx3) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky2 - ky3) / ((double) this->NbrSiteY)))
										   + (cos (2.0 * M_PI * ((double) kx1 - kx3) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky1 - ky3) / ((double) this->NbrSiteY))));
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  this->OneBodyInteractionFactorsupup = new double [this->LzMax + 1];
  this->OneBodyInteractionFactorsdowndown = new double [this->LzMax + 1];
  this->OneBodyInteractionFactorsupdown = new Complex [this->LzMax + 1];
  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
      {
	int Index = (kx1 * this->NbrSiteY) + ky1;
	double d1 = sin (2.0 * M_PI * ((double) kx1) / ((double) this->NbrSiteX));
	double d2 = sin (2.0 * M_PI * ((double) ky1) / ((double) this->NbrSiteY));
	double d3 = (this->BandParameter - cos (2.0 * M_PI * ((double) ky1) / ((double) this->NbrSiteY))
		     - cos (2.0 * M_PI * ((double) kx1) / ((double) this->NbrSiteX)));
	HermitianMatrix TmpOneBobyHamiltonian(2, true);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 0, d3);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 1, Complex(d1, -d2));
	TmpOneBobyHamiltonian.SetMatrixElement(1, 1, -d3);
	ComplexMatrix TmpMatrix(2, 2, true);
	TmpMatrix[0][0] = 1.0;
	TmpMatrix[1][1] = 1.0;
	RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	TmpOneBobyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBobyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << endl;
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactorsupup[Index] = this->KineticFactor * ((TmpDiag(0, 0) * SqrNorm(TmpMatrix[0][0])) 
										+ (TmpDiag(1, 1) * SqrNorm(TmpMatrix[1][0])));
	    this->OneBodyInteractionFactorsdowndown[Index] = this->KineticFactor * ((TmpDiag(0, 0) * SqrNorm(TmpMatrix[0][1])) 
										    + (TmpDiag(1, 1) * SqrNorm(TmpMatrix[1][1])));
	    this->OneBodyInteractionFactorsupdown[Index] = this->KineticFactor * ((TmpDiag(0, 0) * TmpMatrix[0][0] * Conj(TmpMatrix[0][1])) 
										  + (TmpDiag(1, 1) * TmpMatrix[1][0] * Conj(TmpMatrix[1][1])));
	  }
	else
	  {
	    this->OneBodyInteractionFactorsupup[Index] = this->KineticFactor * (SqrNorm(TmpMatrix[1][0]) - SqrNorm(TmpMatrix[0][0]));
	    this->OneBodyInteractionFactorsdowndown[Index] = this->KineticFactor * (SqrNorm(TmpMatrix[1][1]) - SqrNorm(TmpMatrix[0][1]));
	    this->OneBodyInteractionFactorsupdown[Index] = this->KineticFactor * ((TmpMatrix[1][0] * Conj(TmpMatrix[1][1])) - (TmpMatrix[0][0] * Conj(TmpMatrix[0][1])));
	  }
      }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long ParticleOnLatticeWithSpinChernInsulatorHamiltonian::FastMultiplicationMemory(long allowedMemory)
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

  if (this->Architecture->GetOptimizedTypicalRange(this->NbrInteractionPerComponent, MinIndex, MaxIndex) == true)
    {
      this->PrecalculationShift = (int) MinIndex;
      EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
      cout << "distributed calculations have been reoptimized" << endl;
    }  
  if (allowedMemory == 0l)
    {
      delete[] this->NbrInteractionPerComponent;
      this->NbrInteractionPerComponent = 0;
      return 0l;
    }

  long Memory = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    Memory += this->NbrInteractionPerComponent[i];

  cout << "nbr interaction = " << Memory << endl;
  long TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(Complex*)) * EffectiveHilbertSpaceDimension;
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(Complex)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(Complex)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if (this->Particles->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(Complex*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	    Memory += this->NbrInteractionPerComponent[i];
	}
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(Complex*)) * ReducedSpaceDimension) + (Memory * (sizeof (int) + sizeof(Complex)));
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
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(Complex*)) * EffectiveHilbertSpaceDimension) + (Memory * (sizeof (int) + sizeof(Complex)));
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

long ParticleOnLatticeWithSpinChernInsulatorHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  long Memory = 0;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int LastComponent = lastComponent + firstComponent;

  this->EvaluateMNOneBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);
  this->EvaluateMNTwoBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);

  delete TmpParticles;

  return Memory;
}

// enable fast multiplication algorithm
//

void ParticleOnLatticeWithSpinChernInsulatorHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficient = new Complex* [ReducedSpaceDimension];
  //ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles;
  int TotalPos = 0;
  

  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
    {
      this->InteractionPerComponentIndex[TotalPos] = new int [this->NbrInteractionPerComponent[TotalPos]];
      this->InteractionPerComponentCoefficient[TotalPos] = new Complex [this->NbrInteractionPerComponent[TotalPos]];      
//       TmpIndexArray = this->InteractionPerComponentIndex[TotalPos];
//       TmpCoefficientArray = this->InteractionPerComponentCoefficient[TotalPos];
//       Pos = 0l;
//       this->EvaluateMNTwoBodyFastMultiplicationComponent(TmpParticles, i, TmpIndexArray, TmpCoefficientArray, Pos);
//       this->EvaluateMNOneBodyFastMultiplicationComponent(TmpParticles, i, TmpIndexArray, TmpCoefficientArray, Pos);
       ++TotalPos;
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

void ParticleOnLatticeWithSpinChernInsulatorHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{  
  int LastComponent = nbrComponent + firstComponent;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();

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
      this->EvaluateMNTwoBodyFastMultiplicationComponent(TmpParticles, i + this->PrecalculationShift, this->InteractionPerComponentIndex[Pos], 
							 this->InteractionPerComponentCoefficient[Pos], TotalPos);
      this->EvaluateMNOneBodyFastMultiplicationComponent(TmpParticles, i + this->PrecalculationShift, this->InteractionPerComponentIndex[Pos], 
      							 this->InteractionPerComponentCoefficient[Pos], TotalPos);
      
//       for (int j = 0; j < this->NbrInteractionPerComponent[Pos]; ++j)
// 	cout << Pos << " " << j << "  " << this->InteractionPerComponentCoefficient[Pos][j] << endl;
      ++Pos;
    }
  delete TmpParticles;
}


// get the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
//
// nbrThreads = number of threads requested
// segmentIndices = array returning the reference to an array of the first index of each of the segments
// return value = true if no error occured

bool ParticleOnLatticeWithSpinChernInsulatorHamiltonian::GetLoadBalancing(int nbrTasks, long* &segmentIndices)
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  if ((this->NbrInteractionPerComponent != 0) && (this->FastMultiplicationStep != 0))
    {
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;

      if ((this->LoadBalancingArray == 0) || (this->NbrBalancedTasks != nbrTasks))
	{
	  if (this->LoadBalancingArray != 0)
	    delete [] this->LoadBalancingArray;
	  long* SegmentSize = new long[nbrTasks];
	  this->LoadBalancingArray = new long[nbrTasks+1];
	  this->NbrBalancedTasks = nbrTasks;
	  long TmpNbrElement = 0l;
	  for (int i = 0; i < ReducedSpaceDimension; ++i)
	    TmpNbrElement += this->NbrInteractionPerComponent[i];
	  long TmpNbrPerSegment = TmpNbrElement / ((long) nbrTasks);
	  TmpNbrElement = 0l;
	  int Pos = 0;
	  this->LoadBalancingArray[0] = MinIndex;
	  for (int i = 0; i < ReducedSpaceDimension; ++i)
	    {
	      TmpNbrElement += this->NbrInteractionPerComponent[i];
	      if (TmpNbrElement > TmpNbrPerSegment)
		{
		  SegmentSize[Pos] = TmpNbrElement;
		  this->LoadBalancingArray[Pos + 1] = MinIndex + (i * this->FastMultiplicationStep);
		  TmpNbrElement = 0l;
		  ++Pos;
		}
	    }
	  while (Pos < (nbrTasks - 1))
	    {
	      LoadBalancingArray[Pos + 1] = MaxIndex + 1;
	      SegmentSize[Pos] = 0;
	      ++Pos;
	    }
	  this->LoadBalancingArray[nbrTasks] = MaxIndex + 1;
	  SegmentSize[nbrTasks - 1] = TmpNbrElement;
	  
	  cout << "LoadBalancingArray=[ (" << (this->LoadBalancingArray[1] - this->LoadBalancingArray[0]) << ", " << SegmentSize[0] <<")";
	  for (int i = 1; i < nbrTasks; ++i)
	    cout << " (" << (this->LoadBalancingArray[i + 1] - this->LoadBalancingArray[i]) << ", " << SegmentSize[i] << ")";
	  cout << "]"<< endl;
	  delete[] SegmentSize;
	}
    }
  else
    {
      return false;
    }
  segmentIndices = this->LoadBalancingArray;
  return true;
}

