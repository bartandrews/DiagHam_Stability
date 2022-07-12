////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                   class of hamiltonian with particles on                   //
//               Chern insulator in the single band approximation             //
//                           and n-body interaction                           //
//                                                                            //
//                        last modification : 03/08/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeWithKyNBodyDeltaHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MathTools/IntegerAlgebraTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;



// default constructor
//

ParticleOnLatticeWithKyNBodyDeltaHamiltonian::ParticleOnLatticeWithKyNBodyDeltaHamiltonian()
{
  this->TwoBodyFlag = false;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lx = length of simulation cell in x-direction
// ly = length of simulation cell in y-direction
// kyMax = maximum value of momentum in y-direction
// nbrFluxQuanta = number of flux quanta piercing the simulation cell
// contactInteractionU = strength of on-site delta interaction
// reverseHopping = flag to indicate if sign of hopping terms should be reversed
// randomPotential = strength of a random on-site potential
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnLatticeWithKyNBodyDeltaHamiltonian::ParticleOnLatticeWithKyNBodyDeltaHamiltonian(ParticleOnLattice* particles, int nbrParticles, int lx, int ly, int kyMax, int nbrFluxQuanta, int nbrBody , double contactInteractionU, double nBodyContactInteraction, bool reverseHopping, double randomPotential, AbstractArchitecture* architecture, unsigned long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->Lx = lx;
  this->Ly = ly;
  this->KyMax = kyMax;  
  this->NbrCells = lx*ly;
  this->NbrSublattices = 1;
  this->NbrSites = NbrCells*this->NbrSublattices;
  this->HamiltonianShift = 0.0;//4.0 * uPotential;
  this->NBodyValue = nbrBody;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->FluxDensity = ((double)nbrFluxQuanta)/NbrCells;
  this->TwoBodyFlag = false;
  this->ContactInteractionU = contactInteractionU;
  if(this->ContactInteractionU != 0.0)
    this->TwoBodyFlag = true;
  this->NBodyContactInteraction = nBodyContactInteraction;
  this->ReverseHopping = reverseHopping;
  this->DeltaPotential = 0.0;
  this->RandomPotential = randomPotential;
  this->Architecture = architecture;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  this->HaveKySymmetry = true;
  this->DiskStorageFlag = false;
  if(this->TwoBodyFlag == false)
    this->NbrInteractionFactors = 0;
  if (precalculationFileName == 0)
    {
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
		cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
	  cout << endl;
	  if (memory > 0)
	    {
	      this->EnableFastMultiplication();
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnLatticeWithKyNBodyDeltaHamiltonian::~ParticleOnLatticeWithKyNBodyDeltaHamiltonian()
{
	
	if ( this->NBodyContactInteraction != 0.0 )
	{
		for (int i = 0; i < this->NbrNBodySectorSums; ++i)
		{
			if (this->NbrNBodySectorIndicesPerSum[i] > 0)
	{
	  delete[] this->NBodySectorIndicesPerSum[i];
	  delete[] this->NBodyInteractionFactors[i];
	}
    }

      delete[] this->NBodyInteractionFactors;
      delete[] this->NbrNBodySectorIndicesPerSum;
      delete[] this->NBodySectorIndicesPerSum;
    }
}
  
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnLatticeWithKyNBodyDeltaHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;

  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      if (this->TwoBodyFlag == true)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);	  
	      this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	    }
	    this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSource, vDestination);
	}
      else
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);	  
	    }
	  this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSource, vDestination);
	}
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->LowLevelAddMultiplyFastMultiply(vSource, vDestination, firstComponent, LastComponent);
	}
      else
	{
	   this->LowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
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

ComplexVector* ParticleOnLatticeWithKyNBodyDeltaHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors,  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      Complex* Coefficient2 = new Complex [nbrVectors];
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
			if (this->TwoBodyFlag == true)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
	      this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
	    }
	}
      else
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
	    }
	}
	this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSources, vDestinations, nbrVectors);
      
      delete[] Coefficient2;
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->LowLevelMultipleAddMultiplyFastMultiply(vSources, vDestinations, nbrVectors , firstComponent , LastComponent);
	}
      else
	{
	  this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	}
    }
  return vDestinations;
}

ComplexVector& ParticleOnLatticeWithKyNBodyDeltaHamiltonian::ConjugateLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination,  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      this->EvaluateMNOneBodyConjugateAddMultiplyComponent(TmpParticles,firstComponent,LastComponent,1, vSource, vDestination);
      if(this->TwoBodyFlag == true)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      this->EvaluateMNTwoBodyConjugateAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	      this->EvaluateMNNBodyConjugateAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	    }
	  
	}
      else
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      this->EvaluateMNNBodyConjugateAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	    }
	}
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->ConjugateLowLevelAddMultiplyFastMultiply(vSource,  vDestination, firstComponent, LastComponent);
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->ConjugateLowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
  return vDestination;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnLatticeWithKyNBodyDeltaHamiltonian::ConjugateLowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
	
	cout <<"calling non degined function ParticleOnLatticeWithKyNBodyDeltaHamiltonian::ConjugateLowLevelAddMultiplyPartialFastMultiply"<<endl;
	exit(1);
	return vDestination;
 /* int Index;
  Complex TmpInteraction;
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  Complex TmpSum=0.0;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  int* TmpIndexArray;
  unsigned short* TmpCoefficientIndexArray;
  Complex TmpC;
  int TmpNbrRealInteraction;
  int TmpNbrComplexInteraction;
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
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[Pos];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[Pos];
      TmpSum=0.0;
      int Pos2=0;
      for (; Pos2 < TmpNbrRealInteraction; ++Pos2)
	TmpSum +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*vSource[TmpIndexArray[Pos2]];
      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos2)
	TmpSum +=  Conj(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos2]])*vSource[TmpIndexArray[Pos2]];
      vDestination[l] += this->HamiltonianShift * vSource[l];
      l += this->FastMultiplicationStep;
      ++Pos;
    }

  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	this->EvaluateMNOneBodyConjugateAddMultiplyComponent(TmpParticles, firstComponent + k, LastComponent, this->FastMultiplicationStep, vSource,vDestination);
	
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	  {
	    this->EvaluateMNTwoBodyConjugateAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	  }
      }
  delete TmpParticles;
  return vDestination;*/
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

ComplexVector* ParticleOnLatticeWithKyNBodyDeltaHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function ParticleOnLatticeWithKyNBodyDeltaHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply"<<endl;
  return vDestinations;
  
  /*	
  int LastComponent = firstComponent + nbrComponent;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
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
	if (this->TwoBodyFlag == true)
	  {
	    for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
		this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
		}
	      }
	      else
	      {
	    for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
	      }
	  }
	  this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSources, vDestinations, nbrVectors);
			}
  delete[] Coefficient2;
  delete TmpParticles;
  return vDestinations;*/
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnLatticeWithKyNBodyDeltaHamiltonian::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination,int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      if (this->TwoBodyFlag == true)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      this->HermitianEvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	      this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	    }
	}
      else
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      this->HermitianEvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	    }
	}
      this->HermitianEvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSource, vDestination);
      
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->HermitianLowLevelAddMultiplyFastMultiply(vSource, vDestination, firstComponent, LastComponent);
	}
      else
	{
	  this->HermitianLowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
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

ComplexVector* ParticleOnLatticeWithKyNBodyDeltaHamiltonian::HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
												  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      Complex* Coefficient2 = new Complex [nbrVectors];
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      if (this->TwoBodyFlag == true)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      this->HermitianEvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
	      this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
	    }
	  this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSources, vDestinations, nbrVectors);
	}
      else
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      this->HermitianEvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
	    }
	  this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSources, vDestinations, nbrVectors);
	}
      delete[] Coefficient2;
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->HermitianLowLevelMultipleAddMultiplyFastMultiply(vSources, vDestinations, nbrVectors ,firstComponent, LastComponent);
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

ComplexVector* ParticleOnLatticeWithKyNBodyDeltaHamiltonian::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors,  int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function ParticleOnLatticeWithKyNBodyDeltaHamiltonian::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply"<<endl;
  return vDestinations;
//   int LastComponent = firstComponent + nbrComponent;
//   ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
//   int* TmpIndexArray;
//   Complex* TmpCoefficientArray; 
//   int j;
//   int TmpNbrInteraction;
//   firstComponent -= this->PrecalculationShift;
//   LastComponent -= this->PrecalculationShift;
//   int Pos = firstComponent / this->FastMultiplicationStep; 
//   int Pos2;
//   Complex Coefficient;
//   int PosMod = firstComponent % this->FastMultiplicationStep;
//   Complex* Coefficient2 = new Complex [nbrVectors];
//   Complex* TmpSum = new Complex [nbrVectors];
//   if (PosMod != 0)
//     {
//       ++Pos;
//       PosMod = this->FastMultiplicationStep - PosMod;
//     }
//   int l =  PosMod + firstComponent + this->PrecalculationShift;
//   for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
//     {
//       TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
//       TmpIndexArray = this->InteractionPerComponentIndex[Pos];
//       TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
//       for (int k = 0; k < nbrVectors; ++k)
// 	{
// 	  TmpSum[k] = 0.0;
// 	  Coefficient2[k] = vSources[k][l];
// 	}
//       for (j = 0; j < TmpNbrInteraction; ++j)
// 	{
// 	  Pos2 = TmpIndexArray[j];
// 	  Coefficient = TmpCoefficientArray[j];
// 	  for (int k = 0; k < nbrVectors; ++k)
// 	    {
// 	      vDestinations[k][TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient2[k];
// 	      TmpSum[k] += Conj(Coefficient) * vSources[k][Pos2];
// 	    }
// 	}
//      for (int k = 0; k < nbrVectors; ++k)
//        {
// 	 TmpSum[k] += this->HamiltonianShift * Coefficient2[k];
// 	 vDestinations[k][l] += TmpSum[k];
//        }
//       l += this->FastMultiplicationStep;
//       ++Pos;
//     }
//   firstComponent += this->PrecalculationShift;
//   LastComponent += this->PrecalculationShift;
//   for (l = 0; l < this->FastMultiplicationStep; ++l)
//     if (PosMod != l)
//       {	
// 	if (this->TwoBodyFlag == true)
// 	  {
// 	    for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
// 	      {
// 		this->HermitianEvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
// 		this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
// 	      }
// 	  }
// 	else
// 	  {
// 	    for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
// 	      {
// 		this->HermitianEvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
// 	      }
// 	  }
// 	  this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSources, vDestinations, nbrVectors);
//       }
//   delete[] TmpSum;
//   delete[] Coefficient2;
//   delete TmpParticles;
//   return vDestinations;
}

// evaluate all interaction factors
//   

void ParticleOnLatticeWithKyNBodyDeltaHamiltonian::EvaluateInteractionFactors()
{
  int MaxNumberTerms=4*this->NbrSites;
  this->NbrDiagonalInteractionFactors = 0;
  this->DiagonalInteractionFactors = 0;
  this->DiagonalQValues = 0;
  if (this->DeltaPotential != 0.0) ++MaxNumberTerms;
  if (this->RandomPotential != 0.0) MaxNumberTerms+=this->NbrSites;  
  this->HoppingTerms = new Complex[MaxNumberTerms];
  this->KineticQi = new int[MaxNumberTerms];
  this->KineticQf = new int[MaxNumberTerms];
  
  this->NbrHoppingTerms = 0;
  
  double HoppingSign = (this->ReverseHopping ? 1.0 : -1.0);
  Complex TranslationPhase;
  int p = Ly/KyMax;
  int compositeKy, compositeKy2;
  for (int i = 0; i<Lx; ++i) 
    {
      Complex Phase=Polar(1.0,2.0*M_PI*this->FluxDensity*(double)i);
      for (int k=0; k<KyMax; ++k)
	for (int s=0; s<p; ++s)
	  {
	    compositeKy = k*p+s;
	    KineticQi[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy, 0, TranslationPhase);
	    KineticQf[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i+1, compositeKy, 0, TranslationPhase);
	    HoppingTerms[this->NbrHoppingTerms] = HoppingSign*TranslationPhase;	    
	    
#ifdef DEBUG_OUTPUT
	    if (TranslationPhase!=1.0)
	      cout << "(i="<<i<<"->"<<i+1<<") Translation ["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="<<TranslationPhase<<endl;
	    cout << "k="<<k<<", "<<"s="<<s<<" x_0="<<i<<endl;
	    cout << "x: H["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="<<HoppingTerms[this->NbrHoppingTerms]<<" tP="<<TranslationPhase<<endl;
#endif
	    ++this->NbrHoppingTerms;
	    
	    KineticQi[this->NbrHoppingTerms] = KineticQi[this->NbrHoppingTerms-1];
	    KineticQf[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i-1, compositeKy, 0, TranslationPhase);
	    HoppingTerms[this->NbrHoppingTerms] = HoppingSign*TranslationPhase;
	    
#ifdef DEBUG_OUTPUT
	    if (TranslationPhase!=1.0)
	      cout << "(i="<<i<<"->"<<i-1<<") Translation ["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="
		   <<TranslationPhase<<endl;	    
	    cout << "x: H["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="<<HoppingTerms[this->NbrHoppingTerms]<<" tP="<<TranslationPhase<<endl;
#endif
	    ++this->NbrHoppingTerms;
	    
	    // coupling in y-direction: distinguish p=0 and p>0!
	    if (p==1)
	      {
		KineticQi[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, k, 0, TranslationPhase);
		KineticQf[this->NbrHoppingTerms] = KineticQi[this->NbrHoppingTerms];
		HoppingTerms[this->NbrHoppingTerms] = HoppingSign*2.0*cos(2.0*M_PI*((double)k/Ly-this->FluxDensity*(double)i));
#ifdef DEBUG_OUTPUT
		cout << "y - p=1: H["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="<<HoppingTerms[this->NbrHoppingTerms]<<endl;
#endif
		++this->NbrHoppingTerms;
	      }
	    else
	      {
		if (s<p-1)
		  {
		    KineticQi[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy, 0, TranslationPhase);
		    compositeKy2 = k*p+s+1;
		    KineticQf[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy2, 0, TranslationPhase);
		    HoppingTerms[this->NbrHoppingTerms] = HoppingSign*Conj(Phase)*TranslationPhase;
		  }
		else		  
		  {
		    KineticQi[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy, 0, TranslationPhase);
		    compositeKy2 = k*p;
		    KineticQf[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy2, 0, TranslationPhase);
		    //HoppingTerms[this->NbrHoppingTerms] = HoppingSign*Phase*TranslationPhase
				//*Polar(1.0,-2.0*M_PI*((double)k/KyMax-this->FluxDensity*(double)i));
				HoppingTerms[this->NbrHoppingTerms] = HoppingSign*Conj(Phase)*TranslationPhase*Polar(1.0,-2.0*M_PI*((double)k/KyMax));
		  }		
#ifdef DEBUG_OUTPUT
		if (TranslationPhase!=1.0)
		  cout << "(sk)="<<compositeKy<<"->"<<compositeKy2<<") Translation ["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="
				 <<TranslationPhase<<endl;
		cout << "y - p>1: H["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="<<HoppingTerms[this->NbrHoppingTerms]<<" tP="<<TranslationPhase<<endl;
#endif
		++this->NbrHoppingTerms;
		if (s==0)
		  {
		    KineticQi[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy, 0, TranslationPhase);
		    compositeKy2 = k*p+p-1;
		    KineticQf[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy2, 0, TranslationPhase);
		    //HoppingTerms[this->NbrHoppingTerms] = HoppingSign*Conj(Phase)*TranslationPhase
		      //*Polar(1.0,2.0*M_PI*((double)k/KyMax-this->FluxDensity*(double)i));
				HoppingTerms[this->NbrHoppingTerms] = HoppingSign*Phase*TranslationPhase*Polar(1.0,2.0*M_PI*((double)k/KyMax));
		  }
		else
		  {
		    KineticQi[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy, 0, TranslationPhase);
		    compositeKy2 = k*p+s-1;
		    KineticQf[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy2, 0, TranslationPhase);
		    HoppingTerms[this->NbrHoppingTerms] = HoppingSign*Phase*TranslationPhase;
		  }
#ifdef DEBUG_OUTPUT
		if (TranslationPhase!=1.0)
		  cout << "(sk)="<<compositeKy<<"->"<<compositeKy2<<") Translation ["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="
				 <<TranslationPhase<<endl;
		cout << "y - p>1: H["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="<<HoppingTerms[this->NbrHoppingTerms]<<" tP="<<TranslationPhase<<endl;
#endif
		++this->NbrHoppingTerms;
	      }
	  }
    }
  
  cout << "nbr interaction One body= " << this->NbrHoppingTerms << endl;
  cout << "====================================" << endl;
  
  long TotalNbrInteractionFactors = 0;
  
  this->NbrNBodySectorSums = 0;

  if ((this->Particles->GetParticleStatistic() == ParticleOnLattice::BosonicStatistic) && (this->NBodyContactInteraction != 0.0))
    {
      this->NbrNBodySectorSums = this->Lx * this->Ly;
      this->NbrNBodySectorIndicesPerSum = new int[this->NbrNBodySectorSums];
      
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
				this->NbrNBodySectorIndicesPerSum[i] = 0;
      
      double Strength = 1;
      for (int i = 2 ; i<= this->NBodyValue ; i++)
				Strength *= i*i;
      
      for (int i = 2 ; i<= this->NBodyValue ; i++)
				Strength /= ((double)KyMax);
      
      Strength *= this->NBodyContactInteraction;
      
      cout <<"Strength "<< Strength<<endl;
      double ** TmpSortedIndicesPerSumSymmetryFactor;
      int * TmpNbrBodySectorIndicesPerSum;
      int ** TmpNBodySectorIndicesPerSum;
      
      int TmpNbrSum = this->Ly;
      GetAllSymmetricIndices (this->KyMax, this->NBodyValue, TmpNbrBodySectorIndicesPerSum, TmpNBodySectorIndicesPerSum, TmpSortedIndicesPerSumSymmetryFactor);
      
			for (int i=0; i<Lx; ++i)
				for (int Sum = 0; Sum < TmpNbrSum; Sum++)
				{
					this->NbrNBodySectorIndicesPerSum[i + this->Lx * Sum] += TmpNbrBodySectorIndicesPerSum[Sum];
				}
      
      this->NBodySectorIndicesPerSum = new int* [this->NbrNBodySectorSums];
      
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
			{
				if (this->NbrNBodySectorIndicesPerSum[i]  > 0)
				{
					this->NBodySectorIndicesPerSum[i] = new int[this->NBodyValue * this->NbrNBodySectorIndicesPerSum[i]];      
					this->NbrNBodySectorIndicesPerSum[i] = 0;
				}
			}
      
      for (int i=0; i<Lx; ++i)
	for (int Sum = 0; Sum < TmpNbrSum; Sum++)
	  {
	    int TmpSum = i + this->Lx * Sum;
	    for (int IndexinSum = 0; IndexinSum < TmpNbrBodySectorIndicesPerSum[Sum] ; IndexinSum ++)
	      {
		for(int TmpPos = 0; TmpPos< this->NBodyValue; TmpPos++)
		  {
		    this->NBodySectorIndicesPerSum[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum] * this->NBodyValue + TmpPos] = i + this->Lx * TmpNBodySectorIndicesPerSum[Sum][IndexinSum * this->NBodyValue +TmpPos];
		  }
		//SortedIndicesPerSumSymmetryFactor[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum]] = TmpSortedIndicesPerSumSymmetryFactor[Sum][IndexinSum];
		++this->NbrNBodySectorIndicesPerSum[TmpSum];
	      }
	    
	  }
	  
	  	/*for (int Sum = 0; Sum < this->NbrNBodySectorSums; Sum++)
		  {
		    cout <<"Sum = " << Sum <<endl;
		    for (int i = 0; i <  this->NbrNBodySectorIndicesPerSum[Sum]; i++)
				{
					for(int TmpPos = 0; TmpPos< this->NBodyValue; TmpPos++)
					{
						cout <<this->NBodySectorIndicesPerSum[Sum][i * this->NBodyValue +TmpPos]<<" ";
					}
					cout <<endl<<SortedIndicesPerSumSymmetryFactor[Sum][i]<<endl<<endl;
				}
		  }*/
      
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
		//cout <<"Sum = "<<i <<endl;
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int k = i / this->Lx;
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1) // annihilation operators
	    { 
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		// creation operators
		{
			//cout << "j1  = " << j1<<"; j2 = "<<j2<<"; Index = "<<Index<<endl;
		  this->NBodyInteractionFactors[i][Index] =  Strength * TmpSortedIndicesPerSumSymmetryFactor[k][j1] * TmpSortedIndicesPerSumSymmetryFactor[k][j2];
			//this->NBodyInteractionFactors[i][Index] =  Strength * SortedIndicesPerSumSymmetryFactor[i][j1] * SortedIndicesPerSumSymmetryFactor[i][j2];
			//cout <<this->NBodyInteractionFactors[i][Index];
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
      
      for (int Sum = 0; Sum < TmpNbrSum; Sum++)
	{
	  delete [] TmpSortedIndicesPerSumSymmetryFactor[Sum];
	}
      
      delete [] TmpSortedIndicesPerSumSymmetryFactor;
      delete [] TmpNbrBodySectorIndicesPerSum;
      delete [] TmpNBodySectorIndicesPerSum;
    }
  else
    {
      TotalNbrInteractionFactors = 0l;
    }
    if ((this->Particles->GetParticleStatistic() == ParticleOnLattice::BosonicStatistic) && (this->TwoBodyFlag))
    {
      // general four-particle interactions:
      int Pos=0;
      int k4;
      double* TmpCoefficient = new double [Lx * p * KyMax * KyMax * KyMax];
      double Strength = this->ContactInteractionU / KyMax;  // Check Amplitude!!
      double MaxCoefficient = 0.0;
      for (int i=0; i<Lx; ++i)
	for (int s=0; s<p; ++s)
	  for (int k1 = 0; k1 < this->KyMax; ++k1)
	    for (int k2 = 0; k2 <= k1; ++k2)
	      for (int k3 = 0; k3 < this->KyMax; ++k3)
		{
		  k4 = k1 + k2 - k3;
		  if (k4 < 0)
		    k4 += this->KyMax;
		  else
		    if (k4 >= this->KyMax)
		      k4 -= this->KyMax;
		  if (k3 > k4)
		    {
		      if (k1 != k2)
			{
			  TmpCoefficient[Pos] = 4.0*Strength;
			}
		      else
			TmpCoefficient[Pos] = 2.0*Strength;
		      if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
			MaxCoefficient = fabs(TmpCoefficient[Pos]);
		      ++Pos;
		    }
		  else
		    if (k3 == k4)
		      {
			if (k1 != k2)
			  TmpCoefficient[Pos] = 2.0*Strength;
			else
			  TmpCoefficient[Pos] = Strength;
			if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
			  MaxCoefficient = fabs(TmpCoefficient[Pos]);
			++Pos;
		      }
		}
      
      this->NbrInteractionFactors = 0;
      this->Q1Value = new int [Pos];
      this->Q2Value = new int [Pos];
      this->Q3Value = new int [Pos];
      this->Q4Value = new int [Pos];
      this->InteractionFactors = new Complex [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;
      
      for (int i=0; i<Lx; ++i)
	for (int s=0; s<p; ++s)
	  for (int k1 = 0; k1 < this->KyMax; ++k1)
	    for (int k2 = 0; k2 <= k1; ++k2)
	      for (int k3 = 0; k3 < this->KyMax; ++k3)
		{
		  k4 = k1 + k2 - k3;
		  if (k4 < 0)
		    k4 += this->KyMax;
		  else
		    if (k4 >= this->KyMax)
		      k4 -= this->KyMax;
		  if (k3 >= k4)
		    {
		      if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
			{
			  this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
			  this->Q1Value[this->NbrInteractionFactors] = Particles->EncodeQuantumNumber(i, k1*p+s, 0, TranslationPhase);
			  this->Q2Value[this->NbrInteractionFactors] = Particles->EncodeQuantumNumber(i, k2*p+s, 0, TranslationPhase);
			  this->Q3Value[this->NbrInteractionFactors] = Particles->EncodeQuantumNumber(i, k3*p+s, 0, TranslationPhase);
			  this->Q4Value[this->NbrInteractionFactors] = Particles->EncodeQuantumNumber(i, k4*p+s, 0, TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "4-body term: "<<this->InteractionFactors[this->NbrInteractionFactors]<< " a^\\dag_"<< this->Q4Value[this->NbrInteractionFactors]<<"a^\\dag_"<<this->Q3Value[this->NbrInteractionFactors]<<"a_"<<this->Q2Value[this->NbrInteractionFactors]<<"a_"<<this->Q1Value[this->NbrInteractionFactors]<<endl;
#endif
			  ++this->NbrInteractionFactors;
			}
		      ++Pos;
		    }
		}
      delete [] TmpCoefficient;
    }
  else // no such interactions
    {
      this->NbrInteractionFactors = 0;
    }
    
  cout << "nbr interaction Twobody= " << this->NbrInteractionFactors << endl;
  cout << "nbr interaction Nbody= " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long ParticleOnLatticeWithKyNBodyDeltaHamiltonian::FastMultiplicationMemory(long allowedMemory)
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  
  this->NbrRealInteractionPerComponent = new ElementIndexType [EffectiveHilbertSpaceDimension];
  this->NbrComplexInteractionPerComponent = new ElementIndexType [EffectiveHilbertSpaceDimension];   
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    {
      this->NbrRealInteractionPerComponent[i] = 0x0;
      this->NbrComplexInteractionPerComponent[i] = 0x0;
    }
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;

  QHEParticlePrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);

  long Memory = 0;
	
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    {
      Memory += this->NbrRealInteractionPerComponent[i];
      Memory += this->NbrComplexInteractionPerComponent[i];
    }
  
  cout << "nbr interaction = " << Memory << endl;
  
  // memory requirement, ignoring the actual storage size of the values of matrix
  // elements, which is assumed small (maybe need to add an estimate, at least)
  long TmpMemory = allowedMemory - (2*sizeof (ElementIndexType) + sizeof (int*) + sizeof(ElementIndexType*)) * EffectiveHilbertSpaceDimension;
  cout << "of which can be stored: "<<(TmpMemory / ((int) (sizeof (int) + sizeof(ElementIndexType))))<<endl;
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(ElementIndexType)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(ElementIndexType)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if (this->Particles->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  // memory requirement, ignoring the actual storage size of the values of matrix
	  // elements, which is assumed small (maybe need to add an estimate, at least, again!)
	  TmpMemory = allowedMemory - (2*sizeof (ElementIndexType) + sizeof (int*) + sizeof(ElementIndexType*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	    {
	      Memory += this->NbrRealInteractionPerComponent[i];
	      Memory += this->NbrComplexInteractionPerComponent[i];
	    }	  
	}
      
      Memory = ((2*sizeof (ElementIndexType) + sizeof (int*) + sizeof(ElementIndexType*)) * ReducedSpaceDimension) + (Memory * (sizeof (int) + sizeof(ElementIndexType)));
      
      if (this->DiskStorageFlag == false)
	{
	  int TotalReducedSpaceDimension = ReducedSpaceDimension;
	  ElementIndexType* TmpNbrRealInteractionPerComponent = new ElementIndexType [TotalReducedSpaceDimension];
	  ElementIndexType* TmpNbrComplexInteractionPerComponent = new ElementIndexType [TotalReducedSpaceDimension];	  
	  int Pos = 0;
	  for (int i = 0; i < ReducedSpaceDimension; ++i)
	    {
	      TmpNbrRealInteractionPerComponent[i] = this->NbrRealInteractionPerComponent[Pos];
	      TmpNbrComplexInteractionPerComponent[i] = this->NbrComplexInteractionPerComponent[Pos];
	      Pos += this->FastMultiplicationStep;
	    }
	  delete[] this->NbrRealInteractionPerComponent;
	  delete[] this->NbrComplexInteractionPerComponent;
	  this->NbrRealInteractionPerComponent = TmpNbrRealInteractionPerComponent;
	  this->NbrComplexInteractionPerComponent = TmpNbrComplexInteractionPerComponent;
	}
    }
  else
    {
      Memory = ((2*sizeof (ElementIndexType) + sizeof (int*) + sizeof(ElementIndexType*)) * EffectiveHilbertSpaceDimension) + (Memory * (sizeof (int) + sizeof(ElementIndexType)));
      this->FastMultiplicationStep = 1;
    }
  
  cout << "reduction factor=" << this->FastMultiplicationStep << endl;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  cout << "final Memory in bytes = " <<Memory<<endl;
  return Memory;    
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long ParticleOnLatticeWithKyNBodyDeltaHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  long Memory = 0;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  int LastComponent = lastComponent + firstComponent;
  
  this->EvaluateMNNBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);
  if (this->TwoBodyFlag == true)
    {
      this->EvaluateMNTwoBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);
    }
  this->EvaluateMNOneBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);
  
  delete TmpParticles;
  return Memory;
}



// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = index of the last component that has to be precalcualted

void ParticleOnLatticeWithKyNBodyDeltaHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{  
  int LastComponent = nbrComponent + firstComponent;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  
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
      int PosR = 0;
      int PosC = this->NbrRealInteractionPerComponent[Pos];
      this->EvaluateMNOneBodyFastMultiplicationComponent(TmpParticles, i, this->InteractionPerComponentIndex[Pos], 
							 this->InteractionPerComponentCoefficientIndex[Pos], PosR,PosC);
      this->EvaluateMNNBodyFastMultiplicationComponent(TmpParticles, i, this->InteractionPerComponentIndex[Pos], 
						       this->InteractionPerComponentCoefficientIndex[Pos], PosR,PosC);
      if (this->TwoBodyFlag == true)
	{
	  this->EvaluateMNTwoBodyFastMultiplicationComponent(TmpParticles, i, this->InteractionPerComponentIndex[Pos], this->InteractionPerComponentCoefficientIndex[Pos],PosR,PosC);
	}
      ++Pos;
    }
    
  delete TmpParticles;
}

// get all indices needed to characterize a completly symmetric tensor, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// sortedIndicesPerSumSymmetryFactor = reference on a array where symmetry factor (aka inverse of the product of the factorial of the number 
//                                      of time each index appears) are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long ParticleOnLatticeWithKyNBodyDeltaHamiltonian::GetAllSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum, double**& sortedIndicesPerSumSymmetryFactor)
{
  int p = Ly/KyMax;
  long** DimensionSymmetricGroup;
  if (nbrValues >= nbrIndices)
    DimensionSymmetricGroup = GetIrreducibleRepresentationDimensionSymmetricGroup(nbrValues);
  else
    DimensionSymmetricGroup = GetIrreducibleRepresentationDimensionSymmetricGroup(nbrIndices);
  long NbrElements = DimensionSymmetricGroup[nbrValues][nbrIndices];
  
	
  int** Indices = new int* [NbrElements];
  int* Sum = new int [NbrElements];
  int Max;
  int Step;
  int Pos = 0;
  int TmpNbrIndices;
  for (int i = nbrValues - 1; i >= 0; --i)
    {
      Step = DimensionSymmetricGroup[i + 1][nbrIndices - 1];
      for (int j = 0; j < Step; ++j)
	{
	  Indices[Pos] = new int [nbrIndices];
	  Indices[Pos][0] = i;
	  Sum[Pos] = i;
	  ++Pos;
	}
    }
  for (int i = 1; i < nbrIndices; ++i)
    {
      int Pos = 0;
      TmpNbrIndices = nbrIndices - i - 1;
      while (Pos < NbrElements)
	{
	  Max = Indices[Pos][i - 1];
	  Step = DimensionSymmetricGroup[Max + 1][TmpNbrIndices];
	  for (int j = 0; j < Step; ++j)
	    {
	      Indices[Pos][i] = Max;
	      Sum[Pos] += Max;
	      ++Pos;
	    }
	  --Max;
	  for (; Max >= 0; --Max)
	    {
	      Step = DimensionSymmetricGroup[Max + 1][TmpNbrIndices];
	      for (int j = 0; j < Step; ++j)
		{
		  Indices[Pos][i] = Max;
		  Sum[Pos] += Max;
		  ++Pos;
		}
	    }
	}
    }
  int * TmpSum =  new int [p*NbrElements];
  
  for (int s = 0 ; s< p ; s++)
    for (int i = 0; i < NbrElements; ++i)
      {
	TmpSum[p*i+s] = ((Sum[i]*p +s) % this->Ly);
      }
  
  delete[] Sum;
  
  int MaxSum = this->Ly-1;
  long* TmpPos = new long [MaxSum + 1];
  nbrSortedIndicesPerSum = new int [MaxSum + 1];
  sortedIndicesPerSum = new int* [MaxSum + 1];
  sortedIndicesPerSumSymmetryFactor = new double* [MaxSum + 1];
  for (int i = 0; i <= MaxSum; ++i)
    nbrSortedIndicesPerSum[i] = 0;
  for (int i = 0; i < p*NbrElements; ++i)
    {
      ++nbrSortedIndicesPerSum[TmpSum[i]];
    }
  for (int i = 0; i <= MaxSum; ++i)
    {
      TmpPos[i] = 0l;
      sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * nbrIndices];
      sortedIndicesPerSumSymmetryFactor[i] = new double [nbrSortedIndicesPerSum[i]];
      nbrSortedIndicesPerSum[i] = 0;
    }
    
   for (int s = 0 ; s< p ; s++)
     for (int i = 0; i < NbrElements; ++i)
       {
	 Pos = TmpSum[p*i+s];
	 Max = nbrSortedIndicesPerSum[Pos];
	 int* TmpIndices = Indices[i];
	 for (int j = 0; j < nbrIndices; ++j)
	   {
	     sortedIndicesPerSum[Pos][TmpPos[Pos]] = (TmpIndices[j]*p +s) % this->Ly;
	     ++TmpPos[Pos];
	   }
	 
	 double& SymmetryFactor = sortedIndicesPerSumSymmetryFactor[Pos][Max];
	 SymmetryFactor = 1.0;
	 for (int j = 1; j < nbrIndices; ++j)
	   {
	     int TmpSymmetryFactor = 1;
	     while ((j < nbrIndices) && (TmpIndices[j - 1] == TmpIndices[j]))
	       {
		 ++TmpSymmetryFactor;
		 ++j;
	       }
	     if (TmpSymmetryFactor != 1)
	       for (int k = 2; k <= TmpSymmetryFactor; ++k)
		 SymmetryFactor *= (double) k;
	   }
	 SymmetryFactor = 1.0 / SymmetryFactor;
	 ++nbrSortedIndicesPerSum[Pos];
       }
   
   for (int i = 0; i < NbrElements; ++i)
     { 
       delete[] Indices[i];	
     }
   delete[] TmpPos;
   delete[] TmpSum;
   delete[] Indices;
   
   for (int i = 0; i <= nbrValues; ++i)
     delete[] DimensionSymmetricGroup[i];
   delete[] DimensionSymmetricGroup;
   return NbrElements;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnLatticeWithKyNBodyDeltaHamiltonian::LowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  double TmpInteractionRe,TmpInteractionIm;
  int LastComponent = firstComponent + nbrComponent;
  int* TmpIndexArray;
  ElementIndexType* TmpCoefficientIndexArray;
  double TmpRe, TmpIm;
  Complex *TmpCPtr;
  int TmpNbrRealInteraction;
  int TmpNbrComplexInteraction;
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
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[Pos];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[Pos];
      TmpRe = vSource[l].Re;
      TmpIm = vSource[l].Im;
      int Pos2=0;
      for (; Pos2 < TmpNbrRealInteraction; ++Pos2)
	{
	  vDestination.Re(TmpIndexArray[Pos2]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*TmpRe;
	  vDestination.Im(TmpIndexArray[Pos2]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*TmpIm;
	}
      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos2)
	{
	  TmpCPtr= &(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]);
	  vDestination.Re(TmpIndexArray[Pos2]) +=  TmpCPtr->Re*TmpRe-TmpCPtr->Im*TmpIm;
	  vDestination.Im(TmpIndexArray[Pos2]) +=  TmpCPtr->Re*TmpIm+TmpCPtr->Im*TmpRe;		  
	}
      vDestination.Re(l) += this->HamiltonianShift * TmpRe;
      vDestination.Im(l) += this->HamiltonianShift * TmpIm;	      
      l += this->FastMultiplicationStep;
      ++Pos;
    }

  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	if (this->TwoBodyFlag == true)
	{
	  for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	    {
	      this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	      this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	    }
	  
	}
	else
	  {
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	      }
	  }
	this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + k, LastComponent, this->FastMultiplicationStep, vSource, vDestination);
      }
  
  delete TmpParticles;
  return vDestination;
}



// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* ParticleOnLatticeWithKyNBodyDeltaHamiltonian::ConjugateLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  
  if (this->FastMultiplicationFlag == false)
    {
      Complex * TmpCoefficients = new Complex [nbrVectors];
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      this->EvaluateMNOneBodyConjugateAddMultiplyComponent(TmpParticles, firstComponent,LastComponent,1, vSources, vDestinations , nbrVectors);

      if (this->TwoBodyFlag == true)
	{
	  cout <<"Not finished fonction ParticleOnLatticeWithKyNBodyDeltaHamiltonian::ConjugateLowLevelMultipleAddMultiply " <<endl;
	  exit(1);
	  
	  for (int i = firstComponent ; i < LastComponent; i++)
	    {
	      this->EvaluateMNTwoBodyConjugateAddMultiplyComponent(TmpParticles, i ,  vSources, vDestinations, nbrVectors,  TmpCoefficients);
	      
	       //this->EvaluateMNNBodyConjugateAddMultiplyComponent(TmpParticles, i ,  vSources, vDestinations, nbrVectors,  TmpCoefficients);
	    }
	}
      else
	{
	  cout <<"Not finished fonction ParticleOnLatticeWithKyNBodyDeltaHamiltonian::ConjugateLowLevelMultipleAddMultiply " <<endl;
	  exit(1);
	  
	  //for (int i =   firstComponent; i <LastComponent; i++)
	  //this->EvaluateMNNBodyConjugateAddMultiplyComponent(TmpParticles, i ,  vSources, vDestinations, nbrVectors,  TmpCoefficients);
	  
	}
      delete [] TmpCoefficients;
      delete TmpParticles;
    }	  
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->ConjugateLowLevelMultipleAddMultiplyFastMultiply(vSources,vDestinations, nbrVectors,firstComponent, LastComponent);
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->ConjugateLowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	}
    }
  return vDestinations;
}


