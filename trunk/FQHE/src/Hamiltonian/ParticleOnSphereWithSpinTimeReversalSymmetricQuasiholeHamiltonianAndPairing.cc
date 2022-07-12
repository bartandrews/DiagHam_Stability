////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                      Class author: Cecile Repellin                         //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to Laughlin qh on a sphere with      //
//        SU(2) spin with opposite magnetic field for each species            //
//                            and pairing                                     //
//                                                                            //
//                        last modification : 11/05/2016                      //
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


#include "Hamiltonian/ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"
#include "GeneralTools/ArrayTools.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;

//default constructor
ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing()
{
  this->HermitianSymmetryFlag = false;
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// lzmax = maximum Lz value reached by a particle in the state
// architecture = architecture to use for precalculation
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyPotentialPairing =  one-body pairing term (sorted from component on the lowest Lz state to component on the highest Lz state), on site, symmetric spin up / spin down
// chargingEnergy = factor in front of the charging energy (i.e 1/(2C))
// averageNumberParticles = average number of particles in the system
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing(QuasiholeOnSphereWithSpinAndPairing* particles, int lzmax, double* onebodyPotentialUpUp, double* onebodyPotentialDownDown, double* onebodyPotentialPairing, double chargingEnergy, double averageNumberParticles, AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  
  this->ChargingEnergy = chargingEnergy;
  this->AverageNumberParticles = averageNumberParticles;
  this->HamiltonianShift =  0.0;
  long MinIndex;
  long MaxIndex;
  
  this->HermitianSymmetryFlag = true;
//   this->HermitianSymmetryFlag = false;
  
  
  this->Architecture = architecture;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->Memory = memory;
  this->FastMultiplicationFlag = false;
  this->DiskStorageFlag = onDiskCacheFlag;

  this->OneBodyInteractionFactorsupup = 0;
  if (onebodyPotentialUpUp != 0)
    {
      this->OneBodyInteractionFactorsupup = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsupup[i] = onebodyPotentialUpUp[i];
    }
  
  this->OneBodyInteractionFactorsdowndown = 0;
  if (onebodyPotentialDownDown != 0)
    {
      this->OneBodyInteractionFactorsdowndown = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsdowndown[i] = onebodyPotentialDownDown[i];
    }
  this->OneBodyInteractionFactorsPairing = 0;
  if (onebodyPotentialPairing != 0)
    {
      this->OneBodyInteractionFactorsPairing = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->OneBodyInteractionFactorsPairing[i] = onebodyPotentialPairing[i];
	}
    }
  
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      if (TmpMemory < 1024l)
	cout  << "fast = " <<  TmpMemory << "b ";
      else
	if (TmpMemory < (1l << 20))
	  cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	else
	  if (TmpMemory < (1l << 30))
	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	  else
	    cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
      cout << endl;
      this->EnableFastMultiplication();
    }
  else
    {
      if (this->Architecture->HasAutoLoadBalancing())
	{
	  this->FastMultiplicationMemory(0l);
	}
    }
}

// destructor
//

ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::~ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing() 
{
  if (this->OneBodyInteractionFactorsupup != 0)
    delete[] this->OneBodyInteractionFactorsupup;
  if (this->OneBodyInteractionFactorsdowndown != 0)
    delete[] this->OneBodyInteractionFactorsdowndown;
  if (this->OneBodyInteractionFactorsPairing != 0)
    delete [] this->OneBodyInteractionFactorsPairing;
  
  
  
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (QuasiholeOnSphereWithSpinAndPairing*) hilbertSpace;
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift += shift;
}
  
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
													     int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  double TmpCoefficient;
    
  if (this->FastMultiplicationFlag == false)
    {
      int MaximalNumberCouplingElements = ((QuasiholeOnSphereWithSpinAndPairing*) (this->Particles))->GetMaximalNumberCouplingElements();
      int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
      double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
      QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
      int index = firstComponent;
      int NbrElements;
      double ChargeContribution;
      double TmpCoef;
      
      for (int i = 0; i < nbrComponent; ++i)
	{
	  TmpCoef = vSource[index];
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if (TmpCoefficient != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
		      NbrElements = TmpParticles->AduAdd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
		    }	
		}
	    }
	      
	  if (this->OneBodyInteractionFactorsupup != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
	      {
		NbrElements = TmpParticles->AduAu(index, lz, TmpLeftIndices, TmpInteractionElements);
		for (int j = 0; j < NbrElements; ++j)
		  vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
	      }
	    }
	  }
	      
	  if (this->OneBodyInteractionFactorsdowndown != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
	      {
		NbrElements = TmpParticles->AddAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		for (int j = 0; j < NbrElements; ++j)
		  vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
	      }
	    }
	  }
	  
	  if (this->ChargingEnergy != 0.0)
	  {
	    int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	    if (TmpTotalNbrParticles != this->AverageNumberParticles)
	    {
	      ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
	      ChargeContribution *= ChargeContribution;
	      ChargeContribution *= (vSource[index] * this->ChargingEnergy);
	      vDestination[index] += ChargeContribution;
	    }
	  }
	  
	  ++index;
	}
      
      
      if (this->HamiltonianShift != 0.0)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    vDestination[i] += this->HamiltonianShift * vSource[i];
	}
      delete TmpParticles;
      delete[] TmpLeftIndices;
      delete[] TmpInteractionElements;
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
	      TmpCoefficient = vSource[k];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * TmpCoefficient;
	      vDestination[k++] += this->HamiltonianShift * TmpCoefficient;
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      cout << "Error: ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelAddMultiplyDiskStorage not implemented" << endl;
// 	      this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
														int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  Complex TmpCoefficient;
    
  if (this->FastMultiplicationFlag == false)
    {
      int MaximalNumberCouplingElements = ((QuasiholeOnSphereWithSpinAndPairing*) (this->Particles))->GetMaximalNumberCouplingElements();
      int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
      double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
      QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
      int index = firstComponent;
      int NbrElements;
      double ChargeContribution;
      Complex TmpCoef;
      
      for (int i = 0; i < nbrComponent; ++i)
	{
	  TmpCoef = vSource[index];
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if (TmpCoefficient != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			vDestination[TmpLeftIndices[j]] += TmpCoef * (TmpInteractionElements[j] * TmpCoefficient);
		      NbrElements = TmpParticles->AduAdd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			vDestination[TmpLeftIndices[j]] += TmpCoef * (TmpInteractionElements[j] * TmpCoefficient);
		    }	
		}
	    }
	      
	  if (this->OneBodyInteractionFactorsupup != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
		    {
		      NbrElements = TmpParticles->AduAu(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			vDestination[TmpLeftIndices[j]] += TmpCoef * (TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
		    }
		}
	    }
	  
	  if (this->OneBodyInteractionFactorsdowndown != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
		    {
		      NbrElements = TmpParticles->AddAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			vDestination[TmpLeftIndices[j]] += TmpCoef * (TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
		    }
		}
	    }
	  
	  if (this->ChargingEnergy != 0.0)
	    {
	      int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	      if (TmpTotalNbrParticles != this->AverageNumberParticles)
		{
		  ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
		  ChargeContribution *= ChargeContribution;
		  ChargeContribution *= this->ChargingEnergy;
		  vDestination[index] += vSource[index] * ChargeContribution;
		}
	    }
	  
	  ++index;
	}
      
      
      if (this->HamiltonianShift != 0.0)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    vDestination[i] += this->HamiltonianShift * vSource[i];
	}
      delete TmpParticles;
      delete[] TmpLeftIndices;
      delete[] TmpInteractionElements;
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
	      TmpCoefficient = vSource[k];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * TmpCoefficient;
	      vDestination[k++] += this->HamiltonianShift * TmpCoefficient;
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      cout << "Error: ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelAddMultiplyDiskStorage not implemented" << endl;
// 	      this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
  return vDestination;
}



// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
														      int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      int MaximalNumberCouplingElements = ((QuasiholeOnSphereWithSpinAndPairing*) (this->Particles))->GetMaximalNumberCouplingElements();
      int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
      double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
      QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
      
      int index = firstComponent;
      int NbrElements;
      double ChargeContribution;
      double TmpCoef;
      double TmpSum;
      double TmpCoefficient;
      for (int i = 0; i < nbrComponent; ++i)
	{
	  TmpSum = 0.0;
	  TmpCoef = vSource[index];
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if (TmpCoefficient != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
		      {
			vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
			TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * TmpCoefficient);
		      }
		    }	
		}
	    }
	      
	  if (this->OneBodyInteractionFactorsupup != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
	      {
		NbrElements = TmpParticles->AduAu(index, lz, TmpLeftIndices, TmpInteractionElements);
		for (int j = 0; j < NbrElements; ++j)
		{
		  if (TmpLeftIndices[j] <= index)
		  {
		    vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
		    if (TmpLeftIndices[j] < index)
		      TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
		  }
		}
	      }
	    }
	  }
	      
	  if (this->OneBodyInteractionFactorsdowndown != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
	      {
		NbrElements = TmpParticles->AddAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		for (int j = 0; j < NbrElements; ++j)
		{
		  if (TmpLeftIndices[j] <= index)
		  {
		    vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
		    if (TmpLeftIndices[j] < index)
		      TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
		  }
		}
	      }
	    }
	  }
	  
	  if (this->ChargingEnergy != 0.0)
	  {
	    int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	    if (TmpTotalNbrParticles != this->AverageNumberParticles)
	    {
	      ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
	      ChargeContribution *= ChargeContribution;
	      ChargeContribution *= (vSource[index] * this->ChargingEnergy);
	      vDestination[index] += ChargeContribution;
	    }
	  }
	  vDestination[index] += TmpSum;
	  ++index;
	}
      
      if (this->HamiltonianShift != 0.0)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    vDestination[i] += this->HamiltonianShift * vSource[i];
	}
      delete TmpParticles;
      delete[] TmpLeftIndices;
      delete[] TmpInteractionElements;
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
		  TmpSum += TmpCoefficientArray[j] * vSource[TmpIndexArray[j]];
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
		}
	      TmpSum += this->HamiltonianShift * Coefficient;
	      vDestination[k++] += TmpSum;
	    }
	}
      else
	{
	  this->HermitianLowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	}
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
															 int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      int MaximalNumberCouplingElements = ((QuasiholeOnSphereWithSpinAndPairing*) (this->Particles))->GetMaximalNumberCouplingElements();
      int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
      double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
      QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
      
      int index = firstComponent;
      int NbrElements;
      double ChargeContribution;
      Complex TmpCoef;
      Complex TmpSum;
      double TmpCoefficient;
      for (int i = 0; i < nbrComponent; ++i)
	{
	  TmpSum = 0.0;
	  TmpCoef = vSource[index];
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if (TmpCoefficient != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
		      {
			vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
			TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * TmpCoefficient);
		      }
		    }	
		}
	    }
	      
	  if (this->OneBodyInteractionFactorsupup != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
		    {
		      NbrElements = TmpParticles->AduAu(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  if (TmpLeftIndices[j] <= index)
			    {
			      vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
			      if (TmpLeftIndices[j] < index)
				TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
			    }
			}
		    }
		}
	    }
	  
	  if (this->OneBodyInteractionFactorsdowndown != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
		    {
		      NbrElements = TmpParticles->AddAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  if (TmpLeftIndices[j] <= index)
			    {
			      vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
			      if (TmpLeftIndices[j] < index)
				TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
			    }
			}
		    }
		}
	    }
	  
	  if (this->ChargingEnergy != 0.0)
	    {
	      int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	      if (TmpTotalNbrParticles != this->AverageNumberParticles)
		{
		  ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
		  ChargeContribution *= ChargeContribution;
		  ChargeContribution *= this->ChargingEnergy;
		  vDestination[index] += vSource[index] * ChargeContribution;
		}
	    }
	  vDestination[index] += TmpSum;
	  ++index;
	}
      
      if (this->HamiltonianShift != 0.0)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    vDestination[i] += this->HamiltonianShift * vSource[i];
	}
      delete TmpParticles;
      delete[] TmpLeftIndices;
      delete[] TmpInteractionElements;
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
		  TmpSum += TmpCoefficientArray[j] * vSource[TmpIndexArray[j]];
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
		}
	      TmpSum += this->HamiltonianShift * Coefficient;
	      vDestination[k++] += TmpSum;
	    }
	}
      else
	{
	  this->HermitianLowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	}
    }
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  
  
  if (this->FastMultiplicationFlag == false)
    {
      int MaximalNumberCouplingElements = ((QuasiholeOnSphereWithSpinAndPairing*) (this->Particles))->GetMaximalNumberCouplingElements();
      int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
      double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
      QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
      int index = firstComponent;
      int NbrElements;
      double ChargeContribution;
      double TmpOneBodyInteraction;
      for (int i = 0; i < nbrComponent; ++i)
	{
	  for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsPairing != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsPairing[lz];    
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			  vDestinations[k][TmpLeftIndices[j]] += (vSources[k][index] * TmpInteractionElements[j] * TmpOneBodyInteraction);
		      NbrElements = TmpParticles->AduAdd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			  vDestinations[k][TmpLeftIndices[j]] += (vSources[k][index] * TmpInteractionElements[j] * TmpOneBodyInteraction);
		    }
		}
	      
	      if (this->OneBodyInteractionFactorsupup != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsupup[lz];
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AduAu(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			  vDestinations[k][TmpLeftIndices[j]] += (vSources[k][index] * TmpInteractionElements[j] * TmpOneBodyInteraction);
		    }
		}
	      
	      if (this->OneBodyInteractionFactorsdowndown != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsdowndown[lz];
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AddAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			  vDestinations[k][TmpLeftIndices[j]] += (vSources[k][index] * TmpInteractionElements[j] * TmpOneBodyInteraction);
		    }
		}
	    }
	    
	  if (this->ChargingEnergy != 0.0)
	  {
	    int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	    if (TmpTotalNbrParticles != this->AverageNumberParticles)
	    {
	      ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
	      ChargeContribution *= ChargeContribution;
	      ChargeContribution *= this->ChargingEnergy;
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][index] += (ChargeContribution * vSources[k][index]);
	    }
	  }
	  
	  ++index;
	}
      
      if (this->HamiltonianShift != 0.0)
	{
	  for (int k= 0; k < nbrVectors; ++k)
	    {
	      RealVector& TmpDestination = vDestinations[k];
	      RealVector& TmpSource = vSources[k];
	      for (int i = firstComponent; i < LastComponent; ++i)
		TmpDestination[i] += this->HamiltonianShift * TmpSource[i];
	    }
	}
       
      delete[] TmpLeftIndices;
      delete[] TmpInteractionElements;
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  double* Coefficient2 = new double [nbrVectors];
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int j;
	  int Pos;
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
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      cout << "Error: ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelMultipleAddMultiplyDiskStorage not implemented" << endl;
// 	      this->LowLevelMultipleAddMultiplyDiskStorage(vSource, vDestination, nbrVectors, firstComponent, nbrComponent);
	    }
	}
    }
   
  return vDestinations;
}


// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  
  
  if (this->FastMultiplicationFlag == false)
    {
      int MaximalNumberCouplingElements = ((QuasiholeOnSphereWithSpinAndPairing*) (this->Particles))->GetMaximalNumberCouplingElements();
      int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
      double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
      QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
      int index = firstComponent;
      int NbrElements;
      double ChargeContribution;
      double TmpOneBodyInteraction;
      for (int i = 0; i < nbrComponent; ++i)
	{
	  for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsPairing != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsPairing[lz];    
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			  vDestinations[k][TmpLeftIndices[j]] += vSources[k][index] * (TmpInteractionElements[j] * TmpOneBodyInteraction);
		      NbrElements = TmpParticles->AduAdd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			  vDestinations[k][TmpLeftIndices[j]] += vSources[k][index] * (TmpInteractionElements[j] * TmpOneBodyInteraction);
		    }
		}
	      
	      if (this->OneBodyInteractionFactorsupup != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsupup[lz];
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AduAu(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			  vDestinations[k][TmpLeftIndices[j]] += vSources[k][index] * (TmpInteractionElements[j] * TmpOneBodyInteraction);
		    }
		}
	      
	      if (this->OneBodyInteractionFactorsdowndown != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsdowndown[lz];
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AddAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			  vDestinations[k][TmpLeftIndices[j]] += vSources[k][index] * (TmpInteractionElements[j] * TmpOneBodyInteraction);
		    }
		}
	    }
	    
	  if (this->ChargingEnergy != 0.0)
	    {
	      int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	      if (TmpTotalNbrParticles != this->AverageNumberParticles)
		{
		  ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
		  ChargeContribution *= ChargeContribution;
		  ChargeContribution *= this->ChargingEnergy;
		  for (int k = 0; k < nbrVectors; ++k)
		    vDestinations[k][index] += (ChargeContribution * vSources[k][index]);
		}
	    }
	  
	  ++index;
	}
      
      if (this->HamiltonianShift != 0.0)
	{
	  for (int k= 0; k < nbrVectors; ++k)
	    {
	      ComplexVector& TmpDestination = vDestinations[k];
	      ComplexVector& TmpSource = vSources[k];
	      for (int i = firstComponent; i < LastComponent; ++i)
		TmpDestination[i] += this->HamiltonianShift * TmpSource[i];
	    }
	}
       
      delete[] TmpLeftIndices;
      delete[] TmpInteractionElements;
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  Complex* Coefficient2 = new Complex [nbrVectors];
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int j;
	  int Pos;
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
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      cout << "Error: ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelMultipleAddMultiplyDiskStorage not implemented" << endl;
// 	      this->LowLevelMultipleAddMultiplyDiskStorage(vSource, vDestination, nbrVectors, firstComponent, nbrComponent);
	    }
	}
    }
   
  return vDestinations;
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

RealVector* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::HermitianLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
														int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      int MaximalNumberCouplingElements = ((QuasiholeOnSphereWithSpinAndPairing*) (this->Particles))->GetMaximalNumberCouplingElements();
      int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
      double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
      QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
      
      int index = firstComponent;
      int NbrElements;
      double ChargeContribution;
      double TmpOneBodyInteraction;
      double* TmpSum = new double[nbrVectors];
      for (int i = 0; i < nbrComponent; ++i)
	{
	  for (int k = 0; k < nbrVectors; ++k)
	    TmpSum[k] = 0.0;
	  for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsPairing != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsPairing[lz];    
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][TmpLeftIndices[j]] += (vSources[k][index] * TmpInteractionElements[j] * TmpOneBodyInteraction);
			  TmpSum[k] += (vSources[k][TmpLeftIndices[j]] * TmpInteractionElements[j] * TmpOneBodyInteraction);
			}
		      
		    }
		}
	      
	      if (this->OneBodyInteractionFactorsupup != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsupup[lz];
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AduAu(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
		      {
			if (TmpLeftIndices[j] <= index)
			{
			  for (int k = 0; k < nbrVectors; ++k)
			  {
			    vDestinations[k][TmpLeftIndices[j]] += (vSources[k][index] * TmpInteractionElements[j] * TmpOneBodyInteraction);
			    if (TmpLeftIndices[j] < index)
			      TmpSum[k] += (vSources[k][TmpLeftIndices[j]] * TmpInteractionElements[j] * TmpOneBodyInteraction);
			  }
			}
		      }
		    }
		}
	      
	      if (this->OneBodyInteractionFactorsdowndown != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsdowndown[lz];
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AddAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
		      {
			if (TmpLeftIndices[j] <= index)
			{
			for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][TmpLeftIndices[j]] += (vSources[k][index] * TmpInteractionElements[j] * TmpOneBodyInteraction);
			  if (TmpLeftIndices[j] < index)
			    TmpSum[k] += (vSources[k][TmpLeftIndices[j]] * TmpInteractionElements[j] * TmpOneBodyInteraction);
			}
			}
		      }
		    }
		}
	    }
	    
	  if (this->ChargingEnergy != 0.0)
	  {
	    int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	    if (TmpTotalNbrParticles != this->AverageNumberParticles)
	    {
	      ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
	      ChargeContribution *= ChargeContribution;
	      ChargeContribution *= this->ChargingEnergy;
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][index] += (ChargeContribution * vSources[k][index]);
	    }
	  }
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][index] += TmpSum[k];
	  ++index;
	}
      
      if (this->HamiltonianShift != 0.0)
	{
	  for (int k= 0; k < nbrVectors; ++k)
	    {
	      RealVector& TmpDestination = vDestinations[k];
	      RealVector& TmpSource = vSources[k];
	      for (int i = firstComponent; i < LastComponent; ++i)
		TmpDestination[i] += this->HamiltonianShift * TmpSource[i];
	    }
	}
      
      delete[] TmpSum;
      delete[] TmpLeftIndices;
      delete[] TmpInteractionElements;
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
		      TmpSum[l] += Coefficient * vSources[l][Pos];
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

// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
																 int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      int MaximalNumberCouplingElements = ((QuasiholeOnSphereWithSpinAndPairing*) (this->Particles))->GetMaximalNumberCouplingElements();
      int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
      double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
      QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
      
      int index = firstComponent;
      int NbrElements;
      double ChargeContribution;
      double TmpOneBodyInteraction;
      Complex* TmpSum = new Complex[nbrVectors];
      for (int i = 0; i < nbrComponent; ++i)
	{
	  for (int k = 0; k < nbrVectors; ++k)
	    TmpSum[k] = 0.0;
	  for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsPairing != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsPairing[lz];    
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][TmpLeftIndices[j]] += vSources[k][index] * (TmpInteractionElements[j] * TmpOneBodyInteraction);
			  TmpSum[k] += vSources[k][TmpLeftIndices[j]] * (TmpInteractionElements[j] * TmpOneBodyInteraction);
			}
		      
		    }
		}
	      
	      if (this->OneBodyInteractionFactorsupup != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsupup[lz];
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AduAu(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  if (TmpLeftIndices[j] <= index)
			    {
			      for (int k = 0; k < nbrVectors; ++k)
				{
				  vDestinations[k][TmpLeftIndices[j]] += vSources[k][index] * (TmpInteractionElements[j] * TmpOneBodyInteraction);
				  if (TmpLeftIndices[j] < index)
				    TmpSum[k] += vSources[k][TmpLeftIndices[j]] * (TmpInteractionElements[j] * TmpOneBodyInteraction);
				}
			    }
			}
		    }
		}
	      
	      if (this->OneBodyInteractionFactorsdowndown != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsdowndown[lz];
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AddAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  if (TmpLeftIndices[j] <= index)
			    {
			      for (int k = 0; k < nbrVectors; ++k)
				{
				  vDestinations[k][TmpLeftIndices[j]] += vSources[k][index] * (TmpInteractionElements[j] * TmpOneBodyInteraction);
				  if (TmpLeftIndices[j] < index)
				    TmpSum[k] += vSources[k][TmpLeftIndices[j]] * (TmpInteractionElements[j] * TmpOneBodyInteraction);
				}
			    }
			}
		    }
		}
	    }
	    
	  if (this->ChargingEnergy != 0.0)
	    {
	      int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	      if (TmpTotalNbrParticles != this->AverageNumberParticles)
		{
		  ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
		  ChargeContribution *= ChargeContribution;
		  ChargeContribution *= this->ChargingEnergy;
		  for (int k = 0; k < nbrVectors; ++k)
		    vDestinations[k][index] += (ChargeContribution * vSources[k][index]);
		}
	    }
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][index] += TmpSum[k];
	  ++index;
	}
      
      if (this->HamiltonianShift != 0.0)
	{
	  for (int k= 0; k < nbrVectors; ++k)
	    {
	      ComplexVector& TmpDestination = vDestinations[k];
	      ComplexVector& TmpSource = vSources[k];
	      for (int i = firstComponent; i < LastComponent; ++i)
		TmpDestination[i] += this->HamiltonianShift * TmpSource[i];
	    }
	}
      
      delete[] TmpSum;
      delete[] TmpLeftIndices;
      delete[] TmpInteractionElements;
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double Coefficient;
	  Complex* Coefficient2 = new Complex [nbrVectors];
	  double* TmpCoefficientArray; 
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
		      TmpSum[l] += Coefficient * vSources[l][Pos];
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

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::GetHilbertSpace ()
{
  return this->Particles;
}


// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that has to be precalcualted
// return value = number of non-zero matrix element

long ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::PartialFastMultiplicationMemory(int firstComponent, int nbrComponent)
{
  int NbrElements;
  double TmpCoefficient;
  long Memory = 0l;
  
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  cout << "MaximalNumberCouplingElements = " << MaximalNumberCouplingElements << endl;
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  int LastComponent = nbrComponent + firstComponent;
  int* TmpCounting = new int [TmpParticles->GetHilbertSpaceDimension()];
  long TmpTotal = 0l;
  int CurrentNbrCounting = 0;
  int Tmp = 0;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      if (this->OneBodyInteractionFactorsPairing != 0)
	{
	  for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
	      if (TmpCoefficient != 0.0)
		{
		  NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		  Memory += NbrElements;
		  TmpTotal += NbrElements;
		  this->NbrInteractionPerComponent[i - this->PrecalculationShift] += NbrElements;
		  
		  if (this->HermitianSymmetryFlag == false)
		  {
		    NbrElements = TmpParticles->AduAdd(i, lz, TmpLeftIndices, TmpInteractionElements);
		    Memory += NbrElements;
		    TmpTotal += NbrElements;
		    this->NbrInteractionPerComponent[i - this->PrecalculationShift] += NbrElements;
		  }
		}	
	    }
	}
      for (int k = 0; k < CurrentNbrCounting; ++k)
	TmpCounting[k] = 0;
      CurrentNbrCounting = 0;
      if (this->OneBodyInteractionFactorsupup != 0)
	{
	  for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
		{
		  NbrElements = TmpParticles->AduAu(i, lz, TmpLeftIndices, TmpInteractionElements);
		  for (int k = 0; k < NbrElements; ++k)
		    {
		      if ((this->HermitianSymmetryFlag == false) || (TmpLeftIndices[k] <= i))
			CurrentNbrCounting += SearchInSortedArrayAndInsert<int>(TmpLeftIndices[k], TmpCounting, CurrentNbrCounting);
		    }
		}
	    }
	}      
      if (this->OneBodyInteractionFactorsdowndown != 0)
	{
	  for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
		{
		  NbrElements = TmpParticles->AddAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		  for (int k = 0; k < NbrElements; ++k)
		    {
		      if ((this->HermitianSymmetryFlag == false) || (TmpLeftIndices[k] <= i))
			CurrentNbrCounting += SearchInSortedArrayAndInsert<int>(TmpLeftIndices[k], TmpCounting, CurrentNbrCounting);
		    }
		}
	    }
	}
      
      int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
      if ((this->ChargingEnergy != 0.0) && (TmpTotalNbrParticles != this->AverageNumberParticles))
	{
	  CurrentNbrCounting += SearchInSortedArrayAndInsert<int>(i, TmpCounting, CurrentNbrCounting);
	}  
      TmpTotal += CurrentNbrCounting;
      Memory += CurrentNbrCounting;
      this->NbrInteractionPerComponent[i - this->PrecalculationShift] += CurrentNbrCounting;
    }
    
  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  delete[] TmpCounting;
  delete TmpParticles;
  return Memory;
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = index of the last component that has to be precalcualted

void ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{  
  double TmpCoefficient;
  double ChargeContribution;
  int NbrElements;
  int LastComponent = nbrComponent + firstComponent;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();

  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  
  long Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
    
  int* TmpCounting = new int [TmpParticles->GetHilbertSpaceDimension()];
  double* TmpCoefficients = new double [TmpParticles->GetHilbertSpaceDimension()];
  int CurrentNbrCounting = 0;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      int TmpShiftedIndex = i + this->PrecalculationShift;
      for (int k = 0; k < CurrentNbrCounting; ++k)
	TmpCounting[k] = 0;
      CurrentNbrCounting = 0;
      int position = 0;
      if (this->OneBodyInteractionFactorsPairing != 0)
	{
	  for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
	      if (TmpCoefficient != 0.0)
		{
		  NbrElements = TmpParticles->AuAd(TmpShiftedIndex, lz, TmpLeftIndices, TmpInteractionElements);
		  for (int j = 0; j < NbrElements; ++j)
		    {
		      TmpCounting[CurrentNbrCounting] = TmpLeftIndices[j];
		      TmpCoefficients[CurrentNbrCounting] = TmpCoefficient * TmpInteractionElements[j];
		      ++CurrentNbrCounting;
		    }
		  if (this->HermitianSymmetryFlag == false)
		    {
		      NbrElements = TmpParticles->AduAdd(TmpShiftedIndex, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  TmpCounting[CurrentNbrCounting] = TmpLeftIndices[j];
			  TmpCoefficients[CurrentNbrCounting] = TmpCoefficient * TmpInteractionElements[j];
			  ++CurrentNbrCounting;
			}
		    }
		}
	    }
	}
      
      SortArrayUpOrdering (TmpCounting, TmpCoefficients, CurrentNbrCounting);
      if (this->OneBodyInteractionFactorsupup != 0)
	{
	  for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
		{
		  NbrElements = TmpParticles->AduAu(TmpShiftedIndex, lz, TmpLeftIndices, TmpInteractionElements);
		  for (int j = 0; j < NbrElements; ++j)
		    {
		      if (this->HermitianSymmetryFlag == false)
			{
			  CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyInteractionFactorsupup[lz] * TmpInteractionElements[j]);
			}
		      else
			{
			  if (TmpLeftIndices[j] < TmpShiftedIndex)
			    CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyInteractionFactorsupup[lz] * TmpInteractionElements[j]);
			  if (TmpLeftIndices[j] == TmpShiftedIndex)
			    CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyInteractionFactorsupup[lz] * TmpInteractionElements[j] * 0.5);
			}
		    }
		}
	    }
	}      
      if (this->OneBodyInteractionFactorsdowndown != 0)
	{
	  for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
		{
		  NbrElements = TmpParticles->AddAd(TmpShiftedIndex, lz, TmpLeftIndices, TmpInteractionElements);
		  for (int j = 0; j < NbrElements; ++j)
		    {
		      if (this->HermitianSymmetryFlag == false)
			{
			  CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyInteractionFactorsdowndown[lz] * TmpInteractionElements[j]);
			}
		      else
			{
			  if (TmpLeftIndices[j] < TmpShiftedIndex)
			    CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyInteractionFactorsdowndown[lz] * TmpInteractionElements[j]);
			  if (TmpLeftIndices[j] == TmpShiftedIndex)
			    CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyInteractionFactorsdowndown[lz] * TmpInteractionElements[j] * 0.5);
			}
		    }
		}
	    }
	}
      
      if (this->ChargingEnergy != 0.0)
	{
	  int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(TmpShiftedIndex);
	  if (TmpTotalNbrParticles != this->AverageNumberParticles)
	    {
	      ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
	      ChargeContribution *= ChargeContribution;
	      ChargeContribution *= this->ChargingEnergy;
	      if (this->HermitianSymmetryFlag)
		ChargeContribution *= 0.5;
	      CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpShiftedIndex, TmpCounting, TmpCoefficients, CurrentNbrCounting, ChargeContribution);
	    }
	}
      for (int j = 0; j < CurrentNbrCounting; ++j)
	{
	  this->InteractionPerComponentIndex[Pos][position] = TmpCounting[j];
	  this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients[j];
	  ++position;
	}
      ++Pos;

    }
  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  delete[] TmpCounting;
  delete[] TmpCoefficients;

  delete TmpParticles;
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

RealVector& ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
																int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  double TmpCoefficient;
  double TmpCoef;
  double ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
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

  int Index;  
  double TmpInteraction;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	{
	  TmpCoef = vSource[i];
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if (TmpCoefficient != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
		      NbrElements = TmpParticles->AduAdd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
		    }	
		}
	    }
	    
	  if (this->OneBodyInteractionFactorsupup != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	      {
		if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
		  {
		    NbrElements = TmpParticles->AduAu(i, lz, TmpLeftIndices, TmpInteractionElements);
		    for (int j = 0; j < NbrElements; ++j)
		      vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
		  }
	      }
	  }
	      
	  if (this->OneBodyInteractionFactorsdowndown != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
	      {
		NbrElements = TmpParticles->AddAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		for (int j = 0; j < NbrElements; ++j)
		  vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
		}
	    }
	  }
	  
	  if (this->ChargingEnergy != 0.0)
	  {
	    int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	    if (TmpTotalNbrParticles != this->AverageNumberParticles)
	    {
	      ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
	      ChargeContribution *= ChargeContribution;
	      ChargeContribution *= (TmpCoef * this->ChargingEnergy);
	      vDestination[i] += ChargeContribution;
	    }
	  }
      }
      }

  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  delete TmpParticles;
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

ComplexVector& ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
																   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  Complex Coefficient;
  double TmpCoefficient;
  Complex TmpCoef;
  double ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
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

  int Index;  
  double TmpInteraction;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	  {
	    TmpCoef = vSource[i];
	    if (this->OneBodyInteractionFactorsPairing != 0)
	      {
		for (int lz = 0; lz <= this->LzMax; ++lz)
		  {
		    TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		    if (TmpCoefficient != 0.0)
		      {
			NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
			for (int j = 0; j < NbrElements; ++j)
			  vDestination[TmpLeftIndices[j]] += TmpCoef * (TmpInteractionElements[j] * TmpCoefficient);
			NbrElements = TmpParticles->AduAdd(i, lz, TmpLeftIndices, TmpInteractionElements);
			for (int j = 0; j < NbrElements; ++j)
			  vDestination[TmpLeftIndices[j]] += TmpCoef * (TmpInteractionElements[j] * TmpCoefficient);
		      }	
		  }
	      }
	    
	    if (this->OneBodyInteractionFactorsupup != 0)
	      {
		for (int lz = 0; lz <= this->LzMax; ++lz)
		  {
		    if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
		      {
			NbrElements = TmpParticles->AduAu(i, lz, TmpLeftIndices, TmpInteractionElements);
			for (int j = 0; j < NbrElements; ++j)
			  vDestination[TmpLeftIndices[j]] += TmpCoef * (TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
		      }
		  }
	      }
	    
	    if (this->OneBodyInteractionFactorsdowndown != 0)
	      {
		for (int lz = 0; lz <= this->LzMax; ++lz)
		  {
		    if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
		      {
			NbrElements = TmpParticles->AddAd(i, lz, TmpLeftIndices, TmpInteractionElements);
			for (int j = 0; j < NbrElements; ++j)
			  vDestination[TmpLeftIndices[j]] += TmpCoef * (TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
		      }
		  }
	      }
	    
	    if (this->ChargingEnergy != 0.0)
	      {
		int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
		if (TmpTotalNbrParticles != this->AverageNumberParticles)
		  {
		    ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
		    ChargeContribution *= ChargeContribution;
		    ChargeContribution *= this->ChargingEnergy;
		    vDestination[i] += TmpCoef * ChargeContribution;
		  }
	      }
	  }
      }
  
  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  delete TmpParticles;
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

RealVector& ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::HermitianLowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
																	 int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  double TmpCoefficient;
  double TmpCoef;
  double ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
  int j;
  double TmpSum;
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
      TmpSum = 0.0;
      for (j = 0; j < TmpNbrInteraction; ++j)
	{
	  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	  TmpSum += TmpCoefficientArray[j] * vSource[TmpIndexArray[j]];
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
	{
	  TmpSum = 0.0;
	  TmpCoef = vSource[i];
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if (TmpCoefficient != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
		      {
			vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
			TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * TmpCoefficient);
		      }
		    }	
		}
	    }
	    
	  if (this->OneBodyInteractionFactorsupup != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	      {
		if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
		  {
		    NbrElements = TmpParticles->AduAu(i, lz, TmpLeftIndices, TmpInteractionElements);
		    for (int j = 0; j < NbrElements; ++j)
		    {
		      if (TmpLeftIndices[j] <= i)
		      {
			vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
			if (TmpLeftIndices[j] < i)
			  TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);		      
		      }
		    }
		  }
	      }
	  }
	      
	  if (this->OneBodyInteractionFactorsdowndown != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
	      {
		NbrElements = TmpParticles->AddAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		for (int j = 0; j < NbrElements; ++j)
		{
		  if (TmpLeftIndices[j] <= i)
		  {
		    vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
		    if (TmpLeftIndices[j] < i)
		      TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
		  }
		}
		}
	    }
	  }
	  
	  if (this->ChargingEnergy != 0.0)
	  {
	    int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	    if (TmpTotalNbrParticles != this->AverageNumberParticles)
	    {
	      ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
	      ChargeContribution *= ChargeContribution;
	      ChargeContribution *= (TmpCoef * this->ChargingEnergy);
	      vDestination[i] += ChargeContribution;
	    }
	  }
	  vDestination[i] += TmpSum;
	}
      }

  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  delete TmpParticles;
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

ComplexVector& ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::HermitianLowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
																	 int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  Complex Coefficient;
  double TmpCoefficient;
  Complex TmpCoef;
  double ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
  int j;
  Complex TmpSum;
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
      TmpSum = 0.0;
      for (j = 0; j < TmpNbrInteraction; ++j)
	{
	  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	  TmpSum += TmpCoefficientArray[j] * vSource[TmpIndexArray[j]];
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
	{
	  TmpSum = 0.0;
	  TmpCoef = vSource[i];
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if (TmpCoefficient != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
		      {
			vDestination[TmpLeftIndices[j]] += TmpCoef * (TmpInteractionElements[j] * TmpCoefficient);
			TmpSum += vSource[TmpLeftIndices[j]] * (TmpInteractionElements[j] * TmpCoefficient);
		      }
		    }	
		}
	    }
	    
	  if (this->OneBodyInteractionFactorsupup != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
		    {
		      NbrElements = TmpParticles->AduAu(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  if (TmpLeftIndices[j] <= i)
			    {
			      vDestination[TmpLeftIndices[j]] += TmpCoef * (TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
			      if (TmpLeftIndices[j] < i)
				TmpSum += vSource[TmpLeftIndices[j]] * (TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);		      
			    }
			}
		    }
		}
	    }
	      
	  if (this->OneBodyInteractionFactorsdowndown != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
		    {
		      NbrElements = TmpParticles->AddAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  if (TmpLeftIndices[j] <= i)
			    {
			      vDestination[TmpLeftIndices[j]] += TmpCoef * (TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
			      if (TmpLeftIndices[j] < i)
				TmpSum += vSource[TmpLeftIndices[j]] * (TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
			    }
			}
		    }
		}
	    }
	  
	  if (this->ChargingEnergy != 0.0)
	    {
	      int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	      if (TmpTotalNbrParticles != this->AverageNumberParticles)
		{
		  ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
		  ChargeContribution *= ChargeContribution;
		  ChargeContribution *= this->ChargingEnergy;
		  vDestination[i] += TmpCoef * ChargeContribution;
		}
	    }
	  vDestination[i] += TmpSum;
	}
      }

  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  delete TmpParticles;
  return vDestination;
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

RealVector* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double* Coefficient2 = new double [nbrVectors];
  
  double Coefficient;
  double TmpCoefficient;
  double ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  int Pos2;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
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
	    {
	      vDestinations[k][Pos2] += Coefficient  * Coefficient2[k];
	    }
	}
      l += this->FastMultiplicationStep;
      ++Pos;
    }

  int Index;  
  double TmpInteraction;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	{
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if (TmpCoefficient != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][TmpLeftIndices[j]] += (vSources[l][i] * TmpInteractionElements[j] * TmpCoefficient);
		      NbrElements = TmpParticles->AduAdd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][TmpLeftIndices[j]] += (vSources[l][i] * TmpInteractionElements[j] * TmpCoefficient);
		    }	
		}
	    }
	    
	  if (this->OneBodyInteractionFactorsupup != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	      {
		if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
		  {
		    NbrElements = TmpParticles->AduAu(i, lz, TmpLeftIndices, TmpInteractionElements);
		    for (int j = 0; j < NbrElements; ++j)
		      for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][TmpLeftIndices[j]] += (vSources[l][i] * TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
		  }
	      }
	  }
	      
	  if (this->OneBodyInteractionFactorsdowndown != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
	      {
		NbrElements = TmpParticles->AddAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		for (int j = 0; j < NbrElements; ++j)
		  for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][TmpLeftIndices[j]] += (vSources[l][i] * TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
		}
	    }
	  }
	  
	  if (this->ChargingEnergy != 0.0)
	  {
	    int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	    if (TmpTotalNbrParticles != this->AverageNumberParticles)
	    {
	      ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
	      ChargeContribution *= ChargeContribution;
	      ChargeContribution *= (this->ChargingEnergy);
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] += ChargeContribution * vSources[l][i];
	    }
	  }
      }
      }

  delete[] Coefficient2;
  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  delete TmpParticles;
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

ComplexVector* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
																	   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  Complex* Coefficient2 = new Complex [nbrVectors];
  
  double Coefficient;
  double TmpCoefficient;
  double ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  int Pos2;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
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
	    {
	      vDestinations[k][Pos2] += Coefficient  * Coefficient2[k];
	    }
	}
      l += this->FastMultiplicationStep;
      ++Pos;
    }

  int Index;  
  double TmpInteraction;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	{
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if (TmpCoefficient != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][TmpLeftIndices[j]] += vSources[l][i] * (TmpInteractionElements[j] * TmpCoefficient);
		      NbrElements = TmpParticles->AduAdd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][TmpLeftIndices[j]] += vSources[l][i] * (TmpInteractionElements[j] * TmpCoefficient);
		    }	
		}
	    }
	    
	  if (this->OneBodyInteractionFactorsupup != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
		    {
		      NbrElements = TmpParticles->AduAu(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][TmpLeftIndices[j]] += vSources[l][i] * (TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
		    }
		}
	    }
	  
	  if (this->OneBodyInteractionFactorsdowndown != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
		    {
		      NbrElements = TmpParticles->AddAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][TmpLeftIndices[j]] += vSources[l][i] * (TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
		    }
		}
	    }
	  
	  if (this->ChargingEnergy != 0.0)
	    {
	      int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	      if (TmpTotalNbrParticles != this->AverageNumberParticles)
		{
		  ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
		  ChargeContribution *= ChargeContribution;
		  ChargeContribution *= (this->ChargingEnergy);
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][i] += ChargeContribution * vSources[l][i];
		}
	    }
	}
      }

  delete[] Coefficient2;
  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  delete TmpParticles;
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

RealVector* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
																		 int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double* Coefficient2 = new double [nbrVectors];
  
  double Coefficient;
  double TmpCoefficient;
  double ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  int Pos2;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
  int j;
  double* TmpSum = new double[nbrVectors];
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
	      TmpSum[k] += Coefficient * vSources[k][Pos2];
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

  int Index;  
  double TmpInteraction;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	{
	  for (int l = 0; l < nbrVectors; ++l)
	    TmpSum[l] = 0.0;
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if (TmpCoefficient != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int l = 0; l < nbrVectors; ++l)
			{
			  vDestinations[l][TmpLeftIndices[j]] += (vSources[l][i] * TmpInteractionElements[j] * TmpCoefficient);
			  TmpSum[l] += (vSources[l][TmpLeftIndices[j]] * TmpInteractionElements[j] * TmpCoefficient);
			}
		    }	
		}
	    }
	    
	  if (this->OneBodyInteractionFactorsupup != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	      {
		if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
		  {
		    NbrElements = TmpParticles->AduAu(i, lz, TmpLeftIndices, TmpInteractionElements);
		    for (int j = 0; j < NbrElements; ++j)
		    {
		      if (TmpLeftIndices[j] <= i)
		      {
			for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][TmpLeftIndices[j]] += (vSources[l][i] * TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
			if (TmpLeftIndices[j] < i)
			  for (int l = 0; l < nbrVectors; ++l)
			    TmpSum[l] += (vSources[l][TmpLeftIndices[j]] * TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
		      }
		    }
		  }
	      }
	  }
	      
	  if (this->OneBodyInteractionFactorsdowndown != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
	      {
		NbrElements = TmpParticles->AddAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		for (int j = 0; j < NbrElements; ++j)
		{
		  if (TmpLeftIndices[j] <= i)
		  {
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][TmpLeftIndices[j]] += (vSources[l][i] * TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
		    if (TmpLeftIndices[j] < i)
		      for (int l = 0; l < nbrVectors; ++l)
			TmpSum[l] += (vSources[l][TmpLeftIndices[j]] * TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
		  }
		}
	      }
	    }
	  }
	  
	  if (this->ChargingEnergy != 0.0)
	  {
	    int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	    if (TmpTotalNbrParticles != this->AverageNumberParticles)
	    {
	      ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
	      ChargeContribution *= ChargeContribution;
	      ChargeContribution *= (this->ChargingEnergy);
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] += ChargeContribution * vSources[l][i];
	    }
	  }
	  for (int l = 0; l < nbrVectors; ++l)
	    vDestinations[l][i] += TmpSum[l];
      }
      }

  delete[] TmpSum;
  delete[] Coefficient2;
  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  delete TmpParticles;
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

ComplexVector* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
																		    int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  Complex* Coefficient2 = new Complex [nbrVectors];
  
  double Coefficient;
  double TmpCoefficient;
  double ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  int Pos2;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
  int j;
  Complex* TmpSum = new Complex[nbrVectors];
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
	      TmpSum[k] += Coefficient * vSources[k][Pos2];
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

  int Index;  
  double TmpInteraction;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	  {
	    for (int l = 0; l < nbrVectors; ++l)
	      TmpSum[l] = 0.0;
	    if (this->OneBodyInteractionFactorsPairing != 0)
	      {
		for (int lz = 0; lz <= this->LzMax; ++lz)
		  {
		    TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		    if (TmpCoefficient != 0.0)
		      {
			NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
			for (int j = 0; j < NbrElements; ++j)
			  for (int l = 0; l < nbrVectors; ++l)
			    {
			      vDestinations[l][TmpLeftIndices[j]] += vSources[l][i] * (TmpInteractionElements[j] * TmpCoefficient);
			      TmpSum[l] += vSources[l][TmpLeftIndices[j]] * (TmpInteractionElements[j] * TmpCoefficient);
			    }
		      }	
		  }
	      }
	    
	    if (this->OneBodyInteractionFactorsupup != 0)
	      {
		for (int lz = 0; lz <= this->LzMax; ++lz)
		  {
		    if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
		      {
			NbrElements = TmpParticles->AduAu(i, lz, TmpLeftIndices, TmpInteractionElements);
			for (int j = 0; j < NbrElements; ++j)
			  {
			    if (TmpLeftIndices[j] <= i)
			      {
				for (int l = 0; l < nbrVectors; ++l)
				  vDestinations[l][TmpLeftIndices[j]] += vSources[l][i] * (TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
				if (TmpLeftIndices[j] < i)
				  for (int l = 0; l < nbrVectors; ++l)
				    TmpSum[l] += vSources[l][TmpLeftIndices[j]] * (TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
			      }
			  }
		      }
		  }
	      }
	    
	    if (this->OneBodyInteractionFactorsdowndown != 0)
	      {
		for (int lz = 0; lz <= this->LzMax; ++lz)
		  {
		    if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
		      {
			NbrElements = TmpParticles->AddAd(i, lz, TmpLeftIndices, TmpInteractionElements);
			for (int j = 0; j < NbrElements; ++j)
			  {
			    if (TmpLeftIndices[j] <= i)
			      {
				for (int l = 0; l < nbrVectors; ++l)
				  vDestinations[l][TmpLeftIndices[j]] += vSources[l][i] * (TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
				if (TmpLeftIndices[j] < i)
				  for (int l = 0; l < nbrVectors; ++l)
				    TmpSum[l] += vSources[l][TmpLeftIndices[j]] * (TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
			      }
			  }
		      }
		  }
	      }
	    
	    if (this->ChargingEnergy != 0.0)
	      {
		int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
		if (TmpTotalNbrParticles != this->AverageNumberParticles)
		  {
		    ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
		    ChargeContribution *= ChargeContribution;
		    ChargeContribution *= (this->ChargingEnergy);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][i] += ChargeContribution * vSources[l][i];
		  }
	      }
	    for (int l = 0; l < nbrVectors; ++l)
	      vDestinations[l][i] += TmpSum[l];
	  }
      }

  delete[] TmpSum;
  delete[] Coefficient2;
  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  delete TmpParticles;
  return vDestinations;
}

// ask if Hamiltonian implements hermitian symmetry operations
//

bool ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::IsHermitian()
{
  return this->HermitianSymmetryFlag;
}

// get the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
//
// nbrThreads = number of threads requested
// segmentIndices = array returning the reference to an array of the first index of each of the segments
// return value = true if no error occured

bool ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::GetLoadBalancing(int nbrTasks, long* &segmentIndices)
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
	  if (LoadBalancingArray!=0)
	    delete [] LoadBalancingArray;
	  long *SegmentSize = new long[nbrTasks];
	  this->LoadBalancingArray = new long[nbrTasks+1];
	  this->NbrBalancedTasks = nbrTasks;
	  long TmpNbrElement = 0;
	  for (int i=0; i < ReducedSpaceDimension; ++i)
	    TmpNbrElement += this->NbrInteractionPerComponent[i];
	  long TmpNbrPerSegment = TmpNbrElement / nbrTasks;
	  TmpNbrElement = 0;
	  int Pos=0;
	  this->LoadBalancingArray[0] = MinIndex;
	  for (int i = 0; i < ReducedSpaceDimension; ++i)
	    {
	      TmpNbrElement += this->NbrInteractionPerComponent[i];
	      if (TmpNbrElement > TmpNbrPerSegment)
		{
		  SegmentSize[Pos] = TmpNbrElement;
		  this->LoadBalancingArray[Pos + 1]= MinIndex + (i * this->FastMultiplicationStep);
		  TmpNbrElement = 0;
		  ++Pos;
		}
	    }
	  while (Pos < (nbrTasks - 1))
	    {
	      LoadBalancingArray[Pos + 1] = MaxIndex + 1;
	      SegmentSize[Pos] = 0;
	      ++Pos;
	    }
	  LoadBalancingArray[nbrTasks] = MaxIndex + 1;
	  SegmentSize[nbrTasks - 1] = TmpNbrElement;
	  
	  cout << "LoadBalancingArray=[ ("<< LoadBalancingArray[1] - LoadBalancingArray[0] <<", "<<SegmentSize[0]<<")";
	  for (int i = 1; i < nbrTasks; ++i)
	    cout <<" ("<< LoadBalancingArray[i+1] - LoadBalancingArray[i] << ", " << SegmentSize[i] << ")";
	  cout << "]"<< endl;
	  delete [] SegmentSize;
	}
    }
  else
    {
      return false;
    }
  segmentIndices = this->LoadBalancingArray;
  return true;
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::FastMultiplicationMemory(long allowedMemory)
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
  long TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(double*)) * EffectiveHilbertSpaceDimension;
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(double)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(double)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if (EffectiveHilbertSpaceDimension != (ReducedSpaceDimension * this->FastMultiplicationStep))
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
	  if (this->DiskStorageFlag == false)
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
	}
      else
	if (this->DiskStorageFlag == false)
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
