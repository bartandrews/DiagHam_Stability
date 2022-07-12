////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to Laughlin qh on a sphere with      //
//        SU(2) spin with opposite magnetic field for each species,           //
//                      pairing and no momentum conservation                  //
//                                                                            //
//                        last modification : 26/08/2016                      //
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


#include "Hamiltonian/ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta.h"
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
//

ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta()
{
  this->NbrSectorSums = 0;
  this->OneBodyInteractionFactors = 0;
  this->InteractionPerComponentIntegerCoefficient = 0;
  this->HermitianSymmetryFlag = false;
  this->LoadBalancingArray = 0;
}

// constructor from default data
//
// particles = Hilbert space associated to the system
// lzmax = maximum Lz value reached by a particle in the state
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyOffDiagonalPotentialUpUp = off-diagonal contribution of the one-body potential for particles with spin up, the first entry is the annihilation index, 
//                                   the second entry is the momentum tranfer
// onebodyOffDiagonalPotentialDownDown = off-diagonal contribution of the one-body potential for particles with spin down, the first entry is the annihilation index, 
//                                       the second entry is the momentum tranfer
// onebodyPotentialPairing =  one-body pairing term (sorted from component on the lowest Lz state to component on the highest Lz state), on site, symmetric spin up / spin down
// onebodyOffDiagonalPotentialPairing = off diagonal contribution of the one-body pairing term, the first entry is the index of the rightmost creation operator, 
//	    				  the second entry is the momentum tranfer
// chargingEnergy = factor in front of the charging energy (i.e 1/(2C))
// averageNumberParticles = average number of particles in the system
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta(QuasiholeOnSphereWithSpinAndPairing* particles, int lzmax, 
																					     int maxMomentumTransfer, double* onebodyPotentialUpUp, 
																					     double* onebodyPotentialDownDown, Complex** onebodyOffDiagonalPotentialUpUp, 
																					     Complex** onebodyOffDiagonalPotentialDownDown, 
																					     Complex* onebodyPotentialPairing, 
																					     Complex** onebodyOffDiagonalPotentialPairing, double chargingEnergy, 
																					     double averageNumberParticles, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->MaximumMomentumTransfer = maxMomentumTransfer;
  
  this->ChargingEnergy = chargingEnergy;
  this->AverageNumberParticles = averageNumberParticles;
  this->HamiltonianShift =  0.0;
  long MinIndex;
  long MaxIndex;
  
  this->HermitianSymmetryFlag = true;
    
  this->LoadBalancingArray = 0;
  this->NbrBalancedTasks  = 0;

  this->Architecture = architecture;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->Memory = memory;
  this->FastMultiplicationFlag = false;

  this->NbrSectorSums = 0;
  this->OneBodyInteractionFactors = 0;
  this->InteractionPerComponentIntegerCoefficient = 0;

  this->OneBodyInteractionFactorsupup = 0;
  if (onebodyPotentialUpUp != 0)
    {
      this->OneBodyInteractionFactorsupup = new double [this->LzMax + 1];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsupup[i] = onebodyPotentialUpUp[i];
    }
  
  this->OneBodyInteractionFactorsdowndown = 0;
  if (onebodyPotentialDownDown != 0)
    {
      this->OneBodyInteractionFactorsdowndown = new double [this->LzMax + 1];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsdowndown[i] = onebodyPotentialDownDown[i];
    }
  this->OneBodyOffDiagonalInteractionFactorsupup = 0;
  if (onebodyOffDiagonalPotentialUpUp != 0)
    {
      this->OneBodyOffDiagonalInteractionFactorsupup = new Complex*[this->LzMax + 1];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->OneBodyOffDiagonalInteractionFactorsupup[i] = new Complex[this->MaximumMomentumTransfer];
	  for (int j = 0; j < this->MaximumMomentumTransfer; ++j)
	    {
	      this->OneBodyOffDiagonalInteractionFactorsupup[i][j] = onebodyOffDiagonalPotentialUpUp[i][j];
	    }
	}
    }
  this->OneBodyOffDiagonalInteractionFactorsdowndown = 0;
  if (onebodyOffDiagonalPotentialDownDown != 0)
    {
      this->OneBodyOffDiagonalInteractionFactorsdowndown = new Complex*[this->LzMax + 1];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->OneBodyOffDiagonalInteractionFactorsdowndown[i] = new Complex[this->MaximumMomentumTransfer];
	  for (int j = 0; j < this->MaximumMomentumTransfer; ++j)
	    {
	      this->OneBodyOffDiagonalInteractionFactorsdowndown[i][j] = onebodyOffDiagonalPotentialDownDown[i][j];
	    }
	}
    }

  this->OneBodyInteractionFactorsPairing = 0;
  if (onebodyPotentialPairing != 0)
    {
      this->OneBodyInteractionFactorsPairing = new Complex [this->LzMax + 1];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->OneBodyInteractionFactorsPairing[i] = onebodyPotentialPairing[i];
	}
    }
  this->OneBodyOffDiagonalInteractionFactorsPairing = 0;
  if (onebodyOffDiagonalPotentialPairing != 0)
    {
      this->OneBodyOffDiagonalInteractionFactorsPairing = new Complex* [this->LzMax + 1];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->OneBodyOffDiagonalInteractionFactorsPairing[i] = new Complex[2 * this->MaximumMomentumTransfer + 1];
	  for (int j = 0; j < this->MaximumMomentumTransfer; ++j)
	    {
	      this->OneBodyOffDiagonalInteractionFactorsPairing[i][j] = onebodyOffDiagonalPotentialPairing[i][j];
	    }
	  for (int j = this->MaximumMomentumTransfer + 1; j <= (2 * this->MaximumMomentumTransfer); ++j)
	    {
	      this->OneBodyOffDiagonalInteractionFactorsPairing[i][j] = onebodyOffDiagonalPotentialPairing[i][j];
	    }
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
}

// destructor
//

ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::~ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta() 
{
  if (this->OneBodyInteractionFactorsupup != 0)
    delete[] this->OneBodyInteractionFactorsupup;
  if (this->OneBodyInteractionFactorsdowndown != 0)
    delete[] this->OneBodyInteractionFactorsdowndown;
  if (this->OneBodyInteractionFactorsPairing != 0)
    delete [] this->OneBodyInteractionFactorsPairing;
  if (this->OneBodyOffDiagonalInteractionFactorsupup != 0)
    {
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  delete[] this->OneBodyOffDiagonalInteractionFactorsupup[i];
	}
      delete[] this->OneBodyOffDiagonalInteractionFactorsupup;
    }
  if (this->OneBodyOffDiagonalInteractionFactorsdowndown != 0)
    {
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  delete[] this->OneBodyOffDiagonalInteractionFactorsdowndown[i];
	}
      delete[] this->OneBodyOffDiagonalInteractionFactorsdowndown;
    }
  if (this->OneBodyOffDiagonalInteractionFactorsPairing != 0)
    {
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  delete[] this->OneBodyOffDiagonalInteractionFactorsPairing[i];
	}
      delete[] this->OneBodyOffDiagonalInteractionFactorsPairing;
    }  
  
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (QuasiholeOnSphereWithSpinAndPairing*) hilbertSpace;
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::ShiftHamiltonian (double shift)
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

ComplexVector& ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
      Complex ChargeContribution;
      Complex TmpCoef;
      
      for (int i = 0; i < nbrComponent; ++i)
	{
	  TmpCoef = vSource[index];
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if ((TmpCoefficient.Re != 0.0) || (TmpCoefficient.Im != 0.0))
		    {
		      NbrElements = TmpParticles->AuAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
		      TmpCoefficient.Im *= -1.0;
		      NbrElements = TmpParticles->AduAdd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
		    }	
		}
	    }

	  if (this->OneBodyOffDiagonalInteractionFactorsPairing != 0)
	    {
	      for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
		{
		  for (int lz = 0; lz < (this->LzMax - k); ++lz)
		    {
		      TmpCoefficient = this->OneBodyOffDiagonalInteractionFactorsPairing[lz][k];
		      if ((TmpCoefficient.Re != 0.0) || (TmpCoefficient.Im != 0.0))
			{
			  NbrElements = TmpParticles->AuAd(index, lz + k + 1, lz, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
			  NbrElements = TmpParticles->AuAd(index, lz, lz + k + 1, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
			  TmpCoefficient.Im *= -1.0;
			  NbrElements = TmpParticles->AduAdd(index, lz + k + 1, lz, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
			  NbrElements = TmpParticles->AduAdd(index, lz, lz + k + 1, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
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
	  
	  if (this->OneBodyOffDiagonalInteractionFactorsupup != 0)
	    {
	      for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
		{
		  for (int lz = 0; lz < (this->LzMax - k); ++lz)
		    {
		      if (this->OneBodyOffDiagonalInteractionFactorsupup[lz][k] != 0.0)
			{
			  NbrElements = TmpParticles->AduAu(index, lz + k + 1, lz, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyOffDiagonalInteractionFactorsupup[lz][k]);
			  NbrElements = TmpParticles->AduAu(index, lz, lz + k + 1, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * Conj(this->OneBodyOffDiagonalInteractionFactorsupup[lz][k]));
			}
		    }
		}
	    }
	  if (this->OneBodyOffDiagonalInteractionFactorsdowndown != 0)
	    {
	      for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
		{
		  for (int lz = 0; lz < (this->LzMax - k); ++lz)
		    {
		      if (this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k] != 0.0)
			{
			  NbrElements = TmpParticles->AddAd(index, lz + k + 1, lz, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k]);
			  NbrElements = TmpParticles->AddAd(index, lz, lz + k + 1, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * Conj(this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k]));
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
		  ChargeContribution.Re *= ChargeContribution.Re;
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
	  Complex* TmpCoefficientArray; 
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
	  this->LowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
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

ComplexVector& ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
      Complex ChargeContribution;
      Complex TmpCoef;
      Complex TmpSum;
      Complex TmpCoefficient;
      for (int i = 0; i < nbrComponent; ++i)
	{
	  TmpSum = 0.0;
	  TmpCoef = vSource[index];
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if ((TmpCoefficient.Re != 0.0) || (TmpCoefficient.Im != 0.0))
		    {
		      NbrElements = TmpParticles->AuAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
			  TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * Conj(TmpCoefficient));
			}
		    }	
		}
	    }
	  
	  if (this->OneBodyOffDiagonalInteractionFactorsPairing != 0)
	    {
	      for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
		{
		  for (int lz = 0; lz < (this->LzMax - k); ++lz)
		    {
		      TmpCoefficient = this->OneBodyOffDiagonalInteractionFactorsPairing[lz][k];
		      if ((TmpCoefficient.Re != 0.0) || (TmpCoefficient.Im != 0.0))
			{
			  NbrElements = TmpParticles->AuAd(index, lz + k + 1, lz, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    {
			      vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
			      TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * Conj(TmpCoefficient));
			    }
			  NbrElements = TmpParticles->AuAd(index, lz, lz + k + 1, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    {
			      vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
			      TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * Conj(TmpCoefficient));
			    }	
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
	  
	  if (this->OneBodyOffDiagonalInteractionFactorsupup != 0)
	    {
	      for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
		{
		  for (int lz = 0; lz < (this->LzMax - k); ++lz)
		    {
		      if (this->OneBodyOffDiagonalInteractionFactorsupup[lz][k] != 0.0)
			{
			  NbrElements = TmpParticles->AduAu(index, lz + k + 1, lz, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    {
			      if (TmpLeftIndices[j] <= index)
				{
				  vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyOffDiagonalInteractionFactorsupup[lz][k]);
				  if (TmpLeftIndices[j] < index)
				    TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * this->OneBodyOffDiagonalInteractionFactorsupup[lz][k]);
				}
			    }
			  NbrElements = TmpParticles->AduAu(index, lz, lz + k + 1, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    {
			      if (TmpLeftIndices[j] <= index)
				{
				  vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * Conj(this->OneBodyOffDiagonalInteractionFactorsupup[lz][k]));
				  if (TmpLeftIndices[j] < index)
				    TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * Conj(this->OneBodyOffDiagonalInteractionFactorsupup[lz][k]));
				}
			    }
			}
		    }
		}
	    }

	  if (this->OneBodyOffDiagonalInteractionFactorsdowndown != 0)
	    {
	      for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
		{
		  for (int lz = 0; lz < (this->LzMax - k); ++lz)
		    {
		      if (this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k] != 0.0)
			{
			  NbrElements = TmpParticles->AddAd(index, lz + k + 1, lz, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    {
			      if (TmpLeftIndices[j] <= index)
				{
				  vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k]);
				  if (TmpLeftIndices[j] < index)
				    TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k]);
				}
			    }
			  NbrElements = TmpParticles->AddAd(index, lz, lz + k + 1, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    {
			      if (TmpLeftIndices[j] <= index)
				{
				  vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * Conj(this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k]));
				  if (TmpLeftIndices[j] < index)
				    TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * Conj(this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k]));
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
		  ChargeContribution.Re *= ChargeContribution.Re;
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
		  TmpSum += TmpCoefficientArray[j] * vSource[TmpIndexArray[j]];
		  vDestination[TmpIndexArray[j]] +=  Conj(TmpCoefficientArray[j]) * Coefficient;
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

ComplexVector* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, 
																  int nbrVectors, int firstComponent, int nbrComponent)
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
      Complex ChargeContribution;
      Complex TmpOneBodyInteraction;
      for (int i = 0; i < nbrComponent; ++i)
	{
	  for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsPairing != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsPairing[lz];    
		  if ((TmpOneBodyInteraction.Re != 0.0) || (TmpOneBodyInteraction.Im != 0.0))
		    {
		      NbrElements = TmpParticles->AuAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			  vDestinations[k][TmpLeftIndices[j]] += (vSources[k][index] * TmpInteractionElements[j] * TmpOneBodyInteraction);
		      TmpOneBodyInteraction.Im *= -1.0;
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
		  ChargeContribution.Re *= ChargeContribution.Re;
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
	  Complex* TmpCoefficientArray; 
	  int j;
	  int Pos;
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


// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
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
      Complex TmpOneBodyInteraction;
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
		  if ((TmpOneBodyInteraction.Re != 0.0) || (TmpOneBodyInteraction.Im != 0.0))
		    {
		      NbrElements = TmpParticles->AuAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  for (int k = 0; k < nbrVectors; ++k)
			    {
			      vDestinations[k][TmpLeftIndices[j]] += (vSources[k][index] * TmpInteractionElements[j] * TmpOneBodyInteraction);
			      TmpSum[k] += (vSources[k][TmpLeftIndices[j]] * TmpInteractionElements[j] * Conj(TmpOneBodyInteraction));
			    }
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

AbstractHilbertSpace* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::GetHilbertSpace ()
{
  return this->Particles;
}


// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that has to be precalcualted
// return value = number of non-zero matrix element

long ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::PartialFastMultiplicationMemory(int firstComponent, int nbrComponent)
{
  int NbrElements;
  Complex TmpCoefficient;
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
	      if ((TmpCoefficient.Re != 0.0) || (TmpCoefficient.Im != 0.0))
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
      if (this->OneBodyOffDiagonalInteractionFactorsPairing != 0)
	{
	  for (int k = 1; k <= this->MaximumMomentumTransfer; ++k)
	    {
	      for (int lz = 0; lz <= (this->LzMax - k); ++lz)
		{
		  TmpCoefficient = this->OneBodyOffDiagonalInteractionFactorsPairing[lz][this->MaximumMomentumTransfer + k];
		  if ((TmpCoefficient.Re != 0.0) || (TmpCoefficient.Im != 0.0))
		    {
		      NbrElements = TmpParticles->AuAd(i, lz, lz + k, TmpLeftIndices, TmpInteractionElements);
		      Memory += NbrElements;
		      TmpTotal += NbrElements;
		      this->NbrInteractionPerComponent[i - this->PrecalculationShift] += NbrElements;
		      if (this->HermitianSymmetryFlag == false)
			{
			  NbrElements = TmpParticles->AduAdd(i, lz, lz + k, TmpLeftIndices, TmpInteractionElements);
			  Memory += NbrElements;
			  TmpTotal += NbrElements;
			  this->NbrInteractionPerComponent[i - this->PrecalculationShift] += NbrElements;
			}
		    }
		}
	    }
	  for (int k = 1; k <= this->MaximumMomentumTransfer; ++k)
	    {
	      for (int lz = k; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyOffDiagonalInteractionFactorsPairing[lz][this->MaximumMomentumTransfer - k];
		  if ((TmpCoefficient.Re != 0.0) || (TmpCoefficient.Im != 0.0))
		    {
		      NbrElements = TmpParticles->AuAd(i, lz, lz - k, TmpLeftIndices, TmpInteractionElements);
		      Memory += NbrElements;
		      TmpTotal += NbrElements;
		      this->NbrInteractionPerComponent[i - this->PrecalculationShift] += NbrElements;
		      if (this->HermitianSymmetryFlag == false)
			{
			  NbrElements = TmpParticles->AduAdd(i, lz, lz - k, TmpLeftIndices, TmpInteractionElements);
			  Memory += NbrElements;
			  TmpTotal += NbrElements;
			  this->NbrInteractionPerComponent[i - this->PrecalculationShift] += NbrElements;
			}
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
      
      if (this->OneBodyOffDiagonalInteractionFactorsupup != 0)
	{
	  for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
	    {
	      for (int lz = 0; lz < (this->LzMax - k); ++lz)
		{
		  if (this->OneBodyOffDiagonalInteractionFactorsupup[lz][k] != 0.0)
		    {
		      NbrElements = TmpParticles->AduAu(i, lz + k + 1, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  if ((this->HermitianSymmetryFlag == false) || (TmpLeftIndices[j] <= i))
			    CurrentNbrCounting += SearchInSortedArrayAndInsert<int>(TmpLeftIndices[j], TmpCounting, CurrentNbrCounting);
			}
		      NbrElements = TmpParticles->AduAu(i, lz, lz + k + 1, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  if ((this->HermitianSymmetryFlag == false) || (TmpLeftIndices[j] <= i))
			    CurrentNbrCounting += SearchInSortedArrayAndInsert<int>(TmpLeftIndices[j], TmpCounting, CurrentNbrCounting);
			}
		    }
		}
	    }
	}
      
      if (this->OneBodyOffDiagonalInteractionFactorsdowndown != 0)
	{
	  for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
	    {
	      for (int lz = 0; lz < (this->LzMax - k); ++lz)
		{
		  if (this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k] != 0.0)
		    {
		      NbrElements = TmpParticles->AddAd(i, lz + k + 1, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  if ((this->HermitianSymmetryFlag == false) || (TmpLeftIndices[j] <= i))
			    CurrentNbrCounting += SearchInSortedArrayAndInsert<int>(TmpLeftIndices[j], TmpCounting, CurrentNbrCounting);
			}
		      NbrElements = TmpParticles->AddAd(i, lz, lz + k + 1, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  if ((this->HermitianSymmetryFlag == false) || (TmpLeftIndices[j] <= i))
			    CurrentNbrCounting += SearchInSortedArrayAndInsert<int>(TmpLeftIndices[j], TmpCounting, CurrentNbrCounting);
			}
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

void ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{  
  Complex TmpCoefficient;
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
  Complex* TmpCoefficients = new Complex [TmpParticles->GetHilbertSpaceDimension()];
  int CurrentNbrCounting = 0;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      for (int k = 0; k < CurrentNbrCounting; ++k)
	TmpCounting[k] = 0;
      CurrentNbrCounting = 0;
      int position = 0;
      if (this->OneBodyInteractionFactorsPairing != 0)
	{
	  for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
	      if ((TmpCoefficient.Re != 0.0) || (TmpCoefficient.Im != 0.0))
		{
		  NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		  for (int j = 0; j < NbrElements; ++j)
		    {
		      TmpCounting[CurrentNbrCounting] = TmpLeftIndices[j];
		      TmpCoefficients[CurrentNbrCounting] = TmpCoefficient * TmpInteractionElements[j];
		      ++CurrentNbrCounting;
		    }
		  if (this->HermitianSymmetryFlag == false)
		    {
		      NbrElements = TmpParticles->AduAdd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  TmpCounting[CurrentNbrCounting] = TmpLeftIndices[j];
			  TmpCoefficients[CurrentNbrCounting] = TmpCoefficient * Conj(TmpInteractionElements[j]);
			  ++CurrentNbrCounting;
			}
		    }
		}
	    }
	}
      if (this->OneBodyOffDiagonalInteractionFactorsPairing != 0)
	{
	  for (int k = 1; k <= this->MaximumMomentumTransfer; ++k)
	    {
	      for (int lz = 0; lz <= (this->LzMax - k); ++lz)
		{
		  TmpCoefficient = this->OneBodyOffDiagonalInteractionFactorsPairing[lz][this->MaximumMomentumTransfer + k];
		  if ((TmpCoefficient.Re != 0.0) || (TmpCoefficient.Im != 0.0))
		    {
		      NbrElements = TmpParticles->AuAd(i, lz, lz + k, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  TmpCounting[CurrentNbrCounting] = TmpLeftIndices[j];
			  TmpCoefficients[CurrentNbrCounting] = TmpCoefficient * TmpInteractionElements[j];
			  ++CurrentNbrCounting;
			}
		      if (this->HermitianSymmetryFlag == false)
			{
			  TmpCoefficient.Im *= -1.0;
			  NbrElements = TmpParticles->AduAdd(i, lz, lz + k, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    {
			      TmpCounting[CurrentNbrCounting] = TmpLeftIndices[j];
			      TmpCoefficients[CurrentNbrCounting] = TmpCoefficient * Conj(TmpInteractionElements[j]);
			      ++CurrentNbrCounting;
			    }
			}
		    }
		}
	    }
	  for (int k = 1; k <= this->MaximumMomentumTransfer; ++k)
	    {
	      for (int lz = k; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyOffDiagonalInteractionFactorsPairing[lz][this->MaximumMomentumTransfer - k];
		  if ((TmpCoefficient.Re != 0.0) || (TmpCoefficient.Im != 0.0))
		    {
		      NbrElements = TmpParticles->AuAd(i, lz, lz - k, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  TmpCounting[CurrentNbrCounting] = TmpLeftIndices[j];
			  TmpCoefficients[CurrentNbrCounting] = TmpCoefficient * TmpInteractionElements[j];
			  ++CurrentNbrCounting;
			}

		      if (this->HermitianSymmetryFlag == false)
			{
			  TmpCoefficient.Im *= -1.0;
			  NbrElements = TmpParticles->AduAdd(i, lz, lz - k, TmpLeftIndices, TmpInteractionElements);
			  for (int j = 0; j < NbrElements; ++j)
			    {
			      TmpCounting[CurrentNbrCounting] = TmpLeftIndices[j];
			      TmpCoefficients[CurrentNbrCounting] = TmpCoefficient * Conj(TmpInteractionElements[j]);
			      ++CurrentNbrCounting;
			    }
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
		  NbrElements = TmpParticles->AduAu(i, lz, TmpLeftIndices, TmpInteractionElements);
		  for (int j = 0; j < NbrElements; ++j)
		    {
		      if (this->HermitianSymmetryFlag == false)
			CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyInteractionFactorsupup[lz] * TmpInteractionElements[j]);
		      else
			{
			  if (TmpLeftIndices[j] < i)
			    CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyInteractionFactorsupup[lz] * TmpInteractionElements[j]);
			  if (TmpLeftIndices[j] == i)
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
		  NbrElements = TmpParticles->AddAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		  for (int j = 0; j < NbrElements; ++j)
		    {
		      if (this->HermitianSymmetryFlag == false)
			CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyInteractionFactorsdowndown[lz] * TmpInteractionElements[j]);
		      else
			{
			  if (TmpLeftIndices[j] < i)
			    CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyInteractionFactorsdowndown[lz] * TmpInteractionElements[j]);
			  if (TmpLeftIndices[j] == i)
			    CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyInteractionFactorsdowndown[lz] * TmpInteractionElements[j] * 0.5);
			}
		    }
		}
	    }
	}

      if (this->OneBodyOffDiagonalInteractionFactorsupup != 0)
	{
	  for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
	    {
	      for (int lz = 0; lz < (this->LzMax - k); ++lz)
		{
		  if (this->OneBodyOffDiagonalInteractionFactorsupup[lz][k] != 0.0)
		    {
		      NbrElements = TmpParticles->AduAu(i, lz + k + 1, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  if (this->HermitianSymmetryFlag == false)
			    {
			      CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyOffDiagonalInteractionFactorsupup[lz][k] * TmpInteractionElements[j]);
			    }
			  else
			    {
			      if (TmpLeftIndices[j] < i)
				CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyOffDiagonalInteractionFactorsupup[lz][k] * TmpInteractionElements[j]);
			      if (TmpLeftIndices[j] == i)
				CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyOffDiagonalInteractionFactorsupup[lz][k] * TmpInteractionElements[j] * 0.5);
			      
			    }
			}
		      NbrElements = TmpParticles->AduAu(i, lz, lz + k + 1, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  if (this->HermitianSymmetryFlag == false)
			    {
			      CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, Conj(this->OneBodyOffDiagonalInteractionFactorsupup[lz][k]) * TmpInteractionElements[j]);
			    }
			  else
			    {
			      if (TmpLeftIndices[j] < i)
				CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, Conj(this->OneBodyOffDiagonalInteractionFactorsupup[lz][k]) * TmpInteractionElements[j]);
			      if (TmpLeftIndices[j] == i)
				CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, Conj(this->OneBodyOffDiagonalInteractionFactorsupup[lz][k]) * TmpInteractionElements[j] * 0.5);
			      
			    }
			}
		    }
		}
	    }
	}
      
      if (this->OneBodyOffDiagonalInteractionFactorsdowndown != 0)
	{
	  for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
	    {
	      for (int lz = 0; lz < (this->LzMax - k); ++lz)
		{
		  if (this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k] != 0.0)
		    {
		      NbrElements = TmpParticles->AddAd(i, lz + k + 1, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  if (this->HermitianSymmetryFlag == false)
			    {
			      CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k] * TmpInteractionElements[j]);
			    }
			  else
			    {
			      if (TmpLeftIndices[j] < i)
				CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k] * TmpInteractionElements[j]);
			      if (TmpLeftIndices[j] == i)
				CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k] * TmpInteractionElements[j] * 0.5);
			      
			    }
			}
		      NbrElements = TmpParticles->AddAd(i, lz, lz + k + 1, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			{
			  if (this->HermitianSymmetryFlag == false)
			    {
			      CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, Conj(this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k]) * TmpInteractionElements[j]);
			    }
			  else
			    {
			      if (TmpLeftIndices[j] < i)
				CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, Conj(this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k]) * TmpInteractionElements[j]);
			      if (TmpLeftIndices[j] == i)
				CurrentNbrCounting += SearchInArrayAndSetWeight<int>(TmpLeftIndices[j], TmpCounting, TmpCoefficients, CurrentNbrCounting, Conj(this->OneBodyOffDiagonalInteractionFactorsdowndown[lz][k]) * TmpInteractionElements[j] * 0.5);
			      
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
	      if (this->HermitianSymmetryFlag)
		ChargeContribution *= 0.5;
	      CurrentNbrCounting += SearchInArrayAndSetWeight<int>(i, TmpCounting, TmpCoefficients, CurrentNbrCounting, ChargeContribution);
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

ComplexVector& ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::LowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
																	     int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  Complex Coefficient;
  Complex TmpCoefficient;
  Complex TmpCoef;
  Complex ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  int* TmpIndexArray;
  Complex* TmpCoefficientArray; 
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
		    if ((TmpCoefficient.Re != 0.0) || (TmpCoefficient.Im != 0.0))
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
		    ChargeContribution.Re *= ChargeContribution.Re;
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

ComplexVector& ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::HermitianLowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
																		      int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  Complex Coefficient;
  Complex TmpCoefficient;
  Complex TmpCoef;
  Complex ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  int* TmpIndexArray;
  Complex* TmpCoefficientArray; 
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
	    if ((TmpCoefficient.Re != 0.0) || (TmpCoefficient.Im != 0.0))
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
			    TmpSum += (vSource[TmpLeftIndices[j]] * TmpInteractionElements[j] * Conj(TmpCoefficient));
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
		    ChargeContribution.Re *= ChargeContribution.Re;
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

ComplexVector* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::LowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
																		     int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  Complex* Coefficient2 = new Complex [nbrVectors];
  
  Complex Coefficient;
  Complex TmpCoefficient;
  double ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  int Pos2;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  int* TmpIndexArray;
  Complex* TmpCoefficientArray; 
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
		    if ((TmpCoefficient.Re != 0.0) || (TmpCoefficient.Im != 0.0))
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

ComplexVector* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, 
																			      ComplexVector* vDestinations, int nbrVectors, 
																			      int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  Complex* Coefficient2 = new Complex [nbrVectors];
  
  Complex Coefficient;
  Complex TmpCoefficient;
  double ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  int Pos2;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  int* TmpIndexArray;
  Complex* TmpCoefficientArray; 
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

  int Index;  
  double TmpInteraction;
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
		  if ((TmpCoefficient.Re != 0.0) || (TmpCoefficient.Im != 0.0))
		    {
		      NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int l = 0; l < nbrVectors; ++l)
			{
			  vDestinations[l][TmpLeftIndices[j]] += (vSources[l][i] * TmpInteractionElements[j] * TmpCoefficient);
			  TmpSum[l] += (vSources[l][TmpLeftIndices[j]] * Conj(TmpInteractionElements[j] * TmpCoefficient));
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

