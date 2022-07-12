////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//        class of generic hamiltonian for interacting spinless particles     //
//             on lattice written in real space and p-wave pairing            //
//                                                                            //
//                        last modification : 11/06/2016                      //
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


#ifndef PARTICLEONLATTICEREALSPACEPAIRINGHAMILTONIAN_H
#define PARTICLEONLATTICEREALSPACEPAIRINGHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceHamiltonian.h"
#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnLatticeRealSpacePairingHamiltonian : public ParticleOnLatticeRealSpaceHamiltonian
{

 protected:

  // array that contains the number of sites connected to each site via the pairing term
  int* OneBodyGenericPairingNbrConnectedSites;
  // array that contains the indices of the site connected to each site via the pairing term
  int** OneBodyGenericPairingConnectedSites;
  // array that contains all one-body pairing interaction factors
  Complex** OneBodyGenericPairingInteractionFactors;

 public:

  // default constructor
  //
  ParticleOnLatticeRealSpacePairingHamiltonian();

  // constructor
  //
  // nbrParticles = number of particles
  // nbrSites = number of sites
  // tightBinding = hamiltonian corresponding to the tight-binding model in real space
  // densityDensity = matrix that gives the amplitude of each density-density interaction term
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeRealSpacePairingHamiltonian(ParticleOnSphere* particles, int nbrSites, 
					       HermitianMatrix& tightBinding, RealSymmetricMatrix& densityDensity,
					       AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  virtual ~ParticleOnLatticeRealSpacePairingHamiltonian();

  
 protected:
 
  // core part of the FastMultiplication method involving the one-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
							    int* indexArray, Complex* coefficientArray, long& position);

  // core part of the AddMultiply method involving the one-body interaction, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent,
						     int step, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the one-body interaction for a set of vectors, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent,
						     int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors);

  // core part of the PartialFastMultiplicationMemory method involving two-body term and one-body terms
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory);

  // evaluate the one body interaction factors from a tight-binding matrix
  //
  // tightBinding = hamiltonian corresponding to the tight-binding model in real space
  virtual void EvaluateOneBodyFactorsFromTightBingding (HermitianMatrix& tightBinding);
  
};

// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeRealSpacePairingHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent,
															int step, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Index;
  double Coefficient;
  
  if (this->OneBodyGenericInteractionFactors != 0)
    {
      for (int i = firstComponent; i < lastComponent; i += step)
	{ 
	  Coefficient = 0.0;
	  for (int j = 0; j < this->NbrSites; ++j)
	    {
	      int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
	      int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
	      Complex* TmpInteractionFactors = this->OneBodyGenericInteractionFactors[j];
	      Coefficient = particles->A(i,j);
	      double TmpCoefficient;
	      if (Coefficient != 0.0 )
		{
		  for (int k = 0; k < TmpNbrConnectedSites; ++k)
		    {
		      TmpCoefficient = Coefficient;
		      Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient);
		      if (Index < particles->GetHilbertSpaceDimension())
			vDestination[Index] += TmpCoefficient * TmpInteractionFactors[k] * vSource[i];
		    }
		}
	    }
	}
    }

  if (this->OneBodyGenericPairingInteractionFactors != 0)
    {
      for (int i = firstComponent; i < lastComponent; i += step)
	{ 
	  Coefficient = 0.0;
	  for (int j = 0; j < this->NbrSites; ++j)
	    {
	      int TmpNbrConnectedSites = this->OneBodyGenericPairingNbrConnectedSites[j];
	      int* TmpConnectedSites = this->OneBodyGenericPairingConnectedSites[j];
	      Complex* TmpInteractionFactors = this->OneBodyGenericPairingInteractionFactors[j];
	      double TmpCoefficient;
	      for (int k = 0; k < TmpNbrConnectedSites; ++k)
		{
		  Index = particles->AdAd(i, j, TmpConnectedSites[k], TmpCoefficient);
		  if (Index < particles->GetHilbertSpaceDimension())
		    {
		      vDestination[Index] += TmpCoefficient * TmpInteractionFactors[k] * vSource[i];
		    }
		  Index = particles->AA(i, TmpConnectedSites[k], j, TmpCoefficient);
		  if (Index < particles->GetHilbertSpaceDimension())
		    {
		      vDestination[Index] += TmpCoefficient * Conj(TmpInteractionFactors[k]) * vSource[i];
		    }
		}
	    }
	}
    }

  if (this->DiagonalElements != 0)
    {
      for (int i = firstComponent; i < lastComponent; i += step)
	vDestination[i] += this->DiagonalElements[i] * vSource[i];
    }
  else
    {
      for (int i = firstComponent; i < lastComponent; i += step)
	vDestination[i] += this->HamiltonianShift * vSource[i];
    }
}

// core part of the AddMultiply method involving the one-body interaction for a set of vectors, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together

inline void ParticleOnLatticeRealSpacePairingHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent,
													 int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  int Index;
  
  if (this->OneBodyGenericInteractionFactors != 0) 
    {
      for (int p = 0; p < nbrVectors; ++p)
	{
	  ComplexVector& TmpSourceVector = vSources[p];
	  ComplexVector& TmpDestinationVector = vDestinations[p];   
	  for (int i = firstComponent; i < lastComponent; i += step)
	    { 
	      if (this->DiagonalElements == 0)
		TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	      else
		TmpDestinationVector[i] += this->DiagonalElements[i] * TmpSourceVector[i];
	      for (int j = 0; j < this->NbrSites; ++j)
		{
		  int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
		  int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
		  Complex* TmpInteractionFactors = this->OneBodyGenericInteractionFactors[j];
		  Coefficient = particles->A(i,j);
		  double TmpCoefficient;
		  if (Coefficient != 0.0 )
		  {
		    for (int k = 0; k < TmpNbrConnectedSites; ++k)
		      {
			TmpCoefficient = Coefficient;
			Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient);
			if (Index < particles->GetHilbertSpaceDimension())
			  TmpDestinationVector[Index] += TmpCoefficient * TmpInteractionFactors[k] * TmpSourceVector[i];
		      }
		  }
		}
	    }
	}
    }
  else
    {
      for (int p = 0; p < nbrVectors; ++p)
	{
	  ComplexVector& TmpSourceVector = vSources[p];
	  ComplexVector& TmpDestinationVector = vDestinations[p];   
	  if (this->DiagonalElements == 0)
	    {
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		}
	    }
	  else
	    {
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  TmpDestinationVector[i] += this->DiagonalElements[i] * TmpSourceVector[i];
		}
	    }
	}
    }
  if (this->OneBodyGenericPairingInteractionFactors != 0) 
    {
      for (int p = 0; p < nbrVectors; ++p)
	{
	  ComplexVector& TmpSourceVector = vSources[p];
	  ComplexVector& TmpDestinationVector = vDestinations[p];   
	  for (int i = firstComponent; i < lastComponent; i += step)
	    { 
	      if (this->DiagonalElements == 0)
		TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	      else
		TmpDestinationVector[i] += this->DiagonalElements[i] * TmpSourceVector[i];
	      for (int j = 0; j < this->NbrSites; ++j)
		{
		  int TmpNbrConnectedSites = this->OneBodyGenericPairingNbrConnectedSites[j];
		  int* TmpConnectedSites = this->OneBodyGenericPairingConnectedSites[j];
		  Complex* TmpInteractionFactors = this->OneBodyGenericPairingInteractionFactors[j];
		  double TmpCoefficient;
		  for (int k = 0; k < TmpNbrConnectedSites; ++k)
		    {
		      Index = particles->AdAd(i, j, TmpConnectedSites[k], TmpCoefficient);
		      if (Index < particles->GetHilbertSpaceDimension())
			{
			  TmpDestinationVector[Index] += TmpCoefficient * TmpInteractionFactors[k] * TmpSourceVector[i];
			}
		      Index = particles->AA(i, TmpConnectedSites[k], j, TmpCoefficient);
		      if (Index < particles->GetHilbertSpaceDimension())
			{
			  TmpDestinationVector[Index] += TmpCoefficient * Conj(TmpInteractionFactors[k]) * TmpSourceVector[i];
			}		      
		    }
		}
	    }
	}
    }
}

// core part of the FastMultiplication method involving the one-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void ParticleOnLatticeRealSpacePairingHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
												       int* indexArray, Complex* coefficientArray, long& position)
{
  if (this->OneBodyGenericInteractionFactors == 0)
    return;

  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  int Index;

  if (this->HermitianSymmetryFlag == false)
    {
      for (int j = 0; j < this->NbrSites; ++j)
	{
	  int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
	  int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
	  Complex* TmpInteractionFactors = this->OneBodyGenericInteractionFactors[j];
	  Coefficient = particles->A(index, j);
	  double TmpCoefficient;
	  if (Coefficient != 0.0)
	    {
	      for (int k = 0; k < TmpNbrConnectedSites; ++k)
		{
		  TmpCoefficient = Coefficient;
		  Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = TmpCoefficient * TmpInteractionFactors[k];
		      ++position;
		    }
		}
	    }
	  TmpNbrConnectedSites = this->OneBodyGenericPairingNbrConnectedSites[j];
	  TmpConnectedSites = this->OneBodyGenericPairingConnectedSites[j];
	  TmpInteractionFactors = this->OneBodyGenericPairingInteractionFactors[j];
	  for (int k = 0; k < TmpNbrConnectedSites; ++k)
	    {
	      Index = particles->AdAd(index, j, TmpConnectedSites[k], TmpCoefficient);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = TmpCoefficient * TmpInteractionFactors[k];
		  ++position;
		}
	      Index = particles->AA(index, TmpConnectedSites[k], j, TmpCoefficient);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = TmpCoefficient * Conj(TmpInteractionFactors[k]);
		  ++position;
		}
	    }
	}
    }
  else
    {
      for (int j = 0; j < this->NbrSites; ++j)
	{
	  int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
	  int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
	  Complex* TmpInteractionFactors = this->OneBodyGenericInteractionFactors[j];
	  Coefficient = particles->A(index, j);
	  double TmpCoefficient;
	  if (Coefficient != 0.0)
	    {
	      for (int k = 0; k < TmpNbrConnectedSites; ++k)
		{
		  TmpCoefficient = Coefficient;
		  Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient);
		  if (Index <= index)
		    {
		      indexArray[position] = Index;
		      if (Index == index)
			{
			  coefficientArray[position] = 0.5 * TmpCoefficient * TmpInteractionFactors[k];
			}
		      else
			{
			  coefficientArray[position] = TmpCoefficient * TmpInteractionFactors[k];
			}
		      ++position;
		    }
		}
	    }
	  TmpNbrConnectedSites = this->OneBodyGenericPairingNbrConnectedSites[j];
	  TmpConnectedSites = this->OneBodyGenericPairingConnectedSites[j];
	  TmpInteractionFactors = this->OneBodyGenericPairingInteractionFactors[j];
	  for (int k = 0; k < TmpNbrConnectedSites; ++k)
	    {
	      Index = particles->AdAd(index, j, TmpConnectedSites[k], TmpCoefficient);
	      if (Index <= index)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = TmpCoefficient * TmpInteractionFactors[k];
		  ++position;
		}
	      Index = particles->AA(index, TmpConnectedSites[k], j, TmpCoefficient);
	      if (Index <= index)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = TmpCoefficient * Conj(TmpInteractionFactors[k]);
		  ++position;
		}
	    }
	}
    }
}


// core part of the PartialFastMultiplicationMemory method involving two-body term and one-body terms
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void ParticleOnLatticeRealSpacePairingHamiltonian::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory)
{
  if (this->OneBodyGenericInteractionFactors == 0) 
    return;
  int Index;
  double Coefficient = 0.0;
  int Dim = particles->GetHilbertSpaceDimension();
  
  if (this->HermitianSymmetryFlag == false)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrSites; ++j)
	    {
	      int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
	      int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
	      Coefficient = particles->A(i, j);
	      double TmpCoefficient;
	      if (Coefficient != 0.0)
		{
		  for (int k = 0; k < TmpNbrConnectedSites; ++k)
		    {
		      TmpCoefficient = Coefficient;
		      Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		}
	      TmpNbrConnectedSites = this->OneBodyGenericPairingNbrConnectedSites[j];
	      TmpConnectedSites = this->OneBodyGenericPairingConnectedSites[j];
	      for (int k = 0; k < TmpNbrConnectedSites; ++k)
		{
		  Index = particles->AdAd(i, j, TmpConnectedSites[k], TmpCoefficient);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		  Index = particles->AA(i, TmpConnectedSites[k], j, TmpCoefficient);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		}
	    }
	}
    }
  else
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrSites; ++j)
	    {
	      int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
	      int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
	      Coefficient = particles->A (i, j);
	      double TmpCoefficient;
	      if (Coefficient != 0.0)
		{
		  for (int k = 0; k < TmpNbrConnectedSites; ++k)
		    {
		      TmpCoefficient = Coefficient;
		      Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient);
		      if (Index <= i)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		}
	      TmpNbrConnectedSites = this->OneBodyGenericPairingNbrConnectedSites[j];
	      TmpConnectedSites = this->OneBodyGenericPairingConnectedSites[j];
	      for (int k = 0; k < TmpNbrConnectedSites; ++k)
		{
		  Index = particles->AdAd(i, j, TmpConnectedSites[k], TmpCoefficient);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		  Index = particles->AA(i, TmpConnectedSites[k], j, TmpCoefficient);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		}
	    }
	}
    }
}


#endif
