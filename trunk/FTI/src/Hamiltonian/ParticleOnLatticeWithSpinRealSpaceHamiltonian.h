////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                     class author: Cecile Repellin                          //
//                                                                            //
//        class of generic hamiltonian for interacting spinuful particles     //
//                       on lattice written in real space                     //
//                                                                            //
//                        last modification : 15/10/2014                      //
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


#ifndef PARTICLEONLATTICEWITHSPINREALSPACEHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINREALSPACEHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinChernInsulatorHamiltonian.h"
#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnLatticeWithSpinRealSpaceHamiltonian : public ParticleOnLatticeWithSpinChernInsulatorHamiltonian
{

 protected:

  // number of sites 
  int NbrSites;
 
  
  // optional array to store the hamiltonian diagonal elements
  double* DiagonalElements;

  // array that contains the number of up sites connected to each up site
  int* OneBodyGenericNbrConnectedSitesupup;
  // array that contains the number of down sites connected to each down site
  int* OneBodyGenericNbrConnectedSitesdowndown;
  // array that contains the indices of the up site connected to each up site
  int** OneBodyGenericConnectedSitesupup;
  // array that contains the indices of the down site connected to each down site
  int** OneBodyGenericConnectedSitesdowndown;
  // array that contains all one-body interaction factors between up spins
  Complex** OneBodyGenericInteractionFactorsupup;
  // array that contains all one-body interaction factors between down spins
  Complex** OneBodyGenericInteractionFactorsdowndown;

 public:

  // default constructor
  //
  ParticleOnLatticeWithSpinRealSpaceHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSites = number of sites
  // tightBindingupup = hamiltonian corresponding to the tight-binding model in real space for particles with spin up
  // tightBindingdowndown = hamiltonian corresponding to the tight-binding model in real space for particles with spin down
  // densityDensityupup = matrix that gives the amplitude of each density-density interaction term between particles with spin up
  // densityDensitydowndown = matrix that gives the amplitude of each density-density interaction term between particles with spin down
  // densityDensityupdown = matrix that gives the amplitude of each density-density interaction term between particles with spin up and down
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeWithSpinRealSpaceHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSites, 
					HermitianMatrix& tightBindingupup, HermitianMatrix& tightBindingdowndown, RealSymmetricMatrix& densityDensityupup, RealSymmetricMatrix& densityDensitydowndown, RealSymmetricMatrix& densityDensityupdown,
					AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  virtual ~ParticleOnLatticeWithSpinRealSpaceHamiltonian();

  // add an additional S^2 term to the Hamiltonian
  //
  // factor = factor in front of the S^2
  // fixedSz = flag indicating whether Sz needs to be evaluated
  // memory = amount of memory that can be used for S^2  precalculations
  virtual void AddS2 (double factor = 1.0, bool fixedSz = true, long memory = 0l);

  
 protected:
 
  // core part of the FastMultiplication method involving the one-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
							    int* indexArray, Complex* coefficientArray, long& position);

  // core part of the AddMultiply method involving the one-body interaction, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
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
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
						     int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors);

  // core part of the PartialFastMultiplicationMemory method involving two-body term and one-body terms
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory);

  // evaluate the one body interaction factors from a tight-binding matrix
  //
  // tightBindingupup = hamiltonian corresponding to the tight-binding model in real space for spin up particles
  // tightBindingdowndown = hamiltonian corresponding to the tight-binding model in real space for spin down particles
  virtual void EvaluateOneBodyFactorsFromTightBingding (HermitianMatrix& tightBindingupup, HermitianMatrix& tightBindingdowndown);

  // evaluate the two body interaction factors from a generic density-density interaction
  //
  // densityDensityupup = matrix that gives the amplitude of each density-density interaction term for particles with spin up
  // densityDensitydowndown = matrix that gives the amplitude of each density-density interaction term between particles with spin down
  // densityDensityupup = matrix that gives the amplitude of each density-density interaction term for particles with opposite spins
  virtual void EvaluateInteractionFactorsFromDensityDensity (RealSymmetricMatrix& densityDensityupup, RealSymmetricMatrix& densityDensitydowndown, RealSymmetricMatrix& densityDensityupdown) ;
  
};

// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeWithSpinRealSpaceHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
															int step, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Index;
  double Coefficient;
  if (this->OneBodyGenericInteractionFactorsupup != 0)
    {
      if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; i += step)
	    { 
	      Coefficient = 0.0;
	      for (int j = 0; j < this->NbrSites; ++j)
		{
		  Coefficient = particles->Au(i,j);
		  double TmpCoefficient;
		  if (Coefficient != 0.0)
		    {
		      int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
		      int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
		      Complex* TmpInteractionFactorsupup = this->OneBodyGenericInteractionFactorsupup[j];
		      for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient);
			  if (Index < particles->GetHilbertSpaceDimension())
			    vDestination[Index] += TmpCoefficient * TmpInteractionFactorsupup[k] * vSource[i];		  
			}
		    }
		  
		  Coefficient = particles->Ad (i,j);
		  if (Coefficient != 0.0)
		    {
		      int TmpNbrConnectedSitesdowndown = this->OneBodyGenericNbrConnectedSitesdowndown[j];
		      int* TmpConnectedSitesdowndown = this->OneBodyGenericConnectedSitesdowndown[j];
		      Complex* TmpInteractionFactorsdowndown = this->OneBodyGenericInteractionFactorsdowndown[j];
		      for (int k = 0; k < TmpNbrConnectedSitesdowndown; ++k)
			{  
			  TmpCoefficient = Coefficient;
		    Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient);
		    if (Index < particles->GetHilbertSpaceDimension())
		      vDestination[Index] += TmpCoefficient * TmpInteractionFactorsdowndown[k] * vSource[i];
			}
		    }
		}
	    }
	}
      else
	{
	  for (int i = firstComponent; i < lastComponent; i += step)
	    { 
	      Coefficient = 0.0;
	      for (int j = 0; j < this->NbrSites; ++j)
		{
		  Coefficient = particles->Au(i,j);
		  double TmpCoefficient;
		  if (Coefficient != 0.0)
		    {
		      int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
		      int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
		      Complex* TmpInteractionFactorsupup = this->OneBodyGenericInteractionFactorsupup[j];	      
		      for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient);
			  if (Index < particles->GetHilbertSpaceDimension())
			    vDestination[Index] += TmpCoefficient * TmpInteractionFactorsupup[k] * vSource[i];		  
			}
		    }
		}
	    }
	}
    }
  else
    {
      if (this->OneBodyGenericInteractionFactorsdowndown != 0) 
	{
	  for (int i = firstComponent; i < lastComponent; i += step)
	    { 
	      Coefficient = 0.0;
	      for (int j = 0; j < this->NbrSites; ++j)
		{
		  Coefficient = particles->Ad(i,j);
		  double TmpCoefficient;
		  if (Coefficient != 0.0)
		    {
		      int TmpNbrConnectedSitesdowndown = this->OneBodyGenericNbrConnectedSitesdowndown[j];
		      int* TmpConnectedSitesdowndown = this->OneBodyGenericConnectedSitesdowndown[j];
		      Complex* TmpInteractionFactorsdowndown = this->OneBodyGenericInteractionFactorsdowndown[j];	      
		      for (int k = 0; k < TmpNbrConnectedSitesdowndown; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient);
			  if (Index < particles->GetHilbertSpaceDimension())
			    vDestination[Index] += TmpCoefficient * TmpInteractionFactorsdowndown[k] * vSource[i];		  
			}
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

inline void ParticleOnLatticeWithSpinRealSpaceHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
													 int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  int Index;
  if (this->OneBodyGenericInteractionFactorsupup != 0) 
    {
      if (this->OneBodyGenericInteractionFactorsdowndown != 0)
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
		      Coefficient = particles->Au(i,j);
		      if (Coefficient != 0.0)
			{
			  double TmpCoefficient;
			  int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
			  int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
			  Complex* TmpInteractionFactorsupup = this->OneBodyGenericInteractionFactorsupup[j];
			  
			  for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			    {
			      TmpCoefficient = Coefficient;
			      Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient);
			      if (Index < Dim)
				TmpDestinationVector[Index] += Coefficient * TmpInteractionFactorsupup[k] * TmpSourceVector[i];
			    }
			}
		      Coefficient = particles->Ad(i,j);
		      if (Coefficient != 0.0)
			{
			  double TmpCoefficient;
			  int TmpNbrConnectedSitesdowndown = this->OneBodyGenericNbrConnectedSitesdowndown[j];
			  int* TmpConnectedSitesdowndown = this->OneBodyGenericConnectedSitesdowndown[j];
			  Complex* TmpInteractionFactorsdowndown = this->OneBodyGenericInteractionFactorsdowndown[j];
			  for (int k = 0; k < TmpNbrConnectedSitesdowndown; ++k)
			    {
			      Index = particles->Add(TmpConnectedSitesdowndown[k], Coefficient);
			      if (Index < Dim)
				TmpDestinationVector[Index] += Coefficient * TmpInteractionFactorsdowndown[k] * TmpSourceVector[i];
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
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  if (this->DiagonalElements == 0)
		    TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		  else
		    TmpDestinationVector[i] += this->DiagonalElements[i] * TmpSourceVector[i];
		  for (int j = 0; j < this->NbrSites; ++j)
		    {
		      Coefficient = particles->Au(i,j);
		      if (Coefficient != 0.0)
			{
			  double TmpCoefficient;
			  int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
			  int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
			  Complex* TmpInteractionFactorsupup = this->OneBodyGenericInteractionFactorsupup[j];
			  for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			    {
			      TmpCoefficient = Coefficient;
			      Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient);
			      if (Index < Dim)
				TmpDestinationVector[Index] += TmpCoefficient * TmpInteractionFactorsupup[k] * TmpSourceVector[i];
			    }
			}
		    }
		}
	    }      
	}
    }
  else
    {
      if (this->OneBodyGenericInteractionFactorsdowndown != 0)
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
		      Coefficient = particles->Ad(i,j);
		      if (Coefficient != 0.0)
			{
			  double TmpCoefficient;
			  int TmpNbrConnectedSitesdowndown = this->OneBodyGenericNbrConnectedSitesdowndown[j];
			  int* TmpConnectedSitesdowndown = this->OneBodyGenericConnectedSitesdowndown[j];
			  Complex* TmpInteractionFactorsdowndown = this->OneBodyGenericInteractionFactorsdowndown[j];
			  for (int k = 0; k < TmpNbrConnectedSitesdowndown; ++k)
			    {
			      TmpCoefficient = Coefficient;
			      Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient);
			      if (Index < Dim)
				TmpDestinationVector[Index] += TmpCoefficient * TmpInteractionFactorsdowndown[k] * TmpSourceVector[i];
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
    }
}

// core part of the FastMultiplication method involving the one-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void ParticleOnLatticeWithSpinRealSpaceHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
													int* indexArray, Complex* coefficientArray, long& position)
{
  if ((this->OneBodyGenericInteractionFactorsupup == 0) && (this->OneBodyGenericInteractionFactorsdowndown == 0))
    return;

  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  int Index;

  if (this->HermitianSymmetryFlag == false)
    {
      if (this->OneBodyGenericInteractionFactorsupup != 0)
	{
	  if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	    {
	      for (int j = 0; j < this->NbrSites; ++j)
		{
		  Coefficient = particles->Au(index, j);
		  if (Coefficient != 0.0)
		    {
		      double TmpCoefficient;
		      int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
		      int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
		      Complex* TmpInteractionFactorsupup = this->OneBodyGenericInteractionFactorsupup[j];
		      for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient);
			  if (Index < Dim)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = TmpCoefficient * TmpInteractionFactorsupup[k];
			      ++position;
			    }
			}
		    }
		  Coefficient = particles->Ad(index, j);
		  if (Coefficient != 0.0)
		    {
		      double TmpCoefficient;
		      int TmpNbrConnectedSitesdowndown = this->OneBodyGenericNbrConnectedSitesdowndown[j];
		      int* TmpConnectedSitesdowndown = this->OneBodyGenericConnectedSitesdowndown[j];
		      Complex* TmpInteractionFactorsdowndown = this->OneBodyGenericInteractionFactorsdowndown[j];
		      for (int k = 0; k < TmpNbrConnectedSitesdowndown; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient);
			  if (Index < Dim)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = TmpCoefficient * TmpInteractionFactorsdowndown[k];
			      ++position;
			    }
			}
		    }
		}
	    }
	  else
	    {
	      for (int j = 0; j < this->NbrSites; ++j)
		{
		  Coefficient = particles->Au (index, j);
		  if (Coefficient != 0.0)
		    {
		      double TmpCoefficient;
		      int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
		      int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
		      Complex* TmpInteractionFactorsupup = this->OneBodyGenericInteractionFactorsupup[j];
		      for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient);
			  if (Index < Dim)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = TmpCoefficient * TmpInteractionFactorsupup[k];
			      ++position;
			    }
			}
		    }
		}	  
	    }
	}
      else
	{
	  if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	    {
	      for (int j = 0; j < this->NbrSites; ++j)
		{
		  Coefficient = particles->Ad (index, j);
		  if (Coefficient != 0.0)
		    {
		      double TmpCoefficient;
		      int TmpNbrConnectedSitesdowndown = this->OneBodyGenericNbrConnectedSitesdowndown[j];
		      int* TmpConnectedSitesdowndown = this->OneBodyGenericConnectedSitesdowndown[j];
		      Complex* TmpInteractionFactorsdowndown = this->OneBodyGenericInteractionFactorsdowndown[j];
		      for (int k = 0; k < TmpNbrConnectedSitesdowndown; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient);
			  if (Index < Dim)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = TmpCoefficient * TmpInteractionFactorsdowndown[k];
			      ++position;
			    }
			}
		    }
		}
	    } 
	}
    }
  else
    {
      if (this->OneBodyGenericInteractionFactorsupup != 0)
	{
	  if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	    {
	      for (int j = 0; j < this->NbrSites; ++j)
		{
		  Coefficient = particles->Au(index, j);
		  if (Coefficient != 0.0)
		    {
		      double TmpCoefficient;
		      int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
		      int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
		      Complex* TmpInteractionFactorsupup = this->OneBodyGenericInteractionFactorsupup[j];
		      for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->Adu( TmpConnectedSitesupup[k], TmpCoefficient);
			  if (Index <= index)
			    {
			      indexArray[position] = Index;
			      if (Index == index)
				{
				  coefficientArray[position] = 0.5 * TmpCoefficient * TmpInteractionFactorsupup[k];
				}
			      else
				{
				  coefficientArray[position] = TmpCoefficient * TmpInteractionFactorsupup[k];
				}
			      ++position;
			    }
			}
		    }
		  
		  Coefficient = particles->Ad(index, j);
		  if (Coefficient != 0.0)
		    {
		      double TmpCoefficient;
		      int TmpNbrConnectedSitesdowndown = this->OneBodyGenericNbrConnectedSitesdowndown[j];
		      int* TmpConnectedSitesdowndown = this->OneBodyGenericConnectedSitesdowndown[j];
		      Complex* TmpInteractionFactorsdowndown = this->OneBodyGenericInteractionFactorsdowndown[j];
		      for (int k = 0; k < TmpNbrConnectedSitesdowndown; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient);
			  if (Index <= index)
			    {
			      indexArray[position] = Index;
			      if (Index == index)
				{
				  coefficientArray[position] = 0.5 * TmpCoefficient * TmpInteractionFactorsdowndown[k];
				}
			      else
				{
				  coefficientArray[position] = TmpCoefficient * TmpInteractionFactorsdowndown[k];
				}
			      ++position;
			    }
			}
		    }
		}
	    }
	  else
	    {
	      for (int j = 0; j < this->NbrSites; ++j)
		{
		  Coefficient = particles->Au(index, j);
		  if (Coefficient != 0.0)
		    {
		      double TmpCoefficient;
		      int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
		      int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
		      Complex* TmpInteractionFactorsupup = this->OneBodyGenericInteractionFactorsupup[j];
		      for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient);
			  if (Index <= index)
			    {
			      indexArray[position] = Index;
			      if (Index == index)
				{
				  coefficientArray[position] = 0.5 * TmpCoefficient * TmpInteractionFactorsupup[k];
				}
			      else
				{
				  coefficientArray[position] = TmpCoefficient * TmpInteractionFactorsupup[k];
				}
			      ++position;
			    }
			}
		    }
		}
	    }
	}
      else
	{
	  for (int j = 0; j < this->NbrSites; ++j)
	    {
	      Coefficient = particles->Ad(index, j);
	      if (Coefficient != 0.0)
		{
		  double TmpCoefficient;
		  int TmpNbrConnectedSitesdowndown = this->OneBodyGenericNbrConnectedSitesdowndown[j];
		  int* TmpConnectedSitesdowndown = this->OneBodyGenericConnectedSitesdowndown[j];
		  Complex* TmpInteractionFactorsdowndown = this->OneBodyGenericInteractionFactorsdowndown[j];
		  for (int k = 0; k < TmpNbrConnectedSitesdowndown; ++k)
		    {
		      TmpCoefficient = Coefficient;
		      Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient);
		      if (Index <= index)
			{
			  indexArray[position] = Index;
			  if (Index == index)
			    {
			      coefficientArray[position] = 0.5 * TmpCoefficient * TmpInteractionFactorsdowndown[k];
			    }
			      else
				{
				  coefficientArray[position] = TmpCoefficient * TmpInteractionFactorsdowndown[k];
				}
			  ++position;
			}
		    }
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

inline void ParticleOnLatticeWithSpinRealSpaceHamiltonian::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  if ((this->OneBodyGenericInteractionFactorsupup == 0) && (this->OneBodyGenericInteractionFactorsdowndown == 0))
    return;
  int Index;
  double Coefficient = 0.0;
  int Dim = particles->GetHilbertSpaceDimension();
  
  if (this->HermitianSymmetryFlag == false)
    {
      if (this->OneBodyGenericInteractionFactorsupup != 0)
	{
	  if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	    {
	      for (int i = firstComponent; i < lastComponent; ++i)
		{
		  for (int j = 0; j < this->NbrSites; ++j)
		    {
		      Coefficient = particles->Au(i,j);
		      if (Coefficient != 0.0)
			{
			  double TmpCoefficient;
			  int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
			  int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
			  for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			    {
			      TmpCoefficient = Coefficient;
			      Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient);
			      if (Index < Dim)
				{
				  ++memory;
				  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
				}
			    }
			}
		      Coefficient = particles->Ad(i,j);
		      if (Coefficient != 0.0)
			{
			  double TmpCoefficient;
			  int TmpNbrConnectedSitesdowndown = this->OneBodyGenericNbrConnectedSitesdowndown[j];
			  int* TmpConnectedSitesdowndown = this->OneBodyGenericConnectedSitesdowndown[j];
			  for (int k = 0; k < TmpNbrConnectedSitesdowndown; ++k)
			    {
			      TmpCoefficient = Coefficient;
			      Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient);
			      if (Index < Dim)
				{
				  ++memory;
				  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
				}
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
		      Coefficient = particles->Au(i,j);
		      if (Coefficient != 0.0)
			{
			  double TmpCoefficient;
			  int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
			  int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
			  for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			    {
			      TmpCoefficient = Coefficient;
			      Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient);
			      if (Index < Dim)
				{
				  ++memory;
				  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
				}
			    }
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
		  Coefficient = particles->Ad(i,j);
		  if (Coefficient != 0.0)
		    {
		      double TmpCoefficient;
		      int TmpNbrConnectedSitesdowndown = this->OneBodyGenericNbrConnectedSitesdowndown[j];
		      int* TmpConnectedSitesdowndown = this->OneBodyGenericConnectedSitesdowndown[j];
		      for (int k = 0; k < TmpNbrConnectedSitesdowndown; ++k)
			{
			  TmpCoefficient = Coefficient;
			      Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient);
			      if (Index < Dim)
				{
				  ++memory;
				  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
				}
			}
		    }
		}
	    }
	}
    }
  else
    {
      if (this->OneBodyGenericInteractionFactorsupup != 0)
	{
	  if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	    {
	      for (int i = firstComponent; i < lastComponent; ++i)
		{
		  for (int j = 0; j < this->NbrSites; ++j)
		    {
		      Coefficient = particles->Au(i,j);
		      if (Coefficient != 0.0)
			{
			  double TmpCoefficient;
			  int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
			  int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
			  for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			    {
			      TmpCoefficient = Coefficient;
			      Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient);
			      if (Index <= i)
				{
				  ++memory;
				  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
				}
			    }
			}
		      
		      Coefficient = particles->Ad(i,j);
		      if (Coefficient != 0.0)
			{
			  double TmpCoefficient;
			  int TmpNbrConnectedSitesdowndown = this->OneBodyGenericNbrConnectedSitesdowndown[j];
			  int* TmpConnectedSitesdowndown = this->OneBodyGenericConnectedSitesdowndown[j];
			  for (int k = 0; k < TmpNbrConnectedSitesdowndown; ++k)
			    {
			      TmpCoefficient = Coefficient;
			      Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient);
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
	  else
	    {
	      for (int i = firstComponent; i < lastComponent; ++i)
		{
		  for (int j = 0; j < this->NbrSites; ++j)
		    {
		      Coefficient = particles->Au(i,j);
		      if (Coefficient != 0.0)
			{
			  double TmpCoefficient;
			  int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
			  int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
			  for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			    {
			      TmpCoefficient = Coefficient;
			      Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient);
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
	}
      else
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j = 0; j < this->NbrSites; ++j)
		{
		  Coefficient = particles->Ad(i,j);
		  if (Coefficient != 0.0)
		    {
		      double TmpCoefficient;
		      int TmpNbrConnectedSitesdowndown = this->OneBodyGenericNbrConnectedSitesdowndown[j];
		      int* TmpConnectedSitesdowndown = this->OneBodyGenericConnectedSitesdowndown[j];
		      for (int k = 0; k < TmpNbrConnectedSitesdowndown; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient);
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
    }
}


#endif
