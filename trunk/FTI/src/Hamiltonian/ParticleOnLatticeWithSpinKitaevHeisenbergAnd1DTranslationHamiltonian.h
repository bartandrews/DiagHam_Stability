////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//                   class of Hubbard hamiltonian associated                  //
//  to particles on a lattice with spin and translation in the x direction    //
//                                                                            //
//                        last modification : 13/07/2014                      //
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


#ifndef PARTICLEONLATTICEWITHSPINKITAEVHEISENBERGAND1DTRANSLATIONHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINKITAEVHEISENBERGAND1DTRANSLATIONHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian : public ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian
{

 protected:
  
  // momentum along the x direction
  int XMomentum;
  // periodicity in the x direction
  int MaxMomentum;
  // periodicity with respect to site numbering
  int XIncrement;
  
  //array containing all the phase factors that are needed when computing matrix elements
  Complex* ExponentialFactors;

 public:

  // default constructor
  //
  ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSites = number of sites
  // momentum = momentum sector
  // periodicity = periodicity with respect to site numbering 
  // kineticFactorIntra = multiplicative factor in front of the intraspin kinetic term
  // uPotential = Hubbard potential strength
  // j1Factor = strength of the isotropic nearest neighbor spin interaction
  // j2Factor = strength of the anisotropic nearest neighbor spin interaction
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSite, int momentum, int periodicity, char* geometryFile, double kineticFactorIsotropic, double kineticFactorAnisotropic, double uPotential, double j1Factor, double j2Factor, AbstractArchitecture* architecture, long memory = -1);

  // constructor from the explicit the bond description
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSites = number of sites
  // nbrBonds = number of bonds
  // sitesA = array of A sites for each bond
  // sitesB = array of B sites for each bond
  // bondTypes = array that describe each type of bond (0 for x, 1 fo y, 2 for z)
  // momentum = momentum sector
  // periodicity = periodicity with respect to site numbering 
  // kineticFactorIntra = multiplicative factor in front of the intraspin kinetic term
  // uPotential = Hubbard potential strength
  // j1Factor = strength of the isotropic nearest neighbor spin interaction
  // j2Factor = strength of the anisotropic nearest neighbor spin interaction
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSite, 
								       int nbrBonds,  int* sitesA, int* sitesB, int* bondTypes,
								       int momentum, int periodicity, 
								       double kineticFactorIsotropic, double kineticFactorAnisotropic, double uPotential, double j1Factor, double j2Factor, 
								       AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  virtual ~ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian();

  
 protected:
 
  // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the two-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
						     ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the two-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  inline void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
							     ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // core part of the FastMultiplication method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray  
  virtual void EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
							    int* indexArray, Complex* coefficientArray, long& position);

  // core part of the PartialFastMultiplicationMemory method involving two-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations  
  virtual void EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory);

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

  // evaluate all exponential factors
  //   
  virtual void EvaluateExponentialFactors();

};

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  int* TmpIndices2;
  Complex* TmpInteractionFactor;
  int Index;
  int NbrTranslations;
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
      TmpIndices = this->IntraSectorIndicesPerSum[j];
      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
      TmpIndices2 = this->InterSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupupupup[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupup[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndowndowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	}
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  Coefficient3 = particles->AuAd(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupupupdown[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupdown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	}
    }
}

// core part of the AddMultiply method involving the two-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
												      ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  int* TmpIndices;
  int* TmpIndices2;
  Complex* TmpInteractionFactor;
  int NbrTranslations;
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
      TmpIndices = this->IntraSectorIndicesPerSum[j];
      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
      TmpIndices2 = this->InterSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupupupup[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupup[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndowndowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	    }
	}
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  Coefficient3 = particles->AuAd(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupupupdown[j][(i1 * Lim2) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupdown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	    }
	}
    }
}

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added  

inline void ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  //int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  int* TmpIndices2;
  Complex* TmpInteractionFactor;
  int Index;
  Complex TmpSum = 0.0;
  int NbrTranslations;
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
      TmpIndices = this->IntraSectorIndicesPerSum[j];
      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
      TmpIndices2 = this->InterSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupupupup[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupup[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndowndowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	}
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  Coefficient3 = particles->AuAd(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupupupdown[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupdown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	}
    }
  vDestination[index] += TmpSum;
}

// core part of the AddMultiply method involving the two-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
																 ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  //int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  int* TmpIndices;
  int* TmpIndices2;
  Complex* TmpInteractionFactor;
  Complex* TmpSum = new Complex[nbrVectors];
  int NbrTranslations;
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
      TmpIndices = this->IntraSectorIndicesPerSum[j];
      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
      TmpIndices2 = this->InterSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupupupup[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupup[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndowndowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	    }
	}
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  Coefficient3 = particles->AuAd(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupupupdown[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupdown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations] * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	    }
	}
    }
  for (int l = 0; l < nbrVectors; ++l)
    vDestinations[l][index] += TmpSum[l];
  delete[] TmpSum;
}

// core part of the FastMultiplication method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
															       int* indexArray, Complex* coefficientArray, long& position)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  int* TmpIndices2;
  Complex* TmpInteractionFactor;
  int Index;
  int NbrTranslations;

  if (this->HermitianSymmetryFlag == false)
    {
      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
	  TmpIndices = this->IntraSectorIndicesPerSum[j];
	  int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
	  TmpIndices2 = this->InterSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupupupup[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsdowndownupup[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsdowndowndowndown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	    }
	  for (int i1 = 0; i1 < Lim2; i1 += 2)
	    {
	      Coefficient3 = particles->AuAd(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupupupdown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsdowndownupdown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	    }
	}
    }
  else
    {
      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
	  TmpIndices = this->IntraSectorIndicesPerSum[j];
	  int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
	  TmpIndices2 = this->InterSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupupupup[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  if (Index == index)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsdowndownupup[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  if (Index == index)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  if (Index == index)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  if (Index == index)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsdowndowndowndown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  if (Index == index)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  if (Index == index)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		}
	    }
	  for (int i1 = 0; i1 < Lim2; i1 += 2)
	    {
	      Coefficient3 = particles->AuAd(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupupupdown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  if (Index == index)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsdowndownupdown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  if (Index == index)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  if (Index == index)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		}
	    }
	}
    }
}

// core part of the PartialFastMultiplicationMemory method involving two-body term
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  int* TmpIndices2;
  Complex* TmpInteractionFactor;
  int Index;
  int NbrTranslations;

  if (this->HermitianSymmetryFlag == false)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
	      TmpIndices = this->IntraSectorIndicesPerSum[j];
	      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
	      TmpIndices2 = this->InterSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient3 = particles->AuAu(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsupupupup[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactorsdowndownupup[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  Coefficient3 = particles->AdAd(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactorsdowndowndowndown[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		    }
		}
	      for (int i1 = 0; i1 < Lim2; i1 += 2)
		{
		  Coefficient3 = particles->AuAd(i, TmpIndices2[i1], TmpIndices2[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsupupupdown[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactorsdowndownupdown[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
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
	  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
	      TmpIndices = this->IntraSectorIndicesPerSum[j];
	      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
	      TmpIndices2 = this->InterSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient3 = particles->AuAu(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsupupupup[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactorsdowndownupup[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  Coefficient3 = particles->AdAd(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactorsdowndowndowndown[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		    }
		}
	      for (int i1 = 0; i1 < Lim2; i1 += 2)
		{
		  Coefficient3 = particles->AuAd(i, TmpIndices2[i1], TmpIndices2[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsupupupdown[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactorsdowndownupdown[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		    }
		}
	    }
	}
    }
}


// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
															int step, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Index;
  double Coefficient;
  int NbrTranslations;
  if (this->HermitianSymmetryFlag == false)
    {
      if (this->OneBodyGenericInteractionFactorsupup != 0)
	if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	  {
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		for (int j = 0; j < this->NbrSite; ++j) 
		  {
		    for (int k = 0; k < 3; ++k)
		      {
			int j2 = this->MapNearestNeighborBonds[j][k];
			if (j2 < this->NbrSite)
			  {
			    Index = particles->AduAu(i , j, j2, Coefficient, NbrTranslations);
			    if (Index < particles->GetHilbertSpaceDimension())
			      vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k] * this->ExponentialFactors[NbrTranslations] * vSource[i];
			    
			    Index = particles->AddAd(i , j, j2, Coefficient, NbrTranslations);
			    if (Index < particles->GetHilbertSpaceDimension())
			      vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k] * this->ExponentialFactors[NbrTranslations] * vSource[i];		  
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
		for (int j = 0; j < this->NbrSite; ++j) 
		  {
		    for (int k = 0; k < 3; ++k)
		      {
			int j2 = this->MapNearestNeighborBonds[j][k];
			if (j2 < this->NbrSite)
			  {
			    Index = particles->AduAu(i, j, j2, Coefficient, NbrTranslations);
			    if (Index < particles->GetHilbertSpaceDimension())
			      vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k] * this->ExponentialFactors[NbrTranslations] * vSource[i];
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
		  for (int j = 0; j < this->NbrSite; ++j) 
		    {
		      for (int k = 0; k < 3; ++k)
			{
			  int j2 = this->MapNearestNeighborBonds[j][k];
			  if (j2 < this->NbrSite)
			    {
			      Index = particles->AddAd(i, j, j2, Coefficient, NbrTranslations);
			      if (Index < particles->GetHilbertSpaceDimension())
				vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k] * this->ExponentialFactors[NbrTranslations] * vSource[i];		  
			    }
			}
		    }
		}
	    }	
	}
      for (int i = firstComponent; i < lastComponent; i += step)
	vDestination[i] += this->HamiltonianShift * vSource[i];
      
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  double Coefficient;
	  Complex Source;
	  int Dim = particles->GetHilbertSpaceDimension();
	  int Index;
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      Source = vSource[i];
	      for (int j = 0; j < this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAu(i, j, j2, Coefficient, NbrTranslations);
			  if (Index < Dim)
			    vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k])) * Source;
			  Index = particles->AduAd(i, j, j2, Coefficient, NbrTranslations);
			  if (Index < Dim)
			    vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupdown[j][k]) * Source;
			}
		    }
		}
	    }
	}
    }
  else
    {
      if (this->OneBodyGenericInteractionFactorsupup != 0)
	if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	  {
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		for (int j = 0; j < this->NbrSite; ++j) 
		  {
		    for (int k = 0; k < 3; ++k)
		      {
			int j2 = this->MapNearestNeighborBonds[j][k];
			if (j2 < this->NbrSite)
			  {
			    Index = particles->AduAu(i , j, j2, Coefficient, NbrTranslations);
			    if (Index <= i)
			      {
				vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k] * this->ExponentialFactors[NbrTranslations] * vSource[i];
				if (Index < i)
				  {
				    vDestination[i] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsupup[j][k] * this->ExponentialFactors[NbrTranslations]) * vSource[Index];
				  }
			      }
			    
			    Index = particles->AddAd(i , j, j2, Coefficient, NbrTranslations);
			    if (Index <= i)
			      {
				vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k] * this->ExponentialFactors[NbrTranslations] * vSource[i];
				if (Index < i)
				  {
				    vDestination[i] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsdowndown[j][k] * this->ExponentialFactors[NbrTranslations]) * vSource[Index];
				  }
			      }
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
		for (int j = 0; j < this->NbrSite; ++j) 
		  {
		    for (int k = 0; k < 3; ++k)
		      {
			int j2 = this->MapNearestNeighborBonds[j][k];
			if (j2 < this->NbrSite)
			  {
			    Index = particles->AduAu(i, j, j2, Coefficient, NbrTranslations);
			    if (Index <= i)
			      {
				vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k] * this->ExponentialFactors[NbrTranslations] * vSource[i];
				if (Index < i)
				  {
				    vDestination[i] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsupup[j][k] * this->ExponentialFactors[NbrTranslations]) * vSource[Index];
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
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  Coefficient = 0.0;
		  for (int j = 0; j < this->NbrSite; ++j) 
		    {
		      for (int k = 0; k < 3; ++k)
			{
			  int j2 = this->MapNearestNeighborBonds[j][k];
			  if (j2 < this->NbrSite)
			    {
			      Index = particles->AddAd(i, j, j2, Coefficient, NbrTranslations);
			      if (Index <= i)
				{
				  vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k] * this->ExponentialFactors[NbrTranslations] * vSource[i];
				  if (Index < i)
				    {
				      vDestination[i] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsdowndown[j][k] * this->ExponentialFactors[NbrTranslations]) * vSource[Index];
				    }
				}
			    }
			}
		    }
		}
	    }	
	}
      for (int i = firstComponent; i < lastComponent; i += step)
	vDestination[i] += this->HamiltonianShift * vSource[i];
      
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  double Coefficient;
	  Complex Source;
	  int Dim = particles->GetHilbertSpaceDimension();
	  int Index;
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      Source = vSource[i];
	      for (int j = 0; j < this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAu(i, j, j2, Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k])) * Source;
			      if (Index < i)
				{
				  vDestination[i] += (Coefficient * Conj(this->ExponentialFactors[NbrTranslations]) * this->OneBodyGenericInteractionFactorsupdown[j][k]) * vSource[Index];
				}
			    }
			  Index = particles->AduAd(i, j, j2, Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupdown[j][k]) * Source;
			      if (Index < i)
				{
				  vDestination[i] += (Coefficient * Conj(this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupdown[j][k])) * vSource[Index];
				}
			    }
			}
		    }
		}
	    }
	}
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

inline void ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
															int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  int Index;
  int NbrTranslations;
  if (this->HermitianSymmetryFlag == false)
    {
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
		      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		      for (int j = 0; j < this->NbrSite; ++j) 
			{
			  for (int k = 0; k < 3; ++k)
			    {
			      int j2 = this->MapNearestNeighborBonds[j][k];
			      if (j2 < this->NbrSite)
				{
				  Index = particles->AduAu(i, j, j2, Coefficient, NbrTranslations);
				  if (Index < Dim)
				    TmpDestinationVector[Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupup[j][k] * TmpSourceVector[i];
				  
				  Index = particles->AddAd(i, j, j2, Coefficient, NbrTranslations);
				  if (Index < Dim)
				    TmpDestinationVector[Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsdowndown[j][k] * TmpSourceVector[i];
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
		      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		      for (int j = 0; j < this->NbrSite; ++j) 
			{
			  for (int k = 0; k < 3; ++k)
			    {
			      int j2 = this->MapNearestNeighborBonds[j][k];
			      if (j2 < this->NbrSite)
				{
				  Index = particles->AduAu(i, j, j2, Coefficient, NbrTranslations);
				  if (Index < Dim)
				    TmpDestinationVector[Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupup[j][k] * TmpSourceVector[i];
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
		      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		      for (int j = 0; j < this->NbrSite; ++j) 
			{
			  for (int k = 0; k < 3; ++k)
			    {
			      int j2 = this->MapNearestNeighborBonds[j][k];
			      if (j2 < this->NbrSite)
				{
				  Index = particles->AddAd(i, j, j2, Coefficient, NbrTranslations);
				  if (Index < Dim)
				    TmpDestinationVector[Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsdowndown[j][k] * TmpSourceVector[i];
				}
			    }
			}
		    }
		}
	    }	
	}
      for (int p = 0; p < nbrVectors; ++p)
	{
	  ComplexVector& TmpSourceVector = vSources[p];
	  ComplexVector& TmpDestinationVector = vDestinations[p];
	  for (int i = firstComponent; i < lastComponent; i += step)
	    TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	}
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      for (int j = 0; j < this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAu(i, j, j2, Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      for (int p = 0; p < nbrVectors; ++p)
				vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k]) * vSources[p][i];
			    }
			  Index = particles->AduAd(i, j, j2, Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      for (int p = 0; p < nbrVectors; ++p)
				vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupdown[j][k] * vSources[p][i];
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
	      for (int p = 0; p < nbrVectors; ++p)
		{
		  ComplexVector& TmpSourceVector = vSources[p];
		  ComplexVector& TmpDestinationVector = vDestinations[p];
		  
		  for (int i = firstComponent; i < lastComponent; i += step)
		    { 
		      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		      for (int j = 0; j < this->NbrSite; ++j) 
			{
			  for (int k = 0; k < 3; ++k)
			    {
			      int j2 = this->MapNearestNeighborBonds[j][k];
			      if (j2 < this->NbrSite)
				{
				  Index = particles->AduAu(i, j, j2, Coefficient, NbrTranslations);
				  if (Index <= i)
				    {
				      TmpDestinationVector[Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupup[j][k] * TmpSourceVector[i];
				      if (Index < i)
					{
					  TmpDestinationVector[i] += Coefficient * Conj(this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupup[j][k]) * TmpSourceVector[Index];
					}
				    }
				  
				  Index = particles->AddAd(i, j, j2, Coefficient, NbrTranslations);
				  if (Index <= i)
				    {
				      TmpDestinationVector[Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsdowndown[j][k] * TmpSourceVector[i];
				      if (Index < i)
					{
					  TmpDestinationVector[i] += Coefficient * Conj(this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsdowndown[j][k]) * TmpSourceVector[Index];
					}
				    }
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
		      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		      for (int j = 0; j < this->NbrSite; ++j) 
			{
			  for (int k = 0; k < 3; ++k)
			    {
			      int j2 = this->MapNearestNeighborBonds[j][k];
			      if (j2 < this->NbrSite)
				{
				  Index = particles->AduAu(i, j, j2, Coefficient, NbrTranslations);
				  if (Index <= i)
				    {
				      TmpDestinationVector[Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupup[j][k] * TmpSourceVector[i];
				      if (Index < i)
					{
					  TmpDestinationVector[i] += Coefficient * Conj(this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupup[j][k]) * TmpSourceVector[Index];
					}
				    }
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
		      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		      for (int j = 0; j < this->NbrSite; ++j) 
			{
			  for (int k = 0; k < 3; ++k)
			    {
			      int j2 = this->MapNearestNeighborBonds[j][k];
			      if (j2 < this->NbrSite)
				{
				  Index = particles->AddAd(i, j, j2, Coefficient, NbrTranslations);
				  if (Index <= i)
				    {
				      TmpDestinationVector[Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsdowndown[j][k] * TmpSourceVector[i];
				      if (Index < i)
					{
					  TmpDestinationVector[i] += Coefficient * Conj(this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsdowndown[j][k]) * TmpSourceVector[Index];
					}
				    }
				}
			    }
			}
		    }
		}
	    }	
	}
      for (int p = 0; p < nbrVectors; ++p)
	{
	  ComplexVector& TmpSourceVector = vSources[p];
	  ComplexVector& TmpDestinationVector = vDestinations[p];
	  for (int i = firstComponent; i < lastComponent; i += step)
	    TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	}
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      for (int j = 0; j < this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAu(i, j, j2, Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      if (Index < i)
				{
				  for (int p = 0; p < nbrVectors; ++p)
				    {
				      vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k]) * vSources[p][i];
				      vDestinations[p][i] += Coefficient * Conj(this->ExponentialFactors[NbrTranslations]) * this->OneBodyGenericInteractionFactorsupdown[j][k] * vSources[p][Index];
				    }
				}
			      else
				{
				  for (int p = 0; p < nbrVectors; ++p)
				    vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k]) * vSources[p][i];
				}
			    }
			  Index = particles->AduAd(i, j, j2, Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      if (Index < i)
				{
				  for (int p = 0; p < nbrVectors; ++p)
				    {
				      vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupdown[j][k] * vSources[p][i];
				      vDestinations[p][i] += Coefficient * Conj(this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupdown[j][k]) * vSources[p][Index];
				    }
				}
			      else
				{
				  for (int p = 0; p < nbrVectors; ++p)
				    vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupdown[j][k] * vSources[p][i];
				}
			    }
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

inline void ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
															       int* indexArray, Complex* coefficientArray, long& position)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  int Index;
  int NbrTranslations;


  if (this->HermitianSymmetryFlag == false)
    {
      if ((this->OneBodyGenericInteractionFactorsdowndown != 0) && (this->OneBodyGenericInteractionFactorsupup != 0))
	{
	  for (int j = 0; j < this->NbrSite; ++j)
	    {
	      for (int k = 0; k < 3; ++k)
		{
		  int j2 = this->MapNearestNeighborBonds[j][k];
		  if (j2 < this->NbrSite)
		    {
		      Index = particles->AduAu(index + this->PrecalculationShift, j, j2, Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupup[j][k];
			  ++position;
			}
		      Index = particles->AddAd(index + this->PrecalculationShift, j, j2, Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsdowndown[j][k];
			  ++position;
			}
		    }
		}
	    }
	}
      else
	{
	  if (this->OneBodyGenericInteractionFactorsupup != 0)
	    {
	      for (int j = 0; j < this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AduAu(index + this->PrecalculationShift, j, j2, Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupup[j][k];
			      ++position;
			    }
			}
		    }
		}
	    }
	  else
	    {
	      if (this->OneBodyGenericInteractionFactorsdowndown != 0)
		{
		  for (int j = 0; j < this->NbrSite; ++j)
		    {
		      for (int k = 0; k < 3; ++k)
			{
			  int j2 = this->MapNearestNeighborBonds[j][k];
			  if (j2 < this->NbrSite)
			    {
			      Index = particles->AddAd(index + this->PrecalculationShift, j, j2, Coefficient, NbrTranslations);
			      if (Index < Dim)
				{
				  indexArray[position] = Index;
				  coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsdowndown[j][k];
				  ++position;
				}
			    }
			}
		    }
		} 
	    }
	}
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  for (int j = 0; j < this->NbrSite; ++j)
	    {
	      for (int k = 0; k < 3; ++k)
		{
		  int j2 = this->MapNearestNeighborBonds[j][k];
		  if (j2 < this->NbrSite)
		    {
		      Index = particles->AddAu(index + this->PrecalculationShift, j, j2, Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k]);
			  ++position;
			}
		      Index = particles->AduAd(index + this->PrecalculationShift, j, j2, Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupdown[j][k];
			  ++position;
			}
		    }
		}
	    }
	}
    }
  else
    {
      if ((this->OneBodyGenericInteractionFactorsdowndown != 0) && (this->OneBodyGenericInteractionFactorsupup != 0))
	{
	  for (int j = 0; j < this->NbrSite; ++j)
	    {
	      for (int k = 0; k < 3; ++k)
		{
		  int j2 = this->MapNearestNeighborBonds[j][k];
		  if (j2 < this->NbrSite)
		    {
		      Index = particles->AduAu(index, j, j2, Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  indexArray[position] = Index;
			  if (Index == index)
			    {
			      coefficientArray[position] = 0.5 * Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupup[j][k];
			    }
			  else
			    {
			      coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupup[j][k];
			    }
			  ++position;
			}
		      Index = particles->AddAd(index, j, j2, Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  indexArray[position] = Index;
			  if (Index == index)
			    {
			      coefficientArray[position] = 0.5 * Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsdowndown[j][k];
			    }
			  else
			    {
			      coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsdowndown[j][k];
			    }
			  ++position;
			}
		    }
		}
	    }
	}
      else
	{
	  if (this->OneBodyGenericInteractionFactorsupup != 0)
	    {
	      for (int j = 0; j < this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AduAu(index, j, j2, Coefficient, NbrTranslations);
			  if (Index <= index)
			    {
			      indexArray[position] = Index;
			      if (Index == index)
				{
				  coefficientArray[position] = 0.5 * Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupup[j][k];
				}
			      else
				{
				  coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupup[j][k];
				}
			      ++position;
			    }
			}
		    }
		}
	    }
	  else
	    {
	      if (this->OneBodyGenericInteractionFactorsdowndown != 0)
		{
		  for (int j = 0; j < this->NbrSite; ++j)
		    {
		      for (int k = 0; k < 3; ++k)
			{
			  int j2 = this->MapNearestNeighborBonds[j][k];
			  if (j2 < this->NbrSite)
			    {
			      Index = particles->AddAd(index, j, j2, Coefficient, NbrTranslations);
			      if (Index <= index)
				{
				  indexArray[position] = Index;
				  if (Index == index)
				    {
				      coefficientArray[position] = 0.5 * Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsdowndown[j][k];
				    }
				  else
				    {
				      coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsdowndown[j][k];
				    }
				  ++position;
				}
			    }
			}
		    }
		} 
	    }
	}
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  for (int j = 0; j < this->NbrSite; ++j)
	    {
	      for (int k = 0; k < 3; ++k)
		{
		  int j2 = this->MapNearestNeighborBonds[j][k];
		  if (j2 < this->NbrSite)
		    {
		      Index = particles->AddAu(index, j, j2, Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  indexArray[position] = Index;
			  if (Index == index)
			    {
			      coefficientArray[position] = 0.5 * Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k]);
			    }
			  else
			    {
			      coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k]);
			    }
			  ++position;
			}
		      Index = particles->AduAd(index, j, j2, Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  indexArray[position] = Index;
			  if (Index == index)
			    {
			      coefficientArray[position] = 0.5 * Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupdown[j][k];
			    }
			  else
			    {
			      coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyGenericInteractionFactorsupdown[j][k];
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

inline void ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  int NbrTranslations;
  int Dim = particles->GetHilbertSpaceDimension();
  
  if (this->HermitianSymmetryFlag == false)
    {
      if (this->OneBodyGenericInteractionFactorsupup != 0) 
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AduAu(i, j, j2, Coefficient, NbrTranslations);
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
      
      if (this->OneBodyGenericInteractionFactorsdowndown != 0) 
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAd(i, j, j2, Coefficient, NbrTranslations);
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
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAu(i, j, j2, Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			  Index = particles->AduAd(i, j, j2, Coefficient, NbrTranslations);
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
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AduAu(i, j, j2, Coefficient, NbrTranslations);
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
      
      if (this->OneBodyGenericInteractionFactorsdowndown != 0) 
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAd(i, j, j2, Coefficient, NbrTranslations);
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
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAu(i, j, j2, Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			    ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			  Index = particles->AduAd(i, j, j2, Coefficient, NbrTranslations);
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
