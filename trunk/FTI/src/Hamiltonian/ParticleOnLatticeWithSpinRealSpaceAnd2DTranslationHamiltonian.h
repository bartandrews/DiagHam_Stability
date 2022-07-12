////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                       Class author: Cecile Repellin                        //
//                                                                            //
//        class of generic hamiltonian for interacting spinful particles      //
//         on lattice written in real space and handling 2d translations      //
//                                                                            //
//                        last modification : 11/09/2014                      //
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


#ifndef PARTICLEONLATTICEWITHSPINREALSPACEAND2DTRANSLATIONHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINREALSPACEAND2DTRANSLATIONHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceHamiltonian.h"
#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian : public ParticleOnLatticeWithSpinRealSpaceHamiltonian
{

 protected:

  // momentum along the x direction
  int XMomentum;
  // periodicity in the x direction
  int MaxXMomentum;

  // momentum along the y direction
  int YMomentum;
  // periodicity in the y direction
  int MaxYMomentum;
  
  //array containing all the phase factors that are needed when computing matrix elements
  Complex** ExponentialFactors;


 public:

  // default constructor
  //
  ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSites = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = number of momentum sectors in the x direction
  // yMomentum = momentum sector in the x direction
  // maxYMomentum = number of momentum sectors in the x direction
  // tightBindingupup = hamiltonian corresponding to the tight-binding model in real space for particles with spin up
  // tightBindingdowndown = hamiltonian corresponding to the tight-binding model in real space for particles with spin down
  // densityDensityupup = matrix that gives the amplitude of each density-density interaction term between particles with spin up
  // densityDensitydowndown = matrix that gives the amplitude of each density-density interaction term between particles with spin down
  // densityDensityupdown = matrix that gives the amplitude of each density-density interaction term between particles with spin up and down
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSites, 
								int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, 
								HermitianMatrix& tightBindingupup, HermitianMatrix& tightBindingdowndown,
								RealSymmetricMatrix& densityDensityupup, RealSymmetricMatrix& densityDensitydowndown, 
								RealSymmetricMatrix& densityDensityupdown,
								AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  virtual ~ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian();

  // add an additional S^2 term to the Hamiltonian
  //
  // factor = factor in front of the S^2
  // fixedSz = flag indicating whether Sz needs to be evaluated
  // memory = amount of memory that can be used for S^2  precalculations
  virtual void AddS2 (double factor = 1.0, bool fixedSz = true, long memory = 0l);

  
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

inline void ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  int* TmpIndices2;
  Complex* TmpInteractionFactor;
  int Index;
  int NbrTranslationsX;
  int NbrTranslationsY;
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
      TmpIndices = this->IntraSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupup[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      
	      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	}
    }
  for (int j = 0; j < this->NbrInterSectorSums; ++j)
    {
      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
      TmpIndices2 = this->InterSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  Coefficient3 = particles->AuAd(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      
	      TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor)) * Coefficient4;
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

inline void ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
												      ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  int* TmpIndices;
  int* TmpIndices2;
  Complex* TmpInteractionFactor;
  int NbrTranslationsX;
  int NbrTranslationsY;
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
      TmpIndices = this->IntraSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupup[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      
	      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	    }
	}
    }
    for (int j = 0; j < this->NbrInterSectorSums; ++j)
    {
      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
      TmpIndices2 = this->InterSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  Coefficient3 = particles->AuAd(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      
	      TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor) * tmpCoefficients[p];
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

inline void ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
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
  int NbrTranslationsX;
  int NbrTranslationsY;
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
      TmpIndices = this->IntraSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupup[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      
	      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	}
    }
    for (int j = 0; j < this->NbrInterSectorSums; ++j)
    {
      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
      TmpIndices2 = this->InterSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  Coefficient3 = particles->AuAd(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      
	      TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * Coefficient4;
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

inline void ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
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
  int NbrTranslationsX;
  int NbrTranslationsY;
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
      TmpIndices = this->IntraSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupup[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      
	      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	    }
	}
    }
    for (int j = 0; j < this->NbrInterSectorSums; ++j)
    {
      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
      TmpIndices2 = this->InterSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  Coefficient3 = particles->AuAd(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];

	      TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * tmpCoefficients[p];
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

inline void ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, int* indexArray, Complex* coefficientArray, long& position)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  int* TmpIndices2;
  Complex* TmpInteractionFactor;
  int Index;
  int NbrTranslationsX;
  int NbrTranslationsY;

  if (this->HermitianSymmetryFlag == false)
    {
      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
	  TmpIndices = this->IntraSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupup[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	    }
	}
      for (int j = 0; j < this->NbrInterSectorSums; ++j)
	{
	  int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
	  TmpIndices2 = this->InterSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim2; i1 += 2)
	    {
	      Coefficient3 = particles->AuAd(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor);
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
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupup[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		      if (Index <= index)
			{
			  if (Index == index)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		      if (Index <= index)
			{
			  if (Index == index)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		}
	    }
	}
	for (int j = 0; j < this->NbrInterSectorSums; ++j)
	{
	  int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
	  TmpIndices2 = this->InterSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim2; i1 += 2)
	    {
	      Coefficient3 = particles->AuAd(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		      if (Index <= index)
			{
			  if (Index == index)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * (*TmpInteractionFactor);
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

inline void ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  int* TmpIndices2;
  Complex* TmpInteractionFactor;
  int Index;
  int NbrTranslationsX;
  int NbrTranslationsY;

  if (this->HermitianSymmetryFlag == false)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
	      TmpIndices = this->IntraSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient3 = particles->AuAu(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsupup[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
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
		      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
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
	  for (int j = 0; j < this->NbrInterSectorSums; ++j)
	    {
	      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
	      TmpIndices2 = this->InterSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim2; i1 += 2)
		{
		  Coefficient3 = particles->AuAd(i, TmpIndices2[i1], TmpIndices2[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
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
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient3 = particles->AuAu(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsupup[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
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
		      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
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
	  for (int j = 0; j < this->NbrInterSectorSums; ++j)
	    {
	      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
	      TmpIndices2 = this->InterSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim2; i1 += 2)
		{
		  Coefficient3 = particles->AuAd(i, TmpIndices2[i1], TmpIndices2[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
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

inline void ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
														 int step, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Index;
  double Coefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;
  if (this->OneBodyGenericInteractionFactorsupup != 0)
    {
      if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	{
	  if (this->HermitianSymmetryFlag == false)
	    {
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  Coefficient = 0.0;
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
			      Coefficient = TmpCoefficient;
			      Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			      if (Index < particles->GetHilbertSpaceDimension())
				vDestination[Index] += TmpCoefficient * TmpInteractionFactorsupup[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * vSource[i];
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
			      TmpCoefficient = Coefficient;
			      Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			      if (Index < particles->GetHilbertSpaceDimension())
				vDestination[Index] += TmpCoefficient * TmpInteractionFactorsdowndown[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * vSource[i];
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
		      if (Coefficient != 0.0)
			{
			  double TmpCoefficient;
			  int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
			  int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
			  Complex* TmpInteractionFactorsupup = this->OneBodyGenericInteractionFactorsupup[j];
			  for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			    { TmpCoefficient = Coefficient;
			      Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			      if (Index <= i)
				{
				  vDestination[Index] += TmpCoefficient * TmpInteractionFactorsupup[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * vSource[i];
				  if (Index < i)
				    vDestination[i] += TmpCoefficient * Conj(TmpInteractionFactorsupup[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * vSource[Index];
				}
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
			      TmpCoefficient = Coefficient;
			      Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			      if (Index <= i)
				{
				  vDestination[Index] += TmpCoefficient * TmpInteractionFactorsdowndown[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * vSource[i];
				  if (Index < i)
				    vDestination[i] += TmpCoefficient * Conj(TmpInteractionFactorsdowndown[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * vSource[Index];
				}
			    }
			}
		    }
		}
	    }
	}
      else
	{
	  if (this->HermitianSymmetryFlag == false)
	    {
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  Coefficient = 0.0;
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
			      Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			      if (Index < particles->GetHilbertSpaceDimension())
				vDestination[Index] += TmpCoefficient * TmpInteractionFactorsupup[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * vSource[i];
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
		      if (Coefficient != 0.0)
			{
			  double TmpCoefficient;
			  int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
			  int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
			  Complex* TmpInteractionFactorsupup = this->OneBodyGenericInteractionFactorsupup[j];
			  for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			    {
			      Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			      if (Index <= i)
				{
				  vDestination[Index] += TmpCoefficient * TmpInteractionFactorsupup[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * vSource[i];
				  if (Index < i)
				    vDestination[i] += TmpCoefficient * Conj(TmpInteractionFactorsupup[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * vSource[Index];
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
	  if (this->HermitianSymmetryFlag == false)
	    {
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  Coefficient = 0.0;
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
			      Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			      if (Index < particles->GetHilbertSpaceDimension())
				vDestination[Index] += TmpCoefficient * TmpInteractionFactorsdowndown[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * vSource[i];
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
			      Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			      if (Index <= i)
				{
				  vDestination[Index] += TmpCoefficient * TmpInteractionFactorsdowndown[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * vSource[i];
				  if (Index < i)
				    vDestination[i] += TmpCoefficient * Conj(TmpInteractionFactorsdowndown[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * vSource[Index];
				}
			    }
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

inline void ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
													 int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  int Index;
  int NbrTranslationsX;
  int NbrTranslationsY;
  if (this->OneBodyGenericInteractionFactorsupup != 0) 
    {
      if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	{
	  if (this->HermitianSymmetryFlag == false)
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
				  Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
				  if (Index < Dim)
				    TmpDestinationVector[Index] += TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsupup[k] * TmpSourceVector[i];
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
				  TmpCoefficient = Coefficient;
				  Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
				  if (Index < Dim)
				    TmpDestinationVector[Index] += TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsdowndown[k] * TmpSourceVector[i];
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
			  Coefficient = particles->Au(i, j);
			  if (Coefficient != 0.0)
			    {
			      double TmpCoefficient;
			      int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
			      int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
			      Complex* TmpInteractionFactorsupup = this->OneBodyGenericInteractionFactorsupup[j];
			      for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
				{
				  TmpCoefficient = Coefficient;
				  Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
				  if (Index <= i)
				    {
				      TmpDestinationVector[Index] += TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsupup[k] * TmpSourceVector[i];
				      if (Index < i)
					TmpDestinationVector[i] += TmpCoefficient * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsupup[k]) * TmpSourceVector[Index];			      
				    }
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
				  TmpCoefficient = Coefficient;
				  Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
				  if (Index <= i)
				    {
				      TmpDestinationVector[Index] += TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsdowndown[k] * TmpSourceVector[i];
				      if (Index < i)
					TmpDestinationVector[i] += TmpCoefficient * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsdowndown[k]) * TmpSourceVector[Index];			      
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
	  if (this->HermitianSymmetryFlag == false)
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
				  Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
				  if (Index < Dim)
				    TmpDestinationVector[Index] += TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsupup[k] * TmpSourceVector[i];
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
				  Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
				  if (Index <= i)
				    {
				      TmpDestinationVector[Index] += TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsupup[k] * TmpSourceVector[i];
				      if (Index < i)
					TmpDestinationVector[i] += TmpCoefficient * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsupup[k]) * TmpSourceVector[Index];			      
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
	  if (this->HermitianSymmetryFlag == false)
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
				  Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
				  if (Index < Dim)
				    TmpDestinationVector[Index] += TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsdowndown[k] * TmpSourceVector[i];
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
			  Coefficient = particles->Ad(i,j);
			  if (Coefficient != 0.0)
			    {
			      double TmpCoefficient;
			      int TmpNbrConnectedSitesdowndown = this->OneBodyGenericNbrConnectedSitesdowndown[j];
			      int* TmpConnectedSitesdowndown = this->OneBodyGenericConnectedSitesdowndown[j];
			      Complex* TmpInteractionFactorsdowndown = this->OneBodyGenericInteractionFactorsdowndown[j];
			      for (int k = 0; k < TmpNbrConnectedSitesdowndown; ++k)
				{
				  Index = particles->Add( TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
				  if (Index <= i)
				    {
				      TmpDestinationVector[Index] += TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsdowndown[k] * TmpSourceVector[i];
				      if (Index < i)
					TmpDestinationVector[i] += TmpCoefficient * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsdowndown[k]) * TmpSourceVector[Index];			      
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

inline void ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
															       int* indexArray, Complex* coefficientArray, long& position)
{
  if ((this->OneBodyGenericInteractionFactorsupup == 0) && (this->OneBodyGenericInteractionFactorsdowndown == 0))
    return;
  
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  int Index;
  int NbrTranslationsX;
  int NbrTranslationsY;
  if (this->OneBodyGenericInteractionFactorsupup != 0)
    {
      if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	{
	  if (this->HermitianSymmetryFlag == false)
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
			  Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			  if (Index < Dim)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsupup[k];
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
			  Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			  if (Index < Dim)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsdowndown[k];
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
		  Coefficient = particles->Au(index,j);
		  if (Coefficient != 0.0)
		    {
		      double TmpCoefficient;
		      int TmpNbrConnectedSitesupup = this->OneBodyGenericNbrConnectedSitesupup[j];
		      int* TmpConnectedSitesupup = this->OneBodyGenericConnectedSitesupup[j];
		      Complex* TmpInteractionFactorsupup = this->OneBodyGenericInteractionFactorsupup[j];
		      for (int k = 0; k < TmpNbrConnectedSitesupup; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			  if (Index <= index)
			    {
			      if (Index == index)
				{
				  indexArray[position] = Index;
				  coefficientArray[position] = 0.5 * TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsupup[k];
				  ++position;
				}
			      else
				{
				  indexArray[position] = Index;
				  coefficientArray[position] = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsupup[k];
				  ++position; 
				}
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
			  Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			  if (Index <= index)
			    {
			      if (Index == index)
				{
				  indexArray[position] = Index;
				  coefficientArray[position] = 0.5 * TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsdowndown[k];
				  ++position;
				}
			      else
				{
				  indexArray[position] = Index;
				  coefficientArray[position] = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsdowndown[k];
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
	  if (this->HermitianSymmetryFlag == false)
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
			  Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			  if (Index < Dim)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsupup[k];
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
			  Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			  if (Index <= index)
			    {
			      indexArray[position] = Index;
			      if (Index == index)
				{
				  coefficientArray[position] = 0.5 * TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsupup[k];
				}
			      else
				{
				  coefficientArray[position] = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsupup[k];
				}
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
      if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	{
	  if (this->HermitianSymmetryFlag == false)
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
			  Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			  if (Index < Dim)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsdowndown[k];
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
			  Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			  if (Index <= index)
			    {
			      indexArray[position] = Index;
			      if (Index == index)
				{
				  coefficientArray[position] = 0.5 * TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsdowndown[k];
				}
			      else
				{
			  coefficientArray[position] = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactorsdowndown[k];
				}
			      ++position;
			    }
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

inline void ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  if ((this->OneBodyGenericInteractionFactorsupup == 0) && (this->OneBodyGenericInteractionFactorsdowndown == 0))
    return;
  int Index;
  double Coefficient = 0.0;
  int NbrTranslationsX;
  int NbrTranslationsY;
  int Dim = particles->GetHilbertSpaceDimension();
  
  if (this->OneBodyGenericInteractionFactorsupup != 0)
    {
      if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	{
	  if (this->HermitianSymmetryFlag == false)
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
			      Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
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
			      Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
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
			      Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
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
			      Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
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
	  if (this->HermitianSymmetryFlag == false)
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
			      Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
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
			      TmpCoefficient= Coefficient;
			      Index = particles->Adu(TmpConnectedSitesupup[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
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
  else
    {
      if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	{
	  if (this->HermitianSymmetryFlag == false)
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
			      Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
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
		      Coefficient = particles->Ad(i,j);
		      if (Coefficient != 0.0)
			{
			  double TmpCoefficient;
			  int TmpNbrConnectedSitesdowndown = this->OneBodyGenericNbrConnectedSitesdowndown[j];
			  int* TmpConnectedSitesdowndown = this->OneBodyGenericConnectedSitesdowndown[j];
			  for (int k = 0; k < TmpNbrConnectedSitesdowndown; ++k)
			    {
			      TmpCoefficient = Coefficient;
			      Index = particles->Add(TmpConnectedSitesdowndown[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
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
}


#endif
