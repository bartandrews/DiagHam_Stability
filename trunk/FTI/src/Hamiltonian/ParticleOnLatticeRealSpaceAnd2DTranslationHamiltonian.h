////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//        class of generic hamiltonian for interacting spinless particles     //
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


#ifndef PARTICLEONLATTICEREALSPACEAND2DTRANSLATIONHAMILTONIAN_H
#define PARTICLEONLATTICEREALSPACEAND2DTRANSLATIONHAMILTONIAN_H


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


class ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian : public ParticleOnLatticeRealSpaceHamiltonian
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
  ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSites = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = number of momentum sectors in the x direction
  // yMomentum = momentum sector in the x direction
  // maxYMomentum = number of momentum sectors in the x direction
  // tightBinding = hamiltonian corresponding to the tight-binding model in real space
  // densityDensity = matrix that gives the amplitude of each density-density interaction term
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSites, 
							int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, 
							HermitianMatrix& tightBinding, RealSymmetricMatrix& densityDensity,
							AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  virtual ~ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian();

 protected:
 
  // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the two-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
						     ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the two-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
							     ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // core part of the FastMultiplication method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray  
  virtual void EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
							    int* indexArray, Complex* coefficientArray, long& position);

  // core part of the PartialFastMultiplicationMemory method involving two-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations  
  virtual void EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory);

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

inline void ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
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
  for (int j = 0; j < this->NbrSectorSums; ++j)
    {
      int Lim = 2 * this->NbrSectorIndicesPerSum[j];
      TmpIndices = this->SectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AA(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
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

inline void ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
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
  for (int j = 0; j < this->NbrSectorSums; ++j)
    {
      int Lim = 2 * this->NbrSectorIndicesPerSum[j];
      TmpIndices = this->SectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AA(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
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

inline void ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
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
  for (int j = 0; j < this->NbrSectorSums; ++j)
    {
      int Lim = 2 * this->NbrSectorIndicesPerSum[j];
      TmpIndices = this->SectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AA(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
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

inline void ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
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
  for (int j = 0; j < this->NbrSectorSums; ++j)
    {
      int Lim = 2 * this->NbrSectorIndicesPerSum[j];
      TmpIndices = this->SectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AA(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor)) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
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

inline void ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
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
  int NbrTranslationsX;
  int NbrTranslationsY;

  if (this->HermitianSymmetryFlag == false)
    {
      for (int j = 0; j < this->NbrSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrSectorIndicesPerSum[j];
	  TmpIndices = this->SectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->AA(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
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
      for (int j = 0; j < this->NbrSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrSectorIndicesPerSum[j];
	  TmpIndices = this->SectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->AA(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
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

inline void ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory)
{
//  cout <<"inline void ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory)"<<endl;
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
	  for (int j = 0; j < this->NbrSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrSectorIndicesPerSum[j];
	      TmpIndices = this->SectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient3 = particles->AA(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
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
	  for (int j = 0; j < this->NbrSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrSectorIndicesPerSum[j];
	      TmpIndices = this->SectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient3 = particles->AA(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslationsX, NbrTranslationsY);
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

inline void ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent,
															int step, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Index;
  double Coefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;
  if (this->OneBodyGenericInteractionFactors != 0)
    {
      if (this->HermitianSymmetryFlag == false)
	{
	  for (int i = firstComponent; i < lastComponent; i += step)
	    { 
	      Coefficient = 0.0;
	      for (int j = 0; j < this->NbrSites; ++j)
		{
		  int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
		  int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
		  Complex* TmpInteractionFactors = this->OneBodyGenericInteractionFactors[j];
		  for (int k = 0; k < TmpNbrConnectedSites; ++k)
		    {
		      Index = particles->AdA(i, j, TmpConnectedSites[k], Coefficient, NbrTranslationsX, NbrTranslationsY);
		      if (Index < particles->GetHilbertSpaceDimension())
			vDestination[Index] += Coefficient * TmpInteractionFactors[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * vSource[i];
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
		  int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
		  int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
		  Complex* TmpInteractionFactors = this->OneBodyGenericInteractionFactors[j];
		  for (int k = 0; k < TmpNbrConnectedSites; ++k)
		    {
		      Index = particles->AdA(i, j, TmpConnectedSites[k], Coefficient, NbrTranslationsX, NbrTranslationsY);
		      if (Index <= i)
			{
			  vDestination[Index] += Coefficient * TmpInteractionFactors[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * vSource[i];
			  if (Index < i)
			    vDestination[i] += Coefficient * Conj(TmpInteractionFactors[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * vSource[Index];
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

inline void ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent,
													 int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  int Index;
  int NbrTranslationsX;
  int NbrTranslationsY;
  if (this->OneBodyGenericInteractionFactors != 0) 
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
		      int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
		      int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
		      Complex* TmpInteractionFactors = this->OneBodyGenericInteractionFactors[j];
		      Coefficient = particles->A(i, j);
		      double TmpCoefficient;
		      if (Coefficient != 0.0)
		      {
			for (int k = 0; k < TmpNbrConnectedSites; ++k)
			  {
			    TmpCoefficient = Coefficient;
			    Index = particles->Ad (TmpConnectedSites[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			    if (Index < Dim)
			      TmpDestinationVector[Index] += TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactors[k] * TmpSourceVector[i];
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
		      int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
		      int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
		      Complex* TmpInteractionFactors = this->OneBodyGenericInteractionFactors[j];
		      Coefficient = particles->A (i,j);
		      double TmpCoefficient;
		      if (Coefficient != 0)
		      {
			for (int k = 0; k < TmpNbrConnectedSites; ++k)
			  {
			    TmpCoefficient = Coefficient;
			    Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
			    if (Index <= i)
			      {
				TmpDestinationVector[Index] += TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactors[k] * TmpSourceVector[i];
				if (Index < i)
				  TmpDestinationVector[i] += TmpCoefficient * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactors[k]) * TmpSourceVector[Index];			      
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

// core part of the FastMultiplication method involving the one-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
															       int* indexArray, Complex* coefficientArray, long& position)
{

//  cout <<"inline void ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 													       int* indexArray, Complex* coefficientArray, long& position)"<<endl;
  if (this->OneBodyGenericInteractionFactors == 0)
    return;

  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  int Index;
  int NbrTranslationsX;
  int NbrTranslationsY;

  if (this->HermitianSymmetryFlag == false)
    {
      for (int j = 0; j < this->NbrSites; ++j)
	{
	  int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
	  int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
	  Complex* TmpInteractionFactors = this->OneBodyGenericInteractionFactors[j];
	  Coefficient = particles->A (index, j);
	  double TmpCoefficient;
	  if (Coefficient != 0.0)
	  {
	    for (int k = 0; k < TmpNbrConnectedSites; ++k)
	      {
		TmpCoefficient = Coefficient;
		Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (Index < Dim)
		  {
		    indexArray[position] = Index;
		    coefficientArray[position] = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactors[k];
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
	  int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
	  int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
	  Complex* TmpInteractionFactors = this->OneBodyGenericInteractionFactors[j];
	  Coefficient = particles->A (index, j);
	  double TmpCoefficient;
	  if (Coefficient != 0.0)
	  {
	    for (int k = 0; k < TmpNbrConnectedSites; ++k)
	      {
		TmpCoefficient = Coefficient;
		Index = particles->Ad (TmpConnectedSites[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (Index <= index)
		  {
		    indexArray[position] = Index;
		    if (Index == index)
		      {
			coefficientArray[position] = 0.5 * TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactors[k];
		      }
		    else
		      {
			coefficientArray[position] = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpInteractionFactors[k];
		      }
		    ++position;
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

inline void ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory)
{
//   cout <<"inline void ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory)"<<endl;
  if (this->OneBodyGenericInteractionFactors == 0) 
 {
//cout <<" in   if (this->OneBodyGenericInteractionFactors == 0)  "<<endl;
    return;
}
  int Index;
  double Coefficient = 0.0;
  int NbrTranslationsX;
  int NbrTranslationsY;
  int Dim = particles->GetHilbertSpaceDimension();
  
  if (this->HermitianSymmetryFlag == false)
    {
  //  cout <<"  if (this->HermitianSymmetryFlag == false)"<<endl;

      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrSites; ++j)
	    {
	      int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
//              cout <<j<< " "<< TmpNbrConnectedSites<<endl;
	      int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
	      double TmpCoefficient ;
              for (int k = 0; k < TmpNbrConnectedSites; ++k)
		  {
		    TmpCoefficient = 1.0;
//                    cout <<k<<" "<< TmpConnectedSites[k] <<endl;
		    Index = particles->AdA(i, j, TmpConnectedSites[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		    if ((Index < Dim) &&(TmpCoefficient != 0.0 ))
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
		    Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
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


#endif
