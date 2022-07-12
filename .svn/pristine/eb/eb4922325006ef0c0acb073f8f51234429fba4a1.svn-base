////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//            to particles on a torus with magnetic translations and          //
//                      an internal SU(3) degree of freedom                   //
//                                                                            //
//                        last modification : 08/08/2014                      //
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


#ifndef ABSTRACTQHEONTORUSWITHSU3SPINANDMAGNETICTRANSLATIONSHAMILTONIAN_H
#define ABSTRACTQHEONTORUSWITHSU3SPINANDMAGNETICTRANSLATIONSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSU3SpinAndMagneticTranslations.h"
#include "Hamiltonian/ParticleOnLatticeWithSU3SpinChernInsulatorHamiltonian.h"


#include <iostream>


using std::ostream;


class AbstractArchitecture;


class AbstractQHEOnTorusWithSU3SpinAndMagneticTranslationsHamiltonian : public ParticleOnLatticeWithSU3SpinChernInsulatorHamiltonian
{

 protected:
  
  // ratio between the width in the x direction and the width in the y direction
  double Ratio;
  // ratio between the width in the y direction and the width in the x direction
  double InvRatio;

  // maximum momentum value reached by a particle in the state
  int MaxMomentum;
  // momentum value in the x direction (modulo GCD of NbParticles and MaxMomentum)
  int XMomentum;
  // number of Lz values in a state
  int NbrLzValue;
  // GCD of MaxMomentum and NbrFermions (momemta are defined modulo MomentumModulo)
  int MomentumModulo;

  //array containing all the phase factors that are needed when computing matrix elements
  Complex* ExponentialFactors;

 public:

  // destructor
  //
  virtual ~AbstractQHEOnTorusWithSU3SpinAndMagneticTranslationsHamiltonian();
  
 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors() = 0;

  // evaluate all exponential factors
  //   
  virtual void EvaluateExponentialFactors();

  // get all the indices that should appear in the annihilation/creation operators
  //
  virtual void GetIndices();

  // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the two-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int index, ComplexVector* vSources, 
						     ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the two-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  inline void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int index, ComplexVector* vSources, 
							     ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // core part of the FastMultiplication method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray  
  virtual void EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSU3Spin* particles, int index, 
							    int* indexArray, Complex* coefficientArray, long& position);

  // core part of the PartialFastMultiplicationMemory method involving two-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations  
  virtual void EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSU3Spin* particles, int firstComponent, int lastComponent, long& memory);

  // core part of the AddMultiply method involving the one-body interaction, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int firstComponent, int lastComponent,
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
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int firstComponent, int lastComponent,
						     int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors);

  // core part of the FastMultiplication method involving the one-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSU3Spin* particles, int index, 
							    int* indexArray, Complex* coefficientArray, long& position);

  // core part of the PartialFastMultiplicationMemory method involving two-body term and one-body terms
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSU3Spin* particles, int firstComponent, int lastComponent, long& memory);

};

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnTorusWithSU3SpinAndMagneticTranslationsHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int Index;
  int NbrTranslations;
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
      TmpIndices = this->IntraSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->A1A1(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors11[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A2A2(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors22[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A3A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors33[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		  ++TmpInteractionFactor;
		}
	    }
	}
    }
  for (int j = 0; j < this->NbrInterSectorSums; ++j)
    {
      int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
      TmpIndices = this->InterSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->A1A2(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors12[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A1A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors13[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A2A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors23[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
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

inline void AbstractQHEOnTorusWithSU3SpinAndMagneticTranslationsHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int index, ComplexVector* vSources, 
														   ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  int NbrTranslations;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
      TmpIndices = this->IntraSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->A1A1(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors11[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A2A2(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors22[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A3A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors33[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	    }
	}
    }
  for (int j = 0; j < this->NbrInterSectorSums; ++j)
    {
      int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
      TmpIndices = this->InterSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->A1A2(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors12[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A1A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors13[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A2A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors23[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
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

inline void AbstractQHEOnTorusWithSU3SpinAndMagneticTranslationsHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  //int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  Complex TmpSum = 0.0;
  int Index;
  int NbrTranslations;
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
      TmpIndices = this->IntraSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->A1A1(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors11[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A2A2(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors22[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A3A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors33[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	}
    }
  for (int j = 0; j < this->NbrInterSectorSums; ++j)
    {
      int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
      TmpIndices = this->InterSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->A1A2(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors12[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A1A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors13[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A2A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors23[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
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

inline void AbstractQHEOnTorusWithSU3SpinAndMagneticTranslationsHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int index, ComplexVector* vSources, 
													       ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  //int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  int NbrTranslations;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  Complex* TmpSum = new Complex[nbrVectors];
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
      TmpIndices = this->IntraSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->A1A1(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors11[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj(this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A2A2(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors22[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj(this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A3A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors33[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj(this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
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
      int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
      TmpIndices = this->InterSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->A1A2(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors12[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj(this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A1A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors13[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj(this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A2A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors23[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj(this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
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

inline void AbstractQHEOnTorusWithSU3SpinAndMagneticTranslationsHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSU3Spin* particles, int index, 
															 int* indexArray, Complex* coefficientArray, long& position)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int Dim = particles->GetHilbertSpaceDimension();
  int NbrTranslations;

  if (this->HermitianSymmetryFlag == false)
    {
      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
	  TmpIndices = this->IntraSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient2 = particles->A1A1(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors11[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
                          indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient2 = particles->A2A2(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors22[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
                          indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient2 = particles->A3A3(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors33[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
                          indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	    }
	}
      for (int j = 0; j < this->NbrInterSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
	  TmpIndices = this->InterSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient2 = particles->A1A2(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors12[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
                          indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient2 = particles->A1A3(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors13[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
                          indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient2 = particles->A2A3(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors23[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
                          indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
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
      int AbsoluteIndex = index + this->PrecalculationShift;
      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
	  TmpIndices = this->IntraSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient2 = particles->A1A1(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors11[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient2 = particles->A2A2(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors22[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient2 = particles->A3A3(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors33[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
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
	  int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
	  TmpIndices = this->InterSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient2 = particles->A1A2(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors12[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient2 = particles->A1A3(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors13[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient2 = particles->A2A3(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors23[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
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

// core part of the PartialFastMultiplicationMemory method involving two-body term and one-body terms
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void AbstractQHEOnTorusWithSU3SpinAndMagneticTranslationsHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSU3Spin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  // double* TmpInteractionFactor;
  int NbrTranslations;
  int Dim = particles->GetHilbertSpaceDimension();
  
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
		  Coefficient2 = particles->A1A1(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			}
		    }
		  Coefficient2 = particles->A2A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			}
		    }
		  Coefficient2 = particles->A3A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			}
		    }
		}
	    }
	  
	  for (int j = 0; j < this->NbrInterSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
	      TmpIndices = this->InterSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient2 = particles->A1A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			}
		    }
		  Coefficient2 = particles->A1A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			}
		    }
		  Coefficient2 = particles->A2A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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
	  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
	      TmpIndices = this->IntraSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient2 = particles->A1A1(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			}
		    }
		  Coefficient2 = particles->A2A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			}
		    }
		  Coefficient2 = particles->A3A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			}
		    }
		}
	    }
	  for (int j = 0; j < this->NbrInterSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
	      TmpIndices = this->InterSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient2 = particles->A1A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			}
		    }
		  Coefficient2 = particles->A1A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			}
		    }
		  Coefficient2 = particles->A2A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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

// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnTorusWithSU3SpinAndMagneticTranslationsHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int firstComponent, int lastComponent,
														  int step, ComplexVector& vSource, ComplexVector& vDestination)
{
  if (this->OneBodyInteractionFactors11 != 0)
    if (this->OneBodyInteractionFactors22 != 0)
      {
	if (this->OneBodyInteractionFactors33 != 0)
	  {
	    double TmpDiagonal = 0.0;
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		TmpDiagonal = 0.0;
		for (int j = 0; j <= this->LzMax; ++j) 
		  {
		    TmpDiagonal += this->OneBodyInteractionFactors11[j] * particles->Ad1A1(i, j);
		    TmpDiagonal += this->OneBodyInteractionFactors22[j] * particles->Ad2A2(i, j);
		    TmpDiagonal += this->OneBodyInteractionFactors33[j] * particles->Ad3A3(i, j);
		  }
		vDestination[i] += (this->HamiltonianShift + TmpDiagonal)* vSource[i];
	      }
	  }
	else
	  {
	    double TmpDiagonal = 0.0;
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		TmpDiagonal = 0.0;
		for (int j = 0; j <= this->LzMax; ++j) 
		  {
		    TmpDiagonal += this->OneBodyInteractionFactors11[j] * particles->Ad1A1(i, j);
		    TmpDiagonal += this->OneBodyInteractionFactors22[j] * particles->Ad2A2(i, j);
		  }
		vDestination[i] += (this->HamiltonianShift + TmpDiagonal)* vSource[i];
	      }
	  }
      }
    else
      {
	if (this->OneBodyInteractionFactors33 != 0)
	  {
	    double TmpDiagonal = 0.0;
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		TmpDiagonal = 0.0;
		for (int j = 0; j <= this->LzMax; ++j) 
		  {
		    TmpDiagonal += this->OneBodyInteractionFactors11[j] * particles->Ad1A1(i, j);
		    TmpDiagonal += this->OneBodyInteractionFactors33[j] * particles->Ad3A3(i, j);
		  }
		vDestination[i] += (this->HamiltonianShift + TmpDiagonal)* vSource[i];
	      }
	  }
	else
	  {
	    double TmpDiagonal = 0.0;
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		TmpDiagonal = 0.0;
		for (int j = 0; j <= this->LzMax; ++j) 
		  TmpDiagonal += this->OneBodyInteractionFactors11[j] * particles->Ad1A1(i, j);
		vDestination[i] += (this->HamiltonianShift + TmpDiagonal)* vSource[i];
	      }
	  }
      }
  else
    {
      if (this->OneBodyInteractionFactors22 != 0)
	{
	  if (this->OneBodyInteractionFactors33 != 0)
	    {
	      double TmpDiagonal = 0.0;
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  TmpDiagonal = 0.0;
		  for (int j = 0; j <= this->LzMax; ++j) 
		    {
		      TmpDiagonal += this->OneBodyInteractionFactors22[j] * particles->Ad2A2(i, j);
		      TmpDiagonal += this->OneBodyInteractionFactors33[j] * particles->Ad3A3(i, j);
		    }
		  vDestination[i] += (this->HamiltonianShift + TmpDiagonal)* vSource[i];
		}
	    }
	  else
	    {
	      double TmpDiagonal = 0.0;
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  TmpDiagonal = 0.0;
		  for (int j = 0; j <= this->LzMax; ++j) 
		    TmpDiagonal += this->OneBodyInteractionFactors22[j] * particles->Ad2A2(i, j);
		  vDestination[i] += (this->HamiltonianShift + TmpDiagonal)* vSource[i];
		}
	    }
	}	
      else
	{
	  if (this->OneBodyInteractionFactors33 != 0)
	    {
	      double TmpDiagonal = 0.0;
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  TmpDiagonal = 0.0;
		  for (int j = 0; j <= this->LzMax; ++j) 
		    TmpDiagonal += this->OneBodyInteractionFactors33[j] * particles->Ad3A3(i, j);
		  vDestination[i] += (this->HamiltonianShift + TmpDiagonal)* vSource[i];
		}
	    }
	  else
	    {
	      for (int i = firstComponent; i < lastComponent; i += step)
		vDestination[i] += this->HamiltonianShift * vSource[i];
	    }
	}
    }
  if (this->OneBodyInteractionFactors12 != 0)
    {
      double Coefficient;
      Complex Source;
      int NbrTranslations;
      int Dim = particles->GetHilbertSpaceDimension();
      int Index;
      for (int i = firstComponent; i < lastComponent; i += step)
	{
	  Source = vSource[i];
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->Ad2A1(i + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(OneBodyInteractionFactors12[j])) * Source;
		}
	      Index = particles->Ad1A2(i + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * OneBodyInteractionFactors12[j]) * Source;
		}
	    }
	}
    }
  if (this->OneBodyInteractionFactors13 != 0)
    {
      double Coefficient;
      Complex Source;
      int NbrTranslations;
      int Dim = particles->GetHilbertSpaceDimension();
      int Index;
      for (int i = firstComponent; i < lastComponent; i += step)
	{
	  Source = vSource[i];
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->Ad3A1(i + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(OneBodyInteractionFactors13[j])) * Source;
		}
	      Index = particles->Ad1A3(i + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * OneBodyInteractionFactors13[j]) * Source;
		}
	    }
	}
    }
  if (this->OneBodyInteractionFactors23 != 0)
    {
      double Coefficient;
      Complex Source;
      int NbrTranslations;
      int Dim = particles->GetHilbertSpaceDimension();
      int Index;
      for (int i = firstComponent; i < lastComponent; i += step)
	{
	  Source = vSource[i];
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->Ad3A2(i + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(OneBodyInteractionFactors23[j])) * Source;
		}
	      Index = particles->Ad2A3(i + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * OneBodyInteractionFactors23[j]) * Source;
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

inline void AbstractQHEOnTorusWithSU3SpinAndMagneticTranslationsHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int firstComponent, int lastComponent,
														  int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  if (this->OneBodyInteractionFactors11 != 0) 
    if (this->OneBodyInteractionFactors22 != 0)
      {
	double TmpDiagonal = 0.0;
	for (int p = 0; p < nbrVectors; ++p)
	  {
	    ComplexVector& TmpSourceVector = vSources[p];
	    ComplexVector& TmpDestinationVector = vDestinations[p];
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		TmpDiagonal = 0.0;
		for (int j = 0; j <= this->LzMax; ++j) 
		  {
		    TmpDiagonal += this->OneBodyInteractionFactors11[j] * particles->Ad1A1(i, j);
		    TmpDiagonal += this->OneBodyInteractionFactors22[j] * particles->Ad2A2(i, j);
		  }
		TmpDestinationVector[i] += (this->HamiltonianShift + TmpDiagonal)* TmpSourceVector[i];
	      }
	  }
      }
    else
      {
	double TmpDiagonal = 0.0;
	for (int p = 0; p < nbrVectors; ++p)
	  {
	    ComplexVector& TmpSourceVector = vSources[p];
	    ComplexVector& TmpDestinationVector = vDestinations[p];
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		TmpDiagonal = 0.0;
		for (int j = 0; j <= this->LzMax; ++j) 
		  TmpDiagonal += this->OneBodyInteractionFactors11[j] * particles->Ad1A1(i, j);
		TmpDestinationVector[i] += (this->HamiltonianShift + TmpDiagonal)* TmpSourceVector[i];
	      }
	  }
      }
  else
    if (this->OneBodyInteractionFactors22 != 0)
      {
	double TmpDiagonal = 0.0;
	for (int p = 0; p < nbrVectors; ++p)
	  {
	    ComplexVector& TmpSourceVector = vSources[p];
	    ComplexVector& TmpDestinationVector = vDestinations[p];
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		TmpDiagonal = 0.0;
		for (int j = 0; j <= this->LzMax; ++j) 
		  TmpDiagonal += this->OneBodyInteractionFactors22[j] * particles->Ad2A2(i, j);
		TmpDestinationVector[i] += (this->HamiltonianShift + TmpDiagonal)* TmpSourceVector[i];
	      }
	  }
      }	
    else
      for (int p = 0; p < nbrVectors; ++p)
	{
	  ComplexVector& TmpSourceVector = vSources[p];
	  ComplexVector& TmpDestinationVector = vDestinations[p];
	  for (int i = firstComponent; i < lastComponent; i += step)
	    TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	}
  for (int p = 0; p < nbrVectors; ++p)
    {
      ComplexVector& TmpSourceVector = vSources[p];
      ComplexVector& TmpDestinationVector = vDestinations[p];
      for (int i = firstComponent; i < lastComponent; i += step)
	TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
    }
  if (this->OneBodyInteractionFactors12 != 0)
    {
      double Coefficient;
      int Dim = particles->GetHilbertSpaceDimension();
      int Index;
      int NbrTranslations;
      for (int i = firstComponent; i < lastComponent; i += step)
	{
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->Ad2A1(i + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  for (int p = 0; p < nbrVectors; ++p)
		    vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(OneBodyInteractionFactors12[j]) * vSources[p][i];
		}
	      Index = particles->Ad1A2(i + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  for (int p = 0; p < nbrVectors; ++p)
		    vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * OneBodyInteractionFactors12[j] * vSources[p][i];
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

inline void AbstractQHEOnTorusWithSU3SpinAndMagneticTranslationsHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSU3Spin* particles, int index, 
															 int* indexArray, Complex* coefficientArray, long& position)
{
  if ((this->OneBodyInteractionFactors22 != 0) || (this->OneBodyInteractionFactors11 != 0) || (this->OneBodyInteractionFactors33 != 0))
    {
      double TmpDiagonal = 0.0;
      if (this->OneBodyInteractionFactors11 != 0)
	for (int j = 0; j <= this->LzMax; ++j)
	  TmpDiagonal += this->OneBodyInteractionFactors11[j] * particles->Ad1A1(index + this->PrecalculationShift, j);
      if (this->OneBodyInteractionFactors22 != 0)
	for (int j = 0; j <= this->LzMax; ++j)
	  TmpDiagonal += this->OneBodyInteractionFactors22[j] * particles->Ad2A2(index + this->PrecalculationShift, j);	  
      if (this->OneBodyInteractionFactors33 != 0)
	for (int j = 0; j <= this->LzMax; ++j)
	  TmpDiagonal += this->OneBodyInteractionFactors33[j] * particles->Ad3A3(index + this->PrecalculationShift, j);	  
      indexArray[position] = index + this->PrecalculationShift;
      if (this->HermitianSymmetryFlag == true)
	TmpDiagonal *= 0.5;
      coefficientArray[position] = TmpDiagonal;
      ++position;
    }
  if (this->OneBodyInteractionFactors12 != 0)
    {
      if (this->HermitianSymmetryFlag == false)
	{
	  int NbrTranslations;
	  int Dim = particles->GetHilbertSpaceDimension();
	  double Coefficient;
	  int Index;
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->Ad2A1(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyInteractionFactors12[j];
		  ++position;
		}
	      Index = particles->Ad1A2(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyInteractionFactors12[j]);
		  ++position;
		}
	    }
	}
      else
	{
	  int NbrTranslations;
	  int Dim = particles->GetHilbertSpaceDimension();
	  double Coefficient;
	  int Index;
	  int AbsoluteIndex = index + this->PrecalculationShift;
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->Ad2A1(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index <= AbsoluteIndex)
		{
		  if (Index == AbsoluteIndex)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = 0.5 * Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyInteractionFactors12[j];
		      ++position;
		    }
		  else
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyInteractionFactors12[j];
		      ++position;
		    }
		}
	      Index = particles->Ad1A2(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index <= AbsoluteIndex)
		{
		  if (Index == AbsoluteIndex)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = 0.5 * Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyInteractionFactors12[j]);
		      ++position;
		    }
		  else
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyInteractionFactors12[j]);
		      ++position;
		    }
		}
	    }
	}
    }
  if (this->OneBodyInteractionFactors13 != 0)
    {
      if (this->HermitianSymmetryFlag == false)
	{
	  int NbrTranslations;
	  int Dim = particles->GetHilbertSpaceDimension();
	  double Coefficient;
	  int Index;
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->Ad3A1(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyInteractionFactors13[j];
		  ++position;
		}
	      Index = particles->Ad1A3(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyInteractionFactors13[j]);
		  ++position;
		}
	    }
	}
      else
	{
	  int NbrTranslations;
	  int Dim = particles->GetHilbertSpaceDimension();
	  double Coefficient;
	  int Index;
	  int AbsoluteIndex = index + this->PrecalculationShift;
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->Ad3A1(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index <= AbsoluteIndex)
		{
		  if (Index == AbsoluteIndex)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = 0.5 * Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyInteractionFactors13[j];
		      ++position;
		    }
		  else
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyInteractionFactors13[j];
		      ++position;
		    }
		}
	      Index = particles->Ad1A3(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index <= AbsoluteIndex)
		{
		  if (Index == AbsoluteIndex)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = 0.5 * Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyInteractionFactors13[j]);
		      ++position;
		    }
		  else
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyInteractionFactors13[j]);
		      ++position;
		    }
		}
	    }
	}
    }
  if (this->OneBodyInteractionFactors23 != 0)
    {
      if (this->HermitianSymmetryFlag == false)
	{
	  int NbrTranslations;
	  int Dim = particles->GetHilbertSpaceDimension();
	  double Coefficient;
	  int Index;
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->Ad3A2(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyInteractionFactors23[j];
		  ++position;
		}
	      Index = particles->Ad2A3(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyInteractionFactors23[j]);
		  ++position;
		}
	    }
	}
      else
	{
	  int NbrTranslations;
	  int Dim = particles->GetHilbertSpaceDimension();
	  double Coefficient;
	  int Index;
	  int AbsoluteIndex = index + this->PrecalculationShift;
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->Ad3A2(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index <= AbsoluteIndex)
		{
		  if (Index == AbsoluteIndex)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = 0.5 * Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyInteractionFactors23[j];
		      ++position;
		    }
		  else
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyInteractionFactors23[j];
		      ++position;
		    }
		}
	      Index = particles->Ad2A3(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index <= AbsoluteIndex)
		{
		  if (Index == AbsoluteIndex)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = 0.5 * Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyInteractionFactors23[j]);
		      ++position;
		    }
		  else
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyInteractionFactors23[j]);
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

inline void AbstractQHEOnTorusWithSU3SpinAndMagneticTranslationsHamiltonian::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSU3Spin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  int NbrTranslations;
  int Dim = particles->GetHilbertSpaceDimension();
  
  if (this->HermitianSymmetryFlag == false)
    {
      if (this->OneBodyInteractionFactors12 != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j<= this->LzMax; ++j)
		{
		  Index = particles->Ad2A1(i, j, j, Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		  Index = particles->Ad1A2(i, j, j, Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		}
	    }
	}
      if (this->OneBodyInteractionFactors13 != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j<= this->LzMax; ++j)
		{
		  Index = particles->Ad3A1(i, j, j, Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		  Index = particles->Ad1A3(i, j, j, Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		}
	    }
	}
      if (this->OneBodyInteractionFactors23 != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j<= this->LzMax; ++j)
		{
		  Index = particles->Ad3A2(i, j, j, Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		  Index = particles->Ad2A3(i, j, j, Coefficient, NbrTranslations);
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
      if (this->OneBodyInteractionFactors12 != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j<= this->LzMax; ++j)
		{
		  Index = particles->Ad2A1(i, j, j, Coefficient, NbrTranslations);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		  Index = particles->Ad1A2(i, j, j, Coefficient, NbrTranslations);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		}
	    }
	}
      if (this->OneBodyInteractionFactors13 != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j<= this->LzMax; ++j)
		{
		  Index = particles->Ad3A1(i, j, j, Coefficient, NbrTranslations);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		  Index = particles->Ad1A3(i, j, j, Coefficient, NbrTranslations);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		}
	    }
	}
      if (this->OneBodyInteractionFactors23 != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j<= this->LzMax; ++j)
		{
		  Index = particles->Ad3A2(i, j, j, Coefficient, NbrTranslations);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		  Index = particles->Ad2A3(i, j, j, Coefficient, NbrTranslations);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		}
	    }
	}
    }

  if ((this->OneBodyInteractionFactors22 != 0) || (this->OneBodyInteractionFactors11 != 0) || (this->OneBodyInteractionFactors33 != 0))
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  ++memory;
	  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
	}
    }
}

#endif
