////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                   class of quatum spin Hall restricted to two band         //
//                                                                            //
//                        last modification : 27/02/2011                      //
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


#ifndef PARTICLEONLATTICEQUANTUMSPINHALLTWOBANDHAMILTONIAN_H
#define PARTICLEONLATTICEQUANTUMSPINHALLTWOBANDHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinChernInsulatorHamiltonian.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class Abstract2DTightBindingModel;


class ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian : public ParticleOnLatticeWithSpinChernInsulatorHamiltonian
{

 protected:
  
  // pointer to the tight binding model
  Abstract2DTightBindingModel* TightBindingModel;

  // szSymmetryBreaking = amplitude of the Sz symmetry breaking term
  double SzSymmetryBreaking;

  // interaction factors : the table first entry is a the sum of annihilation/creation indices, the second entry linearized index (annihilation_index * nbr_creation_index + creation_index)
  // array containing all interaction factors for spin up and spin up (a^+_u a^+_u a_u a_u)
  Complex** InteractionFactorsupupupup;
  // array containing all interaction factors for spin up and spin up (a^+_u a^+_u a_u a_d)
  Complex** InteractionFactorsupupupdown;
  // array containing all interaction factors for spin up and spin up (a^+_d a^+_d a_u a_u)
  Complex** InteractionFactorsupupdowndown;
  // array containing all interaction factors for spin down and spin down (a^+_d a^+_d a_u a_u)
  Complex** InteractionFactorsdowndownupup;
  // array containing all interaction factors for spin down and spin down (a^+_d a^+_d a_u a_d)
  Complex** InteractionFactorsdowndownupdown;
  // array containing all interaction factors for spin down and spin down (a^+_d a^+_d a_d a_d)
  Complex** InteractionFactorsdowndowndowndown;
  // array containing all interaction factors for spin up and spin down  (a^+_u a^+_d a_u a_u)
  Complex** InteractionFactorsupdownupup;
  // array containing all interaction factors for spin up and spin down  (a^+_u a^+_d a_u a_d)
  Complex** InteractionFactorsupdownupdown;
  // array containing all interaction factors for spin up and spin down  (a^+_u a^+_d a_d a_d)
  Complex** InteractionFactorsupdowndowndown;


 public:

  // default constructor
  //
  ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // bandParameter = band parameter
  // szSymmetryBreaking = amplitude of the Sz symmetry breaking term
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double bandParameter, double szSymmetryBreaking, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian();
  

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
  virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
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

  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the part of the interaction coefficient coming from the two band truncated basis
  //
  // basisMatrices = array of basis matrices
  // momentumIndex1 = momentum index for the first particle
  // momentumIndex2 = momentum index for the second particle
  // bandIndex1 = band index of the first particle
  // bandIndex2 = band index of the second particle
  // return value = coefficient
  Complex ComputeBasisContribution(ComplexMatrix* basisMatrices, int momentumIndex1, int momentumIndex2, int bandIndex1, int bandIndex2);

};

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  int* TmpIndices2;
  Complex* TmpInteractionFactor;
  int Index;
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
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupup[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
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
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndowndowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
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
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupdown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
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

inline void ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
												      ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  int* TmpIndices;
  int* TmpIndices2;
  Complex* TmpInteractionFactor;
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
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupup[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
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
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndowndowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
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
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupdown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
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

inline void ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
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
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupup[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
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
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndowndowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupdown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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

inline void ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
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
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupup[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactorsdowndowndowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactorsdowndownupdown[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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

inline void ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
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
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsdowndownupup[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
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
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsdowndowndowndown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
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
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsdowndownupdown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
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
	  int Lim2 = 0;
	  if (j < this->NbrInterSectorSums)
	  {
	    Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
	    TmpIndices2 = this->InterSectorIndicesPerSum[j];
	  }
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupupupup[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsdowndownupup[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
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
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsdowndowndowndown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		  if (j < this->NbrInterSectorSums)
		  {
		    TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim2) >> 2]);
		  }
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
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
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactorsdowndownupdown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		    }
		  if (j < this->NbrInterSectorSums)
		  {
		    TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
		  }
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
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

inline void ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  int* TmpIndices2;
  Complex* TmpInteractionFactor;
  int Index;

  if (this->HermitianSymmetryFlag == false)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
	      TmpIndices = this->IntraSectorIndicesPerSum[j];
	      int Lim2 = 0;
	      if (j < this->NbrInterSectorSums)
	      {
		Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
		TmpIndices2 = this->InterSectorIndicesPerSum[j];
	      }
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient3 = particles->AuAu(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsupupupup[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      if (j < this->NbrInterSectorSums)
		      {
			TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim2) >> 2]);
		      }
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      if (j < this->NbrInterSectorSums)
		      {
			TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim2) >> 2]);
		      }
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      if (j < this->NbrInterSectorSums)
		      {
		      TmpInteractionFactor = &(this->InteractionFactorsupdownupdown[j][(i1 * Lim2) >> 2]);
		      }
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	      int Lim2 = 0;
	      if (j < this->NbrInterSectorSums)
	      {
		Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
		TmpIndices2 = this->InterSectorIndicesPerSum[j];
	      }
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient3 = particles->AuAu(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsupupupup[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      
		      if (j < this->NbrInterSectorSums)
		      {
			TmpInteractionFactor = &(this->InteractionFactorsupdownupup[j][(i1 * Lim2) >> 2]);
			for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      }
		    }
		  Coefficient3 = particles->AdAd(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      if (j < this->NbrInterSectorSums)
		      {
			TmpInteractionFactor = &(this->InteractionFactorsupdowndowndown[j][(i1 * Lim2) >> 2]);
			for (int i2 = 0; i2 < Lim2; i2 += 2)
			  {
// 			    cout << i << " " ;
// 			    particles->PrintState(cout, i);
// 			    cout << " " << TmpIndices[i1] << " " << TmpIndices[i1 + 1] << " "  << TmpIndices2[i2] << " " << TmpIndices2[i2 + 1] << endl;
			    Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	      for (int i1 = 0; i1 < Lim2; i1 += 2)
		{
		  Coefficient3 = particles->AuAd(i, TmpIndices2[i1], TmpIndices2[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsupupupdown[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
			  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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

#endif
