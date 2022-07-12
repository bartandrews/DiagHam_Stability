////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                class of quantum spin Hall restricted to three bands        //
//                  with a fully SU(3) symmetry breaking interaction          //
//                                                                            //
//                        last modification : 30/06/2012                      //
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


#ifndef PARTICLEONLATTICEQUANTUMSPINHALLFULLTHREEBANDHAMILTONIAN_H
#define PARTICLEONLATTICEQUANTUMSPINHALLFULLTHREEBANDHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSU3Spin.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallThreeBandHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian : public ParticleOnLatticeQuantumSpinHallThreeBandHamiltonian
{

 protected:
  
  // array containing all interaction factors for 11 and 12 (a^+_1 a^+_1 a_1 a_2)
  Complex** InteractionFactors1112;
  // array containing all interaction factors for 11 and 13 (a^+_1 a^+_1 a_1 a_3)
  Complex** InteractionFactors1113;
  // array containing all interaction factors for 11 and 23 (a^+_1 a^+_1 a_2 a_2)
  Complex** InteractionFactors1122;
  // array containing all interaction factors for 11 and 23 (a^+_1 a^+_1 a_2 a_3)
  Complex** InteractionFactors1123;
  // array containing all interaction factors for 11 and 23 (a^+_1 a^+_1 a_3 a_3)
  Complex** InteractionFactors1133;

  // array containing all interaction factors for 11 and 12 (a^+_1 a^+_2 a_1 a_1)
  Complex** InteractionFactors1211;
  // array containing all interaction factors for 11 and 12 (a^+_1 a^+_2 a_1 a_3)
  Complex** InteractionFactors1213;
  // array containing all interaction factors for 11 and 13 (a^+_1 a^+_2 a_2 a_2)
  Complex** InteractionFactors1222;
  // array containing all interaction factors for 11 and 13 (a^+_1 a^+_2 a_2 a_3)
  Complex** InteractionFactors1223;
  // array containing all interaction factors for 11 and 23 (a^+_1 a^+_2 a_3 a_3)
  Complex** InteractionFactors1233;

  // array containing all interaction factors for 11 and 12 (a^+_1 a^+_3 a_1 a_1)
  Complex** InteractionFactors1311;
  // array containing all interaction factors for 11 and 12 (a^+_1 a^+_3 a_1 a_2)
  Complex** InteractionFactors1312;
  // array containing all interaction factors for 11 and 13 (a^+_1 a^+_3 a_2 a_2)
  Complex** InteractionFactors1322;
  // array containing all interaction factors for 11 and 13 (a^+_1 a^+_3 a_2 a_3)
  Complex** InteractionFactors1323;
  // array containing all interaction factors for 11 and 23 (a^+_1 a^+_3 a_3 a_3)
  Complex** InteractionFactors1333;

  // array containing all interaction factors for 11 and 12 (a^+_2 a^+_2 a_1 a_1)
  Complex** InteractionFactors2211;
  // array containing all interaction factors for 11 and 12 (a^+_2 a^+_2 a_1 a_2)
  Complex** InteractionFactors2212;
  // array containing all interaction factors for 11 and 13 (a^+_2 a^+_2 a_1 a_3)
  Complex** InteractionFactors2213;
  // array containing all interaction factors for 11 and 23 (a^+_2 a^+_2 a_2 a_3)
  Complex** InteractionFactors2223;
  // array containing all interaction factors for 11 and 23 (a^+_2 a^+_2 a_3 a_3)
  Complex** InteractionFactors2233;

  // array containing all interaction factors for 11 and 12 (a^+_2 a^+_3 a_1 a_1)
  Complex** InteractionFactors2311;
  // array containing all interaction factors for 11 and 12 (a^+_2 a^+_3 a_1 a_2)
  Complex** InteractionFactors2312;
  // array containing all interaction factors for 11 and 12 (a^+_2 a^+_3 a_1 a_3)
  Complex** InteractionFactors2313;
  // array containing all interaction factors for 11 and 13 (a^+_2 a^+_3 a_2 a_2)
  Complex** InteractionFactors2322;
  // array containing all interaction factors for 11 and 23 (a^+_2 a^+_3 a_3 a_3)
  Complex** InteractionFactors2333;

  // array containing all interaction factors for 11 and 12 (a^+_3 a^+_3 a_1 a_1)
  Complex** InteractionFactors3311;
  // array containing all interaction factors for 11 and 12 (a^+_3 a^+_3 a_1 a_2)
  Complex** InteractionFactors3312;
  // array containing all interaction factors for 11 and 13 (a^+_3 a^+_3 a_1 a_3)
  Complex** InteractionFactors3313;
  // array containing all interaction factors for 11 and 23 (a^+_3 a^+_3 a_2 a_2)
  Complex** InteractionFactors3322;
  // array containing all interaction factors for 11 and 23 (a^+_3 a^+_3 a_2 a_3)
  Complex** InteractionFactors3323;


 public:

  // default constructor
  //
  ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian();

  // destructor
  //
  ~ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian();
  

 protected:
 
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
  virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int index, ComplexVector* vSources, 
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

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  //virtual ComplexVector* LowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  //virtual ComplexVector* HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent);

};

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
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
	  Coefficient3 = particles->A1A1(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1111[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors2211[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors1211[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A2A2(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1122[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors2222[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors1222[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	  Coefficient3 = particles->A1A2(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1112[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors2212[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors1212[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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

inline void ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int index, ComplexVector* vSources, 
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
	  Coefficient3 = particles->A1A1(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1111[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors2211[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors1211[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A2A2(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1122[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors2222[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors1222[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	    }
	}
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  Coefficient3 = particles->A1A2(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1112[j][(i1 * Lim2) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors2212[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors1212[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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

inline void ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
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
	  Coefficient3 = particles->A1A1(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1111[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors2211[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors1211[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->A2A2(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1122[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors2222[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors1222[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	  Coefficient3 = particles->A1A2(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1112[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors2212[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpInteractionFactor = &(this->InteractionFactors1212[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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

inline void ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU3Spin* particles, int index, ComplexVector* vSources, 
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
	  Coefficient3 = particles->A1A1(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1111[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors2211[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors3311[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors1211[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors1311[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors2311[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	  Coefficient3 = particles->A2A2(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1122[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors2222[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors3322[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors1222[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors1322[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors2322[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	  Coefficient3 = particles->A3A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1133[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors2233[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors3333[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors1233[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors1333[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors2333[j][(i1 * Lim) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	  Coefficient3 = particles->A1A2(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1112[j][(i1 * Lim2) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors2212[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors3312[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors1212[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors1312[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors2312[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	  Coefficient3 = particles->A1A3(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1113[j][(i1 * Lim2) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors2213[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors3313[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors1213[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors1313[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors2313[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	  Coefficient3 = particles->A2A3(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactors1123[j][(i1 * Lim2) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors2223[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors3323[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors1223[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors1323[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
	      TmpInteractionFactor = &(this->InteractionFactors2323[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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

inline void ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSU3Spin* particles, int index, 
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
	      Coefficient3 = particles->A1A1(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors1111[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors2211[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors3311[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors1211[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors1311[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors2311[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient3 = particles->A2A2(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors1122[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors2222[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors3322[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors1222[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors1322[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors2322[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient3 = particles->A3A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors1133[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors2233[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors3333[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors1233[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors1333[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors2333[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	      Coefficient3 = particles->A1A2(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors1112[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors2212[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors3312[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors1212[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors1312[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors2312[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient3 = particles->A1A3(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors1113[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors2213[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors3313[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors1213[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors1313[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors2313[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient3 = particles->A2A3(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors1123[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors2223[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors3323[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors1223[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors1323[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		  TmpInteractionFactor = &(this->InteractionFactors2323[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	  int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
	  TmpIndices2 = this->InterSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->A1A1(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors1111[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors2211[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors3311[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors1211[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors2311[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors1211[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	      Coefficient3 = particles->A2A2(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors1122[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors2222[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors3322[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors1222[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors1322[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors2322[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	      Coefficient3 = particles->A3A3(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors1133[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors2233[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors3333[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors1233[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors1333[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors2333[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	      Coefficient3 = particles->A1A2(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors1112[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors2212[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors3312[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors1212[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors1312[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors2312[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	      Coefficient3 = particles->A1A3(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors1113[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors2213[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors3313[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors1213[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors1313[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors2313[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
	      Coefficient3 = particles->A2A3(index, TmpIndices2[i1], TmpIndices2[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors1123[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors2223[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors3323[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors1223[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors1323[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  TmpInteractionFactor = &(this->InteractionFactors2323[j][(i1 * Lim2) >> 2]);
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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

// CONTINUE HERE


// core part of the PartialFastMultiplicationMemory method involving two-body term
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSU3Spin* particles, int firstComponent, int lastComponent, long& memory)
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
	      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
	      TmpIndices2 = this->InterSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient3 = particles->A1A1(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors1111[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2211[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors3311[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1211[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1311[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2311[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  Coefficient3 = particles->A2A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors1122[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2222[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors3322[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1222[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1322[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2322[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  Coefficient3 = particles->A3A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors1133[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2233[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors3333[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1233[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1333[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2333[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  Coefficient3 = particles->A1A2(i, TmpIndices2[i1], TmpIndices2[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors1112[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2212[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors3312[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1212[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1312[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2312[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  Coefficient3 = particles->A1A3(i, TmpIndices2[i1], TmpIndices2[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors1113[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2213[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors3313[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1213[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1313[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2313[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  Coefficient3 = particles->A2A3(i, TmpIndices2[i1], TmpIndices2[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors1123[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2223[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors3323[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1223[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1323[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2323[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  Coefficient3 = particles->A1A1(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors1111[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2211[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors3311[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1211[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1311[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2311[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  Coefficient3 = particles->A2A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors1122[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2222[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors3322[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1222[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1322[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2322[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  Coefficient3 = particles->A3A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors1133[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2233[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors3333[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1233[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1333[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2333[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
		  Coefficient3 = particles->A1A2(i, TmpIndices2[i1], TmpIndices2[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors1112[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2212[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors3312[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1212[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1312[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2312[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  Coefficient3 = particles->A1A3(i, TmpIndices2[i1], TmpIndices2[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors1113[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2213[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors3313[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1213[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1313[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2313[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  Coefficient3 = particles->A2A3(i, TmpIndices2[i1], TmpIndices2[i1 + 1]);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors1123[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2223[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors3323[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1223[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad2(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors1323[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad1Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			  ++TmpInteractionFactor;
			}
		      TmpInteractionFactor = &(this->InteractionFactors2323[j][(i1 * Lim2) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->Ad2Ad3(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
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
