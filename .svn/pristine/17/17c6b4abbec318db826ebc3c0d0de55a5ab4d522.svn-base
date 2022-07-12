////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2013 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of abstract fractional quantum Hall hamiltonian          //
//              associated to particles with SU(2) spin on a sphere           //
//                including all possible two-body coupling terms              //
//                                                                            //
//                        last modification : --/--/2013                      //
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


#ifndef ABSTRACTQHEONSPHEREWITHSPINFULLHAMILTONIAN_H
#define ABSTRACTQHEONSPHEREWITHSPINFULLHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/AbstractQHEOnSphereWithSpinHamiltonian.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractQHEOnSphereWithSpinFullHamiltonian : public AbstractQHEOnSphereWithSpinHamiltonian
{

 protected:

  // number of index sum for the creation (or annhilation) operators in a given spin sector
  int NbrUpUpSectorSums;
  int NbrUpDownSectorSums;
  int NbrDownDownSectorSums;

  // array containing the number of index group per index sum for the creation (or annhilation) operators in a given spin sector
  int* NbrUpUpSectorIndicesPerSum;
  int* NbrUpDownSectorIndicesPerSum;
  int* NbrDownDownSectorIndicesPerSum;
  // array containing the (m1,m2) indices per index sum for the creation (or annhilation) operators in a given spin sector
  int** UpUpSectorIndicesPerSum;
  int** UpDownSectorIndicesPerSum;
  int** DownDownSectorIndicesPerSum;

  // arrays containing all interaction factors, the first index correspond correspond to index sum for the creation (or annhilation) operators
  // the second index is a linearized index (m1,m2) + (n1,n2) * (nbr element in current index sum) (m for creation operators, n for annhilation operators)
  double** InteractionFactorsUpUpUpUp;
  double** InteractionFactorsUpDownUpUp;
  double** InteractionFactorsDownDownUpUp;
  double** InteractionFactorsUpUpUpDown;
  double** InteractionFactorsUpDownUpDown;
 // double** InteractionFactorsDownUpUpDown;
  double** InteractionFactorsDownDownUpDown;
  double** InteractionFactorsUpUpDownDown;
  double** InteractionFactorsUpDownDownDown;
  double** InteractionFactorsDownDownDownDown;


  // number of one body terms in each spin sector
  int NbrOneBodyInteractionFactorsUpUp;
  int NbrOneBodyInteractionFactorsUpDown;
  int NbrOneBodyInteractionFactorsDownUp;
  int NbrOneBodyInteractionFactorsDownDown;

  // arrays containing all one=body interaction factors
  double* OneBodyInteractionFactorsUpUp;
  double* OneBodyInteractionFactorsUpDown;
  double* OneBodyInteractionFactorsDownUp;
  double* OneBodyInteractionFactorsDownDown;

  // m indices appearing in each one body term
  int* OneBodyMValuesUpUp;
  int* OneBodyMValuesUpDown;
  int* OneBodyMValuesDownUp;
  int* OneBodyMValuesDownDown;

 public:

  // destructor
  //
  virtual ~AbstractQHEOnSphereWithSpinFullHamiltonian() = 0;

 protected:
  
  // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector& vSource, RealVector& vDestination);

  // core part of the AddMultiply method involving the two-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector* vSources, 
						     RealVector* vDestinations, int nbrVectors, double* tmpCoefficients);

  // core part of the FastMultiplication method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
							  int* indexArray, double* coefficientArray, long& position);

  // core part of the AddMultiply method involving the one-body interaction, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
						     int step, RealVector& vSource, RealVector& vDestination);

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
						     int step, RealVector* vSources, RealVector* vDestinations, int nbrVectors);

  // core part of the FastMultiplication method involving the one-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  void EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
						    int* indexArray, double* coefficientArray, long& position);

  // core part of the PartialFastMultiplicationMemory method involving two-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory);

  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors() = 0;

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // return value = number of non-zero matrix element
  virtual long PartialFastMultiplicationMemory(int firstComponent, int lastComponent);

};


// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnSphereWithSpinFullHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector& vSource, RealVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int* TmpIndices;
  int* TmpIndices2;
  double* TmpInteractionFactor;
  int Index;
  int Lim2;

  // UpUp sector
  for (int j = 0; j < this->NbrUpUpSectorSums; ++j)
    {
      int Lim = 2 * this->NbrUpUpSectorIndicesPerSum[j];
      TmpIndices = this->UpUpSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      Coefficient3 *= vSource[index];
	      TmpIndices2 = this->UpUpSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpUpSectorIndicesPerSum[j];
	      TmpInteractionFactor = this->InteractionFactorsUpUpUpUp[j] + ((i1 * Lim2) >> 2);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->DownDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrDownDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = this->InteractionFactorsDownDownUpUp[j] + ((i1 * Lim2) >> 2);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->UpDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = this->InteractionFactorsUpDownUpUp[j] + ((i1 * Lim2) >> 2);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		  ++TmpInteractionFactor;
		}
	    }
	}
    }
  // DownDown sector
  for (int j = 0; j < this->NbrDownDownSectorSums; ++j)
    {
      int Lim = 2 * this->NbrDownDownSectorIndicesPerSum[j];
      TmpIndices = this->DownDownSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      Coefficient3 *= vSource[index];
	      TmpIndices2 = this->UpUpSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpUpSectorIndicesPerSum[j];
	      TmpInteractionFactor = this->InteractionFactorsUpUpDownDown[j] + ((i1 * Lim2) >> 2);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->DownDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrDownDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = this->InteractionFactorsDownDownDownDown[j] + ((i1 * Lim2) >> 2);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->UpDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = this->InteractionFactorsUpDownDownDown[j] + ((i1 * Lim2) >> 2);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		  ++TmpInteractionFactor;		  
		}
	    }
	}
    }
  // UpDown sector
 for (int j = 0; j < this->NbrUpDownSectorSums; ++j)
    {
      int Lim = 2 * this->NbrUpDownSectorIndicesPerSum[j];
      TmpIndices = this->UpDownSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AuAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      Coefficient3 *= vSource[index];
	      Lim2 = 2 * this->NbrUpUpSectorIndicesPerSum[j];
	      TmpInteractionFactor = this->InteractionFactorsUpUpUpDown[j] + ((i1 * Lim2) >> 2);
	      TmpIndices2 = this->UpUpSectorIndicesPerSum[j];
// This part is only for tests
//              int Index4 =  (i1 * Lim2  >> 2);
//                  for (int i2 = 0; i2 < Lim2; i2 += 2)
//                     {
//                          Index = particles->AduAdu(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
//                          if (Index < Dim)
//              {                   vDestination[Index] += Coefficient * this->InteractionFactorsUpUpUpDown[j][Index4] *  Coefficient3;
//              cout << "Index4 = " << Index4 << endl;}
//                         ++Index4;
//                     }
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);

		  if (Index < Dim)
                    vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		  ++TmpInteractionFactor;
		}

	      Lim2 = 2 * this->NbrDownDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = this->InteractionFactorsDownDownUpDown[j] + ((i1 * Lim2) >> 2);
	      TmpIndices2 = this->DownDownSectorIndicesPerSum[j];
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->UpDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = this->InteractionFactorsUpDownUpDown[j] + ((i1 * Lim2) >> 2);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
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

inline void AbstractQHEOnSphereWithSpinFullHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector* vSources, RealVector* vDestinations, int nbrVectors, double* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  
  int* TmpIndices;
  int* TmpIndices2;
  double* TmpInteractionFactor;

  int Lim2;

// UpUp sector
  for (int j = 0; j < this->NbrUpUpSectorSums; ++j)
    {
      int Lim = 2 * this->NbrUpUpSectorIndicesPerSum[j];
      TmpIndices = this->UpUpSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      TmpIndices2 = this->UpUpSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpUpSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsUpUpUpUp[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->DownDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrDownDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsDownDownUpUp[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->UpDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsUpDownUpUp[j][(i1 * Lim2) >> 2]);
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
  // DownDown sector
  for (int j = 0; j < this->NbrDownDownSectorSums; ++j)
    {
      int Lim = 2 * this->NbrDownDownSectorIndicesPerSum[j];
      TmpIndices = this->DownDownSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      TmpIndices2 = this->UpUpSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpUpSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsUpUpDownDown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->DownDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrDownDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsDownDownDownDown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->UpDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsUpDownDownDown[j][(i1 * Lim2) >> 2]);
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
  // UpDown sector
  for (int j = 0; j < this->NbrUpDownSectorSums; ++j)
    {
      int Lim = 2 * this->NbrUpDownSectorIndicesPerSum[j];
      TmpIndices = this->UpDownSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AuAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      TmpIndices2 = this->UpUpSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpUpSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsUpUpUpDown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->DownDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrDownDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsDownDownUpDown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->UpDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsUpDownUpDown[j][(i1 * Lim2) >> 2]);
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

// core part of the FastMultiplication method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void AbstractQHEOnSphereWithSpinFullHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, int* indexArray, double* coefficientArray, long& position)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient3 = 0.0;
  int* TmpIndices;
  int* TmpIndices2;
  double* TmpInteractionFactor;
  int Dim = particles->GetHilbertSpaceDimension();
  int Lim2;

// UpUp sector
  for (int j = 0; j < this->NbrUpUpSectorSums; ++j)
    {
      int Lim = 2 * this->NbrUpUpSectorIndicesPerSum[j];
      TmpIndices = this->UpUpSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpIndices2 = this->UpUpSectorIndicesPerSum[j];
              Lim2 = 2 * this->NbrUpUpSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsUpUpUpUp[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
		      ++position;
		    }
		  ++TmpInteractionFactor;		  
		}
	      TmpIndices2 = this->DownDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrDownDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsDownDownUpUp[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
		      ++position;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->UpDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsUpDownUpUp[j][(i1 * Lim2) >> 2]);
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
  // DownDown sector
  for (int j = 0; j < this->NbrDownDownSectorSums; ++j)
    {
      int Lim = 2 * this->NbrDownDownSectorIndicesPerSum[j];
      TmpIndices = this->DownDownSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpIndices2 = this->UpUpSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpUpSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsUpUpDownDown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
		      ++position;
		    }
		  ++TmpInteractionFactor;		  
		}
	      TmpIndices2 = this->DownDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrDownDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsDownDownDownDown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
		      ++position;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->UpDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsUpDownDownDown[j][(i1 * Lim2) >> 2]);
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
  // UpDown sector
  for (int j = 0; j < this->NbrUpDownSectorSums; ++j)
    {
      int Lim = 2 * this->NbrUpDownSectorIndicesPerSum[j];
      TmpIndices = this->UpDownSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  Coefficient3 = particles->AuAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpIndices2 = this->UpUpSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpUpSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsUpUpUpDown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
		      ++position;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->DownDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrDownDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsDownDownUpDown[j][(i1 * Lim2) >> 2]);
	      for (int i2 = 0; i2 < Lim2; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
		      ++position;
		    }
		  ++TmpInteractionFactor;
		}
	      TmpIndices2 = this->UpDownSectorIndicesPerSum[j];
	      Lim2 = 2 * this->NbrUpDownSectorIndicesPerSum[j];
	      TmpInteractionFactor = &(this->InteractionFactorsUpDownUpDown[j][(i1 * Lim2) >> 2]);
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


// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnSphereWithSpinFullHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, int step, RealVector& vSource, RealVector& vDestination)
{

  if ((this->OneBodyInteractionFactorsDownDown != 0) || (this->OneBodyInteractionFactorsUpUp != 0))
  {
    double TmpDiagonal = 0.0;
    for (int i = firstComponent; i < lastComponent; i += step)
      {
        TmpDiagonal = 0.0;
        if (this->OneBodyInteractionFactorsUpUp != 0)
          for (int j = 0; j < this->NbrOneBodyInteractionFactorsUpUp; ++j)
            TmpDiagonal += this->OneBodyInteractionFactorsUpUp[j] * particles->AduAu(i, this->OneBodyMValuesUpUp[j]);
        if (this->OneBodyInteractionFactorsDownDown != 0) 
          for (int j = 0; j < this->NbrOneBodyInteractionFactorsDownDown; ++j) 
            TmpDiagonal += this->OneBodyInteractionFactorsDownDown[j] * particles->AddAd(i, this->OneBodyMValuesDownDown[j]);
        vDestination[i] += (this->HamiltonianShift + TmpDiagonal)* vSource[i];
      }
    }

//  cout << "NbrOneBodyInteractionFactorsDownDown = " << this->NbrOneBodyInteractionFactorsDownDown << endl;
//  for (int i=0; i < this->NbrOneBodyInteractionFactorsDownDown; i++)
//    cout << this->OneBodyInteractionFactorsDownDown[i] << " ";
//  cout << endl;
//  cout << "OneBodyMValuesDownDown " << endl;
//  for (int i=0; i < this->NbrOneBodyInteractionFactorsDownDown; i++)
//    cout << this->OneBodyMValuesDownDown[i] << " ";
//  cout << endl;

  if ((this->OneBodyInteractionFactorsUpDown != 0) || (this->OneBodyInteractionFactorsDownUp != 0))
    {
      double Coefficient;
      double Source;
      int Dim = particles->GetHilbertSpaceDimension();
      int Index;
      for (int i = firstComponent; i < lastComponent; i += step)
	{
	  Source = vSource[i];
	  for (int j = 0; j < this->NbrOneBodyInteractionFactorsDownUp; ++j)
	    {
	      Index = particles->AddAu(i + this->PrecalculationShift, this->OneBodyMValuesDownUp[j], Coefficient);
	      if (Index < Dim)
		  vDestination[Index] += this->OneBodyInteractionFactorsDownUp[j] * OneBodyInteractionFactorsUpDown[j] * Source;
	    }
	  for (int j = 0; j < this->NbrOneBodyInteractionFactorsUpDown; ++j)
	    {
	      Index = particles->AduAd(i + this->PrecalculationShift, this->OneBodyMValuesUpDown[j], Coefficient);
	      if (Index < Dim)
		{
		  vDestination[Index] += this->OneBodyInteractionFactorsUpDown[j] * OneBodyInteractionFactorsUpDown[j] * Source;
		}
	    }
	}
    }

//  for (int i = firstComponent; i < lastComponent; i += step)
//    vDestination[i] += this->HamiltonianShift * vSource[i];

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

inline void AbstractQHEOnSphereWithSpinFullHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, int step, RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  // UpUp and DownDown One-body operators
  if ((this->OneBodyInteractionFactorsDownDown != 0) || (this->OneBodyInteractionFactorsUpUp != 0))
  {
    double TmpDiagonal = 0.0;
    for (int p = 0; p < nbrVectors; ++p)
      {
        RealVector& TmpSourceVector = vSources[p];
        RealVector& TmpDestinationVector = vDestinations[p];
        for (int i = firstComponent; i < lastComponent; i += step)
        { 
          TmpDiagonal = 0.0;
          if (this->OneBodyInteractionFactorsUpUp != 0) 
            for (int j = 0; j < this->NbrOneBodyInteractionFactorsUpUp; ++j) 
              TmpDiagonal += this->OneBodyInteractionFactorsUpUp[j] * particles->AduAu(i, this->OneBodyMValuesUpUp[j]);
          if (this->OneBodyInteractionFactorsDownDown != 0) 
            for (int j = 0; j < this->NbrOneBodyInteractionFactorsDownDown; ++j) 
              TmpDiagonal += this->OneBodyInteractionFactorsDownDown[j] * particles->AddAd(i, this->OneBodyMValuesDownDown[j]);
          TmpDestinationVector[i] += (this->HamiltonianShift + TmpDiagonal ) * TmpSourceVector[i];
        }
      }
    }

  // Hamiltonian Shift
//  for (int p = 0; p < nbrVectors; ++p)
//    {
//      RealVector& TmpSourceVector = vSources[p];
//      RealVector& TmpDestinationVector = vDestinations[p];
//      for (int i = firstComponent; i < lastComponent; i += step)
//        TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
//    }

  if (this->OneBodyInteractionFactorsDownUp != 0)
    {
      double Coefficient;
      int Dim = particles->GetHilbertSpaceDimension();
      int Index;
      for (int i = firstComponent; i < lastComponent; i += step)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactorsDownUp; ++j)
	    {
	      Index = particles->AddAu(i + this->PrecalculationShift, this->OneBodyMValuesDownUp[j], Coefficient);
	      if (Index < Dim)
		{
		  for (int p = 0; p < nbrVectors; ++p)
		    vDestinations[p][Index] += Coefficient * this->OneBodyInteractionFactorsDownUp[j] * vSources[p][i];
		}
	    }
	}
    }
  if (this->OneBodyInteractionFactorsUpDown != 0)
    {
      double Coefficient;
      int Dim = particles->GetHilbertSpaceDimension();
      int Index;
      for (int i = firstComponent; i < lastComponent; i += step)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactorsUpDown; ++j)
	    {	  
	      Index = particles->AduAd(i + this->PrecalculationShift, this->OneBodyMValuesUpDown[j], Coefficient);
	      if (Index < Dim)
		{
		  for (int p = 0; p < nbrVectors; ++p)
		    vDestinations[p][Index] += Coefficient * this->OneBodyInteractionFactorsUpDown[j] * vSources[p][i];
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

inline void AbstractQHEOnSphereWithSpinFullHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
												     int* indexArray, double* coefficientArray, long& position)
{
  if ((this->OneBodyInteractionFactorsDownDown != 0) || (this->OneBodyInteractionFactorsUpUp != 0))
    {
      double TmpDiagonal = 0.0;
      if (this->OneBodyInteractionFactorsUpUp != 0)
        for (int j = 0; j < this->NbrOneBodyInteractionFactorsUpUp; ++j)
          TmpDiagonal += this->OneBodyInteractionFactorsUpUp[j] * particles->AduAu(index + this->PrecalculationShift, this->OneBodyMValuesUpUp[j]);
      if (this->OneBodyInteractionFactorsDownDown != 0)
        for (int j = 0; j < this->NbrOneBodyInteractionFactorsDownDown; ++j)
          TmpDiagonal += this->OneBodyInteractionFactorsDownDown[j] * particles->AddAd(index + this->PrecalculationShift, this->OneBodyMValuesDownDown[j]);	  
	  indexArray[position] = index + this->PrecalculationShift;
	  coefficientArray[position] = TmpDiagonal;
	  ++position;	  
    }
  if (this->OneBodyInteractionFactorsUpDown != 0)
    {
      int Dim = particles->GetHilbertSpaceDimension();
      double Coefficient;
      int Index;
      for (int j = 0; j < this->NbrOneBodyInteractionFactorsDownUp; ++j)
	{
	  Index = particles->AddAu(index + this->PrecalculationShift, this->OneBodyMValuesDownUp[j], Coefficient);
	  if (Index < Dim)
	    {
	      indexArray[position] = Index;
	      coefficientArray[position] = Coefficient * this->OneBodyInteractionFactorsDownUp[j];
	      ++position;
	    }
	}
    }
  if (this->OneBodyInteractionFactorsDownUp != 0)
    {
      int Dim = particles->GetHilbertSpaceDimension();
      double Coefficient;
      int Index;
      for (int j = 0; j < this->NbrOneBodyInteractionFactorsUpDown; ++j)
	{
	  Index = particles->AduAd(index + this->PrecalculationShift, this->OneBodyMValuesUpDown[j], Coefficient);
	  if (Index < Dim)
	    {
	      indexArray[position] = Index;
	      coefficientArray[position] = Coefficient * this->OneBodyInteractionFactorsUpDown[j];
	      ++position;
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

inline void AbstractQHEOnSphereWithSpinFullHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient3 = 0.0;
  int* TmpIndices;
  int* TmpIndices2;
  // double* TmpInteractionFactor;
  int Dim = particles->GetHilbertSpaceDimension();
  int Lim2;
  //int SumIndices;
  //int TmpNbrM3Values;
  //int* TmpM3Values;

  for (int i = firstComponent; i < lastComponent; ++i)
    {
// UpUp sector
      for (int j = 0; j < this->NbrUpUpSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrUpUpSectorIndicesPerSum[j];
	  TmpIndices = this->UpUpSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->AuAu(i, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpIndices2 = this->UpUpSectorIndicesPerSum[j];
		  Lim2 = 2 * this->NbrUpUpSectorIndicesPerSum[j];
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		  TmpIndices2 = this->DownDownSectorIndicesPerSum[j];
		  Lim2 = 2 * this->NbrDownDownSectorIndicesPerSum[j];
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		  TmpIndices2 = this->UpDownSectorIndicesPerSum[j];
		  Lim2 = 2 * this->NbrUpDownSectorIndicesPerSum[j];
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		}
	    }
	}
	//cout << " UpUp : NbrInteractionPerComponent[" << i << "] = " << this->NbrInteractionPerComponent[i] << endl;
      // DownDown sector
      for (int j = 0; j < this->NbrDownDownSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrDownDownSectorIndicesPerSum[j];
	  TmpIndices = this->DownDownSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->AdAd(i, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpIndices2 = this->UpUpSectorIndicesPerSum[j];
		  Lim2 = 2 * this->NbrUpUpSectorIndicesPerSum[j];
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		  TmpIndices2 = this->DownDownSectorIndicesPerSum[j];
		  Lim2 = 2 * this->NbrDownDownSectorIndicesPerSum[j];
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		  TmpIndices2 = this->UpDownSectorIndicesPerSum[j];
		  Lim2 = 2 * this->NbrUpDownSectorIndicesPerSum[j];
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		}
	    }
	}
	//cout << " DownDown : NbrInteractionPerComponent[" << i << "] = " << this->NbrInteractionPerComponent[i] << endl;
      // UpDown sector
      for (int j = 0; j < this->NbrUpDownSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrUpDownSectorIndicesPerSum[j];
      //cout << "j = " << j << endl;
      //cout << "Lim = " << Lim << endl;
	  TmpIndices = this->UpDownSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->AuAd(i, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpIndices2 = this->UpUpSectorIndicesPerSum[j];
		  Lim2 = 2 * this->NbrUpUpSectorIndicesPerSum[j];
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      //cout << "Dim = " << Dim << endl;
		      //cout << "Index = " << Index << endl;
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		  TmpIndices2 = this->DownDownSectorIndicesPerSum[j];
		  Lim2 = 2 * this->NbrDownDownSectorIndicesPerSum[j];
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		  TmpIndices2 = this->UpDownSectorIndicesPerSum[j];
		  Lim2 = 2 * this->NbrUpDownSectorIndicesPerSum[j];
		  for (int i2 = 0; i2 < Lim2; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices2[i2], TmpIndices2[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		}
	    }
	}
	//cout << " UpDown : NbrInteractionPerComponent[" << i << "] = " << this->NbrInteractionPerComponent[i] << endl;
    }
  
  if ((this->OneBodyInteractionFactorsDownDown != 0) || (this->OneBodyInteractionFactorsUpUp != 0))
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  ++memory;
	  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
	}
    }
  if (this->OneBodyInteractionFactorsUpDown != 0)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactorsUpDown; ++j)
	    {
	      Index = particles->AduAd(i, this->OneBodyMValuesUpDown[j], Coefficient);
	      if (Index < Dim)
		{
		  ++memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	    }
	}
    }
  if (this->OneBodyInteractionFactorsDownUp != 0)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactorsDownUp; ++j)
	    {
	      Index = particles->AddAu(i, this->OneBodyMValuesDownUp[j], Coefficient);
	      if (Index < Dim)
		{
		  ++memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	    }
	}
    }
}


#endif
