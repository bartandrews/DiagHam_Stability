////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                 class of quantum spin Hall restricted to four bands        //
//                  with a fully SU(4) symmetry breaking interaction          //
//                                                                            //
//                        last modification : 01/08/2012                      //
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


#ifndef PARTICLEONLATTICEQUANTUMSPINHALLFULLFOURBANDHAMILTONIAN_H
#define PARTICLEONLATTICEQUANTUMSPINHALLFULLFOURBANDHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFourBandHamiltonian.h"
#include "Hamiltonian/AbstractQHEHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian : public ParticleOnLatticeQuantumSpinHallFourBandHamiltonian
{

 protected:
  
  // interaction factors : the table first entry , 
  // first entry is the sigma index for the first creation operator
  // second entry is the sigma index for the second creation operator
  // third entry is the sigma index for the first annihilation operator
  // fourth entry is the sigma index for the second annihilation operator
  // fifth entry is a the sum of annihilation/creation indices
  // sixth entry is the linearized index (annihilation_index * nbr_creation_index + creation_index)
  Complex****** InteractionFactorsSigma;

 public:

  // default constructor
  //
  ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian();

  // destructor
  //
  ~ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian();
  

 protected:
 
  // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the two-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int index, ComplexVector* vSources, 
						     ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the two-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int index, ComplexVector* vSources, 
							      ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // core part of the FastMultiplication method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSU4Spin* particles, int index, 
							    int* indexArray, Complex* coefficientArray, long& position);

  // core part of the PartialFastMultiplicationMemory method involving two-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSU4Spin* particles, int firstComponent, int lastComponent, long& memory);


};

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
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
	  for (int sigma1 = 0; sigma1 < 4; ++sigma1)
	    {
	      Coefficient3 = particles->AsigmaAsigma(index, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
	      if (Coefficient3 != 0.0)
		{
		  Coefficient4 = vSource[index];
		  Coefficient4 *= Coefficient3;
		  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
			  if (Index < Dim)
			    {
			      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
		    {
		      for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
			{			  
			  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j][(i1 * Lim2) >> 2]);
			  for (int i2 = 0; i2 < Lim2; i2 += 2)
			    {
			      Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
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
	}
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  for (int sigma1 = 0; sigma1 < 4; ++sigma1)
	    {
	      for (int sigma2 = sigma1 + 1; sigma2 < 4; ++sigma2)
		{
		  Coefficient3 = particles->AsigmaAsigma(index, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
		  if (Coefficient3 != 0.0)
		    {
		      Coefficient4 = vSource[index];
		      Coefficient4 *= Coefficient3;
		      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			{
			  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j][(i1 * Lim) >> 2]);
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
			      if (Index < Dim)
				{
				  vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
				}
			      ++TmpInteractionFactor;
			    }
			}
		      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			{
			  for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
				{			  
				  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
				  for (int i2 = 0; i2 < Lim2; i2 += 2)
				    {
				      Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
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

inline void ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int index, ComplexVector* vSources, 
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
	  for (int sigma1 = 0; sigma1 < 4; ++sigma1)
	    {
	      Coefficient3 = particles->AsigmaAsigma(index, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
	      if (Coefficient3 != 0.0)
		{
		  for (int p = 0; p < nbrVectors; ++p)
		    tmpCoefficients[p] = Coefficient3 * vSources[p][index];
		  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
			  if (Index < Dim)
			    for (int p = 0; p < nbrVectors; ++p)
			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
			  ++TmpInteractionFactor;
			}
		    }
		  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
		    {
		      for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
			{			  
			  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j][(i1 * Lim) >> 2]);
			  for (int i2 = 0; i2 < Lim2; i2 += 2)
			    {
			      Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
			      if (Index < Dim)
				for (int p = 0; p < nbrVectors; ++p)
				  vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
			      ++TmpInteractionFactor;
			    }
			}
		    }
		}
	    }
	}
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  for (int sigma1 = 0; sigma1 < 4; ++sigma1)
	    {
	      for (int sigma2 = sigma1 + 1; sigma2 < 4; ++sigma2)
		{
		  Coefficient3 = particles->AsigmaAsigma(index, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
		  if (Coefficient3 != 0.0)
		    {
		      for (int p = 0; p < nbrVectors; ++p)
			tmpCoefficients[p] = Coefficient3 * vSources[p][index];
		      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			{
			  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
			      if (Index < Dim)
				for (int p = 0; p < nbrVectors; ++p)
				  vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
			      ++TmpInteractionFactor;
			    }
			}
		      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			{
			  for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
			    {			  
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
			      for (int i2 = 0; i2 < Lim2; i2 += 2)
				{
				  Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
				  if (Index < Dim)
				    for (int p = 0; p < nbrVectors; ++p)
				      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
				  ++TmpInteractionFactor;
				}
			    }
			}
		    }
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

inline void ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
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
	  for (int sigma1 = 0; sigma1 < 4; ++sigma1)
	    {
	      Coefficient3 = particles->AsigmaAsigma(index, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
	      if (Coefficient3 != 0.0)
		{
		  Coefficient4 = vSource[index];
		  Coefficient4 *= Coefficient3;
		  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
			  if (Index <= index)
			    {
			      if (Index < index)
				TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
			      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
		    {
		      for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
			{			  
			  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j][(i1 * Lim2) >> 2]);
			  for (int i2 = 0; i2 < Lim2; i2 += 2)
			    {
			      Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
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
	    }
	}
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  for (int sigma1 = 0; sigma1 < 4; ++sigma1)
	    {
	      for (int sigma2 = sigma1 + 1; sigma2 < 4; ++sigma2)
		{
		  Coefficient3 = particles->AsigmaAsigma(index, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
		  if (Coefficient3 != 0.0)
		    {
		      Coefficient4 = vSource[index];
		      Coefficient4 *= Coefficient3;
		      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			{
			  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j][(i1 * Lim) >> 2]);
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
			      if (Index <= index)
				{
				  if (Index < index)
				    TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
				  vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
				}
			      ++TmpInteractionFactor;
			    }
			}
		      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			{
			  for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
			    {			  
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
			      for (int i2 = 0; i2 < Lim2; i2 += 2)
				{
				  Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
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

inline void ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int index, ComplexVector* vSources, 
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
	  for (int sigma1 = 0; sigma1 < 4; ++sigma1)
	    {
	      Coefficient3 = particles->AsigmaAsigma(index, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
	      if (Coefficient3 != 0.0)
		{
		  for (int p = 0; p < nbrVectors; ++p)
		    tmpCoefficients[p] = Coefficient3 * vSources[p][index];
		  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
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
	      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
		{
		  for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
		    {			  
		      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim2; i2 += 2)
			{
			  Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
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
	}
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  for (int sigma1 = 0; sigma1 < 4; ++sigma1)
	    {
	      for (int sigma2 = sigma1 + 1; sigma2 < 4; ++sigma2)
		{
		  Coefficient3 = particles->AsigmaAsigma(index, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
		  if (Coefficient3 != 0.0)
		    {
		      for (int p = 0; p < nbrVectors; ++p)
			tmpCoefficients[p] = Coefficient3 * vSources[p][index];
		      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			{
			  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
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
		      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			{
			  for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
			    {			  
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
			      for (int i2 = 0; i2 < Lim2; i2 += 2)
				{
				  Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma4, Coefficient);
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

inline void ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSU4Spin* particles, int index, 
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
  int AbsoluteIndex = index + this->PrecalculationShift;

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
	      for (int sigma1 = 0; sigma1 < 4; ++sigma1)
		{
		  Coefficient3 = particles->AsigmaAsigma(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
		  if (Coefficient3 != 0.0)
		    {
		      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			{
			  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j][(i1 * Lim) >> 2]);
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
			      if (Index < Dim)
				{
				  indexArray[position] = Index;
				  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
				  ++position;
				}
			      ++TmpInteractionFactor;
			    }
			}
		      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			{
			  for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
			    {			  
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j][(i1 * Lim2) >> 2]);
			      for (int i2 = 0; i2 < Lim2; i2 += 2)
				{
				  Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
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
	    }
	  for (int i1 = 0; i1 < Lim2; i1 += 2)
	    {
	      for (int sigma1 = 0; sigma1 < 4; ++sigma1)
		{
		  for (int sigma2 = sigma1 + 1; sigma2 < 4; ++sigma2)
		    {
		      Coefficient3 = particles->AsigmaAsigma(AbsoluteIndex, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
		      if (Coefficient3 != 0.0)
			{
			  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			    {
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j][(i1 * Lim) >> 2]);
			      for (int i2 = 0; i2 < Lim; i2 += 2)
				{
				  Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
				  if (Index < Dim)
				    {
				      indexArray[position] = Index;
				      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
				      ++position;
				    }
				  ++TmpInteractionFactor;
				}
			    }
			  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			    {
			      for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
				{			  
				  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
				  for (int i2 = 0; i2 < Lim2; i2 += 2)
				    {
				      Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
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
	  int Count = 0;
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      for (int sigma1 = 0; sigma1 < 4; ++sigma1)
		{
		  Coefficient3 = particles->AsigmaAsigma(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
		  if (Coefficient3 != 0.0)
		    {
		      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			{
			  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j][(i1 * Lim) >> 2]);
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
			      if (Index <= AbsoluteIndex)
				{
				  if (Index == AbsoluteIndex)
				    {
				      indexArray[position] = Index;
				      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * (*TmpInteractionFactor);
				      ++position;
				      ++Count;
				    }
				  else
				    {
				      indexArray[position] = Index;
				      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
				      ++position;
				      ++Count;
				    }
				}
			      ++TmpInteractionFactor;
			    }
			}
		      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			{
			  for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
			    {			  
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j][(i1 * Lim2) >> 2]);
			      for (int i2 = 0; i2 < Lim2; i2 += 2)
				{
				  Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
				  if (Index <= AbsoluteIndex)
				    {
				      if (Index == AbsoluteIndex)
					{
					  indexArray[position] = Index;
					  coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * (*TmpInteractionFactor);
					  ++position;
					  ++Count;
					}
				      else
					{
					  indexArray[position] = Index;
					  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
					  ++position;
					  ++Count;
					}
				    }
				  ++TmpInteractionFactor;
				}
			    }
			}
		    }
		}
	    }
	  for (int i1 = 0; i1 < Lim2; i1 += 2)
	    {
	      for (int sigma1 = 0; sigma1 < 4; ++sigma1)
		{
		  for (int sigma2 = sigma1 + 1; sigma2 < 4; ++sigma2)
		    {
		      Coefficient3 = particles->AsigmaAsigma(AbsoluteIndex, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
		      if (Coefficient3 != 0.0)
			{
			  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			    {
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j][(i1 * Lim) >> 2]);
			      for (int i2 = 0; i2 < Lim; i2 += 2)
				{
				  Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
				  if (Index <= AbsoluteIndex)
				    {
				      if (Index == AbsoluteIndex)
					{
					  indexArray[position] = Index;
					  coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * (*TmpInteractionFactor);
					  ++position;
					  ++Count;
					}
				      else
					{
					  indexArray[position] = Index;
					  coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
					  ++position;
					  ++Count;
					}
				    }
				  ++TmpInteractionFactor;
				}
			    }
			  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			    {
			      for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
				{			  
				  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
				  for (int i2 = 0; i2 < Lim2; i2 += 2)
				    {
				      Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
				      if (Index <= AbsoluteIndex)
					{
					  if (Index == AbsoluteIndex)
					    {
					      indexArray[position] = Index;
					      coefficientArray[position] = 0.5 * Coefficient * Coefficient3 * (*TmpInteractionFactor);
					      ++position;
					      ++Count;
					    }
					  else
					    {
					      indexArray[position] = Index;
					      coefficientArray[position] = Coefficient * Coefficient3 * (*TmpInteractionFactor);
					      ++position;
					      ++Count;
					    }
					}
				      ++TmpInteractionFactor;
				    }
				}
			    }
			}
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

inline void ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSU4Spin* particles, int firstComponent, int lastComponent, long& memory)
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
		  for (int sigma1 = 0; sigma1 < 4; ++sigma1)
		    {
		      Coefficient3 = particles->AsigmaAsigma(i, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
		      if (Coefficient3 != 0.0)
			{
			  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			    {
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j][(i1 * Lim) >> 2]);
			      for (int i2 = 0; i2 < Lim; i2 += 2)
				{
				  Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
				  if (Index < Dim)
				    {
				      ++memory;
				      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
				    }
				  ++TmpInteractionFactor;
				}
			    }
			  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			    {
			      for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
				{			  
				  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j][(i1 * Lim2) >> 2]);
				  for (int i2 = 0; i2 < Lim2; i2 += 2)
				    {
				      Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
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
	      
	      for (int i1 = 0; i1 < Lim2; i1 += 2)
		{
		  for (int sigma1 = 0; sigma1 < 4; ++sigma1)
		    {
		      for (int sigma2 = sigma1 + 1; sigma2 < 4; ++sigma2)
			{
			  Coefficient3 = particles->AsigmaAsigma(i, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
			  if (Coefficient3 != 0.0)
			    {
			      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
				{
				  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j][(i1 * Lim) >> 2]);
				  for (int i2 = 0; i2 < Lim; i2 += 2)
				    {
				      Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
				      if (Index < Dim)
					{
					  ++memory;
					  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
					}
				      ++TmpInteractionFactor;
				    }
				}
			      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
				{
				  for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
				    {			  
				      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
				      for (int i2 = 0; i2 < Lim2; i2 += 2)
					{
					  Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
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
		  for (int sigma1 = 0; sigma1 < 4; ++sigma1)
		    {

		      Coefficient3 = particles->AsigmaAsigma(i, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
		      if (Coefficient3 != 0.0)
			{
			  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			    {
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j][(i1 * Lim) >> 2]);
			      for (int i2 = 0; i2 < Lim; i2 += 2)
				{
				  Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
				  if (Index <= i)
				    {
				      ++memory;
				      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
				    }
				  ++TmpInteractionFactor;
				}
			    }
			  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
			    {
			      for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
				{			  
				  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j][(i1 * Lim2) >> 2]);
				  for (int i2 = 0; i2 < Lim2; i2 += 2)
				    {
				      Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
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
	      for (int i1 = 0; i1 < Lim2; i1 += 2)
		{
		  for (int sigma1 = 0; sigma1 < 4; ++sigma1)
		    {
		      for (int sigma2 = sigma1 + 1; sigma2 < 4; ++sigma2)
			{
			  Coefficient3 = particles->AsigmaAsigma(i, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
			  if (Coefficient3 != 0.0)
			    {
			      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
				{
				  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j][(i1 * Lim) >> 2]);
				  for (int i2 = 0; i2 < Lim; i2 += 2)
				    {
				      Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
				      if (Index <= i)
					{
					  ++memory;
					  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
					}
				      ++TmpInteractionFactor;
				    }
				}
			      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
				{
				  for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
				    {			  
				      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
				      for (int i2 = 0; i2 < Lim2; i2 += 2)
					{
					  Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
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
	    }
	}
    }
}

#endif
