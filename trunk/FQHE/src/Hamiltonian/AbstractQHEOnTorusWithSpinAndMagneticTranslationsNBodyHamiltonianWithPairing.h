////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                       class author: Nicolas Regnault                       //
//                                                                            //
//              class of quantum Hall many-body hamiltonian associated        //
//           to particles on a torus with spin and magnetic translations      //
//                                   and pairing                              //
//                                                                            //
//                        last modification : 08/12/2015                      //
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


#ifndef ABSTRACTQHEONTORUSWITHSPINANDMAGNETICTRANSLATIONSNBODYHAMILTONIANWITHPAIRING_H
#define ABSTRACTQHEONTORUSWITHSPINANDMAGNETICTRANSLATIONSNBODYHAMILTONIANWITHPAIRING_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"
#include "Hamiltonian/AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian.h"


#include <iostream>


using std::ostream;


class AbstractArchitecture;


class AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonianWithPairing : public AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian
{

 protected:
  
  // array containing all interaction factors for the pairing (a^+_d a^+_d a_u a_u or a^+_u a^+_u a_d a_d)
  Complex** InteractionFactorsupupdowndown;

 public:

  // destructor
  //
  virtual ~AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonianWithPairing();
  
 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors() = 0;

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

};

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonianWithPairing::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, 
														     ComplexVector& vSource, ComplexVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  Complex* TmpInteractionFactorPairing;
  int Index;
  int NbrTranslations;
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
	      TmpInteractionFactorPairing = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		    }
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing))) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		  ++TmpInteractionFactorPairing;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
	      TmpInteractionFactorPairing = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing))) * Coefficient4;
		  ++TmpInteractionFactor;
		  ++TmpInteractionFactorPairing;
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
	  Coefficient3 = particles->AuAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonianWithPairing::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
												      ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  int NbrTranslations;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  Complex* TmpInteractionFactorPairing;
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
	      TmpInteractionFactorPairing = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing)) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		  ++TmpInteractionFactorPairing;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
	      TmpInteractionFactorPairing = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing)) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		  ++TmpInteractionFactorPairing;
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
	  Coefficient3 = particles->AuAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonianWithPairing::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  //int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  Complex* TmpInteractionFactorPairing;
  Complex TmpSum = 0.0;
  int Index;
  int NbrTranslations;
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
	      TmpInteractionFactorPairing = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		    }
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing))) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		  ++TmpInteractionFactorPairing;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
	      TmpInteractionFactorPairing = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
		      vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		    }
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
		      vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing))) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		  ++TmpInteractionFactorPairing;
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
	  Coefficient3 = particles->AuAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
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

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonianWithPairing::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
													       ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  //int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  int NbrTranslations;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  Complex* TmpInteractionFactorPairing;
  Complex* TmpSum = new Complex[nbrVectors];
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
	      TmpInteractionFactorPairing = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
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
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing)) * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj(this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing)) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing)) * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		  ++TmpInteractionFactorPairing;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
	      TmpInteractionFactorPairing = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index <= index)
		    {
		      if (Index < index)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing)) * tmpCoefficients[p];
			      TmpSum[p] += (Coefficient * Coefficient3) * Conj(this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing)) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing)) * tmpCoefficients[p];
			    }
			}
		    }
		  ++TmpInteractionFactor;
		  ++TmpInteractionFactorPairing;
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
	  Coefficient3 = particles->AuAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonianWithPairing::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
																       int* indexArray, Complex* coefficientArray, long& position)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  Complex* TmpInteractionFactorPairing;
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
	      Coefficient2 = particles->AuAu(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupup[j][(i1 * Lim) >> 2]);
		  TmpInteractionFactorPairing = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
                          indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			  ++position;
			}
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
                          indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing));
			  ++position;
			}
		      ++TmpInteractionFactor;
		      ++TmpInteractionFactorPairing;
		    }
		}
	      Coefficient2 = particles->AdAd(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
		  TmpInteractionFactorPairing = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
                          indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			  ++position;
			}
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
                          indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing));
			  ++position;
			}
		      ++TmpInteractionFactor;
		      ++TmpInteractionFactorPairing;
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
	      Coefficient2 = particles->AuAd(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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
	      Coefficient2 = particles->AuAu(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupup[j][(i1 * Lim) >> 2]);
		  TmpInteractionFactorPairing = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing));
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing));
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		      ++TmpInteractionFactorPairing;
		    }
		}
	      Coefficient2 = particles->AdAd(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
		  TmpInteractionFactorPairing = &(this->InteractionFactorsupupdowndown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = 0.5 * Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing));
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactorPairing));
			      ++position;
			    }
			}
		      ++TmpInteractionFactor;
		      ++TmpInteractionFactorPairing;
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
	      Coefficient2 = particles->AuAd(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonianWithPairing::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
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
		  Coefficient2 = particles->AuAu(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			}
		    }
		  Coefficient2 = particles->AdAd(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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
		  Coefficient2 = particles->AuAd(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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
		  Coefficient2 = particles->AuAu(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			}
		    }
		  Coefficient2 = particles->AdAd(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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
		  Coefficient2 = particles->AuAd(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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
