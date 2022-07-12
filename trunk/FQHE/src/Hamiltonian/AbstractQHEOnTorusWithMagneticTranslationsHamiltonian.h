////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                to particles on a torus with magnetic translations          //
//                                                                            //
//                        last modification : 18/11/2003                      //
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


#ifndef ABSTRACTQHEONTORUSWITHMAGNETICTRANSLATIONSHAMILTONIAN_H
#define ABSTRACTQHEONTORUSWITHMAGNETICTRANSLATIONSHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian.h"

#include <iostream>


using std::ostream;


class AbstractArchitecture;


class AbstractQHEOnTorusWithMagneticTranslationsHamiltonian : public ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian
{

  friend class QHEParticlePrecalculationOperation;

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
  virtual ~AbstractQHEOnTorusWithMagneticTranslationsHamiltonian();

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
  inline void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
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

};

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int Index;
  int NbrTranslations;
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
		  if ((*TmpInteractionFactor) != 0.0)
		    {
		      Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * (this->ExponentialFactors[NbrTranslations] * Coefficient4);
			}
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

inline void AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
													      ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  int NbrTranslations;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  Complex Tmp;
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
		  if ((*TmpInteractionFactor) != 0.0)
		    {
		      Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  Tmp = Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations];
			  for (int p = 0; p < nbrVectors; ++p)
			    vDestinations[p][Index] += Tmp * tmpCoefficients[p];
			}
		    }
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

inline void AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int Index;
  int NbrTranslations;
  Complex TmpSum = 0.0;
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
		  if ((*TmpInteractionFactor) != 0.0)
		    {
		      Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  if (Index < index)
			    TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(this->ExponentialFactors[NbrTranslations] *  (*TmpInteractionFactor));
			  vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * (this->ExponentialFactors[NbrTranslations] * Coefficient4);
			}
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

inline void AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
														       ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  double Coefficient;
  double Coefficient3;
  int Index;
  int NbrTranslations;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  Complex* TmpSum = new Complex[nbrVectors];
  Complex Tmp;
  Complex Tmp2;
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
		  if ((*TmpInteractionFactor) != 0.0)
		    {
		      Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  if (Index < index)
			    {
			      Tmp = Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations];
			      Tmp2 = (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations]);
			      for (int p = 0; p < nbrVectors; ++p)
				{
				  vDestinations[p][Index] += Tmp * tmpCoefficients[p];
				  TmpSum[p] += Tmp2 * vSources[p][Index];
				}
			    }
			  else
			    {
			      Tmp = Coefficient * (*TmpInteractionFactor) * this->ExponentialFactors[NbrTranslations];
			      for (int p = 0; p < nbrVectors; ++p)
				{
				  vDestinations[p][Index] += Tmp * tmpCoefficients[p];
				}
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

inline void AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
														  int* indexArray, Complex* coefficientArray, long& position)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int NbrTranslations;
  int Dim = particles->GetHilbertSpaceDimension();

  if (this->HermitianSymmetryFlag == false)
    {
      for (int j = 0; j < this->NbrSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrSectorIndicesPerSum[j];
	  TmpIndices = this->SectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient2 = particles->AA(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      if ((*TmpInteractionFactor) != 0.0)
			{
			  Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index < Dim)
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
  else
    {
      for (int j = 0; j < this->NbrSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrSectorIndicesPerSum[j];
	  TmpIndices = this->SectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      int AbsoluteIndex = index + this->PrecalculationShift;
	      Coefficient2 = particles->AA(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactors[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      if ((*TmpInteractionFactor) != 0.0)
			{
			  Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			  if (Index <= AbsoluteIndex)
			    {
			      if (Index == AbsoluteIndex)
				{
				  indexArray[position] = Index;
				  coefficientArray[position] = Coefficient * Coefficient2 * 0.5 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
				  ++position;
				}
			      else
				{
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			      ++position;
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

// core part of the PartialFastMultiplicationMemory method involving two-body term
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  int NbrTranslations;
  int Dim = particles->GetHilbertSpaceDimension();

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
		  Coefficient2 = particles->AA(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      Complex* TmpInteractionFactor = &(this->InteractionFactors[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  if ((*TmpInteractionFactor) != 0.0)
			    {
			      Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			      if (Index < Dim)
				{
				  ++memory;
				  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
				}
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
		  Coefficient2 = particles->AA(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      Complex* TmpInteractionFactor = &(this->InteractionFactors[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  if ((*TmpInteractionFactor) != 0.0)
			    {
			      Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
			      if (Index <= i)
				{
				  ++memory;
				  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
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

#endif
