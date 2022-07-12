////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                     spin and generic 3-body interaction                    //
//                                                                            //
//                        last modification : 26/08/2008                      //
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


#ifndef ABSTRACTQHEONSPHEREWITHSPINNBODYINTERACTIONHAMILTONIANWITHPAIRING_H
#define ABSTRACTQHEONSPHEREWITHSPINNBODYINTERACTIONHAMILTONIANWITHPAIRING_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian.h"

#include <iostream>


using std::ostream;


class ClebschGordanCoefficients;


class AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing : public AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian
{

 protected:

  // this class implements additional fields for pair tunnelling terms a^\dagger_\up a^\dagger_\up a_\down a_\down
  // use conventional storage (slower)
  // number of index sum for the creation (or annhilation) operators for the mixed spin sector
  int NbrPairingSectorSums;
  // array containing the number of index group per index sum for the creation (or annhilation) operators for the mixed spin sector
  int* NbrPairingSectorIndicesPerSum;
  // array containing the (m1,m2) indices per index sum for the creation (or annhilation) operators for the mixed spin sector
  int** PairingSectorIndicesPerSum;

  // array containing all interaction factors for mixed spin (up-up-down-down)
  double** InteractionFactorsPairing;

  
 public:

  // default constructor
  //
  virtual ~AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing() = 0;

  // need to overload multiplication routines to add new terms:

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
					  int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual  RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						   int firstComponent, int nbrComponent);


 protected:
 
  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // return value = number of non-zero matrix element
  virtual long PartialFastMultiplicationMemory(int firstComponent, int lastComponent);

  // enable fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  virtual void PartialEnableFastMultiplication(int firstComponent, int lastComponent);


  // core part of the AddMultiply method involving pairing term
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // nbodyIndex = value of n
  virtual void EvaluatePairingAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector& vSource, RealVector& vDestination);

  // core part of the AddMultiply method involving pairing term for a set ofe vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // nbodyIndex = value of n
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void EvaluatePairingAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector* vSources, 
					   RealVector* vDestinations, int nbrVectors, double* tmpCoefficients);

  // core part of the FastMultiplication method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluatePairingFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
							  int* indexArray, double* coefficientArray, long& position);

    // core part of the PartialFastMultiplicationMemory method involving two-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluatePairingFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory);


};



// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing::EvaluatePairingAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector& vSource, RealVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int* TmpIndices;
  double* TmpInteractionFactor;
  int Index;
  if ( this->NbrPairingSectorSums != 0)
    {
      
      for (int j = 0; j < this->NbrPairingSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrPairingSectorIndicesPerSum[j];
	  TmpIndices = this->PairingSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsPairing[j][(i1 * Lim) >> 2]);
		  Coefficient3 *= vSource[index];
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsPairing[j][(i1 * Lim) >> 2]);
		  Coefficient3 *= vSource[index];
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		      ++TmpInteractionFactor;
		    }
		}

	    }
	}

    }
}


// core part of the AddMultiply method involving pairing term for a set ofe vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// nbodyIndex = value of n
// tmpCoefficients = a temporary array whose size is nbrVectors
inline void AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing::EvaluatePairingAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector* vSources, 
					   RealVector* vDestinations, int nbrVectors, double* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  int* TmpIndices;
  double* TmpInteractionFactor;
  if ( this->NbrPairingSectorSums != 0)
    {
      for (int j = 0; j < this->NbrPairingSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrPairingSectorIndicesPerSum[j];
	  TmpIndices = this->PairingSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsPairing[j][(i1 * Lim) >> 2]);
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
		}
	      Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsPairing[j][(i1 * Lim) >> 2]);
		  for (int p = 0; p < nbrVectors; ++p)
		    tmpCoefficients[p] = Coefficient3 * vSources[p][index];
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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


// core part of the FastMultiplication method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray
inline void AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing::EvaluatePairingFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
							  int* indexArray, double* coefficientArray, long& position)
{
    int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  double* TmpInteractionFactor;
  int Dim = particles->GetHilbertSpaceDimension();
  if (this->NbrPairingSectorSums != 0)
    {
     for (int j = 0; j < this->NbrPairingSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrPairingSectorIndicesPerSum[j];
	  TmpIndices = this->PairingSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient2 = particles->AdAd(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsPairing[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}

	      Coefficient2 = particles->AuAu(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsPairing[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
			  ++position;
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
inline void AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing::EvaluatePairingFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  int Dim = particles->GetHilbertSpaceDimension();

  if (this->NbrPairingSectorSums != 0)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrPairingSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrPairingSectorIndicesPerSum[j];
	      TmpIndices = this->PairingSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient2 = particles->AdAd(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			}
		    }

		  Coefficient2 = particles->AuAu(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
}


#endif
