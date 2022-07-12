////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                 class of quantum spin Hall restricted to two bands         //
//                  with a fully SU(2) symmetry breaking interaction          //
//          (this code is an alternative, sipler but slower version of        //
//                 ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian  )      //
//                                                                            //
//                        last modification : 21/02/2013                      //
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


#ifndef PARTICLEONLATTICEQUANTUMSPINHALLFULLTWOBANDHAMILTONIAN_H
#define PARTICLEONLATTICEQUANTUMSPINHALLFULLTWOBANDHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian.h"
#include "Hamiltonian/AbstractQHEHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian : public ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian
{

 protected:
  
  // interaction factors
  // first entry is the sigma index for the first creation operator
  // second entry is the sigma index for the second creation operator
  // third entry is the sigma index for the first annihilation operator
  // fourth entry is the sigma index for the second annihilation operator
  // fifth entry is a the sum of annihilation/creation indices
  // sixth entry is the linearized index (annihilation_index * nbr_creation_index + creation_index)
  Complex****** InteractionFactorsSigma;

  // number of internal indices (i.e number of entries for indices 0, 1, 2 and 3 of InteractionFactorsSigma)
  int NbrInternalIndices;

  // one-body interaction factor 
  // first entry is the sigma index for the creation operator
  // second entry is the sigma index for the annihilation operator
  // third entry is the linearized momentum index
  Complex*** OneBodyInteractionFactorsSigma;
  
 public:

  // default constructor
  //
  ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian();

  // destructor
  //
  ~ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian();
  

 protected:
 
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


};

// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
												      int step, ComplexVector& vSource, ComplexVector& vDestination)
{
  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
    {
      if (this->OneBodyInteractionFactorsSigma[sigma1][sigma1] != 0)
	{
	  double TmpDiagonal = 0.0;
	  for (int i = firstComponent; i < lastComponent; i += step)
	    { 
	      TmpDiagonal = 0.0;
	      for (int j = 0; j <= this->LzMax; ++j) 
		{
		  TmpDiagonal += this->OneBodyInteractionFactorsSigma[sigma1][sigma1][j].Re * particles->AdsigmaAsigma(i, j, sigma1);
		}
	      vDestination[i] += TmpDiagonal * vSource[i];
	    }
	}
    }
  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
    {
      for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
	{
	  if (this->OneBodyInteractionFactorsSigma[sigma1][sigma2] != 0)
	    {
	      double Coefficient;
	      Complex Source;
	      int Dim = particles->GetHilbertSpaceDimension();
	      int Index;
	      if (this->HermitianSymmetryFlag == false)
		{
		  for (int i = firstComponent; i < lastComponent; i += step)
		    {
		      Source = vSource[i];
		      for (int j = 0; j <= this->LzMax; ++j)
			{
			  Index = particles->AdsigmaAsigma(i, j, sigma1, j, sigma2, Coefficient);
			  if (Index < Dim)
			    {
			      vDestination[Index] += (Coefficient * this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]) * Source;
			    }
			  Index = particles->AdsigmaAsigma(i, j, sigma2, j, sigma1, Coefficient);
			  if (Index < Dim)
			    {
			      vDestination[Index] += (Coefficient * Conj(this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j])) * Source;
			    }
			}
		    }
		}
	      else
		{
		  for (int i = firstComponent; i < lastComponent; i += step)
		    {
		      Source = vSource[i];
		      for (int j = 0; j <= this->LzMax; ++j)
			{
			  Index = particles->AdsigmaAsigma(i, j, sigma1, j, sigma2, Coefficient);
			  if (Index <= i)
			    {
			      if (Index < i)
				{
				  vDestination[Index] += (Coefficient * this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]) * Source;
				  vDestination[i] += (Coefficient * Conj(this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j])) * vSource[Index];
				}
			      else
				{
				  vDestination[Index] += (Coefficient * this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]) * Source;
				}
			    }
			  Index = particles->AdsigmaAsigma(i, j, sigma2, j, sigma1, Coefficient);
			  if (Index < Dim)
			    {
			      if (Index <= i)
				{
				  if (Index < i)
				    {
				      vDestination[Index] += (Coefficient * Conj(this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j])) * Source;
				      vDestination[i] += (Coefficient * this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]) * vSource[Index];
				    }
				  else
				    {
				      vDestination[Index] += (Coefficient * Conj(this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j])) * Source;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  
  for (int i = firstComponent; i < lastComponent; i += step)
    {
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

inline void ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
												      int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
    {
      if (this->OneBodyInteractionFactorsSigma[sigma1][sigma1] != 0)
	{
	  double TmpDiagonal = 0.0;
	  for (int i = firstComponent; i < lastComponent; i += step)
	    { 
	      TmpDiagonal = 0.0;
	      for (int j = 0; j <= this->LzMax; ++j) 
		{
		  TmpDiagonal += this->OneBodyInteractionFactorsSigma[sigma1][sigma1][j].Re * particles->AdsigmaAsigma(i, j, sigma1);
		}
	      for (int p = 0; p < nbrVectors; ++p)
		{
		  vDestinations[p][i] += TmpDiagonal * vSources[p][i];
		}
	    }
	}
    }
  for (int p = 0; p < nbrVectors; ++p)
    {
      ComplexVector& TmpSourceVector = vSources[p];
      ComplexVector& TmpDestinationVector = vDestinations[p];
      for (int i = firstComponent; i < lastComponent; i += step)
	TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
    }
  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
    {
      for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
	{
	  if (this->OneBodyInteractionFactorsSigma[sigma1][sigma2] != 0)
	    {
	      double Coefficient;
	      int Dim = particles->GetHilbertSpaceDimension();
	      int Index;
	      if (this->HermitianSymmetryFlag == false)
		{
		  for (int i = firstComponent; i < lastComponent; i += step)
		    {
		      for (int j = 0; j <= this->LzMax; ++j)
			{
			  Index = particles->AdsigmaAsigma(i, j, sigma1, j, sigma2, Coefficient);
			  if (Index < Dim)
			    {
			      for (int p = 0; p < nbrVectors; ++p)
				vDestinations[p][Index] += Coefficient * (this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]) * vSources[p][i];
			    }
			  Index = particles->AdsigmaAsigma(i, j, sigma2, j, sigma1, Coefficient);
			  if (Index < Dim)
			    {
			      for (int p = 0; p < nbrVectors; ++p)
				vDestinations[p][Index] += Coefficient * Conj(this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]) * vSources[p][i];
			    }
			}
		    }
		}
	      else
		{
		  for (int i = firstComponent; i < lastComponent; i += step)
		    {
		      for (int j = 0; j <= this->LzMax; ++j)
			{
			  Index = particles->AdsigmaAsigma(i, j, sigma1, j, sigma2, Coefficient);
			  if (Index <= i)
			    {
			      if (Index < i)
				{
				  for (int p = 0; p < nbrVectors; ++p)
				    {
				      vDestinations[p][Index] += Coefficient * (this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]) * vSources[p][i];
				      vDestinations[p][i] += Coefficient * Conj(this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]) * vSources[p][Index];
				    }
				}
			      else
				{
				  for (int p = 0; p < nbrVectors; ++p)
				    {
				      vDestinations[p][Index] += Coefficient * (this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]) * vSources[p][i];
				    }
				}
			    }
			  Index = particles->AdsigmaAsigma(i, j, sigma2, j, sigma1, Coefficient);
			  if (Index <= i)
			    {
			      if (Index < i)
				{
				  for (int p = 0; p < nbrVectors; ++p)
				    {
				      vDestinations[p][Index] += Coefficient * Conj(this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]) * vSources[p][i];
				      vDestinations[p][i] += Coefficient * (this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]) * vSources[p][Index];
				    }
				}
			      else
				{
				  for (int p = 0; p < nbrVectors; ++p)
				    vDestinations[p][Index] += Coefficient * Conj(this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]) * vSources[p][i];
				}
			    }
			}
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

inline void ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
														     int* indexArray, Complex* coefficientArray, long& position)
{
  double TmpDiagonal = 0.0;
  bool TmpFlag = false;
  int AbsoluteIndex = index + this->PrecalculationShift;
  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
    {
      if (this->OneBodyInteractionFactorsSigma[sigma1][sigma1] != 0)
	{
	  TmpFlag = true;
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      TmpDiagonal += this->OneBodyInteractionFactorsSigma[sigma1][sigma1][j].Re * particles->AdsigmaAsigma(AbsoluteIndex, j, sigma1);	      
	    }
	}
    }
  if (TmpFlag == true)
    {
      indexArray[position] = AbsoluteIndex;
      if (this->HermitianSymmetryFlag == true)
	TmpDiagonal *= 0.5;
      coefficientArray[position] = TmpDiagonal;
      ++position;
    }

  if (this->HermitianSymmetryFlag == false)
    {
      for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	{
	  for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
	    {
	      if (this->OneBodyInteractionFactorsSigma[sigma1][sigma2] != 0)
		{
		  int Dim = particles->GetHilbertSpaceDimension();
		  double Coefficient;
		  int Index;
		  for (int j = 0; j <= this->LzMax; ++j)
		    {
		      Index = particles->AdsigmaAsigma(AbsoluteIndex, j, sigma1, j, sigma2, Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j];
		      ++position;
			}
		      Index = particles->AdsigmaAsigma(AbsoluteIndex, j, sigma2, j, sigma1, Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Conj(this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]);
			  ++position;
			}
		    }
		}
	    }
	}
    }
  else
    {
      for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	{
	  for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
	    {
	      if (this->OneBodyInteractionFactorsSigma[sigma1][sigma2] != 0)
		{
		  int Dim = particles->GetHilbertSpaceDimension();
		  double Coefficient;
		  int Index;
		  for (int j = 0; j <= this->LzMax; ++j)
		    {
		      Index = particles->AdsigmaAsigma(AbsoluteIndex, j, sigma1, j, sigma2, Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  indexArray[position] = Index;
			  if (Index == AbsoluteIndex)
			    {
			      coefficientArray[position] = 0.5 * Coefficient * this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j];
			    }
			  else
			    {
			      coefficientArray[position] = Coefficient * this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j];
			    }
			  ++position;
			}
		      Index = particles->AdsigmaAsigma(AbsoluteIndex, j, sigma2, j, sigma1, Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  indexArray[position] = Index;
			  if (Index == AbsoluteIndex)
			    {
			      coefficientArray[position] = 0.5 * Coefficient * Conj(this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]);
			    }
			  else
			    {
			      coefficientArray[position] = Coefficient * Conj(this->OneBodyInteractionFactorsSigma[sigma1][sigma2][j]);
			    }			    
			  ++position;
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

inline void ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;

  int Dim = particles->GetHilbertSpaceDimension();
  
  bool TmpFlag = false;
  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
    {
      if (this->OneBodyInteractionFactorsSigma[sigma1][sigma1] != 0)
	{
	  TmpFlag = true;
	}
   }
  if (TmpFlag == true)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  ++memory;
	  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
	}
    }

  if (this->HermitianSymmetryFlag == false)
    {
      for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	{
	  for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
	    {
	      if (this->OneBodyInteractionFactorsSigma[sigma1][sigma2] != 0)
		{
		  for (int i = firstComponent; i < lastComponent; ++i)
		    {
		      for (int j = 0; j <= this->LzMax; ++j)
			{
			  Index = particles->AdsigmaAsigma(i, j, sigma1, j, sigma2, Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			  Index = particles->AdsigmaAsigma(i, j, sigma2, j, sigma1, Coefficient);
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
      for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	{
	  for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
	    {
	      if (this->OneBodyInteractionFactorsSigma[sigma1][sigma2] != 0)
		{
		  for (int i = firstComponent; i < lastComponent; ++i)
		    {
		      for (int j=0; j <= this->LzMax; ++j)
			{
			  Index = particles->AdsigmaAsigma(i, j, sigma1, j, sigma2, Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			  Index = particles->AdsigmaAsigma(i, j, sigma2, j, sigma1, Coefficient);
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

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
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
	  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	    {
	      Coefficient3 = particles->AsigmaAsigma(index, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
	      if (Coefficient3 != 0.0)
		{
		  Coefficient4 = vSource[index];
		  Coefficient4 *= Coefficient3;
		  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
		    {
		      if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1] != 0)
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
		    }
		  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
		    {
		      for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
			{			  
			  if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1] != 0)
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
	}
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	    {
	      for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
		{
		  Coefficient3 = particles->AsigmaAsigma(index, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
		  if (Coefficient3 != 0.0)
		    {
		      Coefficient4 = vSource[index];
		      Coefficient4 *= Coefficient3;
		      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			{
			  if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2] != 0)
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
			}
		      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			{
			  for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
			    {			  
			      if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] != 0)
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
}

// core part of the AddMultiply method involving the two-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
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
	  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	    {
	      Coefficient3 = particles->AsigmaAsigma(index, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
	      if (Coefficient3 != 0.0)
		{
		  for (int p = 0; p < nbrVectors; ++p)
		    tmpCoefficients[p] = Coefficient3 * vSources[p][index];
		  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
		    {
		      if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1] != 0)
			{
			  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j][(i1 * Lim) >> 2]);
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
			      if (Index < Dim)
				{
				  for (int p = 0; p < nbrVectors; ++p)
				    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
				}
			      ++TmpInteractionFactor;
			    }
			}
		    }
		  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
		    {
		      for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
			{			  
			  if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1] != 0)
			    {
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j][(i1 * Lim2) >> 2]);
			      for (int i2 = 0; i2 < Lim2; i2 += 2)
				{
				  Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
				  if (Index < Dim)
				    {
				      for (int p = 0; p < nbrVectors; ++p)
					vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
				    }
				  ++TmpInteractionFactor;
				}
			    }
			}
		    }
		}
	    }
	}
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	    {
	      for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
		{
		  Coefficient3 = particles->AsigmaAsigma(index, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
		  if (Coefficient3 != 0.0)
		    {
		      for (int p = 0; p < nbrVectors; ++p)
			tmpCoefficients[p] = Coefficient3 * vSources[p][index];
		      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			{
			  if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2] != 0)
			    {
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j][(i1 * Lim) >> 2]);
			      for (int i2 = 0; i2 < Lim; i2 += 2)
				{
				  Index = particles->AdsigmaAdsigma(TmpIndices[i2], TmpIndices[i2 + 1], sigma3, sigma3, Coefficient);
				  if (Index < Dim)
				    {
				      for (int p = 0; p < nbrVectors; ++p)
					vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
				    }
				  ++TmpInteractionFactor;
				}
			    }
			}
		      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			{
			  for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
			    {			  
			      if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] != 0)
				{
				  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
				  for (int i2 = 0; i2 < Lim2; i2 += 2)
				    {
				      Index = particles->AdsigmaAdsigma(TmpIndices2[i2], TmpIndices2[i2 + 1], sigma3, sigma4, Coefficient);
				      if (Index < Dim)
					{
					  for (int p = 0; p < nbrVectors; ++p)
					    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
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

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added  

inline void ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
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
	  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	    {
	      Coefficient3 = particles->AsigmaAsigma(index, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
	      if (Coefficient3 != 0.0)
		{
		  Coefficient4 = vSource[index];
		  Coefficient4 *= Coefficient3;
		  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
		    {
		      if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1] != 0)
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
		    }
		  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
		    {
		      for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
			{			  
			  if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1] != 0)
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
	}
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	    {
	      for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
		{
		  Coefficient3 = particles->AsigmaAsigma(index, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
		  if (Coefficient3 != 0.0)
		    {
		      Coefficient4 = vSource[index];
		      Coefficient4 *= Coefficient3;
		      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			{
			  if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2] != 0)
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
			}
		      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			{
			  for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
			    {			  
			      if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] != 0)
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

inline void ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
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
  for (int p = 0; p < nbrVectors; ++p)
    TmpSum[p] = 0.0;
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
      TmpIndices = this->IntraSectorIndicesPerSum[j];
      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
      TmpIndices2 = this->InterSectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 2)
	{
	  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	    {
	      Coefficient3 = particles->AsigmaAsigma(index, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
	      if (Coefficient3 != 0.0)
		{
		  for (int p = 0; p < nbrVectors; ++p)
		    tmpCoefficients[p] = Coefficient3 * vSources[p][index];
		  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
		    {
		      if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1] != 0)
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
		  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
		    {
		      for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
			{			  
			  if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1] != 0)
			    {
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j][(i1 * Lim2) >> 2]);
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
	    }
	}
      for (int i1 = 0; i1 < Lim2; i1 += 2)
	{
	  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	    {
	      for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
		{
		  Coefficient3 = particles->AsigmaAsigma(index, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
		  if (Coefficient3 != 0.0)
		    {
		      for (int p = 0; p < nbrVectors; ++p)
			tmpCoefficients[p] = Coefficient3 * vSources[p][index];
		      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			{
			  if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2] != 0)
			    {
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j][(i1 * Lim) >> 2]);
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
		      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			{
			  for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
			    {			  
			      if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] != 0)
				{
				  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
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

inline void ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
													     int* indexArray, Complex* coefficientArray, long& position)
{
//   cout << endl;
//   this->Particles->PrintState(cout, index);
//   cout << endl;
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
	      for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
		{
		  Coefficient3 = particles->AsigmaAsigma(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
		  if (Coefficient3 != 0.0)
		    {
		      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			{
			  if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1] != 0)
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
			}
		      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			{
			  for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
			    {			  
			      if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1] != 0)
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
	    }
	  for (int i1 = 0; i1 < Lim2; i1 += 2)
	    {
	      for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
		{
		  for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
		    {
		      Coefficient3 = particles->AsigmaAsigma(AbsoluteIndex, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
		      if (Coefficient3 != 0.0)
			{
			  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			    {
			      if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2] != 0)
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
			    }
			  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			    {
			      for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
				{			  
				  if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] != 0)
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
	      for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
		{
		  Coefficient3 = particles->AsigmaAsigma(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
		  if (Coefficient3 != 0.0)
		    {
		      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			{
			  if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1] != 0)
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
			}
		      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			{
			  for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
			    {			  
			      if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1] != 0)
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
	    }
	  for (int i1 = 0; i1 < Lim2; i1 += 2)
	    {
	      for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
		{
		  for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
		    {
		      Coefficient3 = particles->AsigmaAsigma(AbsoluteIndex, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
		      if (Coefficient3 != 0.0)
			{
			  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			    {
			      if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2] != 0)
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
			    }
			  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			    {
			      for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
				{			  
				  if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] != 0)
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
}

// core part of the PartialFastMultiplicationMemory method involving two-body term
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
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
		  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
		    {
		      Coefficient3 = particles->AsigmaAsigma(i, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
		      if (Coefficient3 != 0.0)
			{
			  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			    {
			      if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1] != 0)
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
			    }
			  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			    {
			      for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
				{			  
				  if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1] != 0)
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
		}
	      
	      for (int i1 = 0; i1 < Lim2; i1 += 2)
		{
		  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
		    {
		      for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
			{
			  Coefficient3 = particles->AsigmaAsigma(i, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
			  if (Coefficient3 != 0.0)
			    {
			      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
				{
				  if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2] != 0)
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
				}
			      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
				{
				  for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
				    {			  
				      if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] != 0)
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
		  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
		    {

		      Coefficient3 = particles->AsigmaAsigma(i, TmpIndices[i1], TmpIndices[i1 + 1], sigma1, sigma1);
		      if (Coefficient3 != 0.0)
			{
			  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			    {
			      if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1] != 0)
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
			    }
			  for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
			    {
			      for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
				{			  
				  if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1] != 0)
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
		}
	      for (int i1 = 0; i1 < Lim2; i1 += 2)
		{
		  for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
		    {
		      for (int sigma2 = sigma1 + 1; sigma2 < this->NbrInternalIndices; ++sigma2)
			{
			  Coefficient3 = particles->AsigmaAsigma(i, TmpIndices2[i1], TmpIndices2[i1 + 1], sigma1, sigma2);
			  if (Coefficient3 != 0.0)
			    {
			      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
				{
				  if (this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2] != 0)
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
				}
			      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
				{
				  for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
				    {			  
				      if (this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] != 0)
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
}

#endif
