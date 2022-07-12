////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                class of generic heisenberg hamiltonian written             //
//                             in fermionic language                          //
//                                                                            //
//                        last modification : 09/08/2015                      //
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


#ifndef PARTICLEONLATTICEREALSPACEFERMIONIZEDGENERICHEISENBERGHAMILTONIAN_H
#define PARTICLEONLATTICEREALSPACEFERMIONIZEDGENERICHEISENBERGHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEOnSphereFullHamiltonian.h"
#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian : public AbstractQHEOnSphereFullHamiltonian
{

 protected:

  // number of sites 
  int NbrSites;
  
  // parity of the particle number (0 for even or 1 for odd)
  int Parity;
  
  // true if periodic bounday conditions should be used
  bool PeriodicBoundaryConditionFlag;

  // value of the magetic field on each site
  double* HField;

  // coupling in front of the SxSx term
  double Jxx;
  // coupling in front of the SySy term
  double Jyy;
  // coupling in front of the SzSz term
  double Jzz;

  // optional array to store the hamiltonian diagonal elements
  double* DiagonalElements;

  // array that contains the number of sites connected to each site
  int* OneBodyGenericNbrConnectedSites;
  // array that contains the indices of the site connected to each site
  int** OneBodyGenericConnectedSites;
  // array that contains all one-body interaction factors
  double** OneBodyGenericInteractionFactors;

  // array that contains the number of sites connected to each site via a pairing term
  int* OneBodyPairingGenericNbrConnectedSites;
  // array that contains the indices of the site connected to each site via a pairing term
  int** OneBodyPairingGenericConnectedSites;
  // array that contains all one-body pairing interaction factors
  double** OneBodyPairingGenericInteractionFactors;

 public:

  // default constructor
  //
  ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // parity = parity of the particle number (0 for even or 1 for odd)
  // jxx = coupling in front of the SxSx term
  // jyy = coupling in front of the SySy term
  // jzz = coupling in front of the SzSz term
  // hField = array that gives the on-site magnetic field (zero pointer if none)
  // periodicBoundaryConditionFlag = true if periodic bounday conditions should be used
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian(ParticleOnSphere* particles, int parity, 
								    double jxx, double jyy, double jzz, 
								    double* hField, bool periodicBoundaryConditionFlag,
								    AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  virtual ~ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian();

  
 protected:
 
  // core part of the FastMultiplication method involving the one-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
							    int* indexArray, double* coefficientArray, long& position);

  // core part of the AddMultiply method involving the one-body interaction, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent,
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
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent,
						     int step, RealVector* vSources, RealVector* vDestinations, int nbrVectors);

  // core part of the PartialFastMultiplicationMemory method involving two-body term and one-body terms
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory);

  // evaluate the one body interaction factors
  //
  virtual void EvaluateOneBodyFactors ();

  // evaluate the two body interaction factors
  //
  virtual void EvaluateInteractionFactors ();
  
};


// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent,
															int step, RealVector& vSource, RealVector& vDestination)
{
  int Index;
  double Coefficient;
  if (this->OneBodyGenericInteractionFactors != 0)
    {
      for (int i = firstComponent; i < lastComponent; i += step)
	{ 
	  Coefficient = 0.0;
	  for (int j = 0; j < this->NbrSites; ++j)
	    {
	      int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
	      int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
	      double* TmpInteractionFactors = this->OneBodyGenericInteractionFactors[j];
	      int TmpPairingNbrConnectedSites = this->OneBodyPairingGenericNbrConnectedSites[j];
	      int* TmpPairingConnectedSites = this->OneBodyPairingGenericConnectedSites[j];
	      double* TmpPairingInteractionFactors = this->OneBodyPairingGenericInteractionFactors[j];
	      double TmpCoefficient;
	      Coefficient = particles->A(i,j);
	      if (Coefficient != 0.0 )
		{
		  for (int k = 0; k < TmpNbrConnectedSites; ++k)
		    {
		      TmpCoefficient = Coefficient;
		      Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient);
		      if (Index < particles->GetHilbertSpaceDimension())
			vDestination[Index] += TmpCoefficient * TmpInteractionFactors[k] * vSource[i];
		    }
		  for (int k = 0; k < TmpPairingNbrConnectedSites; ++k)
		    {
		      TmpCoefficient = Coefficient;
		      Index = particles->A(TmpPairingConnectedSites[k], TmpCoefficient);
		      if (Index < particles->GetHilbertSpaceDimension())
			vDestination[Index] += TmpCoefficient * TmpPairingInteractionFactors[k] * vSource[i];
		    }
		}
	      Coefficient = particles->Ad(i,j);
	      if (Coefficient != 0.0 )
		{
		  for (int k = 0; k < TmpPairingNbrConnectedSites; ++k)
		    {
		      TmpCoefficient = Coefficient;
		      Index = particles->Ad(TmpPairingConnectedSites[k], TmpCoefficient);
		      if (Index < particles->GetHilbertSpaceDimension())
			vDestination[Index] += TmpCoefficient * TmpPairingInteractionFactors[k] * vSource[i];
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

inline void ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent,
													 int step, RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  int Index;
  if (this->OneBodyGenericInteractionFactors != 0) 
    {
      for (int p = 0; p < nbrVectors; ++p)
	{
	  RealVector& TmpSourceVector = vSources[p];
	  RealVector& TmpDestinationVector = vDestinations[p];   
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
		  double* TmpInteractionFactors = this->OneBodyGenericInteractionFactors[j];
		  int TmpPairingNbrConnectedSites = this->OneBodyPairingGenericNbrConnectedSites[j];
		  int* TmpPairingConnectedSites = this->OneBodyPairingGenericConnectedSites[j];
		  double* TmpPairingInteractionFactors = this->OneBodyPairingGenericInteractionFactors[j];
		  double TmpCoefficient;
		  Coefficient = particles->A(i,j);
		  if (Coefficient != 0.0 )
		    {
		      for (int k = 0; k < TmpNbrConnectedSites; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient);
			  if (Index < particles->GetHilbertSpaceDimension())
			    TmpDestinationVector[Index] += TmpCoefficient * TmpInteractionFactors[k] * TmpSourceVector[i];
			}
		      for (int k = 0; k < TmpPairingNbrConnectedSites; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->A(TmpPairingConnectedSites[k], TmpCoefficient);
			  if (Index < particles->GetHilbertSpaceDimension())
			    TmpDestinationVector[Index] += TmpCoefficient * TmpPairingInteractionFactors[k] * TmpSourceVector[i];
			}
		    }
		  Coefficient = particles->Ad(i,j);
		  if (Coefficient != 0.0 )
		    {
		      for (int k = 0; k < TmpPairingNbrConnectedSites; ++k)
			{
			  TmpCoefficient = Coefficient;
			  Index = particles->Ad(TmpPairingConnectedSites[k], TmpCoefficient);
			  if (Index < particles->GetHilbertSpaceDimension())
			    TmpDestinationVector[Index] += TmpCoefficient * TmpPairingInteractionFactors[k] * TmpSourceVector[i];
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
	  RealVector& TmpSourceVector = vSources[p];
	  RealVector& TmpDestinationVector = vDestinations[p];   
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

inline void ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
															       int* indexArray, double* coefficientArray, long& position)
{
  if (this->OneBodyGenericInteractionFactors == 0)
    return;

  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  int Index;

  if (this->HermitianSymmetryFlag == false)
    {
      for (int j = 0; j < this->NbrSites; ++j)
	{
	  int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
	  int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
	  double* TmpInteractionFactors = this->OneBodyGenericInteractionFactors[j];
	  int TmpPairingNbrConnectedSites = this->OneBodyPairingGenericNbrConnectedSites[j];
	  int* TmpPairingConnectedSites = this->OneBodyPairingGenericConnectedSites[j];
	  double* TmpPairingInteractionFactors = this->OneBodyPairingGenericInteractionFactors[j];
	  double TmpCoefficient;
	  Coefficient = particles->A(index, j);
	  if (Coefficient != 0.0)
	    {
	      for (int k = 0; k < TmpNbrConnectedSites; ++k)
		{
		  TmpCoefficient = Coefficient;
		  Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = TmpCoefficient * TmpInteractionFactors[k];
		      ++position;
		    }
		}
	      for (int k = 0; k < TmpPairingNbrConnectedSites; ++k)
		{
		  TmpCoefficient = Coefficient;
		  Index = particles->A(TmpPairingConnectedSites[k], TmpCoefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = TmpCoefficient * TmpPairingInteractionFactors[k];
		      ++position;
		    }
		}
	    }
	  Coefficient = particles->Ad(index, j);
	  if (Coefficient != 0.0)
	    {
	      for (int k = 0; k < TmpPairingNbrConnectedSites; ++k)
		{
		  TmpCoefficient = Coefficient;
		  Index = particles->Ad(TmpPairingConnectedSites[k], TmpCoefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = TmpCoefficient * TmpPairingInteractionFactors[k];
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
	  double* TmpInteractionFactors = this->OneBodyGenericInteractionFactors[j];
	  int TmpPairingNbrConnectedSites = this->OneBodyPairingGenericNbrConnectedSites[j];
	  int* TmpPairingConnectedSites = this->OneBodyPairingGenericConnectedSites[j];
	  double* TmpPairingInteractionFactors = this->OneBodyPairingGenericInteractionFactors[j];
	  double TmpCoefficient;
	  Coefficient = particles->A(index, j);
	  if (Coefficient != 0.0)
	    {
	      for (int k = 0; k < TmpNbrConnectedSites; ++k)
		{
		  TmpCoefficient = Coefficient;
		  Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient);
		  if (Index <= index)
		    {
		      indexArray[position] = Index;
		      if (Index == index)
			{
			  coefficientArray[position] = 0.5 * TmpCoefficient * TmpInteractionFactors[k];
			}
		      else
			{
			  coefficientArray[position] = TmpCoefficient * TmpInteractionFactors[k];
			}
		      ++position;
		    }
		}
	      for (int k = 0; k < TmpPairingNbrConnectedSites; ++k)
		{
		  TmpCoefficient = Coefficient;
		  Index = particles->A(TmpPairingConnectedSites[k], TmpCoefficient);
		  if (Index <= index)
		    {
		      indexArray[position] = Index;
		      if (Index == index)
			{
			  coefficientArray[position] = 0.5 * TmpCoefficient * TmpPairingInteractionFactors[k];
			}
		      else
			{
			  coefficientArray[position] = TmpCoefficient * TmpPairingInteractionFactors[k];
			}
		      ++position;
		    }
		}
	    }
	  Coefficient = particles->Ad(index, j);
	  if (Coefficient != 0.0)
	    {
	      for (int k = 0; k < TmpPairingNbrConnectedSites; ++k)
		{
		  TmpCoefficient = Coefficient;
		  Index = particles->Ad(TmpPairingConnectedSites[k], TmpCoefficient);
		  if (Index <= index)
		    {
		      indexArray[position] = Index;
		      if (Index == index)
			{
			  coefficientArray[position] = 0.5 * TmpCoefficient * TmpPairingInteractionFactors[k];
			}
		      else
			{
			  coefficientArray[position] = TmpCoefficient * TmpPairingInteractionFactors[k];
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

inline void ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory)
{
  if (this->OneBodyGenericInteractionFactors == 0) 
    return;
  int Index;
  double Coefficient = 0.0;
  int Dim = particles->GetHilbertSpaceDimension();
  
  if (this->HermitianSymmetryFlag == false)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrSites; ++j)
	    {
	      int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
	      int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
	      int TmpPairingNbrConnectedSites = this->OneBodyPairingGenericNbrConnectedSites[j];
	      int* TmpPairingConnectedSites = this->OneBodyPairingGenericConnectedSites[j];
	      double TmpCoefficient;
	      Coefficient = particles->A(i, j);
	      if (Coefficient != 0.0)
		{
		  for (int k = 0; k < TmpNbrConnectedSites; ++k)
		    {
		      TmpCoefficient = Coefficient;
		      Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		  for (int k = 0; k < TmpPairingNbrConnectedSites; ++k)
		    {
		      TmpCoefficient = Coefficient;
		      Index = particles->A(TmpPairingConnectedSites[k], TmpCoefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		}
	      Coefficient = particles->Ad(i, j);
	      if (Coefficient != 0.0)
		{
		  for (int k = 0; k < TmpPairingNbrConnectedSites; ++k)
		    {
		      TmpCoefficient = Coefficient;
		      Index = particles->Ad(TmpPairingConnectedSites[k], TmpCoefficient);
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
  else
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrSites; ++j)
	    {
	      int TmpNbrConnectedSites = this->OneBodyGenericNbrConnectedSites[j];
	      int* TmpConnectedSites = this->OneBodyGenericConnectedSites[j];
	      int TmpPairingNbrConnectedSites = this->OneBodyPairingGenericNbrConnectedSites[j];
	      int* TmpPairingConnectedSites = this->OneBodyPairingGenericConnectedSites[j];
	      double TmpCoefficient;
	      Coefficient = particles->A (i, j);
	      if (Coefficient != 0.0)
		{
		  for (int k = 0; k < TmpNbrConnectedSites; ++k)
		    {
		      TmpCoefficient = Coefficient;
		      Index = particles->Ad(TmpConnectedSites[k], TmpCoefficient);
		      if (Index <= i)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		  for (int k = 0; k < TmpPairingNbrConnectedSites; ++k)
		    {
		      TmpCoefficient = Coefficient;
		      Index = particles->A(TmpPairingConnectedSites[k], TmpCoefficient);
		      if (Index <= i)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		}
	      Coefficient = particles->Ad (i, j);
	      if (Coefficient != 0.0)
		{
		  for (int k = 0; k < TmpPairingNbrConnectedSites; ++k)
		    {
		      TmpCoefficient = Coefficient;
		      Index = particles->Ad(TmpPairingConnectedSites[k], TmpCoefficient);
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
