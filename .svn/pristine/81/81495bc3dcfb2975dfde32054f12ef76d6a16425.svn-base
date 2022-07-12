////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//                   class of Hubbard hamiltonian associated                  //
//                to particles on a lattice with spin                         //
//                                                                            //
//                        last modification : 19/06/2014                      //
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


#ifndef PARTICLEONLATTICEWITHSPINKITAEVHEISENBERGHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINKITAEVHEISENBERGHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian : public ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian
{

 protected:
  
  // number of sites
  int NbrSite;
  
  //number of bonds in the file describing the bonds
  int NbrBonds;
  //sites in first column in file describing the bonds
  int* SitesA;
  //sites in second column in file describing the bonds
  int* SitesB;
  //third column in file describing the bonds: nature of each bond
  int* Bonds;
  
  //array that contains the information about the lattice geometry
  //first index = index of the site under consideration
  //second index 0 = index of the nearest neighbor linked with x link (= NbrSite if no x link)
  //second index 1 = index of the nearest neighbor linked with y link (= NbrSite if no y link)
  //second index 2 = index of the nearest neighbor linked with z link (= NbrSite if no z link)
  int** MapNearestNeighborBonds;
    
  // multiplicative factor in front of the isotropic nearest neighbor hopping term
  double KineticFactorIsotropic;
  // multiplicative factor in front of the anisotropic nearest neighbor hopping term
  double KineticFactorAnisotropic;
  
  //strength of the isotropic nearest neighbor spin interaction
  double J1Factor;
  //strength of the anisotropic nearest neighbor spin interaction
  double J2Factor;

  // array that contains all one-body interaction factors for particles with spin up
  double** OneBodyGenericInteractionFactorsupup;
  // array that contains all one-body interaction factors for particles with spin down
  double** OneBodyGenericInteractionFactorsdowndown;
  // array that contains all one-body interaction factors for tunnelling terms for particles with different spin
  Complex** OneBodyGenericInteractionFactorsupdown;
  
 public:

  // default constructor
  //
  ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSites = number of sites
  // kineticFactorIntra = multiplicative factor in front of the intraspin kinetic term
  // uPotential = Hubbard potential strength
  // j1Factor = strength of the isotropic nearest neighbor spin interaction
  // j2Factor = strength of the anisotropic nearest neighbor spin interaction
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSite, int** geometryDescription, int nbrBonds, double kineticFactorIsotropic, double kineticFactorAnisotropic, double uPotential, double j1Factor, double j2Factor, AbstractArchitecture* architecture, long memory = -1);
  
  // constructor from the explicit bond description
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSite = number of sites
  // nbrBonds = number of bonds
  // sitesA = array of A sites for each bond
  // sitesB = array of B sites for each bond
  // bondTypes = array that describe each type of bond (0 for x, 1 fo y, 2 for z)
  // kineticFactor = multiplicative factor in front of the kinetic term
  // uPotential = Hubbard potential strength
  // j1Factor = strength of the isotropic nearest neighbor spin interaction
  // j2Factor = strength of the anisotropic nearest neighbor spin interaction
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)  
  ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSite, 
						       int nbrBonds,  int* sitesA, int* sitesB, int* bondTypes, 
						       double kineticFactorIsotropic, double kineticFactorAnisotropic, 
						       double uPotential, double j1Factor, double j2Factor, 
						       AbstractArchitecture* architecture, long memory = -1l);

  // destructor
  //
  virtual ~ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian();

  
  //get the index of the nearest neighbor linked with bond x
  //
  // i = index of the site under consideration
  // return value = index of the site linked to i with x bond
  int GetIndexNearestNeighborXBond(int i);
  
  //get the index of the nearest neighbor linked with bond y
  //
  // i = index of the site under consideration
  // return value = index of the site linked to i with y bond
  int GetIndexNearestNeighborYBond(int i);
  
  //get the index of the nearest neighbor linked with bond z
  //
  // i = index of the site under consideration
  // return value = index of the site linked to i with z bond
  int GetIndexNearestNeighborZBond(int i);
  
  
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
//   virtual void EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory);

  // core part of the PartialFastMultiplicationMemory method involving two-body term and one-body terms
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory);

  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

    
  //fill up the positions of the three types of nearest neighbor bonds
  //
  virtual void PlotMapNearestNeighborBonds();
  
  //find the type of bond that links neighboring sites i and j
  //
  // i = index of the first site
  // j = index of the second site
  //return value = 0 (x bond), 1 (y bond), 2 (z bond)
  virtual int FindBondType(int i, int j);
  

};

// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
												      int step, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Index;
  double Coefficient;
  if (this->HermitianSymmetryFlag == false)
    {
      if (this->OneBodyGenericInteractionFactorsupup != 0)
	if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	  {
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		for (int j = 0; j < this->NbrSite; ++j) 
		  {
		    for (int k = 0; k < 3; ++k)
		      {
			int j2 = this->MapNearestNeighborBonds[j][k];
			if (j2 < this->NbrSite)
			  {
			    Index = particles->AduAu(i , j, j2, Coefficient);
			    if (Index < particles->GetHilbertSpaceDimension())
			      vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k] * vSource[i];
			    
			    Index = particles->AddAd(i , j, j2, Coefficient);
			    if (Index < particles->GetHilbertSpaceDimension())
			      vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k] * vSource[i];		  
			  }
		      }
		  }
	      }
	  }
	else
	  {
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		Coefficient = 0.0;
		for (int j = 0; j < this->NbrSite; ++j) 
		  {
		    for (int k = 0; k < 3; ++k)
		      {
			int j2 = this->MapNearestNeighborBonds[j][k];
			if (j2 < this->NbrSite)
			  {
			    Index = particles->AduAu(i, j, j2, Coefficient);
			    if (Index < particles->GetHilbertSpaceDimension())
			      vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k] * vSource[i];
			  }
		      }
		  }
	      }
	  }
      else
	{
	  if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	    {
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  Coefficient = 0.0;
		  for (int j = 0; j < this->NbrSite; ++j) 
		    {
		      for (int k = 0; k < 3; ++k)
			{
			  int j2 = this->MapNearestNeighborBonds[j][k];
			  if (j2 < this->NbrSite)
			    {
			      Index = particles->AddAd(i, j, j2, Coefficient);
			      if (Index < particles->GetHilbertSpaceDimension())
				vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k] * vSource[i];		  
			    }
			}
		    }
		}
	    }	
	}
      for (int i = firstComponent; i < lastComponent; i += step)
	vDestination[i] += this->HamiltonianShift * vSource[i];
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  double Coefficient;
	  Complex Source;
	  int Dim = particles->GetHilbertSpaceDimension();
	  int Index;
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      Source = vSource[i];
	      for (int j = 0; j < this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAu(i, j, j2, Coefficient);
			  if (Index < Dim)
			    vDestination[Index] += (Coefficient * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k])) * Source;
			  Index = particles->AduAd(i, j, j2, Coefficient);
			  if (Index < Dim)
			    vDestination[Index] += (Coefficient * this->OneBodyGenericInteractionFactorsupdown[j][k]) * Source;
			}
		    }
		}
	    }
	}
    }
  else
    {
      if (this->OneBodyGenericInteractionFactorsupup != 0)
	if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	  {
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		for (int j = 0; j < this->NbrSite; ++j) 
		  {
		    for (int k = 0; k < 3; ++k)
		      {
			int j2 = this->MapNearestNeighborBonds[j][k];
			if (j2 < this->NbrSite)
			  {
			    Index = particles->AduAu(i, j, j2, Coefficient);
			    if (Index <= i)
			      {
				vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k] * vSource[i];
				if (Index < i)
				  {
				    vDestination[i] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsupup[j][k]) * vSource[Index];
				  }
			      }			    
			    Index = particles->AddAd(i , j, j2, Coefficient);
			    if (Index <= i)
			      {
				vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k] * vSource[i];		  
				if (Index < i)
				  {
				    vDestination[i] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsdowndown[j][k]) * vSource[Index];		  
				  }
			      }
			  }
		      }
		  }
	      }
	  }
	else
	  {
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		Coefficient = 0.0;
		for (int j = 0; j < this->NbrSite; ++j) 
		  {
		    for (int k = 0; k < 3; ++k)
		      {
			int j2 = this->MapNearestNeighborBonds[j][k];
			if (j2 < this->NbrSite)
			  {
			    Index = particles->AduAu(i, j, j2, Coefficient);
			    if (Index <= i)
			      {
				vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k] * vSource[i];
				if (Index < i)
				  {
				    vDestination[i] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsupup[j][k]) * vSource[Index];
				  }
			      }
			  }
		      }
		  }
	      }
	  }
      else
	{
	  if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	    {
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  Coefficient = 0.0;
		  for (int j = 0; j < this->NbrSite; ++j) 
		    {
		      for (int k = 0; k < 3; ++k)
			{
			  int j2 = this->MapNearestNeighborBonds[j][k];
			  if (j2 < this->NbrSite)
			    {
			      Index = particles->AddAd(i, j, j2, Coefficient);
			      if (Index <= i)
				{
				  vDestination[Index] += Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k] * vSource[i];		  
				  if (Index < i)
				    {
				      vDestination[i] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsdowndown[j][k]) * vSource[Index];		  
				    }
				}
			    }
			}
		    }
		}
	    }	
	}
      for (int i = firstComponent; i < lastComponent; i += step)
	vDestination[i] += this->HamiltonianShift * vSource[i];
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  double Coefficient;
	  Complex Source;
	  int Dim = particles->GetHilbertSpaceDimension();
	  int Index;
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      Source = vSource[i];
	      for (int j = 0; j < this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAu(i, j, j2, Coefficient);
			  if (Index <= i)
			    {
			      vDestination[Index] += (Coefficient * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k])) * Source;
			      if (Index < i)
				{
				  vDestination[i] += (Coefficient * this->OneBodyGenericInteractionFactorsupdown[j][k]) * vSource[Index];
				}
			    }
			  Index = particles->AduAd(i, j, j2, Coefficient);
			  if (Index <= i)
			    {
			      vDestination[Index] += (Coefficient * this->OneBodyGenericInteractionFactorsupdown[j][k]) * Source;
			      if (Index < i)
				{
				  vDestination[i] += (Coefficient * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k])) * vSource[Index];
				}
			    }
			}
		    }
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

inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
												      int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  int Index;
  if (this->HermitianSymmetryFlag == false)
    {
      if (this->OneBodyGenericInteractionFactorsupup != 0)
	{
	  if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		{
		  ComplexVector& TmpSourceVector = vSources[p];
		  ComplexVector& TmpDestinationVector = vDestinations[p];
		  
		  for (int i = firstComponent; i < lastComponent; i += step)
		    {
		      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		      for (int j = 0; j < this->NbrSite; ++j)
			{
			  for (int k = 0; k < 3; ++k)
			    {
			      int j2 = this->MapNearestNeighborBonds[j][k];
			      if (j2 < this->NbrSite)
				{
				  Index = particles->AduAu(i, j, j2, Coefficient);
				  if (Index < Dim)
				    TmpDestinationVector[Index] += Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k] * TmpSourceVector[i];
				  
				  Index = particles->AddAd(i, j, j2, Coefficient);
				  if (Index < Dim)
				    TmpDestinationVector[Index] += Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k] * TmpSourceVector[i];
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
		  ComplexVector& TmpSourceVector = vSources[p];
		  ComplexVector& TmpDestinationVector = vDestinations[p];  
		  for (int i = firstComponent; i < lastComponent; i += step)
		    {
		      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		      for (int j = 0; j < this->NbrSite; ++j)
			{
			  for (int k = 0; k < 3; ++k)
			    {
			      int j2 = this->MapNearestNeighborBonds[j][k];
			      if (j2 < this->NbrSite)
				{
				  Index = particles->AduAu(i, j, j2, Coefficient);
				  if (Index < Dim)
				    TmpDestinationVector[Index] += Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k] * TmpSourceVector[i];
				}
			    }
			}
		    }
		}
	    }
	}
      else
	{
	  if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		{
		  ComplexVector& TmpSourceVector = vSources[p];
		  ComplexVector& TmpDestinationVector = vDestinations[p];  
		  for (int i = firstComponent; i < lastComponent; i += step)
		    {
		      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		      for (int j = 0; j < this->NbrSite; ++j)
			{
			  for (int k = 0; k < 3; ++k)
			    {
			      int j2 = this->MapNearestNeighborBonds[j][k];
			      if (j2 < this->NbrSite)
				{
				  Index = particles->AddAd(i, j, j2, Coefficient);
				  if (Index < Dim)
				    TmpDestinationVector[Index] += Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k] * TmpSourceVector[i];
				}
			    }
			}
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
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      for (int j = 0; j < this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAu(i + this->PrecalculationShift, j, j2, Coefficient);
			  if (Index < Dim)
			    {
			      for (int p = 0; p < nbrVectors; ++p)
				vDestinations[p][Index] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k]) * vSources[p][i];
			    }
			  Index = particles->AduAd(i + this->PrecalculationShift, j, j2, Coefficient);
			  if (Index < Dim)
			    {
			      for (int p = 0; p < nbrVectors; ++p)
				vDestinations[p][Index] += Coefficient * this->OneBodyGenericInteractionFactorsupdown[j][k] * vSources[p][i];
			    }
			}
		    }
		}
	    }
	}
    }
  else
    {
      if (this->OneBodyGenericInteractionFactorsupup != 0) 
	{
	  if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		{
		  ComplexVector& TmpSourceVector = vSources[p];
		  ComplexVector& TmpDestinationVector = vDestinations[p];
		  
		  for (int i = firstComponent; i < lastComponent; i += step)
		    { 
		      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		      for (int j = 0; j < this->NbrSite; ++j) 
			{
			  for (int k = 0; k < 3; ++k)
			    {
			      int j2 = this->MapNearestNeighborBonds[j][k];
			      if (j2 < this->NbrSite)
				{
				  Index = particles->AduAu(i, j, j2, Coefficient);
				  if (Index <= i)
				    {
				      TmpDestinationVector[Index] += Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k] * TmpSourceVector[i];
				      if (Index < i)
					{
					  TmpDestinationVector[i] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsupup[j][k]) * TmpSourceVector[Index];
					}
				    }
				  
				  Index = particles->AddAd(i, j, j2, Coefficient);
				  if (Index <= i)
				    {
				      TmpDestinationVector[Index] += Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k] * TmpSourceVector[i];
				      if (Index < i)
					{
					  TmpDestinationVector[i] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsdowndown[j][k]) * TmpSourceVector[Index];
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
	      for (int p = 0; p < nbrVectors; ++p)
		{
		  ComplexVector& TmpSourceVector = vSources[p];
		  ComplexVector& TmpDestinationVector = vDestinations[p];   
		  for (int i = firstComponent; i < lastComponent; i += step)
		    { 
		      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		      for (int j = 0; j < this->NbrSite; ++j) 
			{
			  for (int k = 0; k < 3; ++k)
			    {
			      int j2 = this->MapNearestNeighborBonds[j][k];
			      if (j2 < this->NbrSite)
				{
				  Index = particles->AduAu(i, j, j2, Coefficient);
				  if (Index <= i)
				    {
				      TmpDestinationVector[Index] += Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k] * TmpSourceVector[i];
				      if (Index < i)
					{
					  TmpDestinationVector[i] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsupup[j][k]) * TmpSourceVector[Index];
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
	  if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		{
		  ComplexVector& TmpSourceVector = vSources[p];
		  ComplexVector& TmpDestinationVector = vDestinations[p];   
		  for (int i = firstComponent; i < lastComponent; i += step)
		    { 
		      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		      for (int j = 0; j < this->NbrSite; ++j) 
			{
			  for (int k = 0; k < 3; ++k)
			    {
			      int j2 = this->MapNearestNeighborBonds[j][k];
			      if (j2 < this->NbrSite)
				{
				  Index = particles->AddAd(i, j, j2, Coefficient);
				  if (Index <= i)
				    {
				      TmpDestinationVector[Index] += Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k] * TmpSourceVector[i];
				      if (Index < i)
					{
					  TmpDestinationVector[i] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsdowndown[j][k]) * TmpSourceVector[Index];
					}
				    }
				}
			    }
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
	}
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      for (int j = 0; j < this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAu(i, j, j2, Coefficient);
			  if (Index <= i)
			    {
			      if (Index < i)
				{
				  for (int p = 0; p < nbrVectors; ++p)
				    {
				      vDestinations[p][Index] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k]) * vSources[p][i];
				      vDestinations[p][i] += Coefficient * this->OneBodyGenericInteractionFactorsupdown[j][k] * vSources[p][Index];
				    }
				}
			      else
				{
				  for (int p = 0; p < nbrVectors; ++p)
				    vDestinations[p][Index] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k]) * vSources[p][i];
				}
			    }
			  Index = particles->AduAd(i, j, j2, Coefficient);
			  if (Index <= i)
			    {
			      if (Index < i)
				{
				  for (int p = 0; p < nbrVectors; ++p)
				    {
				      vDestinations[p][Index] += Coefficient * this->OneBodyGenericInteractionFactorsupdown[j][k] * vSources[p][i];
				      vDestinations[p][i] += Coefficient * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k]) * vSources[p][Index];
				    }
				}
			      else
				{
				  for (int p = 0; p < nbrVectors; ++p)
				    vDestinations[p][Index] += Coefficient * this->OneBodyGenericInteractionFactorsupdown[j][k] * vSources[p][i];
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

inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
													     int* indexArray, Complex* coefficientArray, long& position)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  int Index;
//   cout << position << endl;
  if (this->HermitianSymmetryFlag == false)
    {
      if ((this->OneBodyGenericInteractionFactorsdowndown != 0) && (this->OneBodyGenericInteractionFactorsupup != 0))
	{
	  for (int j = 0; j < this->NbrSite; ++j)
	    {
	      for (int k = 0; k < 3; ++k)
		{
		  int j2 = this->MapNearestNeighborBonds[j][k];
		  if (j2 < this->NbrSite)
		    {
		      Index = particles->AduAu(index + this->PrecalculationShift, j, j2, Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k];
			  ++position;
			}
		      Index = particles->AddAd(index + this->PrecalculationShift, j, j2, Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k];
			  ++position;
			}
		    }
		}
	    }
	}
      else
	{
	  if (this->OneBodyGenericInteractionFactorsupup != 0)
	    {
	      for (int j = 0; j < this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AduAu(index + this->PrecalculationShift, j, j2, Coefficient);
			  if (Index < Dim)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k];
			      ++position;
			    }
			}
		    }
		}
	    }
	  else
	    if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	      {
		for (int j = 0; j < this->NbrSite; ++j)
		  {
		    for (int k = 0; k < 3; ++k)
		      {
			int j2 = this->MapNearestNeighborBonds[j][k];
			if (j2 < this->NbrSite)
			  {
			    Index = particles->AddAd(index + this->PrecalculationShift, j, j2, Coefficient);
			    if (Index < Dim)
			      {
				indexArray[position] = Index;
				coefficientArray[position] = Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k];
				++position;
			      }
			  }
		      }
		  }
	      } 
	}
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  for (int j = 0; j < this->NbrSite; ++j)
	    {
	      for (int k = 0; k < 3; ++k)
		{
		  int j2 = this->MapNearestNeighborBonds[j][k];
		  if (j2 < this->NbrSite)
		    {
		      Index = particles->AddAu(index + this->PrecalculationShift, j, j2, Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k]);
			  ++position;
			}
		      Index = particles->AduAd(index + this->PrecalculationShift, j, j2, Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * this->OneBodyGenericInteractionFactorsupdown[j][k];
			  ++position;
			}
		    }
		}
	    }
	}
    }
  else
    {
      int AbsoluteIndex = index + this->PrecalculationShift;
      if ((this->OneBodyGenericInteractionFactorsdowndown != 0) && (this->OneBodyGenericInteractionFactorsupup != 0))
	{
	  for (int j = 0; j < this->NbrSite; ++j)
	    {
	      for (int k = 0; k < 3; ++k)
		{
		  int j2 = this->MapNearestNeighborBonds[j][k];
		  if (j2 < this->NbrSite)
		    {
		      Index = particles->AduAu(AbsoluteIndex, j, j2, Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k];
			  ++position;
			}
		      Index = particles->AddAd(AbsoluteIndex, j, j2, Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k];
			  ++position;
			}
		    }
		}
	    }
	}
      else
	{
	  if (this->OneBodyGenericInteractionFactorsupup != 0)
	    {
	      for (int j = 0; j < this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AduAu(AbsoluteIndex, j, j2, Coefficient);
			  if (Index <= AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * this->OneBodyGenericInteractionFactorsupup[j][k];
			      ++position;
			    }
			}
		    }
		}
	    }
	  else
	    if (this->OneBodyGenericInteractionFactorsdowndown != 0)
	      {
		for (int j = 0; j < this->NbrSite; ++j)
		  {
		    for (int k = 0; k < 3; ++k)
		      {
			int j2 = this->MapNearestNeighborBonds[j][k];
			if (j2 < this->NbrSite)
			  {
			    Index = particles->AddAd(AbsoluteIndex, j, j2, Coefficient);
			    if (Index <= AbsoluteIndex)
			      {
				indexArray[position] = Index;
				coefficientArray[position] = Coefficient * this->OneBodyGenericInteractionFactorsdowndown[j][k];
				++position;
			      }
			  }
		      }
		  }
	      } 
	}
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  for (int j = 0; j < this->NbrSite; ++j)
	    {
	      for (int k = 0; k < 3; ++k)
		{
		  int j2 = this->MapNearestNeighborBonds[j][k];
		  if (j2 < this->NbrSite)
		    {
		      Index = particles->AddAu(AbsoluteIndex, j, j2, Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Conj(this->OneBodyGenericInteractionFactorsupdown[j][k]);
			  ++position;
			}
		      Index = particles->AduAd(AbsoluteIndex, j, j2, Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * this->OneBodyGenericInteractionFactorsupdown[j][k];
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

inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;

  int Dim = particles->GetHilbertSpaceDimension();
  
  if (this->HermitianSymmetryFlag == false)
    {
      if (this->OneBodyGenericInteractionFactorsupup != 0) 
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AduAu(i, j, j2, Coefficient);
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
      
      if (this->OneBodyGenericInteractionFactorsdowndown != 0) 
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAd(i, j, j2, Coefficient);
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
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAu(i, j, j2, Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			  Index = particles->AduAd(i, j, j2, Coefficient);
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
      if (this->OneBodyGenericInteractionFactorsupup != 0) 
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AduAu(i, j, j2, Coefficient);
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
      
      if (this->OneBodyGenericInteractionFactorsdowndown != 0) 
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAd(i, j, j2, Coefficient);
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
      if (this->OneBodyGenericInteractionFactorsupdown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
			{
			  Index = particles->AddAu(i, j, j2, Coefficient);
			  if (Index <= i)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			    }
			  Index = particles->AduAd(i, j, j2, Coefficient);
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

//find the type of bond that links neighboring sites i and j
//
// i = index of the first site
// j = index of the second site
//return value = 0 (x bond), 1 (y bond), 2 (z bond)

inline int ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::FindBondType(int i, int j)
{
//  cout << i << " " << j << endl;
  int diff = i - j;
  if (this->NbrBonds == 0)
  {
  //open chain of hexagons
  if ((this->NbrSite % 4) == 2)
    {
      if (diff < 0)
	{
	  diff = j;
	  j = i;
	  i = diff;
	  diff = i - j;
	}
      if (diff == 1)
	return 2;
      switch (i & 3)
	{
	case 0:
	  return 0;
	  break;
	case 3:
	  return 0;
	  break;
	case 1:
	  return 1;
	  break;
	case 2:
	  return 1;
	  break;
	}
      return -1;
    }
  else
    {
      //periodic chain of hexagons
      if (((this->NbrSite % 4) == 0) && (this->NbrSite > 4))
      {
	if (diff < 0)
	{
	  diff = j;
	  j = i;
	  i = diff;
	  diff = i - j;
	}
      if (diff == 1)
	return 2;
      if (diff == this->NbrSite - 2)
      {
	return j;
      }
      switch (i & 3)
	{
	case 0:
	  return 0;
	  break;
	case 3:
	  return 0;
	  break;
	case 1:
	  return 1;
	  break;
	case 2:
	  return 1;
	  break;
	}
      return -1;
      }
      else
      {
      if (abs(diff) != 1)
	return -1;
      else
	{
	  if ((i % 2) == 0)
	    if (diff == -1)
	      return 1;
	    else
	      return 0;
	  
	  else
	    if (diff == 1)
	      return 1;
	    else
	      return 0;
	}      
      }
    }
  }
  else
  {
    for (int k = 0; k < this->NbrBonds; ++k)
      if (((this->SitesA[k] == i) && (this->SitesB[k] == j)) || ((this->SitesA[k] == j) && (this->SitesB[k] == i)))
	return this->Bonds[k];
    return-1;
  }
}

//get the index of the nearest neighbor linked with bond x
//
// i = index of the site under consideration
// return value = index of the site linked to i with x bond
inline int ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::GetIndexNearestNeighborXBond(int i)
{
  return this->MapNearestNeighborBonds[i][0];
}


//get the index of the nearest neighbor linked with bond y
//
// i = index of the site under consideration
// return value = index of the site linked to i with y bond
inline int ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::GetIndexNearestNeighborYBond(int i)
{
  return this->MapNearestNeighborBonds[i][1];
}


//get the index of the nearest neighbor linked with bond z
//
// i = index of the site under consideration
// return value = index of the site linked to i with z bond
inline int ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::GetIndexNearestNeighborZBond(int i)
{
  return this->MapNearestNeighborBonds[i][2];
}
#endif
