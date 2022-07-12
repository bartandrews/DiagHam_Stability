////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class of hamiltonian associated to particles on a cylinder with        //
//        SU(2) spin with opposite magnetic field for each species            //
//      a generic interaction defined by its pseudopotential, pairing         //
//               and translation breaking one body potentials                 //
//                                                                            //
//                        last modification : 19/09/2016                      //
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


#ifndef PARTICLEONCYLINDERWITHSPINTIMEREVERSALSYMMETRYGENERICHAMILTONIANANDPAIRINGALLMOMENTA_H
#define PARTICLEONCYLINDERWITHSPINTIMEREVERSALSYMMETRYGENERICHAMILTONIANANDPAIRINGALLMOMENTA_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinChernInsulatorHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class AbstractArchitecture;


class ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta : public ParticleOnLatticeWithSpinChernInsulatorHamiltonian
{

  friend class QHEParticlePrecalculationOperation;

 protected:

  // Number of Pseudopotential for up-up interaction
  int NbrPseudopotentialsUpUp;
  // pseudopotential coefficients for up-up interaction
  double* PseudopotentialsUpUp;
  // Number of Pseudopotential for down-down interaction
  int NbrPseudopotentialsDownDown;
  // pseudopotential coefficients for down-down interaction
  double* PseudopotentialsDownDown;
  // Number of Pseudopotential for up-down interaction
  int NbrPseudopotentialsUpDown;
  // pseudopotential coefficients for up-down interaction
  double* PseudopotentialsUpDown;

  // ratio between the width in the x direction and the width in the y direction
  double Ratio;
  // ratio between the width in the y direction and the width in the x direction
  double InvRatio;

   // maxixum monentum transfer that can appear in a one body operator
   int MaximumMomentumTransfer;

  // off-diagonal contribution of the one-body potential for particles with spin up, the first entry is the annihilation index, the second entry is the momentum tranfer
  Complex** OneBodyOffDiagonalInteractionFactorsupup;
  // off-diagonal contribution of the one-body potential for particles with spin down, the first entry is the annihilation index, the second entry is the momentum tranfer
  Complex** OneBodyOffDiagonalInteractionFactorsdowndown;
  // array that contains all one-body interaction factors for the pairing term
  Complex* OneBodyInteractionFactorsPairing;
  // off diagonal contribution of the one-body pairing term, the first entry is the index of the rightmost creation operator, the second entry is the momentum tranfer
  Complex** OneBodyOffDiagonalInteractionFactorsPairing; 

  // factor in front of the charging energy (i.e 1/(2C))
  double ChargingEnergy;
  // avearge number of particles in the system
  double AverageNumberParticles;

 public:

  ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // lzmax = maximum Lz value reached by a particle in the state
  // ratio = ratio between the width in the x direction and the width in the y direction
  // maxMomentumTransfer = maxixum monentum transfer that can appear in a one body operator
  // pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
  //                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
  // onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
  // onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
  // onebodyOffDiagonalPotentialUpUp = off-diagonal contribution of the one-body potential for particles with spin up, the first entry is the annihilation index, 
  //                                   the second entry is the momentum tranfer
  // onebodyOffDiagonalPotentialDownDown = off-diagonal contribution of the one-body potential for particles with spin down, the first entry is the annihilation index, 
  //                                       the second entry is the momentum tranfer
  // onebodyPotentialPairing =  one-body pairing term (sorted from component on the lowest Lz state to component on the highest Lz state), on site, symmetric spin up / spin down
  // onebodyOffDiagonalPotentialPairing = off diagonal contribution of the one-body pairing term, the first entry is the index of the rightmost creation operator, 
  //	    				  the second entry is the momentum tranfer
  // chargingEnergy = factor in front of the charging energy (i.e 1/(2C))
  // averageNumberParticles = average number of particles in the system
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta(ParticleOnSphereWithSpin* particles, int lzmax, double ratio, int maxMomentumTransfer, 
											double** pseudoPotential,
											double* onebodyPotentialUpUp, double* onebodyPotentialDownDown,
											Complex** onebodyOffDiagonalPotentialUpUp, Complex** onebodyOffDiagonalPotentialDownDown,
											Complex* onebodyPotentialPairing, Complex** onebodyOffDiagonalPotentialPairing,
											double chargingEnergy, double averageNumberParticles,
											AbstractArchitecture* architecture, long memory = -1);
  
  // destructor
  //
  ~ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta();

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // nbrPseudopotentials = number of pseudopotentials
  // pseudopotentials = pseudopotential coefficients
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials);

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

};

// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
												      int step, ComplexVector& vSource, ComplexVector& vDestination)
{
  if (this->OneBodyInteractionFactorsupup != 0)
    if (this->OneBodyInteractionFactorsdowndown != 0)
      {
	double TmpDiagonal = 0.0;
	for (int i = firstComponent; i < lastComponent; i += step)
	  { 
	    TmpDiagonal = 0.0;
	    for (int j = 0; j <= this->LzMax; ++j) 
	      {
		TmpDiagonal += this->OneBodyInteractionFactorsupup[j] * particles->AduAu(i, j);
		TmpDiagonal += this->OneBodyInteractionFactorsdowndown[j] * particles->AddAd(i, j);
	      }
	    vDestination[i] += (this->HamiltonianShift + TmpDiagonal)* vSource[i];
	  }
      }
    else
      {
	double TmpDiagonal = 0.0;
	for (int i = firstComponent; i < lastComponent; i += step)
	  { 
	    TmpDiagonal = 0.0;
	    for (int j = 0; j <= this->LzMax; ++j) 
	      TmpDiagonal += this->OneBodyInteractionFactorsupup[j] * particles->AduAu(i, j);
	    vDestination[i] += (this->HamiltonianShift + TmpDiagonal)* vSource[i];
	  }
      }
  else
    if (this->OneBodyInteractionFactorsdowndown != 0)
      {
	double TmpDiagonal = 0.0;
	for (int i = firstComponent; i < lastComponent; i += step)
	  { 
	    TmpDiagonal = 0.0;
	    for (int j = 0; j <= this->LzMax; ++j) 
	      TmpDiagonal += this->OneBodyInteractionFactorsdowndown[j] * particles->AddAd(i, j);
	    vDestination[i] += (this->HamiltonianShift + TmpDiagonal)* vSource[i];
	  }
      }	
    else
      for (int i = firstComponent; i < lastComponent; i += step)
	vDestination[i] += this->HamiltonianShift * vSource[i];
  if (this->OneBodyInteractionFactorsupdown != 0)
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
		  Index = particles->AddAu(i + this->PrecalculationShift, j, j, Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * Conj(OneBodyInteractionFactorsupdown[j])) * Source;
		    }
		  Index = particles->AduAd(i + this->PrecalculationShift, j, j, Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * OneBodyInteractionFactorsupdown[j]) * Source;
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
		   Index = particles->AddAu(i + this->PrecalculationShift, j, j, Coefficient);
		   if (Index <= i)
		     {
		       if (Index < i)
			 {
			   vDestination[Index] += (Coefficient * Conj(OneBodyInteractionFactorsupdown[j])) * Source;
			   vDestination[i] += (Coefficient * OneBodyInteractionFactorsupdown[j]) * vSource[Index];
			 }
		       else
			 {
			   vDestination[Index] += (Coefficient * Conj(OneBodyInteractionFactorsupdown[j])) * Source;
			 }
		     }
		   Index = particles->AduAd(i + this->PrecalculationShift, j, j, Coefficient);
		   if (Index < Dim)
		     {
		       if (Index <= i)
			 {
			   if (Index < i)
			     {
			       vDestination[Index] += (Coefficient * OneBodyInteractionFactorsupdown[j]) * Source;
			       vDestination[i] += (Coefficient * Conj(OneBodyInteractionFactorsupdown[j])) * vSource[Index];
			     }
			   else
			     {
			       vDestination[Index] += (Coefficient * OneBodyInteractionFactorsupdown[j]) * Source;
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

inline void ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
																	 int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  if (this->OneBodyInteractionFactorsupup != 0) 
    {
      if (this->OneBodyInteractionFactorsdowndown != 0)
	{
	  double TmpDiagonal = 0.0;
	  for (int p = 0; p < nbrVectors; ++p)
	    {
	      ComplexVector& TmpSourceVector = vSources[p];
	      ComplexVector& TmpDestinationVector = vDestinations[p];
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  TmpDiagonal = 0.0;
		  for (int j = 0; j <= this->LzMax; ++j) 
		    {
		      TmpDiagonal += this->OneBodyInteractionFactorsupup[j] * particles->AduAu(i, j);
		      TmpDiagonal += this->OneBodyInteractionFactorsdowndown[j] * particles->AddAd(i, j);
		    }
		  TmpDestinationVector[i] += (this->HamiltonianShift + TmpDiagonal)* TmpSourceVector[i];
		}
	    }
	}
      else
	{
	  double TmpDiagonal = 0.0;
	  for (int p = 0; p < nbrVectors; ++p)
	    {
	      ComplexVector& TmpSourceVector = vSources[p];
	      ComplexVector& TmpDestinationVector = vDestinations[p];
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  TmpDiagonal = 0.0;
		  for (int j = 0; j <= this->LzMax; ++j) 
		    TmpDiagonal += this->OneBodyInteractionFactorsupup[j] * particles->AduAu(i, j);
		TmpDestinationVector[i] += (this->HamiltonianShift + TmpDiagonal)* TmpSourceVector[i];
		}
	    }
	}
    }
  else
    {
      if (this->OneBodyInteractionFactorsdowndown != 0)
	{
	  double TmpDiagonal = 0.0;
	  for (int p = 0; p < nbrVectors; ++p)
	    {
	      ComplexVector& TmpSourceVector = vSources[p];
	      ComplexVector& TmpDestinationVector = vDestinations[p];
	      for (int i = firstComponent; i < lastComponent; i += step)
		{ 
		  TmpDiagonal = 0.0;
		  for (int j = 0; j <= this->LzMax; ++j) 
		    TmpDiagonal += this->OneBodyInteractionFactorsdowndown[j] * particles->AddAd(i, j);
		  TmpDestinationVector[i] += (this->HamiltonianShift + TmpDiagonal)* TmpSourceVector[i];
		}
	    }
	}	
      else
	for (int p = 0; p < nbrVectors; ++p)
	  {
	    ComplexVector& TmpSourceVector = vSources[p];
	    ComplexVector& TmpDestinationVector = vDestinations[p];
	    for (int i = firstComponent; i < lastComponent; i += step)
	      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	  }
    }
  for (int p = 0; p < nbrVectors; ++p)
    {
      ComplexVector& TmpSourceVector = vSources[p];
      ComplexVector& TmpDestinationVector = vDestinations[p];
      for (int i = firstComponent; i < lastComponent; i += step)
	TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
    }
  if (this->OneBodyInteractionFactorsupdown != 0)
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
		  Index = particles->AddAu(i, j, j, Coefficient);
		  if (Index < Dim)
		    {
		      for (int p = 0; p < nbrVectors; ++p)
			vDestinations[p][Index] += Coefficient * Conj(OneBodyInteractionFactorsupdown[j]) * vSources[p][i];
		    }
		  Index = particles->AduAd(i, j, j, Coefficient);
		  if (Index < Dim)
		    {
		      for (int p = 0; p < nbrVectors; ++p)
			vDestinations[p][Index] += Coefficient * OneBodyInteractionFactorsupdown[j] * vSources[p][i];
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
		  Index = particles->AddAu(i, j, j, Coefficient);
		  if (Index <= i)
		    {
		      if (Index < i)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * Conj(OneBodyInteractionFactorsupdown[j]) * vSources[p][i];
			      vDestinations[p][i] += Coefficient * OneBodyInteractionFactorsupdown[j] * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * Conj(OneBodyInteractionFactorsupdown[j]) * vSources[p][i];
			    }
			}
		    }
		  Index = particles->AduAd(i, j, j, Coefficient);
		  if (Index <= i)
		    {
		      if (Index < i)
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    {
			      vDestinations[p][Index] += Coefficient * OneBodyInteractionFactorsupdown[j] * vSources[p][i];
			      vDestinations[p][i] += Coefficient * Conj(OneBodyInteractionFactorsupdown[j]) * vSources[p][Index];
			    }
			}
		      else
			{
			  for (int p = 0; p < nbrVectors; ++p)
			    vDestinations[p][Index] += Coefficient * OneBodyInteractionFactorsupdown[j] * vSources[p][i];
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

inline void ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
																		int* indexArray, Complex* coefficientArray, long& position)
{
  if ((this->OneBodyInteractionFactorsdowndown != 0) || (this->OneBodyInteractionFactorsupup != 0))
    {
      double TmpDiagonal = 0.0;
      if (this->OneBodyInteractionFactorsupup != 0)
	for (int j = 0; j <= this->LzMax; ++j)
	  TmpDiagonal += this->OneBodyInteractionFactorsupup[j] * particles->AduAu(index + this->PrecalculationShift, j);
      if (this->OneBodyInteractionFactorsdowndown != 0)
	for (int j = 0; j <= this->LzMax; ++j)
	  TmpDiagonal += this->OneBodyInteractionFactorsdowndown[j] * particles->AddAd(index + this->PrecalculationShift, j);	  
      indexArray[position] = index + this->PrecalculationShift;
      if (this->HermitianSymmetryFlag == true)
	TmpDiagonal *= 0.5;
      coefficientArray[position] = TmpDiagonal;
      ++position;
    }
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  int Index;
  if (this->OneBodyInteractionFactorsupdown != 0)
    {
      for (int j = 0; j <= this->LzMax; ++j)
	{
	  Index = particles->AddAu(index + this->PrecalculationShift, j, j, Coefficient);
	  if (Index < Dim)
	    {
	      indexArray[position] = Index;
	      coefficientArray[position] = Coefficient * this->OneBodyInteractionFactorsupdown[j];
	      ++position;
	    }
	  Index = particles->AduAd(index + this->PrecalculationShift, j, j, Coefficient);
	  if (Index < Dim)
	    {
	      indexArray[position] = Index;
	      coefficientArray[position] = Coefficient * Conj(this->OneBodyInteractionFactorsupdown[j]);
	      ++position;
	    }
	}
    }
  if (this->OneBodyInteractionFactorsPairing != 0)
    {
      for (int j = 0; j <= this->LzMax; ++j)
	{
	  Index = particles->AuAd(index + this->PrecalculationShift, j, this->LzMax - j, Coefficient);
	  if (Index < Dim)
	    {
	      indexArray[position] = Index;
	      coefficientArray[position] = Coefficient * this->OneBodyInteractionFactorsPairing[j];
	      ++position;
	    }
	  Index = particles->AduAdd(index + this->PrecalculationShift, j, this->LzMax - j, Coefficient);
	  if (Index < Dim)
	    {
	      indexArray[position] = Index;
	      coefficientArray[position] = -Coefficient * Conj(this->OneBodyInteractionFactorsPairing[j]);
	      ++position;
	    }
	}
    }
  if (this->OneBodyOffDiagonalInteractionFactorsPairing != 0)
    {
      for (int k = 1; k <= this->MaximumMomentumTransfer; ++k)
	{
	  for (int j = 0; j <= (this->LzMax - k); ++j)
	    {
	      Index = particles->AuAd(index + this->PrecalculationShift, j, this->LzMax - (j + k), Coefficient);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = Coefficient * this->OneBodyOffDiagonalInteractionFactorsPairing[j][this->MaximumMomentumTransfer + k];
		  ++position;
		}
	      Index = particles->AduAdd(index + this->PrecalculationShift, j, this->LzMax - (j + k), Coefficient);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = -Coefficient * Conj(this->OneBodyOffDiagonalInteractionFactorsPairing[j][this->MaximumMomentumTransfer + k]);
		  ++position;
		}
	    }
	}
      for (int k = 1; k <= this->MaximumMomentumTransfer; ++k)
	{
	  for (int j = k; j <= this->LzMax; ++j)
	    {
	      Index = particles->AuAd(index + this->PrecalculationShift, j, this->LzMax - (j - k), Coefficient);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = Coefficient * this->OneBodyOffDiagonalInteractionFactorsPairing[j][this->MaximumMomentumTransfer - k];
		  ++position;
		}
	      Index = particles->AduAdd(index + this->PrecalculationShift, j, this->LzMax - (j - k), Coefficient);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = -Coefficient * Conj(this->OneBodyOffDiagonalInteractionFactorsPairing[j][this->MaximumMomentumTransfer - k]);
		  ++position;
		}
	    }
	}
    }
  if (this->OneBodyOffDiagonalInteractionFactorsupup != 0)
    {
      for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
	{
	  for (int j = 0; j < (this->LzMax - k); ++j)
	    {
	      Index = particles->AduAu(index + this->PrecalculationShift, j + k + 1, j, Coefficient);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = Coefficient * this->OneBodyOffDiagonalInteractionFactorsupup[j][k];
		  ++position;
		}
	      Index = particles->AduAu(index + this->PrecalculationShift, j, j + k + 1, Coefficient);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = Coefficient * Conj(this->OneBodyOffDiagonalInteractionFactorsupup[j][k]);
		  ++position;
		}
	    }
	}
    }
  if (this->OneBodyOffDiagonalInteractionFactorsdowndown != 0)
    {
      for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
	{
	  for (int j = 0; j < (this->LzMax - k); ++j)
	    {
	      Index = particles->AddAd(index + this->PrecalculationShift, j + k + 1, j, Coefficient);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = Coefficient * this->OneBodyOffDiagonalInteractionFactorsdowndown[j][k];
		  ++position;
		}
	      Index = particles->AddAd(index + this->PrecalculationShift, j, j + k + 1, Coefficient);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = Coefficient * Conj(this->OneBodyOffDiagonalInteractionFactorsdowndown[j][k]);
		  ++position;
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

inline void ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;

  int Dim = particles->GetHilbertSpaceDimension();
  
  if (this->HermitianSymmetryFlag == false)
    {
      if (this->OneBodyInteractionFactorsupdown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j = 0; j <= this->LzMax; ++j)
		{
		  Index = particles->AddAu(i, j, j, Coefficient);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		  Index = particles->AduAd(i, j, j, Coefficient);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		}
	    }
	}
      if (this->OneBodyInteractionFactorsPairing != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j = 0; j <= this->LzMax; ++j)
		{
		  Index = particles->AuAd(i, j, this->LzMax - j, Coefficient);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		  Index = particles->AduAdd(i, j, this->LzMax - j, Coefficient);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		}
	    }
	}
      if (this->OneBodyOffDiagonalInteractionFactorsPairing != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int k = 1; k <= this->MaximumMomentumTransfer; ++k)
		{
		  for (int j = 0; j <= (this->LzMax - k); ++j)
		    {
		      Index = particles->AuAd(i, j, this->LzMax - (j + k), Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		      Index = particles->AduAdd(i, j, this->LzMax - (j + k), Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		}
	      for (int k = 1; k <= this->MaximumMomentumTransfer; ++k)
		{
		  for (int j = k; j <= this->LzMax; ++j)
		    {
		      Index = particles->AuAd(i, j, this->LzMax - (j - k), Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		      Index = particles->AduAdd(i, j, this->LzMax - (j - k), Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		}
	    }
	}
      if (this->OneBodyOffDiagonalInteractionFactorsupup != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
		{
		  for (int j = 0; j < (this->LzMax - k); ++j)
		    {
		      Index = particles->AduAu(i, j + k + 1, j, Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		      Index = particles->AduAu(i, j, j + k + 1, Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		}
	    }
	}
      if (this->OneBodyOffDiagonalInteractionFactorsdowndown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
		{
		  for (int j = 0; j < (this->LzMax - k); ++j)
		    {
		      Index = particles->AddAd(i, j + k + 1, j, Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		      Index = particles->AddAd(i, j, j + k + 1, Coefficient);
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
      if (this->OneBodyInteractionFactorsupdown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j <= this->LzMax; ++j)
		{
		  Index = particles->AddAu(i, j, j, Coefficient);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		  Index = particles->AduAd(i, j, j, Coefficient);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		}
	    }
	}
      if (this->OneBodyInteractionFactorsPairing != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j = 0; j <= this->LzMax; ++j)
		{
		  Index = particles->AuAd(i, j, this->LzMax - j, Coefficient);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		  Index = particles->AduAdd(i, j, this->LzMax - j, Coefficient);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		}
	    }
	}
      if (this->OneBodyOffDiagonalInteractionFactorsPairing != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int k = 1; k <= this->MaximumMomentumTransfer; ++k)
		{
		  for (int j = 0; j <= (this->LzMax - k); ++j)
		    {
		      Index = particles->AuAd(i, j, this->LzMax - (j + k), Coefficient);
		      if (Index <= i)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		      Index = particles->AduAdd(i, j, this->LzMax - (j + k), Coefficient);
		      if (Index <= i)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		}
	      for (int k = 1; k <= this->MaximumMomentumTransfer; ++k)
		{
		  for (int j = k; j <= this->LzMax; ++j)
		    {
		      Index = particles->AuAd(i, j, this->LzMax - (j - k), Coefficient);
		      if (Index <= i)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		      Index = particles->AduAdd(i, j, this->LzMax - (j - k), Coefficient);
		      if (Index <= i)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		}
	    }
	}
      if (this->OneBodyOffDiagonalInteractionFactorsupup != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
		{
		  for (int j = 0; j < (this->LzMax - k); ++j)
		    {
		      Index = particles->AduAu(i, j + k + 1, j, Coefficient);
		      if (Index <= i)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		      Index = particles->AduAu(i, j, j + k + 1, Coefficient);
		      if (Index <= i)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		}
	    }
	}
      if (this->OneBodyOffDiagonalInteractionFactorsdowndown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int k = 0; k < this->MaximumMomentumTransfer; ++k)
		{
		  for (int j = 0; j < (this->LzMax - k); ++j)
		    {
		      Index = particles->AddAd(i, j + k + 1, j, Coefficient);
		      if (Index <= i)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		      Index = particles->AddAd(i, j, j + k + 1, Coefficient);
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

  if ((this->OneBodyInteractionFactorsdowndown != 0) || (this->OneBodyInteractionFactorsupup != 0))
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  ++memory;
	  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
	}
    }
}


#endif
