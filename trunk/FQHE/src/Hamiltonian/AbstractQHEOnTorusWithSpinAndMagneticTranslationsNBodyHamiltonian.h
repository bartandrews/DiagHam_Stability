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
//                                                                            //
//                        last modification : 23/10/2015                      //
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


#ifndef ABSTRACTQHEONTORUSWITHSPINANDMAGNETICTRANSLATIONSNBODYHAMILTONIAN_H
#define ABSTRACTQHEONTORUSWITHSPINANDMAGNETICTRANSLATIONSNBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian.h"


#include <iostream>


using std::ostream;


class AbstractArchitecture;


class AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian : public ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian
{

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
  virtual ~AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian();
  
 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors() = 0;

  // evaluate all exponential factors
  //   
  virtual void EvaluateExponentialFactors();

  // get all the indices that should appear in the annihilation/creation operators
  //
  // onlyIntraFlag = true if only the intra component indices have to be computed
  virtual void GetIndices(bool onlyIntraFlag = false);

  // get all indices needed to characterize a completly skew symmetric tensor, ordered by the sum of the indices
  //
  // nbrValues = number of different values an index can have
  // nbrIndices = number of indices 
  // nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
  // sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
  // return value = total number of index groups
  virtual long GetAllSkewSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum);

  // get all indices needed to characterize a  tensor made of two completly skew symmetric  sets of indices, sorted by the sum of the indices
  //
  // nbrValues = number of different values an index can have
  // nbrIndicesUp = number of indices for the first set of indices (i.e. spin up)
  // nbrIndicesDown = number of indices for the first set of indices (i.e. spin down), warning nbrIndicesDown should lower of equal to nbrIndicesUp
  // nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
  // sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
  // return value = total number of index groups
  virtual long GetAllTwoSetSkewSymmetricIndices (int nbrValues, int nbrIndicesUp, int nbrIndicesDown, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum);

  // get all indices needed to characterize a completly symmetric tensor, sorted by the sum of the indices
  //
  // nbrValues = number of different values an index can have
  // nbrIndices = number of indices 
  // nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
  // sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
  // sortedIndicesPerSumSymmetryFactor = reference on a array where symmetry factor (aka inverse of the product of the factorial of the number 
  //                                      of time each index appears) are stored (first array dimension corresponding to sum of the indices)
  // return value = total number of index groups
  virtual long GetAllSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum,
				       double**& sortedIndicesPerSumSymmetryFactor);

  // get all indices needed to characterize a  tensor made of two completly symmetric  sets of indices, sorted by the sum of the indices
  //
  // nbrValues = number of different values an index can have
  // nbrIndicesUp = number of indices for the first set of indices (i.e. spin up)
  // nbrIndicesDown = number of indices for the first set of indices (i.e. spin down), warning nbrIndicesDown should lower of equal to nbrIndicesUp
  // nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
  // sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
  // sortedIndicesPerSumSymmetryFactor = reference on a array where symmetry factor (aka inverse of the product of the factorial of the number 
  //                                      of time each index appears) are stored (first array dimension corresponding to sum of the indices)
  // return value = total number of index groups
  virtual long GetAllTwoSetSymmetricIndices (int nbrValues, int nbrIndicesUp, int nbrIndicesDown, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum,
					     double**& sortedIndicesPerSumSymmetryFactor);

  // get all indices (without assuming any symmetry), sorted by the sum of the indices
  //
  // nbrValues = number of different values an index can have
  // nbrIndices = number of indices 
  // nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
  // sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
  // return value = total number of index groups
  virtual long GetAllIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum);

  // core part of the AddMultiply method involving n-body term
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // nbodyIndex = value of n
  void EvaluateMNNBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination, int nbodyIndex);

  // core part of the AddMultiply method involving n-body term for a set ofe vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // nbodyIndex = value of n
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void EvaluateMNNBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
						   ComplexVector* vDestinations, int nbrVectors, int nbodyIndex, Complex* tmpCoefficients);

  // core part of the AddMultiply method involving the n-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void HermitianEvaluateMNNBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination, int nbodyIndex);

  // core part of the AddMultiply method involving the n-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void HermitianEvaluateMNNBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
							    ComplexVector* vDestinations, int nbrVectors, int nbodyIndex,
							    Complex* tmpCoefficients, Complex* tmpSum);

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

  // core part of the FastMultiplication method involving n-body term
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // nbodyIndex = value of n
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNNBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, int nbodyIndex, 
							  int* indexArray, Complex* coefficientArray, long& position);

  // core part of the PartialFastMultiplicationMemory method involving n-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // nbodyIndex = value of n
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNNBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, int nbodyIndex, long& memory);

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

  // core part of the FastMultiplication method involving the one-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
							    int* indexArray, Complex* coefficientArray, long& position);

  // core part of the PartialFastMultiplicationMemory method involving two-body term and one-body terms
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory);

};

// core part of the AddMultiply method involving n-body term
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// nbodyIndex = value of n

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::EvaluateMNNBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, 
														   ComplexVector& vDestination, int nbodyIndex)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int NbrTranslations;

  int SqrNBodyValue = nbodyIndex*nbodyIndex;

  for (int spinSector = 0; spinSector < NbrSpinSectors[nbodyIndex]; ++spinSector)
    {
      int *TmpSpinIndicesA = SpinIndices[nbodyIndex][spinSector];
      int *TmpSpinIndicesAd = SpinIndices[nbodyIndex][spinSector]+nbodyIndex;
      int* TmpIndices;
      Complex* TmpInteractionFactor;
      int Index;
      for (int j = 0; j < this->NbrNBodySpinMomentumSectorSum[nbodyIndex][spinSector]; ++j)
	{
	  int Lim = nbodyIndex * this->NbrNBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
	  TmpIndices = this->NBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
	  for (int i1 = 0; i1 < Lim; i1 += nbodyIndex)
	    {
	      Coefficient3 = particles->ProdA(index, TmpIndices + i1, TmpSpinIndicesA, nbodyIndex);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->NBodyInteractionFactors[nbodyIndex][spinSector][j][(i1 * Lim) / SqrNBodyValue]);
		  Coefficient4 = vSource[index];
		  Coefficient4 *= Coefficient3;
		  for (int i2 = 0; i2 < Lim; i2 += nbodyIndex)
		    {
		      Index = particles->ProdAd(TmpIndices + i2, TmpSpinIndicesAd, nbodyIndex, Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * Coefficient4;
			}
		      ++TmpInteractionFactor;
		    }
		}
	    }
	}
    }
}



// core part of the AddMultiply method involving n-body term for a set ofe vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// nbodyIndex = value of n
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::EvaluateMNNBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, 
														   ComplexVector* vSources, ComplexVector* vDestinations, 
														   int nbrVectors, int nbodyIndex, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int NbrTranslations;

  int SqrNBodyValue = nbodyIndex*nbodyIndex;

  for (int spinSector = 0; spinSector < NbrSpinSectors[nbodyIndex]; ++spinSector)
    {
      int *TmpSpinIndicesA = SpinIndices[nbodyIndex][spinSector];
      int *TmpSpinIndicesAd = SpinIndices[nbodyIndex][spinSector]+nbodyIndex;
      int* TmpIndices;
      Complex* TmpInteractionFactor;
      int Index;
      for (int j = 0; j < this->NbrNBodySpinMomentumSectorSum[nbodyIndex][spinSector]; ++j)
	{
	  int Lim = nbodyIndex * this->NbrNBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
	  TmpIndices = this->NBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
	  for (int i1 = 0; i1 < Lim; i1 += nbodyIndex)
	    {
	      Coefficient3 = particles->ProdA(index, TmpIndices + i1, TmpSpinIndicesA, nbodyIndex);
	      if (Coefficient3 != 0.0)
		{
		  for (int l = 0; l < nbrVectors; ++l)
		    tmpCoefficients[l] = Coefficient3 * vSources[l][index];
		  TmpInteractionFactor = &(this->NBodyInteractionFactors[nbodyIndex][spinSector][j][(i1 * Lim) / SqrNBodyValue]);
		  for (int i2 = 0; i2 < Lim; i2 += nbodyIndex)
		    {
		      Index = particles->ProdAd(TmpIndices + i2, TmpSpinIndicesAd, nbodyIndex, Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
			  Coefficient4 = Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor);
			  for (int l = 0; l < nbrVectors; ++l)
			    vDestinations[l][Index] += Coefficient4 * tmpCoefficients[l];
			}
		      ++TmpInteractionFactor;
		    }
		}
	    }
	}
    }
}


// core part of the AddMultiply method involving the n-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::HermitianEvaluateMNNBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, 
															    ComplexVector& vSource, ComplexVector& vDestination, 
															    int nbodyIndex)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex TmpSum=0.0;
  Complex Coefficient4;
  int NbrTranslations;

  int SqrNBodyValue = nbodyIndex*nbodyIndex;

  for (int spinSector = 0; spinSector < NbrSpinSectors[nbodyIndex]; ++spinSector)
    {
      int *TmpSpinIndicesA = SpinIndices[nbodyIndex][spinSector];
      int *TmpSpinIndicesAd = SpinIndices[nbodyIndex][spinSector]+nbodyIndex;
      int* TmpIndices;
      Complex* TmpInteractionFactor;
      int Index;
      for (int j = 0; j < this->NbrNBodySpinMomentumSectorSum[nbodyIndex][spinSector]; ++j)
	{
	  int Lim = nbodyIndex * this->NbrNBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
	  TmpIndices = this->NBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
	  for (int i1 = 0; i1 < Lim; i1 += nbodyIndex)
	    {
	      Coefficient3 = particles->ProdA(index, TmpIndices + i1, TmpSpinIndicesA, nbodyIndex);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->NBodyInteractionFactors[nbodyIndex][spinSector][j][(i1 * Lim) / SqrNBodyValue]);
		  Coefficient4 = vSource[index];
		  Coefficient4 *= Coefficient3;
		  for (int i2 = 0; i2 < Lim; i2 += nbodyIndex)
		    {
		      Index = particles->ProdAd(TmpIndices + i2, TmpSpinIndicesAd, nbodyIndex, Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  if (Index < index)
			    TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			  vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * Coefficient4;
			}
		      ++TmpInteractionFactor;
		    }
		}
	    }
	}
    }
  vDestination[index] += TmpSum;
}


// core part of the AddMultiply method involving the n-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::HermitianEvaluateMNNBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, 
															    ComplexVector* vSources, ComplexVector* vDestinations, 
															    int nbrVectors, int nbodyIndex, Complex* tmpCoefficients, 
															    Complex* tmpSum)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int NbrTranslations;

  int SqrNBodyValue = nbodyIndex*nbodyIndex;
  
  for (int i=0; i<nbrVectors; ++i)
    tmpSum[i]=0.0;

  for (int spinSector = 0; spinSector < NbrSpinSectors[nbodyIndex]; ++spinSector)
    {
      int *TmpSpinIndicesA = SpinIndices[nbodyIndex][spinSector];
      int *TmpSpinIndicesAd = SpinIndices[nbodyIndex][spinSector]+nbodyIndex;
      int* TmpIndices;
      Complex* TmpInteractionFactor;
      int Index;
      for (int j = 0; j < this->NbrNBodySpinMomentumSectorSum[nbodyIndex][spinSector]; ++j)
	{
	  int Lim = nbodyIndex * this->NbrNBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
	  TmpIndices = this->NBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
	  for (int i1 = 0; i1 < Lim; i1 += nbodyIndex)
	    {
	      Coefficient3 = particles->ProdA(index, TmpIndices + i1, TmpSpinIndicesA, nbodyIndex);
	      if (Coefficient3 != 0.0)
		{
		  for (int l = 0; l < nbrVectors; ++l)
		    tmpCoefficients[l] = Coefficient3 * vSources[l][index];
		  TmpInteractionFactor = &(this->NBodyInteractionFactors[nbodyIndex][spinSector][j][(i1 * Lim) / SqrNBodyValue]);
		  for (int i2 = 0; i2 < Lim; i2 += nbodyIndex)
		    {
		      Index = particles->ProdAd(TmpIndices + i2, TmpSpinIndicesAd, nbodyIndex, Coefficient, NbrTranslations);
		      if (Index <= index)
			{
			  Coefficient4 = this->ExponentialFactors[NbrTranslations] * Coefficient*(*TmpInteractionFactor);
			  for (int l = 0; l < nbrVectors; ++l)
			    {
			      vDestinations[l][Index] += Coefficient4 * tmpCoefficients[l];
			      if (Index < index)
				tmpSum[l] += (Coefficient * Coefficient3) * Conj(this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * vSources[l][Index];
			    }
			}
		      ++TmpInteractionFactor;
		    }
		}
	    }
	}
    }
  for (int l = 0; l < nbrVectors; ++l)
    vDestinations[l][index] += tmpSum[l];
}


// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, 
														     ComplexVector& vSource, ComplexVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
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
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		    }
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    vDestination[Index] += (Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor))) * Coefficient4;
		  ++TmpInteractionFactor;
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

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
												      ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  int NbrTranslations;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
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
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 2)
		{
		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
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

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  //int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
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
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
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
		  ++TmpInteractionFactor;
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

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
													       ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  //int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  int NbrTranslations;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
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
		  ++TmpInteractionFactor;
		}
	    }
	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
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
		  ++TmpInteractionFactor;
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

// core part of the FastMultiplication method involving n-body term
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// nbodyIndex = value of n
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::EvaluateMNNBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, int nbodyIndex, 
															  int* indexArray, Complex* coefficientArray, long& position)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int NbrTranslations;

  int SqrNBodyValue = nbodyIndex*nbodyIndex;

  if (this->HermitianSymmetryFlag == false)
    {
      for (int spinSector = 0; spinSector < NbrSpinSectors[nbodyIndex]; ++spinSector)
	{
	  int *TmpSpinIndicesA = SpinIndices[nbodyIndex][spinSector];
	  int *TmpSpinIndicesAd = SpinIndices[nbodyIndex][spinSector]+nbodyIndex;
	  int* TmpIndices;
	  Complex* TmpInteractionFactor;
	  int Index;
	  for (int j = 0; j < this->NbrNBodySpinMomentumSectorSum[nbodyIndex][spinSector]; ++j)
	    {
	      int Lim = nbodyIndex * this->NbrNBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
	      TmpIndices = this->NBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
	      for (int i1 = 0; i1 < Lim; i1 += nbodyIndex)
		{
		  Coefficient3 = particles->ProdA(index, TmpIndices + i1, TmpSpinIndicesA, nbodyIndex);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->NBodyInteractionFactors[nbodyIndex][spinSector][j][(i1 * Lim) / SqrNBodyValue]);
		      for (int i2 = 0; i2 < Lim; i2 += nbodyIndex)
			{
			  Index = particles->ProdAd(TmpIndices + i2, TmpSpinIndicesAd, nbodyIndex, Coefficient, NbrTranslations);
			  if (Index < Dim)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = (Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * Coefficient3;
			      ++position;
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
      int AbsoluteIndex = index + this->PrecalculationShift;
      for (int spinSector = 0; spinSector < NbrSpinSectors[nbodyIndex]; ++spinSector)
	{
	  int *TmpSpinIndicesA = SpinIndices[nbodyIndex][spinSector];
	  int *TmpSpinIndicesAd = SpinIndices[nbodyIndex][spinSector]+nbodyIndex;
	  int* TmpIndices;
	  Complex* TmpInteractionFactor;
	  int Index;
	  for (int j = 0; j < this->NbrNBodySpinMomentumSectorSum[nbodyIndex][spinSector]; ++j)
	    {
	      int Lim = nbodyIndex * this->NbrNBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
	      TmpIndices = this->NBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
	      for (int i1 = 0; i1 < Lim; i1 += nbodyIndex)
		{
		  Coefficient3 = particles->ProdA(index, TmpIndices + i1, TmpSpinIndicesA, nbodyIndex);
		  if (Coefficient3 != 0.0)
		    {
		      TmpInteractionFactor = &(this->NBodyInteractionFactors[nbodyIndex][spinSector][j][(i1 * Lim) / SqrNBodyValue]);
		      for (int i2 = 0; i2 < Lim; i2 += nbodyIndex)
			{
			  Index = particles->ProdAd(TmpIndices + i2, TmpSpinIndicesAd, nbodyIndex, Coefficient, NbrTranslations);
			  if (Index <= AbsoluteIndex)
			    {
			      if (Index == AbsoluteIndex)
				{
				  indexArray[position] = Index;
				  coefficientArray[position] = 0.5 * (Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * Coefficient3;
				  ++position;
				}
			      else
				{
				  indexArray[position] = Index;
				  coefficientArray[position] = (Coefficient * this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor)) * Coefficient3;
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
}

// core part of the PartialFastMultiplicationMemory method involving n-body term
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// nbodyIndex = value of n
// memory = reference on the amount of memory required for precalculations

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::EvaluateMNNBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, 
																int lastComponent, int nbodyIndex, long& memory)
{
  int Dim = particles->GetHilbertSpaceDimension();
  int Index;
  double Coefficient = 0.0;
  int* TmpIndices;
  double Coefficient3;
  int NbrTranslations;
  if (this->HermitianSymmetryFlag == false)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int spinSector = 0; spinSector < NbrSpinSectors[nbodyIndex]; ++spinSector)
	    {
	      int* TmpSpinIndicesA = SpinIndices[nbodyIndex][spinSector];
	      int* TmpSpinIndicesAd = SpinIndices[nbodyIndex][spinSector] + nbodyIndex;
	      for (int j = 0; j < this->NbrNBodySpinMomentumSectorSum[nbodyIndex][spinSector]; ++j)
		{
		  int Lim = nbodyIndex * this->NbrNBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
		  TmpIndices = this->NBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
		  for (int i1 = 0; i1 < Lim; i1 += nbodyIndex)
		    {
		      Coefficient3 = particles->ProdA(i, TmpIndices + i1, TmpSpinIndicesA, nbodyIndex);
		      if (Coefficient3 != 0.0)
			{
			  for (int i2 = 0; i2 < Lim; i2 += nbodyIndex)
			    {
			      Index = particles->ProdAd(TmpIndices + i2, TmpSpinIndicesAd, nbodyIndex, Coefficient, NbrTranslations);
			      if (Index < Dim)
				{
				  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
				  ++memory;
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
	  for (int spinSector = 0; spinSector < NbrSpinSectors[nbodyIndex]; ++spinSector)
	    {
	      int* TmpSpinIndicesA = SpinIndices[nbodyIndex][spinSector];
	      int* TmpSpinIndicesAd = SpinIndices[nbodyIndex][spinSector] + nbodyIndex;
	      for (int j = 0; j < this->NbrNBodySpinMomentumSectorSum[nbodyIndex][spinSector]; ++j)
		{
		  int Lim = nbodyIndex * this->NbrNBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
		  TmpIndices = this->NBodySpinMomentumSectorIndicesPerSum[nbodyIndex][spinSector][j];
		  for (int i1 = 0; i1 < Lim; i1 += nbodyIndex)
		    {
		      Coefficient3 = particles->ProdA(i, TmpIndices + i1, TmpSpinIndicesA, nbodyIndex);
		      if (Coefficient3 != 0.0)
			{
			  for (int i2 = 0; i2 < Lim; i2 += nbodyIndex)
			    {
			      Index = particles->ProdAd(TmpIndices + i2, TmpSpinIndicesAd, nbodyIndex, Coefficient, NbrTranslations);
			      if (Index <= i)
				{
				  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
				  ++memory;
				}
			    }
			}
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

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
															 int* indexArray, Complex* coefficientArray, long& position)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
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
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
		      if (Index < Dim)
			{
                          indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (this->ExponentialFactors[NbrTranslations] * (*TmpInteractionFactor));
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient2 = particles->AdAd(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient, NbrTranslations);
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
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient2 = particles->AdAd(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
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
		      ++TmpInteractionFactor;
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

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
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

// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
														  int step, ComplexVector& vSource, ComplexVector& vDestination)
{
  if (this->OneBodyInteractionFactorsupup != 0)
    {
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
    }
  else
    {
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
    }
  if (this->OneBodyInteractionFactorsupdown != 0)
    {
      double Coefficient;
      Complex Source;
      int NbrTranslations;
      int Dim = particles->GetHilbertSpaceDimension();
      int Index;
      for (int i = firstComponent; i < lastComponent; i += step)
	{
	  Source = vSource[i];
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->AddAu(i + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(OneBodyInteractionFactorsupdown[j])) * Source;
		}
	      Index = particles->AduAd(i + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  vDestination[Index] += (Coefficient * this->ExponentialFactors[NbrTranslations] * OneBodyInteractionFactorsupdown[j]) * Source;
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

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
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
	{
	  for (int p = 0; p < nbrVectors; ++p)
	    {
	      ComplexVector& TmpSourceVector = vSources[p];
	      ComplexVector& TmpDestinationVector = vDestinations[p];
	      for (int i = firstComponent; i < lastComponent; i += step)
		TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	    }
	}
    }
  if (this->OneBodyInteractionFactorsupdown != 0)
    {
      double Coefficient;
      int Dim = particles->GetHilbertSpaceDimension();
      int Index;
      int NbrTranslations;
      for (int i = firstComponent; i < lastComponent; i += step)
	{
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->AddAu(i + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  for (int p = 0; p < nbrVectors; ++p)
		    vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(OneBodyInteractionFactorsupdown[j]) * vSources[p][i];
		}
	      Index = particles->AduAd(i + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  for (int p = 0; p < nbrVectors; ++p)
		    vDestinations[p][Index] += Coefficient * this->ExponentialFactors[NbrTranslations] * OneBodyInteractionFactorsupdown[j] * vSources[p][i];
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

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
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
  if (this->OneBodyInteractionFactorsupdown != 0)
    {
      if (this->HermitianSymmetryFlag == false)
	{
	  int NbrTranslations;
	  int Dim = particles->GetHilbertSpaceDimension();
	  double Coefficient;
	  int Index;
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->AddAu(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyInteractionFactorsupdown[j];
		  ++position;
		}
	      Index = particles->AduAd(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index < Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyInteractionFactorsupdown[j]);
		  ++position;
		}
	    }
	}
      else
	{
	  int NbrTranslations;
	  int Dim = particles->GetHilbertSpaceDimension();
	  double Coefficient;
	  int Index;
	  int AbsoluteIndex = index + this->PrecalculationShift;
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->AddAu(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index <= AbsoluteIndex)
		{
		  if (Index == AbsoluteIndex)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = 0.5 * Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyInteractionFactorsupdown[j];
		      ++position;
		    }
		  else
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * this->OneBodyInteractionFactorsupdown[j];
		      ++position;
		    }
		}
	      Index = particles->AduAd(index + this->PrecalculationShift, j, j, Coefficient, NbrTranslations);
	      if (Index <= AbsoluteIndex)
		{
		  if (Index == AbsoluteIndex)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = 0.5 * Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyInteractionFactorsupdown[j]);
		      ++position;
		    }
		  else
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * this->ExponentialFactors[NbrTranslations] * Conj(this->OneBodyInteractionFactorsupdown[j]);
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

inline void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  int NbrTranslations;
  int Dim = particles->GetHilbertSpaceDimension();
  
  if (this->HermitianSymmetryFlag == false)
    {
      if (this->OneBodyInteractionFactorsupdown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j<= this->LzMax; ++j)
		{
		  Index = particles->AddAu(i, j, j, Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		  Index = particles->AduAd(i, j, j, Coefficient, NbrTranslations);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
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
	      for (int j=0; j<= this->LzMax; ++j)
		{
		  Index = particles->AddAu(i, j, j, Coefficient, NbrTranslations);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		    }
		  Index = particles->AduAd(i, j, j, Coefficient, NbrTranslations);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
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
