////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of abstract quantum Hall hamiltonian associated            //
//     to particles on a sphere with spin and n-body interaction terms        //
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


#ifndef PARTICLEONLATTICEWITHSPINCHERNINSULATORNBODYHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINCHERNINSULATORNBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinChernInsulatorHamiltonian.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian : public ParticleOnLatticeWithSpinChernInsulatorHamiltonian
{

 protected:
  
  // maximum number of particles that can interact together
  int MaxNBody;
  // indicates which n-body interaction terms are present in the Hamiltonian (only entries beyond three-body are considered!)
  bool* NBodyFlags;
  // number of spin indices - array[NbrNBody] format (annihilation1,...,annihilationNBody)(creation1,...,creationNBody)
  int* NbrSpinSectors;
  // sign in front of each n-body interaction term (including the one coming from the statistics, per spin indices
  double** NBodySign;
  // spin indices for each annihilation/creation opertor set per n-body interaction
  int*** SpinIndices;
  // spin indices for each annihilation/creation opertor set per n-body interaction, encoded in a unique integer
  int** SpinIndicesShort;


  // number of index sum for the creation / annihilation operators for each NBody interaction and spin sector
  int** NbrNBodySpinMomentumSectorSum;
  // number of index groups per index sum for the creation / annihilation operators for each NBody interaction and spin sector
  int*** NbrNBodySpinMomentumSectorIndicesPerSum;
  // array containing the momentum indices per sum index
  int**** NBodySpinMomentumSectorIndicesPerSum;

  // array of the corresponding interaction factors for the above sectors
  Complex**** NBodyInteractionFactors;
  
  // flag to indicate a fully defined (i.e with pseudo-potentials) two body interaction
  bool FullTwoBodyFlag;

 public:

  // destructor
  //
  virtual ~ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian();


  // symmetrize interaction factors to enable hermitian matrix multiplication
  // return = true upon success
  virtual bool HermitianSymmetrizeInteractionFactors();


  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
  virtual  ComplexVector* LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						      int firstComponent, int nbrComponent);
  

 protected:

 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors() = 0;

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
  //  virtual void EvaluateMNNBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, int nbodyIndex, long& memory);

  // core part of the PartialFastMultiplicationMemory method involving n-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // nbodyIndex = value of n
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNNBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, int nbodyIndex, long& memory);

  // compute the permutation array for the interaction indices
  // 
  // permutations = arrays where permuted indices will be stored
  // permutationSign = array with the sign of each permutation (initialized only when dealing with fermionic statistics)
  // nbodyIndex = order of interaction to consider
  // return value = number of permutations
  virtual int ComputePermutations(int**& permutations, double*& permutationSign, int nbodyIndex);

};

// core part of the AddMultiply method involving n-body term
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// nbodyIndex = value of n

inline void ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::EvaluateMNNBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination, int nbodyIndex)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;

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
		      Index = particles->ProdAd(TmpIndices + i2, TmpSpinIndicesAd, nbodyIndex, Coefficient);
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



// core part of the AddMultiply method involving n-body term for a set ofe vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// nbodyIndex = value of n
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::EvaluateMNNBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
													ComplexVector* vDestinations, int nbrVectors, int nbodyIndex, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;

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
		      Index = particles->ProdAd(TmpIndices + i2, TmpSpinIndicesAd, nbodyIndex, Coefficient);
		      if (Index < Dim)
			{
			  Coefficient4 = Coefficient*(*TmpInteractionFactor);
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

inline void ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::HermitianEvaluateMNNBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination, int nbodyIndex)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex TmpSum=0.0;
  Complex Coefficient4;

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
		      Index = particles->ProdAd(TmpIndices + i2, TmpSpinIndicesAd, nbodyIndex, Coefficient);
		      if (Index < Dim)
			{
			  TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
			  vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
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

inline void ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::HermitianEvaluateMNNBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
														    ComplexVector* vDestinations, int nbrVectors, int nbodyIndex, Complex* tmpCoefficients, Complex* tmpSum)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;

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
		      Index = particles->ProdAd(TmpIndices + i2, TmpSpinIndicesAd, nbodyIndex, Coefficient);
		      if (Index < Dim)
			{
			  Coefficient4 = Coefficient*(*TmpInteractionFactor);
			  for (int l = 0; l < nbrVectors; ++l)
			    {
			      vDestinations[l][Index] += Coefficient4 * tmpCoefficients[l];
			      tmpSum[l] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor)) * vSources[l][Index];
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


// core part of the FastMultiplication method involving n-body term
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// nbodyIndex = value of n
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::EvaluateMNNBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, int nbodyIndex, 
														int* indexArray, Complex* coefficientArray, long& position)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;

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
		  for (int i2 = 0; i2 < Lim; i2 += nbodyIndex)
		    {
		      Index = particles->ProdAd(TmpIndices + i2, TmpSpinIndicesAd, nbodyIndex, Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = (Coefficient * (*TmpInteractionFactor)) * Coefficient3;
			  //cout << "(*TmpInteractionFactor["<<position<<"])="<<(*TmpInteractionFactor)<<endl;
			  ++position;
			}
		      ++TmpInteractionFactor;
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

inline void ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::EvaluateMNNBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, 
														     int nbodyIndex, long& memory)
{
  int Dim = particles->GetHilbertSpaceDimension();
  int Index;
  double Coefficient = 0.0;
  int* TmpIndices;
  double Coefficient3;
  for (int i = firstComponent; i < lastComponent; ++i)
    {

      for (int spinSector = 0; spinSector < NbrSpinSectors[nbodyIndex]; ++spinSector)
	{
	  int *TmpSpinIndicesA = SpinIndices[nbodyIndex][spinSector];
	  int *TmpSpinIndicesAd = SpinIndices[nbodyIndex][spinSector]+nbodyIndex;
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
			  Index = particles->ProdAd(TmpIndices + i2, TmpSpinIndicesAd, nbodyIndex, Coefficient);
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



#endif
