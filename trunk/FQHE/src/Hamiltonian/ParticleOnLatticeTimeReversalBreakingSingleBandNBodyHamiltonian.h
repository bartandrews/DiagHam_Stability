////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                class of hamiltonian with particles on lattice              //
//       with time reversal breaking in the single band approximation         //
//                           and n-body interaction                           //
//                                                                            //
//                        last modification : 03/08/2011                      //
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


#ifndef PARTICLEONLATTICETIMEREVERSALBREAKINGSINGLEBANDNBODYHAMILTONIAN_H
#define PARTICLEONLATTICETIMEREVERSALBREAKINGSINGLEBANDNBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnLatticeTimeReversalBreakingSingleBandNBodyHamiltonian : public ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian
{

 protected:
  
  // number of index sum for the creation (or annhilation) operators for the intra spin sector
  int NbrNBodySectorSums;
  // array containing the number of index group per index sum for the creation (or annhilation) operators for the intra spin sector
  int* NbrNBodySectorIndicesPerSum;
  // array containing the (m1,m2,m3) indices per index sum for the creation (or annhilation) operators for the intra spin sector
  int** NBodySectorIndicesPerSum;

  // arrays containing all interaction factors, the first index correspond correspond to index sum for the creation (or annhilation) operators
  // the second index is a linearized index (m1,m2,m3) + (n1,n2,n3) * (nbr element in current index sum) (m for creation operators, n for annhilation operators)
  Complex** NBodyInteractionFactors;

  // true if one/two body terms are set
  bool TwoBodyFlag;  

  // number of interacting particles
  int NBodyValue;
  // square value of NBodyValue
  int SqrNBodyValue;


 public:

  // default constructor
  //
  ParticleOnLatticeTimeReversalBreakingSingleBandNBodyHamiltonian();

  // destructor
  //
  ~ParticleOnLatticeTimeReversalBreakingSingleBandNBodyHamiltonian();
  
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
  virtual ComplexVector* LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						     int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
  virtual ComplexVector* HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							      int firstComponent, int nbrComponent);


 protected:
 
  // core part of the AddMultiply method involving the n-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  virtual void EvaluateMNNBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the n-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void EvaluateMNNBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
						     ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // core part of the AddMultiply method involving the n-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  virtual void HermitianEvaluateMNNBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the n-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void HermitianEvaluateMNNBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
							      ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // core part of the FastMultiplication method involving the n-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNNBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
							      int* indexArray, Complex* coefficientArray, long& position);

  // core part of the PartialFastMultiplicationMemory method involving n-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNNBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory);

  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										 int firstComponent, int nbrComponent);

  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // test the amount of memory needed for fast multiplication algorithm
  //
  // allowedMemory = amount of memory that cam be allocated for fast multiplication
  // return value = amount of memory needed
  virtual long FastMultiplicationMemory(long allowedMemory);

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // return value = number of non-zero matrix element
  virtual long PartialFastMultiplicationMemory(int firstComponent, int lastComponent);

  // enable fast multiplication algorithm
  //
  virtual void EnableFastMultiplication();

  // enable fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // nbrComponent  = index of the last component that has to be precalcualted
  virtual void PartialEnableFastMultiplication(int firstComponent, int nbrComponent);

  // compute the permutation array for the interaction indices
  // 
  // permutations = arrays where permuted indices will be stored
  // permutationSign = array with the sign of each permutation (initialized only when dealing with fermionic statistics)
  // return value = number of permutations
  int ComputePermutations(int**& permutations, double*& permutationSign);

};

// core part of the AddMultiply method involving the n-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeTimeReversalBreakingSingleBandNBodyHamiltonian::EvaluateMNNBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int Index;
 
  for (int j = 0; j < this->NbrNBodySectorSums; ++j)
    {
      int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
      TmpIndices = this->NBodySectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	{
	  Coefficient3 = particles->ProdA(index, TmpIndices + i1, this->NBodyValue);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->NBodyInteractionFactors[j][(i1 * Lim) / this->SqrNBodyValue]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		{
		  Index = particles->ProdAd(TmpIndices + i2, this->NBodyValue, Coefficient);
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


// core part of the AddMultiply method involving the n-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void ParticleOnLatticeTimeReversalBreakingSingleBandNBodyHamiltonian::EvaluateMNNBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
														   ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  for (int j = 0; j < this->NbrNBodySectorSums; ++j)
    {
      int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
      TmpIndices = this->NBodySectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	{
	  Coefficient3 = particles->ProdA(index, TmpIndices + i1, this->NBodyValue);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->NBodyInteractionFactors[j][(i1 * Lim) / this->SqrNBodyValue]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		{
		  Index = particles->ProdAd(TmpIndices + i2, this->NBodyValue, Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
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

inline void ParticleOnLatticeTimeReversalBreakingSingleBandNBodyHamiltonian::HermitianEvaluateMNNBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  //int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int Index;
  Complex TmpSum = 0.0;
  for (int j = 0; j < this->NbrNBodySectorSums; ++j)
    {
      int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
      TmpIndices = this->NBodySectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	{
	  Coefficient3 = particles->ProdA(index, TmpIndices + i1, this->NBodyValue);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->NBodyInteractionFactors[j][(i1 * Lim) / this->SqrNBodyValue]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		{
		  Index = particles->ProdAd(TmpIndices + i2, this->NBodyValue, Coefficient);
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

inline void ParticleOnLatticeTimeReversalBreakingSingleBandNBodyHamiltonian::HermitianEvaluateMNNBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
														 ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  //int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  Complex* TmpSum = new Complex[nbrVectors];
  for (int j = 0; j < this->NbrNBodySectorSums; ++j)
    {
      int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
      TmpIndices = this->NBodySectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	{
	  Coefficient3 = particles->ProdA(index, TmpIndices + i1, this->NBodyValue);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->NBodyInteractionFactors[j][(i1 * Lim) / this->SqrNBodyValue]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		{
		  Index = particles->ProdAd(TmpIndices + i2, this->NBodyValue, Coefficient);
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
  for (int l = 0; l < nbrVectors; ++l)
    vDestinations[l][index] += TmpSum[l];
  delete[] TmpSum;
}

// core part of the FastMultiplication method involving the n-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void ParticleOnLatticeTimeReversalBreakingSingleBandNBodyHamiltonian::EvaluateMNNBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index,
															int* indexArray, Complex* coefficientArray, long& position)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int Dim = particles->GetHilbertSpaceDimension();

  if (this->HermitianSymmetryFlag == false)
    {
      for (int j = 0; j < this->NbrNBodySectorSums; ++j)
	{
	  int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
	  TmpIndices = this->NBodySectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	    {
	      int AbsoluteIndex = index + this->PrecalculationShift;
	      Coefficient2 = particles->ProdA(AbsoluteIndex, TmpIndices + i1, this->NBodyValue);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->NBodyInteractionFactors[j][(i1 * Lim) / this->SqrNBodyValue]);
		  for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		    {
		      Index = particles->ProdAd(TmpIndices + i2, this->NBodyValue, Coefficient);
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
  else
    {
      for (int j = 0; j < this->NbrNBodySectorSums; ++j)
	{
	  int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
	  TmpIndices = this->NBodySectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	    {
	      int AbsoluteIndex = index + this->PrecalculationShift;
	      Coefficient2 = particles->ProdA(AbsoluteIndex, TmpIndices + i1, this->NBodyValue);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->NBodyInteractionFactors[j][(i1 * Lim) / this->SqrNBodyValue]);
		  for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		    {
		      Index = particles->ProdAd(TmpIndices + i2, this->NBodyValue, Coefficient);
		      if (Index <= AbsoluteIndex)
			{
			  if (Index == AbsoluteIndex)
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient2 * 0.5 * (*TmpInteractionFactor);
			      ++position;
			    }
			  else
			    {
			      indexArray[position] = Index;
			      coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
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


// core part of the PartialFastMultiplicationMemory method involving n-body term
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void ParticleOnLatticeTimeReversalBreakingSingleBandNBodyHamiltonian::EvaluateMNNBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  int Dim = particles->GetHilbertSpaceDimension();

  if (this->HermitianSymmetryFlag == true)
   {
     for (int i = firstComponent; i < lastComponent; ++i)
       {
         for (int j = 0; j < this->NbrNBodySectorSums; ++j)
	   {
	     int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
	     TmpIndices = this->NBodySectorIndicesPerSum[j];
	     for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	       {
	         Coefficient2 = particles->ProdA(i, TmpIndices + i1, this->NBodyValue);
	         if (Coefficient2 != 0.0)
		  {
		    for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		      {
		        Index = particles->ProdAd(TmpIndices + i2, this->NBodyValue, Coefficient);
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
  else
  {
   for (int i = firstComponent; i < lastComponent; ++i)
       {
         for (int j = 0; j < this->NbrNBodySectorSums; ++j)
	   {
	     int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
	     TmpIndices = this->NBodySectorIndicesPerSum[j];
	     for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	       {
	         Coefficient2 = particles->ProdA(i, TmpIndices + i1, this->NBodyValue);
	         if (Coefficient2 != 0.0)
		  {
		    for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		      {
		        Index = particles->ProdAd(TmpIndices + i2, this->NBodyValue, Coefficient);
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
