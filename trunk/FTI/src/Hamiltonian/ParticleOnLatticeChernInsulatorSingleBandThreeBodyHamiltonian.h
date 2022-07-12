////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                   class of hamiltonian with particles on                   //
//               Chern insulator in the single band approximation             //
//                         and three-body interaction                         //
//                                                                            //
//                        last modification : 13/07/2011                      //
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


#ifndef PARTICLEONLATTICECHERNINSULATORSINGLEBANDTHREEBODYHAMILTONIAN_H
#define PARTICLEONLATTICECHERNINSULATORSINGLEBANDTHREEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnLatticeChernInsulatorSingleBandThreeBodyHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandHamiltonian
{

 protected:
  
  // number of index sum for the creation (or annhilation) operators for the intra spin sector
  int NbrThreeBodySectorSums;
  // array containing the number of index group per index sum for the creation (or annhilation) operators for the intra spin sector
  int* NbrThreeBodySectorIndicesPerSum;
  // array containing the (m1,m2,m3) indices per index sum for the creation (or annhilation) operators for the intra spin sector
  int** ThreeBodySectorIndicesPerSum;

  // arrays containing all interaction factors, the first index correspond correspond to index sum for the creation (or annhilation) operators
  // the second index is a linearized index (m1,m2,m3) + (n1,n2,n3) * (nbr element in current index sum) (m for creation operators, n for annhilation operators)
  Complex** ThreeBodyInteractionFactors;

  // true if one/two body terms are set
  bool TwoBodyFlag;  
  
 public:

  // default constructor
  //
  ParticleOnLatticeChernInsulatorSingleBandThreeBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // tightBindingModel = pointer to the tight binding model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeChernInsulatorSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY,
          Abstract2DTightBindingModel* tightBindingModel, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeChernInsulatorSingleBandThreeBodyHamiltonian();
  
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
 
  // core part of the AddMultiply method involving the three-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  virtual void EvaluateMNThreeBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the three-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void EvaluateMNThreeBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
						     ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // core part of the AddMultiply method involving the three-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  virtual void HermitianEvaluateMNThreeBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the three-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void HermitianEvaluateMNThreeBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
							      ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // core part of the FastMultiplication method involving the three-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNThreeBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
							      int* indexArray, Complex* coefficientArray, long& position);

  // core part of the PartialFastMultiplicationMemory method involving three-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNThreeBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory);

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

};

// core part of the AddMultiply method involving the three-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeChernInsulatorSingleBandThreeBodyHamiltonian::EvaluateMNThreeBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int Index;
  for (int j = 0; j < this->NbrThreeBodySectorSums; ++j)
    {
      int Lim = 3 * this->NbrThreeBodySectorIndicesPerSum[j];
      TmpIndices = this->ThreeBodySectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 3)
	{
	  Coefficient3 = particles->ProdA(index, TmpIndices + i1, 3);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->ThreeBodyInteractionFactors[j][(i1 * Lim) / 9]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 3)
		{
		  Index = particles->ProdAd(TmpIndices + i2, 3, Coefficient);
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


// core part of the AddMultiply method involving the three-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void ParticleOnLatticeChernInsulatorSingleBandThreeBodyHamiltonian::EvaluateMNThreeBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
														   ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  for (int j = 0; j < this->NbrThreeBodySectorSums; ++j)
    {
      int Lim = 3 * this->NbrThreeBodySectorIndicesPerSum[j];
      TmpIndices = this->ThreeBodySectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 3)
	{
	  Coefficient3 = particles->ProdA(index, TmpIndices + i1, 3);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->ThreeBodyInteractionFactors[j][(i1 * Lim) / 9]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 3)
		{
		  Index = particles->ProdAd(TmpIndices + i2, 3, Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		  ++TmpInteractionFactor;
		}
	    }
	}
    }
}

// core part of the AddMultiply method involving the three-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeChernInsulatorSingleBandThreeBodyHamiltonian::HermitianEvaluateMNThreeBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int Index;
  Complex TmpSum = 0.0;
  for (int j = 0; j < this->NbrThreeBodySectorSums; ++j)
    {
      int Lim = 3 * this->NbrThreeBodySectorIndicesPerSum[j];
      TmpIndices = this->ThreeBodySectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 3)
	{
	  Coefficient3 = particles->ProdA(index, TmpIndices + i1, 3);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->ThreeBodyInteractionFactors[j][(i1 * Lim) / 9]);
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      for (int i2 = 0; i2 < Lim; i2 += 3)
		{
		  Index = particles->ProdAd(TmpIndices + i2, 3, Coefficient);
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


// core part of the AddMultiply method involving the three-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void ParticleOnLatticeChernInsulatorSingleBandThreeBodyHamiltonian::HermitianEvaluateMNThreeBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
														 ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  Complex* TmpSum = new Complex[nbrVectors];
  for (int j = 0; j < this->NbrThreeBodySectorSums; ++j)
    {
      int Lim = 3 * this->NbrThreeBodySectorIndicesPerSum[j];
      TmpIndices = this->ThreeBodySectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += 3)
	{
	  Coefficient3 = particles->ProdA(index, TmpIndices + i1, 3);
	  if (Coefficient3 != 0.0)
	    {
	      TmpInteractionFactor = &(this->ThreeBodyInteractionFactors[j][(i1 * Lim) / 9]);
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      for (int i2 = 0; i2 < Lim; i2 += 3)
		{
		  Index = particles->ProdAd(TmpIndices + i2, 3, Coefficient);
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

// core part of the FastMultiplication method involving the three-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void ParticleOnLatticeChernInsulatorSingleBandThreeBodyHamiltonian::EvaluateMNThreeBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index,
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
      for (int j = 0; j < this->NbrThreeBodySectorSums; ++j)
	{
	  int Lim = 3 * this->NbrThreeBodySectorIndicesPerSum[j];
	  TmpIndices = this->ThreeBodySectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 3)
	    {
	      int AbsoluteIndex = index + this->PrecalculationShift;
	      Coefficient2 = particles->ProdA(AbsoluteIndex, TmpIndices + i1, 3);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->ThreeBodyInteractionFactors[j][(i1 * Lim) / 9]);
		  for (int i2 = 0; i2 < Lim; i2 += 3)
		    {
		      Index = particles->ProdAd(TmpIndices + i2, 3, Coefficient);
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
      for (int j = 0; j < this->NbrThreeBodySectorSums; ++j)
	{
	  int Lim = 3 * this->NbrThreeBodySectorIndicesPerSum[j];
	  TmpIndices = this->ThreeBodySectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 3)
	    {
	      int AbsoluteIndex = index + this->PrecalculationShift;
	      Coefficient2 = particles->ProdA(AbsoluteIndex, TmpIndices + i1, 3);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->ThreeBodyInteractionFactors[j][(i1 * Lim) / 9]);
		  for (int i2 = 0; i2 < Lim; i2 += 3)
		    {
		      Index = particles->ProdAd(TmpIndices + i2, 3, Coefficient);
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


// core part of the PartialFastMultiplicationMemory method involving three-body term
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void ParticleOnLatticeChernInsulatorSingleBandThreeBodyHamiltonian::EvaluateMNThreeBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  int Dim = particles->GetHilbertSpaceDimension();
  int SumIndices;
  int TmpNbrM3Values;
  int* TmpM3Values;


  for (int i = firstComponent; i < lastComponent; ++i)
    {
      for (int j = 0; j < this->NbrThreeBodySectorSums; ++j)
	{
	  int Lim = 3 * this->NbrThreeBodySectorIndicesPerSum[j];
	  TmpIndices = this->ThreeBodySectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 3)
	    {
	      Coefficient2 = particles->ProdA(i, TmpIndices + i1, 3);
	      if (Coefficient2 != 0.0)
		{
		  for (int i2 = 0; i2 < Lim; i2 += 3)
		    {
		      Index = particles->ProdAd(TmpIndices + i2, 3, Coefficient);
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




#endif
