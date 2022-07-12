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
//                                                                            //
//                        last modification : 23/02/2011                      //
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


#ifndef PARTICLEONLATTICETIMEREVERSALBREAKINGSINGLEBANDHAMILTONIAN_H
#define PARTICLEONLATTICETIMEREVERSALBREAKINGSINGLEBANDHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian : public AbstractQHEHamiltonian
{

 protected:
  
  // Hilbert space associated to the system
  ParticleOnSphere* Particles;

  // pointer to the tight binding model
  // number of particles
  int NbrParticles;

  // Max momentum that can be reached written as a linearized index
  int LzMax;

  // shift to apply to go from precalculation index to the corresponding index in the HilbertSpace
  int PrecalculationShift;

  // shift to apply to the Hamiltonian diagonal elements
  double HamiltonianShift;

  // amount of memory (in bytes) that can be used to store precalculated matrix elements
  long Memory;


  // number of index sum for the creation (or annhilation) operators for the intra spin sector
  int NbrSectorSums;
  // array containing the number of index group per index sum for the creation (or annhilation) operators for the intra spin sector
  int* NbrSectorIndicesPerSum;
  // array containing the (m1,m2) indices per index sum for the creation (or annhilation) operators for the intra spin sector
  int** SectorIndicesPerSum;

  // arrays containing all interaction factors, the first index correspond correspond to index sum for the creation (or annhilation) operators
  // the second index is a linearized index (m1,m2) + (n1,n2) * (nbr element in current index sum) (m for creation operators, n for annhilation operators)
  Complex** InteractionFactors;

  // array that contains all one-body interaction factors for particles with spin up
  double* OneBodyInteractionFactors;
  
  // flag for fast multiplication algorithm
  bool FastMultiplicationFlag;
  // step between each precalculated index
  int FastMultiplicationStep;

  // stored interactions per component
  int *NbrInteractionPerComponent;

  // number of tasks for load balancing
  int NbrBalancedTasks;
  // load balancing array for parallelisation, indicating starting indices
  long* LoadBalancingArray;

  // indices of matrix elements per component
  int **InteractionPerComponentIndex;
  // coefficients of matrix elements per component
  Complex** InteractionPerComponentCoefficient;
  // coefficients of matrix elements per component, stored as an unsigned integer array
  unsigned int** InteractionPerComponentIntegerCoefficient;
  // array where the coefficients are stored when using InteractionPerComponentIntegerCoefficient
  Complex* InteractionIntegerMatrixElementCoefficientTable;
  double* InteractionIntegerCreationCoefficientTable;
  double* InteractionIntegerAnnihiliationCoefficientTable;

  // translations of matrix elements per component
  int** InteractionPerComponentNbrTranslation;
  
  // flag for implementation of hermitian symmetry
  bool HermitianSymmetryFlag;
  
 public:

  // default constructor
  //
  ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian();

  // destructor
  //
  ~ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian();
  
  // ask if Hamiltonian implements hermitian symmetry operations
  //
  virtual bool IsHermitian();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  virtual void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  virtual AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  virtual int GetHilbertSpaceDimension ();
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  virtual void ShiftHamiltonian (double shift);

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

  // get the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
  //
  // nbrThreads = number of threads requested
  // segmentIndices = array returning the reference to an array of the first index of each of the segments
  // return value = true if no error occured
  virtual bool GetLoadBalancing(int nbrTasks, long* &segmentIndices);
    
 protected:
 
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
  virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
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

  // core part of the FastMultiplication method involving the one-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
							    int* indexArray, Complex* coefficientArray, long& position);

  // core part of the AddMultiply method involving the one-body interaction, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent,
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
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent,
						     int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors);

  // core part of the PartialFastMultiplicationMemory method involving two-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory);

  // core part of the PartialFastMultiplicationMemory method involving one-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory);

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

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int Index;
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
		  Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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


// core part of the AddMultiply method involving the two-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
													      ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
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
		  Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		  if (Index < Dim)
		    for (int p = 0; p < nbrVectors; ++p)
		      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
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

inline void ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int Index;
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
		  Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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


// core part of the AddMultiply method involving the two-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
														       ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  double Coefficient;
  double Coefficient3;
  int Index;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  Complex* TmpSum = new Complex[nbrVectors];
  for (int l = 0; l < nbrVectors; ++l)
    TmpSum[l] = 0.0;
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
		  Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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

// core part of the FastMultiplication method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
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
		      Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		      Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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


// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent,
													      int step, ComplexVector& vSource, ComplexVector& vDestination)
{
  if (this->OneBodyInteractionFactors != 0)
    {
      double TmpDiagonal = 0.0;
      for (int i = firstComponent; i < lastComponent; i += step)
	{ 
	  TmpDiagonal = 0.0;
	  for (int j = 0; j <= this->LzMax; ++j) 
	    {
	      TmpDiagonal += this->OneBodyInteractionFactors[j] * particles->AdA(i, j);
	    }
	  vDestination[i] += (this->HamiltonianShift + TmpDiagonal)* vSource[i];
	}
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

inline void ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent,
													      int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  if (this->OneBodyInteractionFactors != 0) 
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
		  TmpDiagonal += this->OneBodyInteractionFactors[j] * particles->AdA(i, j);
		}
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

// core part of the FastMultiplication method involving the one-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
														     int* indexArray, Complex* coefficientArray, long& position)
{
  if (this->OneBodyInteractionFactors != 0)
    {
      double TmpDiagonal = 0.0;
      for (int j = 0; j <= this->LzMax; ++j)
	TmpDiagonal += this->OneBodyInteractionFactors[j] * particles->AdA(index + this->PrecalculationShift, j);
      if (this->HermitianSymmetryFlag == true)
	TmpDiagonal *= 0.5;
      indexArray[position] = index + this->PrecalculationShift;
      coefficientArray[position] = TmpDiagonal;
      ++position;	  
    }
}

// core part of the PartialFastMultiplicationMemory method involving two-body term
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;

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
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AdAd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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

// core part of the PartialFastMultiplicationMemory method involving two-body term
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory)
{
  if (this->OneBodyInteractionFactors != 0)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  ++memory;
	  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
	}
    }
}



#endif
