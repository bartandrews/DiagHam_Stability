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


#ifndef PARTICLEONCYLINDERNBODYHAMILTONIAN_H
#define PARTICLEONCYLINDERNBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeTimeReversalBreakingSingleBandNBodyHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnCylinderNBodyHamiltonian : public ParticleOnLatticeTimeReversalBreakingSingleBandNBodyHamiltonian
{

  // ratio between the width in the x direction and the width in the y direction
  double Ratio;

 protected:
  
  int MinSumIndices;
  int MaxSumIndices;

 public:

  // default constructor
  //
  ParticleOnCylinderNBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum momentum
  // nBody = order of n-body interaction
  // ratio = aspect ratio of the cylinder
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnCylinderNBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int maxMomentum, int nBody, double ratio, AbstractArchitecture* architecture, long memory = -1);

 protected:

  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // get all indices needed to characterize a completly symmetric tensor, sorted by the sum of the indices
  //
  // nbrValues = number of different values an index can have
  // nbrIndices = number of indices 
  // nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
  // sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
  // sortedIndicesPerSumSymmetryFactor = reference on a array where symmetry factor (aka product of the factorial of the number 
  //                                      of time each index appears) are stored (first array dimension corresponding to sum of the indices)
  // return value = total number of index groups

  long GetAllSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum, int**& sortedIndicesPerSumSymmetryFactor);

  // get all indices needed to characterize a completly skew symmetric tensor, sorted by the sum of the indices
  //
  // nbrValues = number of different values an index can have
  // nbrIndices = number of indices 
  // nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
  // sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
  // return value = total number of index groups

  long  GetAllSkewSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
										  int**& sortedIndicesPerSum);

  // evaluate the numerical coefficient  in front of the Prod a_m^+ Prod a_m coupling term for bosons
  // Indices1,Indices2 = arrays containing indices
  // return value = numerical coefficient

  Complex EvaluateInteractionCoefficient(int* Indices1, int* Indices2);

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


};

// core part of the AddMultiply method involving the n-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnCylinderNBodyHamiltonian::EvaluateMNNBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int Index;
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
   {
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
  else //Bosons... use ProdAL etc.
   {
    for (int j = 0; j < this->NbrNBodySectorSums; ++j)
      {
        int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
        TmpIndices = this->NBodySectorIndicesPerSum[j];
        for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	  {
	    Coefficient3 = particles->ProdAL(index, TmpIndices + i1, this->NBodyValue);
	    if (Coefficient3 != 0.0)
	      {
	        TmpInteractionFactor = &(this->NBodyInteractionFactors[j][(i1 * Lim) / this->SqrNBodyValue]);
	        Coefficient4 = vSource[index];
	        Coefficient4 *= Coefficient3;
	        for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		  {
		    Index = particles->ProdAdL(TmpIndices + i2, this->NBodyValue, Coefficient);
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


// core part of the AddMultiply method involving the n-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void ParticleOnCylinderNBodyHamiltonian::EvaluateMNNBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
														   ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
   {
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
  else //Bosons...use ProdAL etc.
  {
    for (int j = 0; j < this->NbrNBodySectorSums; ++j)
      {
        int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
        TmpIndices = this->NBodySectorIndicesPerSum[j];
        for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	  {
	    Coefficient3 = particles->ProdAL(index, TmpIndices + i1, this->NBodyValue);
	    if (Coefficient3 != 0.0)
	      {
	        TmpInteractionFactor = &(this->NBodyInteractionFactors[j][(i1 * Lim) / this->SqrNBodyValue]);
	        for (int p = 0; p < nbrVectors; ++p)
		  tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	        for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		  {
		    Index = particles->ProdAdL(TmpIndices + i2, this->NBodyValue, Coefficient);
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

// core part of the AddMultiply method involving the n-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnCylinderNBodyHamiltonian::HermitianEvaluateMNNBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  //int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int Index;
  Complex TmpSum = 0.0;

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
   {
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
  else //Bosons...use ProdAL etc.
   {
    for (int j = 0; j < this->NbrNBodySectorSums; ++j)
      {
        int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
        TmpIndices = this->NBodySectorIndicesPerSum[j];
        for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	  {
	    Coefficient3 = particles->ProdAL(index, TmpIndices + i1, this->NBodyValue);
	    if (Coefficient3 != 0.0)
	      {
	        TmpInteractionFactor = &(this->NBodyInteractionFactors[j][(i1 * Lim) / this->SqrNBodyValue]);
	        Coefficient4 = vSource[index];
	        Coefficient4 *= Coefficient3;
	        for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		  {
		    Index = particles->ProdAdL(TmpIndices + i2, this->NBodyValue, Coefficient);
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
}


// core part of the AddMultiply method involving the n-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void ParticleOnCylinderNBodyHamiltonian::HermitianEvaluateMNNBodyAddMultiplyComponent(ParticleOnSphere* particles, int index, ComplexVector* vSources, 
														 ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  //int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  Complex* TmpSum = new Complex[nbrVectors];

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
   {
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
   }
  else //Bosons...use ProdAL etc.
   {
    for (int j = 0; j < this->NbrNBodySectorSums; ++j)
      {
        int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
        TmpIndices = this->NBodySectorIndicesPerSum[j];
        for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	  {
	    Coefficient3 = particles->ProdAL(index, TmpIndices + i1, this->NBodyValue);
	    if (Coefficient3 != 0.0)
	      {
	        TmpInteractionFactor = &(this->NBodyInteractionFactors[j][(i1 * Lim) / this->SqrNBodyValue]);
	        for (int p = 0; p < nbrVectors; ++p)
		  tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	        for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		  {
		    Index = particles->ProdAdL(TmpIndices + i2, this->NBodyValue, Coefficient);
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
   }
  delete[] TmpSum;
}

// core part of the FastMultiplication method involving the n-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void ParticleOnCylinderNBodyHamiltonian::EvaluateMNNBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index,
														  int* indexArray, Complex* coefficientArray, long& position)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  Complex* TmpInteractionFactor;
  int Dim = particles->GetHilbertSpaceDimension();

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
   {
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
  else //Bosons...use ProdAL etc.
   {
    if (this->HermitianSymmetryFlag == false)
      {
        for (int j = 0; j < this->NbrNBodySectorSums; ++j)
	  {
	    int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
	    TmpIndices = this->NBodySectorIndicesPerSum[j];
	    for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	      {
	        int AbsoluteIndex = index + this->PrecalculationShift;
	        Coefficient2 = particles->ProdAL(AbsoluteIndex, TmpIndices + i1, this->NBodyValue);
	        if (Coefficient2 != 0.0)
		  {
		    TmpInteractionFactor = &(this->NBodyInteractionFactors[j][(i1 * Lim) / this->SqrNBodyValue]);
		    for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		      {
		        Index = particles->ProdAdL(TmpIndices + i2, this->NBodyValue, Coefficient);
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
	        Coefficient2 = particles->ProdAL(AbsoluteIndex, TmpIndices + i1, this->NBodyValue);
	        if (Coefficient2 != 0.0)
		  {
		    TmpInteractionFactor = &(this->NBodyInteractionFactors[j][(i1 * Lim) / this->SqrNBodyValue]);
		    for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		      {
		        Index = particles->ProdAdL(TmpIndices + i2, this->NBodyValue, Coefficient);
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
}


// core part of the PartialFastMultiplicationMemory method involving n-body term
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void ParticleOnCylinderNBodyHamiltonian::EvaluateMNNBodyFastMultiplicationMemoryComponent(ParticleOnSphere* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  int Dim = particles->GetHilbertSpaceDimension();

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
   {
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
  else //Bosons...use ProdAL etc.
   {
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
	           Coefficient2 = particles->ProdAL(i, TmpIndices + i1, this->NBodyValue);
	           if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		        {
		          Index = particles->ProdAdL(TmpIndices + i2, this->NBodyValue, Coefficient);
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
	           Coefficient2 = particles->ProdAL(i, TmpIndices + i1, this->NBodyValue);
	           if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		        {
		          Index = particles->ProdAdL(TmpIndices + i2, this->NBodyValue, Coefficient);
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
}

#endif
