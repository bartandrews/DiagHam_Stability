////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of abstract fractional quantum Hall hamiltonian          //
//              associated to particles with SU(2) spin on a sphere           //
//                                                                            //
//                        last modification : 07/06/2007                      //
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


#ifndef ABSTRACTQHEONSPHEREQUANTUMWELLHAMILTONIAN_H
#define ABSTRACTQHEONSPHEREQUANTUMWELLHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/AbstractQHEOnSphereHamiltonian.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnSphereWithSpinL2Hamiltonian;
class ParticleOnSphereWithSpinS2Hamiltonian;


class AbstractQHEOnSphereQuantumWellHamiltonian : public AbstractQHEOnSphereHamiltonian
{

 protected:
  
  // number of index sum for the creation (or annhilation) operators for the intra spin sector
  int NbrIntraSectorSums;
  // array containing the number of index group per index sum for the creation (or annhilation) operators for the intra spin sector
  int* NbrIntraSectorIndicesPerSum;
  // array containing the (m1,m2) indices per index sum for the creation (or annhilation) operators for the intra spin sector
  int** IntraSectorIndicesPerSum;

  // number of index sum for the creation (or annhilation) operators for the inter spin sector
  int NbrInterSectorSums;
  // array containing the number of index group per index sum for the creation (or annhilation) operators for the inter spin sector
  int* NbrInterSectorIndicesPerSum;
  // array containing the (m1,m2) indices per index sum for the creation (or annhilation) operators for the inter spin sector
  int** InterSectorIndicesPerSum;

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//^^^^^^^^^^^^^^^^^^^^^ M I X E D ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // number of index sum for the creation (or annhilation) operators for the mixed spin sector
  int NbrMixedIntraSectorSums;
  // array containing the number of index group per index sum for the creation (or annhilation) operators for the mixed spin sector
  int* NbrMixedIntraSectorIndicesPerSum;
  // array containing the (m1,m2) indices per index sum for the creation (or annhilation) operators for the mixed spin sector
  int** MixedIntraSectorIndicesPerSum;

  // number of index sum for the creation (or annhilation) operators for the mixed spin sector
  int NbrMixedInterSectorSums;
  // array containing the number of index group per index sum for the creation (or annhilation) operators for the mixed spin sector
  int* NbrMixedInterSectorIndicesPerSum;
  // array containing the (m1,m2) indices per index sum for the creation (or annhilation) operators for the mixed spin sector
  int** MixedInterSectorIndicesPerSum;

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//^^^^^^^^^^^^^^^^^^^^^ M I X E D ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


  // alternative method used to store indices attached to each interaction factor
  //  arrays of m1 and m2 intra-sector indices (i.e. indices for the a_m1_s a_m2_s factors)
  int* M1IntraValue;
  int* M2IntraValue;
  //  number of element in each m1 or m2 intra-sector index array
  int NbrM12IntraIndices;
  //  the m3 intra-sector index array (i.e. index for the ad_m3_s a_(m1+m3-m3)_s factors) for each set of (m1,m2)
  int** M3IntraValues;
  // number of possible m3 intra-sector index for each set of (m1,m2)
  int* NbrM3IntraValues;

  // alternative method used to store indices attached to each interaction factor
  //  arrays of m1 and m2 inter-sector indices (i.e. indices for the a_m1_s a_m2_-s factors)
  int* M1InterValue;
  int* M2InterValue;  
  //  number of element in each m1 or m2 inter-sector index array
  int NbrM12InterIndices;
  //  the m3 inter-sector index array (i.e. index for the ad_m3_s a_(m1+m3-m3)_-s factors) for each set of (m1,m2)
  int** M3InterValues;
  // number of possible m3 inter-sector index for each set of (m1,m2)
  int* NbrM3InterValues;

  // arrays containing all interaction factors, the first index correspond correspond to index sum for the creation (or annhilation) operators
  // the second index is a linearized index (m1,m2) + (n1,n2) * (nbr element in current index sum) (m for creation operators, n for annhilation operators)
  // array containing all interaction factors for spin up and spin up
  double** InteractionFactorsupup;
  // array containing all interaction factors for spin down and spin down
  double** InteractionFactorsdowndown;
  // array containing all interaction factors for spin up and spin down
  double** InteractionFactorsupdown;


//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//^^^^^^^^^^^^^^^^^^^^^ M I X E D ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // array containing all interaction factors for mixed spin (up-up-down-down)
  double** InteractionFactorsmixedintra;
  // array containing all interaction factors for mixed spin (up-up-down-down)
  double** InteractionFactorsmixedinter;
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//^^^^^^^^^^^^^^^^^^^^^ M I X E D ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


  // arrays containing all interaction factors needed for the alternative method
  // array containing all interaction factors for spin up and spin up
  double* M12InteractionFactorsupup;
  // array containing all interaction factors for spin down and spin down
  double* M12InteractionFactorsdowndown;
  // array containing all interaction factors for spin up and spin down
  double* M12InteractionFactorsupdown;

  // array that contains all one-body interaction factors for particles with spin up
  double* OneBodyInteractionFactorsupup;
  // array that contains all one-body interaction factors for particles with spin down
  double* OneBodyInteractionFactorsdowndown;
  // array that contains all one-body interaction factors for tunnelling terms for particles with different spin
  double* OneBodyInteractionFactorsupdown;
  
  // pointer to an optional L^2 operator in the Hamiltonian 
  ParticleOnSphereWithSpinL2Hamiltonian* L2Hamiltonian;
  // pointer to an optional S^2 operator in the Hamiltonian 
  ParticleOnSphereWithSpinS2Hamiltonian* S2Hamiltonian;

 public:

  // destructor
  //
  virtual ~AbstractQHEOnSphereQuantumWellHamiltonian() = 0;

  // ask if Hamiltonian implements methods using hermitian symmetry 
  //
  virtual bool IsHermitian();

  // ask if Hamiltonian implements methods applying the conjugate of the Hamiltonian
  //
  virtual bool IsConjugate();

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
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
  virtual  RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						   int firstComponent, int nbrComponent);

  // add an additional S^2 term to the Hamiltonian
  //
  // totalLz = twice the projected momentum total value
  // totalSz = twice the projected spin total value
  // factor = factor in front of the S^2
  // memory = amount of memory that can be used for S^2  precalculations 
  void AddS2 (int totalLz, int totalSz, double factor = 1.0, long memory = 0l);

  // add an additional L^2 term to the Hamiltonian
  //
  // totalLz = twice the projected momentum total value
  // totalSz = twice the projected spin total value
  // factor = factor in front of the L^2
  // memory = amount of memory that can be used for L^2  precalculations   
  void AddL2 (int totalLz, int totalSz, double factor = 1.0, long memory = 0l);


 protected:
  
  // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector& vSource, RealVector& vDestination);

  // core part of the AddMultiply method involving the two-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector* vSources, 
						     RealVector* vDestinations, int nbrVectors, double* tmpCoefficients);

  // core part of the FastMultiplication method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
							  int* indexArray, double* coefficientArray, long& position);

  // core part of the AddMultiply method involving the one-body interaction, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
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
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
						     int step, RealVector* vSources, RealVector* vDestinations, int nbrVectors);

  // core part of the PartialFastMultiplicationMemory method involving two-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory);

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
  RealVector* LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							     int firstComponent, int nbrComponent);

  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors() = 0;

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
  // lastComponent  = index of the last component that has to be precalcualted
  virtual void PartialEnableFastMultiplication(int firstComponent, int lastComponent);

  // enable fast multiplication algorithm
  //
  // fileName = prefix of the name of the file where temporary matrix elements will be stored
  virtual void EnableFastMultiplicationWithDiskStorage(char* fileName);

};


// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnSphereQuantumWellHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector& vSource, RealVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int* TmpIndices;
  double* TmpInteractionFactor;
  int Index;
  if ( (this->NbrIntraSectorSums != 0) || (this->NbrInterSectorSums != 0) || (this->NbrMixedIntraSectorSums != 0) || (this->NbrMixedInterSectorSums != 0))
    {
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
		  Coefficient3 *= vSource[index];
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
		  Coefficient3 *= vSource[index];
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
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
		  Coefficient3 *= vSource[index];
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		      ++TmpInteractionFactor;
		    }
		}
	    } 
	}


//*********************************************************************************************************
//********************************     M I X E D     T E R M S ***********************************************
//*********************************************************************************************************

     for (int j = 0; j < this->NbrMixedIntraSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrMixedIntraSectorIndicesPerSum[j];
	  TmpIndices = this->MixedIntraSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsmixedintra[j][(i1 * Lim) >> 2]);
		  Coefficient3 *= vSource[index];
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		      ++TmpInteractionFactor;
		    }
		  }
	      Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsmixedintra[j][(i1 * Lim) >> 2]);
		  Coefficient3 *= vSource[index];
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		      ++TmpInteractionFactor;
		    }
		  }

	    }
	}

      for (int j = 0; j < this->NbrMixedInterSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrMixedInterSectorIndicesPerSum[j];
	  TmpIndices = this->MixedInterSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      //Change the order of m3 and m4!
	      Coefficient3 = particles->AuAd(index, TmpIndices[i1 + 1], TmpIndices[i1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsmixedinter[j][(i1 * Lim) >> 2]);
		  Coefficient3 *= vSource[index];
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      //The order of m1 and m2 remains the same
		      Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		      ++TmpInteractionFactor;
		    }
		}
	    } 
	}

//*********************************************************************************************************
//********************************     M I X E D     T E R M S ***********************************************
//*********************************************************************************************************
    }
  else
    {
      double Coefficient2;
      int SumIndices;
      int TmpNbrM3Values;
      int* TmpM3Values;
      int ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  Coefficient = particles->AuAu(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      Coefficient *= vSource[index];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AduAdu(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)			
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsupup[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
      ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  Coefficient = particles->AdAd(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      Coefficient *= vSource[index];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AddAdd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)			
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsdowndown[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
      ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AuAd(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      Coefficient *= vSource[index];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AduAdd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)			
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsupdown[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
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

inline void AbstractQHEOnSphereQuantumWellHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector* vSources, RealVector* vDestinations, int nbrVectors, double* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  
  if (	(this->NbrIntraSectorSums != 0) || (this->NbrInterSectorSums != 0) || (this->NbrMixedIntraSectorSums != 0)|| (this->NbrMixedInterSectorSums != 0)	)
    {
      int* TmpIndices;
      double* TmpInteractionFactor;
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
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			for (int p = 0; p < nbrVectors; ++p)
			  vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
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
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			for (int p = 0; p < nbrVectors; ++p)
			  vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
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
		      Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			for (int p = 0; p < nbrVectors; ++p)
			  vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		      ++TmpInteractionFactor;
		    }
		}
	    }
	}

//*********************************************************************************************************
//********************************     M I X E D     T E R M S ***********************************************
//*********************************************************************************************************
    for (int j = 0; j < this->NbrMixedIntraSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrMixedIntraSectorIndicesPerSum[j];
	  TmpIndices = this->MixedIntraSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsmixedintra[j][(i1 * Lim) >> 2]);
		  for (int p = 0; p < nbrVectors; ++p)
		    tmpCoefficients[p] = Coefficient3 * vSources[p][index];
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			for (int p = 0; p < nbrVectors; ++p)
			  vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		      ++TmpInteractionFactor;
		    }
		}
	      Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsmixedintra[j][(i1 * Lim) >> 2]);
		  for (int p = 0; p < nbrVectors; ++p)
		    tmpCoefficients[p] = Coefficient3 * vSources[p][index];
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			for (int p = 0; p < nbrVectors; ++p)
			  vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		      ++TmpInteractionFactor;
		    }
		}

	    }
	  }


      for (int j = 0; j < this->NbrMixedInterSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrMixedInterSectorIndicesPerSum[j];
	  TmpIndices = this->MixedInterSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      //Change the order of m3 and m4!
	      Coefficient3 = particles->AuAd(index, TmpIndices[i1 + 1], TmpIndices[i1]);
	      if (Coefficient3 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsmixedinter[j][(i1 * Lim) >> 2]);
		  for (int p = 0; p < nbrVectors; ++p)
		    tmpCoefficients[p] = Coefficient3 * vSources[p][index];
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      //Order of m1 and m2 remains the same.
		      Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			for (int p = 0; p < nbrVectors; ++p)
			  vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		      ++TmpInteractionFactor;
		    }
		}
	    }
	}

//*********************************************************************************************************
//********************************     M I X E D     T E R M S ***********************************************
//*********************************************************************************************************
    }
  else
    {
      double Coefficient2;
      int SumIndices;
      int TmpNbrM3Values;
      int* TmpM3Values;
      int ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  Coefficient = particles->AuAu(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      for (int l = 0; l < nbrVectors; ++l)
		tmpCoefficients[l] = Coefficient * vSources[l][index];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AduAdu(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][Index] += tmpCoefficients[l] * this->M12InteractionFactorsupup[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
      ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  Coefficient = particles->AdAd(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      for (int l = 0; l < nbrVectors; ++l)
		tmpCoefficients[l] = Coefficient * vSources[l][index];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AddAdd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][Index] += tmpCoefficients[l] * this->M12InteractionFactorsdowndown[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
      ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AuAd(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      for (int l = 0; l < nbrVectors; ++l)
		tmpCoefficients[l] = Coefficient * vSources[l][index];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AduAdd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][Index] += tmpCoefficients[l] * this->M12InteractionFactorsupdown[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
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

inline void AbstractQHEOnSphereQuantumWellHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, int* indexArray, double* coefficientArray, long& position)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  double* TmpInteractionFactor;
  int Dim = particles->GetHilbertSpaceDimension();
  if ((this->NbrIntraSectorSums != 0) || (this->NbrInterSectorSums != 0) || (this->NbrMixedIntraSectorSums != 0) || (this->NbrMixedInterSectorSums != 0)	)
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
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
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
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
		      Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
//*********************************************************************************************************
//********************************     M I X E D     T E R M S ***********************************************
//*********************************************************************************************************
     for (int j = 0; j < this->NbrMixedIntraSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrMixedIntraSectorIndicesPerSum[j];
	  TmpIndices = this->MixedIntraSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      Coefficient2 = particles->AdAd(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsmixedintra[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}

	      Coefficient2 = particles->AuAu(index + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsmixedintra[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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


     for (int j = 0; j < this->NbrMixedInterSectorSums; ++j)
	{
	  int Lim = 2 * this->NbrMixedInterSectorIndicesPerSum[j];
	  TmpIndices = this->MixedInterSectorIndicesPerSum[j];
	  for (int i1 = 0; i1 < Lim; i1 += 2)
	    {
	      //Change order of m3 and m4!
	      Coefficient2 = particles->AuAd(index + this->PrecalculationShift,  TmpIndices[i1 + 1], TmpIndices[i1]);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsmixedinter[j][(i1 * Lim) >> 2]);
		  for (int i2 = 0; i2 < Lim; i2 += 2)
		    {
		      //Order of m1 and m2 is the same.
		      Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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


//*********************************************************************************************************
//********************************     M I X E D     T E R M S ***********************************************
//*********************************************************************************************************
    }
  else
    {
      int SumIndices;
      int TmpNbrM3Values;
      int* TmpM3Values;
      int ReducedNbrInteractionFactorsupup;
      int ReducedNbrInteractionFactorsupdown;
      int ReducedNbrInteractionFactorsdowndown;      
      ReducedNbrInteractionFactorsupup = 0;
      ReducedNbrInteractionFactorsdowndown = 0;
      
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  Coefficient = particles->AuAu(index + this->PrecalculationShift, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AduAdu(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsupup[ReducedNbrInteractionFactorsupup];
		      ++position;
		    }		      
		  ++ReducedNbrInteractionFactorsupup;
		}    
	    }
	  else
	    ReducedNbrInteractionFactorsupup += this->NbrM3IntraValues[m1];
	  Coefficient = particles->AdAd(index + this->PrecalculationShift, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AddAdd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsdowndown[ReducedNbrInteractionFactorsdowndown];
		      ++position;
		    }		      
		  ++ReducedNbrInteractionFactorsdowndown;
		}    
	    }
	  else
	    ReducedNbrInteractionFactorsdowndown += this->NbrM3IntraValues[m1];
	}
      ReducedNbrInteractionFactorsupdown = 0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AuAd(index + this->PrecalculationShift, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AduAdd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsupdown[ReducedNbrInteractionFactorsupdown];
		      ++position;
		    }		      
		  ++ReducedNbrInteractionFactorsupdown;
		}    
	    }
	  else
	    ReducedNbrInteractionFactorsupdown += this->NbrM3InterValues[m1];
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

inline void AbstractQHEOnSphereQuantumWellHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, int step, RealVector& vSource, RealVector& vDestination)
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
      double Source;
      int Dim = particles->GetHilbertSpaceDimension();
      int Index;
      for (int i = firstComponent; i < lastComponent; i += step)
	{
	  Source = vSource[i];
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->AddAu(i + this->PrecalculationShift, j, j, Coefficient);
	      if (Index < Dim)
		{
		  vDestination[Index] += Coefficient * OneBodyInteractionFactorsupdown[j] * Source;
		}
	      Index = particles->AduAd(i + this->PrecalculationShift, j, j, Coefficient);
	      if (Index < Dim)
		{
		  vDestination[Index] += Coefficient * OneBodyInteractionFactorsupdown[j] * Source;
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

inline void AbstractQHEOnSphereQuantumWellHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, int step, RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  if (this->OneBodyInteractionFactorsupup != 0) 
    if (this->OneBodyInteractionFactorsdowndown != 0)
      {
	double TmpDiagonal = 0.0;
	for (int p = 0; p < nbrVectors; ++p)
	  {
	    RealVector& TmpSourceVector = vSources[p];
	    RealVector& TmpDestinationVector = vDestinations[p];
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
	    RealVector& TmpSourceVector = vSources[p];
	    RealVector& TmpDestinationVector = vDestinations[p];
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		TmpDiagonal = 0.0;
		for (int j = 0; j <= this->LzMax; ++j) 
		  TmpDiagonal += this->OneBodyInteractionFactorsupup[j] * particles->AduAu(i, j);
		TmpDestinationVector[i] += (this->HamiltonianShift + TmpDiagonal)* TmpSourceVector[i];
	      }
	  }
      }
  else
    if (this->OneBodyInteractionFactorsdowndown != 0)
      {
	double TmpDiagonal = 0.0;
	for (int p = 0; p < nbrVectors; ++p)
	  {
	    RealVector& TmpSourceVector = vSources[p];
	    RealVector& TmpDestinationVector = vDestinations[p];
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
	  RealVector& TmpSourceVector = vSources[p];
	  RealVector& TmpDestinationVector = vDestinations[p];
	  for (int i = firstComponent; i < lastComponent; i += step)
	    TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	}
  for (int p = 0; p < nbrVectors; ++p)
    {
      RealVector& TmpSourceVector = vSources[p];
      RealVector& TmpDestinationVector = vDestinations[p];
      for (int i = firstComponent; i < lastComponent; i += step)
	TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
    }
  if (this->OneBodyInteractionFactorsupdown != 0)
    {
      double Coefficient;
      int Dim = particles->GetHilbertSpaceDimension();
      int Index;
      for (int i = firstComponent; i < lastComponent; i += step)
	{
	  for (int j = 0; j <= this->LzMax; ++j)
	    {
	      Index = particles->AddAu(i + this->PrecalculationShift, j, j, Coefficient);
	      if (Index < Dim)
		{
		  for (int p = 0; p < nbrVectors; ++p)
		    vDestinations[p][Index] += Coefficient * OneBodyInteractionFactorsupdown[j] * vSources[p][i];
		}
	      Index = particles->AduAd(i + this->PrecalculationShift, j, j, Coefficient);
	      if (Index < Dim)
		{
		  for (int p = 0; p < nbrVectors; ++p)
		    vDestinations[p][Index] += Coefficient * OneBodyInteractionFactorsupdown[j] * vSources[p][i];
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

inline void AbstractQHEOnSphereQuantumWellHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  // double* TmpInteractionFactor;
  int Dim = particles->GetHilbertSpaceDimension();
  int SumIndices;
  int TmpNbrM3Values;
  int* TmpM3Values;

  for (int i = firstComponent; i < lastComponent; ++i)
    {
      if ((this->NbrIntraSectorSums != 0) || (this->NbrInterSectorSums != 0)|| (this->NbrMixedIntraSectorSums != 0)|| (this->NbrMixedInterSectorSums != 0) )
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
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
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
			  Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			}
		    }
		}
	    }	

//*********************************************************************************************************
//********************************     M I X E D     T E R M S ***********************************************
//*********************************************************************************************************
	  for (int j = 0; j < this->NbrMixedIntraSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrMixedIntraSectorIndicesPerSum[j];
	      TmpIndices = this->MixedIntraSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient2 = particles->AdAd(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			}
		    }

		  Coefficient2 = particles->AuAu(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			}
		    }


		}
	    }


	  for (int j = 0; j < this->NbrMixedInterSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrMixedInterSectorIndicesPerSum[j];
	      TmpIndices = this->MixedInterSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  //Change order of m3 and m4!
		  Coefficient2 = particles->AuAd(i, TmpIndices[i1 + 1], TmpIndices[i1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  //Order of m1 and m2 is the same.
			  Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			}
		    }
		}
	    }	


//*********************************************************************************************************
//********************************     M I X E D     T E R M S ***********************************************
//*********************************************************************************************************
	}    
      else
	{
	  for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	    {
	      Coefficient = particles->AuAu(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		  TmpM3Values = this->M3IntraValues[m1];
		  TmpNbrM3Values = this->NbrM3IntraValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      if (particles->AduAdu(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient) < this->Particles->GetHilbertSpaceDimension())
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }    
		}
	      Coefficient = particles->AdAd(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		  TmpM3Values = this->M3IntraValues[m1];
		  TmpNbrM3Values = this->NbrM3IntraValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      if (particles->AddAdd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient) < this->Particles->GetHilbertSpaceDimension())
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }    
		}
	    }
	  for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	    {
	      Coefficient = particles->AuAd(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      if (particles->AduAdd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient) < this->Particles->GetHilbertSpaceDimension())
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
