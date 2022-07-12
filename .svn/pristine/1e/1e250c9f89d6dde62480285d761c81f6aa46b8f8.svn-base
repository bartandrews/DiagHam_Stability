////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                          hardcore n-body interaction                       //
//                                                                            //
//                        last modification : 29/07/2014                      //
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


#ifndef PARTICLEONTORUSGENERICNBODYWITHMAGNETICTRANSLATIONSHAMILTONIAN_H
#define PARTICLEONTORUSGENERICNBODYWITHMAGNETICTRANSLATIONSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "Hamiltonian/AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian.h"
#include "MathTools/FactorialCoefficient.h" 

#include "Polynomial/Polynomial.h"

#include <iostream>
#include <algorithm>


using std::ostream;



class ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian : public AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian
{

 protected:

  // array where the interaction coefficients are stored
  double** PrecalculatedInteractionCoefficients;
  // number of rows in PrecalculatedInteractionCoefficients
  int NbrEntryPrecalculatedInteractionCoefficients1;
  // number of columns in PrecalculatedInteractionCoefficients
  int NbrEntryPrecalculatedInteractionCoefficients2;

  // numerical prefactor in front of the many-body interaction
  double NBodyPrefactor;

  // Number of Pseudopotential
  int NbrPseudopotentials;
  // pseudopotential coefficients
  double* Pseudopotentials;
  // Laguerre polynomial for the pseudopotentials
  Polynomial* LaguerrePolynomials;

 public:

  // default constructor
  //
  ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // nbrNBody = value of the n (i.e. the n-body interaction)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio, int nbrNBody,
								  AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
								  char* precalculationFileName = 0);

  // destructor
  //
  virtual ~ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian();
  	
  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  virtual void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  virtual void ShiftHamiltonian (double shift);
  
 protected:
  
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();
	
  // evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_nj coupling term
  //
  // creationCoefficients = array that contains the creation coefficients
  // annihilationCoefficients = array that contains the annihilation coefficients
  // return value = numerical coefficient  
  virtual double EvaluateInteractionCoefficient(double* creationCoefficients, double* annihilationCoefficients);
  
  // evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_n, assuming bosonic operators
  //
  // mIndices = array contaning the creation indices
  // nIndices = array contaning the annihilation indices
  // return value = numerical coefficient
  virtual double EvaluateInteractionNIndexSymmetrizedCoefficient(int* mIndices, int* nIndices);

  // evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_n, assuming fermionic operators
  //
  // mIndices = array contaning the creation indices
  // nIndices = array contaning the annihilation indices
  // return value = numerical coefficient
  virtual double EvaluateInteractionNIndexAntiSymmetrizedCoefficient(int* mIndices, int* nIndices);

  // evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_nj coupling term (factor corresponding to the creation operators only) for each integer modulo the NBodyValue
  //
  // mIndices = array that contains the creation indices
  // momentumTransfer = momentum transfer operated by the \prod_i a+_mi \prod_j a_nj operator, in units of the number of flux quanta
  // return value = array of numerical coefficients 
  virtual double* EvaluateInteractionCoefficientCreation(int* mIndices, int momentumTransfer);
  
  
  // evaluate the two nested Gaussian sum for a three body interaction (for test purposes)
  //
  // momFactor = array of indices that contains the information of the creation (or annihilation) indices
  // TmpIndices = array of indices that gives the initial indices that will be incremented in the sum
  // return value = value of the sum
  virtual double EvaluateGaussianSum(int* momFactor, int* TmpIndices);
  
  
  // evaluate the N nested infinite sums of EvaluateInteractionCoefficientCreation
  //
  // nBodyValue = current index being incremented 
  //TmpIndices = array containing the current value of all indices
  // Sum = current value of the sum
  // countIter = array of integers giving, for each TmpIndices[i], the number of times it has been incremented
  // momFactor = array of indices that contains the information of the creation (or annihilation) indices
  //return value = value of the coefficient
  virtual double EvaluateGaussianSum(int nBodyValue, int* TmpIndices, double Sum, int* countIter, int* momFactor);
  
   
  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  virtual double EvaluateTwoBodyInteractionCoefficient(int m1, int m2, int m3, int m4);
  
  // get Fourier transform of the interaction
  //
  // Q2_half = one half of q² value
  // retrun value =Fourier transform of ithe nteraction
  virtual double GetVofQ(double Q2_half);
 
  //evaluate the linearized index corresponding to a series of creation indices, momentum transfer, and modulo
  //
  // TmpIndices = array containing the indexes of the creation operators
  // g = modulo of the index on which the infinite sum will be performed
  // momentumTransfer = transfer of momentum between the creation and annihilation operators, in units of the number of flux quanta
  // return value = linearized index
  virtual int EvaluateLinearizedIndex (int* TmpIndices, int g, int momentumTransfer);
  
  
  // retrieve the separate indices starting from the linearized index
  //
  // linearizedIndex = linearized index
  // g = reference on the modulo of the index on which the infinite sum will be performed
  // momentumTransfer = reference on the transfer of momentum between the creation and annihilation operators, in units of the number of flux quanta
  // indices = reference on a temporary array to store the full index
  // nbrPermutations  = temporary array to store the numebr of permutations
  // return value = pointer to the array containing the creation indexes
  virtual void GetIndicesFromLinearizedIndex (int linearizedIndex, int& g, int& momentumTransfer, int* indices);

  // find the canonical index corresponding to a giving index
  //
  // index = index of the state whose canonical index has to be determined
  // sign = reference on an optional sign resulting when finiding the canonical index
  // indices = reference on a temporary array to store the full index
  // nbrPermutations  = temporary array to store the numebr of permutations
  //return value = canonical index
  virtual int GetCanonicalIndex (int index, int& sign, int* indices, int* nbrPermutations);
  
   // test if creation/annihilation are in canonical form
  //
  // mIndices = array that contains the creation indices (will be modified)
  // nIndices = array that contains the annihilation indices (will be modified)
  // linearizedMIndices = linearized creation index
  // linearizedNIndices = linearized annihilation index
  // totalMomentum = momentum sector of the creation/annihilation indices
  // canonicalMIndices = reference on the linearized creation index of the canonical form
  // canonicalMIndices = reference on the linearized annihilation index of the canonical form
  // canonicalTotalMomentum = reference on the momentum sector of the creation/annihilation indices of the canonical form
  // canonicalPermutationCoefficient = additional permutation coefficient to get the canonical form
  // return value = true if the indices are in canonical form
  virtual bool FindCanonicalIndices(int* mIndices, int* nIndices, int linearizedMIndices, int linearizedNIndices, int totalMomentum, int& canonicalMIndices, int& canonicalNIndices, int& canonicalTotalMomentum, double& canonicalPermutationCoefficient);  

};

// evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_n, assuming bosonic operators
//
// mIndices = array contaning the creation indices
// nIndices = array contaning the annihilation indices
// return value = numerical coefficient

inline double ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::EvaluateInteractionNIndexSymmetrizedCoefficient(int* mIndices, int* nIndices)
{  
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  FactorialCoefficient FactorialNBody;
  FactorialNBody.SetToOne();
  FactorialNBody.FactorialMultiply(this->NBodyValue);
  double FinalFactor = FactorialNBody.GetNumericalValue();
  for (int i = 0; i < this->NBodyValue; ++i)
    FinalFactor *= (M_PI * DoubleNbrLzValue);
  
  
  int* mIndices2 = new int[this->NBodyValue];
  int* nIndices2 = new int[this->NBodyValue];
  
  int momentumTransfer = 0;
  for (int i = 0; i < this->NBodyValue; ++i)
    momentumTransfer += (nIndices[i] - mIndices[i]);
  momentumTransfer /= this->NbrLzValue;
  int shiftedMomentumTransfer = momentumTransfer + this->NBodyValue - 1;
  int shift =  (this->NBodyValue - 1)*this->NBodyValue;
  
  double TmpInteraction; 
  int m1 = mIndices[0];
  int m2 = nIndices[0];
  
  int Factor = this->NbrLzValue;
  for (int k = 1; k < this->NBodyValue; ++k)
    {
      m1 += mIndices[k]*Factor;
      m2 += nIndices[k]*Factor;
      Factor *= this->NbrLzValue;
    }
  int m2Initial = m2;
  
  TmpInteraction = 0.0;
  for (int g = 0; g < this->NBodyValue; ++g)
    {
      TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][shift + g] * this->PrecalculatedInteractionCoefficients[m2][shiftedMomentumTransfer*this->NBodyValue + g];
    }
  
  for (int k = 0 ; k < this->NBodyValue; ++k)
    nIndices2[k] = nIndices[k];
  while (std::prev_permutation(nIndices2, nIndices2 + this->NBodyValue))
    {
      m2 = nIndices2[0];
      Factor = this->NbrLzValue;
      for (int k = 1 ; k < this->NBodyValue; ++k)
	{
	  m2 += nIndices2[k]*Factor;
	  Factor *= this->NbrLzValue;
	}
      
      for (int g = 0; g < this->NBodyValue; ++g)
	TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][shift + g] * this->PrecalculatedInteractionCoefficients[m2][shiftedMomentumTransfer*this->NBodyValue + g];
    }
  
  
  for (int k = 0 ; k < this->NBodyValue; ++k)
    mIndices2[k] = mIndices[k];
  while (std::prev_permutation(mIndices2, mIndices2 + this->NBodyValue))
    {
      m1 = mIndices2[0];
      Factor = this->NbrLzValue;
      for (int k = 1 ; k < this->NBodyValue; ++k)
	{
	  m1 += mIndices2[k]*Factor;
	  Factor *= this->NbrLzValue;
	}
      for (int g = 0; g < this->NBodyValue; ++g)
	TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][shift + g] * this->PrecalculatedInteractionCoefficients[m2][shiftedMomentumTransfer*this->NBodyValue + g];
      
      for (int k = 0 ; k < this->NBodyValue; ++k)
	nIndices2[k] = nIndices[k];
      while (std::prev_permutation(nIndices2, nIndices2 + this->NBodyValue))
	{
	  
	  m2 = nIndices2[0];
	  Factor = this->NbrLzValue;
	  for (int k = 1 ; k < this->NBodyValue; ++k)
	    {
	      m2 += nIndices2[k]*Factor;
	      Factor *= this->NbrLzValue;
	    }
	  
	  for (int g = 0; g < this->NBodyValue; ++g)
	    TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][shift + g] * this->PrecalculatedInteractionCoefficients[m2][shiftedMomentumTransfer*this->NBodyValue + g];
	}
    }    
  
  delete[] mIndices2;
  delete[] nIndices2;
  
  
  TmpInteraction /= FinalFactor;
  return TmpInteraction;
}

// evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_n, assuming fermionic operators
//
// mIndices = array contaning the creation indices
// nIndices = array contaning the annihilation indices
// return value = numerical coefficient

inline double ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::EvaluateInteractionNIndexAntiSymmetrizedCoefficient(int* mIndices, int* nIndices)
{  
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  FactorialCoefficient FactorialNBody;
  FactorialNBody.SetToOne();
  FactorialNBody.FactorialMultiply(this->NBodyValue);
  double FinalFactor = FactorialNBody.GetNumericalValue();
  for (int i = 0; i < this->NBodyValue; ++i)
    FinalFactor *= (M_PI * DoubleNbrLzValue);
  
  
  int* mIndices2 = new int[this->NBodyValue];
  int* nIndices2 = new int[this->NBodyValue];
  
  int momentumTransfer = 0;
  for (int i = 0; i < this->NBodyValue; ++i)
    momentumTransfer += (nIndices[i] - mIndices[i]);
  momentumTransfer /= this->NbrLzValue;
  int shiftedMomentumTransfer = momentumTransfer + this->NBodyValue - 1;
  int shift =  (this->NBodyValue - 1)*this->NBodyValue;
  
  double TmpInteraction; 
  int m1 = mIndices[0];
  int m2 = nIndices[0];
  
  int Factor = this->NbrLzValue;
  for (int k = 1; k < this->NBodyValue; ++k)
    {
      m1 += mIndices[k]*Factor;
      m2 += nIndices[k]*Factor;
      Factor *= this->NbrLzValue;
    }
  int m2Initial = m2;
  
  TmpInteraction = 0.0;
  for (int g = 0; g < this->NBodyValue; ++g)
    {
      TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][shift + g] * this->PrecalculatedInteractionCoefficients[m2][shiftedMomentumTransfer*this->NBodyValue + g];
    }
  
  int signature = 1;
  for (int k = 0 ; k < this->NBodyValue; ++k)
    nIndices2[k] = nIndices[k];
  while (std::prev_permutation(nIndices2, nIndices2 + this->NBodyValue))
    {
      signature *= -1;
      m2 = nIndices2[0];
      Factor = this->NbrLzValue;
      for (int k = 1 ; k < this->NBodyValue; ++k)
	{
	  m2 += nIndices2[k]*Factor;
	  Factor *= this->NbrLzValue;
	}
      
      for (int g = 0; g < this->NBodyValue; ++g)
	TmpInteraction += signature * this->PrecalculatedInteractionCoefficients[m1][shift + g] * this->PrecalculatedInteractionCoefficients[m2][shiftedMomentumTransfer*this->NBodyValue + g];
    }
  
  
  for (int k = 0 ; k < this->NBodyValue; ++k)
    mIndices2[k] = mIndices[k];
  signature = 1;
  while (std::prev_permutation(mIndices2, mIndices2 + this->NBodyValue))
    {
      signature *= -1;
      m1 = mIndices2[0];
      Factor = this->NbrLzValue;
      for (int k = 1 ; k < this->NBodyValue; ++k)
	{
	  m1 += mIndices2[k]*Factor;
	  Factor *= this->NbrLzValue;
	}
      for (int g = 0; g < this->NBodyValue; ++g)
	TmpInteraction += signature * this->PrecalculatedInteractionCoefficients[m1][shift + g] * this->PrecalculatedInteractionCoefficients[m2][shiftedMomentumTransfer*this->NBodyValue + g];
      
      for (int k = 0 ; k < this->NBodyValue; ++k)
	nIndices2[k] = nIndices[k];
      while (std::prev_permutation(nIndices2, nIndices2 + this->NBodyValue))
	{
	  signature *= -1;
	  m2 = nIndices2[0];
	  Factor = this->NbrLzValue;
	  for (int k = 1 ; k < this->NBodyValue; ++k)
	    {
	      m2 += nIndices2[k]*Factor;
	      Factor *= this->NbrLzValue;
	    }
	  
	  for (int g = 0; g < this->NBodyValue; ++g)
	    TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][shift + g] * this->PrecalculatedInteractionCoefficients[m2][shiftedMomentumTransfer*this->NBodyValue + g];
	}
    }    
  
  delete[] mIndices2;
  delete[] nIndices2;
  
  
  TmpInteraction /= FinalFactor;
  return TmpInteraction;
}


//evaluate the linearized index corresponding to a series of creation indices, momentum transfer, and modulo
//
// TmpIndices = array containing the indexes of the creation operators
// g = modulo of the index on which the infinite sum will be performed
// momentumTransfer = transfer of momentum between the creation and annihilation operators, in units of the number of flux quanta
// return value = linearized index

inline int ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::EvaluateLinearizedIndex (int* TmpIndices, int g, int momentumTransfer)
{
  int LinearizedIndex = (g + this->NBodyValue * (this->NBodyValue - 1 + momentumTransfer));
  for (int i = this->NBodyValue - 1; i >= 0; --i)
    {
      LinearizedIndex *= this->NbrLzValue;
      LinearizedIndex += TmpIndices[i];
    }
  return LinearizedIndex;  
}

// retrieve the separate indices starting from the linearized index
//
// linearizedIndex = linearized index
// g = reference on the modulo of the index on which the infinite sum will be performed
// momentumTransfer = reference on the transfer of momentum between the creation and annihilation operators, in units of the number of flux quanta
// indices = temporary array to store the full index
// return value = pointer to the array containing the creation indexes

inline void ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::GetIndicesFromLinearizedIndex (int linearizedIndex, int& g, int& momentumTransfer, int* indices)
{
  int TmpCoefficient = pow(this->NbrLzValue, this->NBodyValue) * this->NBodyValue;
  int TmpCoefficient1;
  momentumTransfer = linearizedIndex / TmpCoefficient;
  linearizedIndex -= momentumTransfer * TmpCoefficient;
  momentumTransfer -= (this->NBodyValue - 1);
  TmpCoefficient /= this->NBodyValue;
  g = linearizedIndex / TmpCoefficient;
  linearizedIndex -= g*TmpCoefficient;
  TmpCoefficient /= this->NbrLzValue;
  
  int i = this->NBodyValue - 1;
  for (; i > 0; --i)
  {
    indices[i] = linearizedIndex / TmpCoefficient;
    linearizedIndex -= indices[i] * TmpCoefficient;
    TmpCoefficient /= this->NbrLzValue;
  }
  indices[0] = linearizedIndex;
}

#endif
