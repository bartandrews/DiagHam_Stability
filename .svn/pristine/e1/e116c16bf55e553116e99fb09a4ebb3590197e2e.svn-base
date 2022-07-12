///////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of abstract fractional quantum Hall hamiltonian          //
//              associated to particles with SU(4) spin on a sphere           //
//                                                                            //
//                        last modification : 28/11/2006                      //
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


#ifndef ABSTRACTQHEONSPHEREWITHSU4SPINCASIMIRHAMILTONIAN_H
#define ABSTRACTQHEONSPHEREWITHSU4SPINCASIMIRHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSU4Spin.h"
#include "Hamiltonian/AbstractQHEOnSphereHamiltonian.h"

class ParticleOnSphereWithSU4SpinL2Hamiltonian;
class ParticleOnSphereWithSU4SpinS2Hamiltonian;

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian : public AbstractQHEOnSphereHamiltonian
{

 protected:
  
  // simple method used to store indices attached to each interaction factor
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

  // arrays storing interaction factors for two-body interactions
  // diagonal in spin and pseudospin
  double *M12InteractionFactorsupup;
  double *M12InteractionFactorsumum;
  double *M12InteractionFactorsdpdp;
  double *M12InteractionFactorsdmdm;
  // crossterm between spin and pseudospin
  double *M12InteractionFactorsupum;
  double *M12InteractionFactorsupdp;
  double *M12InteractionFactorsupdm;
  double *M12InteractionFactorsumdp;
  double *M12InteractionFactorsumdm;
  double *M12InteractionFactorsdpdm;
  // entangled interaction terms
  double *M12InteractionFactorsupdmEnt;
  double *M12InteractionFactorsumdpEnt;

  // array that contains all one-body interaction factors for particles with spin up/isospin plus
  double* OneBodyInteractionFactorsupup;
  // array that contains all one-body interaction factors for particles with spin up/isospin minus
  double* OneBodyInteractionFactorsumum;
  // array that contains all one-body interaction factors for particles with spin down/isospin plus
  double* OneBodyInteractionFactorsdpdp;
  // array that contains all one-body interaction factors for particles with spin down/isospin minus
  double* OneBodyInteractionFactorsdmdm;

  // pointer to an optional L^2 operator in the Hamiltonian 
  ParticleOnSphereWithSU4SpinL2Hamiltonian* L2Hamiltonian;
  // pointer to an optional S^2 operator in the Hamiltonian 
  ParticleOnSphereWithSU4SpinS2Hamiltonian* S2Hamiltonian;


 public:

  // default constructor
  //
  AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian();

  // destructor
  //
  virtual ~AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian();

  // ask if Hamiltonian implements methods using hermitian symmetry 
  //
  virtual bool IsHermitian();

  // ask if Hamiltonian implements methods applying the conjugate of the Hamiltonian
  //
  virtual bool IsConjugate();

  // add an additional S^2 term to the Hamiltonian
  //
  // totalLz = twice the projected momentum total value
  // totalSz = twice the projected spin total value
  // factor = factor in front of the S^2
  // memory = amount of memory that can be used for S^2  precalculations 
  void AddS2 (int totalLz, int totalSz, double factor, long memory);


  // add an additional S^2 term to the Hamiltonian
  //
  // totalLz = twice the projected momentum total value
  // totalSz = twice the projected spin total value
  // factor = factor in front of the S^2
  // memory = amount of memory that can be used for S^2  precalculations 
  void AddL2 (int totalLz, int totalSz, double factor, long memory);
  
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
  virtual RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						  int firstComponent, int nbrComponent);

 protected:
  
  // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int index, RealVector& vSource, RealVector& vDestination);

  // core part of the AddMultiply method involving the two-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int index, RealVector* vSources, 
						     RealVector* vDestinations, int nbrVectors, double* tmpCoefficients);

  // core part of the FastMultiplication method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSU4Spin* particles, int index, 
							    int* indexArray, double* coefficientArray, long& position);

  // core part of the FastMultiplication method involving the one-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSU4Spin* particles, int index, 
							    int* indexArray, double* coefficientArray, long& position);

  // core part of the AddMultiply method involving the one-body interaction, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int firstComponent, int lastComponent,
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
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int firstComponent, int lastComponent,
						     int step, RealVector* vSources, RealVector* vDestinations, int nbrVectors);

  // core part of the PartialFastMultiplicationMemory method involving two-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSU4Spin* particles, int firstComponent, int lastComponent, long& memory);
 
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

inline void AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int index, RealVector& vSource, RealVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();

  double Coefficient, Coefficient2;
  int SumIndices;
  int TmpNbrM3Values;
  int* TmpM3Values;
  int ReducedNbrInteractionFactors = 0;
  int Index;
  if (this->M12InteractionFactorsupup!=NULL)
    {
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  Coefficient = particles->AupAup(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);
	  if (Coefficient != 0.0)
	    {
	      Coefficient *= vSource[index];
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAdup(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsupup[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
    }
  if (this->M12InteractionFactorsumum!=NULL)
    {
      ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  Coefficient = particles->AumAum(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);
	  if (Coefficient != 0.0)
	    {
	      Coefficient *= vSource[index];
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdumAdum(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsumum[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
    }
  if (this->M12InteractionFactorsdpdp!=NULL)
    {
      ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  Coefficient = particles->AdpAdp(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      Coefficient *= vSource[index];
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AddpAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsdpdp[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
    }
  if (this->M12InteractionFactorsdmdm!=NULL)
    {
      ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{ 
	  Coefficient = particles->AdmAdm(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      Coefficient *= vSource[index];
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AddmAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsdmdm[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
    }
  if (this->M12InteractionFactorsupum!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AupAum(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      Coefficient *= vSource[index];
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAdum(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsupum[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsupdp!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AupAdp(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      Coefficient *= vSource[index];		  
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsupdp[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsupdm!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AupAdm(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      Coefficient *= vSource[index];
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsupdm[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsumdp!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AumAdp(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      Coefficient *= vSource[index];
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdumAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsumdp[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsumdm!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AumAdm(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      Coefficient *= vSource[index];
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdumAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsumdm[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsdpdm!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AdpAdm(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      Coefficient *= vSource[index];
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AddpAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsdpdm[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}      
    }
  if (this->M12InteractionFactorsupdmEnt!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AupAdm(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      Coefficient *= vSource[index];
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdumAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsupdmEnt[ReducedNbrInteractionFactors] * Coefficient2;
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
    if (this->M12InteractionFactorsumdpEnt!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AumAdp(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      Coefficient *= vSource[index];
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * this->M12InteractionFactorsumdpEnt[ReducedNbrInteractionFactors] * Coefficient2;
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

inline void AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int index, RealVector* vSources, 
												    RealVector* vDestinations, int nbrVectors, double* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();

  double Coefficient, Coefficient2;
  int Index, SumIndices;
  int TmpNbrM3Values;
  int* TmpM3Values;
  int ReducedNbrInteractionFactors = 0;
  if (this->M12InteractionFactorsupup!=NULL)
    {
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  Coefficient = particles->AupAup(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);
	  if (Coefficient != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient * vSources[p][index];
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAdup(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      Coefficient2 *= this->M12InteractionFactorsupup[ReducedNbrInteractionFactors];
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += tmpCoefficients[l] * Coefficient2;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
    }
  if (this->M12InteractionFactorsumum!=NULL)
    {
      ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  Coefficient = particles->AumAum(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);
	  if (Coefficient != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient * vSources[p][index];
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdumAdum(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      Coefficient2 *= this->M12InteractionFactorsumum[ReducedNbrInteractionFactors];
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += tmpCoefficients[l] * Coefficient2;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
    }
  if (this->M12InteractionFactorsdpdp!=NULL)
    {
      ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  Coefficient = particles->AdpAdp(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient * vSources[p][index];
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AddpAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      Coefficient2 *= this->M12InteractionFactorsdpdp[ReducedNbrInteractionFactors];
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += tmpCoefficients[l] * Coefficient2;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
    }
  if (this->M12InteractionFactorsdmdm!=NULL)
    {
      ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{ 
	  Coefficient = particles->AdmAdm(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient * vSources[p][index];
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AddmAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      Coefficient2 *= this->M12InteractionFactorsdmdm[ReducedNbrInteractionFactors];
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += tmpCoefficients[l] * Coefficient2;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
    }
  if (this->M12InteractionFactorsupum!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AupAum(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient * vSources[p][index];
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAdum(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      Coefficient2 *= this->M12InteractionFactorsupum[ReducedNbrInteractionFactors];
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += tmpCoefficients[l] * Coefficient2;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsupdp!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AupAdp(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient * vSources[p][index];		  
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      Coefficient2 *= this->M12InteractionFactorsupdp[ReducedNbrInteractionFactors];
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += tmpCoefficients[l] * Coefficient2;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsupdm!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AupAdm(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient * vSources[p][index];
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      Coefficient2 *= this->M12InteractionFactorsupdm[ReducedNbrInteractionFactors];
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += tmpCoefficients[l] * Coefficient2;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsumdp!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AumAdp(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient * vSources[p][index];
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdumAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      Coefficient2 *= this->M12InteractionFactorsumdp[ReducedNbrInteractionFactors];
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += tmpCoefficients[l] * Coefficient2;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsumdm!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AumAdm(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient * vSources[p][index];
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdumAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      Coefficient2 *= this->M12InteractionFactorsumdm[ReducedNbrInteractionFactors];
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += tmpCoefficients[l] * Coefficient2;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsdpdm!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AdpAdm(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient * vSources[p][index];
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AddpAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      Coefficient2 *= this->M12InteractionFactorsdpdm[ReducedNbrInteractionFactors];
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += tmpCoefficients[l] * Coefficient2;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsupdmEnt!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AupAdm(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient * vSources[p][index];
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdumAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      Coefficient2 *= this->M12InteractionFactorsupdmEnt[ReducedNbrInteractionFactors];
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += tmpCoefficients[l] * Coefficient2;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsumdpEnt!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AumAdp(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      for (int p = 0; p < nbrVectors; ++p)
		tmpCoefficients[p] = Coefficient * vSources[p][index];
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      Coefficient2 *= this->M12InteractionFactorsumdpEnt[ReducedNbrInteractionFactors];
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += tmpCoefficients[l] * Coefficient2;
		    }
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

inline void AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSU4Spin* particles, int index, 
													   int* indexArray, double* coefficientArray, long& position)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int Dim = particles->GetHilbertSpaceDimension();
  int SumIndices;
  int TmpNbrM3Values;
  int* TmpM3Values;
  int ReducedNbrInteractionFactors = 0;
  if (this->M12InteractionFactorsupup!=NULL)
    {
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  Coefficient = particles->AupAup(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAdup(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsupup[ReducedNbrInteractionFactors];
		      ++position;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
    }
  if (this->M12InteractionFactorsumum!=NULL)
    {
      ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  Coefficient = particles->AumAum(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdumAdum(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsumum[ReducedNbrInteractionFactors];
		      ++position;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
    }
  if (this->M12InteractionFactorsdpdp!=NULL)
    {
      ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  Coefficient = particles->AdpAdp(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AddpAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsdpdp[ReducedNbrInteractionFactors];
		      ++position;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
    }
  if (this->M12InteractionFactorsdmdm!=NULL)
    {
      ReducedNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{ 
	  Coefficient = particles->AdmAdm(index, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
	      TmpM3Values = this->M3IntraValues[m1];
	      TmpNbrM3Values = this->NbrM3IntraValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AddmAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsdmdm[ReducedNbrInteractionFactors];
		      ++position;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
	}
    }
  if (this->M12InteractionFactorsupum!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AupAum(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAdum(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsupum[ReducedNbrInteractionFactors];
		      ++position;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsupdp!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AupAdp(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsupdp[ReducedNbrInteractionFactors];
		      ++position;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsupdm!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AupAdm(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsupdm[ReducedNbrInteractionFactors];
		      ++position;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsumdp!=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AumAdp(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdumAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsumdp[ReducedNbrInteractionFactors];
		      ++position;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsumdm !=NULL)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AumAdm(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdumAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsumdm[ReducedNbrInteractionFactors];
		      ++position;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsdpdm!=0)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AdpAdm(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AddpAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsdpdm[ReducedNbrInteractionFactors];
		      ++position;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsupdmEnt!=0)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AupAdm(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdumAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsupdmEnt[ReducedNbrInteractionFactors];
		      ++position;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
	}
    }
  if (this->M12InteractionFactorsumdpEnt!=0)
    {
      ReducedNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  Coefficient = particles->AumAdp(index, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
	      TmpM3Values = this->M3InterValues[m1];
	      TmpNbrM3Values = this->NbrM3InterValues[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdupAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * this->M12InteractionFactorsumdpEnt[ReducedNbrInteractionFactors];
		      ++position;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
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

inline void AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int firstComponent, int lastComponent,
												    int step, RealVector& vSource, RealVector& vDestination)
{
  if ((this->OneBodyInteractionFactorsupup != 0)&&(this->OneBodyInteractionFactorsumum != 0)&&(this->OneBodyInteractionFactorsdpdp != 0) && (this->OneBodyInteractionFactorsdmdm != 0))
    {
      double TmpDiagonal = 0.0;
      for (int i = firstComponent; i < lastComponent; i += step)
	{ 
	  TmpDiagonal = 0.0;
	  for (int j = 0; j <= this->LzMax; ++j) 
	    {
	      TmpDiagonal += this->OneBodyInteractionFactorsupup[j] * particles->AdupAup(i, j);
	      TmpDiagonal += this->OneBodyInteractionFactorsumum[j] * particles->AdumAum(i, j);
	      TmpDiagonal += this->OneBodyInteractionFactorsdpdp[j] * particles->AddpAdp(i, j);
	      TmpDiagonal += this->OneBodyInteractionFactorsdmdm[j] * particles->AddmAdm(i, j);
	    }
	  vDestination[i] += (this->HamiltonianShift + TmpDiagonal) * vSource[i];
	}
    }
  else
    {
      if ((this->OneBodyInteractionFactorsupup == 0)&&(this->OneBodyInteractionFactorsumum == 0)&&(this->OneBodyInteractionFactorsdpdp == 0) && (this->OneBodyInteractionFactorsdmdm == 0))
	{
	  for (int i = firstComponent; i < lastComponent; i += step)
	    vDestination[i] += this->HamiltonianShift * vSource[i];
	}
      else
	{
	  cout << "Cases with partial one-body interactions not coded, yet in AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::EvaluateMNOneBodyAddMultiplyComponent"<<endl;
	  exit(1);
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

inline void AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSU4Spin* particles, int firstComponent, int lastComponent,
												    int step, RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  if ((this->OneBodyInteractionFactorsupup != 0)&&(this->OneBodyInteractionFactorsumum != 0)&&(this->OneBodyInteractionFactorsdpdp != 0) && (this->OneBodyInteractionFactorsdmdm != 0))
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
		  TmpDiagonal += this->OneBodyInteractionFactorsupup[j] * particles->AdupAup(i, j);
		  TmpDiagonal += this->OneBodyInteractionFactorsumum[j] * particles->AdumAum(i, j);
		  TmpDiagonal += this->OneBodyInteractionFactorsdpdp[j] * particles->AddpAdp(i, j);
		  TmpDiagonal += this->OneBodyInteractionFactorsdmdm[j] * particles->AddmAdm(i, j);
		}
	      TmpDestinationVector[i] += (this->HamiltonianShift + TmpDiagonal)* TmpSourceVector[i];
	    }
	}
    }
  else
    {
      if ((this->OneBodyInteractionFactorsupup == 0)&&(this->OneBodyInteractionFactorsumum == 0)&&(this->OneBodyInteractionFactorsdpdp == 0) && (this->OneBodyInteractionFactorsdmdm == 0))
	{
	  for (int p = 0; p < nbrVectors; ++p)
	    {
	      RealVector& TmpSourceVector = vSources[p];
	      RealVector& TmpDestinationVector = vDestinations[p];
	      for (int i = firstComponent; i < lastComponent; i += step)
		TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	    }
	}
      else
	{
	  cout << "Cases with partial one-body interactions not coded, yet in AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::EvaluateMNOneBodyAddMultiplyComponent"<<endl;
	  exit(1);
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

inline void AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSU4Spin* particles, int index, 
													   int* indexArray, double* coefficientArray, long& position)
{
  if ((this->OneBodyInteractionFactorsupup != 0)||(this->OneBodyInteractionFactorsumum != 0)||(this->OneBodyInteractionFactorsdpdp != 0) || (this->OneBodyInteractionFactorsdmdm != 0))
    {
      double TmpDiagonal = 0.0;
      for (int j = 0; j <= this->LzMax; ++j) 
	{
	  if (this->OneBodyInteractionFactorsupup != 0)
	    TmpDiagonal += this->OneBodyInteractionFactorsupup[j] * particles->AdupAup(index, j);
	  if (this->OneBodyInteractionFactorsumum != 0)
	    TmpDiagonal += this->OneBodyInteractionFactorsumum[j] * particles->AdumAum(index, j);
	  if (this->OneBodyInteractionFactorsdpdp != 0)
	    TmpDiagonal += this->OneBodyInteractionFactorsdpdp[j] * particles->AddpAdp(index, j);
	  if (this->OneBodyInteractionFactorsdmdm != 0)
	    TmpDiagonal += this->OneBodyInteractionFactorsdmdm[j] * particles->AddmAdm(index, j);
	}
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

inline void AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSU4Spin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  // double* TmpInteractionFactor;
  int Dim = particles->GetHilbertSpaceDimension();
  int SumIndices;
  int TmpNbrM3Values;
  int* TmpM3Values;
  
  for (int i = firstComponent; i < lastComponent; ++i)
    {
      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	{
	  if (this->M12InteractionFactorsupup!=0)
	    {
	      Coefficient = particles->AupAup(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		  TmpM3Values = this->M3IntraValues[m1];
		  TmpNbrM3Values = this->NbrM3IntraValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = particles->AdupAdup(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		}
	    }
	  if (this->M12InteractionFactorsumum!=0)
	    {
	      Coefficient = particles->AumAum(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		  TmpM3Values = this->M3IntraValues[m1];
		  TmpNbrM3Values = this->NbrM3IntraValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = particles->AdumAdum(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		}
	    }
	  if (this->M12InteractionFactorsdpdp!=0)
	    {
	      Coefficient = particles->AdpAdp(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		  TmpM3Values = this->M3IntraValues[m1];
		  TmpNbrM3Values = this->NbrM3IntraValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = particles->AddpAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		}
	    }
	  if (this->M12InteractionFactorsdmdm!=0)
	    {
	      Coefficient = particles->AdmAdm(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		  TmpM3Values = this->M3IntraValues[m1];
		  TmpNbrM3Values = this->NbrM3IntraValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = particles->AddmAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		}
	    }
	}
      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	{
	  if (this->M12InteractionFactorsupum!=0)
	    {
	      Coefficient = particles->AupAum(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = particles->AdupAdum(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		}
	    }
	  if (this->M12InteractionFactorsupdp!=0)
	    {
	      Coefficient = particles->AupAdp(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = particles->AdupAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		}
	    }
	  if (this->M12InteractionFactorsupdm!=0)
	    {
	      Coefficient = particles->AupAdm(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = particles->AdupAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		}
	    }
	  if (this->M12InteractionFactorsumdp!=0)
	    {
	      Coefficient = particles->AumAdp(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      //cout << "Adum("<<TmpM3Values[m3]<<")Addp("<<SumIndices - TmpM3Values[m3]<<")Aum("<<this->M1InterValue[m1]<<")Adp("<< this->M2InterValue[m1]<<")"<<endl;
		      Index = particles->AdumAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		}
	    }
	  if (this->M12InteractionFactorsumdm!=0)
	    {
	      Coefficient = particles->AumAdm(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = particles->AdumAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		}
	    }
	  if (this->M12InteractionFactorsdpdm!=0)
	    {
	      Coefficient = particles->AdpAdm(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = particles->AddpAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		}
	    }
	  if (this->M12InteractionFactorsupdmEnt!=0)
	    {
	      Coefficient = particles->AupAdm(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = particles->AdumAddp(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }
		}
	    }
	  if (this->M12InteractionFactorsumdpEnt!=0)
	    {
	      Coefficient = particles->AumAdp(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = particles->AdupAddm(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
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
  if ((this->OneBodyInteractionFactorsupup != 0)||(this->OneBodyInteractionFactorsumum != 0)||(this->OneBodyInteractionFactorsdpdp != 0) || (this->OneBodyInteractionFactorsdmdm != 0))
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  ++memory;
	  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
	}	  
    }
}




#endif
