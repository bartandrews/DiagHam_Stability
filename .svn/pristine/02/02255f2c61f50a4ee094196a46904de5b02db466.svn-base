////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//                class of quantum Hall Hamiltonian associated                //
//                          to particles on a lattice                         //
//                                                                            //
//                        last modification : 12/02/2008                      //
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


#include "config.h"
#include "Hamiltonian/AbstractQHEOnLatticeHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>
#include <fstream>


using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;
using std::ostream;


// default constructor
//
AbstractQHEOnLatticeHamiltonian::AbstractQHEOnLatticeHamiltonian()
{
  this->NbrQ12Indices=0;
}

// destructor
//

AbstractQHEOnLatticeHamiltonian::~AbstractQHEOnLatticeHamiltonian()
{
  if (FastMultiplicationFlag)
    {
      delete [] this->NbrRealInteractionPerComponent;
      delete [] this->NbrComplexInteractionPerComponent;
      long MinIndex, MaxIndex;
      this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
      int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
      for (int i=0; i<EffectiveHilbertSpaceDimension; ++i)
	{
	  delete [] this->InteractionPerComponentIndex[i];
	  delete [] this->InteractionPerComponentCoefficientIndex[i];
	}
      delete [] this->InteractionPerComponentIndex;
      delete [] this->InteractionPerComponentCoefficientIndex;
      this->RealInteractionCoefficients.Empty();
      this->ComplexInteractionCoefficients.Empty();
      this->FastMultiplicationFlag=false;
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void AbstractQHEOnLatticeHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
    if (NbrHoppingTerms>0)
    {
      delete [] this->HoppingTerms;
      delete [] this->KineticQi;
      delete [] this->KineticQf;
    }
  if (NbrInteractionFactors>0)
    {
      delete [] this->InteractionFactors;
      delete [] this->Q1Value;
      delete [] this->Q2Value;
      delete [] this->Q3Value;
      delete [] this->Q4Value;
    }
  if (NbrQ12Indices>0)
    {
      for (int i=0; i<NbrQ12Indices; ++i)
	{
	  delete [] this->Q3PerQ12[i];
	  delete [] this->Q4PerQ12[i];
	}
      delete [] this->NbrQ34Values;
      delete [] this->InteractionFactors;
      delete [] this->Q1Value;
      delete [] this->Q2Value;
      
    }
  if (NbrDiagonalInteractionFactors>0)
    {
      delete [] this->DiagonalInteractionFactors;
      delete [] this->DiagonalQValues;
    }
  this->Particles = (ParticleOnLattice*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// set flux density in units of flux quanta through the lattice
//
// nbrFluxQuanta = flux quantua piercing the lattice
void AbstractQHEOnLatticeHamiltonian::SetNbrFluxQuanta(int nbrFluxQuanta)
{
  this->NbrFluxQuanta=nbrFluxQuanta;
  this->FluxDensity=((double)nbrFluxQuanta)/NbrCells;
  this->Particles->SetNbrFluxQuanta(nbrFluxQuanta);
  if (NbrHoppingTerms>0)
    {
      delete [] this->HoppingTerms;
      delete [] this->KineticQi;
      delete [] this->KineticQf;
    }
  if (NbrInteractionFactors>0)
    {
      delete [] this->InteractionFactors;
      delete [] this->Q1Value;
      delete [] this->Q2Value;
      delete [] this->Q3Value;
      delete [] this->Q4Value;
    }
  if (NbrQ12Indices>0)
    {
      for (int i=0; i<NbrQ12Indices; ++i)
	{
	  delete [] this->Q3PerQ12[i];
	  delete [] this->Q4PerQ12[i];
	}
      delete [] this->NbrQ34Values;
      delete [] this->InteractionFactors;
      delete [] this->Q1Value;
      delete [] this->Q2Value;
      
    }
  if (NbrDiagonalInteractionFactors>0)
    {
      delete [] this->DiagonalInteractionFactors;
      delete [] this->DiagonalQValues;
    }  
  bool EnableFastCalculation=FastMultiplicationFlag;
  if (FastMultiplicationFlag)
    {
      delete [] this->NbrRealInteractionPerComponent;
      delete [] this->NbrComplexInteractionPerComponent;
      long MinIndex, MaxIndex;
      this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
      int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
      for (int i=0; i<EffectiveHilbertSpaceDimension; ++i)
	{
	  delete [] this->InteractionPerComponentIndex[i];
	  delete [] this->InteractionPerComponentCoefficientIndex[i];
	}
      delete [] this->InteractionPerComponentIndex;
      delete [] this->InteractionPerComponentCoefficientIndex;
      this->RealInteractionCoefficients.Empty();
      this->ComplexInteractionCoefficients.Empty();
      this->FastMultiplicationFlag=false;
    }    
  this->EvaluateInteractionFactors();
  if (this->LoadedPrecalculation)
    {
      cout << "Cannot re-enable fast calculation when precalculation data was loaded from file!"<<endl;
      cout << "Reverting to slow calculation"<<endl;
    }
  else if (EnableFastCalculation)
    {
      int TmpMemory = this->FastMultiplicationMemory(0);
      if (TmpMemory < 1024)
	cout  << "fast = " <<  TmpMemory << "b ";
      else
	if (TmpMemory < (1 << 20))
	  cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	else
	  if (TmpMemory < (1 << 30))
	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	  else
	    cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
      cout << endl;
      if (AllowedMemory > 0)
	{
	  this->EnableFastMultiplication();
	}
    }
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* AbstractQHEOnLatticeHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int AbstractQHEOnLatticeHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void AbstractQHEOnLatticeHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractQHEOnLatticeHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  int dim = this->Particles->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractQHEOnLatticeHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& AbstractQHEOnLatticeHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnLatticeHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
							       int firstComponent, int nbrComponent) 
{
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnLatticeHamiltonian::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
								     int firstComponent, int nbrComponent)
{
  return vDestinations;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return vDestination;
}

// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored
RealVector* AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						  int firstComponent, int nbrComponent)
{
  return vDestinations;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiply(ComplexVector& vSource,
								 ComplexVector& vDestination, 
								int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();  
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      // deal with kinetic energy terms first!      
      int qi;
      int qf;
      double TmpInteractionRe,TmpInteractionIm;
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
      for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	{
	  qi = this->KineticQi[j];
	  qf = this->KineticQf[j];
	  TmpInteractionRe = this->HoppingTerms[j].Re;
	  TmpInteractionIm = this->HoppingTerms[j].Im;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      //cout << "element "<<qi<<"->"<<qf<<" on "<<i<<": "; 
	      Index = TmpParticles->AdA(i, qf, qi, Coefficient);
	      //cout << "target: "<<Index<<endl;
	      if (Index < Dim)
		{
		  vDestination.Re(Index) += Coefficient * (TmpInteractionRe*vSource[i].Re - TmpInteractionIm*vSource[i].Im);
		  vDestination.Im(Index) += Coefficient * (TmpInteractionRe*vSource[i].Im + TmpInteractionIm*vSource[i].Re);
		}
	    }
	}
      qi = this->KineticQi[ReducedNbrHoppingTerms];
      qf = this->KineticQf[ReducedNbrHoppingTerms];
      TmpInteractionRe = this->HoppingTerms[ReducedNbrHoppingTerms].Re;
      TmpInteractionIm = this->HoppingTerms[ReducedNbrHoppingTerms].Im;
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  //cout << "element "<<qi<<"->"<<qf<<" on "<<i<<": "; 
	  Index = TmpParticles->AdA(i, qf, qi, Coefficient);
	  //cout << "target: "<<Index<<endl;
	  if (Index < Dim)
	    {
	      vDestination.Re(Index) += Coefficient * (TmpInteractionRe*vSource[i].Re - TmpInteractionIm*vSource[i].Im);
	      vDestination.Im(Index) += Coefficient * (TmpInteractionRe*vSource[i].Im + TmpInteractionIm*vSource[i].Re);
	    }
	  vDestination.Re(i) += this->HamiltonianShift * vSource[i].Re;
	  vDestination.Im(i) += this->HamiltonianShift * vSource[i].Im;
	}         
      
      // four-fermion interactions:
      if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
	      TmpInteractionRe = this->InteractionFactors[j].Re;
	      TmpInteractionIm = this->InteractionFactors[j].Im;
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      vDestination.Re(Index) += Coefficient * (TmpInteractionRe*vSource[i].Re - TmpInteractionIm*vSource[i].Im);
		      vDestination.Im(Index) += Coefficient * (TmpInteractionRe*vSource[i].Im + TmpInteractionIm*vSource[i].Re);
		    }
		}
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2, TmpRe, TmpIm;
	  int ProcessedNbrInteractionFactors;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      ProcessedNbrInteractionFactors = 0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = TmpParticles->AA(i, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      TmpRe = vSource[i].Re*Coefficient;
		      TmpIm = vSource[i].Im*Coefficient;
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = TmpParticles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
			      TmpInteractionRe = this->InteractionFactors[ProcessedNbrInteractionFactors].Re;
			      TmpInteractionIm = this->InteractionFactors[ProcessedNbrInteractionFactors].Im;
			      vDestination.Re(Index) += Coefficient2 * (TmpRe*TmpInteractionRe-TmpIm*TmpInteractionIm);
			      vDestination.Im(Index) += Coefficient2 * (TmpRe*TmpInteractionIm+TmpIm*TmpInteractionRe);
			    }
			  ++ProcessedNbrInteractionFactors;
			}
		    }
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
	    }
	}	  

      // separated diagonal terms as these will be the general rule for contact interactions
      if (NbrDiagonalInteractionFactors>0)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient = TmpParticles->AdAdAADiagonal(i, NbrDiagonalInteractionFactors,
							 DiagonalInteractionFactors, DiagonalQValues);
	      vDestination.Re(i) +=  Coefficient * vSource[i].Re;
	      vDestination.Im(i) +=  Coefficient * vSource[i].Im;
	    }
	}
      
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  unsigned short* TmpCoefficientIndexArray;
	  double TmpRe, TmpIm;
	  unsigned short TmpNbrRealInteraction;
	  unsigned short TmpNbrComplexInteraction;
	  Complex *TmpCPtr;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
	      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
	      TmpRe = vSource[k].Re;
	      TmpIm = vSource[k].Im;
	      int Pos=0;
	      for (; Pos < TmpNbrRealInteraction; ++Pos)
		{
		  vDestination.Re(TmpIndexArray[Pos]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpRe;
		  vDestination.Im(TmpIndexArray[Pos]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpIm;
		}
	      for (int j=0; j < TmpNbrComplexInteraction; ++j)
		{
		  TmpCPtr= &(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]]);
		  vDestination.Re(TmpIndexArray[Pos]) +=  TmpCPtr->Re*TmpRe-TmpCPtr->Im*TmpIm;
		  vDestination.Im(TmpIndexArray[Pos]) +=  TmpCPtr->Re*TmpIm+TmpCPtr->Im*TmpRe;
		  ++Pos;
		}
	      vDestination.Re(k) += this->HamiltonianShift * TmpRe;
	      vDestination.Im(k) += this->HamiltonianShift * TmpIm;	      
	      ++k;
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }

  //cout << "vDestination:" <<endl<<vDestination<<endl;
  
  return vDestination;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										   int firstComponent, int nbrComponent)
{

  cout << "AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiplyPartialFastMultiply must be defined"<<endl;
  exit(-1);
//   int LastComponent = firstComponent + nbrComponent;
//   int Dim = this->Particles->GetHilbertSpaceDimension();
//   double Coefficient;
//   ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
//   int* TmpIndexArray;
//   double* TmpCoefficientArray; 
//   int j;
//   int TmpNbrInteraction;
//   firstComponent -= this->PrecalculationShift;
//   LastComponent -= this->PrecalculationShift;
//   int Pos = firstComponent / this->FastMultiplicationStep; 
//   int PosMod = firstComponent % this->FastMultiplicationStep;
//   if (PosMod != 0)
//     {
//       ++Pos;
//       PosMod = this->FastMultiplicationStep - PosMod;
//     }
//   int l =  PosMod + firstComponent + this->PrecalculationShift;
//   for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
//     {
//       TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
//       TmpIndexArray = this->InteractionPerComponentIndex[Pos];
//       TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
//       Coefficient = vSource[l];
//       for (j = 0; j < TmpNbrInteraction; ++j)
// 	vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
//       vDestination[l] += this->HamiltonianShift * Coefficient;
//       l += this->FastMultiplicationStep;
//       ++Pos;
//     }
//   int Index;
//   int m1;
//   int m2;
//   int m3;
//   int m4;
//   double TmpInteraction;
//   int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
//   firstComponent += this->PrecalculationShift;
//   LastComponent += this->PrecalculationShift;
//   for (int k = 0; k < this->FastMultiplicationStep; ++k)
//     if (PosMod != k)
//       {		
// 	if (this->NbrM12Indices == 0)
// 	  {
// 	    for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
// 	      {
// 		m1 = this->M1Value[j];
// 		m2 = this->M2Value[j];
// 		m3 = this->M3Value[j];
// 		TmpInteraction = this->InteractionFactors[j];
// 		m4 = m1 + m2 - m3;
// 		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
// 		  {
// 		    Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
// 		    if (Index < Dim)
// 		      vDestination.Re(Index) += Real( Coefficient * TmpInteraction * vSource[i];
// 		  }
// 	      }
// 	    m1 = this->M1Value[ReducedNbrInteractionFactors];
// 	    m2 = this->M2Value[ReducedNbrInteractionFactors];
// 	    m3 = this->M3Value[ReducedNbrInteractionFactors];
// 	    TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
// 	    m4 = m1 + m2 - m3;
// 	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
// 	      {
// 		Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
// 		if (Index < Dim)
// 		  vDestination.Re(Index) += Real( Coefficient * TmpInteraction * vSource[i];
// 		vDestination[i] += this->HamiltonianShift * vSource[i];
// 	      }
// 	  }
// 	else
// 	  {
// 	    double Coefficient2;
// 	    int SumIndices;
// 	    int TmpNbrQ34Values;
// 	    int* TmpM3Values;
// 	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
// 	      {
// 		ReducedNbrInteractionFactors = 0;
// 		for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
// 		  {
// 		    Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
// 		    if (Coefficient != 0.0)
// 		      {
// 			SumIndices = this->M1Value[m1] + this->M2Value[m1];
// 			Coefficient *= vSource[i];
// 			TmpNbrQ34Values = this->NbrQ34Values[m1];
// 			TmpM3Values = this->M3Values[m1];
// 			for (m3 = 0; m3 < TmpNbrQ34Values; ++m3)
// 			  {
// 			    Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
// 			    if (Index < Dim)			
// 			      vDestination.Re(Index) += Real( Coefficient * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
// 			    ++ReducedNbrInteractionFactors;
// 			  }
// 		      }
// 		    else
// 		      ReducedNbrInteractionFactors += this->NbrQ34Values[m1];
// 		  }
// 		vDestination[i] += this->HamiltonianShift * vSource[i];
// 	      }
	    
// 	  }
// 	if (this->OneBodyTermFlag == true)
// 	  {
// 	    for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
// 	      {
// 		m1 = this->OneBodyMValues[j];
// 		m2 = this->OneBodyNValues[j];
// 		TmpInteraction = this->OneBodyInteractionFactors[j];
// 		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
// 		  {
// 		    Index = this->Particles->AdA(i, m1, m2, Coefficient);
// 		    if (Index < Dim)
// 		      vDestination.Re(Index) += Real( Coefficient * TmpInteraction * vSource[i];		  
// 		  }
// 	      }
// 	  }
//       }    

//   delete TmpParticles;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using disk storage option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiplyDiskStorage(ComplexVector& vSource, ComplexVector& vDestination, 
									   int firstComponent, int nbrComponent)
{
  cout << "AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiplyDiskStorage must be defined" << endl;
  exit(1);
//   double Coefficient;
//   int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
//   double* BufferCoefficientArray  = new double [this->BufferSize * this->MaxNbrInteractionPerComponent];
//   int TmpNbrIteration = nbrComponent / this->BufferSize;
//   int* TmpIndexArray;
//   double* TmpCoefficientArray;
//   int TmpNbrInteraction;
//   int k = firstComponent;
//   int EffectiveHilbertSpaceDimension;
//   firstComponent -= this->PrecalculationShift;
  
//   ifstream File;
//   File.open(this->DiskStorageFileName, ios::binary | ios::in);
//   File.read ((char*) &EffectiveHilbertSpaceDimension, sizeof(int));
//   long FileJump = 0;
//   for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
//     FileJump += (long) this->NbrInteractionPerComponent[i];
//   FileJump *= sizeof(int);
//   long FileOffset = 0;
//   for (int i = this->DiskStorageStart; i < firstComponent; ++i)
//     FileOffset += this->NbrInteractionPerComponent[i];
//   File.seekg (((FileOffset + EffectiveHilbertSpaceDimension + 1) * sizeof(int)), ios::cur);
//   FileJump += (sizeof(double) - sizeof(int)) * FileOffset;
  
//   for (int i = 0; i < TmpNbrIteration; ++i)
//     {
//       int TmpPos = firstComponent;
//       long ReadBlockSize = 0;
//       for (int j = 0; j < this->BufferSize; ++j)
// 	{
// 	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
// 	  ++TmpPos;
// 	}		  
//       File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
//       FileJump -= sizeof(int) * ReadBlockSize;
//       File.seekg (FileJump, ios::cur);
//       File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
//       FileJump += sizeof(double) * ReadBlockSize;
//       File.seekg (-FileJump, ios::cur);
      
//       TmpIndexArray = BufferIndexArray;
//       TmpCoefficientArray = BufferCoefficientArray;
//       for (int l = 0; l < this->BufferSize; ++l)
// 	{
// 	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
// 	  Coefficient = vSource[k];
// 	  if (TmpNbrInteraction > 0)
// 	    {
// 	      for (int j = 0; j < TmpNbrInteraction; ++j)
// 		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
// 	      TmpIndexArray += TmpNbrInteraction;
// 	      TmpCoefficientArray += TmpNbrInteraction;
// 	    }
// 	  vDestination[k] += this->HamiltonianShift * Coefficient;
// 	  ++k;
// 	  ++firstComponent;
// 	}
//     }
  
//   if ((TmpNbrIteration * this->BufferSize) != nbrComponent)
//     {
//       int TmpPos = firstComponent;
//       int Lim =  nbrComponent % this->BufferSize;
//       long ReadBlockSize = 0;
//       for (int j = 0; j < Lim; ++j)
// 	{
// 	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
// 	  ++TmpPos;
// 	}		  
//       File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
//       FileJump -= sizeof(int) * ReadBlockSize;
//       File.seekg (FileJump, ios::cur);
//       File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
//       FileJump += sizeof(double) * ReadBlockSize;
//       File.seekg (-FileJump, ios::cur);
      
//       TmpIndexArray = BufferIndexArray;
//       TmpCoefficientArray = BufferCoefficientArray;
//       for (int i = 0; i < Lim; ++i)
// 	{
// 	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
// 	  Coefficient = vSource[k];
// 	  if (TmpNbrInteraction > 0)
// 	    {
// 	      for (int j = 0; j < TmpNbrInteraction; ++j)
// 		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
// 	      TmpIndexArray += TmpNbrInteraction;
// 	      TmpCoefficientArray += TmpNbrInteraction;
// 	    }
// 	  vDestination[k] += this->HamiltonianShift * Coefficient;
// 	  ++k;
// 	  ++firstComponent;
// 	}
//     }
  
//   File.close();
//   delete[] BufferIndexArray;
//   delete[] BufferCoefficientArray;
  return vDestination;
}


// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									int firstComponent, int nbrComponent)
{
  cout << "Using non-defined function LowLevelMultipleAddMultiply!"<<endl;
  int LastComponent = firstComponent + nbrComponent - 1; 
  if (true) // test for fast multiplication
    {
      Complex Coefficient;
      // copy single vector occurrence after testing and replace things like
      int Index=0, i=0;
      for (int l = 0; l < nbrVectors; ++l)
	vDestinations[l][Index] += Coefficient * vSources[l][i];  
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  unsigned short* TmpCoefficientIndexArray;
	  double TmpRe, TmpIm;
	  unsigned short TmpNbrRealInteraction;
	  unsigned short TmpNbrComplexInteraction;
	  Complex *TmpCPtr;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
	      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
	      for (int l = 0; l < nbrVectors; ++l)
		{		  
		  TmpRe = vSources[l][k].Re;
		  TmpIm = vSources[l][k].Im;
		  int Pos=0;
		  for (; Pos < TmpNbrRealInteraction; ++Pos)
		    {
		      vDestinations[l].Re(TmpIndexArray[Pos]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpRe;
		      vDestinations[l].Im(TmpIndexArray[Pos]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpIm;
		    }
		  for (int j=0; j < TmpNbrComplexInteraction; ++j)
		    {
		      TmpCPtr= &(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]]);
		      vDestinations[l].Re(TmpIndexArray[Pos]) +=  TmpCPtr->Re*TmpRe-TmpCPtr->Im*TmpIm;
		      vDestinations[l].Im(TmpIndexArray[Pos]) +=  TmpCPtr->Re*TmpIm+TmpCPtr->Im*TmpRe;
		      ++Pos;
		    }
		  vDestinations[l].Re(k) += this->HamiltonianShift * TmpRe;
		  vDestinations[l].Im(k) += this->HamiltonianShift * TmpIm;	      
		  ++k;
		  
		}
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	}
    }
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
											   int firstComponent, int nbrComponent)
{
//   int LastComponent = firstComponent + nbrComponent;
//   int Dim = this->Particles->GetHilbertSpaceDimension();
//   double Coefficient;
//   ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
//   int* TmpIndexArray;
//   Complex* TmpCoefficientArray; 
//   Complex* Coefficient2 = new Complex [nbrVectors];
//   int j;
//   int TmpNbrInteraction;
//   firstComponent -= this->PrecalculationShift;
//   LastComponent -= this->PrecalculationShift;
//   int Pos2;
//   int Pos = firstComponent / this->FastMultiplicationStep; 
//   int PosMod = firstComponent % this->FastMultiplicationStep;
//   if (PosMod != 0)
//     {
//       ++Pos;
//       PosMod = this->FastMultiplicationStep - PosMod;
//     }
//   int l =  PosMod + firstComponent + this->PrecalculationShift;
//   for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
//     {
//       TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
//       TmpIndexArray = this->InteractionPerComponentIndex[Pos];
//       TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
//       for (int k = 0; k < nbrVectors; ++k)
// 	{
// 	  Coefficient2[k] = vSources[k][l];
// 	  vDestinations[k][l] += this->HamiltonianShift * Coefficient2[k];
// 	}
//       for (j = 0; j < TmpNbrInteraction; ++j)
// 	{
// 	  Pos2 = TmpIndexArray[j];
// 	  Coefficient = TmpCoefficientArray[j];
// 	  for (int k = 0; k < nbrVectors; ++k)
// 	    {
// 	      vDestinations[k][Pos2] += Coefficient  * Coefficient2[k];
// 	    }
// 	}
//       l += this->FastMultiplicationStep;
//       ++Pos;
//     }
//   int Index;
//   int m1;
//   int m2;
//   int m3;
//   int m4;
//   double TmpInteraction;
//   firstComponent += this->PrecalculationShift;
//   LastComponent += this->PrecalculationShift;
//   for (int k = 0; k < this->FastMultiplicationStep; ++k)
//     if (PosMod != k)
//       {	
// 	if (this->NbrM12Indices == 0)
// 	  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
// 	    {
// 	      m1 = this->M1Value[j];
// 	      m2 = this->M2Value[j];
// 	      m3 = this->M3Value[j];
// 	      TmpInteraction = this->InteractionFactors[j];
// 	      m4 = m1 + m2 - m3;
// 	      for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
// 		{
// 		  Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
// 		  if (Index < Dim)
// 		    {
// 		      Coefficient *= TmpInteraction;
// 		      for (int l = 0; l < nbrVectors; ++l)
// 			vDestinations[l][Index] += Coefficient * vSources[l][i];
// 		    }
// 		}
// 	    }
// 	else
// 	  {
// 	    double Coefficient2;
// 	    int SumIndices;
// 	    int TmpNbrQ34Values;
// 	    int* TmpM3Values;
// 	    double* TmpCoefficients = new double[nbrVectors];
// 	    int ReducedNbrInteractionFactors;
// 	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
// 	      {
// 		ReducedNbrInteractionFactors = 0;
// 		for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
// 		  {
// 		    Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
// 		    if (Coefficient != 0.0)
// 		      {
// 			SumIndices = this->M1Value[m1] + this->M2Value[m1];
// 			TmpNbrQ34Values = this->NbrQ34Values[m1];
// 			TmpM3Values = this->M3Values[m1];
// 			for (int l = 0; l < nbrVectors; ++l)
// 			  TmpCoefficients[l] = Coefficient * vSources[l][i];
// 			for (m3 = 0; m3 < TmpNbrQ34Values; ++m3)
// 			  {
// 			    Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
// 			    if (Index < Dim)
// 			      for (int l = 0; l < nbrVectors; ++l)
// 				vDestinations[l][Index] += TmpCoefficients[l] * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
// 			    ++ReducedNbrInteractionFactors;
// 			  }
// 		      }
// 		    else
// 		      ReducedNbrInteractionFactors += this->NbrQ34Values[m1];
// 		  }
// 	      }
// 	    delete[] TmpCoefficients;
// 	  }
// 	for (int l = 0; l < nbrVectors; ++l)
// 	  {
// 	    ComplexVector& TmpSourceVector = vSources[l];
// 	    ComplexVector& TmpDestinationVector = vDestinations[l];
// 	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
// 	      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
// 	  }
	
//       }
//   delete[] Coefficient2;
//   delete TmpParticles;
  return vDestinations;
}

// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using disk storage option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyDiskStorage(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										   int firstComponent, int nbrComponent)
{
  cout << "Using non-defined function LowLevelMultipleAddMultiplyDiskStorage!"<<endl;
  // Complex Coefficient;
//   int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
//   Complex* BufferCoefficientArray  = new Complex[this->BufferSize * this->MaxNbrInteractionPerComponent];
//   Complex* Coefficient2 = new Complex [nbrVectors];
//   int TmpNbrIteration = nbrComponent / this->BufferSize;
//   int* TmpIndexArray;
//   Complex* TmpCoefficientArray;
//   int TmpNbrInteraction;
//   int k = firstComponent;
//   int EffectiveHilbertSpaceDimension;
//   int Pos;
//   firstComponent -= this->PrecalculationShift;
  
//   ifstream File;
//   File.open(this->DiskStorageFileName, ios::binary | ios::in);
//   File.read ((char*) &EffectiveHilbertSpaceDimension, sizeof(int));
//   long FileJump = 0;
//   for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
//     FileJump += (long) this->NbrInteractionPerComponent[i];
//   FileJump *= sizeof(int);
//   long FileOffset = 0;
//   for (int i = this->DiskStorageStart; i < firstComponent; ++i)
//     FileOffset += this->NbrInteractionPerComponent[i];
//   File.seekg (((FileOffset + EffectiveHilbertSpaceDimension + 1) * sizeof(int)), ios::cur);
//   FileJump += (sizeof(double) - sizeof(int)) * FileOffset;
  
//   for (int i = 0; i < TmpNbrIteration; ++i)
//     {
//       int TmpPos = firstComponent;
//       long ReadBlockSize = 0;
//       for (int j = 0; j < this->BufferSize; ++j)
// 	{
// 	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
// 	  ++TmpPos;
// 	}		  
//       File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
//       FileJump -= sizeof(int) * ReadBlockSize;
//       File.seekg (FileJump, ios::cur);
//       File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
//       FileJump += sizeof(double) * ReadBlockSize;
//       File.seekg (-FileJump, ios::cur);
      
//       TmpIndexArray = BufferIndexArray;
//       TmpCoefficientArray = BufferCoefficientArray;
//       for (int m = 0; m < this->BufferSize; ++m)
// 	{
// 	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
// 	  for (int l = 0; l < nbrVectors; ++l)
// 	    {
// 	      Coefficient2[l] = vSources[l][k];
// 	      vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
// 	    }
// 	  if (TmpNbrInteraction > 0)
// 	    {
// 	      for (int j = 0; j < TmpNbrInteraction; ++j)
// 		{
// 		  Pos = TmpIndexArray[j];
// 		  Coefficient = TmpCoefficientArray[j];
// 		  for (int l = 0; l < nbrVectors; ++l)
// 		    {
// 		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
// 		    }
// 		}
// 	      TmpIndexArray += TmpNbrInteraction;
// 	      TmpCoefficientArray += TmpNbrInteraction;
// 	    }
// 	  ++k;
// 	  ++firstComponent;
// 	}
//     }
  
//   if ((TmpNbrIteration * this->BufferSize) != nbrComponent)
//     {
//       int TmpPos = firstComponent;
//       int Lim =  nbrComponent % this->BufferSize;
//       long ReadBlockSize = 0;
//       for (int j = 0; j < Lim; ++j)
// 	{
// 	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
// 	  ++TmpPos;
// 	}		  
//       File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
//       FileJump -= sizeof(int) * ReadBlockSize;
//       File.seekg (FileJump, ios::cur);
//       File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
//       FileJump += sizeof(double) * ReadBlockSize;
//       File.seekg (-FileJump, ios::cur);
      
//       TmpIndexArray = BufferIndexArray;
//       TmpCoefficientArray = BufferCoefficientArray;
//       for (int m = 0; m < Lim; ++m)
// 	{
// 	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
// 	  for (int l = 0; l < nbrVectors; ++l)
// 	    {
// 	      Coefficient2[l] = vSources[l][k];
// 	      vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
// 	    }
// 	  if (TmpNbrInteraction > 0)
// 	    {
// 	      for (int j = 0; j < TmpNbrInteraction; ++j)
// 		{
// 		  Pos = TmpIndexArray[j];
// 		  Coefficient = TmpCoefficientArray[j];
// 		  for (int l = 0; l < nbrVectors; ++l)
// 		    {
// 		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
// 		    }
// 		}
// 	      TmpIndexArray += TmpNbrInteraction;
// 	      TmpCoefficientArray += TmpNbrInteraction;
// 	    }
// 	  ++k;
// 	  ++firstComponent;
// 	}
//     }
  
//   File.close();
//   delete[] BufferIndexArray;
//   delete[] BufferCoefficientArray;
//   delete[] Coefficient2;
  return vDestinations;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  vDestination.ClearVector();
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							 int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination.Re(i) = 0.0;
      vDestination.Im(i) = 0.0;
    }
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored
RealVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								      int firstComponent, int nbrComponent)
{
  return vDestination;
}
 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> AbstractQHEOnLatticeHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> AbstractQHEOnLatticeHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// if allowedMemory == 0, the value from preceding calls is used
// return value = amount of memory needed

long AbstractQHEOnLatticeHamiltonian::FastMultiplicationMemory(long allowedMemory)
{
  this->LoadedPrecalculation=false;
  if (allowedMemory>0)
    this->AllowedMemory = allowedMemory;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  this->NbrRealInteractionPerComponent = new unsigned short [EffectiveHilbertSpaceDimension];
  this->NbrComplexInteractionPerComponent = new unsigned short [EffectiveHilbertSpaceDimension];   
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    {
      this->NbrRealInteractionPerComponent[i] = 0x0;
      this->NbrComplexInteractionPerComponent[i] = 0x0;
    }
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start memory" << endl;

  QHEParticlePrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);
  long Memory = 0;
   
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    {
      Memory += this->NbrRealInteractionPerComponent[i];
      Memory += this->NbrComplexInteractionPerComponent[i];
    }
  
  cout << "nbr interaction = " << Memory << endl;
  
  // memory requirement, ignoring the actual storage size of the values of matrix
  // elements, which is assumed small (maybe need to add an estimate, at least)
  long TmpMemory = AllowedMemory - (2*sizeof (unsigned short) + sizeof (int*) + sizeof(unsigned short*)) * EffectiveHilbertSpaceDimension;
  cout << "of which can be stored: "<<(TmpMemory / ((int) (sizeof (int) + sizeof(unsigned short))))<<endl;
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(unsigned short)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(unsigned short)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if (this->Particles->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  // memory requirement, ignoring the actual storage size of the values of matrix
	  // elements, which is assumed small (maybe need to add an estimate, at least, again!)
	  TmpMemory = AllowedMemory - (2*sizeof (unsigned short) + sizeof (int*) + sizeof(unsigned short*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	    {
	      Memory += this->NbrRealInteractionPerComponent[i];
	      Memory += this->NbrComplexInteractionPerComponent[i];
	    }	  
	}
      
      Memory = ((2*sizeof (unsigned short) + sizeof (int*) + sizeof(unsigned short*)) * ReducedSpaceDimension) + (Memory * (sizeof (int) + sizeof(unsigned short)));
      
      if (this->DiskStorageFlag == false)
	{
	  int TotalReducedSpaceDimension = ReducedSpaceDimension;
	  unsigned short* TmpNbrRealInteractionPerComponent = new unsigned short [TotalReducedSpaceDimension];
	  unsigned short* TmpNbrComplexInteractionPerComponent = new unsigned short [TotalReducedSpaceDimension];	  
	  int Pos = 0;
	  for (int i = 0; i < ReducedSpaceDimension; ++i)
	    {
	      TmpNbrRealInteractionPerComponent[i] = this->NbrRealInteractionPerComponent[Pos];
	      TmpNbrComplexInteractionPerComponent[i] = this->NbrComplexInteractionPerComponent[Pos];
	      Pos += this->FastMultiplicationStep;
	    }
	  delete[] this->NbrRealInteractionPerComponent;
	  delete[] this->NbrComplexInteractionPerComponent;
	  this->NbrRealInteractionPerComponent = TmpNbrRealInteractionPerComponent;
	  this->NbrComplexInteractionPerComponent = TmpNbrComplexInteractionPerComponent;
	}      
    }
  else
    {
      Memory = ((2*sizeof (unsigned short) + sizeof (int*) + sizeof(unsigned short*)) * EffectiveHilbertSpaceDimension) + (Memory * (sizeof (int) + sizeof(unsigned short)));
      this->FastMultiplicationStep = 1;
    }

  cout << "reduction factor=" << this->FastMultiplicationStep << endl;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  cout << "final Memory in bytes = " <<Memory<<endl;
  return Memory;    
}



// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that have to be precalcualted
// return value = number of non-zero matrix elements
//
long AbstractQHEOnLatticeHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Index;
  double Coefficient;
  long Memory = 0;
  // deal with kinetic energy terms first!      
  int qi;
  int qf;
  Complex TmpInteraction;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();  
  for (int j = 0; j < NbrHoppingTerms; ++j)
    {
      qi = this->KineticQi[j];
      qf = this->KineticQf[j];
      TmpInteraction = this->HoppingTerms[j];
      if (fabs(TmpInteraction.Im)<1e-14)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpParticles->AdA(i, qf, qi, Coefficient);
	      if (Index < this->Particles->GetHilbertSpaceDimension())
		{
		  ++Memory;		
		  ++this->NbrRealInteractionPerComponent[i - this->PrecalculationShift];		  
		}
	    }
	}
      else
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpParticles->AdA(i, qf, qi, Coefficient);
	      if (Index < this->Particles->GetHilbertSpaceDimension())
		{
		  ++Memory;
		  ++this->NbrComplexInteractionPerComponent[i - this->PrecalculationShift];
		}
	    }
	}
    }
  
  // four-fermion interactions:
  if (this->NbrQ12Indices == 0) // full storage
    {
      for (int j = 0; j < NbrInteractionFactors; ++j) 
	{
	  int q1 = this->Q1Value[j];
	  int q2 = this->Q2Value[j];
	  int q3 = this->Q3Value[j];
	  int q4 = this->Q4Value[j];
	  TmpInteraction = this->InteractionFactors[j];
	  if (fabs(TmpInteraction.Im)<1e-14)
	    {
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		  if (Index < this->Particles->GetHilbertSpaceDimension())
		    {
		      ++Memory;
		      ++this->NbrRealInteractionPerComponent[i - this->PrecalculationShift];
		    }
		}
	    }
	  else
	    {
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		  if (Index < this->Particles->GetHilbertSpaceDimension())
		    {
		      ++Memory;
		      ++this->NbrComplexInteractionPerComponent[i - this->PrecalculationShift];
		    }
		}
	    }	 	  
	}
    }
  else // intelligent storage
    {
      double Coefficient2;
      int ProcessedNbrInteractionFactors;
      int TmpNbrQ34Values;
      int* TmpQ3Values;
      int* TmpQ4Values;
      for (int i = firstComponent; i < LastComponent; ++i)
	{	  
	  ProcessedNbrInteractionFactors = 0;
	  for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
	    {
	      Coefficient = TmpParticles->AA(i, this->Q1Value[i12], this->Q2Value[i12]);
	      if (Coefficient != 0.0)
		{
		  TmpNbrQ34Values = this->NbrQ34Values[i12];
		  TmpQ3Values = this->Q3PerQ12[i12];
		  TmpQ4Values = this->Q4PerQ12[i12];
		  for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
		    {
		      Index = TmpParticles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
		      if (Index < this->Particles->GetHilbertSpaceDimension())
			{
			  ++Memory;
			  if (fabs(this->InteractionFactors[ProcessedNbrInteractionFactors].Im)<1e-14)
			    ++this->NbrRealInteractionPerComponent[i - this->PrecalculationShift];
			  else
			    ++this->NbrComplexInteractionPerComponent[i - this->PrecalculationShift];
			}
		      ++ProcessedNbrInteractionFactors;
		    }
		}
	      else
		ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
	    }
	}
    }	  
  
  // separated diagonal terms as these will be the general rule for contact interactions
  if (NbrDiagonalInteractionFactors>0)
    {	  
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  Coefficient = TmpParticles->AdAdAADiagonal(i, NbrDiagonalInteractionFactors,
						     DiagonalInteractionFactors, DiagonalQValues);
	  if (fabs(Coefficient)>1e-14)
	    {
	      ++Memory;
	      ++this->NbrRealInteractionPerComponent[i - this->PrecalculationShift];
	    }
	}
    }
  
  delete TmpParticles;

  return Memory;
}

// enable fast multiplication algorithm
//

void AbstractQHEOnLatticeHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  int Index;
  int qi;
  int qf;
  int tmpElementPos;
  double Coefficient;
  int* TmpIndexArray;
  unsigned short* TmpCoefficientIndexArray;
  int PosR, PosC;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficientIndex = new unsigned short* [ReducedSpaceDimension];

  int TotalPos = 0;
  int TmpInteractionPerComponent;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
    {
      TmpInteractionPerComponent = this->NbrRealInteractionPerComponent[TotalPos]+this->NbrComplexInteractionPerComponent[TotalPos];
      this->InteractionPerComponentIndex[TotalPos] = new int [TmpInteractionPerComponent];
      this->InteractionPerComponentCoefficientIndex[TotalPos] = new unsigned short [TmpInteractionPerComponent];      
      TmpIndexArray = this->InteractionPerComponentIndex[TotalPos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[TotalPos];
      PosR = 0;  // counter for position of real matrix elements
      PosC = this->NbrRealInteractionPerComponent[TotalPos];  // counter for position of complex matrix elements

      // deal with kinetic energy terms first!             
      for (int j = 0; j < NbrHoppingTerms; ++j)
	{
	  qi = this->KineticQi[j];
	  qf = this->KineticQf[j];
	  // considering: this->HoppingTerms[j]
	  Index = TmpParticles->AdA(i, qf, qi, Coefficient);	  
	  if (Index < Dim)
	    {
	      //cout << "Element ("<<qi<<"->"<<qf<<"): "<<Coefficient<<endl;
	      if (fabs(this->HoppingTerms[j].Im)<1e-14) // real element
		{
		  TmpIndexArray[PosR] = Index;
		  tmpElementPos = RealInteractionCoefficients.InsertElement
		    (Coefficient*this->HoppingTerms[j].Re);
		  if (tmpElementPos > USHRT_MAX )
		    {
		      cout << "Error: too many different real matrix elements for fast storage"<<endl;
		      exit(1);
		    }
		  TmpCoefficientIndexArray[PosR] = (unsigned short) tmpElementPos;
 		  ++PosR;
		}
	      else
		{
		  TmpIndexArray[PosC] = Index;
		  tmpElementPos = ComplexInteractionCoefficients.InsertElement
		    (Coefficient*this->HoppingTerms[j]);
		  if (tmpElementPos > USHRT_MAX )
		    {
		      cout << "Error: too many different complex matrix elements for fast storage"<<endl;
		      exit(1);
		    }
		  TmpCoefficientIndexArray[PosC] = (unsigned short) tmpElementPos;
		  ++PosC;
		}
	      // cout << "connecting :"<<Index<<", "<<i<<": "<<Coefficient*this->HoppingTerms[j]<<endl;
	    }
	}

      // four-fermion interactions:
      if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];	       
	      Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
	      if (Index < Dim)
		{
		  if (fabs(this->InteractionFactors[j].Im)<1e-14) // real element
		    {
		      TmpIndexArray[PosR] = Index;
		      tmpElementPos = RealInteractionCoefficients.InsertElement
			(Coefficient*this->InteractionFactors[j].Re);
		      if (tmpElementPos > USHRT_MAX )
			{
			  cout << "Error: too many different real matrix elements for fast storage"<<endl;
			  exit(1);
			}
		      TmpCoefficientIndexArray[PosR] = (unsigned short) tmpElementPos;
		      ++PosR;
		    }
		  else
		    {
		      TmpIndexArray[PosC] = Index;
		      tmpElementPos = ComplexInteractionCoefficients.InsertElement
			(Coefficient*this->InteractionFactors[j]);
		      if (tmpElementPos > USHRT_MAX )
			{
			  cout << "Error: too many different complex matrix elements for fast storage"<<endl;
			  exit(1);
			}
		      TmpCoefficientIndexArray[PosC] = (unsigned short) tmpElementPos;
		      ++PosC;
		    }
		  //cout << "4b - connecting :"<<Index<<", "<<i<<": "<<Coefficient*this->InteractionFactors[j]<< " (q's=["<<q1<<","<<q2<<","<<q3<<","<<q4<<"])"<<endl;
		}
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2;
	  int ProcessedNbrInteractionFactors = 0;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	  for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
	    {
	      Coefficient = TmpParticles->AA(i, this->Q1Value[i12], this->Q2Value[i12]);
	      if (Coefficient != 0.0)
		{
		  TmpNbrQ34Values = this->NbrQ34Values[i12];
		  TmpQ3Values = this->Q3PerQ12[i12];
		  TmpQ4Values = this->Q4PerQ12[i12];
		  for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
		    {
		      Index = TmpParticles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
		      if (Index < Dim)
			{
			  if (fabs(this->InteractionFactors[ProcessedNbrInteractionFactors].Im)<1e-14) 
			    {
			      TmpIndexArray[PosR] = Index;
			      tmpElementPos = RealInteractionCoefficients.InsertElement
				(Coefficient*this->InteractionFactors[ProcessedNbrInteractionFactors].Re);
			      if (tmpElementPos > USHRT_MAX )
				{
				  cout << "Error: too many different real matrix elements for fast storage"<<endl;
				  exit(1);
				}
			      TmpCoefficientIndexArray[PosR] = (unsigned short) tmpElementPos;
			      ++PosR;
			    }
			  else
			    {
			      TmpIndexArray[PosC] = Index;
			      tmpElementPos = ComplexInteractionCoefficients.InsertElement
				(Coefficient*this->InteractionFactors[ProcessedNbrInteractionFactors]);
			      if (tmpElementPos > USHRT_MAX )
				{
				  cout << "Error: too many different complex matrix elements for fast storage"<<endl;
				  exit(1);
				}
			      TmpCoefficientIndexArray[PosC] = (unsigned short) tmpElementPos;
			      ++PosC;
			    }
			  ++ProcessedNbrInteractionFactors;
			}		     
		      else
			ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		    }
		}
	    }
	}

      // separated diagonal terms as these will be the general rule for contact interactions
      if (NbrDiagonalInteractionFactors>0)
	{
	  Coefficient = TmpParticles->AdAdAADiagonal(i, NbrDiagonalInteractionFactors,
						     DiagonalInteractionFactors, DiagonalQValues);
	  if (fabs(Coefficient)>1e-14)
	    {
	      TmpIndexArray[PosR] = i;
	      tmpElementPos = RealInteractionCoefficients.InsertElement(Coefficient);
	      if (tmpElementPos > USHRT_MAX )
		{
		  cout << "Error: too many different real matrix elements for fast storage"<<endl;
		  exit(1);
		}
	      TmpCoefficientIndexArray[PosR] = (unsigned short) tmpElementPos;
	      ++PosR;
	      // cout << "diag - connecting :"<<i<<", "<<i<<": "<<Coefficient<<endl;
	    }	   
	}
      ++TotalPos;
    }
      
  delete TmpParticles;

  cout << "Nbr distinct matrix elements: "<<RealInteractionCoefficients.GetNbrElements()<<" real, "
       << ComplexInteractionCoefficients.GetNbrElements()<<" complex"<<endl;
   
  this->FastMultiplicationFlag = true;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted

void AbstractQHEOnLatticeHamiltonian::PartialEnableFastMultiplication(int firstComponent, int lastComponent)
{
  cout << "Using non-defined function PartialEnableFastMultiplication!"<<endl;  
//   int Index;
//   double Coefficient;
//   int m1;
//   int m2;
//   int m3;
//   int m4;
//   int* TmpIndexArray;
//   double* TmpCoefficientArray;
//   int Pos;
//   int Min = firstComponent / this->FastMultiplicationStep;
//   int Max = lastComponent / this->FastMultiplicationStep;
//   ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
 
//   for (int i = Min; i < Max; ++i)
//     {
//       this->InteractionPerComponentIndex[i] = new int [this->NbrInteractionPerComponent[i]];
//       this->InteractionPerComponentCoefficient[i] = new double [this->NbrInteractionPerComponent[i]];      
//       TmpIndexArray = this->InteractionPerComponentIndex[i];
//       TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
//       Pos = 0;
//       for (int j = 0; j < this->NbrInteractionFactors; ++j) 
// 	{
// 	  m1 = this->M1Value[j];
// 	  m2 = this->M2Value[j];
// 	  m3 = this->M3Value[j];
// 	  m4 = m1 + m2 - m3;
// 	  Index = TmpParticles->AdAdAA(i * this->FastMultiplicationStep, m1, m2, m3, m4, Coefficient);
// 	  if (Index < this->Particles->GetHilbertSpaceDimension())
// 	    {
// 	      TmpIndexArray[Pos] = Index;
// 	      TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[j];
// 	      ++Pos;
// 	    }
// 	}
//       if (this->OneBodyTermFlag == true)
// 	{
// 	  for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
// 	    {
// 	      m1 = this->OneBodyMValues[j];
// 	      m2 = this->OneBodyNValues[j];
// 	      Index = this->Particles->AdA(i + this->PrecalculationShift, m1, m2, Coefficient);
// 	      if (Index < this->Particles->GetHilbertSpaceDimension())
// 		{
// 		  TmpIndexArray[Pos] = Index;
// 		  TmpCoefficientArray[Pos] = Coefficient * this->OneBodyInteractionFactors[j];
// 		  ++Pos;
// 		}
// 	    }
// 	}
//    }
//   delete TmpParticles;
}

// enable fast multiplication algorithm using on disk cache 
//
// fileName = prefix of the name of the file where temporary matrix elements will be stored

void AbstractQHEOnLatticeHamiltonian::EnableFastMultiplicationWithDiskStorage(char* fileName)
{
  cout << "Using non-defined function EnableFastMultiplicationWithDiskStorage!"<<endl;
//   if (this->FastMultiplicationStep == 1)
//     {
//       this->DiskStorageFlag = false;
//       this->DiskStorageFileName = 0;
//       this->EnableFastMultiplication();
//       return;
//     }
//   this->DiskStorageFlag = true;
//   this->DiskStorageFileName = new char [strlen(fileName) + 8];
//   sprintf (this->DiskStorageFileName, "%s.ham", fileName);

//   long MinIndex;
//   long MaxIndex;
//   this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
//   int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
//   this->DiskStorageStart = (int) MinIndex;
//   int DiskStorageEnd = 1 + (int) MaxIndex;

//   int Index;
//   int m1;
//   int m2;
//   int m3;
//   int m4;
//   double Coefficient;
//   int* TmpIndexArray;
//   double* TmpCoefficientArray;
//   int Pos;
//   timeval TotalStartingTime2;
//   timeval TotalEndingTime2;
//   double Dt2;
//   gettimeofday (&(TotalStartingTime2), 0);
//   cout << "start" << endl;
//   this->InteractionPerComponentIndex = 0;
//   this->InteractionPerComponentCoefficient = 0;
//   this->MaxNbrInteractionPerComponent = 0;

//   int TotalPos = 0;
//   ofstream File;
//   File.open(this->DiskStorageFileName, ios::binary | ios::out);
 
//   File.write((char*) &(EffectiveHilbertSpaceDimension), sizeof(int));
//   File.write((char*) &(this->FastMultiplicationStep), sizeof(int));
//   File.write((char*) this->NbrInteractionPerComponent, sizeof(int) * EffectiveHilbertSpaceDimension);

//   long FileJump = 0;
//   for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
//     {
//       FileJump += (long) this->NbrInteractionPerComponent[i];
//       if (this->MaxNbrInteractionPerComponent < this->NbrInteractionPerComponent[i])
// 	this->MaxNbrInteractionPerComponent = this->NbrInteractionPerComponent[i];
//     }
//   FileJump *= sizeof(int);

//   TmpIndexArray = new int [this->MaxNbrInteractionPerComponent];
//   TmpCoefficientArray = new double [this->MaxNbrInteractionPerComponent];      
//   double Coefficient2;
//   int SumIndices;
//   int TmpNbrQ34Values;
//   int* TmpM3Values;
//   int ReducedNbrInteractionFactors;

//   for (int i = this->DiskStorageStart; i < DiskStorageEnd; ++i)
//     {
//       if (this->NbrInteractionPerComponent[TotalPos] > 0)
// 	{
// 	  Pos = 0;
// 	  if (this->NbrM12Indices == 0)
// 	    {
// 	      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
// 		{
// 		  m1 = this->M1Value[j];
// 		  m2 = this->M2Value[j];
// 		  m3 = this->M3Value[j];
// 		  m4 = m1 + m2 - m3;
// 		  Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
// 		  if (Index < this->Particles->GetHilbertSpaceDimension())
// 		    {
// 		      TmpIndexArray[Pos] = Index;
// 		      TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[j];
// 		      ++Pos;
// 		    }
// 		}
// 	    }
// 	  else
// 	    {
// 	      ReducedNbrInteractionFactors = 0;
// 	      for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
// 		{
// 		  Coefficient = this->Particles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
// 		  if (Coefficient != 0.0)
// 		    {
// 		      SumIndices = this->M1Value[m1] + this->M2Value[m1];
// 		      TmpM3Values = this->M3Values[m1];
// 		      TmpNbrQ34Values = this->NbrQ34Values[m1];
// 		      for (m3 = 0; m3 < TmpNbrQ34Values; ++m3)
// 			{
// 			  Index = this->Particles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
// 			  if (Index < this->Particles->GetHilbertSpaceDimension())
// 			    {
// 			      TmpIndexArray[Pos] = Index;
// 			      TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
// 			      ++Pos;
// 			    }
// 			  ++ReducedNbrInteractionFactors;
// 			}    
// 		    }
// 		  else
// 		    ReducedNbrInteractionFactors += this->NbrQ34Values[m1];
// 		}	      
// 	    }
// 	  if (this->OneBodyTermFlag == true)
// 	    {
// 	      for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
// 		{
// 		  m1 = this->OneBodyMValues[j];
// 		  m2 = this->OneBodyNValues[j];
// 		  Index = this->Particles->AdA(i, m1, m2, Coefficient);
// 		  if (Index < this->Particles->GetHilbertSpaceDimension())
// 		    {
// 		      TmpIndexArray[Pos] = Index;
// 		      TmpCoefficientArray[Pos] = Coefficient * this->OneBodyInteractionFactors[j];
// 		      ++Pos;
// 		    }
// 		}
// 	    }
// 	  File.write((char*) TmpIndexArray, sizeof(int) * this->NbrInteractionPerComponent[TotalPos]);
// 	  FileJump -= sizeof(int) * this->NbrInteractionPerComponent[TotalPos];
// 	  File.seekp(FileJump, ios::cur);
// 	  File.write((char*) TmpCoefficientArray, sizeof(double) * this->NbrInteractionPerComponent[TotalPos]);
// 	  FileJump += sizeof(double) * this->NbrInteractionPerComponent[TotalPos];
// 	  File.seekp(-FileJump, ios::cur);	  
// 	}
//       ++TotalPos;
//     }
//   delete[] TmpIndexArray;
//   delete[] TmpCoefficientArray;
//   File.close();

//   this->FastMultiplicationFlag = true;
//   this->BufferSize = this->Memory / ((this->MaxNbrInteractionPerComponent * (sizeof(int) + sizeof(double))) + sizeof(int*) + sizeof(double*));

//   gettimeofday (&(TotalEndingTime2), 0);
//   cout << "------------------------------------------------------------------" << endl << endl;;
//   Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
//     ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
//   cout << "time = " << Dt2 << endl;
}

// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs

bool AbstractQHEOnLatticeHamiltonian::SavePrecalculation (char* fileName)
{
  if (this->FastMultiplicationFlag)
    {
      ofstream File;
      File.open(fileName, ios::binary | ios::out);
      int Tmp = this->Particles->GetHilbertSpaceDimension();
      File.write((char*) &(Tmp), sizeof(int));
      File.write((char*) &(this->FastMultiplicationStep), sizeof(int));
      Tmp /= this->FastMultiplicationStep;
      if ((Tmp * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
	++Tmp;
      File.write((char*) this->NbrRealInteractionPerComponent, sizeof(unsigned short) * Tmp);
      File.write((char*) this->NbrComplexInteractionPerComponent, sizeof(unsigned short) * Tmp);
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * (this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i]));
	}
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentCoefficientIndex[i]), sizeof(unsigned short) * (this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i]));	  
	}
      RealInteractionCoefficients.WriteArray(File);
      ComplexInteractionCoefficients.WriteArray(File);
      File.close();
      return true;
    }
  else
    {
      return false;
    }
}

// load precalculations from a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be read
// return value = true if no error occurs

bool AbstractQHEOnLatticeHamiltonian::LoadPrecalculation (char* fileName)
{
  this->LoadedPrecalculation=true;
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  int Tmp;
  File.read((char*) &(Tmp), sizeof(int));
  if (Tmp != this->Particles->GetHilbertSpaceDimension())
    {
      File.close();
      return false;
    }
  File.read((char*) &(this->FastMultiplicationStep), sizeof(int));
  Tmp /= this->FastMultiplicationStep;
  if ((Tmp * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
    ++Tmp;
  this->NbrRealInteractionPerComponent = new unsigned short [Tmp];
  this->NbrComplexInteractionPerComponent = new unsigned short [Tmp];
  File.read((char*) this->NbrRealInteractionPerComponent, sizeof(unsigned short) * Tmp);
  File.read((char*) this->NbrComplexInteractionPerComponent, sizeof(unsigned short) * Tmp);

  this->InteractionPerComponentIndex = new int* [Tmp];
  this->InteractionPerComponentCoefficientIndex = new unsigned short* [Tmp];
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [(this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i])];
      File.read((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * (this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i]));	  
    }
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentCoefficientIndex[i]=new unsigned short[(this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i])];
      File.read((char*) (this->InteractionPerComponentCoefficientIndex[i]), sizeof(unsigned short) * (this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i]));	  
    }
  RealInteractionCoefficients.ReadArray(File);
  ComplexInteractionCoefficients.ReadArray(File);
   
  File.close();
  this->FastMultiplicationFlag = true;
  return true;
}

