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
#include "GeneralTools/StringTools.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperationWithMatrixElements.h"

#include <iostream>
#include <sys/time.h>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <cstring>

using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;
using std::ostream;



// default constructor
//
AbstractQHEOnLatticeHamiltonian::AbstractQHEOnLatticeHamiltonian() :
  MaxElementIndex (std::numeric_limits<ElementIndexType>::max())
{
  this->NbrQ12Indices=0;
  this->NbrRealInteractionPerComponent=0;
  this->NbrComplexInteractionPerComponent=0;
  this->NbrDiagonalInteractionFactors=0;
  this->NbrRhoRhoInteractionFactors=0;
  this->LoadBalancingArray=0;
  this->NbrBalancedTasks=0;
  this->FastMultiplicationStep=0;
  this->HermitianSymmetryFlag=false;
  this->HaveComplexMatrixElements=true;
  this->HaveTestedForComplexMatrixElement=false;
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

      int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
	++ReducedSpaceDimension;
      for (int i=0; i<ReducedSpaceDimension; ++i)
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
  if (this->NbrRhoRhoInteractionFactors>0)  
    {
      delete [] RhoRhoInteractionFactors;
      delete [] RhoRhoQ12Values;
    }
  this->Particles = (ParticleOnLattice*) hilbertSpace;
  this->HaveTestedForComplexMatrixElement=false;
  this->EvaluateInteractionFactors();
}

// set flux density in units of flux quanta through the lattice
//
// nbrFluxQuanta = flux quantua piercing the lattice
void AbstractQHEOnLatticeHamiltonian::SetNbrFluxQuanta(int nbrFluxQuanta)
{
  this->NbrFluxQuanta=nbrFluxQuanta;
  this->FluxDensity=((double)nbrFluxQuanta)/NbrSites;
  this->Particles->SetNbrFluxQuanta(nbrFluxQuanta);
  if (NbrHoppingTerms>0)
    {
      delete [] this->HoppingTerms;
      delete [] this->KineticQi;
      delete [] this->KineticQf;
    }
  if (NbrInteractionFactors>0 && NbrQ12Indices==0)
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
  HaveTestedForComplexMatrixElement=false;
  this->EvaluateInteractionFactors();
  if (this->LoadedPrecalculation)
    {
      cout << "Cannot re-enable fast calculation when precalculation data was loaded from file!"<<endl;
      cout << "Reverting to slow calculation"<<endl;
    }
  else if (EnableFastCalculation)
    {
      int TmpMemory = this->FastMultiplicationMemory(0);
      cout  << "fast = ";
      PrintMemorySize(cout, TmpMemory)<<endl;
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

// ask if Hamiltonian implements hermitian symmetry operations
//
bool AbstractQHEOnLatticeHamiltonian::IsHermitian()
{
  return HermitianSymmetryFlag;
}

// ask if Hamiltonian implements conjugate methods
//
bool AbstractQHEOnLatticeHamiltonian::IsConjugate()
{
  return true;
}

// ask if Hamiltonian has complex matrix elements
//
bool AbstractQHEOnLatticeHamiltonian::IsComplex()
{
  if (this->HaveTestedForComplexMatrixElement==true)
    {
      return this->HaveComplexMatrixElements;
    }
  else
    {
      this->HaveComplexMatrixElements=false;
      // deal with kinetic energy terms first!             
      for (int j = 0; j < NbrHoppingTerms; ++j)
	{
	  if (fabs(this->HoppingTerms[j].Im)>LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD) // real element
	    {
	      this->HaveComplexMatrixElements=true;
	      return true;
	    }
	}
      // four-fermion interactions:
      if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      if (fabs(this->InteractionFactors[j].Im)>LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD) // real element
		{
		  this->HaveComplexMatrixElements=true;
		  return true;
		}
	    }
	}
      else // intelligent storage
	{
	  int ProcessedNbrInteractionFactors=0;
	  for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
	    {
	      int TmpNbrQ34Values=NbrQ34Values[i12];
	      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
		{
		  if (fabs(this->InteractionFactors[ProcessedNbrInteractionFactors].Im)>LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD)
		    {
		      this->HaveComplexMatrixElements=true;
		      return true;
		    }
		  ++ProcessedNbrInteractionFactors;
		}
	    }
	}
    }
  return this->HaveComplexMatrixElements;
}


// count interaction terms
// return = number of interaction terms
long AbstractQHEOnLatticeHamiltonian::CountTwoBodyInteractionTerms()
{

  if (this->NbrQ12Indices == 0)
    {
      return this->NbrInteractionFactors; // full storage, so simple to answer this.
    }
  else
    {
      // count of interaction factors
      int CurrentNbrInteractionFactors=0;
      for (int q12 = 0; q12 < this->NbrQ12Indices; ++q12)
	CurrentNbrInteractionFactors+=this->NbrQ34Values[q12];
      return CurrentNbrInteractionFactors;
    }
}

// symmetrize interaction factors to enable hermitian matrix multiplication
// return = true upon success
bool AbstractQHEOnLatticeHamiltonian::HermitianSymmetrizeInteractionFactors()
{
  if (HermitianSymmetryFlag)
    return true;

  if (this->Particles->HaveOrder()==false)
    {
      cout << "Hamiltonian tried to use hermitian symmetry, but this is not implemented in HilbertSpace!"<<endl;
      HermitianSymmetryFlag=false;
      return false;
    }

  cout << "Using hermitian symmetry"<<endl;

  int *M = new int[2];
  int *N = new int[2];

  // single particle terms
  if (NbrHoppingTerms>0)
    {
      cout << "One-body hopping terms before hermitian symmetry: "<<this->NbrHoppingTerms<<endl;
      int TmpNbrHoppingTerms = 0;
      int *Flags = new int[this->NbrHoppingTerms];
      for (int j = 0; j < NbrHoppingTerms; ++j) 
	{
	  M[0] = this->KineticQi[j];
	  N[0] = this->KineticQf[j];
	  Flags[j] = this->Particles->CheckOrder(M, N, 1);
	  // cout << "M="<<M[0]<<", N="<<N[0]<<", order: "<<Flags[j]<<" element: "<<HoppingTerms[j]<<endl;
	  if (Flags[j]>0)
	    ++TmpNbrHoppingTerms;
	  else if (Flags[j]==0)
	    {
	      ++TmpNbrHoppingTerms;
	      HoppingTerms[j]*=0.5;
	    }
	}
      Complex *TmpHoppingTerms = new Complex[TmpNbrHoppingTerms];
      int *TmpQi = new int[TmpNbrHoppingTerms];
      int *TmpQf = new int[TmpNbrHoppingTerms];
      int Pos=0;
      for (int j = 0; j < this->NbrHoppingTerms; ++j)
	if (Flags[j]>=0)
	  {
	    TmpHoppingTerms[Pos]=this->HoppingTerms[j];
	    TmpQi[Pos]=this->KineticQi[j];
	    TmpQf[Pos]=this->KineticQf[j];
	    ++Pos;
	  }
      delete [] this->HoppingTerms;
      delete [] this->KineticQi;
      delete [] this->KineticQf;
      delete [] Flags;
      this->HoppingTerms = TmpHoppingTerms;
      this->KineticQi = TmpQi;
      this->KineticQf = TmpQf; 
      this->NbrHoppingTerms = TmpNbrHoppingTerms;
      cout << "One-body hopping terms after hermitian symmetry: "<<this->NbrHoppingTerms<<endl;
    }

  cout << "Two-body hopping terms before hermitian symmetry: "<<this->CountTwoBodyInteractionTerms()<<endl;

  if (this->NbrQ12Indices == 0)
    {
      if (NbrInteractionFactors>0)
	{
	  int TmpNbrInteractionFactors = 0;
	  int *Flags = new int[NbrInteractionFactors];
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      M[0] = this->Q1Value[j];
	      M[1] = this->Q2Value[j];
	      N[0] = this->Q3Value[j];
	      N[1] = this->Q4Value[j];
	      Flags[j] = this->Particles->CheckOrder (M, N, 2);
	      cout << "Flag("<<this->Q1Value[j]<<", "<<this->Q2Value[j]<<", "<<this->Q3Value[j]<<", "<<this->Q4Value[j]<<")="<<Flags[j]<<endl;
	      if (Flags[j]>0)
		++TmpNbrInteractionFactors;
	      else
		{
		  if (Flags[j]==0)
		    {
		      ++TmpNbrInteractionFactors;
		      this->InteractionFactors[j]*=0.5; // diagonal term: make up for double counting
		    }
		  else
		    {
		      cout << "Discarding element "<<this->Q1Value[j]<<", "<<this->Q2Value[j]<<", "<<this->Q3Value[j]<<", "<<this->Q4Value[j]<<endl;
		    }
		}
	    }
	  Complex* TmpInteractionFactors = new Complex[TmpNbrInteractionFactors];
	  int* TmpQ1Value = new int[TmpNbrInteractionFactors];
	  int* TmpQ2Value = new int[TmpNbrInteractionFactors];
	  int* TmpQ3Value = new int[TmpNbrInteractionFactors];
	  int* TmpQ4Value = new int[TmpNbrInteractionFactors];
	  int Pos=0;
	  for (int j = 0; j < NbrInteractionFactors; ++j)
	    {
	      if (Flags[j]>=0)
		{
		  TmpInteractionFactors[Pos]=InteractionFactors[j];
		  TmpQ1Value[Pos]=Q1Value[j];
		  TmpQ2Value[Pos]=Q2Value[j];
		  TmpQ3Value[Pos]=Q3Value[j];
		  TmpQ4Value[Pos]=Q4Value[j];
		  ++Pos;
		}
	    }
	  delete [] InteractionFactors;
	  delete [] Q1Value;
	  delete [] Q2Value;
	  delete [] Q3Value;
	  delete [] Q4Value;
	  this->InteractionFactors = TmpInteractionFactors;
	  this->NbrInteractionFactors = TmpNbrInteractionFactors;
	  this->Q1Value = TmpQ1Value;
	  this->Q2Value = TmpQ2Value;
	  this->Q3Value = TmpQ3Value;
	  this->Q4Value = TmpQ4Value;
	  delete [] Flags;
	}
    }
  else
    {
      int OldNbrQ34Values;
      int* OldQ3PerQ12;
      int* OldQ4PerQ12;
      int TmpNbrQ12Values = 0;
      int* Q12Flags = new int[this->NbrQ12Indices];
      int TmpNbrQ34Values;
      int* TmpQ3PerQ12;
      int* TmpQ4PerQ12;
      int* Q34Flags;
      // quick 5count of interaction factors
      int OldNbrInteractionFactors=0;
      for (int q12 = 0; q12 < this->NbrQ12Indices; ++q12)
	OldNbrInteractionFactors+=this->NbrQ34Values[q12];      
      int TmpNbrInteractionFactors=0;
      Complex *TmpInteractionFactors=new Complex[OldNbrInteractionFactors];
      int Pos=0;
      for (int q12 = 0; q12 < this->NbrQ12Indices; ++q12)
	{
	  M[0]=this->Q1Value[q12];
	  M[1]=this->Q2Value[q12];
	  OldNbrQ34Values = this->NbrQ34Values[q12];
	  OldQ3PerQ12 = this->Q3PerQ12[q12];
	  OldQ4PerQ12 = this->Q4PerQ12[q12];
	  TmpNbrQ34Values = 0;
	  Q34Flags = new int[OldNbrQ34Values];
	  for (int q34 = 0; q34 < OldNbrQ34Values; ++q34)
	    {
	      N[0]=OldQ3PerQ12[q34];
	      N[1]=OldQ4PerQ12[q34];
	      Q34Flags[q34] = this->Particles->CheckOrder(M, N, 2);
	      //cout << "Flag("<<M[0]<<", "<<M[1]<<", "<<N[0]<<", "<<N[1]<<")="<<Q34Flags[q34]<<endl;
	      if (Q34Flags[q34]>0)
		{
		  ++TmpNbrQ34Values;
		  TmpInteractionFactors[TmpNbrInteractionFactors++]=this->InteractionFactors[Pos];
		}
	      else if (Q34Flags[q34]==0)
		{
		  ++TmpNbrQ34Values;
		  TmpInteractionFactors[TmpNbrInteractionFactors++]=0.5*this->InteractionFactors[Pos];
		}
// 	      else
// 		{
// 		  cout << "Discarding element "<<M[0]<<", "<<M[1]<<", "<<N[0]<<", "<<N[1]<<endl;
// 		}
	      
	      ++Pos;
	    }
	  if (TmpNbrQ34Values>0)
	    {
	      //cout << "Q1="<<M[0]<<", Q2="<<M[1]<<": ";
	      ++TmpNbrQ12Values;
	      Q12Flags[q12]=1;
	      TmpQ3PerQ12 = new int[TmpNbrQ34Values];
	      TmpQ4PerQ12 = new int[TmpNbrQ34Values];
	      int Pos2=0;
	      for (int q34 = 0; q34 < OldNbrQ34Values; ++q34)
		if (Q34Flags[q34]>=0)
		  {
		    TmpQ3PerQ12[Pos2]=OldQ3PerQ12[q34];
		    TmpQ4PerQ12[Pos2]=OldQ4PerQ12[q34];
		    //cout << " Q3: " << TmpQ3PerQ12[Pos2];
		    //cout << " Q4: " << TmpQ4PerQ12[Pos2];
		    Pos2++;
		  }
	      //cout << endl;
	      delete [] OldQ3PerQ12;
	      delete [] OldQ4PerQ12;
	      this->Q3PerQ12[q12] = TmpQ3PerQ12;
	      this->Q4PerQ12[q12] = TmpQ4PerQ12;
	      this->NbrQ34Values[q12] = TmpNbrQ34Values;
 	    }
	  else
	    {
	      Q12Flags[q12]=-1;
	      delete [] OldQ3PerQ12;
	      delete [] OldQ4PerQ12;
	    }
	  delete [] Q34Flags;
	}
      if (this->NbrQ12Indices!=TmpNbrQ12Values)
	{
	  int *NewQ1Value=new int[TmpNbrQ12Values];
	  int *NewQ2Value=new int[TmpNbrQ12Values];
	  int **NewQ3PerQ12=new int*[TmpNbrQ12Values];
	  int **NewQ4PerQ12=new int*[TmpNbrQ12Values];
	  int *NewNbrQ34Values=new int[TmpNbrQ12Values];
	  Pos = 0;
	  for (int q12 = 0; q12 < this->NbrQ12Indices; ++q12)
	    if (Q12Flags[q12]>0)
	      {
		NewQ1Value[Pos]=this->Q1Value[q12];
		NewQ2Value[Pos]=this->Q2Value[q12];
		NewQ3PerQ12[Pos]=this->Q3PerQ12[q12];
		NewQ4PerQ12[Pos]=this->Q4PerQ12[q12];
		NewNbrQ34Values[Pos]=this->NbrQ34Values[q12];
		++Pos;
	      }
	  delete [] this->Q1Value;
	  delete [] this->Q2Value;
	  delete [] this->Q3PerQ12;
	  delete [] this->Q4PerQ12;
	  delete [] this->NbrQ34Values;
	  this->Q1Value=NewQ1Value;
	  this->Q2Value=NewQ2Value;
	  this->Q3PerQ12=NewQ3PerQ12;
	  this->Q4PerQ12=NewQ4PerQ12;
	  this->NbrQ34Values=NewNbrQ34Values;
	  this->NbrQ12Indices=TmpNbrQ12Values;
	}
      // reduce size of table InteractionFactors to match new size, and copy contents
      delete [] this->InteractionFactors;
      this->InteractionFactors = new Complex[TmpNbrInteractionFactors];
      for (int i=0; i<TmpNbrInteractionFactors; ++i)
	this->InteractionFactors[i]=TmpInteractionFactors[i];
      delete [] TmpInteractionFactors;
    }

  cout << "Two-body hopping terms after hermitian symmetry: "<<this->CountTwoBodyInteractionTerms()<<endl;


  // diagonal terms are always the same... so we're done

  delete [] M;
  delete [] N;
  this->HermitianSymmetryFlag=true;
  return true;
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


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSource, vDestination);
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->LowLevelAddMultiplyFastMultiply(vSource, vDestination, firstComponent, LastComponent);
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
  return vDestination;
}

ComplexVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiplyFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int lastComponent)
{
  int* TmpIndexArray;
  ElementIndexType* TmpCoefficientIndexArray;
  double TmpRe, TmpIm;
  ElementIndexType TmpNbrRealInteraction;
  ElementIndexType TmpNbrComplexInteraction;
  Complex *TmpCPtr;
  int k = firstComponent;
  firstComponent -= this->PrecalculationShift;
  lastComponent -= this->PrecalculationShift;
  for (int i = firstComponent; i < lastComponent; ++i, ++k)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
      TmpRe = vSource[k].Re;
      TmpIm = vSource[k].Im;
      int Pos=0;
			//cout <<"i = "<< i<<endl;
			//cout <<"TmpNbrRealInteraction = "<<TmpNbrRealInteraction<<endl;
      for (; Pos < TmpNbrRealInteraction; ++Pos)
	{
		/*cout <<"TmpIndexArray[Pos] = "<<TmpIndexArray[Pos]<<endl;
		cout <<"vDestination.Re(TmpIndexArray[Pos]) = " << vDestination.Re(TmpIndexArray[Pos])<<endl;
		cout <<"TmpCoefficientIndexArray[Pos] = "<<TmpCoefficientIndexArray[Pos]<<endl;
		cout <<"TmpRe = " <<TmpRe<<endl;
		cout <<" RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]] = "<<RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]<<endl;*/
	  vDestination.Re(TmpIndexArray[Pos]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpRe;
	  vDestination.Im(TmpIndexArray[Pos]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpIm;
	}
      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos)
	{
	  TmpCPtr= &(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]]);
	  vDestination.Re(TmpIndexArray[Pos]) +=  TmpCPtr->Re*TmpRe-TmpCPtr->Im*TmpIm;
	  vDestination.Im(TmpIndexArray[Pos]) +=  TmpCPtr->Re*TmpIm+TmpCPtr->Im*TmpRe;		  
	}
      vDestination.Re(k) += this->HamiltonianShift * TmpRe;
      vDestination.Im(k) += this->HamiltonianShift * TmpIm;	      
    }
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

ComplexVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int* TmpIndexArray;
  ElementIndexType* TmpCoefficientIndexArray;
  double TmpRe, TmpIm;
  Complex *TmpCPtr;
  int TmpNbrRealInteraction;
  int TmpNbrComplexInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  int l =  PosMod + firstComponent + this->PrecalculationShift;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[Pos];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[Pos];
      TmpRe = vSource[l].Re;
      TmpIm = vSource[l].Im;
      int Pos2=0;
      for (; Pos2 < TmpNbrRealInteraction; ++Pos2)
	{
	  vDestination.Re(TmpIndexArray[Pos2]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*TmpRe;
	  vDestination.Im(TmpIndexArray[Pos2]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*TmpIm;
	}
      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos2)
	{
	  TmpCPtr= &(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]);
	  vDestination.Re(TmpIndexArray[Pos2]) +=  TmpCPtr->Re*TmpRe-TmpCPtr->Im*TmpIm;
	  vDestination.Im(TmpIndexArray[Pos2]) +=  TmpCPtr->Re*TmpIm+TmpCPtr->Im*TmpRe;		  
	}
      vDestination.Re(l) += this->HamiltonianShift * TmpRe;
      vDestination.Im(l) += this->HamiltonianShift * TmpIm;	      
      l += this->FastMultiplicationStep;
      ++Pos;
    }

  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	  this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + k, LastComponent, this->FastMultiplicationStep, vSource, vDestination);
      }
  delete TmpParticles;
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
  cout << "Attention: AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiplyDiskStorage must be defined" << endl;
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
  cout << "todo!" << endl;
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      Complex* Coefficient2 = new Complex [nbrVectors];
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles , firstComponent , LastComponent , 1 , vSources , vDestinations,nbrVectors);
      delete [] Coefficient2;
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->LowLevelMultipleAddMultiplyFastMultiply(vSources, vDestinations, nbrVectors , firstComponent , LastComponent);
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

ComplexVector* AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyFastMultiply(ComplexVector* vSources, ComplexVector * vDestinations,  int nbrVectors, int firstComponent, int lastComponent)
{
  int* TmpIndexArray;
  int Index;
  double TmpRealCoefficient;
  ElementIndexType* TmpCoefficientIndexArray;
  ElementIndexType TmpNbrRealInteraction;
  ElementIndexType TmpNbrComplexInteraction;
  Complex TmpCoefficient;
  Complex* Coefficient2 = new Complex [nbrVectors];
  int k = firstComponent;
  firstComponent -= this->PrecalculationShift;
  lastComponent -= this->PrecalculationShift;
  for (int i = firstComponent; i < lastComponent; ++i, ++k)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
      for (int l = 0; l < nbrVectors; ++l)
	{
	  Coefficient2[l] = vSources[l][k];
	  vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
	}
      int Pos=0;
      for (; Pos < TmpNbrRealInteraction; ++Pos)
	{
	  Index = TmpIndexArray[Pos];
	  TmpRealCoefficient = RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]];
	  for (int l = 0; l < nbrVectors; ++l)
	    vDestinations[l][Index] +=  TmpRealCoefficient * Coefficient2[l];
	}
      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos)
	{
	  Index = TmpIndexArray[Pos];
	  TmpCoefficient = ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]];
	  for (int l = 0; l < nbrVectors; ++l)
	    vDestinations[l][Index] +=  TmpCoefficient * Coefficient2[l];
	}
    }
  delete [] Coefficient2;
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

ComplexVector* AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors,  int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply"<<endl;
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
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyDiskStorage!"<<endl;
  return vDestinations;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination,  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      this->EvaluateMNOneBodyConjugateAddMultiplyComponent(TmpParticles,firstComponent,LastComponent,1, vSource, vDestination);
      for (int i = firstComponent; i < LastComponent; ++i)
	this->EvaluateMNTwoBodyConjugateAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->ConjugateLowLevelAddMultiplyFastMultiply(vSource,  vDestination, firstComponent, LastComponent);
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->ConjugateLowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
  return vDestination;
}

ComplexVector& AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelAddMultiplyFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										   int firstComponent, int lastComponent)
{
  int* TmpIndexArray;
  ElementIndexType* TmpCoefficientIndexArray;
  Complex TmpSum;
  ElementIndexType TmpNbrRealInteraction;
  ElementIndexType TmpNbrComplexInteraction;
  int k = firstComponent;
  firstComponent -= this->PrecalculationShift;
  lastComponent -= this->PrecalculationShift;
  for (int i = firstComponent; i < lastComponent; ++i, ++k)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
      TmpSum=0.0;
      int Pos=0;
      for (; Pos < TmpNbrRealInteraction; ++Pos)
	TmpSum +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*vSource[TmpIndexArray[Pos]];
      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos)
	TmpSum +=  Conj(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]])*vSource[TmpIndexArray[Pos]];
      vDestination[k] += TmpSum + this->HamiltonianShift * vSource[k];
    }
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

ComplexVector& AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
												int firstComponent, int nbrComponent)
{
  Complex TmpInteraction;
  int LastComponent = firstComponent + nbrComponent;
  Complex TmpSum=0.0;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  int* TmpIndexArray;
  ElementIndexType* TmpCoefficientIndexArray;
  Complex TmpC;
  int TmpNbrRealInteraction;
  int TmpNbrComplexInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  int l =  PosMod + firstComponent + this->PrecalculationShift;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[Pos];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[Pos];
      TmpSum=0.0;
      int Pos2=0;
      for (; Pos2 < TmpNbrRealInteraction; ++Pos2)
	TmpSum +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*vSource[TmpIndexArray[Pos2]];
      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos2)
	TmpSum +=  Conj(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos2]])*vSource[TmpIndexArray[Pos2]];
      vDestination[l] += this->HamiltonianShift * vSource[l];
      l += this->FastMultiplicationStep;
      ++Pos;
    }

  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	this->EvaluateMNOneBodyConjugateAddMultiplyComponent(TmpParticles, firstComponent + k, LastComponent, this->FastMultiplicationStep, vSource,vDestination);
	
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	  {
	    this->EvaluateMNTwoBodyConjugateAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	  }
      }
  delete TmpParticles;
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

ComplexVector& AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelAddMultiplyDiskStorage(ComplexVector& vSource, ComplexVector& vDestination, 
									   int firstComponent, int nbrComponent)
{
  cout << "Attention: AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelAddMultiplyDiskStorage must be defined" << endl;
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

ComplexVector* AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  
  if (this->FastMultiplicationFlag == false)
    {
      Complex * TmpCoefficients = new Complex [nbrVectors];
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      this->EvaluateMNOneBodyConjugateAddMultiplyComponent(TmpParticles, firstComponent,LastComponent,1, vSources, vDestinations , nbrVectors);
      for (int i =   firstComponent; i <LastComponent; i++)
	this->EvaluateMNTwoBodyConjugateAddMultiplyComponent(TmpParticles, i ,  vSources, vDestinations, nbrVectors,  TmpCoefficients);
      delete [] TmpCoefficients;
      delete TmpParticles;
    }	  
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->ConjugateLowLevelMultipleAddMultiplyFastMultiply(vSources,vDestinations, nbrVectors,firstComponent, LastComponent);
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->ConjugateLowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
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

ComplexVector* AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiplyFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors,  int firstComponent, int lastComponent)
{
  Complex TmpCoefficient;
  int Index;
  double TmpRealCoefficient;
	
	int* TmpIndexArray;
  ElementIndexType* TmpCoefficientIndexArray;
  ElementIndexType TmpNbrRealInteraction;
  ElementIndexType TmpNbrComplexInteraction;
  Complex* TmpSum = new Complex [nbrVectors];
	
  int k = firstComponent;
  firstComponent -= this->PrecalculationShift;
  lastComponent -= this->PrecalculationShift;
  for (int i = firstComponent; i < lastComponent; ++i, ++k)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
			
      for (int l = 0; l < nbrVectors; ++l)
				TmpSum[l] = 0.0;
			
      int Pos=0;
      for (; Pos < TmpNbrRealInteraction; ++Pos)
	{
	  Index = TmpIndexArray[Pos];
	  TmpRealCoefficient = RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]];
	  for (int l = 0; l < nbrVectors; ++l)
	    TmpSum[l] +=  TmpRealCoefficient * vSources[l][Index];
	}
      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos)
	{
	  Index = TmpIndexArray[Pos];
	  TmpCoefficient = Conj(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]]);
	  for (int l = 0; l < nbrVectors; ++l)
	    TmpSum[l] +=  TmpCoefficient * vSources[l][Index];
	}
      for (int l = 0; l < nbrVectors; ++l)
	vDestinations[l][i] += TmpSum[l] + this->HamiltonianShift * vSources[l][i];
    }
  delete [] TmpSum;  
	
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

ComplexVector* AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
													int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply"<<endl;
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

ComplexVector* AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiplyDiskStorage(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
												int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiplyDiskStorage!"<<endl;
  return vDestinations;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
      this->HermitianEvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSource, vDestination);
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->HermitianLowLevelAddMultiplyFastMultiply(vSource, vDestination, firstComponent, LastComponent);
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->HermitianLowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->HermitianLowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
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

ComplexVector& AbstractQHEOnLatticeHamiltonian::HermitianLowLevelAddMultiplyFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										   int firstComponent, int lastComponent)
{
  int* TmpIndexArray;
  ElementIndexType* TmpCoefficientIndexArray;
  Complex TmpElement;
  ElementIndexType TmpNbrRealInteraction;
  ElementIndexType TmpNbrComplexInteraction;
  Complex TmpSum;
  int k = firstComponent;
  firstComponent -= this->PrecalculationShift;
  lastComponent -= this->PrecalculationShift;
  // cout << "firstComponent="<<firstComponent<<", lastComponent="<< lastComponent<<endl;
  for (int i = firstComponent; i < lastComponent; ++i, ++k)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
      TmpElement = vSource[k];
      TmpSum = 0.0;
      int Pos=0;
      for (; Pos < TmpNbrRealInteraction; ++Pos)
	{
	  vDestination[TmpIndexArray[Pos]] +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpElement;
	  TmpSum +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]] * vSource[TmpIndexArray[Pos]];
	}
      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos)
	{
	  vDestination[TmpIndexArray[Pos]] +=  ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpElement;
	  TmpSum +=  Conj(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]]) * vSource[TmpIndexArray[Pos]];
	}
      vDestination[k] += TmpSum + this->HamiltonianShift * TmpElement;
    }
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

ComplexVector& AbstractQHEOnLatticeHamiltonian::HermitianLowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
												int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  int* TmpIndexArray;
  ElementIndexType* TmpCoefficientIndexArray;
  Complex TmpInteraction, TmpC;
  Complex TmpSum;
  int TmpNbrRealInteraction;
  int TmpNbrComplexInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  int l =  PosMod + firstComponent + this->PrecalculationShift;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[Pos];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[Pos];
      TmpSum=0.0;
      TmpC = vSource[l];
      int Pos2=0;
      for (; Pos2 < TmpNbrRealInteraction; ++Pos2)
	{
	  vDestination[TmpIndexArray[Pos2]] +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*TmpC;
	  TmpSum += RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*vSource[TmpIndexArray[Pos2]];
	}
      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos2)
	{
	  TmpInteraction= ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos2]];
	  vDestination[TmpIndexArray[Pos2]] +=  TmpInteraction*TmpC;
	  TmpSum +=  Conj(TmpInteraction) * vSource[TmpIndexArray[Pos2]];
	}
      vDestination[l] += TmpSum + this->HamiltonianShift * TmpC;
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	this->HermitianEvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + k, LastComponent,this->FastMultiplicationStep, vSource, vDestination);
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	  {
	    this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	  }
      }
  delete TmpParticles;
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

ComplexVector& AbstractQHEOnLatticeHamiltonian::HermitianLowLevelAddMultiplyDiskStorage(ComplexVector& vSource, ComplexVector& vDestination, 
											int firstComponent, int nbrComponent)
{
  cout << "Attention: AbstractQHEOnLatticeHamiltonian::HermitianLowLevelAddMultiplyDiskStorage must be defined" << endl;
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

ComplexVector* AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      Complex* TmpCoefficients = new Complex [nbrVectors];
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      this->HermitianEvaluateMNOneBodyAddMultiplyComponent(TmpParticles,firstComponent, LastComponent, 1, vSources, vDestinations,nbrVectors);
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles,i , vSources,  vDestinations,  nbrVectors, TmpCoefficients);
	}
      delete [] TmpCoefficients;
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->HermitianLowLevelMultipleAddMultiplyFastMultiply( vSources, vDestinations, nbrVectors, firstComponent, LastComponent);
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->HermitianLowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
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

ComplexVector* AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiplyFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int lastComponent)
{
  int* TmpIndexArray;
  int Index;
  double TmpRealCoefficient;
  Complex TmpCoefficient;
  ElementIndexType* TmpCoefficientIndexArray;
  ElementIndexType TmpNbrRealInteraction;
  ElementIndexType TmpNbrComplexInteraction;
  Complex* Coefficient2 = new Complex [nbrVectors];
  Complex* TmpSum = new Complex[nbrVectors];
  for (int l = 0; l < nbrVectors; ++l)
    TmpSum[l] = 0.0;
  int k = firstComponent;
  firstComponent -= this->PrecalculationShift;
  lastComponent -= this->PrecalculationShift;
  for (int i = firstComponent; i < lastComponent; ++i, ++k)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
      for (int l = 0; l < nbrVectors; ++l)
	{
	  TmpSum[l] = 0.0;
	  Coefficient2[l] = vSources[l][k];
	}
      int Pos=0;
      for (; Pos < TmpNbrRealInteraction; ++Pos)
	{
	  Index = TmpIndexArray[Pos];
	  TmpRealCoefficient = RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]];
	  for (int l = 0; l < nbrVectors; ++l)
	    {
	      vDestinations[l][Index] +=  TmpRealCoefficient * Coefficient2[l];
	      TmpSum[l] += TmpRealCoefficient * vSources[l][Index];
	    }
	}
      for (int j = 0; j < TmpNbrComplexInteraction; ++j, ++Pos)
	{
	  Index = TmpIndexArray[Pos];
	  TmpCoefficient = ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]];
	  for (int l = 0; l < nbrVectors; ++l)
	    vDestinations[l][Index] +=  TmpCoefficient * Coefficient2[l];
	  TmpCoefficient.Conjugate();
	  for (int l = 0; l < nbrVectors; ++l)
	    TmpSum[l] += TmpCoefficient * vSources[l][Index];
	}
      for (int l = 0; l < nbrVectors; ++l)
	vDestinations[l][k] += TmpSum[l] + this->HamiltonianShift * Coefficient2[l];
    }
  delete [] Coefficient2;
  delete [] TmpSum;
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

ComplexVector* AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
													int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply"<<endl;
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

ComplexVector* AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiplyDiskStorage(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
												int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiplyDiskStorage!"<<endl;
  return vDestinations;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSource, vDestination);
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->LowLevelAddMultiplyFastMultiply(vSource, vDestination, firstComponent, LastComponent);
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

RealVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiplyFastMultiply(RealVector& vSource, RealVector& vDestination, 
										   int firstComponent, int lastComponent)
{
  int* TmpIndexArray;
  ElementIndexType* TmpCoefficientIndexArray;
  double TmpRe;
  ElementIndexType TmpNbrRealInteraction;
  int k = firstComponent;
  firstComponent -= this->PrecalculationShift;
  lastComponent -= this->PrecalculationShift;
  for (int i = firstComponent; i < lastComponent; ++i, ++k)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
      TmpRe = vSource[k];
      int Pos=0;
      for (; Pos < TmpNbrRealInteraction; ++Pos)
	{
	  vDestination[TmpIndexArray[Pos]] +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpRe;
	}
      vDestination[k] += this->HamiltonianShift * TmpRe;
    }
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

RealVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
										   int firstComponent, int nbrComponent)
{
  double TmpInteraction;
  int LastComponent = firstComponent + nbrComponent;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  int* TmpIndexArray;
  ElementIndexType* TmpCoefficientIndexArray;
  double TmpRe;
  int TmpNbrRealInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  int l =  PosMod + firstComponent + this->PrecalculationShift;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[Pos];
      TmpRe = vSource[l];
      int Pos2=0;
      for (; Pos2 < TmpNbrRealInteraction; ++Pos2)
	{
	  vDestination[TmpIndexArray[Pos2]] +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*TmpRe;
	}
      vDestination[l] += this->HamiltonianShift * TmpRe;
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	  this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + k, LastComponent, this->FastMultiplicationStep, vSource, vDestination);
      }
  delete TmpParticles;
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

RealVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
									    int firstComponent, int nbrComponent)
{
  cout << "Attention: AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiplyDiskStorage must be defined" << endl;
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
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      double* Coefficient2 = new double [nbrVectors];
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles , firstComponent , LastComponent , 1 , vSources , vDestinations,nbrVectors);
      delete [] Coefficient2;
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, LastComponent);
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


RealVector* AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										     int firstComponent, int lastComponent)
{
  int* TmpIndexArray;
  int Index;
  double TmpRealCoefficient;
  ElementIndexType* TmpCoefficientIndexArray;
  ElementIndexType TmpNbrRealInteraction;
  double* Coefficient2 = new double [nbrVectors];
  int k = firstComponent;
  firstComponent -= this->PrecalculationShift;
  lastComponent -= this->PrecalculationShift;
  for (int i = firstComponent; i < lastComponent; ++i, ++k)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
      for (int l = 0; l < nbrVectors; ++l)
	{
	  Coefficient2[l] = vSources[l][k];
	  vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
	}
      int Pos=0;
      for (; Pos < TmpNbrRealInteraction; ++Pos)
	{
	  Index = TmpIndexArray[Pos];
	  TmpRealCoefficient = RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]];
	  for (int l = 0; l < nbrVectors; ++l)
	    vDestinations[l][Index] +=  TmpRealCoefficient * Coefficient2[l];
	}
    }
  delete [] Coefficient2;
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

RealVector* AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											   int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply"<<endl;
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

RealVector* AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										    int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyDiskStorage!"<<endl;
  return vDestinations;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelAddMultiply(RealVector& vSource,RealVector& vDestination,  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      this->EvaluateMNOneBodyConjugateAddMultiplyComponent(TmpParticles,firstComponent,LastComponent,1, vSource, vDestination);
      for (int i = firstComponent; i < LastComponent; ++i)
	this->EvaluateMNTwoBodyConjugateAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->ConjugateLowLevelAddMultiplyFastMultiply(vSource, vDestination,firstComponent, LastComponent);
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->ConjugateLowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
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

RealVector& AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelAddMultiplyFastMultiply(RealVector& vSource, RealVector& vDestination, 
										   int firstComponent, int lastComponent)
{
  int* TmpIndexArray;
  ElementIndexType* TmpCoefficientIndexArray;
  double TmpSum;
  ElementIndexType TmpNbrRealInteraction;
  int k = firstComponent;
  firstComponent -= this->PrecalculationShift;
  lastComponent -= this->PrecalculationShift;
  for (int i = firstComponent; i < lastComponent; ++i, ++k)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
      TmpSum=0.0;
      int Pos=0;
      for (; Pos < TmpNbrRealInteraction; ++Pos)
	TmpSum +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*vSource[TmpIndexArray[Pos]];
      vDestination[k] += TmpSum + this->HamiltonianShift * vSource[k];
    }
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

RealVector& AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
											     int firstComponent, int nbrComponent)
{
  double TmpInteraction;
  int LastComponent = firstComponent + nbrComponent;
  double TmpSum=0.0;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  int* TmpIndexArray;
  ElementIndexType* TmpCoefficientIndexArray;
  int TmpNbrRealInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  int l =  PosMod + firstComponent + this->PrecalculationShift;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[Pos];
      TmpSum=0.0;
      int Pos2=0;
      for (; Pos2 < TmpNbrRealInteraction; ++Pos2)
	TmpSum +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*vSource[TmpIndexArray[Pos2]];
      vDestination[l] += this->HamiltonianShift * vSource[l];
      l += this->FastMultiplicationStep;
      ++Pos;
    }

  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	this->EvaluateMNOneBodyConjugateAddMultiplyComponent(TmpParticles, firstComponent + k, LastComponent, this->FastMultiplicationStep, vSource,vDestination);
	
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	  {
	    this->EvaluateMNTwoBodyConjugateAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	  }
      }
  delete TmpParticles;
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

RealVector& AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
									   int firstComponent, int nbrComponent)
{
  cout << "Attention: AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelAddMultiplyDiskStorage must be defined" << endl;
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

RealVector* AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      double * TmpCoefficients = new double [nbrVectors];
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      this->EvaluateMNOneBodyConjugateAddMultiplyComponent(TmpParticles, firstComponent,LastComponent,1, vSources, vDestinations , nbrVectors);
      for (int i =   firstComponent; i <LastComponent; i++)
	this->EvaluateMNTwoBodyConjugateAddMultiplyComponent(TmpParticles, i ,  vSources, vDestinations, nbrVectors,  TmpCoefficients);
      delete [] TmpCoefficients;
      delete TmpParticles;
    }	  
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, LastComponent);
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->ConjugateLowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
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

RealVector* AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiplyFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											      int firstComponent, int lastComponent)
{
  int* TmpIndexArray;
  int Index;
  double TmpRealCoefficient;
  ElementIndexType* TmpCoefficientIndexArray;
  ElementIndexType TmpNbrRealInteraction;
  double* TmpSum = new double [nbrVectors];
  int k = firstComponent;
  firstComponent -= this->PrecalculationShift;
  lastComponent -= this->PrecalculationShift;
  for (int i = firstComponent; i < lastComponent; ++i, ++k)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
      for (int l = 0; l < nbrVectors; ++l)
	TmpSum[l] = 0.0;
      int Pos=0;
      for (; Pos < TmpNbrRealInteraction; ++Pos)
	{
	  Index = TmpIndexArray[Pos];
	  TmpRealCoefficient = RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]];
	  for (int l = 0; l < nbrVectors; ++l)
	    TmpSum[l] +=  TmpRealCoefficient * vSources[l][Index];
	}
      for (int l = 0; l < nbrVectors; ++l)
	TmpSum[l] += TmpSum[l] + this->HamiltonianShift * vSources[l][k];
    }
  delete [] TmpSum;
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

RealVector* AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											   int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply"<<endl;
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

RealVector* AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										   int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiplyDiskStorage!"<<endl;
  return vDestinations;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnLatticeHamiltonian::HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
      this->HermitianEvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSource, vDestination);
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->HermitianLowLevelAddMultiplyFastMultiply(vSource, vDestination, firstComponent, LastComponent);
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->HermitianLowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->HermitianLowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
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

RealVector& AbstractQHEOnLatticeHamiltonian::HermitianLowLevelAddMultiplyFastMultiply(RealVector& vSource, RealVector& vDestination, 
										   int firstComponent, int lastComponent)
{
  int* TmpIndexArray;
  ElementIndexType* TmpCoefficientIndexArray;
  double TmpElement;
  ElementIndexType TmpNbrRealInteraction;
  ElementIndexType TmpNbrComplexInteraction;
  double TmpSum;
  int k = firstComponent;
  firstComponent -= this->PrecalculationShift;
  lastComponent -= this->PrecalculationShift;
  for (int i = firstComponent; i < lastComponent; ++i, ++k)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
      TmpElement = vSource[k];
      TmpSum = 0.0;
      int Pos=0;
      for (; Pos < TmpNbrRealInteraction; ++Pos)
	{
	  vDestination[TmpIndexArray[Pos]] +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpElement;
	  TmpSum +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]] * vSource[TmpIndexArray[Pos]];
	}
      vDestination[k] += TmpSum + this->HamiltonianShift * TmpElement;
    }
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

RealVector& AbstractQHEOnLatticeHamiltonian::HermitianLowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
											     int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  int* TmpIndexArray;
  ElementIndexType* TmpCoefficientIndexArray;
  double TmpC;
  double TmpSum;
  int TmpNbrRealInteraction;
  int TmpNbrComplexInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  int l =  PosMod + firstComponent + this->PrecalculationShift;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[Pos];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[Pos];
      TmpSum=0.0;
      TmpC = vSource[l];
      int Pos2=0;
      for (; Pos2 < TmpNbrRealInteraction; ++Pos2)
	{
	  vDestination[TmpIndexArray[Pos2]] +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*TmpC;
	  TmpSum += RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*vSource[TmpIndexArray[Pos2]];
	}
      vDestination[l] += TmpSum + this->HamiltonianShift * TmpC;
      l += this->FastMultiplicationStep;
      ++Pos;
    }

  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
	
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	this->HermitianEvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + k, LastComponent,this->FastMultiplicationStep, vSource, vDestination);
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	  {
	    this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
	  }
      }
  delete TmpParticles;
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

RealVector& AbstractQHEOnLatticeHamiltonian::HermitianLowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
									   int firstComponent, int nbrComponent)
  {
  cout << "Attention: AbstractQHEOnLatticeHamiltonian::HermitianLowLevelAddMultiplyDiskStorage must be defined" << endl;
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

RealVector* AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      double* TmpCoefficients = new double [nbrVectors];
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      this->HermitianEvaluateMNOneBodyAddMultiplyComponent(TmpParticles,firstComponent, LastComponent, 1, vSources, vDestinations,nbrVectors);
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles,i , vSources,  vDestinations,  nbrVectors, TmpCoefficients);
	}
      delete [] TmpCoefficients;
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  this->HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, LastComponent);
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->HermitianLowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
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

RealVector* AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiplyFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											   int firstComponent, int lastComponent)
{
  int* TmpIndexArray;
  int Index;
  double TmpRealCoefficient;
  ElementIndexType* TmpCoefficientIndexArray;
  ElementIndexType TmpNbrRealInteraction;
  double* Coefficient2 = new double [nbrVectors];
  double* TmpSum = new double[nbrVectors];
  for (int l = 0; l < nbrVectors; ++l)
    TmpSum[l] = 0.0;
  int k = firstComponent;
  firstComponent -= this->PrecalculationShift;
  lastComponent -= this->PrecalculationShift;
  for (int i = firstComponent; i < lastComponent; ++i, ++k)
    {
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
      for (int l = 0; l < nbrVectors; ++l)
	{
	  TmpSum[l] = 0.0;
	  Coefficient2[l] = vSources[l][k];
	}
      int Pos=0;
      for (; Pos < TmpNbrRealInteraction; ++Pos)
	{
	  Index = TmpIndexArray[Pos];
	  TmpRealCoefficient = RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]];
	  for (int l = 0; l < nbrVectors; ++l)
	    {
	      vDestinations[l][Index] +=  TmpRealCoefficient * Coefficient2[l];
	      TmpSum[l] += TmpRealCoefficient * vSources[l][Index];
	    }
	}
      for (int l = 0; l < nbrVectors; ++l)
	vDestinations[l][k] += TmpSum[l] + this->HamiltonianShift * Coefficient2[l];
    }
  delete [] Coefficient2;	
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

RealVector* AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											   int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply"<<endl;
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

RealVector* AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											     int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiplyDiskStorage!"<<endl;
  return vDestinations;
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

// get the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
// nbrThreads = number of threads requested
// segmentIndices = array returning the reference to an array of the first index of each of the segments
//
bool AbstractQHEOnLatticeHamiltonian::GetLoadBalancing(int nbrTasks, long* &segmentIndices)
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;

  if (((NbrRealInteractionPerComponent!=0)||(NbrComplexInteractionPerComponent!=0))&&(this->FastMultiplicationStep!=0))
    {
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      if ((LoadBalancingArray==0)||(NbrBalancedTasks!=nbrTasks))
	{
	  if (LoadBalancingArray!=0)
	    delete [] LoadBalancingArray;
	  this->LoadBalancingArray = new long[nbrTasks+1];
	  long *SegmentSize = new long[nbrTasks];
	  this->NbrBalancedTasks=nbrTasks;
	  long TmpNbrElement=0;
	  for (int i=0; i<ReducedSpaceDimension; ++i)
	    TmpNbrElement+=NbrRealInteractionPerComponent[i]+NbrComplexInteractionPerComponent[i];
	  long TmpNbrPerSegment = TmpNbrElement/nbrTasks;
	  TmpNbrElement=0;
	  int Pos=0;
	  this->LoadBalancingArray[0]=MinIndex;
	  for (int i=0; i<ReducedSpaceDimension; ++i)
	    {
	      TmpNbrElement+=NbrRealInteractionPerComponent[i]+NbrComplexInteractionPerComponent[i];
	      if (TmpNbrElement>TmpNbrPerSegment)
		{
		  SegmentSize[Pos]=TmpNbrElement;
		  LoadBalancingArray[Pos+1]=MinIndex+i*this->FastMultiplicationStep;
		  TmpNbrElement=0;
		  ++Pos;
		}
	    }
	  LoadBalancingArray[nbrTasks]=MaxIndex+1;
	  SegmentSize[nbrTasks-1]=TmpNbrElement;
	  
	  cout << "LoadBalancingArray=[ ("<<LoadBalancingArray[1]-LoadBalancingArray[0]<<", "<<SegmentSize[0]<<")";
	  for (int i=1; i<nbrTasks; ++i)
	    cout <<" ("<<LoadBalancingArray[i+1]-LoadBalancingArray[i]<<", "<<SegmentSize[i]<<")";
	  cout << "]"<< endl;
	  delete [] SegmentSize;
	}
    }
  else
    {
      if ((LoadBalancingArray==0)||(NbrBalancedTasks!=nbrTasks))
	{
	  if (LoadBalancingArray!=0)
	    delete [] LoadBalancingArray;
	  this->LoadBalancingArray = new long[nbrTasks+1];
	  
	  int Step = EffectiveHilbertSpaceDimension / nbrTasks;
	  this->LoadBalancingArray[0]=MinIndex;
	  for (int i=0; i<nbrTasks; ++i)
	    LoadBalancingArray[i]=MinIndex+i*Step;
	  LoadBalancingArray[nbrTasks]=MaxIndex+1;
	}
    }
  segmentIndices=LoadBalancingArray;
  return true;
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
  this->NbrRealInteractionPerComponent = new ElementIndexType [EffectiveHilbertSpaceDimension];
  this->NbrComplexInteractionPerComponent = new ElementIndexType [EffectiveHilbertSpaceDimension];   
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

  const char* ComplexFileName = "complex-nbr-interaction.dat";
  const char* RealFileName = "real-nbr-interaction.dat";

  int *TmpInteractionPerComponentIndex=NULL;
  if (this->Architecture->LoadOptimizedTypicalRange (TmpInteractionPerComponentIndex, MinIndex, MaxIndex))
    {
      delete [] TmpInteractionPerComponentIndex;
      this->PrecalculationShift = (int) MinIndex;
      EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;

#ifdef ABSTRACTQHEONLATTICEHAMILTONIAN_LONGINDEX
      this->Architecture->LoadArray(ComplexFileName, NbrComplexInteractionPerComponent);
      this->Architecture->LoadArray(RealFileName, NbrRealInteractionPerComponent);
      this->RealInteractionCoefficients.ReadArray("real-entries.dat");
      this->ComplexInteractionCoefficients.ReadArray("complex-entries.dat");
#else
      cout << "Error: loading of load balancing not implemented for short arrays"<<endl;
      exit(1);
#endif
    }
  else
    {  
#ifdef ABSTRACTQHEONLATTICEHAMILTONIAN_SORTED
      QHEParticlePrecalculationOperationWithMatrixElements Operation(this, true, /* tolerance = */ LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD);
      Operation.ApplyOperation(this->Architecture);
      Operation.GetMatrixElements(this->RealInteractionCoefficients, this->ComplexInteractionCoefficients);
#else
      QHEParticlePrecalculationOperation Operation(this);
      Operation.ApplyOperation(this->Architecture);
#endif
  
      int *TmpInteractionPerComponentIndex = new int[EffectiveHilbertSpaceDimension];
      for (long i=0; i<EffectiveHilbertSpaceDimension; ++i)
	TmpInteractionPerComponentIndex[i] = this->NbrRealInteractionPerComponent[i] + this->NbrComplexInteractionPerComponent[i];
      // adapt load balancing for memory:
      if (this->Architecture->GetOptimizedTypicalRange(TmpInteractionPerComponentIndex, MinIndex, MaxIndex, this->RealInteractionCoefficients, this->ComplexInteractionCoefficients) == true)
	{
	  long OldEffectiveHilbertSpaceDimension = EffectiveHilbertSpaceDimension;
	  this->PrecalculationShift = (int) MinIndex;
	  EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
#ifdef ABSTRACTQHEONLATTICEHAMILTONIAN_LONGINDEX
	  // balance original arrays
	  this->Architecture->RebalanceArray(NbrComplexInteractionPerComponent, ComplexFileName);
	  this->Architecture->RebalanceArray(NbrRealInteractionPerComponent, RealFileName);
#else
	  delete [] TmpInteractionPerComponentIndex;
	  TmpInteractionPerComponentIndex = new int[OldEffectiveHilbertSpaceDimension];
	  for (long i=0; i<OldEffectiveHilbertSpaceDimension; ++i)
	    TmpInteractionPerComponentIndex[i] = (int)NbrComplexInteractionPerComponent[i];
	  this->Architecture->RebalanceArray(TmpInteractionPerComponentIndex, ComplexFileName);
	  delete [] this->NbrComplexInteractionPerComponent;
	  this->NbrComplexInteractionPerComponent = new ElementIndexType[EffectiveHilbertSpaceDimension];
	  for (long i=0; i<EffectiveHilbertSpaceDimension; ++i)
	    NbrComplexInteractionPerComponent[i] = (ElementIndexType)TmpInteractionPerComponentIndex[i];
      
	  delete [] TmpInteractionPerComponentIndex;
	  TmpInteractionPerComponentIndex = new int[OldEffectiveHilbertSpaceDimension];
	  for (long i=0; i<OldEffectiveHilbertSpaceDimension; ++i)
	    TmpInteractionPerComponentIndex[i] = (int)NbrRealInteractionPerComponent[i];
	  this->Architecture->RebalanceArray(TmpInteractionPerComponentIndex, RealFileName);
	  delete [] this->NbrRealInteractionPerComponent;
	  this->NbrRealInteractionPerComponent = new ElementIndexType[EffectiveHilbertSpaceDimension];
	  for (long i=0; i<EffectiveHilbertSpaceDimension; ++i)
	    NbrRealInteractionPerComponent[i] = (ElementIndexType)TmpInteractionPerComponentIndex[i];
#endif
	  // cout << "distributed calculations successfully reoptimized" << endl;
	}
      delete [] TmpInteractionPerComponentIndex;

    } // done precalculating memory size
  
  // in case of no memory, delete allocated memory.
  if (allowedMemory == 0l)
    {
      delete[] this->NbrRealInteractionPerComponent;
      delete[] this->NbrComplexInteractionPerComponent;
      this->NbrRealInteractionPerComponent = 0;
      this->NbrComplexInteractionPerComponent = 0;
      this->RealInteractionCoefficients.Empty();
      this->ComplexInteractionCoefficients.Empty();
      return 0l;
    }

#ifdef ABSTRACTQHEONLATTICEHAMILTONIAN_LONGINDEX
  // fix order of matrix elements, in case any missing entries are found.
  this->RealInteractionCoefficients.FixOrder();
  this->ComplexInteractionCoefficients.FixOrder();
#endif


  long Memory = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    {
      Memory += this->NbrRealInteractionPerComponent[i];
      Memory += this->NbrComplexInteractionPerComponent[i];
    }
  
  cout << "nbr interaction = " << Memory << endl;

  cout << "Nbr unique real elements: "<<this->RealInteractionCoefficients.GetNbrElements()<<endl;
  cout << "Nbr unique complex elements: "<<this->ComplexInteractionCoefficients.GetNbrElements()<<endl;
  cout << "Number of interaction factors: "<<this->NbrInteractionFactors<<" & diagonal terms: "<< this->NbrDiagonalInteractionFactors <<endl;

  // memory requirement, ignoring the actual storage size of the values of matrix
  // elements, which is assumed small (maybe need to add an estimate, at least)
  long TmpMemory = AllowedMemory - (2*sizeof (ElementIndexType) + sizeof (int*) + sizeof(ElementIndexType*)) * EffectiveHilbertSpaceDimension;
  cout << "of which can be stored: "<<(TmpMemory / ((int) (sizeof (int) + sizeof(ElementIndexType))))<<endl;
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(ElementIndexType)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(ElementIndexType)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if (this->Particles->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  // memory requirement, ignoring the actual storage size of the values of matrix
	  // elements, which is assumed small (maybe need to add an estimate, at least, again!)
	  TmpMemory = AllowedMemory - (2*sizeof (ElementIndexType) + sizeof (int*) + sizeof(ElementIndexType*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	    {
	      Memory += this->NbrRealInteractionPerComponent[i];
	      Memory += this->NbrComplexInteractionPerComponent[i];
	    }	  
	}
      
      Memory = ((2*sizeof (ElementIndexType) + sizeof (int*) + sizeof(ElementIndexType*)) * ReducedSpaceDimension) + (Memory * (sizeof (int) + sizeof(ElementIndexType)));
      
      if (this->DiskStorageFlag == false)
	{
	  int TotalReducedSpaceDimension = ReducedSpaceDimension;
	  ElementIndexType* TmpNbrRealInteractionPerComponent = new ElementIndexType [TotalReducedSpaceDimension];
	  ElementIndexType* TmpNbrComplexInteractionPerComponent = new ElementIndexType [TotalReducedSpaceDimension];	  
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
      Memory = ((2*sizeof (ElementIndexType) + sizeof (int*) + sizeof(ElementIndexType*)) * EffectiveHilbertSpaceDimension) + (Memory * (sizeof (int) + sizeof(ElementIndexType)));
      this->FastMultiplicationStep = 1;
      // we could safely delete all interaction factors at this point.

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
  long Memory = 0;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  int LastComponent =  nbrComponent + firstComponent;
  this->EvaluateMNOneBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);
  this->EvaluateMNTwoBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);
  
  delete TmpParticles;
  return Memory;
}


// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that have to be precalcualted
// realInteractionCoefficients = reference on an object collecting unique real matrix elements for this thread
// complexInteractionCoefficients = reference on an object collecting unique complex matrix elements for this thread
// return value = number of non-zero matrix elements
//
long AbstractQHEOnLatticeHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int nbrComponent, SortedRealUniqueArray &realInteractionCoefficients, SortedComplexUniqueArray &complexInteractionCoefficients)
{
  long Memory = 0;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  int LastComponent =  nbrComponent + firstComponent;
  this->EvaluateMNOneBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory, realInteractionCoefficients, complexInteractionCoefficients);
  this->EvaluateMNTwoBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory, realInteractionCoefficients, complexInteractionCoefficients);
  // sort all entries when done.
  realInteractionCoefficients.SortEntries();
  complexInteractionCoefficients.SortEntries();
  
  delete TmpParticles;
  return Memory;
}

// enable fast multiplication algorithm
//

void AbstractQHEOnLatticeHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficientIndex = new ElementIndexType* [ReducedSpaceDimension];
  
  // allocate all memory at the outset:
  for (int i = 0; i < ReducedSpaceDimension; ++i)
    {
      //cout <<"i = "<< i << this->NbrRealInteractionPerComponent[i]<<" "<<this->NbrComplexInteractionPerComponent[i]<<endl;
      this->InteractionPerComponentIndex[i] = new int [this->NbrRealInteractionPerComponent[i] + this->NbrComplexInteractionPerComponent[i]];
      this->InteractionPerComponentCoefficientIndex[i] = new ElementIndexType [this->NbrRealInteractionPerComponent[i]
									     +this->NbrComplexInteractionPerComponent[i]];
    }
  
  QHEParticlePrecalculationOperation Operation(this, false);
  Operation.ApplyOperation(this->Architecture);
      
  cout << "Nbr distinct matrix elements: "<<RealInteractionCoefficients.GetNbrElements()<<" real, "
       << ComplexInteractionCoefficients.GetNbrElements()<<" complex"<<endl;
   
  this->FastMultiplicationFlag = true;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
}

// enable fast multiplication algorithm using on disk cache 
//
// fileName = prefix of the name of the file where temporary matrix elements will be stored

void AbstractQHEOnLatticeHamiltonian::EnableFastMultiplicationWithDiskStorage(char* fileName)
{
  cout << "Using non-defined function EnableFastMultiplicationWithDiskStorage!"<<endl;
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted

void AbstractQHEOnLatticeHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{
  int LastComponent = nbrComponent + firstComponent;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  long Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      int PosR = 0;
      int PosC = this->NbrRealInteractionPerComponent[Pos];
      this->EvaluateMNOneBodyFastMultiplicationComponent(TmpParticles, i, this->InteractionPerComponentIndex[Pos], this->InteractionPerComponentCoefficientIndex[Pos], PosR,PosC);
      this->EvaluateMNTwoBodyFastMultiplicationComponent(TmpParticles, i, this->InteractionPerComponentIndex[Pos], 
							 this->InteractionPerComponentCoefficientIndex[Pos], PosR,PosC);
      ++Pos;
    }
  
  
  //   long TotalPos = ((firstComponent - this->PrecalculationShift - 1) / this->FastMultiplicationStep) + 1;
  //   int InitalPos = ((firstComponent - 1) / this->FastMultiplicationStep) + 1;
  //   InitalPos *= this->FastMultiplicationStep;
  //   for (int i = InitalPos; i < LastComponent; i += this->FastMultiplicationStep)
  //     {
  //       this->EvaluateFastMultiplicationComponent(TmpParticles, i, this->InteractionPerComponentIndex[TotalPos], 
  // 						this->InteractionPerComponentCoefficientIndex[TotalPos], TotalPos);
  //     }

  // cout << "Distinct matrix elements in PartialEnableFastMultiplication: real = "<< RealInteractionCoefficients.GetNbrElements() << " in array ("<<&RealInteractionCoefficients<<"), complex = " << ComplexInteractionCoefficients.GetNbrElements() << " in array ("<<&ComplexInteractionCoefficients<<")"<<endl;
  
  delete TmpParticles;
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
      File.write((char*) this->NbrRealInteractionPerComponent, sizeof(ElementIndexType) * Tmp);
      File.write((char*) this->NbrComplexInteractionPerComponent, sizeof(ElementIndexType) * Tmp);
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * (this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i]));
	}
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentCoefficientIndex[i]), sizeof(ElementIndexType) * (this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i]));	  
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
  this->NbrRealInteractionPerComponent = new ElementIndexType [Tmp];
  this->NbrComplexInteractionPerComponent = new ElementIndexType [Tmp];
  File.read((char*) this->NbrRealInteractionPerComponent, sizeof(ElementIndexType) * Tmp);
  File.read((char*) this->NbrComplexInteractionPerComponent, sizeof(ElementIndexType) * Tmp);

  this->InteractionPerComponentIndex = new int* [Tmp];
  this->InteractionPerComponentCoefficientIndex = new ElementIndexType* [Tmp];
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [(this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i])];
      File.read((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * (this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i]));	  
    }
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentCoefficientIndex[i]=new ElementIndexType[(this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i])];
      File.read((char*) (this->InteractionPerComponentCoefficientIndex[i]), sizeof(ElementIndexType) * (this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i]));	  
    }
  RealInteractionCoefficients.ReadArray(File);
  ComplexInteractionCoefficients.ReadArray(File);
   
  File.close();
  this->FastMultiplicationFlag = true;
  return true;
}

