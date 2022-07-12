////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of two dimension Heisenberg model                  //
//                             and 2d translations                            //
//                                                                            //
//                        last modification : 13/05/2018                      //
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


#include "Hamiltonian/TwoDimensionalHeisenbergAnd2DTranslationHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/RealVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Architecture/ArchitectureOperation/GenericHamiltonianPrecalculationOperation.h"
#include "Output/MathematicaOutput.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::TwoDimensionalHeisenbergAnd2DTranslationHamiltonian()
{
  this->Chain = 0;
  this->XMomentum = 0;
  this->YMomentum = 0;
  this->NbrSpinX = 0;
  this->NbrSpinY = 0;
  this->NbrSpin = 0;
  this->JFactor = 0.0;
  this->JzFactor = 0.0;
  this->HalfJFactor = 0.0;
  this->SzSzContributions = 0;
  this->FastMultiplicationFlag = false;
  this->FastMultiplicationStep = 0;
  this->NbrInteractionPerComponent = 0;
  this->NbrBalancedTasks = 0;
  this->LoadBalancingArray = 0;
  this->InteractionPerComponentIndex = 0;
  this->InteractionPerComponentCoefficient = 0;
  this->HermitianSymmetryFlag = true;
  this->Architecture = 0;
  this->PrecalculationShift = 0;
}

// constructor from default data
//
// chain = pointer to Hilbert space of the associated system
// xMomentum = momentum along the x direction
// nbrSpinX = number of spin along the x direction
// yMomentum = momentum along the y direction
// nbrSpinY = number of spin along the y direction
// jFactor = Heisenberg XX coupling constant between nearest neighbors
// jzFactor = Heisenberg Z coupling constant between nearest neighbors
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::TwoDimensionalHeisenbergAnd2DTranslationHamiltonian(AbstractSpinChain* chain, int xMomentum, int nbrSpinX, 
													 int yMomentum, int nbrSpinY, double jFactor, double jzFactor,
													 AbstractArchitecture* architecture, long memory)
{
  this->Chain = chain;
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->JFactor = jFactor;
  this->JzFactor = jzFactor;
  this->HalfJFactor = 0.5 * this->JFactor;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  //  this->HermitianSymmetryFlag = false;
  this->HermitianSymmetryFlag = true;
  this->SzSzContributions = new double [MaxIndex - MinIndex + 1];
  this->EvaluateExponentialFactors();
  this->EvaluateDiagonalMatrixElements();
  if (memory == 0l)
    {
      this->FastMultiplicationFlag = false;
    }
  else
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      if (TmpMemory < 1024l)
	cout  << "fast = " <<  TmpMemory << "b ";
      else
	if (TmpMemory < (1l << 20))
	  cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	else
	  if (TmpMemory < (1l << 30))
	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	  else
	    cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
      cout << endl;
      this->EnableFastMultiplication();
    }
}

// destructor
//

TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::~TwoDimensionalHeisenbergAnd2DTranslationHamiltonian() 
{
  if (this->FastMultiplicationFlag == true)
    {
      long MinIndex;
      long MaxIndex;
      this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
      int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
      int ReducedDim = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      if ((ReducedDim * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
	++ReducedDim;
      for (int i = 0; i < ReducedDim; ++i)
	{
	  if (this->NbrInteractionPerComponent[i] > 0)
	    {
	      delete[] this->InteractionPerComponentIndex[i];
	      delete[] this->InteractionPerComponentCoefficient[i];
	    }
	}
      delete[] this->InteractionPerComponentCoefficient;
    }
}

// ask if Hamiltonian implements hermitian symmetry operations
//

bool TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::IsHermitian()
{
  return this->HermitianSymmetryFlag;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::ShiftHamiltonian (double shift)
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int TmpNbrComponents = ((int) (MaxIndex - MinIndex)) + 1;
  for (int i = 0; i < TmpNbrComponents; i ++)
    this->SzSzContributions[i] += shift;
}
  
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
											int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  int NbrTranslationsX;
  int NbrTranslationsY;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination[i] += this->SzSzContributions[i - this->PrecalculationShift] * vSource[i];
      Complex& TmpValue = vSource[i];
      for (int j = 0; j < this->NbrSpinX; j++)
 	{
 	  for (int k = 0; k < this->NbrSpinY; k++)
 	    {
	      int Index0 = this->GetLinearizedIndex(j, k);
	      int IndexX0 = this->GetSafeLinearizedIndex(j + 1, k);
	      int IndexY0 = this->GetSafeLinearizedIndex(j, k + 1);
 	      pos = this->Chain->SmiSpj(Index0, IndexX0, i, coef, NbrTranslationsX, NbrTranslationsY);
 	      if (pos != dim)
 		{
 		  vDestination[pos] += (coef * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValue;
 		}
 	      pos = this->Chain->SmiSpj(IndexX0, Index0, i, coef, NbrTranslationsX, NbrTranslationsY);
 	      if (pos != dim)
 		{
 		  vDestination[pos] += (coef * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValue;
 		}
 	      pos = this->Chain->SmiSpj(Index0, IndexY0, i, coef, NbrTranslationsX, NbrTranslationsY);
 	      if (pos != dim)
 		{
 		  vDestination[pos] += (coef * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValue;
 		}
 	      pos = this->Chain->SmiSpj(IndexY0, Index0, i, coef, NbrTranslationsX, NbrTranslationsY);
 	      if (pos != dim)
 		{
 		  vDestination[pos] += (coef * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValue;
 		}
	    }
	}
    }
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
												int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  Complex* TmpValues = new Complex[nbrVectors];
  int NbrTranslationsX;
  int NbrTranslationsY;
  for (int k = 0; k < nbrVectors; ++k)
    {
      ComplexVector& TmpSource = vSources[k];
      ComplexVector& TmpDestination = vDestinations[k];
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  TmpDestination[i] += this->SzSzContributions[i - this->PrecalculationShift] * TmpSource[i];
	}
    }
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int k = 0; k < nbrVectors; ++k)
	{
	  TmpValues[k] = vSources[k][i];
	}
//       for (int j = 0; j < this->NbrSpinX; j++)
// 	{
// 	  for (int k = 0; k < this->NbrSpinY; k++)
// 	    {
// 	      pos = this->Chain->Spi(this->GetLinearizedIndex(j, k), i, coef, NbrTranslationsX, NbrTranslationsY);
// 	      if (pos != dim)
// 		{
// 		  for (int l = 0; l < nbrVectors; ++l)
// 		    vDestinations[l][pos] += (coef * this->HxFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValues[l];
// 		}
// 	      pos = this->Chain->Smi(this->GetLinearizedIndex(j, k), i, coef, NbrTranslationsX, NbrTranslationsY);
// 	      if (pos != dim)
// 		{
// 		  for (int l = 0; l < nbrVectors; ++l)
// 		    vDestinations[l][pos] += (coef * this->HxFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValues[l];
// 		}
// 	    }
// 	}
    }
  delete[] TmpValues;
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

ComplexVector& TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
												 int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      AbstractSpinChain* TmpParticles = (AbstractSpinChain*) this->Chain->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->HermitianEvaluateAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  Complex Coefficient;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      Coefficient = vSource[k];
	      Complex TmpSum = 0.0;
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  TmpSum += Conj(TmpCoefficientArray[j]) * vSource[TmpIndexArray[j]];
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
		}
	      TmpSum += this->SzSzContributions[i - this->PrecalculationShift] * Coefficient;
	      vDestination[k++] += TmpSum;
	    }
	}
      else
	{
	  //	  this->HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	}
    }
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

ComplexVector* TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
													 int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      Complex* Coefficient2 = new Complex [nbrVectors];
      AbstractSpinChain* TmpParticles = (AbstractSpinChain*) this->Chain->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->HermitianEvaluateAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
      delete[] Coefficient2;
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  Complex Coefficient;
	  Complex* Coefficient2 = new Complex [nbrVectors];
	  Complex* TmpCoefficientArray; 
	  int j;
	  int Pos;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  Complex* TmpSum = new Complex [nbrVectors];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  TmpSum[l] = 0.0;
		  Coefficient2[l] = vSources[l][k];
		}
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		      TmpSum[l] += Conj(Coefficient) * vSources[l][Pos];
		    }
		}
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  TmpSum[l] += this->SzSzContributions[i - this->PrecalculationShift] * Coefficient2[l];
		  vDestinations[l][k] += TmpSum[l];
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	  delete[] TmpSum;
	}
      else
	{
	  //	  this->HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	}
    }
  return vDestinations;
}

// evaluate diagonal matrix elements
// 

void TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::EvaluateDiagonalMatrixElements()
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int dim = this->Chain->GetHilbertSpaceDimension();
  int TmpNbrComponents = ((int) (MaxIndex - MinIndex)) + 1;
  // SzSz part
  for (int i = 0; i < TmpNbrComponents; i++)
    {
      // SzSz part
      this->SzSzContributions[i] = 0.0;
      double Tmp = 0.0;
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      Tmp += this->Chain->SziSzj(this->GetLinearizedIndex(j, k), 
					 this->GetSafeLinearizedIndex(j + 1, k), i + this->PrecalculationShift);
	      Tmp += this->Chain->SziSzj(this->GetLinearizedIndex(j, k), 
					 this->GetSafeLinearizedIndex(j, k + 1), i + this->PrecalculationShift);
	    }
	}
      Tmp *= this->JzFactor; 
      this->SzSzContributions[i] += Tmp;
   }
}

// evaluate all exponential factors
//   

void TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::EvaluateExponentialFactors()
{
  this->ExponentialFactors = new Complex*[this->NbrSpinX];
  for (int i = 0; i < this->NbrSpinX; ++i)
    { 
      this->ExponentialFactors[i] = new Complex[this->NbrSpinY];
      for (int j = 0; j < this->NbrSpinY; ++j)
	{ 
	  this->ExponentialFactors[i][j] = Phase(2.0 * M_PI * ((this->XMomentum * ((double) i) / ((double) this->NbrSpinX))
							       + (this->YMomentum * ((double) j) / ((double) this->NbrSpinY))));
	}
    }
}

// core part of the AddMultiply method 
// 
// chain = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

void TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::HermitianEvaluateAddMultiplyComponent(AbstractSpinChain* chain, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  int NbrTranslationsX;
  int NbrTranslationsY;
  double Coefficient;
  int Index;
  Complex& TmpValue = vSource[index];
  Complex TmpSum = this->SzSzContributions[index - this->PrecalculationShift] * TmpValue;
  for (int j = 0; j < this->NbrSpinX; j++)
    {
      for (int k = 0; k < this->NbrSpinY; k++)
	{
	  int Index0 = this->GetLinearizedIndex(j, k);
	  int IndexX0 = this->GetSafeLinearizedIndex(j + 1, k);
	  int IndexY0 = this->GetSafeLinearizedIndex(j, k + 1);
	  Index = this->Chain->SmiSpj(Index0, IndexX0, index, Coefficient, NbrTranslationsX, NbrTranslationsY);
	  if (Index <= index)
	    {
	      vDestination[Index] += (Coefficient * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValue;
	      if (Index < index)
		{
		  TmpSum += (Coefficient * this->HalfJFactor) * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * vSource[Index];
		}
	    }
	  Index = this->Chain->SmiSpj(IndexX0, Index0, index, Coefficient, NbrTranslationsX, NbrTranslationsY);
	  if (Index <= index)
	    {
	      vDestination[Index] += (Coefficient * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValue;
	      if (Index < index)
		{
		  TmpSum += (Coefficient * this->HalfJFactor) * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * vSource[Index];
		}
	    }
	  Index = this->Chain->SmiSpj(Index0, IndexY0, index, Coefficient, NbrTranslationsX, NbrTranslationsY);
	  if (Index <= index)
	    {
	      vDestination[Index] += (Coefficient * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValue;
	      if (Index < index)
		{
		  TmpSum += (Coefficient * this->HalfJFactor) * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * vSource[Index];
		}
	    }
	  Index = this->Chain->SmiSpj(IndexY0, Index0, index, Coefficient, NbrTranslationsX, NbrTranslationsY);
	  if (Index <= index)
	    {
	      vDestination[Index] += (Coefficient * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValue;
	      if (Index < index)
		{
		  TmpSum += (Coefficient * this->HalfJFactor) * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]) * vSource[Index];
		}
	    }
	}
    }
  vDestination[index] += TmpSum;
}


// core part of the AddMultiply method for a set of vectors
// 
// chain = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

void TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::HermitianEvaluateAddMultiplyComponent(AbstractSpinChain* chain, int index, ComplexVector* vSources, 
												ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  Complex* TmpSum = new Complex[nbrVectors];
  for (int l = 0; l < nbrVectors; ++l)
    {
      tmpCoefficients[index] = vSources[l][index];
      TmpSum[l] = this->SzSzContributions[index - this->PrecalculationShift] * tmpCoefficients[index];
    }

  for (int l = 0; l < nbrVectors; ++l)
    vDestinations[l][index] += TmpSum[l];
  delete[] TmpSum;
}

// core part of the FastMultiplication method
// 
// chain = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray  

void TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::EvaluateFastMultiplicationComponent(AbstractSpinChain* chain, int index, 
											      int* indexArray, Complex* coefficientArray, long& position)
{
  int Index;
  int Dim = chain->GetHilbertSpaceDimension();
  double TmpCoefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;
  int AbsoluteIndex = index + this->PrecalculationShift;

  if (this->HermitianSymmetryFlag == false)
    {
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      int Index0 = this->GetLinearizedIndex(j, k);
	      int IndexX0 = this->GetSafeLinearizedIndex(j + 1, k);
	      int IndexY0 = this->GetSafeLinearizedIndex(j, k + 1);
	      Index = this->Chain->SmiSpj(Index0, IndexX0, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (Index != Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = (TmpCoefficient * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  ++position;
		}
	      Index = this->Chain->SmiSpj(IndexX0, Index0, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (Index != Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = (TmpCoefficient * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  ++position;
		}
	      Index = this->Chain->SmiSpj(Index0, IndexY0, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (Index != Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = (TmpCoefficient * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  ++position;
		}
	      Index = this->Chain->SmiSpj(IndexY0, Index0, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (Index != Dim)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = (TmpCoefficient * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  ++position;
		}
	    }
	}
    }
  else
    {
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      int Index0 = this->GetLinearizedIndex(j, k);
	      int IndexX0 = this->GetSafeLinearizedIndex(j + 1, k);
	      int IndexY0 = this->GetSafeLinearizedIndex(j, k + 1);
	      Index = this->Chain->SmiSpj(Index0, IndexX0, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (Index <= AbsoluteIndex)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = (TmpCoefficient * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  if (Index == AbsoluteIndex)
		    {
		      coefficientArray[position] *= 0.5;
		    }
		  ++position;
		}
	      Index = this->Chain->SmiSpj(IndexX0, Index0, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (Index <= AbsoluteIndex)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = (TmpCoefficient * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  if (Index == AbsoluteIndex)
		    {
		      coefficientArray[position] *= 0.5;
		    }
		  ++position;
		}
	      Index = this->Chain->SmiSpj(Index0, IndexY0, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (Index <= AbsoluteIndex)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = (TmpCoefficient * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  if (Index == AbsoluteIndex)
		    {
		      coefficientArray[position] *= 0.5;
		    }
		  ++position;
		}
	      Index = this->Chain->SmiSpj(IndexY0, Index0, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (Index <= AbsoluteIndex)
		{
		  indexArray[position] = Index;
		  coefficientArray[position] = (TmpCoefficient * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  if (Index == AbsoluteIndex)
		    {
		      coefficientArray[position] *= 0.5;
		    }
		  ++position;
		}
	    }
	}
    }
}

// core part of the PartialFastMultiplicationMemory
// 
// chain = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations  

void TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::EvaluateFastMultiplicationMemoryComponent(AbstractSpinChain* chain, int firstComponent, int lastComponent, long& memory)
{
  int Dim = chain->GetHilbertSpaceDimension();
  int Index;
  double TmpCoefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;

  if (this->HermitianSymmetryFlag == false)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrSpinX; j++)
	    {
	      for (int k = 0; k < this->NbrSpinY; k++)
		{
		  int Index0 = this->GetLinearizedIndex(j, k);
		  int IndexX0 = this->GetSafeLinearizedIndex(j + 1, k);
		  int IndexY0 = this->GetSafeLinearizedIndex(j, k + 1);
		  Index = this->Chain->SmiSpj(Index0, IndexX0, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index != Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		    }
		  Index = this->Chain->SmiSpj(IndexX0, Index0, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index != Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		    }
		  Index = this->Chain->SmiSpj(Index0, IndexY0, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index != Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		    }
		  Index = this->Chain->SmiSpj(IndexY0, Index0, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index != Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		    }
		}
	    }
	}
    }
  else
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrSpinX; j++)
	    {
	      for (int k = 0; k < this->NbrSpinY; k++)
		{
		  int Index0 = this->GetLinearizedIndex(j, k);
		  int IndexX0 = this->GetSafeLinearizedIndex(j + 1, k);
		  int IndexY0 = this->GetSafeLinearizedIndex(j, k + 1);
		  Index = this->Chain->SmiSpj(Index0, IndexX0, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		    }
		  Index = this->Chain->SmiSpj(IndexX0, Index0, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		    }
		  Index = this->Chain->SmiSpj(Index0, IndexY0, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (Index <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		    }
		  Index = this->Chain->SmiSpj(IndexY0, Index0, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
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

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::FastMultiplicationMemory(long allowedMemory)
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  this->NbrInteractionPerComponent = new int [EffectiveHilbertSpaceDimension];
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    this->NbrInteractionPerComponent[i] = 0;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;

  GenericHamiltonianPrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);
  //  this->PartialFastMultiplicationMemory(0, this->Chain->GetHilbertSpaceDimension());

 
  if (this->Architecture->GetOptimizedTypicalRange(this->NbrInteractionPerComponent, MinIndex, MaxIndex) == true)
    {
      this->PrecalculationShift = (int) MinIndex;
      EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
      cout << "distributed calculations have been reoptimized" << endl;
    }  
  if (allowedMemory == 0l)
    {
      delete[] this->NbrInteractionPerComponent;
      this->NbrInteractionPerComponent = 0;
      return 0l;
    }
  long Memory = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    {
      Memory += this->NbrInteractionPerComponent[i];
    }
  cout << "nbr interaction = " << Memory << endl;
  long TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(Complex*)) * EffectiveHilbertSpaceDimension;
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(Complex)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(Complex)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if (this->Chain->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(Complex*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	    Memory += this->NbrInteractionPerComponent[i];
	}
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(Complex*)) * ReducedSpaceDimension) + (Memory * (sizeof (int) + sizeof(Complex)));
      long ResidualMemory = allowedMemory - Memory;
      if (ResidualMemory > 0)
	{
	  int TotalReducedSpaceDimension = ReducedSpaceDimension;
	  int* TmpNbrInteractionPerComponent = new int [TotalReducedSpaceDimension];
	  int i = 0;
	  int Pos = 0;
	  for (; i < ReducedSpaceDimension; ++i)
	    {
	      TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[Pos];
	      Pos += this->FastMultiplicationStep;
	    }
	  delete[] this->NbrInteractionPerComponent;
	  this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
	}
      else
	{
	  int* TmpNbrInteractionPerComponent = new int [ReducedSpaceDimension];
	  for (int i = 0; i < ReducedSpaceDimension; ++i)
	    TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[i * this->FastMultiplicationStep];
	  delete[] this->NbrInteractionPerComponent;
	  this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
	}
    }
  else
    {
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(Complex*)) * EffectiveHilbertSpaceDimension) + (Memory * (sizeof (int) + sizeof(Complex)));
      this->FastMultiplicationStep = 1;
    }

  cout << "reduction factor=" << this->FastMultiplicationStep << endl;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  return Memory;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  long Memory = 0l;
  AbstractSpinChain* TmpSpinChain = (AbstractSpinChain*) this->Chain->Clone();
  int LastComponent = lastComponent + firstComponent;
  this->EvaluateFastMultiplicationMemoryComponent(TmpSpinChain, firstComponent, LastComponent, Memory);
  delete TmpSpinChain;
  return Memory;
}

// enable fast multiplication algorithm
//

void TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  int* TmpIndexArray;
  Complex* TmpCoefficientArray;
  long Pos;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficient = new Complex* [ReducedSpaceDimension];

  for (int i = 0; i < ReducedSpaceDimension; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [this->NbrInteractionPerComponent[i]];
      this->InteractionPerComponentCoefficient[i] = new Complex [this->NbrInteractionPerComponent[i]];
    }

  GenericHamiltonianPrecalculationOperation Operation(this, false);
  Operation.ApplyOperation(this->Architecture);
  //  this->PartialEnableFastMultiplication(0, this->Chain->GetHilbertSpaceDimension());

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
// nbrComponent  = index of the last component that has to be precalcualted

void TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{  
  int LastComponent = nbrComponent + firstComponent;
  AbstractSpinChain* TmpSpinChain = (AbstractSpinChain*) this->Chain->Clone();

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
      long TotalPos = 0;
      this->EvaluateFastMultiplicationComponent(TmpSpinChain, i, this->InteractionPerComponentIndex[Pos], 
						this->InteractionPerComponentCoefficient[Pos], TotalPos);
      ++Pos;
    }
  delete TmpSpinChain;
}

