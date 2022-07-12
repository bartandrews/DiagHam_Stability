////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class of pair-hopping hanmiltonian written in the spin-1 language      //
//                               with translations                            //
//                                                                            //
//                        last modification : 13/03/2019                      //
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


#include "Hamiltonian/PairHoppingHamiltonianWithTranslations.h"
#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslations.h"
#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslationsLong.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"

#include <iostream>


using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

PairHoppingHamiltonianWithTranslations::PairHoppingHamiltonianWithTranslations()
{
}

// constructor from default data
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// pValue = value that defines the filling factor p/(2p+1)
// electrostaticPerturbation1 = first electrostatic perturbation 
// electrostaticPerturbation2 = second electrostatic perturbation 
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication

PairHoppingHamiltonianWithTranslations::PairHoppingHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin, int pValue,
									       double electrostaticPerturbation1, double electrostaticPerturbation2,
									       AbstractArchitecture* architecture, long memory)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->PValue = pValue;
  this->J = 1.0;
  this->HalfJ = this->J * 0.5;
  this->Jz = this->J;
  this->ElectrostaticPerturbation1 = electrostaticPerturbation1;
  this->ElectrostaticPerturbation2 = electrostaticPerturbation2;
  this->EvaluateCosinusTable();
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->SzSzContributions = new double [MaxIndex - MinIndex + 1];
  this->PrecalculationShift = (int) MinIndex;  
  this->HermitianSymmetryFlag = true;
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

PairHoppingHamiltonianWithTranslations::~PairHoppingHamiltonianWithTranslations() 
{
}

// evaluate all cosinus/sinus that are needed when computing matrix elements
//

void PairHoppingHamiltonianWithTranslations::EvaluateCosinusTable()
{
  int TmpNbrUnitCells = this->NbrSpin / this->PValue;
  this->CosinusTable = new double [2 * TmpNbrUnitCells];
  this->SinusTable = new double [2 * TmpNbrUnitCells];
  this->ExponentialTable = new Complex [2 * TmpNbrUnitCells];
  double Coef = 2.0 * M_PI * ((double) this->Chain->GetMomentum()) / ((double) (TmpNbrUnitCells));
  for (int i = 0; i < (2 * (TmpNbrUnitCells)); ++i)
    {
      this->CosinusTable[i] = cos(Coef * ((double) i));
      this->SinusTable[i] = sin(Coef * ((double) i));
      this->ExponentialTable[i] = Phase(Coef * ((double) i));
    }
}

// evaluate diagonal matrix elements
// 

void PairHoppingHamiltonianWithTranslations::EvaluateDiagonalMatrixElements()
{
  int Dim = this->Chain->GetHilbertSpaceDimension();

  if ((this->ElectrostaticPerturbation1 == 0.0) && (this->ElectrostaticPerturbation2 == 0.0))
    {
      for (int i = 0; i < Dim; i++)
	{
	  this->SzSzContributions[i] = 0.0;
	}
    }
  else
    {
      PairHoppingP1AsSpin1ChainWithTranslations* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslations*) this->Chain->Clone();
      int NbrUnitCells = this->NbrSpin / this->PValue;
      for (int i = 0; i < Dim; i++)
	{
	  double Tmp = 0.0;
	  for (int j = 1; j < NbrUnitCells; ++j)
	    {
	      Tmp += this->ElectrostaticPerturbation1 * TmpSpace->ZPlus(j * this->PValue - 1, i) * TmpSpace->ZMinus(j * this->PValue, i);
	    }	
	  Tmp += this->ElectrostaticPerturbation1 * TmpSpace->ZPlus(this->NbrSpin - 1, i) * TmpSpace->ZMinus(0, i);
	  for (int j = 0; j < NbrUnitCells; ++j)
	    {
	      for (int p = 1; p < this->PValue; ++p)
		{
		  Tmp += this->ElectrostaticPerturbation1 * TmpSpace->ZPlus(j * this->PValue + p - 1, i) * TmpSpace->Z0(j * this->PValue + p, i);
		  Tmp += this->ElectrostaticPerturbation1 * TmpSpace->Z0(j * this->PValue + p - 1, i) * TmpSpace->ZMinus(j * this->PValue + p, i);
		  Tmp += this->ElectrostaticPerturbation2 * TmpSpace->ZPlus(j * this->PValue + p - 1, i) * TmpSpace->ZPlus(j * this->PValue + p, i);
		  Tmp += this->ElectrostaticPerturbation2 * TmpSpace->Z0(j * this->PValue + p - 1, i) * TmpSpace->Z0(j * this->PValue + p, i);
		  Tmp += this->ElectrostaticPerturbation2 * TmpSpace->ZMinus(j * this->PValue + p - 1, i) * TmpSpace->ZMinus(j * this->PValue + p, i);
		}
	    }
	  this->SzSzContributions[i] = Tmp;
	}
      delete TmpSpace;
    }
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void PairHoppingHamiltonianWithTranslations::ShiftHamiltonian (double shift)
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  for (long i = MinIndex; i <= MaxIndex; ++i)
    {
      this->SzSzContributions[i - MinIndex] += shift;
    }
}
  
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& PairHoppingHamiltonianWithTranslations::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									     int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      int dim = this->Chain->GetHilbertSpaceDimension();
      double Coef;
      Complex Coef2;
      double Coef3;
      int NbrTranslation;
      int NbrTranslation2;
      int pos;
      int pos2;
      int MaxPos = this->NbrSpin - 1;
      int NbrUnitCells = this->NbrSpin / this->PValue;
      if (this->NbrSpin <= 32)
	{
	  PairHoppingP1AsSpin1ChainWithTranslations* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslations*) this->Chain->Clone();
	  for (int i = firstComponent; i < LastComponent; i++)
	    {
	      Complex TmpValue = vSource[i];
	      vDestination[i] += this->SzSzContributions[i] * TmpValue;
	      for (int j = 1; j < NbrUnitCells; ++j)
		{
		  pos = TmpSpace->PlusMinusOperator(j - 1, j, i, Coef, NbrTranslation);
		  if (pos != dim)
		    {
		      vDestination[pos] += Coef * this->ExponentialTable[NbrTranslation]  * TmpValue;
		    }
		}
	      pos = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, i, Coef, NbrTranslation);
	      if (pos != dim)
		{
		  vDestination[pos] += Coef * this->ExponentialTable[NbrTranslation]  * TmpValue;
		}
	      for (int j = 0; j < NbrUnitCells; ++j)
		{
		  for (int p = 1; p < this->PValue; ++p)
		    {
		      pos = TmpSpace->SwapOperator(j, p - 1, i, Coef, NbrTranslation);
		      if (pos != dim)
			{
			  vDestination[pos] += Coef * this->ExponentialTable[NbrTranslation]  * TmpValue;
			}
		    }
		}
	    }
	  delete TmpSpace;
	}
      else
	{
	  PairHoppingP1AsSpin1ChainWithTranslationsLong* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslationsLong*) this->Chain->Clone();
	  for (int i = firstComponent; i < LastComponent; i++)
	    {
	      Complex TmpValue = vSource[i];
	      vDestination[i] += this->SzSzContributions[i] * TmpValue;
	      for (int j = 1; j < NbrUnitCells; ++j)
		{
		  pos = TmpSpace->PlusMinusOperator(j - 1, j, i, Coef, NbrTranslation);
		  if (pos != dim)
		    {
		      vDestination[pos] += Coef * this->ExponentialTable[NbrTranslation]  * TmpValue;
		    }
		}
	      pos = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, i, Coef, NbrTranslation);
	      if (pos != dim)
		{
		  vDestination[pos] += Coef * this->ExponentialTable[NbrTranslation]  * TmpValue;
		}
	      for (int j = 0; j < NbrUnitCells; ++j)
		{
		  for (int p = 1; p < this->PValue; ++p)
		    {
		      pos = TmpSpace->SwapOperator(j, p - 1, i, Coef, NbrTranslation);
		      if (pos != dim)
			{
			  vDestination[pos] += Coef * this->ExponentialTable[NbrTranslation]  * TmpValue;
			}
		    }
		}
	    }
	  delete TmpSpace;
	}
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
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
		}
	      vDestination[k++] += this->SzSzContributions[k] * Coefficient;
	    }
	}
    }
  return vDestination;
}
 
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& PairHoppingHamiltonianWithTranslations::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										    int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      int dim = this->Chain->GetHilbertSpaceDimension();
      double Coef;
      Complex Coef2;
      double Coef3;
      int NbrTranslation;
      int NbrTranslation2;
      int pos;
      int pos2;
      int MaxPos = this->NbrSpin - 1;
      int NbrUnitCells = this->NbrSpin / this->PValue;
      if (this->NbrSpin <= 32)
	{
	  PairHoppingP1AsSpin1ChainWithTranslations* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslations*) this->Chain->Clone();
	  for (int i = firstComponent; i < LastComponent; i++)
	    {
	      Complex TmpValue = vSource[i];
	      Complex TmpSum = 0.0;
	      vDestination[i] += this->SzSzContributions[i] * TmpValue;
	      for (int j = 1; j < NbrUnitCells; ++j)
		{
		  pos = TmpSpace->PlusMinusOperator(j - 1, j, i, Coef, NbrTranslation);
		  if (pos <= i)
		    {
		      vDestination[pos] += Coef * this->ExponentialTable[NbrTranslation]  * TmpValue;
		      if (pos < i)
			{
			  TmpSum += Coef * Conj(this->ExponentialTable[NbrTranslation]) * vSource[pos];
			}
		    }
		}
	      pos = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, i, Coef, NbrTranslation);
	      if (pos <= i)
		{
		  vDestination[pos] += Coef * this->ExponentialTable[NbrTranslation]  * TmpValue;
		  if (pos < i)
		    {
		      TmpSum += Coef * Conj(this->ExponentialTable[NbrTranslation]) * vSource[pos];
		    }
		}
	      for (int j = 0; j < NbrUnitCells; ++j)
		{
		  for (int p = 1; p < this->PValue; ++p)
		    {
		      pos = TmpSpace->SwapOperator(j, p - 1, i, Coef, NbrTranslation);
		      if (pos <= i)
			{
			  vDestination[pos] += Coef * this->ExponentialTable[NbrTranslation]  * TmpValue;
			  if (pos < i)
			    {
			      TmpSum += Coef * Conj(this->ExponentialTable[NbrTranslation])  * vSource[pos];
			    }
			}
		    }
		}
	      vDestination[i] += TmpSum;
	    }
	  delete TmpSpace;
	}
      else
	{
	  PairHoppingP1AsSpin1ChainWithTranslationsLong* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslationsLong*) this->Chain->Clone();
	  for (int i = firstComponent; i < LastComponent; i++)
	    {
	      Complex TmpValue = vSource[i];
	      Complex TmpSum = 0.0;
	      vDestination[i] += this->SzSzContributions[i] * TmpValue;
	      for (int j = 1; j < NbrUnitCells; ++j)
		{
		  pos = TmpSpace->PlusMinusOperator(j - 1, j, i, Coef, NbrTranslation);
		  if (pos <= i)
		    {
		      vDestination[pos] += Coef * this->ExponentialTable[NbrTranslation]  * TmpValue;
		      if (pos < i)
			{
			  TmpSum += Coef * Conj(this->ExponentialTable[NbrTranslation])  *  vSource[pos];
			}
		    }
		}
	      pos = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, i, Coef, NbrTranslation);
	      if (pos <= i)
		{
		  vDestination[pos] += Coef * this->ExponentialTable[NbrTranslation]  * TmpValue;
		  if (pos < i)
		    {
		      TmpSum += Coef * Conj(this->ExponentialTable[NbrTranslation])  * vSource[pos];
		    }
		}
	      for (int j = 0; j < NbrUnitCells; ++j)
		{
		  for (int p = 1; p < this->PValue; ++p)
		    {
		      pos = TmpSpace->SwapOperator(j, p - 1, i, Coef, NbrTranslation);
		      if (pos <= i)
			{
			  vDestination[pos] += Coef * this->ExponentialTable[NbrTranslation]  * TmpValue;
			  if (pos < i)
			    {
			      TmpSum += Coef * Conj(this->ExponentialTable[NbrTranslation])  * vSource[pos];
			    }
			}
		    }
		}
	      vDestination[i] += TmpSum;
	    }
	  delete TmpSpace;
	}
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
	      TmpSum += this->SzSzContributions[k] * Coefficient;
	      vDestination[k++] += TmpSum;
	    }
	}
    }
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

ComplexVector* PairHoppingHamiltonianWithTranslations::LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      double Coef;
      double Coef3;
      int NbrTranslation;
      int NbrTranslation2;
      int pos;
      int pos2;
      int MaxPos = this->NbrSpin - 1;
      int dim = this->Chain->GetHilbertSpaceDimension();
      Complex* TmpValues = new Complex[nbrVectors];
      Complex TmpCoef;
      int NbrUnitCells = this->NbrSpin / this->PValue;
      if (this->NbrSpin <= 32)
	{
	  PairHoppingP1AsSpin1ChainWithTranslations* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslations*) this->Chain->Clone();
	  for (int i = firstComponent; i < LastComponent; i++)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  TmpValues[k] = vSources[k][i];
		  vDestinations[k][i] += this->SzSzContributions[i] * TmpValues[k];
		}
	      for (int j = 1; j < NbrUnitCells; ++j)
		{
		  pos = TmpSpace->PlusMinusOperator(j - 1, j, i, Coef, NbrTranslation);
		  if (pos != dim)
		    {
		      TmpCoef = Coef * this->ExponentialTable[NbrTranslation];
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] += TmpCoef * TmpValues[k];
			}
		    }
		}
	      pos = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, i, Coef, NbrTranslation);
	      if (pos != dim)
		{
		  TmpCoef = Coef * this->ExponentialTable[NbrTranslation];
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      vDestinations[k][pos] += TmpCoef * TmpValues[k];
		    }
		}
	      for (int j = 0; j < NbrUnitCells; ++j)
		{
		  for (int p = 1; p < this->PValue; ++p)
		    {
		      pos = TmpSpace->SwapOperator(j, p - 1, i, Coef, NbrTranslation);
		      if (pos != dim)
			{
			  TmpCoef = Coef * this->ExponentialTable[NbrTranslation];
			  for (int k = 0; k < nbrVectors; ++k)
			    {
			      vDestinations[k][pos] += TmpCoef * TmpValues[k];
			    }
			}
		    }
		}
	    }      
	  delete TmpSpace;
	}
      else
	{
	  PairHoppingP1AsSpin1ChainWithTranslationsLong* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslationsLong*) this->Chain->Clone();
	  for (int i = firstComponent; i < LastComponent; i++)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  TmpValues[k] = vSources[k][i];
		  vDestinations[k][i] += this->SzSzContributions[i] * TmpValues[k];
		}
	      for (int j = 1; j < NbrUnitCells; ++j)
		{
		  pos = TmpSpace->PlusMinusOperator(j - 1, j, i, Coef, NbrTranslation);
		  if (pos != dim)
		    {
		      TmpCoef = Coef * this->ExponentialTable[NbrTranslation];
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] += TmpCoef * TmpValues[k];
		    }
		    }
		}
	      pos = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, i, Coef, NbrTranslation);
	      if (pos != dim)
		{
		  TmpCoef = Coef * this->ExponentialTable[NbrTranslation];
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      vDestinations[k][pos] += TmpCoef * TmpValues[k];
		    }
		}
	      for (int j = 0; j < NbrUnitCells; ++j)
		{
		  for (int p = 1; p < this->PValue; ++p)
		    {
		      pos = TmpSpace->SwapOperator(j, p - 1, i, Coef, NbrTranslation);
		      if (pos != dim)
			{
			  TmpCoef = Coef * this->ExponentialTable[NbrTranslation];
			  for (int k = 0; k < nbrVectors; ++k)
			    {
			      vDestinations[k][pos] += TmpCoef * TmpValues[k];
			    }
			}
		    }
		}
	    }      
	  delete TmpSpace;
	}    
      delete[] TmpValues;
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
		  TmpSum[l] += this->SzSzContributions[k] * Coefficient2[l];
		  vDestinations[l][k] += TmpSum[l];
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	  delete[] TmpSum;
	}
    }
  return vDestinations;
}
 
// core part of the FastMultiplication method
// 
// chain = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray  

void PairHoppingHamiltonianWithTranslations::EvaluateFastMultiplicationComponent(AbstractSpinChain* chain, int index, 
										 int* indexArray, Complex* coefficientArray, long& position)
{
  double Coef;
  Complex Coef2;
  double Coef3;
  int NbrTranslation;
  int NbrTranslation2;
  int TmpIndex;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  int Dim = chain->GetHilbertSpaceDimension();
  int NbrUnitCells = this->NbrSpin / this->PValue;
  int AbsoluteIndex = index + this->PrecalculationShift;
  if (this->HermitianSymmetryFlag == false)
    {
      if (this->NbrSpin <= 32)
	{
	  PairHoppingP1AsSpin1ChainWithTranslations* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslations*) chain;
	  for (int j = 1; j < NbrUnitCells; ++j)
	    {
	      TmpIndex = TmpSpace->PlusMinusOperator(j - 1, j, AbsoluteIndex, Coef, NbrTranslation);
	      if (TmpIndex != Dim)
		{
		  indexArray[position] = TmpIndex;
		  coefficientArray[position] = Coef * this->ExponentialTable[NbrTranslation];
		  ++position;
		}
	    }
	  TmpIndex = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, AbsoluteIndex, Coef, NbrTranslation);
	  if (TmpIndex != Dim)
	    {
	      indexArray[position] = TmpIndex;
	      coefficientArray[position] = Coef * this->ExponentialTable[NbrTranslation];
	      ++position;
	    }
	  for (int j = 0; j < NbrUnitCells; ++j)
	    {
	      for (int p = 1; p < this->PValue; ++p)
		{
		  TmpIndex = TmpSpace->SwapOperator(j, p - 1, AbsoluteIndex, Coef, NbrTranslation);
		  if (TmpIndex != Dim)
		    {
		      indexArray[position] = TmpIndex;
		      coefficientArray[position] = Coef * this->ExponentialTable[NbrTranslation];
		      ++position;
		    }
		}
	    }
	}
      else
	{
	  PairHoppingP1AsSpin1ChainWithTranslationsLong* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslationsLong*) this->Chain->Clone();
	  for (int j = 1; j < NbrUnitCells; ++j)
	    {
	      TmpIndex = TmpSpace->PlusMinusOperator(j - 1, j, AbsoluteIndex, Coef, NbrTranslation);
	      if (TmpIndex != Dim)
		{
		  indexArray[position] = TmpIndex;
		  coefficientArray[position] = Coef * this->ExponentialTable[NbrTranslation];
		  ++position;
		}
	    }
	  TmpIndex = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, AbsoluteIndex, Coef, NbrTranslation);
	  if (TmpIndex != Dim)
	    {
	      indexArray[position] = TmpIndex;
	      coefficientArray[position] = Coef * this->ExponentialTable[NbrTranslation];
	      ++position;
	    }
	  for (int j = 0; j < NbrUnitCells; ++j)
	    {
	      for (int p = 1; p < this->PValue; ++p)
		{
		  TmpIndex = TmpSpace->SwapOperator(j, p - 1, AbsoluteIndex, Coef, NbrTranslation);
		  if (TmpIndex != Dim)
		    {
		      indexArray[position] = TmpIndex;
		      coefficientArray[position] = Coef * this->ExponentialTable[NbrTranslation];
		      ++position;
		    }
		}
	    }
	}
    }
  else
     {
      if (this->NbrSpin <= 32)
	{
	  PairHoppingP1AsSpin1ChainWithTranslations* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslations*) chain;
	  for (int j = 1; j < NbrUnitCells; ++j)
	    {
	      TmpIndex = TmpSpace->PlusMinusOperator(j - 1, j, AbsoluteIndex, Coef, NbrTranslation);
	      if (TmpIndex <= AbsoluteIndex)
		{
		  indexArray[position] = TmpIndex;
		  coefficientArray[position] = Coef * this->ExponentialTable[NbrTranslation];
		  if (TmpIndex == AbsoluteIndex)
		    {
		      coefficientArray[position] *= 0.5;
		    }
		  ++position;
		}
	    }
	  TmpIndex = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, AbsoluteIndex, Coef, NbrTranslation);
	  if (TmpIndex <= AbsoluteIndex)
	    {
	      indexArray[position] = TmpIndex;
	      coefficientArray[position] = Coef * this->ExponentialTable[NbrTranslation];
	      if (TmpIndex == AbsoluteIndex)
		{
		  coefficientArray[position] *= 0.5;
		}
	      ++position;
	    }
	  for (int j = 0; j < NbrUnitCells; ++j)
	    {
	      for (int p = 1; p < this->PValue; ++p)
		{
		  TmpIndex = TmpSpace->SwapOperator(j, p - 1, AbsoluteIndex, Coef, NbrTranslation);
		  if (TmpIndex <= AbsoluteIndex)
		    {
		      indexArray[position] = TmpIndex;
		      coefficientArray[position] = Coef * this->ExponentialTable[NbrTranslation];
		      if (TmpIndex == AbsoluteIndex)
			{
			  coefficientArray[position] *= 0.5;
			}
		      ++position;
		    }
		}
	    }
	}
      else
	{
	  PairHoppingP1AsSpin1ChainWithTranslationsLong* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslationsLong*) this->Chain->Clone();
	  for (int j = 1; j < NbrUnitCells; ++j)
	    {
	      TmpIndex = TmpSpace->PlusMinusOperator(j - 1, j, AbsoluteIndex, Coef, NbrTranslation);
	      if (TmpIndex <= AbsoluteIndex)
		{
		  indexArray[position] = TmpIndex;
		  coefficientArray[position] = Coef * this->ExponentialTable[NbrTranslation];
		  if (TmpIndex == AbsoluteIndex)
		    {
		      coefficientArray[position] *= 0.5;
		    }
		  ++position;
		}
	    }
	  TmpIndex = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, AbsoluteIndex, Coef, NbrTranslation);
	  if (TmpIndex <= AbsoluteIndex)
	    {
	      indexArray[position] = TmpIndex;
	      coefficientArray[position] = Coef * this->ExponentialTable[NbrTranslation];
	      if (TmpIndex == AbsoluteIndex)
		{
		  coefficientArray[position] *= 0.5;
		}
	      ++position;
	    }
	  for (int j = 0; j < NbrUnitCells; ++j)
	    {
	      for (int p = 1; p < this->PValue; ++p)
		{
		  TmpIndex = TmpSpace->SwapOperator(j, p - 1, AbsoluteIndex, Coef, NbrTranslation);
		  if (TmpIndex <= AbsoluteIndex)
		    {
		      indexArray[position] = TmpIndex;
		      coefficientArray[position] = Coef * this->ExponentialTable[NbrTranslation];
		      if (TmpIndex == AbsoluteIndex)
			{
			  coefficientArray[position] *= 0.5;
			}
		      ++position;
		    }
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

void PairHoppingHamiltonianWithTranslations::EvaluateFastMultiplicationMemoryComponent(AbstractSpinChain* chain, int firstComponent, int lastComponent, long& memory)
{
  double Coef;
  Complex Coef2;
  double Coef3;
  int NbrTranslation;
  int NbrTranslation2;
  int TmpIndex;
  int TmpIndex2;
  int MaxPos = this->NbrSpin - 1;
  int Dim = chain->GetHilbertSpaceDimension();
  int NbrUnitCells = this->NbrSpin / this->PValue;
  if (this->HermitianSymmetryFlag == false)
    {
      if (this->NbrSpin <= 32)
	{
	  PairHoppingP1AsSpin1ChainWithTranslations* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslations*) chain;
	  for (int i = firstComponent; i < lastComponent; i++)
	    {
	      for (int j = 1; j < NbrUnitCells; ++j)
		{
		  TmpIndex = TmpSpace->PlusMinusOperator(j - 1, j, i, Coef, NbrTranslation);
		  if (TmpIndex < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];		     
		    }
		}
	      TmpIndex = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, i, Coef, NbrTranslation);
	      if (TmpIndex < Dim)
		{
		  ++memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];		     
		}
	      for (int j = 0; j < NbrUnitCells; ++j)
		{
		  for (int p = 1; p < this->PValue; ++p)
		    {
		      TmpIndex = TmpSpace->SwapOperator(j, p - 1, i, Coef, NbrTranslation);
		      if (TmpIndex < Dim)
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
	  PairHoppingP1AsSpin1ChainWithTranslationsLong* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslationsLong*) this->Chain->Clone();
	  for (int i = firstComponent; i < lastComponent; i++)
	    {
	      for (int j = 1; j < NbrUnitCells; ++j)
		{
		  TmpIndex = TmpSpace->PlusMinusOperator(j - 1, j, i, Coef, NbrTranslation);
		  if (TmpIndex < Dim)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];		     
		    }
		}
	      TmpIndex = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, i, Coef, NbrTranslation);
	      if (TmpIndex < Dim)
		{
		  ++memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];		     
		}
	      for (int j = 0; j < NbrUnitCells; ++j)
		{
		  for (int p = 1; p < this->PValue; ++p)
		    {
		      TmpIndex = TmpSpace->SwapOperator(j, p - 1, i, Coef, NbrTranslation);
		      if (TmpIndex < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];		     
			}
		    }
		}
	    }
	}
    }
  else
    {
     if (this->NbrSpin <= 32)
	{
	  PairHoppingP1AsSpin1ChainWithTranslations* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslations*) chain;
	  for (int i = firstComponent; i < lastComponent; i++)
	    {
	      for (int j = 1; j < NbrUnitCells; ++j)
		{
		  TmpIndex = TmpSpace->PlusMinusOperator(j - 1, j, i, Coef, NbrTranslation);
		  if (TmpIndex <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];		     
		    }
		}
	      TmpIndex = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, i, Coef, NbrTranslation);
	      if (TmpIndex <= i)
		{
		  ++memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];		     
		}
	      for (int j = 0; j < NbrUnitCells; ++j)
		{
		  for (int p = 1; p < this->PValue; ++p)
		    {
		      TmpIndex = TmpSpace->SwapOperator(j, p - 1, i, Coef, NbrTranslation);
		      if (TmpIndex <= i)
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
	  PairHoppingP1AsSpin1ChainWithTranslationsLong* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslationsLong*) this->Chain->Clone();
	  for (int i = firstComponent; i < lastComponent; i++)
	    {
	      for (int j = 1; j < NbrUnitCells; ++j)
		{
		  TmpIndex = TmpSpace->PlusMinusOperator(j - 1, j, i, Coef, NbrTranslation);
		  if (TmpIndex <= i)
		    {
		      ++memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];		     
		    }
		}
	      TmpIndex = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, i, Coef, NbrTranslation);
	      if (TmpIndex <= i)
		{
		  ++memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];		     
		}
	      for (int j = 0; j < NbrUnitCells; ++j)
		{
		  for (int p = 1; p < this->PValue; ++p)
		    {
		      TmpIndex = TmpSpace->SwapOperator(j, p - 1, i, Coef, NbrTranslation);
		      if (TmpIndex <= i)
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

