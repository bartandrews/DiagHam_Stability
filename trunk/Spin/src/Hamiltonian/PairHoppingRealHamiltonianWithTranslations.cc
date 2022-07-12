////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class of pair-hopping hanmiltonian written in the spin-1 language      //
//                        with translations at k=0 or pi                      //
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


#include "Hamiltonian/PairHoppingRealHamiltonianWithTranslations.h"
#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslations.h"
#include "Vector/RealVector.h"
#include "Vector/RealVector.h"
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

PairHoppingRealHamiltonianWithTranslations::PairHoppingRealHamiltonianWithTranslations()
{
}

// constructor from default data
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// pValue = value that defines the filling factor p/(2p+1)

PairHoppingRealHamiltonianWithTranslations::PairHoppingRealHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin, int pValue)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->PValue = pValue;
  this->GlobalEnergyShift = 0.0;
  this->J = 1.0;
  this->HalfJ = this->J * 0.5;
  this->Jz = this->J;
  this->SzSzContributions = new double [1];
  this->EvaluateCosinusTable();
}

// destructor
//

PairHoppingRealHamiltonianWithTranslations::~PairHoppingRealHamiltonianWithTranslations() 
{
}

void PairHoppingRealHamiltonianWithTranslations::ShiftHamiltonian (double shift)
{
  this->GlobalEnergyShift = shift;
}

// evaluate all cosinus/sinus that are needed when computing matrix elements
//

void PairHoppingRealHamiltonianWithTranslations::EvaluateCosinusTable()
{
  int TmpNbrUnitCells = this->NbrSpin / this->PValue;
  this->CosinusTable = new double [2 * TmpNbrUnitCells];
  this->SinusTable = new double [2 * TmpNbrUnitCells];
  this->ExponentialTable = new double [2 * TmpNbrUnitCells];
  double Coef = 2.0 * M_PI * ((double) this->Chain->GetMomentum()) / ((double) TmpNbrUnitCells);
  for (int i = 0; i < (2 * TmpNbrUnitCells); ++i)
    {
      this->CosinusTable[i] = cos(Coef * ((double) i));
      this->SinusTable[i] = 0.0;
      this->ExponentialTable[i] = this->CosinusTable[i];
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

RealVector& PairHoppingRealHamiltonianWithTranslations::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									    int firstComponent, int nbrComponent)
{
  double Coef;
  double Coef2;
  double Coef3;
  int NbrTranslation;
  int NbrTranslation2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  int Last = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  int NbrUnitCells = this->NbrSpin / this->PValue;
  PairHoppingP1AsSpin1ChainWithTranslations* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslations*) this->Chain->Clone();
  for (int i = firstComponent; i < Last; i++)
    {
      double TmpValue = vSource[i];
      vDestination[i] += this-> GlobalEnergyShift * TmpValue;
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

RealVector* PairHoppingRealHamiltonianWithTranslations::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									       int firstComponent, int nbrComponent)
{
  double Coef;
  double Coef3;
  int NbrTranslation;
  int NbrTranslation2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  int Last = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double* TmpValues = new double[nbrVectors];
  double TmpCoef;
  int NbrUnitCells = this->NbrSpin / this->PValue;
  PairHoppingP1AsSpin1ChainWithTranslations* TmpSpace = (PairHoppingP1AsSpin1ChainWithTranslations*) this->Chain->Clone();
  for (int i = firstComponent; i < Last; i++)
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
  delete[] TmpValues;
  return vDestinations;
}
 
