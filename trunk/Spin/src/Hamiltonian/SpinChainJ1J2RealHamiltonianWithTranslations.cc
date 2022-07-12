////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of spin chain J1-J2hamiltonian with translations           //
//                         at inversion symmetric points                      //
//                                                                            //
//                        last modification : 13/07/2016                      //
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


#include "Hamiltonian/SpinChainJ1J2RealHamiltonianWithTranslations.h"
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

SpinChainJ1J2RealHamiltonianWithTranslations::SpinChainJ1J2RealHamiltonianWithTranslations()
{
}

// constructor from default data
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// j1 = coupling constants between nearest neighbor spins
// j2 = coupling constants between second nearest neighbor spins

SpinChainJ1J2RealHamiltonianWithTranslations::SpinChainJ1J2RealHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin, double j1, double j2)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->J1 = j1;
  this->HalfJ1 = this->J1 * 0.5;
  this->J1z = this->J1;
  this->J2 = j2;
  this->HalfJ2 = this->J2 * 0.5;
  this->J2z = this->J2;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  this->EvaluateCosinusTable();
}

// destructor
//

SpinChainJ1J2RealHamiltonianWithTranslations::~SpinChainJ1J2RealHamiltonianWithTranslations() 
{
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& SpinChainJ1J2RealHamiltonianWithTranslations::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									     int firstComponent, int nbrComponent)
{
  double Coef;
  int NbrTranslation;
  int pos;
  int MaxPos = this->NbrSpin - 2;
  int Last = firstComponent + nbrComponent;
  for (int i = firstComponent; i < Last; i++)
    {
      vDestination[i] += this->SzSzContributions[i] * vSource[i];
      for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ1;
	      vDestination[pos] += Coef * (vSource[i] * this->ExponentialTable[NbrTranslation]);
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ1;
	      vDestination[pos] += Coef * (vSource[i] * this->ExponentialTable[NbrTranslation]);
	    }
	  pos = this->Chain->SmiSpj(j, j + 2, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ2;
	      vDestination[pos] += Coef * (vSource[i] * this->ExponentialTable[NbrTranslation]);
	    }
	  pos = this->Chain->SmiSpj(j + 2, j, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ2;
	      vDestination[pos] += Coef * (vSource[i] * this->ExponentialTable[NbrTranslation]);
	    }
	}    
      pos = this->Chain->SmiSpj(MaxPos, MaxPos + 1, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ1;
	  vDestination[pos] += Coef * (vSource[i] * this->ExponentialTable[NbrTranslation]);
	}
      pos = this->Chain->SmiSpj(MaxPos + 1, MaxPos, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ1;
	  vDestination[pos] += Coef * (vSource[i] * this->ExponentialTable[NbrTranslation]);
	}      
      pos = this->Chain->SmiSpj(MaxPos + 1, 0, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ1;
	  vDestination[pos] += Coef * (vSource[i] * this->ExponentialTable[NbrTranslation]);
	}
      pos = this->Chain->SmiSpj(0, MaxPos + 1, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ1;
	  vDestination[pos] += Coef * (vSource[i] * this->ExponentialTable[NbrTranslation]);
	}      
      pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ2;
	  vDestination[pos] += Coef * (vSource[i] * this->ExponentialTable[NbrTranslation]);
	}
      pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ2;
	  vDestination[pos] += Coef * (vSource[i] * this->ExponentialTable[NbrTranslation]);
	}      
      pos = this->Chain->SmiSpj(MaxPos + 1, 1, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ2;
	  vDestination[pos] += Coef * (vSource[i] * this->ExponentialTable[NbrTranslation]);
	}
      pos = this->Chain->SmiSpj(1, MaxPos + 1, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ2;
	  vDestination[pos] += Coef * (vSource[i] * this->ExponentialTable[NbrTranslation]);
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

RealVector* SpinChainJ1J2RealHamiltonianWithTranslations::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										  int firstComponent, int nbrComponent)
{
  double Coef;
  double Coef3;
  int NbrTranslation;
  int NbrTranslation2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 2;
  int Last = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double* TmpValues = new double[nbrVectors];
  double TmpCoef;
  for (int i = firstComponent; i < Last; i++)
    {
      for (int k = 0; k < nbrVectors; ++k)
	{
	  TmpValues[k] = vSources[k][i];
	  vDestinations[k][i] += this->SzSzContributions[i] * TmpValues[k];
	}
      for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      TmpCoef = this->HalfJ1 * Coef * this->ExponentialTable[NbrTranslation];
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][pos] += TmpCoef * TmpValues[k];
		}
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
	  if (pos != dim)
	    { 
	      TmpCoef = this->HalfJ1 * Coef * this->ExponentialTable[NbrTranslation];
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][pos] += TmpCoef * TmpValues[k];
		}
	    }
	  pos = this->Chain->SmiSpj(j, j + 2, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      TmpCoef = this->HalfJ2 * Coef * this->ExponentialTable[NbrTranslation];
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][pos] += TmpCoef * TmpValues[k];
		}
	    }
	  pos = this->Chain->SmiSpj(j + 2, j, i, Coef, NbrTranslation);
	  if (pos != dim)
	    { 
	      TmpCoef = this->HalfJ2 * Coef * this->ExponentialTable[NbrTranslation];
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][pos] += TmpCoef * TmpValues[k];
		}
	    }
	}    
      pos = this->Chain->SmiSpj(MaxPos, MaxPos + 1, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = this->HalfJ1 * Coef * this->ExponentialTable[NbrTranslation];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];
	    }
	}
      pos = this->Chain->SmiSpj(MaxPos + 1, MaxPos, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = this->HalfJ1 * Coef * this->ExponentialTable[NbrTranslation];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];
	    }
	}
      pos = this->Chain->SmiSpj(MaxPos + 1, 0, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = this->HalfJ1 * Coef * this->ExponentialTable[NbrTranslation];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];
	    }
	}
      pos = this->Chain->SmiSpj(0, MaxPos + 1, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = this->HalfJ1 * Coef * this->ExponentialTable[NbrTranslation];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];
	    }
	}

      pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = this->HalfJ2 * Coef * this->ExponentialTable[NbrTranslation];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];
	    }
	}
      pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = this->HalfJ2 * Coef * this->ExponentialTable[NbrTranslation];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];
	    }
	}
      pos = this->Chain->SmiSpj(MaxPos + 1, 1, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = this->HalfJ2 * Coef * this->ExponentialTable[NbrTranslation];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];
	    }
	}
      pos = this->Chain->SmiSpj(1, MaxPos + 1, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = this->HalfJ2 * Coef * this->ExponentialTable[NbrTranslation];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];
	    }
	}
    }
  delete[] TmpValues;
  return vDestinations;
}
 
// evaluate diagonal matrix elements
// 

void SpinChainJ1J2RealHamiltonianWithTranslations::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  int ReducedNbrSites = this->NbrSpin - 2;
  for (int i = 0; i < dim; i++)
    {
      this->SzSzContributions[i] = 0.0;
      for (int j = 0; j < ReducedNbrSites; j++)
	{
	  this->SzSzContributions[i] += this->Chain->SziSzj(j, j + 1, i) * this->J1z;
	  this->SzSzContributions[i] += this->Chain->SziSzj(j, j + 2, i) * this->J2z;
	}
      this->SzSzContributions[i] += this->Chain->SziSzj(this->NbrSpin - 2, this->NbrSpin - 1, i) * this->J1z;
      this->SzSzContributions[i] += this->Chain->SziSzj(this->NbrSpin - 1, 0, i) * this->J1z;
      this->SzSzContributions[i] += this->Chain->SziSzj(this->NbrSpin - 2, 0, i) * this->J2z;
      this->SzSzContributions[i] += this->Chain->SziSzj(this->NbrSpin - 1, 1, i) * this->J2z;
    }
}

