////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of Potts 3 chain hamiltonian with translations             //
//                     for the dual O'Brien-Fendley model                     //
//                                                                            //
//                        last modification : 02/12/2019                      //
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


#include "Hamiltonian/Potts3ChainDualOBrienFendleyRealHamiltonianWithTranslations.h"
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


// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin

Potts3ChainDualOBrienFendleyRealHamiltonianWithTranslations::Potts3ChainDualOBrienFendleyRealHamiltonianWithTranslations(Potts3ChainWithTranslations* chain, int nbrSpin)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->ReducedNbrSpin = this->NbrSpin -1;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  this->EvaluateCosinusTable();
}

// destructor
//

Potts3ChainDualOBrienFendleyRealHamiltonianWithTranslations::~Potts3ChainDualOBrienFendleyRealHamiltonianWithTranslations() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void Potts3ChainDualOBrienFendleyRealHamiltonianWithTranslations::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (Potts3ChainWithTranslations*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* Potts3ChainDualOBrienFendleyRealHamiltonianWithTranslations::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int Potts3ChainDualOBrienFendleyRealHamiltonianWithTranslations::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& Potts3ChainDualOBrienFendleyRealHamiltonianWithTranslations::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
											     int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double Coefficient;
  double CoefficientTranslations;
  double SzkCoefficient;
  int pos;
  int NbrTranslations;
  Potts3ChainWithTranslations* TmpChain = (Potts3ChainWithTranslations*) (this->Chain);
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      double TmpValue = vSource[i];
      vDestination[i] += this->SzSzContributions[i] * TmpValue;
      for (int j = 0; j < this->ReducedNbrSpin; ++j)
	{
	  pos = TmpChain->SzkSmiSpj(j, j + 1, j, i, CoefficientTranslations, SzkCoefficient, NbrTranslations);
	  if (pos != dim)
	    {
	      TmpChain->Szi(j + 1, i, Coefficient);
	      if (Coefficient == 1.0)
		{
		  vDestination[pos] += (2.0 * CoefficientTranslations) * this->ExponentialTable[NbrTranslations] * TmpValue;
		}
	      else
		{
		  vDestination[pos] -= CoefficientTranslations * this->ExponentialTable[NbrTranslations] * TmpValue;
		}
	      if (SzkCoefficient == 1.0)
		{
		  vDestination[pos] += (2.0 * CoefficientTranslations) * this->ExponentialTable[NbrTranslations] * TmpValue;
		}
	      else
		{
		  vDestination[pos] -= CoefficientTranslations * this->ExponentialTable[NbrTranslations] * TmpValue;
		}
	    }
	  pos = TmpChain->SzkSmiSpj(j + 1, j, j, i, CoefficientTranslations, SzkCoefficient, NbrTranslations);
	  if (pos != dim)
	    {
	      TmpChain->Szi(j + 1, i, Coefficient);
	      if (Coefficient == -1.0)
		{
		  vDestination[pos] += (2.0 * CoefficientTranslations) * this->ExponentialTable[NbrTranslations] * TmpValue;
		}
	      else
		{
		  vDestination[pos] -= CoefficientTranslations * this->ExponentialTable[NbrTranslations] * TmpValue;
		}
	      if (SzkCoefficient == -1.0)
		{
		  vDestination[pos] += (2.0 * CoefficientTranslations) * this->ExponentialTable[NbrTranslations] * TmpValue;
		}
	      else
		{
		  vDestination[pos] -= CoefficientTranslations * this->ExponentialTable[NbrTranslations] * TmpValue;
		}
	    }
	}
      pos = TmpChain->SzkSmiSpj(this->ReducedNbrSpin, 0, this->ReducedNbrSpin, i, CoefficientTranslations, SzkCoefficient, NbrTranslations);
      if (pos != dim)
	{
	  TmpChain->Szi(0, i, Coefficient);
	  if (Coefficient == 1.0)
	    {
	      vDestination[pos] += (2.0 * CoefficientTranslations) * this->ExponentialTable[NbrTranslations] * TmpValue;
	    }
	  else
	    {
	      vDestination[pos] -= CoefficientTranslations * this->ExponentialTable[NbrTranslations] * TmpValue;
	    }
	  if (SzkCoefficient == 1.0)
	    {
	      vDestination[pos] += (2.0 * CoefficientTranslations) * this->ExponentialTable[NbrTranslations] * TmpValue;
	    }
	  else
	    {
	      vDestination[pos] -= CoefficientTranslations * this->ExponentialTable[NbrTranslations] * TmpValue;
	    }
	}
      pos = TmpChain->SzkSmiSpj(0, this->ReducedNbrSpin, this->ReducedNbrSpin, i, CoefficientTranslations, SzkCoefficient, NbrTranslations);
      if (pos != dim)
	{
	  TmpChain->Szi(0, i, Coefficient);
	  if (Coefficient == -1.0)
	    {
	      vDestination[pos] += (2.0 * CoefficientTranslations) * this->ExponentialTable[NbrTranslations] * TmpValue;
	    }
	  else
	    {
	      vDestination[pos] -= CoefficientTranslations * this->ExponentialTable[NbrTranslations] * TmpValue;
	    }
	  if (SzkCoefficient == -1.0)
	    {
	      vDestination[pos] += (2.0 * CoefficientTranslations) * this->ExponentialTable[NbrTranslations] * TmpValue;
	    }
	  else
	    {
	      vDestination[pos] -= CoefficientTranslations * this->ExponentialTable[NbrTranslations] * TmpValue;
	    }
	}
    }
  return vDestination;
}

// evaluate diagonal matrix elements
// 

void Potts3ChainDualOBrienFendleyRealHamiltonianWithTranslations::EvaluateDiagonalMatrixElements()
{
  int Dimension = this->Chain->GetHilbertSpaceDimension();
  double Coefficient;
  for (int i = 0; i < Dimension; i++)
    {
      double Tmp = 0.0;
      for (int j = 0; j < this->NbrSpin; j++)
	{
	  this->Chain->Szi(j, i, Coefficient);
	  Tmp += -2.0 + 3.0 * (Coefficient * Coefficient);
	}
      this->SzSzContributions[i] = Tmp;
    }
}

