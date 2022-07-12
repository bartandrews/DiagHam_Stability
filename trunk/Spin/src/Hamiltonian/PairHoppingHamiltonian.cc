////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class of pair-hopping hanmiltonian written in the spin-1 language      //
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


#include "Hamiltonian/PairHoppingHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "HilbertSpace/PairHoppingP1AsSpin1Chain.h"
#include "HilbertSpace/PairHoppingP1AsSpin1ChainLong.h"

#include <iostream>


using std::cout;
using std::endl;
using std::ostream;


// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// pValue = value that defines the filling factor p/(2p+1)
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

PairHoppingHamiltonian::PairHoppingHamiltonian(AbstractSpinChain* chain, int nbrSpin, int pValue, bool periodicBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->PValue = pValue;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  this->GlobalEnergyShift = 0.0;
}

// destructor
//

PairHoppingHamiltonian::~PairHoppingHamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void PairHoppingHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Chain = (PairHoppingP1AsSpin1Chain*) hilbertSpace;
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* PairHoppingHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

PairHoppingHamiltonian& PairHoppingHamiltonian::SetChain(AbstractSpinChain* chain)
{  
  this->Chain = (PairHoppingP1AsSpin1Chain*) chain;
  return *this;  
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int PairHoppingHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void PairHoppingHamiltonian::ShiftHamiltonian (double shift)
{
  this->GlobalEnergyShift = shift;
}
  
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& PairHoppingHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
							  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double Coef;
  double Coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  int NbrUnitCells = this->NbrSpin / this->PValue;
  if (this->NbrSpin <= 32)
    {
      PairHoppingP1AsSpin1Chain* TmpSpace = (PairHoppingP1AsSpin1Chain*) this->Chain->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  double TmpValue = vSource[i];
	  vDestination[i] += this->GlobalEnergyShift * TmpValue;
	  for (int j = 1; j < NbrUnitCells; ++j)
	    {
	      pos = TmpSpace->PlusMinusOperator(j - 1, j, i);
	      if (pos != dim)
		{
		  vDestination[pos] += TmpValue;
		}
	    }
	  if (this->PeriodicBoundaryConditions == true)
	    {
	      pos = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, i);
	      if (pos != dim)
		{
		  vDestination[pos] += TmpValue;
		}
	    }      
	  for (int j = 0; j < NbrUnitCells; ++j)
	    {
	      for (int p = 1; p < this->PValue; ++p)
		{
		  pos = TmpSpace->SwapOperator(j, p - 1, i);
		  if (pos != dim)
		    {
		      vDestination[pos] += TmpValue;
		    }
		}
	    }
	}
      delete TmpSpace;
   }
  else
    {
      PairHoppingP1AsSpin1ChainLong* TmpSpace = (PairHoppingP1AsSpin1ChainLong*) this->Chain->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  double TmpValue = vSource[i];
	  vDestination[i] += this->GlobalEnergyShift * TmpValue;
	  for (int j = 1; j < NbrUnitCells; ++j)
	    {
	      pos = TmpSpace->PlusMinusOperator(j - 1, j, i);
	      if (pos != dim)
		{
		  vDestination[pos] += TmpValue;
		}
	    }
	  if (this->PeriodicBoundaryConditions == true)
	    {
	      pos = TmpSpace->PlusMinusOperator(NbrUnitCells - 1, 0, i);
	      if (pos != dim)
		{
		  vDestination[pos] += TmpValue;
		}
	    }      
	  for (int j = 0; j < NbrUnitCells; ++j)
	    {
	      for (int p = 1; p < this->PValue; ++p)
		{
		  pos = TmpSpace->SwapOperator(j, p - 1, i);
		  if (pos != dim)
		    {
		      vDestination[pos] += TmpValue;
		    }
		}
	    }
	}
      delete TmpSpace;
   }
  return vDestination;
}

