////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of spin chain hamiltonian implementing the              //
//                             O'Brien-Fendley model                          //
//                                                                            //
//                        last modification : 06/01/2022                      //
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


#include "Hamiltonian/SpinChainOBrienFendleyHamiltonian.h"
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
// squareFactor = numerical factor in front of the 1/3 (S_i S_i+1)^2 term
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

SpinChainOBrienFendleyHamiltonian::SpinChainOBrienFendleyHamiltonian(AbstractSpinChain* chain, int nbrSpin, bool periodicBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->SquareFactor = 1.0;
  this->LinearFactor = -1.0;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

SpinChainOBrienFendleyHamiltonian::~SpinChainOBrienFendleyHamiltonian() 
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

RealVector& SpinChainOBrienFendleyHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double Coef;
  double Coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      double TmpValue = vSource[i];
      vDestination[i] += this->SzSzContributions[i] * TmpValue;
      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; ++j)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, Coef);
	  if (pos != dim)
	    {
	      Coef2 = 0.5 * Coef * TmpValue;
	      vDestination[pos] += this->LinearFactor * Coef2;
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, Coef);
	  if (pos != dim)
	    {
	      Coef2 = 0.5 * Coef * TmpValue;
	      vDestination[pos] += this->LinearFactor * Coef2;
	    }
	  pos = this->Chain->SmiSpjSmiSpj(j, j + 1, j, j + 1, i, Coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += (0.25 * this->SquareFactor * Coef) * TmpValue;	      
	    }	  
	  pos = this->Chain->SmiSpjSmiSpj(j + 1, j, j + 1, j, i, Coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += (0.25 * this->SquareFactor * Coef) * TmpValue;	      
	    }	  
	}
      if (this->PeriodicBoundaryConditions == true)
	{
	  pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef);
	  if (pos != dim)
	    {
	      Coef2 = 0.5 * Coef * TmpValue;
	      vDestination[pos] += this->LinearFactor * Coef2;
	    }
	  pos = this->Chain->SmiSpj(0, MaxPos, i, Coef);
	  if (pos != dim)
	    {
	      Coef2 = 0.5 * Coef * TmpValue;
	      vDestination[pos] += this->LinearFactor * Coef2;
	    }
	  pos = this->Chain->SmiSpjSmiSpj(MaxPos, 0, MaxPos, 0, i, Coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += (0.25 * this->SquareFactor * Coef) * TmpValue;	      
	    }	  
	  pos = this->Chain->SmiSpjSmiSpj(0, MaxPos, 0, MaxPos, i, Coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += (0.25 * this->SquareFactor * Coef) * TmpValue;	      
	    }	  
	}
    }
  return vDestination;
}

// evaluate diagonal matrix elements
// 

void SpinChainOBrienFendleyHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; ++i)
    {
      this->SzSzContributions[i] = 0.0;
    }
}

