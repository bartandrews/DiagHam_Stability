////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of spin chain hamiltonian for the AKLT model            //
//                         written as a stabilizer code                       //
//                                                                            //
//                        last modification : 20/12/2017                      //
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


#include "Hamiltonian/SpinChainAKLTStabilizerHamiltonian.h"
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

SpinChainAKLTStabilizerHamiltonian::SpinChainAKLTStabilizerHamiltonian()
{
  this->Chain = 0;
  this->NbrSpin = 0;
  this->PeriodicBoundaryConditions = false;
  this->HamiltonianShift = 0.0;  
}

// constructor from default data
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

SpinChainAKLTStabilizerHamiltonian::SpinChainAKLTStabilizerHamiltonian(AbstractSpinChain* chain, int nbrSpin, bool periodicBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  this->HamiltonianShift = 0.0;
}

// destructor
//

SpinChainAKLTStabilizerHamiltonian::~SpinChainAKLTStabilizerHamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinChainAKLTStabilizerHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Chain = (AbstractSpinChain*) hilbertSpace;
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinChainAKLTStabilizerHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

SpinChainAKLTStabilizerHamiltonian& SpinChainAKLTStabilizerHamiltonian::SetChain(AbstractSpinChain* chain)
{  
  this->Chain = chain;
  return *this;  
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int SpinChainAKLTStabilizerHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void SpinChainAKLTStabilizerHamiltonian::ShiftHamiltonian (double shift)
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

RealVector& SpinChainAKLTStabilizerHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								    int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double TmpCoefficient1;
  double TmpCoefficient2;
  int pos;
  int MaxPos = this->NbrSpin - 2;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      double TmpValue = -2.0 * vSource[i];
      vDestination[i] += this->HamiltonianShift;
      for (int j = 0; j < MaxPos; ++j)
	{
	  TmpCoefficient1 = this->Chain->SziSzj(j, j + 2, i);
 	  pos = this->Chain->Spi(j + 1, i, TmpCoefficient2);
 	  if (pos != dim)
 	    {
 	      vDestination[pos] += TmpValue * TmpCoefficient1 * TmpCoefficient2;
 	    }
 	  pos = this->Chain->Smi(j + 1, i, TmpCoefficient2);
 	  if (pos != dim)
 	    {
 	      vDestination[pos] += TmpValue * TmpCoefficient1 * TmpCoefficient2;
 	    }
	}
      if (this->PeriodicBoundaryConditions == true)
	{
	  TmpCoefficient1 = this->Chain->SziSzj(this->NbrSpin - 2, 0, i);
 	  pos = this->Chain->Spi(this->NbrSpin - 1, i, TmpCoefficient2);
 	  if (pos != dim)
 	    {
 	      vDestination[pos] += TmpValue * TmpCoefficient1 * TmpCoefficient2;
 	    }
 	  pos = this->Chain->Smi(this->NbrSpin - 1, i, TmpCoefficient2);
 	  if (pos != dim)
 	    {
 	      vDestination[pos] += TmpValue * TmpCoefficient1 * TmpCoefficient2;
 	    }
	  TmpCoefficient1 = this->Chain->SziSzj(this->NbrSpin - 1, 1, i);
 	  pos = this->Chain->Spi(0, i, TmpCoefficient2);
 	  if (pos != dim)
 	    {
 	      vDestination[pos] += TmpValue * TmpCoefficient1 * TmpCoefficient2;
 	    }
 	  pos = this->Chain->Smi(0, i, TmpCoefficient2);
 	  if (pos != dim)
 	    {
 	      vDestination[pos] += TmpValue * TmpCoefficient1 * TmpCoefficient2;
 	    }
	}
    }
  return vDestination;
}

