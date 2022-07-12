////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of spin chain hamiltonian for the                   //
//                            4 qbits stabilizer code                         //
//                                                                            //
//                        last modification : 28/01/2018                      //
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


#include "Hamiltonian/SpinChain4QbitsStabilizerHamiltonian.h"
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


// constructor from default data
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

SpinChain4QbitsStabilizerHamiltonian::SpinChain4QbitsStabilizerHamiltonian(AbstractSpinChain* chain, int nbrSpin, bool periodicBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  this->HamiltonianShift = 0.0;
}

// destructor
//

SpinChain4QbitsStabilizerHamiltonian::~SpinChain4QbitsStabilizerHamiltonian() 
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

RealVector& SpinChain4QbitsStabilizerHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								      int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double TmpCoefficient1;
  double TmpCoefficient2;
  double TmpCoefficient3;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 3;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      double TmpValue = -2.0 * vSource[i];
      vDestination[i] += this->HamiltonianShift;
      for (int j = 0; j < MaxPos; ++j)
	{
	  TmpCoefficient1 = this->Chain->SziSzj(j, j + 3, i);
 	  pos = this->Chain->Spi(j + 1, i, TmpCoefficient2);
 	  if (pos != dim)
 	    {
	      pos2 = this->Chain->Spi(j + 2, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
	      pos2 = this->Chain->Smi(j + 2, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
 	    }
 	  pos = this->Chain->Smi(j + 1, i, TmpCoefficient2);
 	  if (pos != dim)
 	    {
	      pos2 = this->Chain->Spi(j + 2, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
	      pos2 = this->Chain->Smi(j + 2, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
 	    }
	}
      if (this->PeriodicBoundaryConditions == true)
	{
	  TmpCoefficient1 = this->Chain->SziSzj(this->NbrSpin - 3, 0, i);
 	  pos = this->Chain->Spi(this->NbrSpin - 2, i, TmpCoefficient2);
 	  if (pos != dim)
 	    {
	      pos2 = this->Chain->Spi(this->NbrSpin - 1, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
	      pos2 = this->Chain->Smi(this->NbrSpin - 1, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
 	    }
 	  pos = this->Chain->Smi(this->NbrSpin - 2, i, TmpCoefficient2);
 	  if (pos != dim)
 	    {
	      pos2 = this->Chain->Spi(this->NbrSpin - 1, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
	      pos2 = this->Chain->Smi(this->NbrSpin - 1, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
 	    }
	  TmpCoefficient1 = this->Chain->SziSzj(this->NbrSpin - 2, 1, i);
 	  pos = this->Chain->Spi(this->NbrSpin - 1, i, TmpCoefficient2);
 	  if (pos != dim)
 	    {
	      pos2 = this->Chain->Spi(0, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
	      pos2 = this->Chain->Smi(0, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
 	    }
 	  pos = this->Chain->Smi(this->NbrSpin - 1, i, TmpCoefficient2);
 	  if (pos != dim)
 	    {
	      pos2 = this->Chain->Spi(0, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
	      pos2 = this->Chain->Smi(0, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
 	    }
	  TmpCoefficient1 = this->Chain->SziSzj(this->NbrSpin - 1, 2, i);
 	  pos = this->Chain->Spi(0, i, TmpCoefficient2);
 	  if (pos != dim)
 	    {
	      pos2 = this->Chain->Spi(1, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
	      pos2 = this->Chain->Smi(1, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
 	    }
 	  pos = this->Chain->Smi(0, i, TmpCoefficient2);
 	  if (pos != dim)
 	    {
	      pos2 = this->Chain->Spi(1, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
	      pos2 = this->Chain->Smi(1, pos, TmpCoefficient3);
	      if (pos2 != dim)
		{
		  vDestination[pos2] += TmpValue * TmpCoefficient1 * TmpCoefficient2 * TmpCoefficient3;
		}
 	    }
	}
    }
  return vDestination;
}

