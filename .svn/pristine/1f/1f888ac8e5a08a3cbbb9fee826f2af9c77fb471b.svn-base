////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                          Class author Cecile Repellin                      //
//                                                                            //
//                 class of wave function for the ground state                //
//                            of bosons on CP2                                //
//                                                                            //
//                        last modification : 08/02/2013                      //
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
#include "Tools/FQHEWaveFunction/FQHECP2GeneralizedLaughlinWaveFunction.h"
#include "Vector/RealVector.h"
#include "HilbertSpace/BosonOnCP2.h"
#include "Matrix/ComplexLapackDeterminant.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
//

FQHECP2GeneralizedLaughlinWaveFunction::FQHECP2GeneralizedLaughlinWaveFunction()
{
}

// constructor
//
// nbrParticles = number of particles
// wavefunctions = array to the wavefunctions that have to be multiplied
// nbrWaveFunctions = number of wavefunctions to multiply
// exponent = exponent of the generalized Laughlin wavefunction
 
FQHECP2GeneralizedLaughlinWaveFunction::FQHECP2GeneralizedLaughlinWaveFunction(ComplexLapackDeterminant* determinant, int nbrParticles, int nbrFluxQuanta, int exponent)
{
  this->NbrParticles = nbrParticles;
  this->OddFlag = false;
  this->LaughlinExponent = exponent;
  if ((nbrFluxQuanta % this->LaughlinExponent) == 0)
    this->NbrFluxQuanta = nbrFluxQuanta / this->LaughlinExponent;
  if ((this->LaughlinExponent == 2) && (nbrFluxQuanta == 3))
  {
    this->NbrFluxQuanta = 1;
    this->OddFlag = true;
  }
  this->quantumNumberR = new int[this->NbrParticles];
  this->quantumNumberS = new int[this->NbrParticles];
  this->GetQuantumNumbersFromLinearizedIndex(this->quantumNumberR, this->quantumNumberS);
  this->TmpDeterminant = determinant;
}


// destructor
//

FQHECP2GeneralizedLaughlinWaveFunction::~FQHECP2GeneralizedLaughlinWaveFunction()
{
  delete[] this->quantumNumberR;
  delete[] this->quantumNumberS;
}


// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex FQHECP2GeneralizedLaughlinWaveFunction::CalculateFromCoordinates(RealVector& x)
{  
  Complex MatrixElement;
  double Arg;
  int r;
  int s;
  int t;
  Complex result;
  if (this->OddFlag == false)
  {
    for (int i = 0; i < this->NbrParticles; ++i)
      {
	for (int j = 0; j < this->NbrParticles; ++j)
	{
	  r = this->quantumNumberR[i];
	  s = this->quantumNumberS[i];
	  t = this->NbrFluxQuanta  - r - s;
	  Arg = x[4*j + 2] * (double) (s) + x[4*j + 3] * (double) (t);
	  MatrixElement.Re = cos(Arg);
	  MatrixElement.Im = sin(Arg);
	  MatrixElement *= pow( cos(x[4*j]) , (double) (r)) * pow( sin(x[4*j]) * cos(x[4*j + 1]), (double) (s)) * pow( sin(x[4*j]) * sin(x[4*j + 1]), (double) (t));
	  this->TmpDeterminant->SetMatrixElement(i, j, MatrixElement);
	}
      }
    result = this->TmpDeterminant->Determinant();
    Complex TmpResult = result;
    for (int i = 1; i < this->LaughlinExponent; ++i)
      TmpResult *= result;
    return TmpResult;
   }
  else
  {
    result = 1.0;
    for (int particle = 0; particle < 4; ++particle)
    {
      for (int i = 0; i < 4; ++i)
	{
	  for (int j = 0; j < 4; ++j)
	      {
		if ( (i != particle) && (j != particle) )
		{
		  int Tmpi = i;
		  int Tmpj = j;
		  if (i == 3)
		    Tmpi = particle;
		  if (j == 3)
		    Tmpj = particle;
		  r = this->quantumNumberR[Tmpi];
		  s = this->quantumNumberS[Tmpi];
		  t = this->NbrFluxQuanta  - r - s;
		  Arg = x[4*j + 2] * (double) (s) + x[4*j + 3] * (double) (t);
		  MatrixElement.Re = cos(Arg);
		  MatrixElement.Im = sin(Arg);
		  MatrixElement *= pow( cos(x[4*j]) , (double) (r)) * pow( sin(x[4*j]) * cos(x[4*j + 1]), (double) (s)) * pow( sin(x[4*j]) * sin(x[4*j + 1]), (double) (t));
		  this->TmpDeterminant->SetMatrixElement(Tmpi, Tmpj, MatrixElement);
		}
	      }
	    }
	result *= this->TmpDeterminant->Determinant();
	}
  return result;
  }
}



