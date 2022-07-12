////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                      Class author Cecile Repellin                          //
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


#ifndef FQHECP2GENERALIZEDLAUGHLINWAVEFUNCTION_H
#define FQHECP2GENERALIZEDLAUGHLINWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "GeneralTools/GarbageFlag.h"
#include "HilbertSpace/BosonOnCP2.h"
#include "Matrix/ComplexLapackDeterminant.h"


class FQHECP2GeneralizedLaughlinWaveFunction
{

 protected:

  // number of particles
  int NbrParticles;
  // number of flux quanta
  int NbrFluxQuanta;
  //One Body Hilbert Space
  BosonOnCP2* Space;
  // array that gives the value of r for one particle corresponding to the linearized index
  int* quantumNumberR;
  // array that gives the value of s for one particle corresponding to the linearized index
  int* quantumNumberS;
  // determinant whose square has to be computed
  ComplexLapackDeterminant* TmpDeterminant;
  // flag that indicates whether the wave function should be computed for NbrFluxQuanta odd (true) or even (false)
  bool OddFlag;
  // exponent of the generalized Laughlin wavefunction
  int LaughlinExponent;

 public:

  // default constructor
  //
  FQHECP2GeneralizedLaughlinWaveFunction();

  // constructor
  //
  // nbrParticles = number of particles
  // nbrFluxQuanta = number of flux quanta
  // exponent = exponent of the generalized Laughlin wavefunction
  FQHECP2GeneralizedLaughlinWaveFunction(ComplexLapackDeterminant* determinant, int nbrParticles, int nbrFluxQuanta, int exponent = 2);

   // destructor
  //
  ~FQHECP2GeneralizedLaughlinWaveFunction();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  virtual Complex CalculateFromCoordinates(RealVector& x);
  
  // get the quantum numbers tz, y of a one particle state 
  //
  //quantumNumberTz = array that gives the quantum number tz for a single particle stateDescription
  //quantumNumberY = array that gives the quantum number y for a single particle state
  inline void GetQuantumNumbersFromLinearizedIndex(int* quantumNumberR, int* quantumNumberS)
  {
    for (int tzMax = 0; tzMax <= this->NbrFluxQuanta; ++tzMax)
    {
      for (int shiftedTz = 0; shiftedTz <= tzMax; ++shiftedTz)
	{
	  int tz = 2*shiftedTz - tzMax;
	  int y = 3*tzMax - 2*this->NbrFluxQuanta;
	  int index = this->GetLinearizedIndex(tz, y, 1);
	  quantumNumberR[index] = (y + 3*tz + 2*this->NbrFluxQuanta)/6;
	  quantumNumberS[index] = (y - 3*tz + 2*this->NbrFluxQuanta)/6;
// 	  cout << tz << " " << y << " " << quantumNumberR[index] << " " << quantumNumberS[index] << endl;
	}
    }
  }

  //compute the linearized index for the single particles quantum numbers (tz,y)
  //
  //tz = integer with the value of quantum number tz
  //y = integer with the value of quantum number y
  //nbrParticles = number of particles involved (1 for HilbertSpace, 2 for two-body interaction...)
  inline int GetLinearizedIndex(int tz, int y, int nbrParticles)
  {
    int tzMax = (y + 2*nbrParticles*this->NbrFluxQuanta)/3;
    int index = tzMax*(tzMax + 1)/2 + (tz + tzMax)/2 ;
//     cout <<index << endl;
//     int index = this->LzMax - (tzMax*(tzMax + 1)/2 + (tz + tzMax)/2) ;
    return index;
  }

  
};

#endif
