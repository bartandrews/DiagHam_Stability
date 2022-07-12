////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                           Class author : Cecile Repellin                   //
//                                                                            //
//                 class of function basis for particle on CP2                //
//                                                                            //
//                        last modification : 07/02/2013                      //
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


#ifndef PARTICLEONCP2FUNCTIONBASIS_H
#define PARTICLEONCP2FUNCTIONBASIS_H


#include "config.h"
#include "FunctionBasis/AbstractFunctionBasis.h"

#define M_SQRT_2  1.414213562373095

class ParticleOnCP2FunctionBasis : public AbstractFunctionBasis
{
  protected :
    // twice the maximum Lz value reached by a particle
    int LzMax;
    // array containing numerical prefactor of each function
    double* Prefactor;
    // Chirality of the phase
    int Chirality;
    //number of flux quanta
    int NbrFluxQuanta;
    // array that gives the value of r for one particle corresponding to the linearized index
    int* quantumNumberR;
    // array that gives the value of s for one particle corresponding to the linearized index
    int* quantumNumberS;
  
  public:
   
  // constructor
  //
  // nbrFluxQuanta = number of flux quanta
  // chirality = flag that allows to choose between either one of two conventions for
  // the phase of the orbitals
  ParticleOnCP2FunctionBasis(int nbrFluxQuanta, int chirality=1);

  
  // destructor
  //
  ~ParticleOnCP2FunctionBasis ();

  // get value of the i-th function at a given point (for functions which take values in C)
  //
  // value = reference on the value where the function has to be evaluated
  // result = reference on the value where the result has to be stored
  // index = function index 
  void GetFunctionValue(RealVector& value, Complex& result, int index);

  // get value of the i-th function at a given point (theta, phi=0), which happens to be real
  //
  // theta = reference on the value where the function has to be evaluated
  // index = function index 

  double GetRealFunctionValue(double theta1, double theta2, int index);

  // calculate them all:
  void GetRealFunctionValues(double theta1, double theta2, double *result);
  
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


