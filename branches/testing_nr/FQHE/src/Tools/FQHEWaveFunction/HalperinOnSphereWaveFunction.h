////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of Halperin wave function on sphere                 //
//                                                                            //
//                        last modification : 11/11/2006                      //
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


#ifndef HALPERINONSPHEREWAVEFUNCTION_H
#define HALPERINONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"


class HalperinOnSphereWaveFunction: public Abstract1DComplexFunctionOnSphere
{

 protected:

  // number of particles with spin up
  int NbrSpinUpParticles;
  // number of particles with spin down
  int NbrSpinDownParticles;
  // total number of particles
  int TotalNbrParticles;

  // m1 index ( (m1, m2, n) Halperin wave function)
  int M1Index;
  // m2 index ( (m1, m2, n) Halperin wave function)
  int M2Index;
  // n index ( (m1, m2, n) Halperin wave function)
  int NIndex;

  // temporary arrays used during wave function evaluation
  Complex* UCoordinates;
  Complex* VCoordinates;

 public:

  // constructor
  //
  // nbrSpinUpParticles = number of particles with spin up
  // nbrSpinDownParticles = number of particles with spin down
  // m1Index = m1 index ( (m1, m2, n) Halperin wave function)
  // m2Index = m2 index ( (m1, m2, n) Halperin wave function)
  // nIndex = n index ( (m1, m2, n) Halperin wave function)
  HalperinOnSphereWaveFunction(int nbrSpinUpParticles, int nbrSpinDownParticles, int m1Index, int m2Index, int nIndex);

  // copy constructor
  //
  // function = reference on the wave function to copy
  HalperinOnSphereWaveFunction(const HalperinOnSphereWaveFunction& function);

  // destructor
  //
   ~HalperinOnSphereWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point (the first 2*nbrSpinUpParticles coordinates correspond to the position of the spin up particles, 
  //                                     the other 2*nbrSpinDownParticles are spin down particle positions)
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

  // evaluate function at a given point
  //
  // uv = ensemble of spinor variables on sphere describing point
  //      where function has to be evaluated
  //      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
  // return value = function value at (uv)
  Complex CalculateFromSpinorVariables(ComplexVector& uv);

};

#endif
