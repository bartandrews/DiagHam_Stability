////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of Laughlin wave function on disk with one quasihole           //
//                         (without the gaussian factor)                      //
//                                                                            //
//                        last modification : 15/11/2008                      //
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


#ifndef FQHEDISKLAUGHLINONEQUASIHOLEWAVEFUNCTION_H
#define FQHEDISKLAUGHLINONEQUASIHOLEWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"


class FQHEDiskLaughlinOneQuasiholeWaveFunction: public Abstract1DComplexFunction
{

 protected:

  // position of the quasihole
  Complex ZHole;

  // number of particles
  int NbrParticles;

  // inverse value of the filling factor
  int InvFillingFactor;

 public:

  // constructor
  //
  // nbrParticles = number of particles
  // zHole = position of the quasihole
  // invFillingFactor = inverse value of the filling factor
  FQHEDiskLaughlinOneQuasiholeWaveFunction(int nbrParticles, Complex zHole, int invFillingFactor);

  // copy constructor
  //
  // function = reference on the wave function to copy
  FQHEDiskLaughlinOneQuasiholeWaveFunction(const FQHEDiskLaughlinOneQuasiholeWaveFunction& function);

  // destructor
  //
   ~FQHEDiskLaughlinOneQuasiholeWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

};

#endif
