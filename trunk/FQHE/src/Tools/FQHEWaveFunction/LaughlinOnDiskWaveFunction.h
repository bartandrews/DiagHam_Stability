////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of Laughlin wave function on disk                  //
//                         (without the gaussian factor)                      //
//                                                                            //
//                        last modification : 10/10/2004                      //
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


#ifndef LAUGHLINONDISKWAVEFUNCTION_H
#define LAUGHLINONDISKWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"


class LaughlinOnDiskWaveFunction: public Abstract1DComplexFunction
{

 protected:

  // number of particles
  int NbrParticles;

  // inverse value of the filling factor
  int InvFillingFactor;

  // invert of the maximum x value
  double InvScale;

  // log - scale (to be added in exponent)
  double LogScale;

  // last value of Sums of Square coordinates on invoking operator()
  double SumSqr;

  // flag for exponential factor
  bool ExponentialFactors;

 public:

  // constructor
  //
  // nbrParticles = number of particles
  // invFillingFactor = inverse value of the filling factor
  // scale = typical sytem size
  // useExponentials = flag whether to use exponential factors
  LaughlinOnDiskWaveFunction(int nbrParticles, int invFillingFactor, double scale = 1.0, bool useExponentials = false);

  // copy constructor
  //
  // function = reference on the wave function to copy
  LaughlinOnDiskWaveFunction(const LaughlinOnDiskWaveFunction& function);

  // destructor
  //
   ~LaughlinOnDiskWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // change the normalization of the funtion by a multiplicative factor
  // factor = factor to be multiplied
  virtual void Renormalize(double factor);

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

};

#endif
