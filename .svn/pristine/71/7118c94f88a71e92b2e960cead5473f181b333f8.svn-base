////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of exact state wave function                    //
//                                                                            //
//                        last modification : 20/04/2005                      //
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


#ifndef EXACTWAVEFUNCTION_H
#define EXACTWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "Vector/RealVector.h"

#ifdef USE_HILBERT_SPACE

class AbstractQHEParticle;
class AbstractFunctionBasis;


class ExactWaveFunction: public Abstract1DComplexFunction
{

 protected:

  // vector that describes exact state components in ExactStateSpace basis
  RealVector ExactState;

  // Hilbert space associated to the exact state
  AbstractQHEParticle* ExactStateSpace;

  // one body real space basis to use 
  AbstractFunctionBasis* OneBodyBasis;

  // index that describes the one component exact state (-1 if  ExactState has to be used instead of the one component exact state)
  int ExactStateIndex;

 public:

  // constructor
  //
  // components = vector that describes exact state components in ExactStateSpace basis
  // exactStateSpace = Hilbert space associated to the exact state
  // oneBodyBasis = one body real space basis to use 
  ExactWaveFunction(RealVector& components, AbstractQHEParticle* exactStateSpace, AbstractFunctionBasis* oneBodyBasis);

  // constructor for a one component exact state
  //
  // index = index that describes the one component exact state
  // exactStateSpace = Hilbert space associated to the exact state
  // oneBodyBasis = one body real space basis to use 
  ExactWaveFunction(int index, AbstractQHEParticle* exactStateSpace, AbstractFunctionBasis* oneBodyBasis);

  // copy constructor
  //
  // function = reference on the wave function to copy
  ExactWaveFunction(const ExactWaveFunction& function);

  // destructor
  //
   ~ExactWaveFunction();

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

#endif

