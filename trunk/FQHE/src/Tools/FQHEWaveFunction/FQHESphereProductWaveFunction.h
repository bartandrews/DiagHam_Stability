////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of wave function obtained from product of            //
//                            wave function on sphere                         //
//                                                                            //
//                        last modification : 28/08/2009                      //
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


#ifndef FQHESPHEREPRODUCTWAVEFUNCTION_H
#define FQHESPHEREPRODUCTWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"
#include "GeneralTools/GarbageFlag.h"


class FQHESphereProductWaveFunction: public Abstract1DComplexFunctionOnSphere
{

 protected:

  // number of particles
  int NbrParticles;

  // additional power of the Jastrow factor in front of the wave function product
  int JastrowFactor;  

  // pointer to the wave functions
  Abstract1DComplexFunctionOnSphere** BaseWaveFunctions;
  // number of wavefunctions to multiply
  int NbrWaveFunctions;

  // temporary vector to store permutation of the input uv vector
  ComplexVector TemporaryUV;

  // if false use the same coordinates to evaluate wavefunctions, if true use coordinates 0 to nbrParticles/nbrWaveFunctions - 1 for the first wavefunctions and so on
  bool SeparateCoordinateFlag;
  // temporary vector to store permutation of the input uv vector in separate coodinate mode
  ComplexVector TemporaryUVSeparateCoordinates;

 public:

  // default constructor
  //
  FQHESphereProductWaveFunction();

  // constructor
  //
  // nbrParticles = number of particles
  // waveFunctions = array to the wavefunctions that have to be multiplied
  // nbrWaveFunctions = number of wavefunctions to multiply
  // jastrowFactor = multiply(if positive) or divide (if negative) the wavefunction by an overall Jastrow factor
  // separateCoordinates = if false use the same coordinates to evaluate wavefunctions, if true use coordinates 0 to nbrParticles/nbrWaveFunctions - 1 for the first wavefunctions and so on
  FQHESphereProductWaveFunction(int nbrParticles, Abstract1DComplexFunctionOnSphere** waveFunctions, int nbrWaveFunctions, int jastrowFactor = 0, bool separateCoordinates = false);

  // copy constructor
  //
  // function = reference on the wave function to copy
  FQHESphereProductWaveFunction(const FQHESphereProductWaveFunction& function);

  // destructor
  //
  ~FQHESphereProductWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  virtual Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
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
  virtual Complex CalculateFromSpinorVariables(ComplexVector& uv);

};

#endif
