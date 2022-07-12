////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    class of unprojected Jain composite fermion wave function on sphere     //
//                                                                            //
//                        last modification : 13/04/2005                      //
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


#ifndef UNPROJECTEDJAINCFONSPHEREWAVEFUNCTION_H
#define UNPROJECTEDJAINCFONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "Tools/FQHEWaveFunction/JainCFOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/LandauSpectrumOnSphere.h"


class ConfigurationParser;


class UnprojectedJainCFOnSphereWaveFunction: public JainCFOnSphereWaveFunction
{

 protected:


 public:

  // constructor
  //
  // filename = name of the file describing the occupation of the pseudo-Landau levels
  UnprojectedJainCFOnSphereWaveFunction(char* filename);

  // copy constructor
  //
  // function = reference on the wave function to copy
  UnprojectedJainCFOnSphereWaveFunction(const UnprojectedJainCFOnSphereWaveFunction& function);

  // destructor
  //
   ~UnprojectedJainCFOnSphereWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

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
  Complex CalculateFromSpinorVariables(ComplexVector& uv);

 protected:

  // evaluate monopole spherical harmonic 
  //
  // coordinate = index of the coordinate
  // momentum = monopole spherical harmonic Lz momentum (plus S shift)
  // landauLevel = index of the Landau level
  // maximumMomentum = maxixum momentum that can be reached in the Landau level
  // return value = value of the monopole spherical harmonic at the given point
  Complex EvaluateMonopoleHarmonic (int coordinate, int momentum, int landauLevel, int maximumMomentum);

  // evaluate constant factors that appears in the sum of projected monopole harmonic (except LLL)
  //
  void EvaluateSumPrefactors();

};

#endif
