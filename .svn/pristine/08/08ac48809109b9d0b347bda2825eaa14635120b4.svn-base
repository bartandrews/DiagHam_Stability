////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of Jain composite fermion wave function on sphere          //
//                                                                            //
//                        last modification : 10/01/2005                      //
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


#ifndef JAINCFONSPHEREWAVEFUNCTION_H
#define JAINCFONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "Tools/FQHEWaveFunction/JainCFFilledLevelOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/LandauSpectrumOnSphere.h"


class ConfigurationParser;


class JainCFOnSphereWaveFunction: public JainCFFilledLevelOnSphereWaveFunction
{

 protected:

  // description of the  occupation of the pseudo-Landau levels
  LandauSpectrumOnSphere* LevelOccupation;

  // number of wave functions that appear in the linear combination
  int NbrLinearCombination;
  // coefficients in front of each wave functions that appear in the linear combination  
  double* LinearCombinationCoefficients;

 public:

  // default constructor
  //
  JainCFOnSphereWaveFunction();

  // constructor
  //
  // filename = name of the file describing the occupation of the pseudo-Landau levels
  JainCFOnSphereWaveFunction(char* filename);

  // copy constructor
  //
  // function = reference on the wave function to copy
  JainCFOnSphereWaveFunction(const JainCFOnSphereWaveFunction& function);

  // destructor
  //
   ~JainCFOnSphereWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  virtual Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  virtual Complex operator ()(RealVector& x);

  // evaluate function at a given point
  //
  // uv = ensemble of spinor variables on sphere describing point
  //      where function has to be evaluated
  //      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
  // return value = function value at (uv)
  Complex CalculateFromSpinorVariables(ComplexVector& uv);

 protected:

  // get occupation information from a formatted string
  //
  // descriptionString = pointer to the string containing the description
  // descriptionArray = reference on the array where description has to be stored
  // return value = number of particles (0 if an error occured)
  int ParseOccupationDescription (char* descriptionString, int**& descriptionArray);

  // parse general informations about the composite fermion state
  // 
  // state = reference on the configuration parser that contains the informations
  // return value = false if an error occured  
  bool ParseGeneralInformation(ConfigurationParser& state);

};

#endif
