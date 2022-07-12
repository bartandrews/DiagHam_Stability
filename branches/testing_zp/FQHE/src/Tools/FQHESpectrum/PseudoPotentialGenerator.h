////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class for calculation of Pseudopotential Coefficients  //
//                                                                            //
//                        last modification : 19/11/2007                      //
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


#ifndef PSEUDOPOTENTIALGENERATOR_H
#define PSEUDOPOTENTIALGENERATOR_H


#include "config.h"
#include "MathTools/ClebschGordanCoefficients.h"


class PseudoPotentialGenerator
{
private:
  // main Vector couplings
  ClebschGordanCoefficients *MainCoefficients;
  
  // individual Vector couplings
  ClebschGordanCoefficients *Coefficients;
  
  // new formfactors for finite thickness/separation:
  double *FormFactors;

  // number of flux
  int NbrFlux;

  // Landau level index
  int LandauLevel;

  // maximum momentum index
  int MaxMomentum;
  
  // separation of layers
  double LayerSeparation;
  
public:

  // default constructor
  PseudoPotentialGenerator();

  // constructor
  // nbrFlux = nbr Flux of sphere
  // landauLevel = LL index
  // layerSeparation = layer separation in magnetic lengths
  PseudoPotentialGenerator(int nbrFlux, int landauLevel, double layerSeparation=0.0);

  // default destructor
  ~PseudoPotentialGenerator();

  // evalute a pseudopotential for coulomb interaction in a given Landau level
  //
  // relativeM = relative angular momentum value
  // layerSeparation = layer separation d in bilayer, or layer thickness d modeled by interaction 1/sqrt(r^2+d^2)
  // return value = V_m(d)
  //
  double Evaluate(int relativeM, double layerSeparation);

  // get a particular pseudopotential without changing d
  // 
  double Evaluate(int relativeM);

 private:

  // auxiliary function to evaluate form-factors
  void EvaluateFormFactors();
  
};

#endif
