////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of function basis for particle on sphere             //
//                                                                            //
//                        last modification : 10/12/2002                      //
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


#ifndef PARTICLEONCYLINDERFUNCTIONBASIS_H
#define PARTICLEONCYLINDERFUNCTIONBASIS_H


#include "config.h"
#include "FunctionBasis/AbstractFunctionBasis.h"


class ParticleOnCylinderFunctionBasis: public AbstractFunctionBasis
{

 protected:

  // the maximum momentum value reached by a particle
  int MaxMomentum;

  //Landau level index
  int LandauLevel;

  //ratio between the length and circuference of a cylinder  
  double Ratio;
  // cylinder perimeter
  double Perimeter;
  // kappa factor (i.e. 2 PI / Perimeter)
  double Kappa;

  // wavefunction normalization factor
  double Normalization;

  // shift that is applied to the wavefunction index when evaluting position
  double IndexShift;

 public:

  // constructor
  //
  // maxMomentum = maximum momentum reached by a particle
  // landauLevel = Landau level index
  // ratio = aspect ratio of the cylinder
  // indexShiftFlag = true if apply a (maxMomentum / 2) shift to state index
  ParticleOnCylinderFunctionBasis(int maxMomentum, int landauLevel, double ratio, bool indexShiftFlag = true);

  // get value of the i-th function at a given point (for functions which take values in C)
  //
  // x, y = coordinates where the function should be evaluated
  // index = the function index 
  // returns value of the wavefunction
  Complex GetFunctionValue(double x, double y, double index);

  // get value of the i-th function at a given point (for functions which take values in C)
  //
  // x, y = coordinates where the function should be evaluated
  // index = the function index 
  // returns value of the wavefunction
  Complex GetFunctionValue(double x, double y, int index);

};

#endif
