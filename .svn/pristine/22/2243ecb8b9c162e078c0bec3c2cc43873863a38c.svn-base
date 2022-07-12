////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of function basis for particle on sphere             //
//                            in a given Landau level                         //
//                                                                            //
//                        last modification : 29/11/2005                      //
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


#ifndef PARTICLEONSPHEREGENERICLLFUNCTIONBASIS_H
#define PARTICLEONSPHEREGENERICLLFUNCTIONBASIS_H


#include "config.h"
#include "FunctionBasis/AbstractFunctionBasis.h"


class ParticleOnSphereGenericLLFunctionBasis: public AbstractFunctionBasis
{

 protected:

  // twice the maximum Lz value reached by a particle in the lowest Landau level (i.e. number of flux quanta)
  int LzMax;

  // index of the Landau level to consider (0 for the lowest Landau level)
  int LandauLevel;

  // array containing normalization prefactor of each monopole harmonic
  double* NormalizationPrefactors;

  // array containing constant factors that appears in the sum of monopole harmonic
  double** SumPrefactors;
  
 public:

  // constructor
  //
  // lzMax = twice the maximum Lz value reached by a particle in the lowest Landau level (i.e. number of flux quanta)
  // landauLevel = index of the Landau level to consider (0 for the lowest Landau level)
  ParticleOnSphereGenericLLFunctionBasis(int lzMax, int landauLevel);

  // destructor
  //
  ~ParticleOnSphereGenericLLFunctionBasis ();

  // get value of the i-th function at a given point (for functions which take values in C)
  //
  // value = reference on the value where the function has to be evaluated
  // result = reference on the value where the result has to be stored
  // index = function index 
  void GetFunctionValue(RealVector& value, Complex& result, int index);

 protected:

  // evaluate normalization factors of monopole harmonics
  //
  void EvaluateNormalizationPrefactors();

  // evaluate constant factors that appears in the sum of projected monopole harmonic
  //
  void EvaluateSumPrefactors();

};

#endif


