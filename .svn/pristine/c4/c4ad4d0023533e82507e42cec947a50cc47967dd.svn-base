////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of function basis for particle on Chern insulator         //
//                        within single band approximation                    //
//                                                                            //
//                        last modification : 02/03/2011                      //
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


#ifndef PARTICLEONCHERNINSULATORSINGLEBANDFUNCTIONBASIS_H
#define PARTICLEONCHERNINSULATORSINGLEBANDFUNCTIONBASIS_H


#include "config.h"
#include "FunctionBasis/AbstractFunctionBasis.h"


class ParticleOnChernInsulatorSingleBandFunctionBasis: public AbstractFunctionBasis
{

 protected:

  // number of sites in the x direction
  int NbrSiteX;
  // number of sites in the y direction
  int NbrSiteY;
  // band parameter
  double BandParameter;

  // cosine/sine argument along x i.e. 2 pi / NbrSiteX
  double InvNbrSiteX;
  // cosine/sine argument along y i.e. 2 pi / NbrSiteY
  double InvNbrSiteY;

  // normalization coeffiicents coming from the truncated basis
  Complex* NormalizationCoefficients;

 public:

  // constructor
  //
  ParticleOnChernInsulatorSingleBandFunctionBasis(int nbrSiteX, int nbrSiteY, double mass);

  // destructor
  //
  ~ParticleOnChernInsulatorSingleBandFunctionBasis ();

  // get value of the index-th function at a given point (for functions which take values in C)
  //
  // value = reference on the value where the function has to be evaluated
  // index = linearized momentum index i.e. kx * Ny + ky
  // result = reference on the value where the result has to be stored
  // index = function index 
  void GetFunctionValue(RealVector& value, Complex& result, int index);

};

#endif


