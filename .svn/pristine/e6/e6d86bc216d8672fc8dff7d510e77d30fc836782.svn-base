////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of function basis for particle on checkerboard lattice      //
//                                                                            //
//                        last modification : 29/06/2011                      //
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


#ifndef PARTICLEONCHECKERBOARDLATTICEFUNCTIONBASIS_H
#define PARTICLEONCHECKERBOARDLATTICEFUNCTIONBASIS_H


#include "config.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "Matrix/ComplexMatrix.h"


class ParticleOnCheckerboardLatticeFunctionBasis: public AbstractFunctionBasis
{

 protected:

  // number of sites in the x direction
  int NbrSiteX;
  // number of sites in the y direction
  int NbrSiteY;

  // hoping amplitude between neareast neighbor sites
  double NNHoping;
  // hoping amplitude between next neareast neighbor sites
  double NextNNHoping;
  // hoping amplitude between second next neareast neighbor sites
  double SecondNextNNHoping;
  // four times the sublattice staggered chemical potential 
  double MuS;

  // cosine/sine argument along x i.e. 2 pi / NbrSiteX
  double InvNbrSiteX;
  // cosine/sine argument along y i.e. 2 pi / NbrSiteY
  double InvNbrSiteY;

  // one body basis transformation of the band structure
  ComplexMatrix* OneBodyBasis;

 public:

  // constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // t2p = hoping amplitude between second next neareast neighbor sites
  // mus = sublattice staggered chemical potential 
  ParticleOnCheckerboardLatticeFunctionBasis(int nbrSiteX, int nbrSiteY, double t1, double t2, double t2p, double mus);
  // destructor
  //
  ~ParticleOnCheckerboardLatticeFunctionBasis ();

  // get value of the index-th function at a given point (for functions which take values in C)
  //
  // value = reference on the value where the function has to be evaluated
  // index = linearized momentum index i.e. kx * Ny + ky
  // result = reference on the value where the result has to be stored
  // index = function index
  void GetFunctionValue(RealVector& value, Complex& result, int index);

};

#endif


