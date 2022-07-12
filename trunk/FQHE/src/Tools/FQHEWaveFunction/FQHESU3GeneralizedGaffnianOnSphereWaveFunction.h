////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of SU(3) generalized Gaffnian wave function on sphere       //
//                                                                            //
//                        last modification : 19/04/2008                      //
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


#ifndef FQHESU3GENERALIZEDGAFFNIANONSPHEREWAVEFUNCTION_H
#define FQHESU3GENERALIZEDGAFFNIANONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "Tools/FQHEWaveFunction/FQHESphereSymmetrizedSUKToU1WaveFunction.h"
#include "Matrix/ComplexMatrix.h"


class ConfigurationParser;


class FQHESU3GeneralizedGaffnianOnSphereWaveFunction: public FQHESphereSymmetrizedSUKToU1WaveFunction
{

 protected:

  // coefficient of the intra-component correlations
  int MIntra;
  // coefficient of the inter-component correlations
  int MInter;

 public:

  // constructor
  //
  // nbrParticles = total number of particles 
  // mIntra =  coefficient of the intra-component correlations
  // mInter = coefficient of the inter-component correlations
  FQHESU3GeneralizedGaffnianOnSphereWaveFunction(int nbrParticles, int mIntra, int mInter);

  // copy constructor
  //
  // function = reference on the wave function to copy
  FQHESU3GeneralizedGaffnianOnSphereWaveFunction(const FQHESU3GeneralizedGaffnianOnSphereWaveFunction& function);

  // constructor from configuration file
  //
  // configuration = reference on the configuration file parser
  // errorFlag = reference on the error flag that is set to false if an error occured while retriving the configuration
  // nbrParticles = reference on the total number of particles (computed from the configuration file datas)
  // lzMax = twice the maximum angular momentum for a particle (computed from the configuration file datas)
  FQHESU3GeneralizedGaffnianOnSphereWaveFunction(ConfigurationParser& configuration, bool& errorFlag, int& nbrParticles, int& lzMax);

  // destructor
  //
   ~FQHESU3GeneralizedGaffnianOnSphereWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

 protected:

  // evaluate function at a given point(the first 2*N1 coordinates correspond to the position of the type 1 particles, 
  //                                     the following 2*N2 coordinates correspond to the position of the type 2 particles,
  //                                     last the 2*N3 coordinates correspond to the position of the type 3 particles)
  //
  // uv = ensemble of spinor variables on sphere describing point
  //      where function has to be evaluated
  //      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
  //      this method is for only for internal class usage
  // return value = function value at (uv)
  virtual Complex LocalCalculateFromSpinorVariables(ComplexVector& uv);


};

#endif
