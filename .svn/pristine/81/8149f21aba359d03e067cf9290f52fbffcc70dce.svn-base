////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of Pfaffian wave function with two quasiholes on disk          //
//                                                                            //
//                        last modification : 23/10/2008                      //
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


#ifndef PFAFFIANONDISKTWOQUASIHOLEWAVEFUNCTION_H
#define PFAFFIANONDISKTWOQUASIHOLEWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"


class PfaffianOnDiskTwoQuasiholeWaveFunction: public Abstract1DComplexFunction
{

 protected:

  // number of particles
  int NbrParticles;

  // position of the first quasihole
  Complex ZHole1;
  // position of the second quasihole
  Complex ZHole2;

  // Flag for bosons/fermions
  bool FermionFlag;

  // temporary array where the Pfaffian has to be stored
  ComplexSkewSymmetricMatrix TmpPfaffian;

 public:

  // constructor
  //
  // nbrParticles = number of particles
  // zHole1 = position of the first quasihole
  // zHole2 = position of the second quasihole (spherical coordinates, theta angle)
  // fermions = flag indicating whether to calculate bosonic or fermionic pfaffian
  PfaffianOnDiskTwoQuasiholeWaveFunction(int nbrParticles, Complex zHole1, Complex zHole2, bool fermions=false);

  // copy constructor
  //
  // function = reference on the wave function to copy
  PfaffianOnDiskTwoQuasiholeWaveFunction(const PfaffianOnDiskTwoQuasiholeWaveFunction& function);

  // destructor
  //
   ~PfaffianOnDiskTwoQuasiholeWaveFunction();

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
