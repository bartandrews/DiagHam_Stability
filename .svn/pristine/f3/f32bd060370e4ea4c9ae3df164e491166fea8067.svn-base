////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of Pfaffian wave function with two quasiholes on sphere         //
//                                                                            //
//                        last modification : 18/07/2006                      //
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


#ifndef PFAFFIANONSPHERETWOQUASIHOLEWAVEFUNCTION_H
#define PFAFFIANONSPHERETWOQUASIHOLEWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"


class PfaffianOnSphereTwoQuasiholeWaveFunction: public Abstract1DComplexFunctionOnSphere
{

 protected:

  // number of particles
  int NbrParticles;

  // position of the first quasihole (spherical coordinates)
  double Theta1;
  double Phi1;
  // position of the second quasihole (spherical coordinates)
  double Theta2;
  double Phi2;

  // position of the first quasihole (spinor coordinates)
  Complex U1;
  Complex V1;

  // position of the second quasihole (spinor coordinates)
  Complex U2;
  Complex V2;
  
  // Flag for bosons/fermions
  bool FermionFlag;

  // temporary array where the Pfaffian has to be stored
  ComplexSkewSymmetricMatrix TmpPfaffian;

 public:

  // constructor
  //
  // nbrParticles = number of particles
  // theta1 = position of the first quasihole (spherical coordinates, theta angle)
  // phi1 = position of the first quasihole (spherical coordinates, phi angle)
  // theta2 = position of the second quasihole (spherical coordinates, theta angle)
  // phi2 = position of the second quasihole (spherical coordinates, phi angle)
  // fermions = flag indicating whether to calculate bosonic or fermionic pfaffian
  PfaffianOnSphereTwoQuasiholeWaveFunction(int nbrParticles, double theta1=0.0, double phi1=0.0, double theta2=M_PI, double phi2=0.0, bool fermions=false);

  // copy constructor
  //
  // function = reference on the wave function to copy
  PfaffianOnSphereTwoQuasiholeWaveFunction(const PfaffianOnSphereTwoQuasiholeWaveFunction& function);

  // destructor
  //
   ~PfaffianOnSphereTwoQuasiholeWaveFunction();

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
  virtual Complex CalculateFromSpinorVariables(ComplexVector& uv);

};

#endif
