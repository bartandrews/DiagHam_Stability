////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class of Laughlin wave function with one quasielectron on sphere       //
//                                                                            //
//                        last modification : 30/10/2008                      //
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


#ifndef FQHESPHERELAUGHLINONEQUASIELECTRONWAVEFUNCTION_H
#define FQHESPHERELAUGHLINONEQUASIELECTRONWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"


class FQHESphereLaughlinOneQuasielectronWaveFunction: public Abstract1DComplexFunctionOnSphere
{

 protected:

  // number of particles
  int NbrParticles;

  // position of the quasielectron (spinor coordinates) and its conjugate
  Complex UElectron;
  Complex VElectron;
  Complex ConjUElectron;
  Complex ConjVElectron;
  
  // Flag for bosons/fermions
  bool FermionFlag;

  // temporary array where the Jastrow coefficients have to be stored
  Complex** TmpJastrow;
  Complex** TmpSqrJastrow;
  // temporary array used to store weights
  Complex* TmpWeights;

 public:

  // constructor
  //
  // nbrParticles = number of particles
  // theta = position of the quasielectron (spherical coordinates, theta angle)
  // phi = position of the quasielectron (spherical coordinates, phi angle)
  // fermions = flag indicating whether to calculate bosonic or fermionic laughlin wave function
  FQHESphereLaughlinOneQuasielectronWaveFunction(int nbrParticles, double theta = 0.0, double phi = 0.0, bool fermions=false);

  // copy constructor
  //
  // function = reference on the wave function to copy
  FQHESphereLaughlinOneQuasielectronWaveFunction(const FQHESphereLaughlinOneQuasielectronWaveFunction& function);

  // destructor
  //
   ~FQHESphereLaughlinOneQuasielectronWaveFunction();

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
