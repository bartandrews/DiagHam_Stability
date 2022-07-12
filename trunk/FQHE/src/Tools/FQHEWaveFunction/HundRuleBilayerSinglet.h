////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2008 Gunnar Moeller                   //
//                                                                            //
//                                                                            //
//           class implementing composite fermion state with partially        //
//             filled highest CF shell for a wave function on sphere          //
//                                                                            //
//                        last modification : 16/01/2008                      //
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


#ifndef HUNDRULEBILAYERSINGLET_H
#define HUNDRULEBILAYERSINGLET_H


#include "config.h"
#include "HundRuleCFStates.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Matrix/ComplexMatrix.h"
#include "MathTools/Complex.h"



class HundRuleBilayerSinglet : public Abstract1DComplexFunctionOnSphere
{

 protected:
  // number of particles per layer
  int NbrParticlesPerLayer;
  
  // angular momentum of the single layer states
  int LPerLayer;

  // number of non-zero Clebsch-Gordan Coefficients
  int NbrCouplings;

  // position in array for states with non-zero Clebsch-Gordan Coefficients
  int *MPositions;

  // Coupling coefficients:
  double *Couplings;

  // vector to take up half of the coordinates:
  RealVector Part;

  // complex vector to take up half of the coordinates (used by CalculateFromSpinorCoordinates)
  ComplexVector PartC;

  // storage space for calculating states
  Complex* ResultsLayer1;
  Complex* ResultsLayer2;

  // the generator for states (used for both layers)
  HundRuleCFStates *CFStates;

  // GarbageFlag
  GarbageFlag Flag;

 public:

  // default constructor
  HundRuleBilayerSinglet();
  
  // standard constructor
  HundRuleBilayerSinglet(int nbrParticlesPerLayer, int nbrEffectiveFlux = 1, int jastrowP = 1);

  // copy constructor
  HundRuleBilayerSinglet(HundRuleBilayerSinglet &toCopy);

  // destructor
  ~HundRuleBilayerSinglet();

  // assignment operator
  HundRuleBilayerSinglet& operator = (HundRuleBilayerSinglet &toCopy);

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
  virtual Complex CalculateFromSpinorVariables(ComplexVector& uv);

  // set wavefunction to one for a given set of particle coordinates
  void AdaptNorm(RealVector& x);

  // utility function to set the right dynamic interval for Monte-Carlo
  void AdaptAverageMCNorm(int thermalize = 500, int average = 1000);

};

#endif
