////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of Jain composite fermion wave function on sphere          //
//                      with filled (pseudo) Landau levels                    //
//                                                                            //
//                        last modification : 16/09/2004                      //
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


#ifndef JAINCFFILLEDLEVELONSPHEREWAVEFUNCTION_H
#define JAINCFFILLEDLEVELONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"


class JainCFFilledLevelOnSphereWaveFunction: public Abstract1DComplexFunctionOnSphere
{

 protected:

  // number of particles
  int NbrParticles;

  // number of Landau levels filled with composite fermions
  int NbrLandauLevels;

  // twice the value of the momentum in the lowest pseudo-Landau level
  int TwiceS;
  // flag to indicate if flux is reversed
  bool ReverseFluxFlag;

  // power to which the Jastrow factor has to be raised
  int ActualJastrowPower;

  // for fermionic case, half the power to which the Jastrow factor has to be raised, for bosonic case, half the power + 1 to which the Jastrow factor has to be raised
  int JastrowPower;

  // array with powers of JastrowPower
  double* JastrowPowerPowers;

  // array containing prefactors of each projected monopole harmonic
  double** NormalizationPrefactors;

  // array containing constant factors that appears in the sum of projected monopole harmonic (except LLL)
  double*** SumPrefactors;

  // garbage flag to avoid duplication of precalculation array
  GarbageFlag Flag;

  // temporary array used to store (u_i v_j - u_j v_i)^-1 factors
  Complex** JastrowFactorElements;
  // temporary array  used to store f(a,b) = S_k' (u_i v_k - u_k v_i)^-(a+b) * (u_i v_k)^a (u_k v_i)^b factors that appear in the CF monopole spherical harmonic 
  Complex*** DerivativeFactors;
  // a duplicate array of DerivativeFactors used for precalculations
  Complex*** DerivativeFactors2;

  // temporary array used to store u spinor coordinates
  Complex* SpinorUCoordinates;
  // temporary array used to store v spinor coordinates
  Complex* SpinorVCoordinates;
  // temporary array used to store powers of the u spinor coordinates
  Complex** SpinorUCoordinatePower;
  // temporary array used to store powers of the v spinor coordinates
  Complex** SpinorVCoordinatePower;

 public:

  // default constructor
  //
  JainCFFilledLevelOnSphereWaveFunction();

  // constructor
  //
  // nbrParticles = number of particles
  // nbrLandauLevel = number of Landau levels filled with composite fermions
  // jastrowPower = power to which the Jastrow factor has to be raised
  JainCFFilledLevelOnSphereWaveFunction(int nbrParticles, int nbrLandauLevels, int jastrowPower);

  // copy constructor
  //
  // function = reference on the wave function to copy
  JainCFFilledLevelOnSphereWaveFunction(const JainCFFilledLevelOnSphereWaveFunction& function);

  // destructor
  //
  ~JainCFFilledLevelOnSphereWaveFunction();

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

 protected:

  // evaluate precalculation tables used during wave function evaluation (called at each evaluation)
  //
  // derivativeFlag = indicate if precalculation tables invloved in derivative evaluation have to be calculated
  // return value = value of the Jastrow factor
  Complex EvaluateTables(bool derivativeFlag = true);

  // evaluate normalization factors of projected monopole harmonics
  //
  void EvaluateNormalizationPrefactors();

  // evaluate constant factors that appears in the sum of projected monopole harmonic (except LLL)
  //
  virtual void EvaluateSumPrefactors();

  // evaluate composite fermion monopole spherical harmonic 
  //
  // coordinate = index of the main coordinate (aka coordinate before project onto the lowest Landau level)
  // momentum = monopole spherical harmonic Lz momentum
  // landauLevel = index of the pseudo Landau level
  // maximumMomentum = maxixum momentum that can be reached in the current pseudo Landau level
  // return value = value of the monopole spherical harmonic at the givne point
  Complex EvaluateCFMonopoleHarmonic (int coordinate, int momentum, int landauLevel, int maximumMomentum);

  // evaluate derivative part of the composite fermion monopole spherical harmonic 
  //
  // index = particle index
  // alpha = number of (d\du) derivates
  // beta = number of (d\dv) derivates
  // return value = derivative contribution to the monopole spherical harmonic  
  Complex EvaluateCFMonopoleHarmonicDerivative(int index, int alpha, int beta);

};

#endif
