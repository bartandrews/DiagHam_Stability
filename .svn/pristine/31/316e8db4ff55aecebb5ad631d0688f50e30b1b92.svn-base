////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                          DiagHam  version 0.01                             //
//                                                                            //
//                     Copyright (C) 2007 Gunnar Möller                       //
//                                                                            //
//                                                                            //
//           class of Jain composite fermion wave function on sphere          //
//                      with filled (pseudo) Landau levels                    //
//                                                                            //
//                        last modification : 18/05/2007                      //
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


#ifndef SLBSWAVEFUNCTIONUNPROJECTED_H
#define SLBSWAVEFUNCTIONUNPROJECTED_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "Tools/FQHEWaveFunction/JainCFOnSphereOrbitals.h"

class SLBSWavefunctionUnprojected: public Abstract1DComplexFunctionOnSphere
{

 protected:

  JainCFOnSphereOrbitals *OrbitalFactory;

  bool NegativeFieldFlag;
  
  int NbrParticles;
  int EffectiveFlux;

  // number of Landau levels filled with composite fermions
  int NbrLandauLevels;

  // twice the value of the momentum in the lowest pseudo-Landau level
  int TwiceS;

  // power to which the Jastrow factor has to be raised
  int ActualJastrowPower;


  // array containing prefactors of each projected monopole harmonic
  double** NormalizationPrefactors;

  // array containing constant factors that appears in the sum of projected monopole harmonic (except LLL)
  double*** SumPrefactors;

    // for temporary storage:
  Complex * SpinorUCoordinates;
  Complex * SpinorVCoordinates;
  // temporary array used to store powers of the u spinor coordinates
  Complex** SpinorUCoordinatePower;
  // temporary array used to store powers of the v spinor coordinates
  Complex** SpinorVCoordinatePower;
  // temporary array  used to store f(a,b) = S_k' (u_i v_k - u_k v_i)^-(a+b) * (u_i v_k)^a (u_k v_i)^b factors that appear in the CF monopole spherical harmonic 
  Complex*** DerivativeFactors;
  // a duplicate array of DerivativeFactors used for precalculations
  Complex*** DerivativeFactors2;


  GarbageFlag Flag;
  
  // factor used to multiply each element of the Slater determinant
  double SlaterElementNorm;

  ComplexMatrix Orbitals;
  ComplexMatrix *Slater;

  // for Pfaffian part
  ComplexSkewSymmetricMatrix *Pfaffian;
  
  
  // For internal communication with AdaptNorm:
  Complex Determinant;
  Complex PfaffianValue;

  // factor of a filled LL;
  Complex Chi1;

  // single-particle Jastrow factors
  Complex **Jij;

  // if particles are very close to each other, interpolation occurs in JainCFOrbitals
  // this variable is used to pass on this value between the different subroutines
  // double Interpolation;

  // temporary array used to store (u_i v_j - u_j v_i)^-1 factors
  Complex** JastrowFactorElements;
  
 public:

  // default constructor
  //
  SLBSWavefunctionUnprojected();

  // constructor
  //
  // nbrParticles = number of particles
  // nbrLandauLevels = number of LL's of CF orbitals to fill (unprojected)
  // negFluxFlag = indicating sign of Jain wavefunction component
  SLBSWavefunctionUnprojected(int nbrParticles, int nbrLandauLevels=2, bool negFluxFlag=false);

  // copy constructor
  //
  // function = reference on the wave function to copy
  SLBSWavefunctionUnprojected(const SLBSWavefunctionUnprojected& function);

  // destructor
  //
  ~SLBSWavefunctionUnprojected();

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


  // normalize the wave-function to one for the given particle positions
  // x = point where the function has to be evaluated
  void AdaptNorm(RealVector& x);

  // normalize the wave-function over an average number of MC positions
  void AdaptAverageMCNorm(int thermalize=100, int average=250);

 private:

  // evaluate monopole spherical harmonic 
  //
  // coordinate = index of the coordinate
  // momentum = monopole spherical harmonic Lz momentum (plus S shift)
  // landauLevel = index of the Landau level
  // maximumMomentum = maxixum momentum that can be reached in the Landau level
  // return value = value of the monopole spherical harmonic at the given point
  Complex EvaluateMonopoleHarmonic (int coordinate, int momentum, int landauLevel, int maximumMomentum);

  // evaluate constant factors that appears in the sum of projected monopole harmonic (except LLL)
  //
  void EvaluateSumPrefactors();

  // evaluate normalization factors of projected monopole harmonics
  //
  void EvaluateNormalizationPrefactors();

  
  // evaluate precalculation tables used during wave function evaluation (called at each evaluation)
  //
  // requires SpinorUCoordinates and SpinorVCoordinates to be initialized prior to call
  // derivativeFlag = indicate if precalculation tables invloved in derivative evaluation have to be calculated
  // return value = value of the Jastrow factor
  
  Complex EvaluateTables(bool derivativeFlag);

};




#endif //SLBSWAVEFUNCTIONUNPROJECTED_H
