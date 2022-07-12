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


#ifndef TWOTHIRDUNPOLARIZEDCF_H
#define TWOTHIRDUNPOLARIZEDCF_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"
#include "Matrix/ComplexMatrix.h"
#include "Tools/FQHEWaveFunction/JainCFOnSphereOrbitals.h"

class TwoThirdUnpolarizedCF: public Abstract1DComplexFunctionOnSphere
{

 protected:

  JainCFOnSphereOrbitals *OrbitalFactoryUp;
  JainCFOnSphereOrbitals *OrbitalFactoryDown;

  int NbrParticles;
  int NbrParticlesPerSpin;
  int EffectiveFlux;

  GarbageFlag Flag;
  
  // factor used to multiply each element of the Slater matrix and Cauchy Permanent
  double SlaterElementNorm;

  ComplexMatrix OrbitalsUp;
  ComplexMatrix OrbitalsDown;
  ComplexMatrix *SlaterUp;
  ComplexMatrix *SlaterDown;
  
  // For internal communication with AdaptNorm:
  Complex DeterminantUpValue;
  Complex DeterminantDownValue;

  // single-particle Jastrow factors
  Complex **Jij;

  // inter-spin Jastrow factor
  Complex InterSpinJastrow;
  
  // if particles are very close to each other, interpolation occurs in JainCFOrbitals
  // this variable is used to pass on this value between the different subroutines
  double InterpolationUp;
  double InterpolationDown;

  // for testing:
  Complex * SpinorUCoordinates;
  Complex * SpinorVCoordinates;
  
 public:

  // default constructor
  //
  TwoThirdUnpolarizedCF();

  // constructor
  //
  // nbrParticles = number of particles
  TwoThirdUnpolarizedCF(int nbrParticles);

  // copy constructor
  //
  // function = reference on the wave function to copy
  TwoThirdUnpolarizedCF(const TwoThirdUnpolarizedCF& function);

  // destructor
  //
  ~TwoThirdUnpolarizedCF();

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

  // do all precalculation operations required for a new set of positions

  void EvaluateTables();

};




#endif //TWOTHIRDUNPOLARIZEDCF
