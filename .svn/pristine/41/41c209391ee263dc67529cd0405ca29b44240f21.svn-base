////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                  //
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


#ifndef PAIREDCFONSPHERESPINSINGLETWAVEFUNCTION_H
#define PAIREDCFONSPHERESPINSINGLETWAVEFUNCTION_H


#ifdef HAVE_LAPACK
// you may opt out of using the LAPACK determinant routines by quoting the definition below
#define USE_LAPACK_CFCB
#endif

#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunctionOnSphere.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Tools/FQHEWaveFunction/JainCFOnSphereOrbitals.h"

class PairedCFOnSphereSpinSingletWaveFunction: public Abstract1DComplexTrialFunctionOnSphere
{

 protected:

  JainCFOnSphereOrbitals *Orbitals1;
  JainCFOnSphereOrbitals *Orbitals2;

  int NbrLandauLevels;
  int NbrUp;
  int AbsEffectiveFlux;
  int JastrowPower;

  bool HaveBosons;

  GarbageFlag Flag;
  
  // factor used to multiply each element of the Slater matrix
  double ElementNorm;
  
#ifdef USE_LAPACK_CFCB
  ComplexLapackDeterminant *Matrix;
#else
  ComplexMatrix *Matrix;
#endif
  
  // matrices to store values of orbitals
  ComplexMatrix OrbitalValues1;
  ComplexMatrix OrbitalValues2;

  // single-particle Jastrow factors
  Complex *J11;
  Complex *J12;
  Complex *J21;
  Complex *J22;

  // precalculated sums over Orbitals
  Complex **gAlpha;

  // precalculated inter-spin Jastrow-factor
  Complex InterSpinJastrow;

  // if particles are very close to each other, interpolation occurs in JainCFOrbitals
  // this variable is used to pass on this value between the different subroutines
  double Interpolation;
  
 public:

  // default constructor
  //
  PairedCFOnSphereSpinSingletWaveFunction();

  // constructor
  //
  // nbrParticles = total number of particles 
  // nbrLandauLevel = number of Landau levels filled with composite fermions
  // nbrEffectiveFlux = number of flux quanta of the magnetic monopole field experienced by CF's
  // cfCoefficients = prefactors of CF orbitals in shells 0, 1, 2, ... , nbrLandauLevel-1
  // jastrowPower = power to which the Jastrow factor has to be raised
  PairedCFOnSphereSpinSingletWaveFunction(int nbrParticles, int nbrLandauLevels, int nbrEffectiveFlux,
					  double * cfCoefficients, int jastrowPower=2);

  // copy constructor
  //
  // function = reference on the wave function to copy
  PairedCFOnSphereSpinSingletWaveFunction(const PairedCFOnSphereSpinSingletWaveFunction& function);

  // destructor
  //
  ~PairedCFOnSphereSpinSingletWaveFunction();

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

  // get a value of the wavefunction for the last set of coordinates, but with different variational coefficients
 Complex GetForOtherParameters( double *coefficients);

  // do many evaluations, storing the result in the vector results given in the call
  // x: positions to evaluate the wavefuntion in
  // format for passing parameters in the matrix coefficients coefficients[nbrSet][LandauLevel],
  // the entry [][NbrLandauLevels] corresponds to the MooreRead Term.
  virtual void GetForManyParameters(ComplexVector &results, RealVector& x, double **coefficients);
  
  // set new values of the trial coefficients (keeping the number of LL's)
  virtual void SetTrialParameters(double * coefficients);

  // normalize the wave-function to one for the given particle positions
  // x = point where the function has to be evaluated
  void AdaptNorm(RealVector& x);

  // normalize the wave-function over an average number of MC positions
  void AdaptAverageMCNorm(int thermalize=100, int average=250);

 private:

  // do all precalculation operations required for a new set of positions

  void EvaluateTables();

  // for calculating (-1)^x
  //
  // x = x value 
  // return value = (-1)^x
  double fsgn(int x);

};

// for calculating (-1)^x
//
// x = x value 
// return value = (-1)^x

inline double PairedCFOnSphereSpinSingletWaveFunction::fsgn(int x)
{
  if (x & 1) 
    return -1.0;
  else 
    return 1.0;
}

#endif //PAIREDCFONSPHERESPINSINGLETWAVEFUNCTION
