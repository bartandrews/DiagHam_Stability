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


#ifndef PAIREDCFONSPHERE2QHWAVEFUNCTION_H
#define PAIREDCFONSPHERE2QHWAVEFUNCTION_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunctionOnSphere.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "Tools/FQHEWaveFunction/JainCFOnSphereOrbitals.h"
#include "Tools/FQHEMonteCarlo/ParticleOnSphereCollection.h"

class PairedCFOnSphere2QHWaveFunction: public Abstract1DComplexTrialFunctionOnSphere
{

 protected:

  JainCFOnSphereOrbitals *Orbitals;

  double MooreReadCoefficient;
  int NbrLandauLevels;
  int NbrParticles;
  int AbsEffectiveFlux;
  int LValue;
  int MValue;

  GarbageFlag Flag;
  
  // factor used to multiply each element of the Slater matrix
  double ElementNorm;

  // slater determinant used for calculating individual function values
  ComplexSkewSymmetricMatrix *Slater;

  // matrix with pair wavefunction
  ComplexSkewSymmetricMatrix *Gij;

  // matrix to store values of orbitals
  ComplexMatrix OrbitalValues;

  // number of terms in Vector coefficients
  int NbrCouplings;

  // Vector Coupling coefficients:
  double *Couplings;

  // According angular momentum values
  int *M1Values;
  int *M2Values;

  // Basis functions for quasiholes
  ParticleOnSphereFunctionBasis *QHBasis;

  // Twice the maximum angular momentum of a quasi-hole:
  int QHLzMax;

  // number of MC steps used for internal integration over QP coordinates
  int QuasiParticleMCSteps;

  // coordinates of quasiparticles as used in Monte Carlo

  ParticleOnSphereCollection * QuasiParticles; 
  
  // single-particle Jastrow factors
  Complex *Ji;

  // quasi-hole Jastrow factor elements
  Complex *Jqh1;
  Complex *Jqh2;

  // values of quasi hole basis functions
  Complex *QH1BasisValues;
  Complex *QH2BasisValues;
  
  // precalculated sums over Orbitals
  Complex **gAlpha;

  // if particles are very close to each other, interpolation occurs in JainCFOrbitals
  // this variable is used to pass on this value between the different subroutines
  double Interpolation;
  
 public:

  // constants

  enum
    {
      DefaultConventions = false,
      AssumeOldConventions = true
    };

  // default constructor
  //
  PairedCFOnSphere2QHWaveFunction();

  // constructor
  //
  // nbrParticles = number of particles
  // nbrLandauLevel = number of Landau levels filled with composite fermions
  // nbrEffectiveFlux = number of flux quanta of the magnetic monopole field experienced by CF's
  // lValue = value of angular momentum for total state
  // mValue = value of z component of total angular momentum
  // MooreReadCoefficient = prefactor of singular 1/z term in pair-wave function
  // CFCoefficients = prefactors of CF orbitals in shells 0, 1, 2, ... , nbrLandauLevel-1
  // quasiParticleMCSteps = number of steps for internal integration over QP coordinates
  // correctPrefactors = flag that enables the correction of prefactors to adopt the conventions of previous code
  // jastrowPower = power to which the Jastrow factor has to be raised
  PairedCFOnSphere2QHWaveFunction(int nbrParticles, int nbrLandauLevels, int nbrEffectiveFlux, int lValue, int mValue, double MooreReadCoefficient,
				  double * CFCoefficients, int quasiParticleMCSteps=1000, bool correctPrefactors=false, int jastrowPower=2);

  // copy constructor
  //
  // function = reference on the wave function to copy
  PairedCFOnSphere2QHWaveFunction(const PairedCFOnSphere2QHWaveFunction& function);

  // destructor
  //
  ~PairedCFOnSphere2QHWaveFunction();

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
  virtual Complex GetForOtherParameters( double *coefficients);
 
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

  
  // finish evaluation of the function after preparing all arrays with preceding function calls
  // integrating over quasihole coordinates by Monte-Carlo
  Complex MonteCarloEvaluate();

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

inline double PairedCFOnSphere2QHWaveFunction::fsgn(int x)
{
  if (x & 1) 
    return -1.0;
  else 
    return 1.0;
}


#endif //PAIREDCFONSPHERE2QHWAVEFUNCTION
