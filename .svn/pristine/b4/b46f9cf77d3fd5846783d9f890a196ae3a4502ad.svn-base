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


#ifndef CFONSPHEREWITHSPINPARTONTUNNELLINGWAVEFUNCTION_H
#define CFONSPHEREWITHSPINPARTONTUNNELLINGWAVEFUNCTION_H


#ifdef HAVE_LAPACK
// you may opt out of using the LAPACK determinant routines by quoting the definition below
#define USE_LAPACK_CFCB
#endif

#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Tools/FQHEWaveFunction/JainCFOnSphereOrbitals.h"

class CFOnSphereWithSpinPartonTunnellingWaveFunction: public Abstract1DComplexFunctionOnSphere
{

 protected:

  JainCFOnSphereOrbitals *Orbitals1;
  JainCFOnSphereOrbitals *Orbitals2;

  // Fermi momentum in Up - layer
  int LFBonding;
  // Fermi momentum in Down - layer
  int LFAntiBonding;
  // maximum Fermi momentum
  int MaxLF;
  // Number of particles in Up - layer
  int NbrParticlesUp;
  // Number of particles in Down - layer
  int NbrParticlesDown;
  // Number of orbitals in bonding FS
  int NbrOrbsBonding;
  // Number of orbitals in anti-bonding FS
  int NbrOrbsAntiBonding;
  // total number of particles;
  int NbrParticles;
  
  int AbsEffectiveFlux;

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
/*   Complex *J12; */
/*   Complex *J21; */
  Complex *J22;

  // precalculated inter-layer distances
  //  Complex **InterLayerDistances;

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
  CFOnSphereWithSpinPartonTunnellingWaveFunction();

  // constructor
  //
  // nbrParticles = total number of particles
  // lFBonding = Fermi momentum in bonding FS
  // lFAntiBonding = Fermi momentum in antibonding FS
  // nbrEffectiveFlux = number of flux quanta of the magnetic monopole field experienced by CF's
  // jastrowPower = power to which the Jastrow factor has to be raised
  // nbrParticlesUp = number of particles in upper layer (-1 if half filling each)
  CFOnSphereWithSpinPartonTunnellingWaveFunction(int nbrParticles, int lFBonding, int lFAntiBonding, int nbrEffectiveFlux,
						 int jastrowPower=2, int nbrParticlesUp=-1);

  // copy constructor
  //
  // function = reference on the wave function to copy
  CFOnSphereWithSpinPartonTunnellingWaveFunction(const CFOnSphereWithSpinPartonTunnellingWaveFunction& function);

  // destructor
  //
  ~CFOnSphereWithSpinPartonTunnellingWaveFunction();

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
  // Complex GetForOtherParameters( double *coefficients);

  // do many evaluations, storing the result in the vector results given in the call
  // x: positions to evaluate the wavefuntion in
  // format for passing parameters in the matrix coefficients coefficients[nbrSet][LandauLevel],
  // the entry [][NbrLandauLevels] corresponds to the MooreRead Term.
  // virtual void GetForManyParameters(ComplexVector &results, RealVector& x, double **coefficients);
  
  // set new values of the trial coefficients (keeping the number of LL's)
  // virtual void SetTrialParameters(double * coefficients);

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

inline double CFOnSphereWithSpinPartonTunnellingWaveFunction::fsgn(int x)
{
  if (x & 1) 
    return -1.0;
  else 
    return 1.0;
}

#endif //CFONSPHEREWITHSPINPARTONTUNNELLINGWAVEFUNCTION
