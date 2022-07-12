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


#ifndef PAIREDCFONSPHEREWAVEFUNCTION_H
#define PAIREDCFONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunctionOnSphere.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "Tools/FQHEWaveFunction/JainCFOnSphereOrbitals.h"

class PairedCFOnSphereWaveFunction: public Abstract1DComplexTrialFunctionOnSphere
{

 protected:

  JainCFOnSphereOrbitals *Orbitals;

  double MooreReadCoefficient;
  int NbrLandauLevels;
  int NbrParticles;
  int AbsEffectiveFlux;

  GarbageFlag Flag;
  
  // factor used to multiply each element of the Slater matrix
  double ElementNorm;
  
  ComplexSkewSymmetricMatrix *Slater;

  // matrix to store values of orbitals
  ComplexMatrix OrbitalValues;

  // single-particle Jastrow factors
  Complex *Ji;

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
  PairedCFOnSphereWaveFunction();

  // constructor
  //
  // nbrParticles = number of particles
  // nbrLandauLevel = number of Landau levels filled with composite fermions
  // nbrEffectiveFlux = number of flux quanta of the magnetic monopole field experienced by CF's
  // MooreReadCoefficient = prefactor of singular 1/z term in pair-wave function
  // CFCoefficients = prefactors of CF orbitals in shells 0, 1, 2, ... , nbrLandauLevel-1
  // correctPrefactors = flag that enables the correction of prefactors to adopt the conventions of previous code
  // jastrowPower = power to which the Jastrow factor has to be raised
  PairedCFOnSphereWaveFunction(int nbrParticles, int nbrLandauLevels, int nbrEffectiveFlux, double MooreReadCoefficient,
			       double * CFCoefficients, bool correctPrefactors=false, int jastrowPower=2);

  // copy constructor
  //
  // function = reference on the wave function to copy
  PairedCFOnSphereWaveFunction(const PairedCFOnSphereWaveFunction& function);

  // destructor
  //
  ~PairedCFOnSphereWaveFunction();

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

inline double PairedCFOnSphereWaveFunction::fsgn(int x)
{
  if (x & 1) 
    return -1.0;
  else 
    return 1.0;
}


#endif //PAIREDCFONSPHEREWAVEFUNCTION
