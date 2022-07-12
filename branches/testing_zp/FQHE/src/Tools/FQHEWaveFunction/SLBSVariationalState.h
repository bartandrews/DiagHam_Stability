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


#ifndef SLBSVARIATIONALSTATE_H
#define SLBSVARIATIONALSTATE_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunctionOnSphere.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "Tools/FQHEWaveFunction/JainCFOnSphereOrbitals.h"
#include "Tools/FQHEWaveFunction/PairedCFOnSphereWaveFunction.h"

class SLBSVariationalState: public Abstract1DComplexTrialFunctionOnSphere
{

 protected:

  JainCFOnSphereOrbitals *OrbitalFactory;

  int NbrParticles;
  int EffectiveFlux;

  GarbageFlag Flag;
  
  // factor used to multiply each element of the Slater determinant
  double SlaterElementNorm;

  ComplexMatrix Orbitals;
  ComplexMatrix *Slater;

  // general paired state for Pfaffian part:
  PairedCFOnSphereWaveFunction *PfaffianPart;
  
  // For internal communication with AdaptNorm:
  Complex Determinant;
  Complex PfaffianValue;

  // factor of a filled LL;
  Complex Chi1;

  // single-particle Jastrow factors
  Complex **Jij;

  // if particles are very close to each other, interpolation occurs in JainCFOrbitals
  // this variable is used to pass on this value between the different subroutines
  double Interpolation;

  // for testing:
  Complex * SpinorUCoordinates;
  Complex * SpinorVCoordinates;
  
 public:

  // default constructor
  //
  SLBSVariationalState();


  // constructor
  //
  // nbrParticles = number of particles
  // nbrVariationalLandauLevels = number of CF LL's to be calculated (= length of array 'variationalCoefficients')
  // MooreReadCoefficient = prefactor of Moore-Read 1/z term
  // variationalCoefficients = variational coefficients of pair wavefunction
  SLBSVariationalState(int nbrParticles, int nbrVariationalLandauLevels, double MooreReadCoefficient, double *variationalCoefficients);

  // copy constructor
  //
  // function = reference on the wave function to copy
  SLBSVariationalState(const SLBSVariationalState& function);

  // destructor
  //
  ~SLBSVariationalState();

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

};




#endif //SLBSVARIATIONALSTATE
