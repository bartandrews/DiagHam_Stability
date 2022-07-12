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


#ifndef HUNDRULECFSTATES_H
#define HUNDRULECFSTATES_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"
#include "JainCFOnSphereOrbitals.h"
#include "GeneralTools/GarbageFlag.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "MathTools/Complex.h"
#include "Vector/ComplexVector.h"
#include "SlaterSuperposition.h"


// switch if using LAPACK or native DiagHam routines for Determinants
#define WANT_LAPACK

#ifdef __LAPACK__
#ifdef WANT_LAPACK
#define  __USE_LAPACK_HERE__
#endif
#endif



class HundRuleCFStates : public Abstract1DComplexFunctionOnSphere
{

 protected:

  // number of particles
  int NbrParticles;

  // number of shells filled
  int NumShells;

  // number of orbitals in highest shell
  int NbrOrbitalsInHighestShell;

  // number of particles in highest shell
  int NbrParticlesInHighestShell;

  // maximal angular momentum of the highest shell
  int HighestShellLzMax;

  // twice the value of the momentum of the angular momentum addition
  int TotalL;

  // flag to indicate if flux is reversed
  bool ReverseFluxFlag;

  // twice the value of the flux through the sphere
  int TwoS;
  
  // twice the value of the effective flux
  int TwoSEff;

  // power of Jastrow factors div 2
  int JastrowP;

  // Generator of projected CF orbitals
  JainCFOnSphereOrbitals *Orbitals;

  // single-particle Jastrow factors
  Complex *Ji;

  // if particles are very close to each other, interpolation occurs in JainCFOrbitals
  // this variable is used to pass on this value between the different subroutines
  double Interpolation;

  // Generator of total angular momentum eigenstates (not working for the moment)
  SlaterSuperposition *LEigenstates;

  // array with final description of slater determinants:
  int *NbrTermsPerLz;
  SlaterComponent **TermsPerLz;

  // state returned when prompted with operator()
  int SelectMPosition;

  // for internal communication:
  ComplexMatrix OrbitalValues;

  // matrix to for calculation of slater determinant
#ifdef __USE_LAPACK_HERE__
  ComplexLapackDeterminant SlaterDeterminant;
#else
  ComplexMatrix SlaterDeterminant;
#endif

  // constant used to normalize function
  double ElementNorm;
  
  // GarbageFlag
  GarbageFlag Flag;
  

 public:

  // default constructor
  HundRuleCFStates();

  // standard constructor
  // nbrParticles = number of particles in system
  // nbrEffectiveFlux = effective flux seen by composite fermions
  // jastrowP = power of jastrow factors div 2
  // overrideK = calculate for value of angular momentum < maximumL (only active for 2 or more particles)
  //
  HundRuleCFStates(int nbrParticles, int nbrEffectiveFlux, int jastrowP = 1, int overrideL=-1);

  // copy constructor
  HundRuleCFStates(HundRuleCFStates &toCopy);

  // destructor
  ~HundRuleCFStates();

  // assignment operator
  HundRuleCFStates& operator = (HundRuleCFStates &toCopy);
  
  // clone function 
  //
  // return value = clone of the function 
  virtual Abstract1DComplexFunction* Clone ();

  // request angular momentum
  int GetTotalL() { return this->TotalL; }
  
  // select Lz of state returned with operator()
  // return value=true upon success (M in bounds)
  bool SelectMValue(int newMValue);

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

  // evaluate all different M states for the given particle coordinates
  //
  // x = point where the function has to be evaluated
  // value returned in vector result
  void GetValues(RealVector& x, Complex* result);

  // evaluate all different M states for the given particle coordinates
  //
  // uv = ensemble of spinor variables on sphere describing point
  //      where function has to be evaluated
  // value returned in vector result
  void GetValuesFromSpinorVariables(ComplexVector& uv, Complex *result);


  // set wavefunction to one for a given set of particle coordinates
  void AdaptNorm(RealVector& x);
  
  // utility function to set the right dynamic interval for Monte-Carlo
  void AdaptAverageMCNorm(int thermalize = 500 , int average = 1000);


 private:

  // add up the Slater determinants for a single state:
  Complex EvaluateAState(int LzPosition);

  // evaluate all precalculations
  void EvaluateTables();

  

};

#endif
