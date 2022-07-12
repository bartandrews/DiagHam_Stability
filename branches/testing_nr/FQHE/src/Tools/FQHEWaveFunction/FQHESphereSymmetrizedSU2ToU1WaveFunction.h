////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of U(1) wave function obtained form symmetrization of a   //
//                        SU(2) wave function on sphere                       //
//                                                                            //
//                        last modification : 20/10/2008                      //
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


#ifndef FQHESPHERESYMMETRIZEDSU2TOU1WAVEFUNCTION_H
#define FQHESPHERESYMMETRIZEDSU2TOU1WAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"
#include "GeneralTools/GarbageFlag.h"


class FQHESphereSymmetrizedSU2ToU1WaveFunction: public Abstract1DComplexFunctionOnSphere
{

 protected:

  // number of particles
  int NbrParticles;

  // number of particles with spin up
  int NbrParticlesUp;

  // number of particles with spin down
  int NbrParticlesDown;

  // true if the final state should be a fermionic state
  bool FermionFlag;  

  // true if the summation has to be done over all the permutations of S_NbrParticles 
  bool FullySymmetrize;

  // internal storage of spinor coordinates
  Complex* SpinorUCoordinates;
  Complex* SpinorVCoordinates;

  // pointer to the base SU(2) wave function
  Abstract1DComplexFunctionOnSphere* SU2WaveFunction;
  
  // array containing description of each permutation that appears in the calculation symmetrization process
  unsigned long* Permutations;
  // number of permutations that appears in the symmetrization process
  unsigned long NbrPermutations;

  // garable flag associated to the Permutations array
  GarbageFlag Flag;

  // temporary vector to store permutation of the input uv vector
  ComplexVector TemporaryUV;

 public:

  // default constructor
  //
  FQHESphereSymmetrizedSU2ToU1WaveFunction();

  // constructor
  //
  // nbrParticles = total number of particles
  // nbrParticlesUp = number of particles with spin up
  // sU2Wavefunction = pointer to the base SU(2) wave function
  // fermionFlag = true if the final state should be a fermionic state
  FQHESphereSymmetrizedSU2ToU1WaveFunction(int nbrParticles, int nbrParticlesUp, Abstract1DComplexFunctionOnSphere* sU2Wavefunction, bool fermionFlag = false);

  // constructor from data file 
  //
  // filename = pointer to the file name that described the symmetrization procedure
  // sU2Wavefunction = pointer to the base SU(2) wave function
  // fermionFlag = true if the final state should be a fermionic state
  FQHESphereSymmetrizedSU2ToU1WaveFunction(char* filename, Abstract1DComplexFunctionOnSphere* sU2Wavefunction, bool fermionFlag = false);

  // copy constructor
  //
  // function = reference on the wave function to copy
  FQHESphereSymmetrizedSU2ToU1WaveFunction(const FQHESphereSymmetrizedSU2ToU1WaveFunction& function);

  // destructor
  //
  ~FQHESphereSymmetrizedSU2ToU1WaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  virtual Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

  // evaluate function at a given point
  //
  // uv = ensemble of spinor variables on sphere describing point
  //      where function has to be evaluated
  //      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
  // return value = function value at (uv)
  virtual Complex CalculateFromSpinorVariables(ComplexVector& uv);

  // write all permutations requested to symmetrize the SU(2) state to data file 
  //
  // filename = pointer to the file name that described the symmetrization procedure
  // return value = true if no error occured
  virtual bool WritePermutations(char* filename);

 protected:

  // get all permutations requested to symmetrize the SU(2) state from data file 
  //
  // filename = pointer to the file name that described the symmetrization procedure
  // return value = true if no error occured
  virtual bool ReadPermutations(char* filename);

  // evaluate all permutations requested to symmetrize the SU(2) state
  //
  virtual void EvaluatePermutations();

  // evaluate function at a given point(the first 2*N1 coordinates correspond to the position of the type 1 particles, 
  //                                     the following 2*N2 coordinates correspond to the position of the type 2 particles,
  //                                     last the 2*N3 coordinates correspond to the position of the type 3 particles)
  //
  // uv = ensemble of spinor variables on sphere describing point
  //      where function has to be evaluated
  //      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
  //      this method is for only for internal class usage
  // return value = function value at (uv)
  virtual Complex LocalCalculateFromSpinorVariables(ComplexVector& uv);

};

#endif
