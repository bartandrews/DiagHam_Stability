////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    class of Pfaffian wave function with two quasielectrons on sphere       //
//                                                                            //
//                        last modification : 23/10/2008                      //
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


#ifndef PFAFFIANONSPHERETWOQUASIELECTRONWAVEFUNCTION_H
#define PFAFFIANONSPHERETWOQUASIELECTRONWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"


class PfaffianOnSphereTwoQuasielectronWaveFunction: public Abstract1DComplexFunctionOnSphere
{

 protected:

  // number of particles
  int NbrParticles;

  // position of the first quasielectron (spinor coordinates) and its conjugate
  Complex UElectron1;
  Complex VElectron1;
  Complex ConjUElectron1;
  Complex ConjVElectron1;

  // position of the second quasielectron (spinor coordinates) and its conjugate
  Complex UElectron2;
  Complex VElectron2;
  Complex ConjUElectron2;
  Complex ConjVElectron2;
  
  // Flag for bosons/fermions
  bool FermionFlag;

  // temporary array where the Pfaffian has to be stored
  Complex** TmpPfaffian;
  // temporary array where the sqyuare elements of the Pfaffian have to be stored
  Complex** TmpSqrPfaffian;
  // temporary array where indices are stored
  int* TmpIndexArray;
  // temporary array used to store weights
  Complex* TmpWeights1;
  Complex* TmpWeights2;

  // array containing description of each permutation that appears in the calculation symmetrization process
  unsigned long* Permutations1;
  unsigned long* Permutations2;
  // number of permutations that appears in the symmetrization process
  unsigned long NbrPermutations;
  // garable flag associated to the Permutations array
  GarbageFlag Flag;

  // index of coordinate that will be changed in the next evaluation of the wavefunction
  int NextCoordinate;

  // temporary arrays and variable use to back-up coefficients involving a previous coordinate
  Complex* TmpPreviousPfaffian;
  Complex* TmpPreviousSqrPfaffian;
  Complex TmpPreviousWeights1;
  Complex TmpPreviousWeights2;
 


 public:

  // constructor
  //
  // nbrParticles = number of particles
  // theta1 = position of the first quasielectron (spherical coordinates, theta angle)
  // phi1 = position of the first quasielectron (spherical coordinates, phi angle)
  // theta2 = position of the second quasielectron (spherical coordinates, theta angle)
  // phi2 = position of the second quasielectron (spherical coordinates, phi angle)
  // fermions = flag indicating whether to calculate bosonic or fermionic pfaffian
  PfaffianOnSphereTwoQuasielectronWaveFunction(int nbrParticles, double theta1=0.0, double phi1=0.0, double theta2=M_PI, double phi2=0.0, bool fermions=false);

  // constructor using permutation description stored in a file
  //
  // filename = pointer to the file name that described the symmetrization procedure
  // theta1 = position of the first quasielectron (spherical coordinates, theta angle)
  // phi1 = position of the first quasielectron (spherical coordinates, phi angle)
  // theta2 = position of the second quasielectron (spherical coordinates, theta angle)
  // phi2 = position of the second quasielectron (spherical coordinates, phi angle)
  // fermions = flag indicating whether to calculate bosonic or fermionic pfaffian
  PfaffianOnSphereTwoQuasielectronWaveFunction(char* filename, double theta1 = 0.0, double phi1 = 0.0, 
					       double theta2 = M_PI, double phi2 = 0.0, bool fermions = false);
  // copy constructor
  //
  // function = reference on the wave function to copy
  PfaffianOnSphereTwoQuasielectronWaveFunction(const PfaffianOnSphereTwoQuasielectronWaveFunction& function);

  // destructor
  //
   ~PfaffianOnSphereTwoQuasielectronWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

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

  // write all permutations requested to symmetrize the state to data file 
  //
  // filename = pointer to the file name that described the symmetrization procedure
  // return value = true if no error occured
  virtual bool WritePermutations(char* filename);

 protected:

  // get all permutations requested to symmetrize the state from data file 
  //
  // filename = pointer to the file name that described the symmetrization procedure
  // return value = true if no error occured
  virtual bool ReadPermutations(char* filename);

  // evaluate all permutations requested to symmetrize the state
  //
  virtual void EvaluatePermutations();

};

#endif
