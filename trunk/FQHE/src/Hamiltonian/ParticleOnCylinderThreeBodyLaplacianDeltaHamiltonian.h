////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                          laplacian delta interaction                       //
//                                                                            //
//                        last modification : 29/06/2010                      //
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


#ifndef PARTICLEONCYLINDERTHREEBODYLAPLACIANDELTAHAMILTONIAN_H
#define PARTICLEONCYLINDERTHREEBODYLAPLACIANDELTAHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnCylinderThreeBodyHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnCylinderThreeBodyLaplacianDeltaHamiltonian : public AbstractQHEOnCylinderThreeBodyHamiltonian
{
 protected:
  
  double HypermetricTheta1;
  double HypermetricTheta2;

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // ratio = ratio between the width in the x direction and the width in the y direction
  // theta1 = hypermetric angle theta1
  // theta2 = hypermetric angle theta2
  // confinement = amplitude of the quadratic confinement potential
  // electricFieldParameter = amplitude of the electric field along the cylinder
  // bFieldfParameter = amplitude of the magnetic field (to set the energy scale)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnCylinderThreeBodyLaplacianDeltaHamiltonian(ParticleOnSphere* particles, int nbrParticles, int maxMomentum, double ratio, double theta1, double theta2, double confinement, double electricFieldParameter, double bFieldParameter,
					   AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnCylinderThreeBodyLaplacianDeltaHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);

 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a+_m3 a_m4 a_m5 a_m6 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // m5 = fifth index
  // m6 = sixth index
  // return value = numerical coefficient
  Complex EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int m5, int m6);

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a^+_m3 a_m4 a_m5 a_m6 coupling term for bosons
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // m5 = fifth index
  // m6 = sixth index
  // return value = numerical coefficient

  Complex EvaluateInteractionCoefficientBosons(int m1, int m2, int m3, int m4, int m5, int m6);

  // Get the number of permutations of annihilation/creation indices c_n1 c_n2 c_n3 for bosons

  int NumberOfPermutations(int n1, int n2, int n3);

};

#endif
