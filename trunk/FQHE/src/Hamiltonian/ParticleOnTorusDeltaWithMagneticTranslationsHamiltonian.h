////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                  delta interaction and magnetic translations               //
//                                                                            //
//                        last modification : 18/11/2003                      //
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


#ifndef PARTICLEONTORUSDELTAWITHMAGNETICTRANSLATIONSHAMILTONIAN_H
#define PARTICLEONTORUSDELTAWITHMAGNETICTRANSLATIONSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorus.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnTorusDeltaWithMagneticTranslationsHamiltonian : public ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian
{

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // ratio = ratio between the width in the x direction and the width in the y direction
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusDeltaWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum,  int xMomentum,
							  double ratio, AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnTorusDeltaWithMagneticTranslationsHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

 private:
 
  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4);

  // evaluate Wigner crystal energy per particle
  //
  // return value = Wigner crystal energy per particle
  double EvaluateWignerCrystalEnergy ();

  // evaluate Misra function (integral of t^n exp (-xt) between 1 and +inf)
  //
  // n = index of the Misra function
  // x = point where the function has to be evaluated (> 0)
  // return value = value of the n-Misra function at x
  double MisraFunction (double n, double x);

  // evaluate part of the integral needed in the Misra function (integral of t^n exp (-xt) between min and max)
  //
  // n = index of the Misra function
  // x = point where the function has to be evaluated (> 0)
  // min = lower bound of the integral
  // max = upper bound of the integral
  // nbrSubdivision = number of subdivision used for the integral
  // return value = value of the integral
  double PartialMisraFunction (double n, double x, double min, double max, int nbrSubdivision);

};

#endif
