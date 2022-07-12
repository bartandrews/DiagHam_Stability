////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//      coulomb interaction, magnetic translations, and mass anisotropy       //
//                                                                            //
//                        last modification : 14/05/2014                      //
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


#ifndef PARTICLEONTORUSCOULOMBMASSANISOTROPYWITHMAGNETICTRANSLATIONSHAMILTONIAN_H
#define PARTICLEONTORUSCOULOMBMASSANISOTROPYWITHMAGNETICTRANSLATIONSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian.h"
#include "Polynomial/Polynomial.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnTorusCoulombMassAnisotropyWithMagneticTranslationsHamiltonian : public ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian
{

 protected:

  // anisotropy parameter alpha (q_g^2 = alpha q_x^2 + q_y^2 / alpha)
  double Anisotropy;
  // invert of Anisotropy
  double InvAnisotropy;

 public:

  // constructor from default data
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // ratio = ratio between the width in the x direction and the width in the y direction
  // anisotropy = anisotropy parameter alpha (q_g^2 = alpha q_x^2 + q_y^2 / alpha)
  // haveCoulomb = flag indicating whether a coulomb term is present
  // landauLevel = landauLevel to be simulated
  // nbrPseudopotentials = number of pseudopotentials indicated
  // pseudopotentials = pseudopotential coefficients
  // noWignerEnergy = do not consider the energy contribution from the Wigner crystal 
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusCoulombMassAnisotropyWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum,
							    double ratio, double anisotropy, bool haveCoulomb, int landauLevel, 
							    int nbrPseudopotentials, double* pseudopotentials, bool noWignerEnergy, AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnTorusCoulombMassAnisotropyWithMagneticTranslationsHamiltonian();

 protected:
 
  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4);

  // get fourier transform of interaction
  //
  // Q2_half = one half of q² value
  // Q2_halfgaussian = one half of q_g² value for the one body part (i.e. including the anisotropy)
  double GetVofQ(double Q2_half, double Q2_halfgaussian);

};

#endif
