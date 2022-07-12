////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//       coulomb interaction, a double gate and magnetic translations         //
//                                                                            //
//                        last modification : 30/10/2019                      //
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


#ifndef PARTICLEONTORUSDOUBLEGATTEDCOULOMBWITHMAGNETICTRANSLATIONSHAMILTONIAN_H
#define PARTICLEONTORUSDOUBLEGATTEDCOULOMBWITHMAGNETICTRANSLATIONSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian.h"
#include "Polynomial/Polynomial.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian : public ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian
{

 protected:

  // distance of any of the two gates to the electron gas
  double GateDistance;

 public:

  // default constructor
  //
  ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // ratio = ratio between the width in the x direction and the width in the y direction
  // gateDistance = distance of any of the two gates to the electron gas (in lb units)
  // haveCoulomb = flag indicating whether a coulomb term is present
  // landauLevel = landauLevel to be simulated
  // nbrPseudopotentials = number of pseudopotentials indicated
  // pseudopotentials = pseudopotential coefficients
  // noWignerEnergy = do not consider the energy contribution from the Wigner crystal 
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum,
								       double ratio, double gateDistance, bool haveCoulomb, int landauLevel, int nbrPseudopotentials, double* pseudopotentials, bool noWignerEnergy, AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

 protected:
 
  // get fourier transform of interaction
  // Q2_half = one half of q² value
  virtual double GetVofQ(double Q2_half);

};

#endif
