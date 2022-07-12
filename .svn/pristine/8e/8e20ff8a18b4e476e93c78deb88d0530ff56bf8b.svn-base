////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//        generic two body interaction with an anisotropic mass term          //
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


#ifndef PARTICLEONTORUSMASSANISOTROPYGENERICHAMILTONIAN_H
#define PARTICLEONTORUSMASSANISOTROPYGENERICHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnTorusGenericHamiltonian.h"

#include <iostream>


using std::ostream;


class ParticleOnTorusMassAnisotropyGenericHamiltonian : public ParticleOnTorusGenericHamiltonian
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
  // ratio = ratio between the width in the x direction and the width in the y direction
  // anisotropy = anisotropy parameter alpha (q_g^2 = alpha q_x^2 + q_y^2 / alpha)
  // nbrPseudopotentials = number of pseudopotentials
  // pseudopotential = pseudopotential coefficients
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusMassAnisotropyGenericHamiltonian(ParticleOnTorus* particles, int nbrParticles, int maxMomentum, double ratio, 
						   double anisotropy, int nbrPseudopotentials, double* pseudopotentials,
						   AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnTorusMassAnisotropyGenericHamiltonian();


 protected:
 
  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4);

};

#endif
