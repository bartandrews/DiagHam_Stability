////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a cylinder with      //
//     SU(2) spin and a generic interaction defined by its pseudopotential    //
//                                                                            //
//                        last modification : 19/08/2016                      //
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


#ifndef PARTICLEONCYLINDERWITHSPINGENERICHAMILTONIAN_H
#define PARTICLEONCYLINDERWITHSPINGENERICHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/ParticleOnTorusWithSpinGenericHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class AbstractArchitecture;


class ParticleOnCylinderWithSpinGenericHamiltonian : public ParticleOnTorusWithSpinGenericHamiltonian
{

 protected:


 public:


  // default constructor
  //
  ParticleOnCylinderWithSpinGenericHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // ratio = ratio between the width in the x direction and the width in the y direction
  // architecture = architecture to use for precalculation
  // nbrPseudopotentialsUpUp = number of pseudopotentials for up-up interaction
  // pseudopotentialsUpUp = pseudopotential coefficients for up-up interaction
  // nbrPseudopotentialsDownDown = number of pseudopotentials for down-down interaction
  // pseudopotentialsDownDown = pseudopotential coefficients for down-down interaction
  // nbrPseudopotentialsUpDown = number of pseudopotentials for up-down interaction
  // pseudopotentialsUpDown = pseudopotential coefficients for up-down interaction
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  // onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
  // onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
  // onebodyPotentialUpDown =  one-body tunnelling potential (sorted from component on the lowest Lz state to component on the highest Lz state), on site, symmetric spin up / spin down
  ParticleOnCylinderWithSpinGenericHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, double ratio,
					       int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp,
					       int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown,
					       int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,
					       double spinFluxUp, double spinFluxDown, 
					       AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0, 
					       double * oneBodyPotentialUpUp = 0, double * oneBodyPotentialDownDown = 0, double * oneBodyPotentialUpDown = 0);


  // destructor
  //
  ~ParticleOnCylinderWithSpinGenericHamiltonian();


 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // nbrPseudopotentials = number of pseudopotentials
  // pseudopotentials = pseudopotential coefficients
  // spinFluxM1 = additional inserted flux for m1
  // spinFluxM2 = additional inserted flux for m2
  // spinFluxM3 = additional inserted flux for m3
  // spinFluxM4 = additional inserted flux for m4
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
						double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4);

};

#endif
