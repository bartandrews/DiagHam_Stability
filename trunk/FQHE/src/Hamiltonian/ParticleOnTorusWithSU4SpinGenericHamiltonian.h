////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                  generic two body interaction and SU(4) spin               //
//                                                                            //
//                        last modification : 26/06/2012                      //
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


#ifndef PARTICLEONTORUSGENERICWITHSU4SPINHAMILTONIAN_H
#define PARTICLEONTORUSGENERICWITHSU4SPINHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSU4Spin.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnTorusWithSU4SpinHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class Polynomial;


class ParticleOnTorusWithSU4SpinGenericHamiltonian : public AbstractQHEOnTorusWithSU4SpinHamiltonian
{

 protected:


  // array with the pseudo-potentials (ordered such that the last element corresponds to the delta interaction)
  // first index refered to the spin sector (sorted as upplus-upplus, upplus-upminus, upplus-downplus, upplus-downminus, 
  //                                                   upminus-upminus, upminus-downplus, upminus-downminus, downplus-downplus, downplus-downminus, downminus-downminus)
  double** Pseudopotentials;
  // nbrPseudoPotentials = array with the number of pseudo-potentials per interaction type
  int* NbrPseudopotentials;
  // maximum number of pseudopotentials
  int MaxNbrPseudopotentials;
  // Laguerre polynomial for the pseudopotentials
  Polynomial* LaguerrePolynomials;

  // additional inserted flux for spin 1
  double SpinFlux1;
  // additional inserted flux for spin 2
  double SpinFlux2;
  // additional inserted flux for spin 3
  double SpinFlux3;
  // additional inserted flux for spin 4
  double SpinFlux4;

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // ratio = ratio between the width in the x direction and the width in the y direction
  // nbrPseudoPotentials = array with the number of pseudo-potentials per interaction type
  // pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
  //                   first index refered to the spin sector (sorted as upplus-upplus, upplus-upminus, upplus-downplus, upplus-downminus, 
  //                                                                     upminus-upminus, upminus-downplus, upminus-downminus, downplus-downplus, downplus-downminus, downminus-downminus)
  // spinFlux1 = additional inserted flux for spin 1
  // spinFlux2 = additional inserted flux for spin 2
  // spinFlux3 = additional inserted flux for spin 3
  // spinFlux4 = additional inserted flux for spin 4
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusWithSU4SpinGenericHamiltonian(ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int maxMomentum, double ratio, 
					       int* nbrPseudopotentials, double** pseudoPotentials,
					       double spinFlux1, double spinFlux2, double spinFlux3, double spinFlux4,
					       AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);
  
  // destructor
  //
  ~ParticleOnTorusWithSU4SpinGenericHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);


 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4
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
  double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4,
					int nbrPseudopotentials, double* pseudopotentials,
					double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4);


};

#endif
