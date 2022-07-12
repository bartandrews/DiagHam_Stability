////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//     SU(4) spin and a generic interaction defined by its pseudopotential    //
//                                                                            //
//                        last modification : 28/11/2006                      //
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


#ifndef PARTICLEONSPHEREWITHSU4SPINGENERICHAMILTONIAN_H
#define PARTICLEONSPHEREWITHSU4SPINGENERICHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSU4Spin.h"
#include "Hamiltonian/AbstractQHEOnSphereWithSU4SpinHamiltonian.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class AbstractArchitecture;


class ParticleOnSphereWithSU4SpinGenericHamiltonian : public AbstractQHEOnSphereWithSU4SpinHamiltonian
{

  friend class QHEParticlePrecalculationOperation;

 protected:

  // array with the pseudo-potentials (ordered such that the last element corresponds to the delta interaction)
  // first index refered to the spin/isospin sector (sorted as up-up, up-um, up-dp, up-dm, um-um, um-dp, um-dm, dp-dp, dp-dm and dm-dm)
  double** Pseudopotentials;

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // architecture = architecture to use for precalculation
  // pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
  //                   first index refered to the spin/isospin sector (sorted as up-up, up-um, up-dp, up-dm, um-um, um-dp, um-dm, dp-dp, dp-dm and dm-dm)
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereWithSU4SpinGenericHamiltonian(ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int lzmax, double** pseudoPotential, double zeeman,
						AbstractArchitecture* architecture, long memory = -1, 
						bool onDiskCacheFlag = false, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnSphereWithSU4SpinGenericHamiltonian();

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
  virtual void EvaluateInteractionFactors();

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // clebsch = reference to the Clebsch-Gordan coefficients
  // pseudopotentials = pseudopotential coefficients
  // return value = numerical coefficient  
  virtual double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, ClebschGordanCoefficients& clebsch, double* pseudopotentials);

};

#endif
