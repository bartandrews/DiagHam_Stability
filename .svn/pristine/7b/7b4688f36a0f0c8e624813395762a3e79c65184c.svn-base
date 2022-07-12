////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class of hamiltonian associated to particles on a cylinder with        //
//        SU(2) spin with opposite magnetic field for each species            //
//     a generic interaction defined by its pseudopotential and pairing       //
//                                                                            //
//                        last modification : 15/09/2016                      //
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


#ifndef PARTICLEONCYLINDERWITHSPINTIMEREVERSALSYMMETRYGENERICHAMILTONIANANDPAIRING_H
#define PARTICLEONCYLINDERWITHSPINTIMEREVERSALSYMMETRYGENERICHAMILTONIANANDPAIRING_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/ParticleOnSphereWithSpinTimeReversalSymmetricGenericHamiltonianAndPairing.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class AbstractArchitecture;


class ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairing : public ParticleOnSphereWithSpinTimeReversalSymmetricGenericHamiltonianAndPairing
{

  friend class QHEParticlePrecalculationOperation;

 protected:

  // Number of Pseudopotential for up-up interaction
  int NbrPseudopotentialsUpUp;
  // pseudopotential coefficients for up-up interaction
  double* PseudopotentialsUpUp;
  // Number of Pseudopotential for down-down interaction
  int NbrPseudopotentialsDownDown;
  // pseudopotential coefficients for down-down interaction
  double* PseudopotentialsDownDown;
  // Number of Pseudopotential for up-down interaction
  int NbrPseudopotentialsUpDown;
  // pseudopotential coefficients for up-down interaction
  double* PseudopotentialsUpDown;

  // ratio between the width in the x direction and the width in the y direction
  double Ratio;
  // ratio between the width in the y direction and the width in the x direction
  double InvRatio;

 public:

  ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairing();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // lzmax = maximum Lz value reached by a particle in the state
  // ratio = ratio between the width in the x direction and the width in the y direction
  // pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
  //                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
  // onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
  // onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
  // onebodyPotentialPairing =  one-body pairing term (sorted from component on the lowest Lz state to component on the highest Lz state), on site, symmetric spin up / spin down
  // chargingEnergy = factor in front of the charging energy (i.e 1/(2C))
  // averageNumberParticles = average number of particles in the system
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairing(ParticleOnSphereWithSpin* particles, int lzmax, double ratio, double** pseudoPotential,
									      double* onebodyPotentialUpUp, double* onebodyPotentialDownDown,
									      double* onebodyPotentialPairing, double chargingEnergy, double averageNumberParticles,
									      AbstractArchitecture* architecture, long memory = -1, 
									      bool onDiskCacheFlag = false, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairing();

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
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials);

};

#endif
