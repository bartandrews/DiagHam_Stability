////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//     SU(2) spin and a generic interaction defined by its pseudopotential    //
//                                                                            //
//                        last modification : 07/06/2007                      //
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


#ifndef PARTICLEONSPHERETWOLANDAULEVELHAMILTONIAN_H
#define PARTICLEONSPHERETWOLANDAULEVELHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/AbstractQHEOnSphereWithSpinFullHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class AbstractArchitecture;


class ParticleOnSphereTwoLandauLevelHamiltonian : public AbstractQHEOnSphereWithSpinFullHamiltonian
{

  friend class QHEParticlePrecalculationOperation;

 protected:

  // array with the pseudo-potentials (ordered such that the last element corresponds to the delta interaction)
  // first index refered to the spin sector (sorted as up-up, down-down, up-down)
  double** PseudoPotentials;

  // Total cyclotron energy difference between the two Landau levels
  double TotalCyclotronEnergy;

  // difference of indices between the lower and higher Landau levels
  int LandauLevelIndexDifference;
  // twice the maximum momentum a particle can reach in the upper Landau level
  int LzMaxUp;
  // twice the maximum momentum a particle can reach in the lower Landau level
  int LzMaxDown;


 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state in the lower Landau level
  // landauLevelIndexDifference = difference of indices between the lower and higher Landau levels
  // pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
  //                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
  // cyclotronEnergy = cyclotron energy in e^2/(epsilon l_b) unit
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereTwoLandauLevelHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, int landauLevelIndexDifference, 
					    double** pseudoPotential, double cyclotronEnergy,					    
					    AbstractArchitecture* architecture, long memory = -1, 
					    bool onDiskCacheFlag = false, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnSphereTwoLandauLevelHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();


};

#endif
