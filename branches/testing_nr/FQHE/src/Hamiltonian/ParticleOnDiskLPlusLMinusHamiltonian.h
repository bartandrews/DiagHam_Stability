////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a disk where         //
//             the hamiltonian is reduced to the L^+ L^- operator             //
//                                                                            //
//                        last modification : 04/07/2008                      //
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


#ifndef PARTICLEONDISKLPLUSLMINUSHAMILTONIAN_H
#define PARTICLEONDISKLPLUSLMINUSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEOnSphereHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class AbstractArchitecture;


class ParticleOnDiskLPlusLMinusHamiltonian : public AbstractQHEOnSphereHamiltonian
{

  friend class QHEParticlePrecalculationOperation;

 protected:

  // number of flux quanta
  int NbrFluxQuanta;

  // multiplicative factor in front of the L^+ L^- operator in the Hamiltonian
  double LFactor;


 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // nbrFluxQuanta = number of flux quanta
  // architecture = architecture to use for precalculation
  // lFactor = multiplicative factor in front of the L^+L^- operator in the Hamiltonian
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnDiskLPlusLMinusHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, int nbrFluxQuanta,
				       AbstractArchitecture* architecture, double lFactor = 1.0, 
				       long memory = -1, bool onDiskCacheFlag = false, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnDiskLPlusLMinusHamiltonian();



 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();


};

#endif
