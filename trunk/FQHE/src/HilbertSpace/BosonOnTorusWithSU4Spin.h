////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of bosons on a torus with SU(4) spin                //
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


#ifndef BOSONONTORUSWITHSU4SPIN_H
#define BOSONONTORUSWITHSU4SPIN_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereWithSU4Spin.h"

#include <iostream>

class BosonOnTorusShort;

class BosonOnTorusWithSU4Spin :  public BosonOnSphereWithSU4Spin
{

  friend class BosonOnTorusShort;

 protected:

  // momentum along the y direction
  int KyMomentum;

 public:

  // constructor with a constraint on total spin momentum and total momentum
  // 
  // nbrBosons = number of bosons
  // totalSpin = twice the total spin value
  // totalIsospin = twice the total isospin value
  // totalEntanglement = twice the total entanglement value
  // maxMomentum = momentum maximum value for a boson
  // kyMomentum = momentum along the y direction
  // memory = amount of memory granted for precalculations
  BosonOnTorusWithSU4Spin (int nbrBosons, int totalSpin, int totalIsospin, int totalEntanglement, int maxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnTorusWithSU4Spin(const BosonOnTorusWithSU4Spin& bosons);

  // destructor
  //
  ~BosonOnTorusWithSU4Spin ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnTorusWithSU4Spin& operator = (const BosonOnTorusWithSU4Spin& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

 protected:

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // currentKyUpPlus = current momentum along y for a single particle with up-plus
  // currentKyUpMinus = current momentum along y for a single particle with up-minus
  // currentKyDownPlus = current momentum along y for a single particle with down-plus
  // currentKyDownMinus = current momentum along y for a single particle with down-minus
  // currentTotalKy = current total momentum along y
  // nbrNUpPlus = number of particles with quantum number up-plus
  // nbrNUpMinus = number of particles with quantum number up-minus
  // nbrNDownPlus = number of particles with quantum number down-plus
  // nbrNDownMinus = number of particles with quantum number down-minus
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int nbrBosons, int currentKyUpPlus, int currentKyUpMinus, int currentKyDownPlus, int currentKyDownMinus, int currentTotalKy, 
		      int nbrNUpPlus, int nbrNUpMinus, int nbrNDownPlus, int nbrNDownMinus, long pos);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // nbrNUpPlus = number of particles with quantum number up-plus
  // nbrNUpMinus = number of particles with quantum number up-minus
  // nbrNDownPlus = number of particles with quantum number down-plus
  // nbrNDownMinus = number of particles with quantum number down-minus
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy, int nbrNUpPlus, int nbrNUpMinus, int nbrNDownPlus, int nbrNDownMinus);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // nbrNUpPlusSector = number of particles with quantum number up-plus that belong to the subsytem 
  // nbrNUpMinusSector = number of particles with quantum number up-minus that belong to the subsytem 
  // nbrNDownPlusSector = number of particles with quantum number down-plus that belong to the subsytem 
  // nbrNDownMinusSector = number of particles with quantum number down-plus that belong to the subsytem 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, 
								     int nbrNUpPlusSector, int nbrNUpMinusSector, int nbrNDownPlusSector, int nbrNDownMinusSector, RealVector& groundState, AbstractArchitecture* architecture);

  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in source Hilbert space
  // nbrIndex = number of indices to consider in source Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
								  RealVector& groundState, RealSymmetricMatrix* densityMatrix);

};


#endif


