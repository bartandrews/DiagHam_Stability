////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of bosons on a torus with SU(3) spin                //
//                                                                            //
//                        last modification : 25/06/2012                      //
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


#ifndef BOSONONTORUSWITHSU3SPIN_H
#define BOSONONTORUSWITHSU3SPIN_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereWithSU3Spin.h"

#include <iostream>

class BosonOnTorusShort;

class BosonOnTorusWithSU3Spin :  public BosonOnSphereWithSU3Spin
{

  friend class BosonOnTorusShort;

 protected:

  // momentum along the y direction
  int KyMomentum;

 public:

  // constructor with a constraint ontotal momentum
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  // kyMomentum = momentum along the y direction
  // memory = amount of memory granted for precalculations
  BosonOnTorusWithSU3Spin (int nbrBosons, int maxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // constructor with a constraint on total spin momentum and total momentum
  // 
  // nbrBosons = number of bosons
  // totalTz = twice the total Tz value
  // totalY = three time the total Y value
  // maxMomentum = momentum maximum value for a boson
  // kyMomentum = momentum along the y direction
  // memory = amount of memory granted for precalculations
  BosonOnTorusWithSU3Spin (int nbrBosons, int totalTz, int totalY, int maxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnTorusWithSU3Spin(const BosonOnTorusWithSU3Spin& bosons);

  // destructor
  //
  ~BosonOnTorusWithSU3Spin ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnTorusWithSU3Spin& operator = (const BosonOnTorusWithSU3Spin& bosons);

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
  // currentKy1 = current momentum along y for a single type 1 particle
  // currentKy2 = current momentum along y for a single type 2 particle
  // currentKy3 = current momentum along y for a single type 3 particle
  // currentTotalKy = current total momentum along y
  // nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
  // nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
  // nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int nbrBosons, int currentKy1, int currentKy2, int currentKy3, int currentTotalKy, 
		      int nbrN1, int nbrN2, int nbrN3, long pos);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle 
  // currentTotalKy = current total momentum along y
  // currentFermionicPositionKy1 = current fermionic position within the state description for the type 1 particles
  // currentFermionicPositionKy2 = current fermionic position within the state description for the type 2 particles
  // currentFermionicPositionKy3 = current fermionic position within the state description for the type 3 particles
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int nbrBosons, int currentKy, int currentTotalKy, int currentFermionicPositionKy1,
		      int currentFermionicPositionKy2, int currentFermionicPositionKy3, long pos);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
  // nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
  // nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy, int nbrN1, int nbrN2, int nbrN3);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // nbrN1Sector = number of type 1 particles  that belong to the subsytem 
  // nbrN2Sector = number of type 1 particles  that belong to the subsytem 
  // nbrN3Sector = number of type 1 particles  that belong to the subsytem 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, 
								     int nbrN1Sector, int nbrN2Sector, int nbrN3Sector, RealVector& groundState, AbstractArchitecture* architecture);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // nbrN1Sector = number of type 1 particles  that belong to the subsytem 
  // nbrN2Sector = number of type 1 particles  that belong to the subsytem 
  // nbrN3Sector = number of type 1 particles  that belong to the subsytem 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, 
								     ComplexVector& groundState, AbstractArchitecture* architecture);
  
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
								  ComplexVector& groundState, HermitianMatrix* densityMatrix);

};


#endif


