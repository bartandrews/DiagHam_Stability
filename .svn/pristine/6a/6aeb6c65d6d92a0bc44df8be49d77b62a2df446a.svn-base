////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of bosons on a torus with spin                   //
//                                                                            //
//                        last modification : 03/04/2012                      //
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


#ifndef BOSONONTORUSWITHSPIN_H
#define BOSONONTORUSWITHSPIN_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"

#include <iostream>

class BosonOnTorusShort;

class BosonOnTorusWithSpin :  public BosonOnSphereWithSU2Spin
{

  friend class BosonOnTorusShort;
  friend class BosonOnTorusWithSpinAndMagneticTranslations;
  friend class BosonOnTorusWithSpinAllSzAndMagneticTranslations;

 protected:

  // momentum along the y direction
  int KyMomentum;

  // flag to indicate that the Hilbert space should preserve Sz
  bool SzFlag;

 public:

  // constructor with a constraint on the total Ky momentum
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  // kyMomentum = momentum along the y direction
  // memory = amount of memory granted for precalculations
  BosonOnTorusWithSpin (int nbrBosons, int maxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // constructor with a constraint on total spin momentum and total momentum
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  // totalSpin = twice the total spin along z
  // kyMomentum = momentum along the y direction
  // memory = amount of memory granted for precalculations
  BosonOnTorusWithSpin (int nbrBosons, int maxMomentum, int totalSpin, int kyMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnTorusWithSpin(const BosonOnTorusWithSpin& bosons);

  // destructor
  //
  ~BosonOnTorusWithSpin ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnTorusWithSpin& operator = (const BosonOnTorusWithSpin& bosons);

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

  // project out any configurations that have particles on levels other than lll
  //
  // inputVector = vector to apply the projection to
  // outputVector = projected vector
  // finalSpace = reference to output vector space
  void ProjectionInTheLowestLevel(RealVector &inputVector, RealVector & outputVector, BosonOnTorusShort *finalSpace);
  void ProjectionInTheLowestLevel(ComplexVector &inputVector, ComplexVector & outputVector, BosonOnTorusShort *finalSpace);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, 
								     RealVector& groundState, AbstractArchitecture* architecture);
 // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, 
								     ComplexVector& groundState, AbstractArchitecture* architecture);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // nbrNUpSector = number of spin up  that belong to the subsytem 
  // nbrNDownSector = number of spin down  that belong to the subsytem 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, int nbrNUpSector, int nbrNDownSector, RealVector& groundState, AbstractArchitecture* architecture);

  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in source Hilbert space
  // nbrIndex = number of indices to consider in source Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
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
  long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
							  ComplexVector& groundState, HermitianMatrix* densityMatrix);

 protected:

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // currentFermionicPositionUp = current fermionic position within the state description for the spin up
  // currentFermionicPositionDown = current fermionic position within the state description for the spin down
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int nbrBosons, int currentKy, int currentTotalKy, int currentFermionicPositionUp, int currentFermionicPositionDown, long pos);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // currentFermionicPositionUp = current fermionic position within the state description for the spin up
  // currentFermionicPositionDown = current fermionic position within the state description for the spin down
  // nbrSpinUp = number of particles with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int nbrBosons, int currentKy, int currentTotalKy, 
		      int currentFermionicPositionUp, int currentFermionicPositionDown, int nbrSpinUp, long pos);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy);

  // evaluate Hilbert space dimension for a given total spin momentum
  //
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // nbrSpinUp = number of particles with spin up
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy, int nbrSpinUp);


};


#endif


