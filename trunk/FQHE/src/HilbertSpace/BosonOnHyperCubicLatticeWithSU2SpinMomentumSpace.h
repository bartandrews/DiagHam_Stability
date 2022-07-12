////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//            class of bosons on a 4D hypercubic lattice with SU(2) spin      //
//                                in momentum space                           //
//                                                                            //
//                        last modification : 04/04/2012                      //
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


#ifndef BOSONONHYPERCUBICLATTICEWITHSU2SPINMOMENTUMSPACE_H
#define BOSONONHYPERCUBICLATTICEWITHSU2SPINMOMENTUMSPACE_H

#include "config.h"
#include "HilbertSpace/BosonOnCubicLatticeWithSU2SpinMomentumSpace.h"

#include <iostream>



class BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace : public BosonOnCubicLatticeWithSU2SpinMomentumSpace
{

 protected:

  // number of sites in the t direction
  int NbrSiteT;
  // number of sites in the direction perpendicular to X
  int NbrSiteYZT;
  // number of sites in the direction perpendicular to XY
  int NbrSiteZT;

  // momentum along the t direction
  int KtMomentum;

  // flag to indicate that the Hilbert space should preserve Sz
  bool SzFlag;

 public:

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSiteZ = number of sites in the z direction
  // nbrSiteT = number of sites in the t direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // kzMomentum = momentum along the z direction
  // ktMomentum = momentum along the t direction
  // memory = amount of memory granted for precalculations
  BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace (int nbrBosons, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int nbrSiteT,
						    int kxMomentum, int kyMomentum, int kzMomentum, int ktMomentum, unsigned long memory = 10000000);

  // basic constructor when Sz is conserved
  // 
  // nbrBosons = number of bosons
  // nbrSpinUp = number of particles with spin up
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSiteZ = number of sites in the z direction
  // nbrSiteT = number of sites in the t direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // kzMomentum = momentum along the z direction
  // ktMomentum = momentum along the t direction
  // memory = amount of memory granted for precalculations
  BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace (int nbrBosons, int nbrSpinUp, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int nbrSiteT,
						    int kxMomentum, int kyMomentum, int kzMomentum, int ktMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace(const BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace& bosons);

  // destructor
  //
  ~BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace& operator = (const BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = kx sector in which the density matrix has to be evaluated 
  // kySector = ky sector in which the density matrix has to be evaluated 
  // kzSector = kz sector in which the density matrix has to be evaluated 
  // ktSector = kt sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, int kzSector, int ktSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate a density matrix of a subsystem of the whole system described by a given sum of projectors, using particle partition. The density matrix is only evaluated in a given momentum sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = kx sector in which the density matrix has to be evaluated 
  // kySector = kx sector in which the density matrix has to be evaluated 
  // kzSector = kx sector in which the density matrix has to be evaluated 
  // ktSector = kt sector in which the density matrix has to be evaluated 
  // nbrGroundStates = number of projectors
  // groundStates = array of degenerate groundstates associated to each projector
  // weights = array of weights in front of each projector
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, int kzSector, int ktSector, 
									 int nbrGroundStates, ComplexVector* groundStates, double* weights, AbstractArchitecture* architecture = 0);

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentKz = current momentum along z for a single particle
  // currentKt = current momentum along z for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalKz = current total momentum along z
  // currentTotalKt = current total momentum along t
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentKz, int currentKt, 
					     int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt);

  // evaluate Hilbert space dimension with a fixed number of bosons with spin up
  //
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentKz = current momentum along z for a single particle
  // currentKt = current momentum along z for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalKz = current total momentum along z
  // currentTotalKt = current total momentum along t
  // nbrSpinUp = number of particles with spin up
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentKz, int currentKt,
					     int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt, int nbrSpinUp);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentKz = current momentum along z for a single particle
  // currentKt = current momentum along z for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalKz = current total momentum along z
  // currentTotalKt = current total momentum along t
  // currentFermionicPositionUp = current fermionic position within the state description for the spin up
  // currentFermionicPositionDown = current fermionic position within the state description for the spin down
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentKz, int currentKt, 
			      int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt, 
			      int currentFermionicPositionUp, int currentFermionicPositionDown, long pos);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentKz = current momentum along z for a single particle
  // currentKt = current momentum along z for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalKz = current total momentum along z
  // currentTotalKt = current total momentum along t
  // currentFermionicPositionUp = current fermionic position within the state description for the spin up
  // currentFermionicPositionDown = current fermionic position within the state description for the spin down
  // nbrSpinUp = number of particles with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentKz, int currentKt, 
			      int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt, 
			      int currentFermionicPositionUp, int currentFermionicPositionDown, int nbrSpinUp, long pos);

};


#endif


