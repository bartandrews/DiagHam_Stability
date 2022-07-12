////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on cubic lattice with spin             //
//                                  in momentum space                         //
//                                                                            //
//                        last modification : 23/06/2011                      //
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


#ifndef FERMIONONCUBICLATTICEWITHSPINMOMENTUMSPACE_H
#define FERMIONONCUBICLATTICEWITHSPINMOMENTUMSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"

#include <iostream>



class FermionOnCubicLatticeWithSpinMomentumSpace : public FermionOnSphereWithSpin
{

 protected:

  // number of sites in the x direction
  int NbrSiteX;
  // number of sites in the y direction
  int NbrSiteY;
  // number of sites in the z direction
  int NbrSiteZ;
  // number of sites in the plane perpendicular to X
  int NbrSiteYZ;

  // momentum along the x direction
  int KxMomentum;
  // momentum along the y direction
  int KyMomentum;
  // momentum along the z direction
  int KzMomentum;

  // flag to indicate that the Hilbert space should preserve Sz
  bool SzFlag;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSiteZ = number of sites in the z direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // kzMomentum = momentum along the z direction
  // memory = amount of memory granted for precalculations
  FermionOnCubicLatticeWithSpinMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int kxMomentum, int kyMomentum, int kzMomentum, unsigned long memory = 10000000);

  // basic constructor when Sz is preserved
  // 
  // nbrFermions = number of fermions
  // nbrSpinUp = number of particles with spin up
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSiteZ = number of sites in the z direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // kzMomentum = momentum along the z direction
  // memory = amount of memory granted for precalculations
  FermionOnCubicLatticeWithSpinMomentumSpace (int nbrFermions, int nbrSpinUp, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int kxMomentum, int kyMomentum, int kzMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnCubicLatticeWithSpinMomentumSpace(const FermionOnCubicLatticeWithSpinMomentumSpace& fermions);

  // destructor
  //
  ~FermionOnCubicLatticeWithSpinMomentumSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnCubicLatticeWithSpinMomentumSpace& operator = (const FermionOnCubicLatticeWithSpinMomentumSpace& fermions);

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
  // kySector = kx sector in which the density matrix has to be evaluated 
  // kzSector = kx sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, int kzSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate a density matrix of a subsystem of the whole system described by a given sum of projectors, using particle partition. The density matrix is only evaluated in a given momentum sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = kx sector in which the density matrix has to be evaluated 
  // kySector = kx sector in which the density matrix has to be evaluated 
  // kzSector = kx sector in which the density matrix has to be evaluated 
  // nbrGroundStates = number of projectors
  // groundStates = array of degenerate groundstates associated to each projector
  // weights = array of weights in front of each projector
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, int kzSector, 
									 int nbrGroundStates, ComplexVector* groundStates, double* weights, AbstractArchitecture* architecture = 0);

  // apply the inversion symmetry i.e (k_x,k_y,k_z)->(-k_x,-k_y,-k_z) to a state 
  //
  // inputstate = reference on the input state
  // inputSpace = pointer to the Hilbert space associated to the input state
  // return value = resulting state 
  ComplexVector InversionSymmetry(ComplexVector& state, FermionOnCubicLatticeWithSpinMomentumSpace* inputSpace);

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentKz = current momentum along z for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalKz = current total momentum along z
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz);

  // evaluate Hilbert space dimension with a fixed number of fermions with spin up
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentKz = current momentum along z for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalKz = current total momentum along z
  // nbrSpinUp = number of fermions with spin up
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz, int nbrSpinUp);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentKz = current momentum along z for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalKz = current total momentum along z
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz, long pos);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentKz = current momentum along z for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalKz = current total momentum along z
  // nbrSpinUp = number of fermions with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz, int nbrSpinUp, long pos);

  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in complementary Hilbert space
  // nbrIndex = number of indices to consider in complementary Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
								  ComplexVector& groundState,  HermitianMatrix* densityMatrix);

  // core part of the evaluation density matrix particle partition calculation involving a sum of projetors 
  // 
  // minIndex = first index to consider in source Hilbert space
  // nbrIndex = number of indices to consider in source Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // nbrGroundStates = number of projectors
  // groundStates = array of degenerate groundstates associated to each projector
  // weights = array of weights in front of each projector
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
								  int nbrGroundStates, ComplexVector* groundStates, double* weights, HermitianMatrix* densityMatrix);

};


#endif


