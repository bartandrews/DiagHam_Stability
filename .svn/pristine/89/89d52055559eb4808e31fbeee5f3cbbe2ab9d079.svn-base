////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                          class of bosons on square lattice                 //
//          in momentum space, forbidding orbital multiple occupancy          //
//                 allowing up to 128 orbitals on 64-bit machines             //
//                                                                            //
//                        last modification : 14/02/2017                      //
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


#ifndef BOSONONSQUARELATTICEMOMENTUMSPACEHARDCORELONG_H
#define BOSONONSQUARELATTICEMOMENTUMSPACEHARDCORELONG_H

#include "config.h"
#include "HilbertSpace/BosonOnSphereLongHardcore.h"

#include <iostream>



class BosonOnSquareLatticeMomentumSpaceHardcoreLong : public BosonOnSphereLongHardcore
{
  friend class BosonOnSquareLatticeWithSU2SpinMomentumSpace;
  friend class FQHESquareLatticeSymmetrizeU1U1StateOperation;


 protected:

  // number of sites in the x direction
  int NbrSiteX;
  // number of sites in the y direction
  int NbrSiteY;
  // momentum along the x direction
  int KxMomentum;
  // momentum along the y direction
  int KyMomentum;

 public:

  // default constructor
  // 
  BosonOnSquareLatticeMomentumSpaceHardcoreLong ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // memory = amount of memory granted for precalculations
  BosonOnSquareLatticeMomentumSpaceHardcoreLong (int nbrBosons, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSquareLatticeMomentumSpaceHardcoreLong(const BosonOnSquareLatticeMomentumSpaceHardcoreLong& bosons);

  // destructor
  //
  ~BosonOnSquareLatticeMomentumSpaceHardcoreLong ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSquareLatticeMomentumSpaceHardcoreLong& operator = (const BosonOnSquareLatticeMomentumSpaceHardcoreLong& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // get the list of occupied orbitals in a given state
  //
  // state = ID of the state
  // orbitals = list of orbitals to be filled
  virtual void GetOccupied(int state, int* orbitals);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate a density matrix of a subsystem of the whole system described by a given sum of projectors, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // nbrGroundStates = number of projectors
  // groundStates = array of degenerate groundstates associated to each projector
  // weights = array of weights in front of each projector
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, 
									 int nbrGroundStates, ComplexVector* groundStates, double* weights, AbstractArchitecture* architecture = 0);

  // get the occupation of a given orbital in a state, not assuming the coordinates are valid
  //
  // index = state index
  // xPosition = orbital x coordinates
  // yPosition = orbital y coordinates
  // return value = orbital occupation
  virtual unsigned long GetSafeOccupation(int index, int xPosition, int yPosition);

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentFermionicPosition = current fermionic position within the state description
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int currentFermionicPosition, long pos);


};


// get the occupation of a given orbital in a state, not assuming the coordinates are valid
//
// index = state index
// xPosition = orbital x coordinates
// yPosition = orbital y coordinates
// return value = orbital occupation

inline unsigned long BosonOnSquareLatticeMomentumSpaceHardcoreLong::GetSafeOccupation(int index, int xPosition, int yPosition)
{
  if ((xPosition < 0) || (yPosition < 0) || (xPosition >= this->NbrSiteX) || (yPosition >= this->NbrSiteY))
    {
      return 0x0ul;
    }
  else
    {
      return (((unsigned long) (this->StateDescription[index] >> ((xPosition * this->NbrSiteY) + yPosition))) & 0x1ul);
    }
}

#endif


