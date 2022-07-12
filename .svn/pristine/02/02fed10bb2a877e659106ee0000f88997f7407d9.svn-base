////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Cecile Repellin                 //
//                                                                            //
//                                                                            //
//                          class of bosons on th 4D sphere                   //
//             with a number of orbitals up to 128 (64 bits system)           //
//                                                                            //
//                        last modification : 9/11/2012                      //
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


#ifndef BOSONON4DSPHERELONG_H
#define BOSONON4DSPHERELONG_H

#include "config.h"
#include "HilbertSpace/BosonOnSphereLong.h"

#include <iostream>



class BosonOn4DSphereLong: public BosonOnSphereLong
{

  friend class FQHESquareLatticeSymmetrizeU1U1StateOperation;
  
 protected:

  // Number of flux quanta
  int NbrFluxQuanta;
  // total value of jz
  int TotalJz;
  // shifted value of total jz
  int ShiftedTotalJz;
  // total value of kz
  int TotalKz;
  // shifted value of total kz
  int ShiftedTotalKz;
  // array that gives the value of j for one particle corresponding to the linearized index
  int* quantumNumberJ;
  // array that gives the value of jz for one particle corresponding to the linearized index
  int* quantumNumberJz;
  // array that gives the value of kz for one particle corresponding to the linearized index
  int* quantumNumberKz;
  
 public:

  // default constructor
  // 
  BosonOn4DSphereLong ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrFluxQuanta = number of flux quanta (p)
  // totalJz = total value of jz
  // totalKz = total value of kz
  // memory = amount of memory granted for precalculations
  BosonOn4DSphereLong (int nbrBosons, int nbrFluxQuanta, int totalJz, int totalKz, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy 
  BosonOn4DSphereLong(const BosonOn4DSphereLong& bosons);

  // destructor
  //
  ~BosonOn4DSphereLong ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOn4DSphereLong& operator = (const BosonOn4DSphereLong& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);
/*
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given momentum sector and fixed number of particles
  // 
  // subsytemSizeX = number of sites along the x direction that belong to the subsytem
  // subsytemSizeY = number of sites along the y direction that belong to the subsytem
  // subsytemStartX = x momentum marking the beginning of the rectangluar subsystem
  // subsytemStartY = y momentum marking the beginning of the rectangluar subsystem
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = Kx momentum sector in which the entanglement matrix has to be evaluated 
  // kySector = Ky momentum sector in which the entanglement matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem
  virtual HermitianMatrix EvaluatePartialDensityMatrixMomentumSpace (int subsytemSizeX, int subsytemSizeY, int subsytemStartX, int subsytemStartY, int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState);

  // core part of the evaluation density matrix momentum space partition calculation
  // 
  // minIndex = first index to consider in complementary Hilbert space
  // nbrIndex = number of indices to consider in complementary Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixMomentumSpaceCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace, ComplexVector& groundState,  HermitianMatrix* densityMatrix);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given momentum sector and fixed number of particles
  // 
  // subsytemSizeX = number of sites along the x direction that belong to the subsytem
  // subsytemSizeY = number of sites along the y direction that belong to the subsytem
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = Kx momentum sector in which the entanglement matrix has to be evaluated 
  // kySector = Ky momentum sector in which the entanglement matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = entanglement matrix of the subsytem
  virtual ComplexMatrix EvaluatePartialEntanglementMatrixMomentumSpace (int subsytemSizeX, int subsytemSizeY, int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = Kx momentum sector in which the entanglement matrix has to be evaluated 
  // kySector = Ky momentum sector in which the entanglement matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, bool removeBinomialCoefficient = false);

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

  // find state index from an array
  //
  // stateDescription = array describing the state (stored as kx1,ky1,kx2,ky2,...)
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndexFromArray(int* stateDescription);

  // symmetrized a product of two uncoupled states 
  //
  // outputVector = reference on the vector which will contain the symmetrozed state
  // leftVector = reference on the vector associated to the first color
  // rightVector = reference on the vector associated to the second color
  // leftSpace = pointer to the Hilbert space of the first color
  // rightSpace = pointer to the Hilbert space of the second color
  // unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
  // return value = symmetrized state

  ComplexVector SymmetrizeU1U1State (ComplexVector& leftVector, ComplexVector& rightVector, BosonOnSquareLatticeMomentumSpaceLong* leftSpace, BosonOnSquareLatticeMomentumSpaceLong* rightSpace, bool unnormalizedBasisFlag, AbstractArchitecture* architecture);*/


 protected:

   // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentJ = current value of j for a single particle
  // currentJz = current value of jz for a single particle
  // currentKz = current value of kz for a single particle
  // currentTotalJz = current total value of Jz
  // currentTotalKz = current total value of Kz
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int currentJ, int currentJz,  int currentKz, int currentTotalJz, int currentTotalKz);
  
   // generate all states corresponding to the constraints
  // 
  // stateDescription = array that gives each state description
  // nbrBosons = number of bosons
  // currentJ = current value of j for a single particle
  // currentJz = current value of jz for a single particle
  // currentKz = current value of kz for a single particle
  // currentTotalJz = current total value of Jz
  // currentTotalKz = current total value of Kz
  // currentFermionicPosition = current fermionic position within the state description
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(ULONGLONG* stateDescription, int nbrBosons, int currentJ, int currentJz, int currentKz, int currentTotalJz, int currentTotalKz, int currentFermionicPosition, long pos);

//   // core part of the evaluation density matrix particle partition calculation
//   // 
//   // minIndex = first index to consider in complementary Hilbert space
//   // nbrIndex = number of indices to consider in complementary Hilbert space
//   // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
//   // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
//   // groundState = reference on the total system ground state
//   // densityMatrix = reference on the density matrix where result has to stored
//   // return value = number of components that have been added to the density matrix
//   virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
// 								  ComplexVector& groundState,  HermitianMatrix* densityMatrix);
// 
//   // core part of the evaluation density matrix particle partition calculation involving a sum of projetors 
//   // 
//   // minIndex = first index to consider in source Hilbert space
//   // nbrIndex = number of indices to consider in source Hilbert space
//   // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
//   // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
//   // nbrGroundStates = number of projectors
//   // groundStates = array of degenerate groundstates associated to each projector
//   // weights = array of weights in front of each projector
//   // densityMatrix = reference on the density matrix where result has to stored
//   // return value = number of components that have been added to the density matrix
//   virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
// 								  int nbrGroundStates, ComplexVector* groundStates, double* weights, HermitianMatrix* densityMatrix);

// get the quantum numbers j, jz, kz of a one particle state 
//
//quantumNumberJ = array that gives the quantum number j for a single particle state
//quantumNumberJz = array that gives the quantum number j for a single particle stateDescription
//quantumNumberKz = array that gives the quantum number j for a single particle state
inline void GetQuantumNumbersFromLinearizedIndex(int* quantumNumberJ, int* quantumNumberJz, int* quantumNumberKz)
  {
   for (int j = 0; j <= this->NbrFluxQuanta; ++j)
     {
      for (int jz = 0; jz <= j; ++jz)
      {
	  for (int kz = 0; kz <= this->NbrFluxQuanta - j; ++kz)
	  {
	    int index = -((j-1)*j*(2*j-1))/6 + this->NbrFluxQuanta*j*(j-1)/2 + (this->NbrFluxQuanta+1)*j+ (this->NbrFluxQuanta + 1 - j)*jz +kz;
	    quantumNumberJ[index] = j;
	    quantumNumberJz[index] = jz;
	    quantumNumberKz[index] = kz;
	    
	  }
	}
      }
  }

};


#endif


