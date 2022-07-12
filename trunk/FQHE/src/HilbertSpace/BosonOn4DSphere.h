////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Cecile Repellin                 //
//                                                                            //
//                                                                            //
//                          class of bosons on 4D sphere                      //
//                                                                            //
//                                                                            //
//                        last modification : 19/10/2012                      //
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


#ifndef BOSONON4DSPHERE_H
#define BOSONON4DSPHERE_H

#include "config.h"
#include "HilbertSpace/BosonOnSphereShort.h"

#include <iostream>



class BosonOn4DSphere : public BosonOnSphereShort
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
  BosonOn4DSphere ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrFluxQuanta = number of flux quanta (p)
  // totalJz = total value of jz
  // totalKz = total value of kz
  // memory = amount of memory granted for precalculations
  BosonOn4DSphere (int nbrBosons, int nbrFluxQuanta, int totalJz, int totalKz, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOn4DSphere(const BosonOn4DSphere& bosons);

  // destructor
  //
  ~BosonOn4DSphere ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOn4DSphere& operator = (const BosonOn4DSphere& bosons);

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
  
    // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given (Jz, Kz) sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // jzSector = Jz sector in which the density matrix has to be evaluated 
  // kzSector = Kz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition(int nbrBosonSector, int jzSector, int kzsector,  RealVector& groundState, AbstractArchitecture* architecture = 0);

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
  virtual long GenerateStates(unsigned long* stateDescription, int nbrBosons, int currentJ, int currentJz, int currentKz, int currentTotalJz, int currentTotalKz, int currentFermionicPosition, long pos);

 
  
  
};


#endif


