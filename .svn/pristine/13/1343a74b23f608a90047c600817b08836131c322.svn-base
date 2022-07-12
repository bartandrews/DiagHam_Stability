////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                             Class author Cecile Repellin                   //
//                                                                            //
//                          class of bosons on CP2                            //
//                                                                            //
//                                                                            //
//                        last modification : 08/01/2013                      //
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


#ifndef BOSONONCP2_H
#define BOSONONCP2_H

#include "config.h"
#include "HilbertSpace/BosonOnSphereShort.h"

#include <iostream>



class BosonOnCP2 : public BosonOnSphereShort
{
  
 protected:

  // Number of flux quanta
  int NbrFluxQuanta;
  // total value of Tz
  int TotalTz;
  // total value of Y
  int TotalY;
  // total value of r
  int TotalR;
  // total value of s
  int TotalS;
  // array that gives the value of tz for one particle corresponding to the linearized index
  int* quantumNumberTz;
  // array that gives the value of y for one particle corresponding to the linearized index
  int* quantumNumberY;
  // array that gives the value of r for one particle corresponding to the linearized index
  int* quantumNumberR;
  // array that gives the value of s for one particle corresponding to the linearized index
  int* quantumNumberS;
  

 public:

  // default constructor
  // 
  BosonOnCP2 ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrFluxQuanta = number of flux quanta (p)
  // totalJz = total value of jz
  // totalKz = total value of kz
  // memory = amount of memory granted for precalculations
  BosonOnCP2 (int nbrBosons, int nbrFluxQuanta, int totalTz, int totalY, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnCP2(const BosonOnCP2& bosons);

  // destructor
  //
  ~BosonOnCP2 ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnCP2& operator = (const BosonOnCP2& bosons);

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
  RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition(int nbrBosonSector, int tzSector, int ysector,  RealVector& groundState, AbstractArchitecture* architecture = 0);

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
  
  // get the quantum numbers tz, y of a one particle state 
  //
  //quantumNumberTz = array that gives the quantum number tz for a single particle stateDescription
  //quantumNumberY = array that gives the quantum number y for a single particle state
  inline void GetQuantumNumbersFromLinearizedIndex(int* quantumNumberTz, int* quantumNumberY, int* quantumNumberR, int* quantumNumberS)
  {
    cout << this->NbrFluxQuanta << endl;
    for (int tzMax = 0; tzMax <= this->NbrFluxQuanta; ++tzMax)
    {
      for (int shiftedTz = 0; shiftedTz <= tzMax; ++shiftedTz)
	{
	  int tz = 2*shiftedTz - tzMax;
	  int y = 3*tzMax - 2*this->NbrFluxQuanta;
	  int index = this->GetLinearizedIndex(tz, y, 1);
	  quantumNumberTz[index] = tz;
	  quantumNumberY[index] = y;
	  quantumNumberR[index] = (y + 3*tz + 2*this->NbrFluxQuanta)/6;
	  quantumNumberS[index] = (y - 3*tz + 2*this->NbrFluxQuanta)/6;
// 	  cout << tz << " " << y << " " << quantumNumberR[index] << " " << quantumNumberS[index] << endl;
	}
    }
  }

  
  //compute the linearized index for the single particles quantum numbers (tz,y)
  //
  //tz = integer with the value of quantum number tz
  //y = integer with the value of quantum number y
  //nbrParticles = number of particles involved (1 for HilbertSpace, 2 for two-body interaction...)
  virtual int GetLinearizedIndex(int tz, int y, int nbrParticles);
  
 protected:

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentTz = current value of Tz for a single particle
  //currentTzMax = maximum value of Tz for a given value of y
  // currentY = current value of y for a single particle
  // currentTotalTz = current total value of tz
  // currentTotalY = current total value of y
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int currentTz, int currentTzMax,  int currentY, int currentTotalTz, int currentTotalY);

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
  virtual long GenerateStates(unsigned long* stateDescription, int nbrBosons, int currentTz, int currentTzMax, int currentY, int currentTotalTz, int currentTotalY, int currentFermionicPosition, long pos);

  // request whether state with given index satisfies a general Pauli exclusion principle
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  virtual bool HasPauliExclusions(int index, int pauliK, int pauliR);
  
  // convert a state such that its components are now expressed in the unnormalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has to be normalized to 1
  // symmetryFactor = if true also remove the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertToUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);
  
  
};


//compute the linearized index for the single particles quantum numbers (tz,y)
  //
  //tz = integer with the value of quantum number tz
  //y = integer with the value of quantum number y
  //nbrParticles = number of particles involved (1 for HilbertSpace, 2 for two-body interaction...)
inline int BosonOnCP2::GetLinearizedIndex(int tz, int y, int nbrParticles)
  {
    int tzMax = (y + 2*nbrParticles*this->NbrFluxQuanta)/3;
    int index = tzMax*(tzMax + 1)/2 + (tz + tzMax)/2 ;
//     cout <<index << endl;
//     int index = this->LzMax - (tzMax*(tzMax + 1)/2 + (tz + tzMax)/2) ;
    return index;
  }

#endif


