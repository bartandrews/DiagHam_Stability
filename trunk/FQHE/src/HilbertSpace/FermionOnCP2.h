////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                     class author : Cecile Repellin                         //
//                                                                            //
//                        class of fermions the CP2                           //
//                                                                            //
//                                                                            //
//                        last modification : 25/02/2013                      //
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


#ifndef FERMIONONCP2_H
#define FERMIONONCP2_H

#include "config.h"
#include "HilbertSpace/FermionOnSphere.h"

#include <iostream>



class FermionOnCP2 : public FermionOnSphere
{

 protected:
  // Number of flux quanta
  int NbrFluxQuanta;
  // Minimum value of Y
  int MinY;
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
  FermionOnCP2 ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrFluxQuanta = number of flux quanta (p)
  // totalJz = total value of jz
  // totalKz = total value of kz
  // memory = amount of memory granted for precalculations
  FermionOnCP2 (int nbrFermions, int nbrFluxQuanta, int totalTz, int totalY, unsigned long memory = 10000000);
  
  // constructor for truncated Hilbert space
  // 
  // nbrBosons = number of bosons
  // nbrFluxQuanta = number of flux quanta (p) of full Hilbert space
  // minY = minimum accessible value of y
  // totalTz = total value of Tz
  // totalY = total value of y
  // memory = amount of memory granted for precalculations
  FermionOnCP2 (int nbrFermions, int nbrFluxQuanta, int minY, int totalTz, int totalY, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  FermionOnCP2(const FermionOnCP2& bosons);

  // destructor
  //
  ~FermionOnCP2 ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnCP2& operator = (const FermionOnCP2& bosons);

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
  // nbrFermionSector = number of particles that belong to the subsytem 
  // jzSector = Jz sector in which the density matrix has to be evaluated 
  // kzSector = Kz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition(int nbrFermionSector, int tzSector, int ysector,  RealVector& groundState, AbstractArchitecture* architecture = 0);
  
  
  // find state index from an array
  //
  // stateDescription = array describing the state (stored as kx1,ky1,kx2,ky2,...)
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndexFromArray(int* stateDescription);
  
  
  // get the quantum numbers tz, y of a one particle state 
  //
  //quantumNumberTz = array that gives the quantum number tz for a single particle stateDescription
  //quantumNumberY = array that gives the quantum number y for a single particle state
  void GetQuantumNumbersFromLinearizedIndex(int* quantumNumberTz, int* quantumNumberY, int* quantumNumberR, int* quantumNumberS);
  
  //compute the linearized index for the single particles quantum numbers (tz,y)
  //
  //tz = integer with the value of quantum number tz
  //y = integer with the value of quantum number y
  //nbrParticles = number of particles involved (1 for HilbertSpace, 2 for two-body interaction...)
  //return value = index of the state
  int GetLinearizedIndex(int tz, int y, int nbrParticles);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given tz, y sector.
  // 
  // nbrFermionSector = number of particles that belong to the subsytem 
  // tzSector = tz sector in which the density matrix has to be evaluated
  // ySector = y sector in which the density matrix has to be evaluated
  // groundState = reference on the total system ground state
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrixParticlePartition (int nbrFermionSector, int tzSector, int ySector, RealVector& groundState, bool removeBinomialCoefficient = false);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
  // The entanglement matrix is only evaluated in a given tz, Y sector and both A and B are resticted to a given number of orbitals
  //  
  // nbrFermionSector = number of particles that belong to the subsytem 
  // tzSector = Tz sector in which the density matrix has to be evaluated 
  // ySector = Y sector in which the density matrix has to be evaluated 
  // maxYA = maximum value of Y in the A partition
  // minYB = minimum value of Y in the B partition
  // groundState = reference on the total system ground state
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrixParticlePartition(int nbrFermionSector, int tzSector, int ySector,
									       int maxYA, int minYB, 
									       RealVector& groundState, bool removeBinomialCoefficient);
  
  // evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
  // The entanglement matrix is only evaluated in a given tz, y sector and computed from precalculated particle entanglement matrix. The cut has to be at a given value of Y
  // 
  // nbrFermionSector = number of particles that belong to the subsytem 
  // tzSector = Tz sector in which the density matrix has to be evaluated 
  // ySector = Y sector in which the density matrix has to be evaluated 
  // maxYA = maximum value of Y in the A partition
  // weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
  // minYB = minimum value of Y that has to be kept in the B partition
  // weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
  // entanglementMatrix = reference on the entanglement matrix (will be overwritten)
  // return value = reference on the entanglement matrix
  virtual RealMatrix& EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int tzSector, int ySector, int maxYA, double* weightOrbitalA, 
														int minYB, double* weightOrbitalB, RealMatrix& entanglementMatrix);
  
  
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
  // nbrBosons = number of bosons
  // currentJ = current value of j for a single particle
  // currentJz = current value of jz for a single particle
  // currentKz = current value of kz for a single particle
  // currentTotalJz = current total value of Jz
  // currentTotalKz = current total value of Kz
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int currentTz, int currentTzMax, int currentY, int currentTotalTz, int currentTotalY, long pos);
  
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

 // get the quantum numbers tz, y of a one particle state 
  //
  //quantumNumberTz = array that gives the quantum number tz for a single particle stateDescription
  //quantumNumberY = array that gives the quantum number y for a single particle state
  inline void FermionOnCP2::GetQuantumNumbersFromLinearizedIndex(int* quantumNumberTz, int* quantumNumberY, int* quantumNumberR, int* quantumNumberS)
  {
    for (int tzMax = (2*this->NbrFluxQuanta + this->MinY) / 3; tzMax <= this->NbrFluxQuanta; ++tzMax)
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
  //return value = index of the state
  inline int FermionOnCP2::GetLinearizedIndex(int tz, int y, int nbrParticles)
  {
    int TruncationShift = (2 * this->NbrFluxQuanta + this->MinY) * (2*this->NbrFluxQuanta + this->MinY + 3) / 18;
    int tzMax = (y + 2*nbrParticles*this->NbrFluxQuanta)/3;
    int index = tzMax*(tzMax + 1)/2 + (tz + tzMax)/2  - TruncationShift;
//     cout <<index << endl;
//     int index = this->LzMax - (tzMax*(tzMax + 1)/2 + (tz + tzMax)/2) ;
    return index;
  }

#endif
