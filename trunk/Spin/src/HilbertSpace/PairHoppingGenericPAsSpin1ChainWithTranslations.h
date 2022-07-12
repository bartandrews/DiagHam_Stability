////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of pair hopping with generic p                  //
//            Hilbert space written as spin 1 chain with translations         //
//                                                                            //
//                        last modification : 16/10/2019                      //
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


#ifndef PAIRHOPPINGGENERICPSPIN1CHAINWITHTRANSLATIONS_H
#define PAIRHOPPINGGENERICPSPIN1CHAINWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainWithTranslations.h"

#include <iostream>


using std::ostream;


class HermitianMatrix;
class RealMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


class PairHoppingGenericPAsSpin1ChainWithTranslations : public PairHoppingP2AsSpin1ChainWithTranslations
{

 protected:

  // number of spin 1 per unit cell
  int UnitCellSize;

  // mask to isolate the first unit cell
  unsigned long FirstUnitCellMask;

  // number of possible unit cells satisfying the pair hopping constraints
  int NbrPossibleUnitCells;
  // possible unit cells satisfying the pair hopping constraints
  unsigned long* PossibleUnitCells;
  // number of plusses for each possible unit cells satisfying the pair hopping constraints
  int* PossibleUnitCellsNbrPlus;
  // number of plusses for each possible unit cells satisfying the pair hopping constraints
  int* PossibleUnitCellsNbrMinus;
  
  // number of possible unit cells satisfying the pair hopping constraints for each possible number of minuses
  int* NbrPossibleUnitCellsPerNbrMinus;
  // possible unit cells satisfying the pair hopping constraints for each possible number of minuses
  unsigned long** PossibleUnitCellsPerNbrMinus;
  // number of plusses for each possible unit cells satisfying the pair hopping constraints for each possible number of minuses
  int** PossibleUnitCellsNbrPlusPerNbrMinus;
  
 public:

  // default constructor
  //
  PairHoppingGenericPAsSpin1ChainWithTranslations ();

  // constructor for complete Hilbert space
  //
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // pValue = p value
  // memory = amount of memory granted for precalculations
  PairHoppingGenericPAsSpin1ChainWithTranslations (int chainLength, int momentum, int pValue, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  PairHoppingGenericPAsSpin1ChainWithTranslations (const PairHoppingGenericPAsSpin1ChainWithTranslations& chain);

  // destructor
  //
  ~PairHoppingGenericPAsSpin1ChainWithTranslations ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  PairHoppingGenericPAsSpin1ChainWithTranslations& operator = (const PairHoppingGenericPAsSpin1ChainWithTranslations& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // apply the swap operator within the unit cell
  //
  // unitCellCoordinate = coordinate of the unit cell
  // siteIndex = index of the leftmost site within the unit cell
  // state = index of the state on which the operator has to be applied
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state 
  virtual int SwapOperator (int unitCellCoordinate, int siteIndex, int state, double& coefficient, int& nbrTranslation);

  // apply the operator coupling neighboring unit cells
  //
  // leftUnitCellCoordinate = coordinate of the left unit cell
  // rightUnitCellCoordinate = coordinate of the right unit cell
  // state = index of the state on which the operator has to be applied
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state 
  virtual int PlusMinusOperator (int leftUnitCellCoordinate, int rightUnitCellCoordinate, int state, double& coefficient, int& nbrTranslation);  
  
  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated  (fixing the boundray unit cells to nbr left minuses= szSector / (pValue + 1), nbr right pluses = szSector % (pValue + 1))
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture = 0);
    
  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated  (fixing the boundray unit cells to nbr left minuses= szSector / (pValue + 1), nbr right pluses = szSector % (pValue + 1))
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);
    
 protected:

  // evaluate Hilbert space dimension
  //
  // previousSpin = value of the previous spin (0 for -1, 1 for 0 and 2 for +1)
  // initialNbrMinus = number of -1 spins in the first unit cell
  // sitePosition = site on chain where spin has to be changed
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int previousSpin, int initialNbrMinus, int sitePosition);

  // generate all states with neither constraint from boundary conditions nor discrete symmtry constraint
  //
  // statePosition = position for the new states
  // previousSpin = value of the previous spin (0 for -1, 1 for 0 and 2 for +1)
  // initialNbrMinus = number of -1 spins in the first unit cell
  // sitePosition = site on chain where spin has to be changed
  // return value = number of generated states
  virtual long RawGenerateStates(long statePosition, int previousSpin, int initialNbrMinus, int sitePosition);

  // generate all states with constraints 
  //
  // return value = number of generated states
  virtual long GenerateStates();

  // generate all the posible unit cells 
  //
  virtual void GenerateAllPossibleUnitCells();

};


#endif


