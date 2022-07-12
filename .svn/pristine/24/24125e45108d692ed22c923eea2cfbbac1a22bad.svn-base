////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of pair hopping p=1 (a.k.a. PXP model)             //
//                     Hilbert space written as spin 1 chain                  //
//                                                                            //
//                        last modification : 10/03/2019                      //
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


#ifndef PAIRHOPPINGP1SPIN1CHAIN_H
#define PAIRHOPPINGP1SPIN1CHAIN_H


#include "config.h"
#include "HilbertSpace/Spin1Chain.h"

#include <iostream>


using std::ostream;


class HermitianMatrix;
class RealMatrix;
class ComplexMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


class PairHoppingP1AsSpin1Chain : public Spin1Chain 
{

  friend class Spin1ChainWithTranslations;
  friend class Spin1ChainWithTranslationsAndSzSymmetry;
  friend class Spin1ChainWithTranslationsAndInversionSymmetry;
  friend class Spin1ChainWithSzSymmetry;
  friend class PairHoppingP1AsSpin1ChainWithTranslations;
  friend class PairHoppingP2AsSpin1ChainWithTranslations;
  
 protected:

  // true if the system uses periodic boundary conditions
  bool PeriodicBoundaryConditions;
  
  // true if the hilbert space has to be generated for the entanglement matrix calculations
  bool UseForEntanglementMatrix;
  
public:

  // default constructor
  //
  PairHoppingP1AsSpin1Chain ();

  // constructor for complete Hilbert space
  //
  // chainLength = number of spin / group of 2p+1 orbitals
  // periodicBoundaryConditions = true if the system uses periodic boundary conditions
  // memorySize = memory size in bytes allowed for look-up table
  // useForEntanglementMatrix = true if the hilbert space has to be generated for the entanglement matrix calculations
  PairHoppingP1AsSpin1Chain (int chainLength, bool periodicBoundaryConditions = false, int memorySize = -1, bool useForEntanglementMatrix = false);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  PairHoppingP1AsSpin1Chain (const PairHoppingP1AsSpin1Chain& chain);

  // destructor
  //
  virtual ~PairHoppingP1AsSpin1Chain ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  PairHoppingP1AsSpin1Chain& operator = (const PairHoppingP1AsSpin1Chain& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // apply the swap operator within the unit cell
  //
  // unitCellCoordinate = coordinate of the unit cell
  // siteIndex = index of the leftmost site within the unit cell
  // state = index of the state on which the operator has to be applied
  // return value = index of resulting state 
  virtual int SwapOperator (int unitCellCoordinate, int siteIndex, int state);

  // apply the swap operator within the unit cell with a constraint of the unit cell parity
  //
  // unitCellCoordinate = coordinate of the unit cell
  // siteIndex = index of the leftmost site within the unit cell
  // state = index of the state on which the operator has to be applied
  // return value = index of resulting state 
  virtual int SwapOperatorPlus (int unitCellCoordinate, int siteIndex, int state);

  // apply the operator coupling neighboring unit cells
  //
  // leftUnitCellCoordinate = coordinate of the left unit cell
  // rightUnitCellCoordinate = coordinate of the right unit cell
  // state = index of the state on which the operator has to be applied
  // return value = index of resulting state 
  virtual int PlusMinusOperator (int leftUnitCellCoordinate, int rightUnitCellCoordinate, int state);  
  
  // apply the operator coupling neighboring unit cells with a constraint of the unit cell parity
  //
  // leftUnitCellCoordinate = coordinate of the left unit cell
  // rightUnitCellCoordinate = coordinate of the right unit cell
  // state = index of the state on which the operator has to be applied
  // return value = index of resulting state 
  virtual int PlusMinusOperatorPlus (int leftUnitCellCoordinate, int rightUnitCellCoordinate, int state);  
  
 protected:

  // evaluate Hilbert space dimension
  //
  // previousSpin = value of the previous spin (0 for -1, 1 for 0 and 2 for +1)
  // sitePosition = site on chain where spin has to be changed
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int previousSpin, int sitePosition);

  // generate all states with neither constraint from boundary conditions nor discrete symmtry constraint
  //
  // statePosition = position for the new states
  // previousSpin = value of the previous spin (0 for -1, 1 for 0 and 2 for +1)
  // sitePosition = site on chain where spin has to be changed
  // return value = number of generated states
  virtual long RawGenerateStates(long statePosition, int previousSpin, int sitePosition);

  // generate all states with constraints 
  //
  // return value = number of generated states
  virtual long GenerateStates();
  
};

#endif


