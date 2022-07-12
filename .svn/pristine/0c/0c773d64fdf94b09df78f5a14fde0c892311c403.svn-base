////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of pair hopping p=1 (a.k.a. PXP model)             //
//            Hilbert space written as spin 1 chain with translations         //
//                                                                            //
//                        last modification : 13/03/2019                      //
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


#ifndef PAIRHOPPINGP1SPIN1CHAINWITHTRANSLATIONS_H
#define PAIRHOPPINGP1SPIN1CHAINWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"

#include <iostream>


using std::ostream;


class HermitianMatrix;
class RealMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


class PairHoppingP1AsSpin1ChainWithTranslations : public Spin1ChainWithTranslations
{

 protected:

 public:

  // default constructor
  //
  PairHoppingP1AsSpin1ChainWithTranslations ();

  // constructor for complete Hilbert space
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // memory = amount of memory granted for precalculations
  PairHoppingP1AsSpin1ChainWithTranslations (int chainLength, int momentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  PairHoppingP1AsSpin1ChainWithTranslations (const PairHoppingP1AsSpin1ChainWithTranslations& chain);

  // destructor
  //
  ~PairHoppingP1AsSpin1ChainWithTranslations ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  PairHoppingP1AsSpin1ChainWithTranslations& operator = (const PairHoppingP1AsSpin1ChainWithTranslations& chain);

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
  
  // apply the Z_0 operator (i.e the projector on |0><0|) on a given site
  //
  // globalSiteIndex = global index of the site
  // state = index of the state on which the operator has to be applied
  // return value = numerical coefficient once the projection is applied
  virtual double Z0 (int globalSiteIndex, int state);

  // apply the Z_+ operator (i.e the projector on |+><+|) on a given site
  //
  // globalSiteIndex = global index of the site
  // state = index of the state on which the operator has to be applied
  // return value = numerical coefficient once the projection is applied
  virtual double ZPlus (int globalSiteIndex, int state);

  // apply the Z_- operator (i.e the projector on |-><-|) on a given site
  //
  // globalSiteIndex = global index of the site
  // state = index of the state on which the operator has to be applied
  // return value = numerical coefficient once the projection is applied
  virtual double ZMinus (int globalSiteIndex, int state);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated (unused here)
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);
  
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

  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
  // szSymmetrySign = reference on the additional sign coming from the Sz<->-Sz symmetry
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslations, double& szSymmetrySign);

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();
  
};

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
// szSymmetrySign = reference on the additional sign coming from the Sz<->-Sz symmetry
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline unsigned long PairHoppingP1AsSpin1ChainWithTranslations::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslations, double& szSymmetrySign)
{
  return this->Spin1ChainWithTranslations::FindCanonicalForm(stateDescription, nbrTranslations);
}

// apply the Z_0 operator (i.e the projector on |0><0|) on a given site
//
// globalSiteIndex = global index of the site
// state = index of the state on which the operator has to be applied
// return value = numerical coefficient once the projection is applied

inline double PairHoppingP1AsSpin1ChainWithTranslations::Z0 (int globalSiteIndex, int state)
{
  if (((this->StateDescription[state] >> (globalSiteIndex * 2)) & 0x3ul) == 0x2ul)
    {
      return 1.0;
    }
  else
    {
      return 0.0;
    }
}

// apply the Z_+ operator (i.e the projector on |+><+|) on a given site
//
// globalSiteIndex = global index of the site
// state = index of the state on which the operator has to be applied
// return value = numerical coefficient once the projection is applied

inline double PairHoppingP1AsSpin1ChainWithTranslations::ZPlus (int globalSiteIndex, int state)
{
  if (((this->StateDescription[state] >> (globalSiteIndex * 2)) & 0x3ul) == 0x3ul)
    {
      return 1.0;
    }
  else
    {
      return 0.0;
    }
}

// apply the Z_- operator (i.e the projector on |-><-|) on a given site
//
// globalSiteIndex = global index of the site
// state = index of the state on which the operator has to be applied
// return value = numerical coefficient once the projection is applied

inline double PairHoppingP1AsSpin1ChainWithTranslations::ZMinus (int globalSiteIndex, int state)
{
  if (((this->StateDescription[state] >> (globalSiteIndex * 2)) & 0x3ul) == 0x0ul)
    {
      return 1.0;
    }
  else
    {
      return 0.0;
    }
}

#endif


