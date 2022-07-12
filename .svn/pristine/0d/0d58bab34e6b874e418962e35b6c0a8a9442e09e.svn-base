////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                        Class author: Cecile Repellin                       //
//                                                                            //
//                                                                            //
//                     class of SU(3) (3 color) spin chain                    //
//                                                                            //
//                        last modification : 13/11/2017                      //
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


#ifndef SU3SPINCHAIN_H
#define SU3SPINCHAIN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Vector/RealVector.h"

#include <iostream>


using std::ostream;


class HermitianMatrix;
class RealMatrix;
class ComplexMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


class SU3SpinChain : public AbstractSpinChain
{

  friend class SU3SpinChainAnd2DTranslation;
//   friend class SU3SpinChainWithTranslationsAndSzSymmetry;
//   friend class SU3SpinChainWithTranslationsAndInversionSymmetry;
//   friend class SU3SpinChainWithSzSymmetry;

 protected:

  // flag to indicate if the total Tz and Y components are fixed
  bool FixedSpinProjectionFlag;
  //  total tz component (if fixed)
  int Tz;
  //  total y component (if fixed)
  int Y;
  
  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given maxMomentum sector
  int LookUpTableMemorySize;
  // shift used in each maxMomentum sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used maxMomentum value of the state an the second 
  int** LookUpTable;

  unsigned long* StateDescription;

 public:

  // default constructor
  //
  SU3SpinChain ();

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // tz = twice the value of total tz component
  // y = value of the total y component
  // memorySize = memory size in bytes allowed for look-up table
  SU3SpinChain (int chainLength, int tz, int y, int memorySize) ;

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  SU3SpinChain (const SU3SpinChain& chain);

  // destructor
  //
  virtual ~SU3SpinChain ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  SU3SpinChain& operator = (const SU3SpinChain& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // get the value of the spin (i.e. S) at a given site
  // 
  // site = site index
  // return value = twice the spin
  virtual int GetLocalSpin(int site);

  
  // return value of spin projection on (Oz) for a given state
  //
  // index = index of the state to test
  // return value = spin projection on (Oz)
  virtual int TotalTz (int index);

 
  // return index of resulting state from application of two-site exchange operator on a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = index of resulting state
  virtual int Pij (int i, int j, int state);
  
  // return index of resulting state from application of four-site exchange operator on a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = index of resulting state
  virtual int Pijkl (int i, int j, int k, int l, int state);
  
  // return eigenvalue of Sz_i associated to a given state
  //
  // i = first position
  // state = index of the state to consider
  // return value = corresponding eigenvalue
  virtual double Szi (int i, int state);


  // translate a state assuming the system have periodic boundary conditions (increasing the site index)
  //
  // nbrTranslations = number of translations to apply
  // state = index of the state to translate 
  // return value = index of resulting state
  virtual int TranslateState (int nbrTranslations, int state);

  // find state index
  //
  // state = state description
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long state);

  // find state index
  //
  // stateDescription = state description
  // maxBitPosition = maximum bit set to one in stateDescription
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int maxBitPosition);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

 protected:
  

  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = state that has been produced with the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state  
  virtual int SymmetrizeResult(unsigned long state, double& coefficient);

  // evaluate Hilbert space dimension
  //
  // sz = twice the Sz value
  // nbrSites = number of sites
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int tz, int y, int nbrSites);

  // generate all states corresponding to a given total Sz and no discrete symmetry constraint
  //
  // statePosition = position for the new states
  // sitePosition = site on chain where spin has to be changed
  // currentSz = total Sz value of current state
  // return value = number of generated states
  virtual long RawGenerateStates(long statePosition, int sitePosition, int currentTz, int currentY); 

  // generate all states with constraints 
  //
  // return value = number of generated states
  virtual long GenerateStates();

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);
  
};

// get the value of the spin (i.e. S) at a given site
// 
// site = site index
// return value = twice the spin

inline int SU3SpinChain::GetLocalSpin(int site)
{
  return 2;
}

// find state index
//
// stateDescription = state description
// return value = corresponding index

inline int SU3SpinChain::FindStateIndex(unsigned long stateDescription)
{
  int TmpMaxBitPosition = this->ChainLength << 1;
  while ((TmpMaxBitPosition > 0) && ((stateDescription >> TmpMaxBitPosition) == 0x0ul))
    {
      --TmpMaxBitPosition;
    }
  return this->FindStateIndex(stateDescription, TmpMaxBitPosition);
}

// factorized code that is used to symmetrize the result of any operator action
//
// state = state that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state  

inline int SU3SpinChain::SymmetrizeResult(unsigned long state, double& coefficient)
{
  return this->FindStateIndex(state);
}

#endif


