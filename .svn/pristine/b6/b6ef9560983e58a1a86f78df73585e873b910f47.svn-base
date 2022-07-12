////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of spin 1 chain with the Sz<->-Sz symmetry             //
//                                                                            //
//                        last modification : 11/05/2017                      //
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


#ifndef SPIN1CHAINWITHSZSYMMETRY_H
#define SPIN1CHAINWITHSZSYMMETRY_H


#include "config.h"
#include "HilbertSpace/Spin1Chain.h"
#include "Vector/RealVector.h"

#include <iostream>


using std::ostream;
using std::hex;
using std::dec;


class HermitianMatrix;
class RealMatrix;
class ComplexMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


//  precalculation table used to apply the Sz<->-Sz symmetry
//
static  unsigned long Spin1ChainWithSzInversionTable[] = {
0xfful, 0xfdul, 0xfeul, 0xfcul, 0xf7ul, 0xf5ul, 0xf6ul, 0xf4ul, 0xfbul, 0xf9ul, 0xfaul, 0xf8ul, 0xf3ul, 0xf1ul, 0xf2ul, 0xf0ul, 
0xdful, 0xddul, 0xdeul, 0xdcul, 0xd7ul, 0xd5ul, 0xd6ul, 0xd4ul, 0xdbul, 0xd9ul, 0xdaul, 0xd8ul, 0xd3ul, 0xd1ul, 0xd2ul, 0xd0ul, 
0xeful, 0xedul, 0xeeul, 0xecul, 0xe7ul, 0xe5ul, 0xe6ul, 0xe4ul, 0xebul, 0xe9ul, 0xeaul, 0xe8ul, 0xe3ul, 0xe1ul, 0xe2ul, 0xe0ul, 
0xcful, 0xcdul, 0xceul, 0xccul, 0xc7ul, 0xc5ul, 0xc6ul, 0xc4ul, 0xcbul, 0xc9ul, 0xcaul, 0xc8ul, 0xc3ul, 0xc1ul, 0xc2ul, 0xc0ul, 
0x7ful, 0x7dul, 0x7eul, 0x7cul, 0x77ul, 0x75ul, 0x76ul, 0x74ul, 0x7bul, 0x79ul, 0x7aul, 0x78ul, 0x73ul, 0x71ul, 0x72ul, 0x70ul, 
0x5ful, 0x5dul, 0x5eul, 0x5cul, 0x57ul, 0x55ul, 0x56ul, 0x54ul, 0x5bul, 0x59ul, 0x5aul, 0x58ul, 0x53ul, 0x51ul, 0x52ul, 0x50ul, 
0x6ful, 0x6dul, 0x6eul, 0x6cul, 0x67ul, 0x65ul, 0x66ul, 0x64ul, 0x6bul, 0x69ul, 0x6aul, 0x68ul, 0x63ul, 0x61ul, 0x62ul, 0x60ul, 
0x4ful, 0x4dul, 0x4eul, 0x4cul, 0x47ul, 0x45ul, 0x46ul, 0x44ul, 0x4bul, 0x49ul, 0x4aul, 0x48ul, 0x43ul, 0x41ul, 0x42ul, 0x40ul, 
0xbful, 0xbdul, 0xbeul, 0xbcul, 0xb7ul, 0xb5ul, 0xb6ul, 0xb4ul, 0xbbul, 0xb9ul, 0xbaul, 0xb8ul, 0xb3ul, 0xb1ul, 0xb2ul, 0xb0ul, 
0x9ful, 0x9dul, 0x9eul, 0x9cul, 0x97ul, 0x95ul, 0x96ul, 0x94ul, 0x9bul, 0x99ul, 0x9aul, 0x98ul, 0x93ul, 0x91ul, 0x92ul, 0x90ul, 
0xaful, 0xadul, 0xaeul, 0xacul, 0xa7ul, 0xa5ul, 0xa6ul, 0xa4ul, 0xabul, 0xa9ul, 0xaaul, 0xa8ul, 0xa3ul, 0xa1ul, 0xa2ul, 0xa0ul, 
0x8ful, 0x8dul, 0x8eul, 0x8cul, 0x87ul, 0x85ul, 0x86ul, 0x84ul, 0x8bul, 0x89ul, 0x8aul, 0x88ul, 0x83ul, 0x81ul, 0x82ul, 0x80ul, 
0x3ful, 0x3dul, 0x3eul, 0x3cul, 0x37ul, 0x35ul, 0x36ul, 0x34ul, 0x3bul, 0x39ul, 0x3aul, 0x38ul, 0x33ul, 0x31ul, 0x32ul, 0x30ul, 
0x1ful, 0x1dul, 0x1eul, 0x1cul, 0x17ul, 0x15ul, 0x16ul, 0x14ul, 0x1bul, 0x19ul, 0x1aul, 0x18ul, 0x13ul, 0x11ul, 0x12ul, 0x10ul, 
0x2ful, 0x2dul, 0x2eul, 0x2cul, 0x27ul, 0x25ul, 0x26ul, 0x24ul, 0x2bul, 0x29ul, 0x2aul, 0x28ul, 0x23ul, 0x21ul, 0x22ul, 0x20ul, 
0x0ful, 0x0dul, 0x0eul, 0x0cul, 0x07ul, 0x05ul, 0x06ul, 0x04ul, 0x0bul, 0x09ul, 0x0aul, 0x08ul, 0x03ul, 0x01ul, 0x02ul, 0x00ul};



class Spin1ChainWithSzSymmetry : public Spin1Chain
{

  friend class Spin1ChainWithTranslations;
  friend class Spin1ChainWithTranslationsAndSzSymmetry;
  friend class Spin1ChainWithTranslationsAndInversionSymmetry;

 protected:

  // sign of the Sz<->-Sz symmetry sector
  double SzSymmetrySector;

  // mask needed for the Sz<->-Sz symmetry application
  unsigned long SzSymmetryMask;

  // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;
  // number of state in each orbit
  int* NbrStateInOrbit;

 public:

  // default constructor
  //
  Spin1ChainWithSzSymmetry ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin
  // szSymmetrySector = Sz<->-Sz symmetry sector (can be either +1 or -1)
  // memorySize = memory size in bytes allowed for look-up table
  Spin1ChainWithSzSymmetry (int chainLength, int szSymmetrySector, unsigned long memorySize);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // szSymmetrySector = Sz<->-Sz symmetry sector (can be either +1 or -1)
  // sz = twice the value of total Sz component (should be equal to zero)
  // memorySize = memory size in bytes allowed for look-up table
  Spin1ChainWithSzSymmetry (int chainLength, int szSymmetrySector, int sz, unsigned long memory) ;

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1ChainWithSzSymmetry (const Spin1ChainWithSzSymmetry& chain);

  // destructor
  //
  ~Spin1ChainWithSzSymmetry ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1ChainWithSzSymmetry& operator = (const Spin1ChainWithSzSymmetry& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // return index of resulting state from application of S+_i operator on a given state
  //
  // i = position of S+ operator
  // state = index of the state to be applied on S+_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int Spi (int i, int state, double& coefficient);

  // return index of resulting state from application of S-_i operator on a given state
  //
  // i = position of S- operator
  // state = index of the state to be applied on S-_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int Smi (int i, int state, double& coefficient);

  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSpj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S+_i S+_j operator on a given state
  //
  // i = position of first S+ operator
  // j = position of second S+ operator
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SpiSpj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S-_i S-_j operator on a given state
  //
  // i = position of first S- operator
  // j = position of second S- operator
  // state = index of the state to be applied on S-_i S-_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSmj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S+_i Sz_j operator on a given state
  //
  // i = position of S+ operator
  // j = position of Sz operator
  // state = index of the state to be applied on S+_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SpiSzj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S-_i Sz_j operator on a given state
  //
  // i = position of S- operator
  // j = position of Sz operator
  // state = index of the state to be applied on S-_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSzj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S-_i1 S+_j1 S-_i2 S+_j2 operator on a given state
  //
  // i1 = position of leftmost S- operator
  // j1 = position of leftmost S+ operator
  // i2 = position of rightmost S- operator
  // j2 = position of rightmost S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSpjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient);

  // return index of resulting state from application of Sz_i1 Sz_j1 S-_i2 S+_j2 operator on a given state
  //
  // i1 = position of first Sz operator
  // j1 = position of second Sz operator
  // i2 = position of S- operator
  // j2 = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SziSzjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture = 0);
	
  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

 protected:

  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = state that has been produced with the operator action
  // nbrStateInOrbit = original number of states in the orbit before the operator action
 // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state  
  virtual int SymmetrizeResult(unsigned long state, int nbrStateInOrbit, double& coefficient);

  // generate all states with constraints 
  //
  // return value = number of generated states
  virtual long GenerateStates();

  // find canonical form of a state description
  //
  // stateDescription = unsigned integer describing the state
  // szSymmetrySign = reference on the additional sign coming from the Sz<->-Sz symmetry
  // return value = canonical form of a state description 
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, double& szSymmetrySign);

  //  test if the state is compatible with the discrete symmetry sector
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the discrete symmetry sector
  virtual bool TestDiscreteSymmetryConstraint(unsigned long stateDescription);

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  virtual int FindOrbitSize(unsigned long stateDescription);

  // apply the Sz<->-Sz symmetry to a state description
  //
  // stateDescription = reference on the state on which the Sz<->-Sz symmetry has to be applied
  virtual void ApplySzSymmetry(unsigned long& stateDescription);

  // compute the rescaling factors
  //
  virtual void ComputeRescalingFactors();

};

// factorized code that is used to symmetrize the result of any operator action
//
// state = state that has been produced with the operator action
// nbrStateInOrbit = original number of states in the orbit before the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state  

inline int Spin1ChainWithSzSymmetry::SymmetrizeResult(unsigned long state, int nbrStateInOrbit, double& coefficient)
{
  double TmpSign;
  state = this->FindCanonicalForm(state, TmpSign);
  int TmpMaxMomentum = 2 * this->ChainLength;
  while (((state >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
    --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= TmpSign;
      coefficient *= this->RescalingFactors[nbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
     }
  return TmpIndex;
}

// apply the Sz<->-Sz symmetry to a state description
//
// stateDescription = reference on the state on which the Sz<->-Sz symmetry has to be applied

inline void Spin1ChainWithSzSymmetry::ApplySzSymmetry(unsigned long& stateDescription)
{
#ifdef __64_BITS__
  unsigned long TmpState = Spin1ChainWithSzInversionTable[stateDescription & 0xfful];
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 8) & 0xfful] << 8;
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 16) & 0xfful] << 16;
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 24) & 0xfful] << 24;
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 32) & 0xfful] << 32;
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 40) & 0xfful] << 40;
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 48) & 0xfful] << 48;
  TmpState |= Spin1ChainWithSzInversionTable[stateDescription >> 56] << 56; 
#else
  unsigned long TmpState = Spin1ChainWithSzInversionTable[stateDescription & 0xfful];
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 8) & 0xfful] << 8;
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 16) & 0xfful] << 16;
  TmpState |= Spin1ChainWithSzInversionTable[stateDescription >> 24] << 24;
#endif	
  TmpState &= this->SzSymmetryMask;
  stateDescription = TmpState;
}

// find canonical form of a state description 
//
// stateDescription = unsigned integer describing the state
// szSymmetrySign = reference on the additional sign coming from the Sz<->-Sz symmetry
// return value = canonical form of a state description

inline unsigned long Spin1ChainWithSzSymmetry::FindCanonicalForm(unsigned long stateDescription, double& szSymmetrySign)
{
  unsigned long CanonicalState = stateDescription;
  this->ApplySzSymmetry(stateDescription);  
  if (stateDescription < CanonicalState)
    {
      CanonicalState = stateDescription;
      szSymmetrySign = this->SzSymmetrySector;      
    }
  else
    {
      szSymmetrySign = 1.0;
    }
  return CanonicalState;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int Spin1ChainWithSzSymmetry::FindOrbitSize(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  this->ApplySzSymmetry(TmpStateDescription);
  if (stateDescription == TmpStateDescription)
    return 1;
  else
    return 2;
}

//  test if the state is compatible with the discrete symmetry sector
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the discrete symmetry sector

inline bool Spin1ChainWithSzSymmetry::TestDiscreteSymmetryConstraint(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  this->ApplySzSymmetry(TmpStateDescription);
  if ((stateDescription == TmpStateDescription) && (this->SzSymmetrySector < 0.0))
    return false;
  else
    return true;
}

#endif


