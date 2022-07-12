////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of spin 1/2 chain with a fixed Sz value                //
//                and a fixed parity under the mirror symmetry                //
//                                                                            //
//                        last modification : 29/06/2015                      //
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


#ifndef SPIN1_2CHAINMIRRORSYMMETRY_H
#define SPIN1_2CHAINMIRRORSYMMETRY_H


#include "config.h"
#include "HilbertSpace/Spin1_2ChainNew.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;
using std::hex;
using std::dec;


#ifdef __64_BITS__
#define SPIN1_2CHAIN_MIRRORSYMMETRIC_BIT  0x8000000000000000ul
#define SPIN1_2CHAIN_MIRRORSYMMETRIC_MASK 0x7ffffffffffffffful
#define SPIN1_2CHAIN_PARITYSYMMETRIC_BIT  0x4000000000000000ul
#define SPIN1_2CHAIN_PARITYSYMMETRIC_MASK 0xbffffffffffffffful
#define SPIN1_2CHAIN_FULLSYMMETRY_BIT 0xc000000000000000ul
#define SPIN1_2CHAIN_FULLSYMMETRY_MASK 0x3ffffffffffffffful
#else
#define SPIN1_2CHAIN_MIRRORSYMMETRIC_BIT  0x80000000ul
#define SPIN1_2CHAIN_MIRRORSYMMETRIC_MASK 0x7ffffffful
#define SPIN1_2CHAIN_PARITYSYMMETRIC_BIT  0x40000000ul
#define SPIN1_2CHAIN_PARITYSYMMETRIC_MASK 0xbffffffful
#define SPIN1_2CHAIN_FULLSYMMETRY_BIT 0xc0000000ul
#define SPIN1_2CHAIN_FULLSYMMETRY_MASK 0x3ffffffful
#endif


static unsigned long Spin1_2ChainMirrorSymmetryMirrorTable[] =  {0x0ul, 0x80ul, 0x40ul, 0xc0ul, 0x20ul, 0xa0ul, 0x60ul, 0xe0ul, 0x10ul, 0x90ul, 0x50ul, 
								 0xd0ul, 0x30ul, 0xb0ul, 0x70ul, 0xf0ul, 0x8ul, 0x88ul, 0x48ul, 0xc8ul, 0x28ul, 0xa8ul, 
								 0x68ul, 0xe8ul, 0x18ul, 0x98ul, 0x58ul, 0xd8ul, 0x38ul, 0xb8ul, 0x78ul, 0xf8ul, 0x4ul, 
								 0x84ul, 0x44ul, 0xc4ul, 0x24ul, 0xa4ul, 0x64ul, 0xe4ul, 0x14ul, 0x94ul, 0x54ul, 0xd4ul, 
								 0x34ul, 0xb4ul, 0x74ul, 0xf4ul, 0xcul, 0x8cul, 0x4cul, 0xccul, 0x2cul, 0xacul, 0x6cul, 
								 0xecul, 0x1cul, 0x9cul, 0x5cul, 0xdcul, 0x3cul, 0xbcul, 0x7cul, 0xfcul, 0x2ul, 0x82ul, 
								 0x42ul, 0xc2ul, 0x22ul, 0xa2ul, 0x62ul, 0xe2ul, 0x12ul, 0x92ul, 0x52ul, 0xd2ul, 0x32ul, 
								 0xb2ul, 0x72ul, 0xf2ul, 0xaul, 0x8aul, 0x4aul, 0xcaul, 0x2aul, 0xaaul, 0x6aul, 0xeaul, 
								 0x1aul, 0x9aul, 0x5aul, 0xdaul, 0x3aul, 0xbaul, 0x7aul, 0xfaul, 0x6ul, 0x86ul, 0x46ul, 
								 0xc6ul, 0x26ul, 0xa6ul, 0x66ul, 0xe6ul, 0x16ul, 0x96ul, 0x56ul, 0xd6ul, 0x36ul, 0xb6ul, 
								 0x76ul, 0xf6ul, 0xeul, 0x8eul, 0x4eul, 0xceul, 0x2eul, 0xaeul, 0x6eul, 0xeeul, 0x1eul, 
								 0x9eul, 0x5eul, 0xdeul, 0x3eul, 0xbeul, 0x7eul, 0xfeul, 0x1ul, 0x81ul, 0x41ul, 0xc1ul, 
								 0x21ul, 0xa1ul, 0x61ul, 0xe1ul, 0x11ul, 0x91ul, 0x51ul, 0xd1ul, 0x31ul, 0xb1ul, 0x71ul, 
								 0xf1ul, 0x9ul, 0x89ul, 0x49ul, 0xc9ul, 0x29ul, 0xa9ul, 0x69ul, 0xe9ul, 0x19ul, 0x99ul, 
								 0x59ul, 0xd9ul, 0x39ul, 0xb9ul, 0x79ul, 0xf9ul, 0x5ul, 0x85ul, 0x45ul, 0xc5ul, 0x25ul, 
								 0xa5ul, 0x65ul, 0xe5ul, 0x15ul, 0x95ul, 0x55ul, 0xd5ul, 0x35ul, 0xb5ul, 0x75ul, 0xf5ul, 
								 0xdul, 0x8dul, 0x4dul, 0xcdul, 0x2dul, 0xadul, 0x6dul, 0xedul, 0x1dul, 0x9dul, 0x5dul, 
								 0xddul, 0x3dul, 0xbdul, 0x7dul, 0xfdul, 0x3ul, 0x83ul, 0x43ul, 0xc3ul, 0x23ul, 0xa3ul, 
								 0x63ul, 0xe3ul, 0x13ul, 0x93ul, 0x53ul, 0xd3ul, 0x33ul, 0xb3ul, 0x73ul, 0xf3ul, 0xbul, 
								 0x8bul, 0x4bul, 0xcbul, 0x2bul, 0xabul, 0x6bul, 0xebul, 0x1bul, 0x9bul, 0x5bul, 0xdbul, 
								 0x3bul, 0xbbul, 0x7bul, 0xfbul, 0x7ul, 0x87ul, 0x47ul, 0xc7ul, 0x27ul, 0xa7ul, 0x67ul, 
								 0xe7ul, 0x17ul, 0x97ul, 0x57ul, 0xd7ul, 0x37ul, 0xb7ul, 0x77ul, 0xf7ul, 0xful, 0x8ful, 
								 0x4ful, 0xcful, 0x2ful, 0xaful, 0x6ful, 0xeful, 0x1ful, 0x9ful, 0x5ful, 0xdful, 0x3ful, 
								 0xbful, 0x7ful, 0xfful};


class Spin1_2ChainMirrorSymmetry : public Spin1_2ChainNew
{

 protected:

  // largest integer associated to a configuration belonging to the Hilbert space
  unsigned long  MaxStateDescription;
  // smallest integer associated to a configuration belonging to the Hilbert space
  unsigned long MinStateDescription;

  // shift to apply to a state before mirroring its expression
  int MirrorShift;
  // shift to apply to a state after mirroring its expression
  int MirrorUnshift;

  // parity under the mirror symmetry
  int MirrorParity;
  // additional sign due to the mirror symmetry 
  double MirrorParitySign;

  // a temporary variable to store the symmetry of a state
  unsigned long SOperatorSignature;

 public:


  // default constructor
  //
  Spin1_2ChainMirrorSymmetry ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin 1/2
  // sz = twice the value of total Sz component
  // parity = parity under the mirror symmetry
  // memorySize = memory size in bytes allowed for look-up table
  Spin1_2ChainMirrorSymmetry (int chainLength, int sz, int parity, int memorySize);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1_2ChainMirrorSymmetry (const Spin1_2ChainMirrorSymmetry& chain);

  // destructor
  //
  ~Spin1_2ChainMirrorSymmetry ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1_2ChainMirrorSymmetry& operator = (const Spin1_2ChainMirrorSymmetry& chain);

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

  // compute the parity (prod_i Sz_i) for a given state
  //
  // state = index of the state to be applied on Sz_i operator
  // return value = total Sz value
  virtual unsigned long GetParity (int state);

  // return index of resulting state from application of P_ij operator on a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to be applied on P_ij operator
  // return value = index of resulting state
  virtual int Pij (int i, int j, int state);

  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSpj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S+_i S-_j operator on a given state
  //
  // i = position of S+ operator
  // j = position of S- operator
  // state = index of the state to be applied on S+_i S-_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SpiSmj (int i, int j, int state, double& coefficient);

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

  // return index of resulting state from application of S+_i S-_j Sz_k operator on a given state
  //
  // i = position of S+ operator
  // j = position of S- operator
  // k = position of Sz operator
  // state = index of the state to be applied on S+_i S-_j Sz_k operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SpiSmjSzk (int i, int j, int k, int state, double& coefficient);

  // translate a state assuming the system have periodic boundary
  // conditions (increasing the site index)
  //
  // nbrTranslations = number of translations to apply
  // state = index of the state to translate 
  // return value = index of resulting state
  virtual int TranslateState (int nbrTranslations, int state);

 protected:

  // find state index
  //
  // stateDescription = state description
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription);

  // factorized code that is used to symmetrize result of any spin operator
  //
  // state = reference on the state that has been produced with any operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state
  virtual int SymmetrizeResult(unsigned long state, double& coefficient);      

  // get canonical expression of a given state and its symmetry
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state (with symmetry bit)
  virtual unsigned long GetSignedCanonicalState (unsigned long initialState);

  // get the mirror symmetric of a given state 
  //
  // initialState = state to which the mirror symmetry has to be applied
  // return value = flipped configuration
  virtual unsigned long ApplyMirrorSymmetry (unsigned long initialState);

  // get the normalization factor in front of each basis state (i.e. 1/sqrt(orbit size))
  //
  // return value = pointer to normalization factors
  virtual double* GetBasisNormalization();
 
};

// factorized code that is used to symmetrize result of any spin operator
//
// state = reference on the state that has been produced with any operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int Spin1_2ChainMirrorSymmetry::SymmetrizeResult(unsigned long state, double& coefficient)
{
  unsigned long TmpState = this->GetSignedCanonicalState(state);
  if (this->SOperatorSignature != (TmpState & SPIN1_2CHAIN_MIRRORSYMMETRIC_BIT))
    {
      if (this->SOperatorSignature == 0x0ul)
	{
	  coefficient *= M_SQRT1_2;
	}
      else
	{
	  coefficient *= M_SQRT2;
	}
    }
  TmpState &= SPIN1_2CHAIN_FULLSYMMETRY_MASK;
  if (state > TmpState)
    coefficient *= this->MirrorParitySign;
  return this->FindStateIndex(TmpState);
}


// get canonical expression of a given state and its symmetry
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state (with symmetry bit)

inline unsigned long Spin1_2ChainMirrorSymmetry::GetSignedCanonicalState (unsigned long initialState)
{
  unsigned long TmpState = this->ApplyMirrorSymmetry(initialState);
  if (TmpState < initialState)
    {
      return (TmpState | SPIN1_2CHAIN_MIRRORSYMMETRIC_BIT);
    }
  else
    {
      if (TmpState != initialState)
	{
	  return (initialState | SPIN1_2CHAIN_MIRRORSYMMETRIC_BIT);
	}
      else
	{
	  return initialState;
	}
    }
}

// get the mirror symmetric of a given state 
//
// initialState = state to which the mirror symmetry has to be applied
// return value = flipped configuration

inline unsigned long Spin1_2ChainMirrorSymmetry::ApplyMirrorSymmetry (unsigned long initialState)
{
  initialState <<= this->MirrorShift;
#ifdef __64_BITS__
  unsigned long TmpState = Spin1_2ChainMirrorSymmetryMirrorTable[initialState & 0xfful] << 56;
  TmpState |= Spin1_2ChainMirrorSymmetryMirrorTable[(initialState >> 8) & 0xfful] << 48;
  TmpState |= Spin1_2ChainMirrorSymmetryMirrorTable[(initialState >> 16) & 0xfful] << 40;
  TmpState |= Spin1_2ChainMirrorSymmetryMirrorTable[(initialState >> 24) & 0xfful] << 32;
  TmpState |= Spin1_2ChainMirrorSymmetryMirrorTable[(initialState >> 32) & 0xfful] << 24;
  TmpState |= Spin1_2ChainMirrorSymmetryMirrorTable[(initialState >> 40) & 0xfful] << 16;
  TmpState |= Spin1_2ChainMirrorSymmetryMirrorTable[(initialState >> 48) & 0xfful] << 8;
  TmpState |= Spin1_2ChainMirrorSymmetryMirrorTable[(initialState >> 56) & 0xfful]; 
#else
  unsigned long TmpState = Spin1_2ChainMirrorSymmetryMirrorTable[initialState & 0xfful] << 24;
  TmpState |= Spin1_2ChainMirrorSymmetryMirrorTable[(initialState >> 8) & 0xfful] << 16;
  TmpState |= Spin1_2ChainMirrorSymmetryMirrorTable[(initialState >> 16) & 0xfful] << 8;
  TmpState |= Spin1_2ChainMirrorSymmetryMirrorTable[(initialState >> 24) & 0xfful];
#endif	
  TmpState >>= this->MirrorUnshift;
  return TmpState;
}

#endif


