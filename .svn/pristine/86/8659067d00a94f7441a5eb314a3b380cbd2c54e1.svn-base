////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of spin 1/2 chain with a fixed Sz value                //
//         and a fixed parity under the parity symmetry (aka Sz<->-Sz)        //
//                                                                            //
//                        last modification : 21/07/2015                      //
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


#ifndef SPIN1_2CHAINPARITYSYMMETRY_H
#define SPIN1_2CHAINPARITYSYMMETRY_H


#include "config.h"
#include "HilbertSpace/Spin1_2ChainMirrorSymmetry.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;
using std::hex;
using std::dec;


class Spin1_2ChainParitySymmetry : public Spin1_2ChainMirrorSymmetry
{

 protected:

  // mask to apply during the parity symmetry
  unsigned long ParityMask;

  // parity under the parity symmetry
  int Parity;
  // additional sign due to the parity symmetry 
  double ParitySign;


 public:


  // default constructor
  //
  Spin1_2ChainParitySymmetry ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin 1/2
  // parity = parity under the parity symmetry
  // memorySize = memory size in bytes allowed for look-up table
  Spin1_2ChainParitySymmetry (int chainLength, int parity, int memorySize);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1_2ChainParitySymmetry (const Spin1_2ChainParitySymmetry& chain);

  // destructor
  //
  ~Spin1_2ChainParitySymmetry ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1_2ChainParitySymmetry& operator = (const Spin1_2ChainParitySymmetry& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

 protected:

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

  // get the parity symmetric of a given state 
  //
  // initialState = state to which the parity symmetry has to be applied
  // return value = flipped configuration
  virtual unsigned long ApplyParitySymmetry (unsigned long initialState);

};

// factorized code that is used to symmetrize result of any spin operator
//
// state = reference on the state that has been produced with any operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int Spin1_2ChainParitySymmetry::SymmetrizeResult(unsigned long state, double& coefficient)
{
  unsigned long TmpState = this->GetSignedCanonicalState(state);
  if (this->SOperatorSignature != (TmpState & SPIN1_2CHAIN_PARITYSYMMETRIC_BIT))
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
    coefficient *= this->ParitySign;
  return this->FindStateIndex(TmpState);
}


// get canonical expression of a given state and its symmetry
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state (with symmetry bit)

inline unsigned long Spin1_2ChainParitySymmetry::GetSignedCanonicalState (unsigned long initialState)
{
  unsigned long TmpState = this->ApplyParitySymmetry(initialState);
  if (TmpState < initialState)
    {
      return (TmpState | SPIN1_2CHAIN_PARITYSYMMETRIC_BIT);
    }
  else
    {
      if (TmpState != initialState)
	{
	  return (initialState | SPIN1_2CHAIN_PARITYSYMMETRIC_BIT);
	}
      else
	{
	  return initialState;
	}
    }
}

// get the parity symmetric of a given state 
//
// initialState = state to which the parity symmetry has to be applied
// return value = flipped configuration

inline unsigned long Spin1_2ChainParitySymmetry::ApplyParitySymmetry (unsigned long initialState)
{
  return ((~initialState) & this->ParityMask);
}

#endif


