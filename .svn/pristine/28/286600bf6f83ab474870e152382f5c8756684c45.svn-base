////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//   class of spin 1/2 chain without any Sz contraint but a fixed Sz parity   //
//                                                                            //
//                        last modification : 06/01/2014                      //
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


#ifndef SPIN1_2CHAINFIXEDPARITY_H
#define SPIN1_2CHAINFIXEDPARITY_H


#include "config.h"
#include "HilbertSpace/Spin1_2ChainFull.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;


class Spin1_2ChainFixedParity : public Spin1_2ChainFull
{

 protected:

  // parity of the total (Sz + 1/2) (can be 0 or 1)
  int SzParity;

 public:


  // default constructor
  //
  Spin1_2ChainFixedParity ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin 1/2
  // parity = parity of the total (Sz + 1/2) (can be 0 or 1)
  Spin1_2ChainFixedParity (int chainLength, int parity);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1_2ChainFixedParity (const Spin1_2ChainFixedParity& chain);

  // destructor
  //
  ~Spin1_2ChainFixedParity ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1_2ChainFixedParity& operator = (const Spin1_2ChainFixedParity& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

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

  // compute the parity (prod_i Sz_i) for a given state
  //
  // state = index of the state to be applied on Sz_i operator
  // return value = 0 if prod_i Sz_i = 1, 1 if prod_i Sz_i = -1
  virtual unsigned long Parity (int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // shift = position of the A part leftmost site within the full system
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int nbrSites, int szSector, int shift, RealVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // shift = position of the A part leftmost site within the full system
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, int shift, RealVector& groundState, AbstractArchitecture* architecture = 0);

 protected:

  // find state index
  //
  // state = state description
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long state);

};

// find state index
//
// state = state description
// return value = corresponding index

inline int Spin1_2ChainFixedParity::FindStateIndex(unsigned long state)
{
  return ((int) (state >> 1));    
}

// compute the parity (prod_i Sz_i) for a given state
//
// state = index of the state to be applied on Sz_i operator
// return value = 0 if prod_i Sz_i = 1, 1 if prod_i Sz_i = -1

inline unsigned long Spin1_2ChainFixedParity::Parity (int state)
{
  return (unsigned long) this->SzParity;
}

#endif


