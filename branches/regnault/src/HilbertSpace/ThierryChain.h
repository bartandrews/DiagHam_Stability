////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class of spin 1 chain                           //
//                                                                            //
//                        last modification : 20/02/2001                      //
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


#ifndef THIERRYCHAIN_H
#define THIERRYCHAIN_H


#include "config.h"

#include <iostream>


using std::ostream;


class ThierryChain
{

 protected:

  int* GarbageFlag;

  int Spin2ChainLength;
  int Spin3_2ChainLength;
  int Spin3_2Start;
  int HilbertSpaceDimension;

  int Sz;
  int MaxSpin3_2ContributionToSz;

  unsigned long* ChainDescription;
  
  int* LookUpTable;
  unsigned long LookUpTableMask;
  int LookUpPosition;
  int LookUpPosition2;

 public:
 
  // default constructor
  //
  ThierryChain ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // spin2ChainLength = number of spin 2
  // spin3_2ChainLength = number of spin 3/2
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  ThierryChain (int spin2ChainLength, int spin3_2ChainLength, int sz, int memorySize);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  ThierryChain (const ThierryChain& chain);

  // destructor
  //
  ~ThierryChain ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  ThierryChain& operator = (const ThierryChain& chain);

  // re-initialize chain with another total Sz component
  //
  // sz = twice the value of total Sz component
  // return value = reference on current chain  
  ThierryChain& Reinitialize(int sz);

  // return Hilbert space dimension
  //
  // return value = Hilbert space dimension
  int GetHilbertSpaceDimension() {return this->HilbertSpaceDimension;};

  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = pointer to double where numerical coefficient has to be stored
  // return value = index of resulting state
  int SmiSpj (int i, int j, int state, double* coefficient);

  // return eigenvalue of Sz_i Sz_j associated to a given state
  //
  // i = position of left Sz operator
  // j = position of right Sz operator
  // state = index of the state to consider
  // return value = corresponding eigenvalue
  double SziSzj (int i, int j, int state);

  // find state index
  //
  // state = state description
  // return value = corresponding index
  int FindStateIndex(unsigned long state);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

 private:

  // generate all states
  //
  // initialStatePosition = position of the state to be duplicated
  // statePosition = position for the new states
  // sitePosition = site on chain where spin has to be changed
  // return value = number of generated states
  int GenerateSpin2States(int initialStatePosition, int statePosition, int sitePosition);

  // generate all states
  //
  // initialStatePosition = position of the state to be duplicated
  // statePosition = position for the new states
  // sitePosition = site on chain where spin has to be changed
  // return value = number of generated states
  int GenerateSpin3_2States(int initialStatePosition, int statePosition, int sitePosition);

  // generate Spin 2 states for a given total spin projection Sz
  //
  // statePosition = position for the new states
  // sitePosition = site on chain where spin has to be changed
  // currentStateDescription = description of current state
  // currentSz = total Sz value of current state
  // return value = number of generated states
  int GenerateSpin2States(int statePosition, int sitePosition, unsigned currentStateDescription, int currentSz);

  // generate Spin 3/2 states for a given total spin projection Sz
  //
  // statePosition = position for the new states
  // sitePosition = site on chain where spin has to be changed
  // currentStateDescription = description of current state
  // currentSz = total Sz value of current state
  // return value = number of generated states  
  int GenerateSpin3_2States(int statePosition, int sitePosition, unsigned currentStateDescription, int currentSz);

};

#endif


