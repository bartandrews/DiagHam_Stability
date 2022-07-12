////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of state of fermion on a sphere without any restriction       //
//    on the number of fermion nor the number of states that can be reached   //
//                                                                            //
//                        last modification : 02/03/2004                      //
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


#ifndef FERMIONONSPHERELONGSTATE_H
#define FERMIONONSPHERELONGSTATE_H


#include "config.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;
using std::dec;
using std::hex;


// get the reduced number of state (aka the number of unsigned long per state) minus 1
// 
// nbrState = number of states in which fermions can lie
// return value = reduced number of state
int FermionOnSphereLongStateGetReducedNbrState (int nbrState);

// get the number of the state in the last unsigned long array describing the whole state
// 
// nbrState = number of states in which fermions can lie
// return value = number of the state in the last unsigned long array describing the whole state
int FermionOnSphereLongStateGetRemainderNbrState (int nbrState);


class FermionOnSphereLongState
{

 private:

  // array describing the fermion state
  unsigned long* StateDescription;

 public:

  // default constructor
  // 
  FermionOnSphereLongState();
  
  // basic constructor
  // 
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  FermionOnSphereLongState(const int& reducedNbrState);
  
  // copy constructor
  // 
  // state = reference on the state to copy
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  FermionOnSphereLongState(FermionOnSphereLongState& state, const int& reducedNbrState);

  // destructor
  // 
  ~FermionOnSphereLongState();

  // resize the current state
  //
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = reference on  the current state
  FermionOnSphereLongState& Resize(const int& reducedNbrState);

  // assign a state to the current one
  //
  // state = reference on the state to assign
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = reference on  the current state
  FermionOnSphereLongState& Assign(FermionOnSphereLongState& state, const int& reducedNbrState);

  // assign a state to the current one and after undefined it (array tranfert)
  //
  // state = reference on the state to assign
  // return value = reference on  the current state
  FermionOnSphereLongState& TransfertState(FermionOnSphereLongState& state);

  // get hash key associated to the state
  //
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // keyMask = mask giving how many bits have to be kept
  // return value = hash key
  unsigned long GetHashKey (const int& reducedNbrState, const unsigned long& keyMask);

  // get the highest state index for which the state is occupied
  //
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = highest state index
  int GetHighestIndex (const int& reducedNbrState);

  // set all ocupations to zero
  //
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = reference on  the current state
  FermionOnSphereLongState& EmptyState(const int& reducedNbrState);

  // set occupation of a state to one
  //
  // stateIndex = index of the state whose occupation has to be set 
  void SetOccupation (const int& stateIndex);

  // set occupation of a state to zero
  //
  // stateIndex = index of the state whose occupation has to be set 
  void UnsetOccupation (const int& stateIndex);

  // get occupation of a state 
  //
  // stateIndex = index of the state whose occupation has to be set 
  // return value = number of fermions in the given state
  unsigned long GetOccupation (const int& stateIndex);

  // increment occupation of a state 
  //
  // stateIndex = index of the state whose occupation has to be set 
  void IncrementOccupation (const int& stateIndex);
  
  // test if the state is full and if it is not, increment its occupation
  //
  // stateIndex = index of the state whose occupation has to be set 
  // return value = false if the state is empty
  bool TestAndIncrementOccupation (const int& stateIndex);
  
  // decrement occupation of a state (without testing if the state es empty)
  //
  // stateIndex = index of the state whose occupation has to be set 
  void DecrementOccupation (const int& stateIndex);
  
  // test if the state is empty and if it is not, decrement its occupation
  //
  // stateIndex = index of the state whose occupation has to be set 
  // return value = false if the state is empty
  bool TestAndDecrementOccupation (const int& stateIndex);

  // stateIndex = index of the state whose occupation has to be set 
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // signLookUpTable = pointer to an array containing the parity of the number of one for each integer ranging from 0 to 65535
  // signLookUpTableMaks = pointer to an array containing the parity of the mask on the bits to keep for each shift that is requested by sign evaluation
  // coefficient = reference on a coefficient which will be multiplied by the sign of the permutation
  void GetPermutationSign(int stateIndex, int reducedNbrState, double* signLookUpTable, unsigned long* signLookUpTableMaks, double& coefficient);
  
  // swap two states
  //
  // state = reference on the state to swap with the current one
  void SwapStates (FermionOnSphereLongState& state);

  // test if the current state is identical to another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the two states are identical
  bool Equal (FermionOnSphereLongState& state, int reducedNbrState);

  // test if the current state is different to another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the two states are different
  bool Different (FermionOnSphereLongState& state, int reducedNbrState);
  
  // test if the current state is greater than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater than the other state
  bool Greater (FermionOnSphereLongState& state, int reducedNbrState);

  // test if the current state is greater or equal than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater or equal than the other state
  bool GreaterOrEqual (FermionOnSphereLongState& state, int reducedNbrState);

  // test if the current state is lesser than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater than the other state
  bool Lesser (FermionOnSphereLongState& state, const int reducedNbrState);

  // test if the current state is lesser or equal than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater or equal than the other state
  bool LesserOrEqual (FermionOnSphereLongState& state, int reducedNbrState);

  // print a given state
  //
  // str = reference on current output stream 
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
  // return value = reference on current output stream 
  ostream& PrintState (ostream& str, const int& reducedNbrState, const int& nbrStateRemainder);

};

// get the reduced number of state (aka the number of unsigned long per state) minus 1
// 
// nbrState = number of states in which fermions can lie
// return value = reduced number of state

inline int FermionOnSphereLongStateGetReducedNbrState (int nbrState)
{
#ifdef __64_BITS__
  if ((nbrState & ((unsigned long) 0x3f)) == 0)
    return ((nbrState >> 6) - 1);
  else
    return (nbrState >> 6);
#else
  if ((nbrState & ((unsigned long) 0x1f)) == 0)
    return ((nbrState >> 5) - 1);
  else
    return (nbrState >> 5);
#endif
}

// get the number of the state in the last unsigned long array describing the whole state
// 
// nbrState = number of states in which fermions can lie
// return value = number of the state in the last unsigned long array describing the whole state

inline int FermionOnSphereLongStateGetRemainderNbrState (int nbrState)
{
#ifdef __64_BITS__
  if ((nbrState & 0x3f) == 0)
    return 64;
  else
    return (nbrState & 0x3f);
#else
  if ((nbrState & 0x1f) == 0)
    return 32;
  else
    return (nbrState & 0x1f);
#endif
}

// resize the current state
//
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// return value = reference on  the current state

inline FermionOnSphereLongState& FermionOnSphereLongState::Resize(const int& reducedNbrState)
{
  if (this->StateDescription != 0)
    delete[] this->StateDescription;
  this->StateDescription = new unsigned long [reducedNbrState + 1];
  for (int i = 0; i <= reducedNbrState; ++i)
    this->StateDescription[i] = 0;  
  return *this;
}

// assign a state to the current one
//
// state = reference on the state to assign
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// return value = reference on  the current state

inline FermionOnSphereLongState& FermionOnSphereLongState::Assign(FermionOnSphereLongState& state, const int& reducedNbrState)
{
  for (int i = 0; i <= reducedNbrState; ++i)
    this->StateDescription[i] = state.StateDescription[i];
  return *this;
}

// assign a state to the current one and after undefined it (array tranfert)
//
// state = reference on the state to assign
// return value = reference on  the current state

inline FermionOnSphereLongState& FermionOnSphereLongState::TransfertState(FermionOnSphereLongState& state)
{
  this->StateDescription = state.StateDescription;
  state.StateDescription = 0;
  return *this;
}

// set all ocupations to zero
//
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// return value = reference on  the current state

inline FermionOnSphereLongState& FermionOnSphereLongState::EmptyState(const int& reducedNbrState)
{
  for (int i = 0; i <= reducedNbrState; ++i)
    this->StateDescription[i] = (unsigned long) 0;
  return *this;
}

// get hash key associated to the state
//
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// keyMask = mask giving how many bits have to be kept
// return value = hash key

inline unsigned long FermionOnSphereLongState::GetHashKey (const int& reducedNbrState, const unsigned long& keyMask)
{
  unsigned long Key = this->StateDescription[0];
  for (int i = 1; i <= reducedNbrState; ++i)
    Key += (this->StateDescription[i] << i);
  return (Key & keyMask);
}

// get the highest state index for which the state is occupied
//
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// return value = highest state index

inline int FermionOnSphereLongState::GetHighestIndex (const int& reducedNbrState)
{
  int Index = 0;
  unsigned long Mask;
  int j;
  for (int i = reducedNbrState; ((i >= 0) && (Index == 0)); --i)
    {
#ifdef __64_BITS__
      Mask = (unsigned long) 0x8000000000000000;
      if ((this->StateDescription[i] & Mask) != 0)
	Index = (i << 6) + 63;
      else
	{
	  j = 62;
	  Mask = (unsigned long) 0x4000000000000000;	  
	  while ((Index == 0) && (j >= 0))	
	    {	      
	      if ((this->StateDescription[i] & Mask) != 0)
		Index = (i << 6) + j;
	      --j;
	      Mask >>= 1;
	    }
	}
#else
      Mask = (unsigned long) 0x80000000;
      if ((this->StateDescription[i] & Mask) != 0)
	Index = (i << 5) + 31;
      else
	{
	  j = 30;
	  Mask = (unsigned long) 0x40000000;	  
	  while ((Index == 0) && (j >= 0))	
	    {	      
	      if ((this->StateDescription[i] & Mask) != 0)
		Index = (i << 5) + j;
	      --j;
	      Mask >>= 1;
	    }
	}
#endif
    }
  return Index;  
}

// set occupation of a state to one
//
// stateIndex = index of the state whose occupation has to be set 

inline void FermionOnSphereLongState::SetOccupation (const int& stateIndex)
{
#ifdef __64_BITS__
  this->StateDescription[stateIndex >> 6] |=  ((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x3f));
#else
  this->StateDescription[stateIndex >> 5] |=  ((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x1f));
#endif
}

// set occupation of a state to zero
//
// stateIndex = index of the state whose occupation has to be set 

inline void FermionOnSphereLongState::UnsetOccupation (const int& stateIndex)
{
#ifdef __64_BITS__
  this->StateDescription[stateIndex >> 6] &=  ~(((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x3f)));
#else
  this->StateDescription[stateIndex >> 5] &=  ~(((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x1f)));
#endif
}

// get occupation of a state 
//
// stateIndex = index of the state whose occupation has to be set 
// return value = number of fermions in the given state

inline unsigned long FermionOnSphereLongState::GetOccupation (const int& stateIndex)
{
#ifdef __64_BITS__
  return (this->StateDescription[stateIndex >> 6] & (((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x3f))));
#else
  return (this->StateDescription[stateIndex >> 5] & (((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x1f))));
#endif
}

// increment occupation of a state 
//
// stateIndex = index of the state whose occupation has to be set 

inline void FermionOnSphereLongState::IncrementOccupation (const int& stateIndex)
{
#ifdef __64_BITS__
  this->StateDescription[stateIndex >> 6] |=  ((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x3f));
#else
  this->StateDescription[stateIndex >> 5] |=  ((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x1f));
#endif
}

// test if the state is full and if it is not, increment its occupation
//
// stateIndex = index of the state whose occupation has to be set 
// return value = false if the state is empty

inline bool FermionOnSphereLongState::TestAndIncrementOccupation (const int& stateIndex)
{
#ifdef __64_BITS__
  if ((this->StateDescription[stateIndex >> 6] & (((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x3f)))) != 0)
    return false;
  else
    {
      this->StateDescription[stateIndex >> 6] |=  (((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x3f)));
      return true;
    }
#else
  if ((this->StateDescription[stateIndex >> 5] & (((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x1f)))) != 0)
    return false;
  else
    {
      this->StateDescription[stateIndex >> 5] |=  (((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x1f)));
      return true;
    }
#endif
}
  
// decrement occupation of a state (without testing if the state is empty)
//
// stateIndex = index of the state whose occupation has to be set 

inline void FermionOnSphereLongState::DecrementOccupation (const int& stateIndex)
{
#ifdef __64_BITS__
  this->StateDescription[stateIndex >> 6] &= ~(((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x3f)));
#else
  this->StateDescription[stateIndex >> 5] &= ~(((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x1f)));
#endif
}

// test if the state is empty an if it is not, decrement its occupation
//
// stateIndex = index of the state whose occupation has to be set 
// return value = false if the state is empty

inline bool FermionOnSphereLongState::TestAndDecrementOccupation (const int& stateIndex)
{
#ifdef __64_BITS__
  if ((this->StateDescription[stateIndex >> 6] & (((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x3f)))) == 0)
    return false;
  else
    {
      this->StateDescription[stateIndex >> 6] &=  ~(((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x3f)));
      return true;
    }
#else
  if ((this->StateDescription[stateIndex >> 5] & (((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x1f)))) == 0)
    return false;
  else
    {
      this->StateDescription[stateIndex >> 5] &=  ~(((unsigned long) 0x1) << (stateIndex & ((unsigned long) 0x1f)));
      return true;
    }
#endif
}

// get the sign resulting from the permutation 
//
// stateIndex = index of the state whose occupation has to be set 
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// signLookUpTable = pointer to an array containing the parity of the number of one for each integer ranging from 0 to 65535
// signLookUpTableMaks = pointer to an array containing the parity of the mask on the bits to keep for each shift that is requested by sign evaluation
// coefficient = reference on a coefficient which will be multiplied by the sign of the permutation

inline void FermionOnSphereLongState::GetPermutationSign(int stateIndex, int reducedNbrState, double* signLookUpTable, unsigned long* signLookUpTableMaks, double& coefficient)
{
#ifdef  __64_BITS__
  int tmp = (stateIndex >> 6);
#else
  int tmp = (stateIndex >> 5);
#endif
  for (reducedNbrState; reducedNbrState > tmp; --reducedNbrState)
    {
#ifdef  __64_BITS__
      coefficient *= signLookUpTable[this->StateDescription[reducedNbrState] & ((unsigned long) 0xffff)];
      coefficient *= signLookUpTable[(this->StateDescription[reducedNbrState] >> 16) & ((unsigned long) 0xffff)];
      coefficient *= signLookUpTable[(this->StateDescription[reducedNbrState] >> 32) & ((unsigned long) 0xffff)];
      coefficient *= signLookUpTable[(this->StateDescription[reducedNbrState] >>48)  & ((unsigned long) 0xffff)];
#else
      coefficient *= signLookUpTable[this->StateDescription[reducedNbrState] & 0xffff];
      coefficient *= signLookUpTable[(this->StateDescription[reducedNbrState] >> 16) & 0xffff];
#endif
   }
#ifdef  __64_BITS__
  stateIndex &= 0x3f;
  coefficient *= signLookUpTable[(this->StateDescription[tmp] >> stateIndex) & signLookUpTableMaks[stateIndex]];
  coefficient *= signLookUpTable[(this->StateDescription[tmp] >> (stateIndex + 16)) & signLookUpTableMaks[stateIndex + 16]];
  coefficient *= signLookUpTable[(this->StateDescription[tmp] >> (stateIndex + 32))  & signLookUpTableMaks[stateIndex + 32]];
  coefficient *= signLookUpTable[(this->StateDescription[tmp] >> (stateIndex + 48))  & signLookUpTableMaks[stateIndex + 48]];
#else
  stateIndex &= 0x1f;
  coefficient *= signLookUpTable[(this->StateDescription[tmp] >> stateIndex) & signLookUpTableMaks[stateIndex]];
  coefficient *= signLookUpTable[(this->StateDescription[tmp] >> (stateIndex + 16)) & signLookUpTableMaks[stateIndex + 16]];
#endif
}
  
// swap two states
//
// state = reference on the state to swap with the current one

inline void FermionOnSphereLongState::SwapStates (FermionOnSphereLongState& state)
{
  unsigned long* TmpStateDescription = this->StateDescription;
  this->StateDescription = state.StateDescription;
  state.StateDescription = TmpStateDescription;
}

// test if the current state is identical to another state
//
// state = reference on the state to compare with
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the two states are identical

inline bool FermionOnSphereLongState::Equal (FermionOnSphereLongState& state, int reducedNbrState)
{
  while (reducedNbrState >= 0)
    if (state.StateDescription[reducedNbrState] != this->StateDescription[reducedNbrState])
      return false;
    else
      --reducedNbrState;
  return true;
}

// test if the current state is different to another state
//
// state = reference on the state to compare with
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the two states are different

inline bool FermionOnSphereLongState::Different (FermionOnSphereLongState& state, int reducedNbrState)
{
  while (reducedNbrState >= 0)
    if (state.StateDescription[reducedNbrState] != this->StateDescription[reducedNbrState])
      return true;
    else
      --reducedNbrState;
  return false;
}

// test if the current state is greater than another state
//
// state = reference on the state to compare with
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the current state is greater than the other state

inline bool FermionOnSphereLongState::Greater (FermionOnSphereLongState& state, int reducedNbrState)
{
  while (reducedNbrState >= 0)
    if (state.StateDescription[reducedNbrState] > this->StateDescription[reducedNbrState])
      return false;
    else
      if (state.StateDescription[reducedNbrState] == this->StateDescription[reducedNbrState])
	--reducedNbrState;
      else
	return true;
  return false;
}


// test if the current state is greater or equal than another state
//
// state = reference on the state to compare with
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the current state is greater or equal than the other state

inline bool FermionOnSphereLongState::GreaterOrEqual (FermionOnSphereLongState& state, int reducedNbrState)
{
  while (reducedNbrState >= 0)
    if (state.StateDescription[reducedNbrState] > this->StateDescription[reducedNbrState])
      return false;
    else
      if (state.StateDescription[reducedNbrState] == this->StateDescription[reducedNbrState])
	--reducedNbrState;
      else
	return true;
  return true;
}

// test if the current state is lesser than another state
//
// state = reference on the state to compare with
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the current state is lower than the other state

inline bool FermionOnSphereLongState::Lesser (FermionOnSphereLongState& state, int reducedNbrState)
{
  while (reducedNbrState >= 0)
    if (state.StateDescription[reducedNbrState] < this->StateDescription[reducedNbrState])
      return false;
    else
      if (state.StateDescription[reducedNbrState] == this->StateDescription[reducedNbrState])
	--reducedNbrState;
      else
	return true;
  return false;
}


// test if the current state is lesser or equal than another state
//
// state = reference on the state to compare with
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the current state is lower or equal than the other state

inline bool FermionOnSphereLongState::LesserOrEqual (FermionOnSphereLongState& state, int reducedNbrState)
{
  while (reducedNbrState >= 0)
    if (state.StateDescription[reducedNbrState] < this->StateDescription[reducedNbrState])
      return false;
    else
      if (state.StateDescription[reducedNbrState] == this->StateDescription[reducedNbrState])
	--reducedNbrState;
      else
	return true;
  return true;
}

// print a given state
//
// str = reference on current output stream 
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
// return value = reference on current output stream 

inline ostream& FermionOnSphereLongState::PrintState (ostream& str, const int& reducedNbrState, const int& nbrStateRemainder)
{
  for (int i = 0; i < reducedNbrState; ++i)
    {
#ifdef __64_BITS__
      for (int j = 0; j < 64; ++j)
#else
      for (int j = 0; j < 32; ++j)
#endif
	str << ((this->StateDescription[i] >> j) & ((unsigned long) 0x1)) << " ";
    }
  for (int i = 0; i < nbrStateRemainder; ++i)
    str << ((this->StateDescription[reducedNbrState] >> i) & ((unsigned long) 0x1)) << " ";
  return str;
}


#endif
