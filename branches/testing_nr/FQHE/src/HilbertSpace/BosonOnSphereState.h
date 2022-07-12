////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of state of boson on a sphere                  //
//                                                                            //
//                        last modification : 03/10/2004                      //
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


#ifndef BOSONONSPHERESTATE_H
#define BOSONONSPHERESTATE_H


#include "config.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;
using std::dec;
using std::hex;


// get the reduced number of state (aka the number of unsigned long per state) minus 1
// 
// nbrState = number of states in which bosons can lie
// return value = reduced number of state
int GetReducedNbrState (int nbrState);

// get the number of the state in the last unsigned long array describing the whole state
// 
// nbrState = number of states in which bosons can lie
// return value = number of the state in the last unsigned long array describing the whole state
int GetRemainderNbrState (int nbrState);


class BosonOnSphereState
{

 private:

  // array describing the boson state
  unsigned long* StateDescription;

 public:

  // default constructor
  // 
  BosonOnSphereState();
  
  // basic constructor
  // 
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  BosonOnSphereState(const int& reducedNbrState);
  
  // copy constructor
  // 
  // state = reference on the state to copy
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  BosonOnSphereState(BosonOnSphereState& state, const int& reducedNbrState);

  // destructor
  // 
  ~BosonOnSphereState();

  // resize the current state
  //
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = reference on  the current state
  BosonOnSphereState& Resize(const int& reducedNbrState);

  // assign a state to the current one
  //
  // state = reference on the state to assign
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = reference on  the current state
  BosonOnSphereState& Assign(BosonOnSphereState& state, const int& reducedNbrState);

  // assign a state to the current one and after undefined it (array tranfert)
  //
  // state = reference on the state to assign
  // return value = reference on  the current state
  BosonOnSphereState& TransfertState(BosonOnSphereState& state);

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
  BosonOnSphereState& EmptyState(const int& reducedNbrState);

  // set occupation of a state 
  //
  // stateIndex = index of the state whose occupation has to be set 
  // nbrBosons = number of bosons in the given state
  void SetOccupation (const int& stateIndex, const unsigned long& nbrBosons);

  // get occupation of a state 
  //
  // stateIndex = index of the state whose occupation has to be set 
  // return value = number of bosons in the given state
  unsigned long GetOccupation (const int& stateIndex);

  // increment occupation of a state 
  //
  // stateIndex = index of the state whose occupation has to be set 
  // tmpShift1 = refence on temporary variable used during incrementation
  // tmpShift2 = refence on temporary variable used to evaluated a bit shift
  // return = occupation of the stateIndex-th state after the increment operation
  unsigned long IncrementOccupation (const int& stateIndex, unsigned long& tmpShift1, int& tmpShift2);
  
  // decrement occupation of a state (without testing if the state es empty)
  //
  // stateIndex = index of the state whose occupation has to be set 
  // tmpShift1 = refence on temporary variable used during decrementation
  // tmpShift2 = refence on temporary variable used to evaluated a bit shift
  // return = occupation of the stateIndex-th state before the decrement operation
  unsigned long DecrementOccupation (const int& stateIndex, unsigned long& tmpShift1, int& tmpShift2);
  
  // test if the state is empty an if it is not, decrement its occupation
  //
  // stateIndex = index of the state whose occupation has to be set 
  // tmpShift1 = refence on temporary variable used to evaluated a bit shift
  // tmpShift2 = refence on temporary variable used to evaluated a bit shift
  // return value = false if the state is empty
  inline bool TestAndDecrementOccupation (const int& stateIndex, int& tmpShift1, int& tmpShift2);
  
  // swap two states
  //
  // state = reference on the state to swap with the current one
  void SwapStates (BosonOnSphereState& state);

  // test if the current state is identical to another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the two states are identical
  bool Equal (BosonOnSphereState& state, int reducedNbrState);

  // test if the current state is different to another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the two states are different
  bool Different (BosonOnSphereState& state, int reducedNbrState);
  
  // test if the current state is greater than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater than the other state
  bool Greater (BosonOnSphereState& state, int reducedNbrState);

  // test if the current state is greater or equal than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater or equal than the other state
  bool GreaterOrEqual (BosonOnSphereState& state, int reducedNbrState);

  // test if the current state is lesser than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater than the other state
  bool Lesser (BosonOnSphereState& state, const int reducedNbrState);

  // test if the current state is lesser or equal than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater or equal than the other state
  bool LesserOrEqual (BosonOnSphereState& state, int reducedNbrState);

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
// nbrState = number of states in which bosons can lie
// return value = reduced number of state

inline int GetReducedNbrState (int nbrState)
{
#ifdef __64_BITS__
  if ((nbrState & ((unsigned long) 0x7)) == 0)
    return ((nbrState >> 3) - 1);
  else
    return (nbrState >> 3);
#else
  if ((nbrState & ((unsigned long) 0x3)) == 0)
    return ((nbrState >> 2) - 1);
  else
    return (nbrState >> 2);
#endif
}

// get the number of the state in the last unsigned long array describing the whole state
// 
// nbrState = number of states in which bosons can lie
// return value = number of the state in the last unsigned long array describing the whole state

inline int GetRemainderNbrState (int nbrState)
{
#ifdef __64_BITS__
  if ((nbrState & ((unsigned long) 0x7)) == 0)
    return 8;
  else
    return (nbrState & ((unsigned long) 0x7));
#else
  if ((nbrState & ((unsigned long) 0x3)) == 0)
    return 4;
  else
    return (nbrState & ((unsigned long) 0x3));
#endif
}

// resize the current state
//
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// return value = reference on  the current state

inline BosonOnSphereState& BosonOnSphereState::Resize(const int& reducedNbrState)
{
  if (this->StateDescription != 0)
    delete[] this->StateDescription;
  this->StateDescription = new unsigned long [reducedNbrState + 1];
  for (int i = 0; i <= reducedNbrState; ++i)
    this->StateDescription[i] = 0l;  
  return *this;
}

// assign a state to the current one
//
// state = reference on the state to assign
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// return value = reference on  the current state

inline BosonOnSphereState& BosonOnSphereState::Assign(BosonOnSphereState& state, const int& reducedNbrState)
{
  for (int i = 0; i <= reducedNbrState; ++i)
    this->StateDescription[i] = state.StateDescription[i];
  return *this;
}

// assign a state to the current one and after undefined it (array tranfert)
//
// state = reference on the state to assign
// return value = reference on  the current state

inline BosonOnSphereState& BosonOnSphereState::TransfertState(BosonOnSphereState& state)
{
  this->StateDescription = state.StateDescription;
  state.StateDescription = 0;
  return *this;
}

// set all ocupations to zero
//
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// return value = reference on  the current state

inline BosonOnSphereState& BosonOnSphereState::EmptyState(const int& reducedNbrState)
{
  for (int i = 0; i <= reducedNbrState; ++i)
    this->StateDescription[i] = 0l;
  return *this;
}

// get hash key associated to the state
//
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// keyMask = mask giving how many bits have to be kept
// return value = hash key

inline unsigned long BosonOnSphereState::GetHashKey (const int& reducedNbrState, const unsigned long& keyMask)
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

inline int BosonOnSphereState::GetHighestIndex (const int& reducedNbrState)
{
  int Index = 0;
  for (int i = reducedNbrState; ((i >= 0) && (Index == 0)); --i)
    {
#ifdef __64_BITS__
      if ((this->StateDescription[i] & ((unsigned long) 0xff00000000000000)) != 0)
	Index = (i << 3) + 7;
      else
	if ((this->StateDescription[i] & ((unsigned long) 0xff000000000000)) != 0)
	  Index = (i << 3) + 6;
	else
	  if ((this->StateDescription[i] & ((unsigned long) 0xff0000000000)) != 0)
	    Index = (i << 3) + 5;
	  else
	    if ((this->StateDescription[i] & ((unsigned long) 0xff00000000)) != 0)
	      Index = (i << 3) + 4;
	    else
	      if ((this->StateDescription[i] & ((unsigned long) 0xff000000)) != 0)
		Index = (i << 3) + 3;
	      else
		if ((this->StateDescription[i] & ((unsigned long) 0xff0000)) != 0)
		  Index = (i << 3) + 2;
		else
		  if ((this->StateDescription[i] & ((unsigned long) 0xff00)) != 0)
		    Index = (i << 3) + 1;
		  else
		    if ((this->StateDescription[i] & ((unsigned long) 0xff)) != 0)
		      Index = (i << 3);
#else
      if ((this->StateDescription[i] & ((unsigned long) 0xff000000)) != 0)
	Index = (i << 2) + 3;
      else
	if ((this->StateDescription[i] & ((unsigned long) 0xff0000)) != 0)
	  Index = (i << 2) + 2;
	else
	  if ((this->StateDescription[i] & ((unsigned long) 0xff00)) != 0)
	    Index = (i << 2) + 1;
	  else
	    if ((this->StateDescription[i] & ((unsigned long) 0xff)) != 0)
	      Index = (i << 2);
#endif
    }
  return Index;  
}

// set occupation of a state 
//
// stateIndex = index of the state whose occupation has to be set 
// nbrBosons = number of bosons in the given state

inline void BosonOnSphereState::SetOccupation (const int& stateIndex, const unsigned long& nbrBosons)
{
#ifdef __64_BITS__
  this->StateDescription[stateIndex >> 3] &= ~(((unsigned long) 0xff) << ((stateIndex & ((unsigned long) 0x7)) << 3));
  this->StateDescription[stateIndex >> 3] |= nbrBosons << ((stateIndex & ((unsigned long) 0x7)) << 3);
#else
  this->StateDescription[stateIndex >> 2] &= ~(((unsigned long) 0xff) << ((stateIndex & ((unsigned long) 0x3)) << 3));
  this->StateDescription[stateIndex >> 2] |= nbrBosons << ((stateIndex & ((unsigned long) 0x3)) << 3);
#endif
}

// get occupation of a state 
//
// stateIndex = index of the state whose occupation has to be set 
// return value = number of bosons in the given state

inline unsigned long BosonOnSphereState::GetOccupation (const int& stateIndex)
{
#ifdef __64_BITS__
  return ((this->StateDescription[stateIndex >> 3] >> ((stateIndex & ((unsigned long) 0x7)) << 3)) & ((unsigned long) 0xff));
#else
  return ((this->StateDescription[stateIndex >> 2] >> ((stateIndex & ((unsigned long) 0x3)) << 3)) & ((unsigned long) 0xff));
#endif
}

// increment occupation of a state 
//
// stateIndex = index of the state whose occupation has to be set 
// tmpShift1 = refence on temporary variable used during incrementation
// tmpShift2 = refence on temporary variable used to evaluated a bit shift
// return = occupation of the stateIndex-th state after the increment operation

inline unsigned long BosonOnSphereState::IncrementOccupation (const int& stateIndex, unsigned long& tmpShift1, int& tmpShift2)
{
#ifdef __64_BITS__
  tmpShift2 = (stateIndex & ((unsigned long) 0x7)) << 3;
  tmpShift1 = ((this->StateDescription[stateIndex >> 3] >> tmpShift2) & ((unsigned long) 0xff)) + 1;
  this->StateDescription[stateIndex >> 3] &= ~(((unsigned long) 0xff) << tmpShift2);
  this->StateDescription[stateIndex >> 3] |= tmpShift1 << tmpShift2;
#else
  tmpShift2 = (stateIndex & ((unsigned long) 0x3)) << 3;
  tmpShift1 = ((this->StateDescription[stateIndex >> 2] >> tmpShift2) & ((unsigned long) 0xff)) + 1;
  this->StateDescription[stateIndex >> 2] &= ~(((unsigned long) 0xff) << tmpShift2);
  this->StateDescription[stateIndex >> 2] |= tmpShift1 << tmpShift2;
#endif
  return tmpShift1;
}

// decrement occupation of a state (without testing if the state es empty)
//
// stateIndex = index of the state whose occupation has to be set 
// tmpShift1 = refence on temporary variable used during decrementation
// tmpShift2 = refence on temporary variable used to evaluated a bit shift

inline unsigned long BosonOnSphereState::DecrementOccupation (const int& stateIndex, unsigned long& tmpShift1, int& tmpShift2)
{
#ifdef __64_BITS__
  tmpShift2 = (stateIndex & ((unsigned long) 0x7)) << 3;
  tmpShift1 = ((this->StateDescription[stateIndex >> 3] >> tmpShift2) & ((unsigned long) 0xff));
  this->StateDescription[stateIndex >> 3] &= ~(((unsigned long) 0xff) << tmpShift2);
  this->StateDescription[stateIndex >> 3] |= (tmpShift1 - 1) << tmpShift2;
#else
  tmpShift2 = (stateIndex & ((unsigned long) 0x3)) << 3;
  tmpShift1 = ((this->StateDescription[stateIndex >> 2] >> tmpShift2) & ((unsigned long) 0xff));
  this->StateDescription[stateIndex >> 2] &= ~(((unsigned long) 0xff) << tmpShift2);
  this->StateDescription[stateIndex >> 2] |= (tmpShift1 - 1) << tmpShift2;
#endif
  return tmpShift1;
}

// test if the state is empty an if it is not, decrement its occupation
//
// stateIndex = index of the state whose occupation has to be set 
// tmpShift1 = refence on temporary variable used to evaluated a bit shift
// tmpShift2 = refence on temporary variable used to evaluated a bit shift
// return value = false if the state is empty

inline bool BosonOnSphereState::TestAndDecrementOccupation (const int& stateIndex, int& tmpShift1, int& tmpShift2)
{
#ifdef __64_BITS__
  tmpShift1 = stateIndex >> 3;
  tmpShift2 = (stateIndex & ((unsigned long) 0x7)) << 3;
#else
  tmpShift1 = stateIndex >> 2;
  tmpShift2 = (stateIndex & ((unsigned long) 0x3)) << 3;
#endif
  if (this->StateDescription[tmpShift1] >> tmpShift2)
    {
      this->StateDescription[tmpShift1] &= ~(((unsigned long) 0xff) << tmpShift2);
      this->StateDescription[tmpShift1] |= (((this->StateDescription[tmpShift1] >> tmpShift2) & ((unsigned long) 0xff)) - 1) << tmpShift2;
      return true;
    }
  else
    return false;
}

// swap two states
//
// state = reference on the state to swap with the current one

inline void BosonOnSphereState::SwapStates (BosonOnSphereState& state)
{
  unsigned long* TmpStateDescription = this->StateDescription;
  this->StateDescription = state.StateDescription;
  state.StateDescription = TmpStateDescription;
}

// test if the current state is identical to another state
//
// state = reference on the state to compare with
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the two states are identical

inline bool BosonOnSphereState::Equal (BosonOnSphereState& state, int reducedNbrState)
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
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the two states are different

inline bool BosonOnSphereState::Different (BosonOnSphereState& state, int reducedNbrState)
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
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the current state is greater than the other state

inline bool BosonOnSphereState::Greater (BosonOnSphereState& state, int reducedNbrState)
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
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the current state is greater or equal than the other state

inline bool BosonOnSphereState::GreaterOrEqual (BosonOnSphereState& state, int reducedNbrState)
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

inline bool BosonOnSphereState::Lesser (BosonOnSphereState& state, int reducedNbrState)
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

inline bool BosonOnSphereState::LesserOrEqual (BosonOnSphereState& state, int reducedNbrState)
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

inline ostream& BosonOnSphereState::PrintState (ostream& str, const int& reducedNbrState, const int& nbrStateRemainder)
{
  for (int i = 0; i < reducedNbrState; ++i)
    {
#ifdef __64_BITS__
      for (int j = 0; j < 8; ++j)
#else
      for (int j = 0; j < 4; ++j)
#endif
	str << ((this->StateDescription[i] >> (j << 3)) & ((unsigned long) 0xff)) << " ";
    }
  for (int i = 0; i < nbrStateRemainder; ++i)
    str << ((this->StateDescription[reducedNbrState] >> (i << 3)) & ((unsigned long) 0xff)) << " ";
  return str;
}


#endif
