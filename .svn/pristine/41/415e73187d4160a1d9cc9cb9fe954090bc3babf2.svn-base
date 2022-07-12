////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of state of boson on a torus                   //
//                                                                            //
//                        last modification : 14/10/2003                      //
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


#ifndef BOSONONTORUSSTATE_H
#define BOSONONTORUSSTATE_H


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


class BosonOnTorusState
{

 private:

  // array describing the boson state
  unsigned long* StateDescription;

 public:

  // default constructor
  // 
  BosonOnTorusState();
  
  // basic constructor
  // 
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  BosonOnTorusState(const int& reducedNbrState);
  
  // copy constructor
  // 
  // state = reference on the state to copy
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  BosonOnTorusState(BosonOnTorusState& state, const int& reducedNbrState);

  // destructor
  // 
  ~BosonOnTorusState();

  // resize the current state
  //
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = reference on  the current state
  BosonOnTorusState& Resize(const int& reducedNbrState);

  // assign a state to the current one
  //
  // state = reference on the state to assign
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = reference on  the current state
  BosonOnTorusState& Assign(BosonOnTorusState& state, const int& reducedNbrState);

  // assign a state to the current one and after undefined it (array tranfert)
  //
  // state = reference on the state to assign
  // return value = reference on  the current state
  BosonOnTorusState& TransfertState(BosonOnTorusState& state);

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
  BosonOnTorusState& EmptyState(const int& reducedNbrState);

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
  void SwapStates (BosonOnTorusState& state);

  // test if the current state is identical to another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the two states are identical
  bool Equal (BosonOnTorusState& state, int reducedNbrState);

  // test if the current state is different to another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the two states are different
  bool Different (BosonOnTorusState& state, int reducedNbrState);
  
  // test if the current state is greater than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater than the other state
  bool Greater (BosonOnTorusState& state, int reducedNbrState);

  // test if the current state is greater or equal than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater or equal than the other state
  bool GreaterOrEqual (BosonOnTorusState& state, int reducedNbrState);

  // test if the current state is lesser than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater than the other state
  bool Lesser (BosonOnTorusState& state, const int reducedNbrState);

  // test if the current state is lesser or equal than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater or equal than the other state
  bool LesserOrEqual (BosonOnTorusState& state, int reducedNbrState);

  // put the state in a canonical form
  // 
  // tmpState = reference temporary state with the same number of states than the current one
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
  // nbrState = number of state in the current boson state
  // nbrTranslation = reference where the number of translations used to obtain the canonical form will be stored
  // nbrTranslationStep = step to used between two translations
  void PutInCanonicalForm(BosonOnTorusState& tmpState, const int& reducedNbrStat, const int& nbrStateRemainder,
			  const int& nbrState, int& nbrTranslation, const int& nbrTranslationStep);

  // put the state in a canonical form and find how many translations are needed to obtain the same state
  // 
  // tmpState = reference temporary state with the same number of states than the current one
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
  // nbrState = number of state in the current boson state
  // nbrTranslation = reference where the number of translations used to obtain the canonical form will be stored
  // nbrTranslationStep = step to used between two translations
  // nbrTranslationToIdentity = reference on the number of translation needed to obtain the same state
  void PutInCanonicalForm(BosonOnTorusState& tmpState, const int& reducedNbrStat, const int& nbrStateRemainder,
			  const int& nbrState, int& nbrTranslation, const int& nbrTranslationStep, int& nbrTranslationToIdentity);

  // get the number of translation to obtain the same state
  // 
  // tmpState = reference temporary state with the same number of states than the current one
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
  // nbrState = number of state in the current boson state
  // nbrTranslationStep = step to used between two translations
  // return value = number of translation to obtain the same state
  int GetStateSymmetry(BosonOnTorusState& tmpState, const int& reducedNbrState, const int& nbrStateRemainder, 
		       const int& nbrState, const int& nbrTranslationStep);

  // shift a state to the left
  //
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
  // nbrTranslation = magnitude of the translation to apply
  void LeftShiftState(const int& reducedNbrState, const int& nbrStateRemainder, const int& nbrTranslation);

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

inline BosonOnTorusState& BosonOnTorusState::Resize(const int& reducedNbrState)
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

inline BosonOnTorusState& BosonOnTorusState::Assign(BosonOnTorusState& state, const int& reducedNbrState)
{
  for (int i = 0; i <= reducedNbrState; ++i)
    this->StateDescription[i] = state.StateDescription[i];
  return *this;
}

// assign a state to the current one and after undefined it (array tranfert)
//
// state = reference on the state to assign
// return value = reference on  the current state

inline BosonOnTorusState& BosonOnTorusState::TransfertState(BosonOnTorusState& state)
{
  this->StateDescription = state.StateDescription;
  state.StateDescription = 0;
  return *this;
}

// set all ocupations to zero
//
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1
// return value = reference on  the current state

inline BosonOnTorusState& BosonOnTorusState::EmptyState(const int& reducedNbrState)
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

inline unsigned long BosonOnTorusState::GetHashKey (const int& reducedNbrState, const unsigned long& keyMask)
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

inline int BosonOnTorusState::GetHighestIndex (const int& reducedNbrState)
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

inline void BosonOnTorusState::SetOccupation (const int& stateIndex, const unsigned long& nbrBosons)
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

inline unsigned long BosonOnTorusState::GetOccupation (const int& stateIndex)
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

inline unsigned long BosonOnTorusState::IncrementOccupation (const int& stateIndex, unsigned long& tmpShift1, int& tmpShift2)
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

inline unsigned long BosonOnTorusState::DecrementOccupation (const int& stateIndex, unsigned long& tmpShift1, int& tmpShift2)
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

inline bool BosonOnTorusState::TestAndDecrementOccupation (const int& stateIndex, int& tmpShift1, int& tmpShift2)
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

inline void BosonOnTorusState::SwapStates (BosonOnTorusState& state)
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

inline bool BosonOnTorusState::Equal (BosonOnTorusState& state, int reducedNbrState)
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

inline bool BosonOnTorusState::Different (BosonOnTorusState& state, int reducedNbrState)
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

inline bool BosonOnTorusState::Greater (BosonOnTorusState& state, int reducedNbrState)
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

inline bool BosonOnTorusState::GreaterOrEqual (BosonOnTorusState& state, int reducedNbrState)
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

inline bool BosonOnTorusState::Lesser (BosonOnTorusState& state, int reducedNbrState)
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

inline bool BosonOnTorusState::LesserOrEqual (BosonOnTorusState& state, int reducedNbrState)
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

// put the state in a canonical form
// 
// tmpState = reference temporary state with the same number of states than the current one
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
// nbrState = number of state in the current boson state
// nbrTranslation = reference where the number of translations used to obtain the canonical form will be stored
// nbrTranslationStep = step to used between two translations

inline void BosonOnTorusState::PutInCanonicalForm(BosonOnTorusState& tmpState, const int& reducedNbrState, 
						  const int& nbrStateRemainder, const int& nbrState, 
						  int& nbrTranslation, const int& nbrTranslationStep)
{
  tmpState.Assign(*this, reducedNbrState);
  nbrTranslation = 0; 
  int TmpNbrTranslation = nbrTranslationStep;
  while (TmpNbrTranslation < nbrState)
    {
      tmpState.LeftShiftState(reducedNbrState, nbrStateRemainder, nbrTranslationStep);
      if (tmpState.Lesser(*this, reducedNbrState))
	{
	  nbrTranslation = TmpNbrTranslation;
	  this->Assign(tmpState, reducedNbrState);	  
	}
      TmpNbrTranslation += nbrTranslationStep;
    }
}

// put the state in a canonical form and find how many translations are needed to obtain the same state
// 
// tmpState = reference temporary state with the same number of states than the current one
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
// nbrState = number of state in the current boson state
// nbrTranslation = reference where the number of translations used to obtain the canonical form will be stored
// nbrTranslationStep = step to used between two translations
// nbrTranslationToIdentity = reference on the number of translation needed to obtain the same state

inline void BosonOnTorusState::PutInCanonicalForm(BosonOnTorusState& tmpState, const int& reducedNbrState, 
						  const int& nbrStateRemainder, const int& nbrState, int& nbrTranslation, 
						  const int& nbrTranslationStep, int& nbrTranslationToIdentity)
{
  tmpState.Assign(*this, reducedNbrState);
  nbrTranslation = 0; 
  nbrTranslationToIdentity = nbrTranslationStep;
  tmpState.LeftShiftState(reducedNbrState, nbrStateRemainder, nbrTranslationStep);
  while ((nbrTranslationToIdentity < nbrState) && (tmpState.Different(*this, reducedNbrState)))
    {
      if (tmpState.Lesser(*this, reducedNbrState))
	{
	  nbrTranslation = nbrTranslationToIdentity;
	  this->Assign(tmpState, reducedNbrState);	  
	}
      tmpState.LeftShiftState(reducedNbrState, nbrStateRemainder, nbrTranslationStep);
      nbrTranslationToIdentity += nbrTranslationStep;
    }
}

// get the number of translation to obtain the same state
// 
// tmpState = reference temporary state with the same number of states than the current one
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
// nbrState = number of state in the current boson state
// nbrTranslationStep = step to used between two translations
// return value = number of translation to obtain the same state

inline int BosonOnTorusState::GetStateSymmetry(BosonOnTorusState& tmpState, const int& reducedNbrState, 
						const int& nbrStateRemainder, const int& nbrState, 
						const int& nbrTranslationStep)
{
  tmpState.Assign(*this, reducedNbrState);
  int TmpNbrTranslation = nbrTranslationStep;
  tmpState.LeftShiftState(reducedNbrState, nbrStateRemainder, nbrTranslationStep);
  while ((TmpNbrTranslation < nbrState) && (tmpState.Different(*this, reducedNbrState))) 
    {
      tmpState.LeftShiftState(reducedNbrState, nbrStateRemainder, nbrTranslationStep);
      TmpNbrTranslation += nbrTranslationStep;
    }
  return TmpNbrTranslation;
}

// shift a state to the left
//
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
// nbrTranslation = magnitude of the translation to apply

inline void BosonOnTorusState::LeftShiftState(const int& reducedNbrState, const int& nbrStateRemainder, 
					      const int& nbrTranslation)
{
#ifdef __64_BITS__
  unsigned long Remainder;
  // WARNING : this part of the code is not speed optimal but it doesn't need too much temporary variables and is fast enough
  // if the translation is lower than 15 (true for most of the system size reachable with a 64 bits architecture)
  if (nbrStateRemainder != 8)
    for (int j = 0; j < (nbrTranslation >> 3); ++j)
      {
	Remainder = this->StateDescription[0];
	for (int i = 0; i < reducedNbrState; ++i)
	  {
	    this->StateDescription[i] = this->StateDescription[i + 1];
	  }      
	this->StateDescription[reducedNbrState - 1] |= Remainder << (nbrStateRemainder << 3);
	this->StateDescription[reducedNbrState] = Remainder >> ((8 - nbrStateRemainder) << 3);
      }
  else
    for (int j = 0; j < (nbrTranslation >> 3); ++j)
      {
	Remainder = this->StateDescription[0];
	for (int i = 0; i < reducedNbrState; ++i)
	  {
	    this->StateDescription[i] = this->StateDescription[i + 1];
	  }      
	this->StateDescription[reducedNbrState] = Remainder;
      }
  switch (nbrTranslation & 0x7)
    {
    case 1:
      {	
	Remainder = this->StateDescription[0] & ((unsigned long) 0xff);
	for (int i = 0; i < reducedNbrState; ++i)
	  {
	    this->StateDescription[i] >>= 8;
	    this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xff)) << 56;
	  }
	this->StateDescription[reducedNbrState] >>= 8;
	this->StateDescription[reducedNbrState] |= Remainder << ((nbrStateRemainder - 1) << 3);
      }
      break;
    case 2:
      {	
	if (nbrStateRemainder >= 2)
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffff);
	    for (int i = 0; i < reducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 16;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffff)) << 48;
	      }
	    this->StateDescription[reducedNbrState] >>= 16;
	    this->StateDescription[reducedNbrState] |= Remainder << ((nbrStateRemainder - 2) << 3);
	  }
	else
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffff);
	    int ReducedNbrState = reducedNbrState - 1;
	    for (int i = 0; i < ReducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 16;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffff)) << 16;
	      }
	    this->StateDescription[ReducedNbrState] >>= 16;
	    this->StateDescription[ReducedNbrState] |= this->StateDescription[reducedNbrState] << 48;
	    this->StateDescription[ReducedNbrState] |= (Remainder & ((unsigned long) 0xff)) << 56;
	    this->StateDescription[reducedNbrState] = Remainder >> 8;
	  }
      }
      break;
    case 3:
      {	
	if (nbrStateRemainder >= 3)
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffffff);
	    for (int i = 0; i < reducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 24;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffffff)) << 40;
	      }
	    this->StateDescription[reducedNbrState] >>= 24;
	    this->StateDescription[reducedNbrState] |= Remainder << ((nbrStateRemainder - 3) << 3);
	  }
	else
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffffff);
	    int ReducedNbrState = reducedNbrState - 1;
	    for (int i = 0; i < ReducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 24;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffffff)) << 40;
	      }
	    this->StateDescription[ReducedNbrState] >>= 24;
	    this->StateDescription[ReducedNbrState] |= this->StateDescription[reducedNbrState] << 40;
	    this->StateDescription[ReducedNbrState] |= (Remainder << ((nbrStateRemainder + 5) << 3));
	    this->StateDescription[reducedNbrState] = Remainder >> ((3 - nbrStateRemainder) << 3);
	  }
      }
      break;
    case 4:
      {	
	if (nbrStateRemainder >= 4)
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffffffff);
	    for (int i = 0; i < reducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 32;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffffffff)) << 32;
	      }
	    this->StateDescription[reducedNbrState] >>= 32;
	    this->StateDescription[reducedNbrState] |= Remainder << ((nbrStateRemainder - 4) << 3);
	  }
	else
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffffffff);
	    int ReducedNbrState = reducedNbrState - 1;
	    for (int i = 0; i < ReducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 32;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffffffff)) << 32;
	      }
	    this->StateDescription[ReducedNbrState] >>= 32;
	    this->StateDescription[ReducedNbrState] |= this->StateDescription[reducedNbrState] << 32;
	    this->StateDescription[ReducedNbrState] |= (Remainder << ((nbrStateRemainder + 4) << 3));
	    this->StateDescription[reducedNbrState] = Remainder >> ((4 - nbrStateRemainder) << 3);
	  }
      }
      break;
    case 5:
      {	
	if (nbrStateRemainder >= 5)
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffffffffff);
	    for (int i = 0; i < reducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 40;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffffffffff)) << 24;
	      }
	    this->StateDescription[reducedNbrState] >>= 40;
	    this->StateDescription[reducedNbrState] |= Remainder << ((nbrStateRemainder - 5) << 3);
	  }
	else
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffffffffff);
	    int ReducedNbrState = reducedNbrState - 1;
	    for (int i = 0; i < ReducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 40;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffffffffff)) << 24;
	      }
	    this->StateDescription[ReducedNbrState] >>= 40;
	    this->StateDescription[ReducedNbrState] |= this->StateDescription[reducedNbrState] << 24;
	    this->StateDescription[ReducedNbrState] |= (Remainder << ((nbrStateRemainder + 3) << 3));
	    this->StateDescription[reducedNbrState] = Remainder >> ((5 - nbrStateRemainder) << 3);
	  }
      }
      break;
    case 6:
      {	
	if (nbrStateRemainder >= 6)
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffffffffffff);
	    for (int i = 0; i < reducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 48;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffffffffffff)) << 16;
	      }
	    this->StateDescription[reducedNbrState] >>= 48;
	    this->StateDescription[reducedNbrState] |= Remainder << ((nbrStateRemainder - 6) << 3);
	  }
	else
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffffffffffff);
	    int ReducedNbrState = reducedNbrState - 1;
	    for (int i = 0; i < ReducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 48;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffffffffffff)) << 16;
	      }
	    this->StateDescription[ReducedNbrState] >>= 48;
	    this->StateDescription[ReducedNbrState] |= this->StateDescription[reducedNbrState] << 16;
	    this->StateDescription[ReducedNbrState] |= (Remainder << ((nbrStateRemainder + 2) << 3));
	    this->StateDescription[reducedNbrState] = Remainder >> ((6 - nbrStateRemainder) << 3);
	  }
      }
      break;
    case 7:
      {	
	if (nbrStateRemainder >= 7)
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffffffffffffff);
	    for (int i = 0; i < reducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 56;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffffffffffffff)) << 8;
	      }
	    this->StateDescription[reducedNbrState] >>= 56;
	    this->StateDescription[reducedNbrState] |= Remainder << ((nbrStateRemainder - 7) << 3);
	  }
	else
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffffffffffffff);
	    int ReducedNbrState = reducedNbrState - 1;
	    for (int i = 0; i < ReducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 56;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffffffffffffff)) << 8;
	      }
	    this->StateDescription[ReducedNbrState] >>= 56;
	    this->StateDescription[ReducedNbrState] |= this->StateDescription[reducedNbrState] << 8;
	    this->StateDescription[ReducedNbrState] |= (Remainder << ((nbrStateRemainder + 1) << 3));
	    this->StateDescription[reducedNbrState] = Remainder >> ((7 - nbrStateRemainder) << 3);
	  }
      }
      break;
    }

#else
  unsigned long Remainder;
  // WARNING : this part of the code is not speed optimal but it doesn't need too much temporary variables and is fast enough
  // if the translation is lower than 8 (true for most of the system size reachable with a 32 bits architecture)
  if (nbrStateRemainder != 4)
    for (int j = 0; j < (nbrTranslation >> 2); ++j)
      {
	Remainder = this->StateDescription[0];
	for (int i = 0; i < reducedNbrState; ++i)
	  {
	    this->StateDescription[i] = this->StateDescription[i + 1];
	  }      
	this->StateDescription[reducedNbrState - 1] |= Remainder << (nbrStateRemainder << 3);
	this->StateDescription[reducedNbrState] = Remainder >> ((4 - nbrStateRemainder) << 3);
      }
  else
    for (int j = 0; j < (nbrTranslation >> 2); ++j)
      {
	Remainder = this->StateDescription[0];
	for (int i = 0; i < reducedNbrState; ++i)
	  {
	    this->StateDescription[i] = this->StateDescription[i + 1];
	  }      
	this->StateDescription[reducedNbrState] = Remainder;
      }
  switch (nbrTranslation & 0x3)
    {
    case 1:
      {	
	Remainder = this->StateDescription[0] & ((unsigned long) 0xff);
	for (int i = 0; i < reducedNbrState; ++i)
	  {
	    this->StateDescription[i] >>= 8;
	    this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xff)) << 24;
	  }
	this->StateDescription[reducedNbrState] >>= 8;
	this->StateDescription[reducedNbrState] |= Remainder << ((nbrStateRemainder - 1) << 3);
      }
      break;
    case 2:
      {	
	if (nbrStateRemainder >= 2)
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffff);
	    for (int i = 0; i < reducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 16;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffff)) << 16;
	      }
	    this->StateDescription[reducedNbrState] >>= 16;
	    this->StateDescription[reducedNbrState] |= Remainder << ((nbrStateRemainder - 2) << 3);
	  }
	else
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffff);
	    int ReducedNbrState = reducedNbrState - 1;
	    for (int i = 0; i < ReducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 16;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffff)) << 16;
	      }
	    this->StateDescription[ReducedNbrState] >>= 16;
	    this->StateDescription[ReducedNbrState] |= this->StateDescription[reducedNbrState] << 16;
	    this->StateDescription[ReducedNbrState] |= (Remainder & ((unsigned long) 0xff)) << 24;
	    this->StateDescription[reducedNbrState] = Remainder >> 8;
	  }
      }
      break;
    case 3:
      {	
	if (nbrStateRemainder >= 3)
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffffff);
	    for (int i = 0; i < reducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 24;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffffff)) << 8;
	      }
	    this->StateDescription[reducedNbrState] >>= 24;
	    this->StateDescription[reducedNbrState] |= Remainder << ((nbrStateRemainder - 3) << 3);
	  }
	else
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffffff);
	    int ReducedNbrState = reducedNbrState - 1;
	    for (int i = 0; i < ReducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 24;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffffff)) << 8;
	      }
	    this->StateDescription[ReducedNbrState] >>= 24;
	    this->StateDescription[ReducedNbrState] |= this->StateDescription[reducedNbrState] << 8;
	    this->StateDescription[ReducedNbrState] |= (Remainder << ((nbrStateRemainder + 1) << 3));
	    this->StateDescription[reducedNbrState] = Remainder >> ((3 - nbrStateRemainder) << 3);
	  }
      }
      break;
    }
#endif
}

// print a given state
//
// str = reference on current output stream 
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
// return value = reference on current output stream 

inline ostream& BosonOnTorusState::PrintState (ostream& str, const int& reducedNbrState, const int& nbrStateRemainder)
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
