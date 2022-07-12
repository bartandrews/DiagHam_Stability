////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2009 Nicolas Regnault                  //
//                         class author: Gunnar Möller                        //
//                                                                            //
//                                                                            //
//                       class of abstract spin collection                    //
//                                                                            //
//                        last modification : 18/04/2001                      //
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


#include "config.h"
#include "GenericSUNSpinCollection.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/StringTools.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"

#include <cstdlib>
#include <limits.h>
#include <iostream>
using std::cout;
using std::endl;


// verbosity flag
//#define VERBOSE

// testing flag;
//#define TESTING


// bitset only needed in testing mode
#ifdef TESTING
#include <bitset>
using std::bitset;
#endif

// create Hilbert-space for a collection of spins
// levelN = level of SU(N) symmetry
// nbrSpins = number of spins in collection
// cartanQuantumNumbers = quantum numbers of the first N-1 operators of the Cartan-algebra
// memory = amount of memory granted for precalculations
GenericSUNSpinCollection::GenericSUNSpinCollection(int levelN, int nbrSpins, int *cartanQuantumNumbers, unsigned long memory)
{
  this->LevelN = levelN;
  this->NbrSpins = nbrSpins;
  this->Flag.Initialize();

  if (LevelN>MAXLEVEL)
    {
      cout << "The class GenericSUNSpinCollection is limited to SU(N) levels N<="<<MAXLEVEL<<endl;
      exit(-1);
    }
#ifdef __64_BITS__
  if (NbrSpins>64/BITS)
    {
      cout << "The class GenericSUNSpinCollection is limited to "<<64/BITS<<" spins in 64 bit mode"<<endl;
      exit(-1);
    }
#else
  if (NbrSpins>32/BITS)
    {
      cout << "The class GenericSUNSpinCollection is limited to "<<32/BITS<<" spins in 32 bit mode"<<endl;
      exit(-1);
    }
#endif

  this->CartanQuantumNumbers = new int[MAXLEVEL];
  this->CartanQuantumNumbers[LevelN-1]=this->NbrSpins;
  for (int i=0; i<LevelN-1; ++i)
    {
      this->CartanQuantumNumbers[i] = cartanQuantumNumbers[i];
      this->CartanQuantumNumbers[LevelN-1]-=cartanQuantumNumbers[i];
    }
  if (this->CartanQuantumNumbers[LevelN-1]<0)
    {
      cout << "The first N-1 Cartan quantum numbers have to add to less than NbrSpins."<<endl;
      exit(-1);
    }
  for (int i=LevelN; i<MAXLEVEL; ++i) this->CartanQuantumNumbers[i]=0;
  
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrSpins, this->LevelN,
								    this->CartanQuantumNumbers);

  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrSpins-1, this->LevelN, 0l, CartanQuantumNumbers[0],
						       CartanQuantumNumbers[1], CartanQuantumNumbers[2],
						       CartanQuantumNumbers[3], CartanQuantumNumbers[4],
						       CartanQuantumNumbers[5], CartanQuantumNumbers[6],
						       CartanQuantumNumbers[7]);
  if (TmpHilbertSpaceDimension != this->HilbertSpaceDimension)
    {
      cout << "Mismatch in State-count and State Generation in GenericSUNSpinCollection!" << endl;
      exit(1);
    }

//   for (int i=0; i<HilbertSpaceDimension; ++i)
//     PrintStateOnly(cout, i)<<endl;  


  this->GenerateLookUpTable(memory);

#ifdef VERBOSE
  for (int i=0; i<HilbertSpaceDimension; ++i)
    PrintState(cout, i)<<endl;
#endif

#ifdef __DEBUG__
  long UsedMemory = 0;
  UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long));
  cout << "memory requested for Hilbert space = ";
  PrintMemorySize(cout, UsedMemory)<<endl;
  UsedMemory = this->LookUpTableDepth*sizeof(int);
  for (int i=1; i<this->LookUpTableDepth; ++i) UsedMemory *= this->LevelN;
  cout << "memory requested for lookup table = ";
  PrintMemorySize(cout, UsedMemory)<<endl;
#endif
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}


// copy constructor
GenericSUNSpinCollection::GenericSUNSpinCollection (GenericSUNSpinCollection &collection)
{
  this->HilbertSpaceDimension = collection.HilbertSpaceDimension;
  this->Flag = collection.Flag;
  this->LevelN = collection.LevelN;
  this->NbrSpins = collection.NbrSpins;
  this->CartanQuantumNumbers = collection.CartanQuantumNumbers;
  this->StateDescription = collection.StateDescription;
  this->LookUpTableDepth = collection.LookUpTableDepth;
  this->LookUpTableShift = collection.LookUpTableShift;
  this->LookUpTable = collection.LookUpTable;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
GenericSUNSpinCollection::~GenericSUNSpinCollection ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->CartanQuantumNumbers;
      delete[] this->LookUpTable;
    }
}


// assignment operator
GenericSUNSpinCollection& GenericSUNSpinCollection::operator = (GenericSUNSpinCollection &collection)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->CartanQuantumNumbers;
      // delete[] this->LookUpTable;
    }
  this->HilbertSpaceDimension = collection.HilbertSpaceDimension;
  this->Flag = collection.Flag;
  this->LevelN = collection.LevelN;
  this->NbrSpins = collection.NbrSpins;
  this->CartanQuantumNumbers = collection.CartanQuantumNumbers;
  this->StateDescription = collection.StateDescription;
  this->LookUpTableDepth = collection.LookUpTableDepth;
  this->LookUpTableShift = collection.LookUpTableShift;
  this->LookUpTable = collection.LookUpTable;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}


// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* GenericSUNSpinCollection::Clone()
{
  return new GenericSUNSpinCollection(*this);
}


// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number
List<AbstractQuantumNumber*> GenericSUNSpinCollection::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  for (int i=0; i<LevelN; ++i)
  TmpList += new SzQuantumNumber (this->CartanQuantumNumbers[i]);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number
AbstractQuantumNumber* GenericSUNSpinCollection::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->CartanQuantumNumbers[index]);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace
AbstractHilbertSpace* GenericSUNSpinCollection::ExtractSubspace (AbstractQuantumNumber& q, 
								 SubspaceSpaceConverter& converter)
{
  return 0;
}


// permutation operator of two spins
int GenericSUNSpinCollection::SpinPermutation(int index, int s1, int s2)
{
  unsigned long State = this->StateDescription[index];
//   bitset<32> b = State;
//   cout << "P["<<s1<<", "<<s2<<", "<<b<<"] = ";
  int S1Pos=BITS*s1;
  unsigned Spin1 = ((State>>(s1*BITS))&MASK);
  State &= (~0ul ^ (MASK<<S1Pos));
  int S2Pos=BITS*s2;
  unsigned Spin2 = ((State>>(s2*BITS))&MASK);
  State &= (~0ul ^ (MASK<<S2Pos));
  State |= (Spin2<<S1Pos);
  State |= (Spin1<<S2Pos);
//   b=State;
//   cout << b << endl;
  return this->FindStateIndex(State);
}

// cyclic permutation of three spins
// index = index of state to perform on
// s1 = index of spin 1
// s2 = index of spin 2
// s3 = index of spin 3
int GenericSUNSpinCollection::CyclicSpinPermutation(int index, int s1, int s2, int s3)
{
  unsigned long State = this->StateDescription[index];
//   bitset<32> b = State;
//   cout << "P["<<s1<<", "<<s2<<", "<<b<<"] = ";
  int LastPos=BITS*s3;
  unsigned BufferedSpin = ((State>>(s3*BITS))&MASK);
  State &= (~0ul ^ (MASK<<LastPos));

  int CurrentPos=BITS*s2;
  unsigned CurrentSpin = ((State>>(s2*BITS))&MASK);
  State &= (~0ul ^ (MASK<<CurrentPos));
  State |= (CurrentSpin<<LastPos);

  LastPos=CurrentPos;
  CurrentPos=BITS*s1;
  CurrentSpin = ((State>>(s1*BITS))&MASK);
  State &= (~0ul ^ (MASK<<CurrentPos));
  State |= (CurrentSpin<<LastPos);

  State |= (BufferedSpin<<CurrentPos);
//   b=State;
//   cout << b << endl;
  return this->FindStateIndex(State);
}



// cyclic permutation of k spins
// index = index of state to perform on
// numS = number of spins to permute
// si = indices of spins
// return = index of final state
int GenericSUNSpinCollection::CyclicSpinPermutation(int index, int numS, int *si)
{
  unsigned long State = this->StateDescription[index];
#ifdef TESTING
  bitset<32> b = State;
  cout << "P["<<si[0];
  for (int i=1; i<numS; ++i)
    cout <<", "<<si[i];
  cout <<"; "<< b<< "] = ?"<<endl;
#endif
  int *spinPt=&(si[numS-1]);
  int LastPos=BITS* (*spinPt);
  unsigned BufferedSpin = ((State>>(*spinPt*BITS))&MASK);
  State &= (~0ul ^ (MASK<<LastPos));
#ifdef TESTING
  cout << "Buffered spin ("<<si[numS-1] <<"): "<< BufferedSpin <<endl;
  b=State; 
#endif
  int CurrentPos=0;
  unsigned CurrentSpin;
  for (int n=numS-2; n>=0; --n)
    {
      --spinPt;
      CurrentPos = BITS* *spinPt;
      CurrentSpin = ((State>>(*spinPt *BITS))&MASK);
#ifdef TESTING
      cout << "now moving spin "<<*spinPt<<endl;
      cout << "value: "<<CurrentSpin<<endl;
#endif
      State &= (~0ul ^ (MASK<<CurrentPos));
      State |= (CurrentSpin<<LastPos);
      LastPos = CurrentPos;
#ifdef TESTING
      b=State; 
      cout << "moved:        "<<b<<endl;
#endif
    }

  State |= (BufferedSpin<<CurrentPos);
  
#ifdef TESTING
  b=State;
  cout <<"final state:  "<< b << endl;
#endif
  return this->FindStateIndex(State);
}

// get diagonal terms for an S*S interaction (counting instances for connections with same prefactor)
// index = number of state to be considered
int GenericSUNSpinCollection::S2DiagonalElements(int index)
{
  return this->HilbertSpaceDimension;
}

// get off-diagonal terms for an S*S interaction (counting instances for connections with same prefactor)
// index = number of state to be considered
// spin = number of spin associated with first spin operator
// targetIndices = states connected to (size: N-1)
//
void GenericSUNSpinCollection::S2OffDiagonalElements(int index, int spin, int *targetIndices)
{
}




// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int GenericSUNSpinCollection::FindStateIndex(unsigned long stateDescription)
{  
  unsigned Key = this->GenerateKey(stateDescription >> this->LookUpTableShift);
  //cout << stateDescription << " key="<<Key<<" LKUP="<<this->LookUpTable[Key] << endl;
  int PosMax = this->LookUpTable[Key]-1;
  int PosMin = this->LookUpTable[Key+1];
  int PosMid = (PosMax + PosMin) >> 1;
  //cout << "PosMax="<<PosMax<<" PosMin="<<PosMin<<endl;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ( (PosMin != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMin = PosMid;
	}
      else
	{
	  PosMax = PosMid;
	} 
      PosMid = (PosMax + PosMin) >> 1;
      CurrentState = this->StateDescription[PosMid];
      //      cout << "loop PosMax = "<<PosMax<<", PosMid = "<<PosMid<<", PosMin = "<<PosMin<<endl;
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    return PosMax;
}


// evaluate the Hilbert space dDimension
// nbrSpins = number of spins
// levelN = SU(N) level
// cartanQuantumNumbers = full set of eigenvalues of the Cartan operators
//
long GenericSUNSpinCollection::EvaluateHilbertSpaceDimension(int nbrSpins, int levelN, int *cartanQuantumNumbers)
{
  FactorialCoefficient Dimension(1);
  Dimension.FactorialMultiply(nbrSpins);
  for (int i=0; i<levelN; ++i)
    {
      Dimension.FactorialDivide(cartanQuantumNumbers[i]);
    }
  return Dimension.GetIntegerValue();
}

// generate all states corresponding to the constraints
//
// nbrSpin = number of spin to be considered next
// levelN = highest remaining SU(N)-level
// pos = position in StateDescription array where to store states
// t_i = number of spins remaining with <T_i> = t_i
// return value = position from which new states have to be stored
long GenericSUNSpinCollection::GenerateStates(int nbrSpin, int levelN, long pos,
					      int t0, int t1, int t2, int t3, int t4, int t5, int t6, int t7)
{
  if (nbrSpin < 0) 
    return pos;
  long TmpPos;
  unsigned long Mask;
  if (nbrSpin==0)
    {      
      this->StateDescription[pos] = (unsigned long) (t1+0x2l*t2+0x3l*t3+0x4l*t4+0x5l*t5+0x6l*t6+0x7l*t7);
      return (pos + 1l);
    }
  if (t7 > 0)
    {
      TmpPos = this->GenerateStates(nbrSpin - 1, levelN, pos, t0, t1, t2, t3, t4, t5, t6, t7-1);
      Mask = 0x7ul << ( nbrSpin * BITS );
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if (t6 > 0)
    {
      TmpPos = this->GenerateStates(nbrSpin - 1, levelN, pos, t0, t1, t2, t3, t4, t5, t6-1, t7);
      Mask = 0x6ul << ( nbrSpin * BITS );
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if (t5 > 0)
    {
      TmpPos = this->GenerateStates(nbrSpin - 1, levelN, pos, t0, t1, t2, t3, t4, t5-1, t6, t7);
      Mask = 0x5ul << ( nbrSpin * BITS );
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if (t4 > 0)
    {
      TmpPos = this->GenerateStates(nbrSpin - 1, levelN, pos, t0, t1, t2, t3, t4-1, t5, t6, t7);
      Mask = 0x4ul << ( nbrSpin * BITS );
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if (t3 > 0)
    {
      TmpPos = this->GenerateStates(nbrSpin - 1, levelN, pos, t0, t1, t2, t3-1, t4, t5, t6, t7);
      Mask = 0x3ul << ( nbrSpin * BITS );
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if (t2 > 0)
    {
      TmpPos = this->GenerateStates(nbrSpin - 1, levelN, pos, t0, t1, t2-1, t3, t4, t5, t6, t7);
      Mask = 0x2ul << ( nbrSpin * BITS );
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if (t1 > 0)
    {
      TmpPos = this->GenerateStates(nbrSpin - 1, levelN, pos, t0, t1-1, t2, t3, t4, t5, t6, t7);
      Mask = 0x1ul << ( nbrSpin * BITS );
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if (t0 > 0)
    {
      TmpPos = this->GenerateStates(nbrSpin - 1, levelN, pos, t0-1, t1, t2, t3, t4, t5, t6, t7);
      pos = TmpPos;
    }
  return pos;
}


// generate look-up table associated to current Hilbert space
// 
// memeory = memory size that can be allocated for the look-up table
void GenericSUNSpinCollection::GenerateLookUpTable(unsigned long memory)
{
  unsigned long MaxNbrIntegers = memory/sizeof(int);
  unsigned long TmpMemory=this->LevelN;
  unsigned MaxKey=this->LevelN;
  this->LookUpTableDepth=1;
  while ((TmpMemory*LevelN<MaxNbrIntegers)&&(LookUpTableDepth<this->NbrSpins)&&(MaxKey<UINT_MAX/LevelN))
    {
      TmpMemory*=this->LevelN;
      MaxKey*=this->LevelN;
      MaxKey+=this->LevelN;
      this->LookUpTableDepth++;
    }
  this->LookUpTableShift=BITS*(this->NbrSpins-this->LookUpTableDepth);
  this->LookUpTable = new int[TmpMemory+1];
//   cout << "LookUpTableDepth="<<LookUpTableDepth<<endl;
//   cout << "LookUpTableShift="<<LookUpTableShift<<endl;
//   cout << "LookUpTableSize="<<TmpMemory<<endl;
  unsigned long CurrentIndexedBits=this->StateDescription[0]>>LookUpTableShift;
  unsigned CurrentKey=this->GenerateKey(CurrentIndexedBits);
  unsigned NewKey;
  for (unsigned t=CurrentKey+1; t<TmpMemory+1; ++t)
    {
      this->LookUpTable[t]=0;
    }
  //  cout << "CurrentKey at start: "<<CurrentKey<<endl;
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentIndexedBits != this->StateDescription[i]>>LookUpTableShift)
	{
	  CurrentIndexedBits = this->StateDescription[i]>>LookUpTableShift;
	  NewKey = this->GenerateKey(CurrentIndexedBits);
	  for (unsigned t=CurrentKey; t>NewKey; --t)
	    this->LookUpTable[t]=i;
	  CurrentKey=NewKey;
	}
    }
  for (unsigned t=CurrentKey; t>0; --t)
    this->LookUpTable[t]=this->HilbertSpaceDimension;
  this->LookUpTable[0]=this->HilbertSpaceDimension;

//   for (int i=0; i<TmpMemory+1;++i)
//     cout << "LookUpTable["<<i<<"]="<<LookUpTable[i]<<endl;
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 
ostream& GenericSUNSpinCollection::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  Str << " ";
  for (int i=this->NbrSpins-1; i>=0; --i)
    Str << ((TmpState>>(i*BITS))&MASK)<<" ";
  Str << " position = " << this->FindStateIndex(TmpState);
  if (state !=  this->FindStateIndex(TmpState))
    Str << " error! ";
  return Str;
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 
ostream& GenericSUNSpinCollection::PrintStateOnly (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  Str << state<<": "<<TmpState<< " = ";
  for (int i=this->NbrSpins-1; i>=0; --i)
    Str << ((TmpState>>(i*BITS))&MASK)<<" ";
  return Str;
}
