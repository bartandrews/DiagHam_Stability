////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//          class of hard-core boson with a single quantum number             //
//                                                                            //
//                        last modification : 11/02/2008                      //
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
#include "HilbertSpace/HardCoreBosonLong.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "QuantumNumber/NumberParticleQuantumNumber.h"

#include <bitset>
#include <iostream>
#include <cstdlib>

using std::bitset;
using std::cout;
using std::endl;

// default constructor
HardCoreBosonLong::HardCoreBosonLong()
{
  this->NbrBosons=0;
  this->HilbertSpaceDimension=0;
}


// constructor
// nbrBosons = number of bosons in system
// nbrStates = number of states available
// memory = amount of memory granted for precalculations
HardCoreBosonLong::HardCoreBosonLong(int nbrBosons, int nbrStates, unsigned long memory)
{
#ifdef __64_BITS__  
  if (nbrStates>127)    
#else
  if (nbrStates>63)    
#endif
    {
      cout<<"HardCoreBosonLong: Cannot represent the "<<nbrStates<<" states requested in a single word"<<endl;
      exit(1);
    }
  this->NbrBosons=nbrBosons;
  this->NbrStates=nbrStates;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(nbrBosons,nbrStates);
  this->StateDescription=new ULONGLONG[this->HilbertSpaceDimension];
  this->StateHighestBit=new int[this->HilbertSpaceDimension];
  this->GenerateStates(nbrBosons,nbrStates);
  this->TargetSpace=this;
  this->GenerateLookUpTable(memory);
  this->Flag.Initialize();
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
#ifdef __DEBUG__
  // for (int i=0; i<HilbertSpaceDimension; ++i)
//     {
//       this->PrintState(cout,i);
//       cout << endl;
//     }
  unsigned long UsedMemory = 0;  
  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
  // review here!
  //  UsedMemory += this->NbrLzValue * sizeof(int);
  // UsedMemory += this->NbrLzValue * ((unsigned long) this->LookUpTableMemorySize) * sizeof(int);
  // UsedMemory +=  (1 << this->MaximumSignLookUp) * sizeof(double);
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
#endif
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension<<endl;
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
HardCoreBosonLong::HardCoreBosonLong(const HardCoreBosonLong& bosons)
{
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->NbrStates = bosons.NbrStates;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateHighestBit = bosons.StateHighestBit;
  this->Flag = bosons.Flag;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}


// virtual destructor
//

HardCoreBosonLong::~HardCoreBosonLong ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrStates; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
}

// assignment (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space
HardCoreBosonLong& HardCoreBosonLong::operator = (const HardCoreBosonLong& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrStates; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
  if (this->TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->NbrStates = bosons.NbrStates;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateHighestBit = bosons.StateHighestBit;
  this->Flag = bosons.Flag;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}


// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space
AbstractHilbertSpace* HardCoreBosonLong::Clone()
{
  return new HardCoreBosonLong(*this);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void HardCoreBosonLong::SetTargetSpace(HardCoreBosonLong* targetSpace)
{
  this->TargetSpace=targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int HardCoreBosonLong::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->GetHilbertSpaceDimension();
}

// get the particle statistic 
//
// return value = particle statistic
int HardCoreBosonLong::GetParticleStatistic()
{
  return AbstractQHEParticle::BosonicStatistic;
}

// get information about any additional symmetry of the Hilbert space
//
// return value = symmetry id

int HardCoreBosonLong::GetHilbertSpaceAdditionalSymmetry()
{
  return HardCoreBosonLong::NoSymmetry;
}



// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number
List<AbstractQuantumNumber*> HardCoreBosonLong::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new NumberParticleQuantumNumber(this->NbrBosons);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number
AbstractQuantumNumber* HardCoreBosonLong::GetQuantumNumber (int index)
{
  return new NumberParticleQuantumNumber(this->NbrBosons);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace
AbstractHilbertSpace* HardCoreBosonLong::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int HardCoreBosonLong::AdAdAA (int index, int m1, int m2, int n1, int n2)
{
  int StateHighestBit = this->StateHighestBit[index];
  ULONGLONG State = this->StateDescription[index];
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & (( (ULONGLONG) (0x1ul)) << ((ULONGLONG) n1) )) == (ULONGLONG)0x0ul ) 
      || ((State & (( (ULONGLONG) (0x1ul)) << ((ULONGLONG) n2) )) == (ULONGLONG)0x0ul) || (n1 == n2) || (m1 == m2))
    {
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  ULONGLONG TmpState = State;
  // perform annihilation operators
  TmpState &= ~(( (ULONGLONG) (0x1ul)) << ((ULONGLONG) n2) );
  if (NewHighestBit == n2)
    while ((TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  TmpState &= ~(( (ULONGLONG) (0x1ul)) << ((ULONGLONG) n1) );
  if (NewHighestBit == n1)
    while ((TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  // create particle at m2
  if ((TmpState & (((ULONGLONG) (0x1ul)) << ((ULONGLONG) m2)))!= ((ULONGLONG) 0x0ul  ))
    {
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewHighestBit)
    {
      NewHighestBit = m2;
    }
  TmpState |= (((ULONGLONG) (0x1ul)) << ((ULONGLONG) m2));
  // create particle at m1
  if ((TmpState & (((ULONGLONG) (0x1ul)) << ((ULONGLONG) m1)))!= ((ULONGLONG) 0x0ul  ))
    {
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewHighestBit)
    {
      NewHighestBit = m1;
    }
  TmpState |= (((ULONGLONG) (0x1ul)) << ((ULONGLONG) m1));
  return this->TargetSpace->FindStateIndex(TmpState, NewHighestBit);
}


// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double HardCoreBosonLong::AA (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];

  if (((ProdATemporaryState & (((ULONGLONG) (0x1ul)) << ((ULONGLONG) n1))) ==  ((ULONGLONG) 0x0ul)) || ((ProdATemporaryState & (((ULONGLONG) (0x1ul)) << ((ULONGLONG) n2))) == ((ULONGLONG) 0x0ul)) || (n1 == n2))
    return 0.0;

  this->ProdAHighestBit = this->StateHighestBit[index];

  this->ProdATemporaryState &= ~(((ULONGLONG) (0x1ul)) << ((ULONGLONG) n2));
  
  this->ProdATemporaryState &= ~(((ULONGLONG) (0x1ul)) << ((ULONGLONG) n1));

  while ((this->ProdATemporaryState >> this->ProdAHighestBit) == 0)
    --this->ProdAHighestBit;
  return 1.0;
}


// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int HardCoreBosonLong::AdAd (int m1, int m2)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  if ((TmpState & (((ULONGLONG) (0x1ul)) << ((ULONGLONG) m2))) !=  ((ULONGLONG) 0x0ul))
    {
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = this->ProdAHighestBit;
  if (m2 > NewHighestBit)
    NewHighestBit = m2;
  TmpState |= (((ULONGLONG) (0x1ul)) << ((ULONGLONG) m2));
  if ((TmpState & (((ULONGLONG) (0x1ul)) << ((ULONGLONG) m1)))!= ((ULONGLONG) 0x0ul))
    {
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewHighestBit)
    NewHighestBit = m1;
  TmpState |= (((ULONGLONG) (0x1ul)) << ((ULONGLONG) m1));
  return this->TargetSpace->FindStateIndex(TmpState, NewHighestBit);
}


// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int HardCoreBosonLong::AdA (int index, int m, int n)
{
  if (m!=n)
    {
      int StateHighestBit = this->StateHighestBit[index];
      ULONGLONG State = this->StateDescription[index];
      if ((n > StateHighestBit) || ((State & (((ULONGLONG) (0x1ul)) << ((ULONGLONG) n ))) == ( (ULONGLONG)0x0ul) ))
	{
	  return this->HilbertSpaceDimension;
	}
      int NewHighestBit = StateHighestBit;
      ULONGLONG TmpState = State;
      // perform annihilation operators
      TmpState &= ~(((ULONGLONG) (0x1ul)) << ((ULONGLONG) n) );
      if (NewHighestBit == n)
	while ((TmpState >> NewHighestBit) == 0)
	  --NewHighestBit;
      // create particle at m
      if ((TmpState & (((ULONGLONG) (0x1ul)) << ((ULONGLONG) m) )) != ((ULONGLONG) 0x0ul))
	{
	  return this->HilbertSpaceDimension;
	}
      if (m > NewHighestBit)
	{
	  NewHighestBit = m;
	}
      TmpState |= (((ULONGLONG) (0x1ul)) << ((ULONGLONG) m));
      return this->TargetSpace->FindStateIndex(TmpState, NewHighestBit);
    }
  else          
    return this->AdA (index, m);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m
int HardCoreBosonLong::AdA (int index, int m)
{
  int StateHighestBit = this->StateHighestBit[index];
  ULONGLONG State = this->StateDescription[index];
  if ((m > StateHighestBit) || ((State & (((ULONGLONG) (0x1ul)) << ((ULONGLONG) m))) ==  ((ULONGLONG) 0x0ul)))
    {
      return this->HilbertSpaceDimension;
    }
  else
    {
      return index;
    }
}

  // evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex HardCoreBosonLong::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
{
  return this->EvaluateWaveFunction(state, position, basis, 0, this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// return value = wave function evaluated at the given location

Complex HardCoreBosonLong::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, int nextCoordinates)
{
  return this->EvaluateWaveFunctionWithTimeCoherence(state, position, basis, nextCoordinates, 0, 
						     this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location
Complex HardCoreBosonLong::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
						int firstComponent, int nbrComponent)
{
  return Complex(0.0, 0.0);
}


// evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex HardCoreBosonLong::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, 
								 int nextCoordinates, int firstComponent, 
								 int nbrComponent)
{
  return Complex(0.0, 0.0);
}
                                                                                                                        
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void HardCoreBosonLong::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& HardCoreBosonLong::PrintState (ostream& Str, int state)
{
  ULONGLONG TmpState = this->StateDescription[state];
  for (int i = 0; i < this->NbrStates; ++i)
  {
    if( ((TmpState >> ((ULONGLONG) i))& ((ULONGLONG) 0x1ul)) ==  ((ULONGLONG) 0x1ul))
      Str<< "1 ";
    else
      Str<< "0 ";
  }
    Str << " position = " << this->FindStateIndex(TmpState, this->StateHighestBit[state]);
  if (state !=  this->FindStateIndex(TmpState, this->StateHighestBit[state]))
    Str << " error! ";
  return Str;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int HardCoreBosonLong::FindStateIndex(ULONGLONG stateDescription, int highestBit)
{
  // bitset <32> b = stateDescription;
//   cout << "in FindStateIndex: desc=" << b<<", highest bit=" << highestBit << endl;
  long PosMax = stateDescription >> this->LookUpTableShift[highestBit];
  long PosMin = this->LookUpTable[highestBit][PosMax];
  PosMax = this->LookUpTable[highestBit][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  ULONGLONG CurrentState = this->StateDescription[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    return PosMin;
}


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// nbrStates = number of elementary states available
//
int HardCoreBosonLong::EvaluateHilbertSpaceDimension(int nbrBosons,int nbrStates)
{  
  FactorialCoefficient F;
  F.SetToOne();
  F.FactorialMultiply(nbrStates);
  F.FactorialDivide(nbrStates-nbrBosons);
  F.FactorialDivide(nbrBosons);
  return (int)F.GetIntegerValue();
}

// generate many-body states without any additional symmetries
// nbrBosons = number of bosons
// nbrStates = number of elementary states available
//
// (using routines from UnsignedIntegerTools)
//
int HardCoreBosonLong::GenerateStates(int nbrBosons,int nbrStates)
{
  int countdown=this->HilbertSpaceDimension-1;
  int presentHighestBit=nbrBosons-1;    
  this->StateHighestBit[countdown] = presentHighestBit;
  this->StateDescription[countdown--]=smallestOne(nbrBosons);
  while (countdown>-1)
    {      
      this->StateDescription[countdown]=nextone(this->StateDescription[countdown+1]);
      if (this->StateDescription[countdown] & (0x1ul << (presentHighestBit+1)))
	++presentHighestBit;
      this->StateHighestBit[countdown]=presentHighestBit;
      --countdown;
    }
  if (this->StateDescription[0]!=biggestOne(nbrBosons, nbrStates))
    {
      bitset <64> b=this->StateDescription[0];
      cout << "Error in GenerateStates: Last state is: " << b << endl;
      exit(1);
    }
  return 0;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void HardCoreBosonLong::GenerateLookUpTable(unsigned long memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * this->NbrStates);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrStates)
    this->MaximumLookUpShift = this->NbrStates;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrStates];
  this->LookUpTableShift = new int [this->NbrStates];
  for (int i = 0; i < this->NbrStates; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentHighestBit = this->StateHighestBit[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
  if (CurrentHighestBit < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentHighestBit] = 0;
  else
    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentHighestBit];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentHighestBit != this->StateHighestBit[i])
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
	  /*	  for (unsigned long j = 0; j <= this->LookUpTableMemorySize; ++j)
	    cout << TmpLookUpTable[j] << " ";
	    cout << endl << "-------------------------------------------" << endl;*/
 	  CurrentHighestBit = this->StateHighestBit[i];
	  TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
	  if (CurrentHighestBit < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentHighestBit] = 0;
	  else
	    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentHighestBit];
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  CurrentLookUpTableValue = this->LookUpTableMemorySize;
	  while (CurrentLookUpTableValue > TmpLookUpTableValue)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      while (CurrentLookUpTableValue > TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  while (CurrentLookUpTableValue > 0)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[0] = this->HilbertSpaceDimension - 1;
  /*  for (unsigned long j = 0; j <= this->LookUpTableMemorySize; ++j)
    cout << TmpLookUpTable[j] << " ";
    cout << endl << "-------------------------------------------" << endl;*/
}

