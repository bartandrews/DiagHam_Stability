////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of fermions on sphere without fixed total Lz             //
//                                                                            //
//                        last modification : 19/04/2010                      //
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
#include "HilbertSpace/FermionOnSphereFull.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/Endian.h"
#include "MathTools/BinomialCoefficients.h"

#include <math.h>
#include <fstream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
//

FermionOnSphereFull::FermionOnSphereFull()
{
  this->HilbertSpaceDimension = 0;
  this->LargeHilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// lzMax = twice the maximum Lz value reached by a fermion
// memory = amount of memory granted for precalculations
// referenceState = array that describes the reference state to start from

FermionOnSphereFull::FermionOnSphereFull (int nbrFermions, int lzMax, unsigned long memory)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->TotalLz = 0;

#ifdef __64_BITS__
  this->InvertShift = 32 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 16 - ((this->LzMax + 1 ) >> 1);
#endif
  if ((this->LzMax & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();

  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateLzMax = new int [this->LargeHilbertSpaceDimension];
  this->TotalLzValues = new int [this->LargeHilbertSpaceDimension];
  if (this->NbrFermions > 0)
    {
      this->GenerateStates(this->NbrFermions, this->LzMax, this->LzMax, 0);
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->StateDescription[i];
	  int TmpTotalLz = 0;
	  for (int j = this->StateLzMax[i]; j >=0; --j)
	    TmpTotalLz += j * ((int) ((TmpState >> j) & 0x1ul));
	  this->TotalLzValues[i] = (TmpTotalLz << 1) - (this->NbrFermions * this->LzMax) ;
	}
    }
  else
    {
      this->StateDescription[0] = 0x0ul; 
      this->StateLzMax[0] = 0;
    }
  this->GenerateLookUpTable(memory);

/*
  double coeff;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Clone();
  for (int k=0; k<this->HilbertSpaceDimension; ++k)
  {
     cout<<"k= "<<k<<" "; this->PrintState(cout, k); cout<<endl;
     for (int i=0; i<=this->LzMax; ++i)
       for (int j=0; j<=this->LzMax; ++j)
         {
           int Index = TmpParticles->AdA(k,i,j,coeff);
           if (Index < this->HilbertSpaceDimension)
             cout<<"i= "<<i<<" j="<<j<<" "<<Index<<endl;
         }
  }
  delete TmpParticles;
  exit(1);
*/

#ifdef __DEBUG__
  unsigned long UsedMemory = 0l;
  UsedMemory += ((unsigned long) this->LargeHilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
  UsedMemory +=  (1 << this->MaximumSignLookUp) * sizeof(double);
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
#endif
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereFull::FermionOnSphereFull(const FermionOnSphereFull& fermions)
{
  if (fermions.TargetSpace != ((FermionOnSphere*)&fermions))
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->TotalLzValues = fermions.TotalLzValues;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
}

// destructor
//

FermionOnSphereFull::~FermionOnSphereFull ()
{
}


// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereFull& FermionOnSphereFull::operator = (const FermionOnSphereFull& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
  if (fermions.TargetSpace != ((FermionOnSphere*)&fermions))
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->TotalLzValues = fermions.TotalLzValues;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->InitializeWaveFunctionEvaluation();
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereFull::Clone()
{
  return new FermionOnSphereFull(*this);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereFull::AdA (int index, int m)
{
  return 0.0;
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereFull::AdA (int index, int m, int n, double& coefficient)
{
  int StateLzMax = this->StateLzMax[index];
  unsigned long State = this->StateDescription[index];
  if ((n > StateLzMax) || ((State & (((unsigned long) (0x1)) << n)) == 0))
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  int NewLzMax = StateLzMax;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  TmpState &= ~(((unsigned long) (0x1)) << n);
  if ((TmpState != 0x0ul))
    {
      while ((TmpState >> NewLzMax) == 0)
	--NewLzMax;
    }
  else
    NewLzMax = 0;
  if ((TmpState & (((unsigned long) (0x1)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  if (m > NewLzMax)
    {
      NewLzMax = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m);
  int NewIndex = this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
  //cout<<"AdA "; this->PrintState(cout,index);cout<<" m= "<<m<<" n="<<n<<" NewLzMax="<<NewLzMax<<" NewIndex= "<<NewIndex<<endl;
  return NewIndex;
}

// carefully test whether state is in Hilbert-space and find corresponding state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
// return value = corresponding index, or dimension of space, if not found
int FermionOnSphereFull::FindStateIndex(unsigned long stateDescription, int lzmax)
{
//      for (int i=0; i<HilbertSpaceDimension; ++i)
//	if (this->StateDescription[i] == stateDescription)
//          return i;
//      return this->HilbertSpaceDimension;

   int PosMax = this->HilbertSpaceDimension - 1;
   int PosMin = 0;

   while (PosMin <= PosMax) 
     {
       int PosMid = (PosMin + PosMax) >> 1;  // compute mid point.
       unsigned long CurrentState = this->StateDescription[PosMid];
       if (CurrentState > stateDescription) 
           PosMin = PosMid + 1;  // repeat search in top half.
       else if (CurrentState < stateDescription) 
           PosMax = PosMid - 1; // repeat search in bottom half.
       else
           return PosMid;     // found it. return position /////
   }
   return this->HilbertSpaceDimension;    // failed to find key

/*
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
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

*/

}

// generate all states (i.e. all possible skew symmetric polynomials with fixed Lz)
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// currentLzMax = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereFull::GenerateStates(int nbrFermions, int lzMax, int currentLzMax, long pos)
{
  if ((nbrFermions == 0) || (currentLzMax < (nbrFermions - 1)))
    return pos;
  if (nbrFermions == 1)
    {
      int currentLzMaxOld = currentLzMax;
      for ( ; currentLzMax >= 0; --currentLzMax)
	{
	  this->StateDescription[pos] = 0x1ul << currentLzMax;
	  this->StateLzMax[pos] = lzMax;
	  ++pos;
	}
      return pos;
    }
  int ReducedCurrentLzMax = currentLzMax - 1;
  long TmpPos = this->GenerateStates(nbrFermions - 1, lzMax, ReducedCurrentLzMax, pos);
  unsigned long Mask = 0x1ul << currentLzMax;
  for (long i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  if (lzMax == currentLzMax)
    return this->GenerateStates(nbrFermions, ReducedCurrentLzMax, ReducedCurrentLzMax, TmpPos);
  else
    return this->GenerateStates(nbrFermions, lzMax, ReducedCurrentLzMax, TmpPos);
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// return value = Hilbert space dimension

long FermionOnSphereFull::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax)
{
  BinomialCoefficients TmpCoefficients (lzMax + 1);
  return TmpCoefficients(lzMax + 1, nbrFermions);
}


// convert the vector with a given Lz to the full space (all Lz components)
// inputState = input vector
// inputSpace = input Hilbert space with given Lz
// return value = vector in the full Hilbert space

void FermionOnSphereFull::ConvertToAllLz (ComplexVector& inputState, ParticleOnSphere* inputSpace, ComplexVector& outputState)
{
  FermionOnSphere* InputSpace = (FermionOnSphere*) inputSpace;

  int Index;
  for (int i = 0; i < InputSpace->GetHilbertSpaceDimension(); i++)
   {
     Index = this->FindStateIndex(InputSpace->StateDescription[i], InputSpace->StateLzMax[i]);
     if (Index < this->HilbertSpaceDimension)
      {
        outputState[Index] = inputState[i];
      }
     else
      {
       cout<<"Component not found!"<<endl;
       exit(-1);
      }
   }
}
