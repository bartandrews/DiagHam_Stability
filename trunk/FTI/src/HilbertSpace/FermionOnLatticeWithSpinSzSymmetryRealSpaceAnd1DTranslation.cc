////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//       in real space with translation invariance in one direction           //
//                          and the Sz<->-Sz symmetry                         //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                        last modification : 30/07/2014                      //
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
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"
#include "GeneralTools/StringTools.h"

#include <math.h>
#include <cstdlib>
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

FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation::FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation ()
{
  this->NbrFermions = 0;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalSpin = 0;
  this->MaxMomentum = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSite = 0;
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = 0;
  this->XMomentum = 0; 
  this->YMomentum = 0;
  this->SzParitySign = 1.0;
  this->StateShift = 2 * (this->MaxMomentum / this->MomentumModulo);
  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = ((unsigned long) 1);
  this->MaximumSignLookUp = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->StateHighestBit = 0;  
  this->LargeHilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// momnetum = momentum sector
// periodicity = periodicity with respect to site numbering 
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation::FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation (int nbrFermions, int nbrSite, int momentum, int periodicity, 
															  bool minusSzParity, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = this->NbrSite / periodicity;
  cout << "MomentumModulo=" << MomentumModulo<<endl;
  this->XMomentum = momentum % this->MomentumModulo;
  this->YMomentum = 0;
  this->StateShift = 2 * periodicity;
  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;
  this->SzParitySign = 1.0;
  if (minusSzParity == true)
    this->SzParitySign = -1.0;

  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions);
  cout << "intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->GenerateSignLookUpTable();
      this->LargeHilbertSpaceDimension  = this->GenerateStates(true, true);
      this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
      cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->HilbertSpaceDimension > 0)
	this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
      UsedMemory = this->NbrMomentum * sizeof(int);
      UsedMemory += this->NbrMomentum * this->LookUpTableMemorySize * sizeof(int);
      cout << "memory requested for lookup table = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
#endif
    }
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation::FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation(const FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;

      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      delete[] this->SignLookUpTable;
      delete[] this->NbrParticleLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->NbrSite = fermions.NbrSite;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->NbrFermionStates = fermions.NbrFermionStates;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;
  this->SzParitySign = fermions.SzParitySign;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation::~FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation& FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation::operator = (const FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;

      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      delete[] this->SignLookUpTable;
      delete[] this->NbrParticleLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->NbrSite = fermions.NbrSite;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->NbrFermionStates = fermions.NbrFermionStates;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;
  this->SzParitySign = fermions.SzParitySign;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation::Clone()
{
  return new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      Tmp = (TmpState >> (i << 1)) & 03ul;
      switch (Tmp)
	{
	case 0x0ul:
	  Str << "0 ";
	  break;
	case 0x1ul:
	  Str << "d ";
	  break;
	case 0x2ul:
	  Str << "u ";
	  break;
	case 0x3ul:
	  Str << "X ";
	  break;
	}
    }
  Str << " " << this->FindStateIndex(this->StateDescription[state], this->StateHighestBit[state]);
  if (this->FindStateIndex(this->StateDescription[state], this->StateHighestBit[state]) != state)
    {
      Str << "  error";
    }
  Str << " orb. size=" << this->NbrStateInOrbit[state];
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation::RawGenerateStates(int nbrFermions, int currentSite, long pos)
{
  if (nbrFermions == 0)
    {
      this->StateDescription[pos] = 0x0ul;	  
      return (pos + 1l);
    }
  if ((currentSite < 0) || (nbrFermions < 0))
    return pos;
  if (nbrFermions == 1)
    {
      for (int j = currentSite; j >= 0; --j)
	{
	  this->StateDescription[pos] = 0x2ul << (j << 1);
	  ++pos;
	  this->StateDescription[pos] = 0x1ul << (j << 1);
	  ++pos;
	}
      return pos;
    }
  long TmpPos = this->RawGenerateStates(nbrFermions - 2, currentSite - 1, pos);
  unsigned long Mask = 0x3ul << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->RawGenerateStates(nbrFermions - 1, currentSite - 1, pos);
  Mask = 0x2ul << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
   TmpPos = this->RawGenerateStates(nbrFermions - 1, currentSite - 1, pos);
   Mask = 0x1ul << (currentSite << 1);
   for (; pos < TmpPos; ++pos)
     this->StateDescription[pos] |= Mask;
  return this->RawGenerateStates(nbrFermions, currentSite - 1, pos);
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// return value = Hilbert space dimension

long FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation::EvaluateHilbertSpaceDimension(int nbrFermions)
{
  BinomialCoefficients binomials(2*this->NbrSite);
  long dimension = binomials(2*this->NbrSite, this->NbrFermions);
  return dimension;
}


// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation::FindStateIndex(char* stateDescription)
{
  char** TmpDescription;
  if (SplitLine(stateDescription, TmpDescription, ' ') != this->MaxMomentum)
    return -1;
  unsigned long TmpState = 0x0ul;
  int TmpNbrParticles = 0;
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      if (TmpDescription[i][0] == 'u')
	{
	  TmpState |= 0x2ul << (2 * i);
	  ++TmpNbrParticles;	  
	}
      else
	{
	  if (TmpDescription[i][0] == 'd')
	    {
	      TmpState |= 0x1ul << (2 * i);
	      ++TmpNbrParticles;	  
	    }
	  else
	    {
	      if (TmpDescription[i][0] == 'X')
		{
		  TmpState |= 0x3ul << (2 * i);
		  TmpNbrParticles += 2;	  
		}
	      else
		{
		  if (TmpDescription[i][0] != '0')
		    {
		      return -1;
		    }
		}
	    }
	}
      delete[] TmpDescription[i];
    }
  delete[] TmpDescription;
  if (TmpNbrParticles != this->NbrFermions)
    return -1;
  int TmpLzMax = 2 * this->MaxMomentum + 1;
  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  return this->FindStateIndex(TmpState, TmpLzMax);
}

// find state index
//
// stateDescription = unsigned longeger describing the state
// maxMomentum = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation::FindStateIndex(unsigned long stateDescription, int maxMomentum)
{
  if ((stateDescription > this->StateDescription[0]) || 
      (stateDescription < this->StateDescription[this->HilbertSpaceDimension - 1]))
    return this->HilbertSpaceDimension;

//   long PosMax = stateDescription >> this->LookUpTableShift[maxMomentum];
//   long PosMin = this->LookUpTable[maxMomentum][PosMax];
//   PosMax = this->LookUpTable[maxMomentum][PosMax + 1];

  long PosMax = 0l;
  long PosMin = this->LargeHilbertSpaceDimension - 1l;
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
  //  cout << hex << CurrentState << " " <<  stateDescription << " " << this->StateDescription[0] << " " << this->StateDescription[this->HilbertSpaceDimension - 1] << dec << " " << PosMid << " " << PosMin << " " << PosMax << endl;
  if (CurrentState == stateDescription)
    return PosMid;
  else
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMin;
//   if ((stateDescription < this->StateDescription[0]) || (stateDescription > this->StateDescription[this->HilbertSpaceDimension - 1]))
//     return this->HilbertSpaceDimension;
//   long PosMax = stateDescription >> this->LookUpTableShift[maxMomentum];
//   long PosMin = this->LookUpTable[maxMomentum][PosMax];
//   PosMax = this->LookUpTable[maxMomentum][PosMax + 1];
//   long PosMid = (PosMin + PosMax) >> 1;
//   unsigned long CurrentState = this->StateDescription[PosMid];
//   while ((PosMin != PosMid) && (CurrentState != stateDescription))
//     {
//       if (CurrentState > stateDescription)
// 	{
// 	  PosMax = PosMid;
// 	}
//       else
// 	{
// 	  PosMin = PosMid;
// 	} 
//       PosMid = (PosMin + PosMax) >> 1;
//       CurrentState = this->StateDescription[PosMid];
//     }
//   if (CurrentState == stateDescription)
//     return PosMid;
//   else
//     if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
//       return this->HilbertSpaceDimension;
//     else
//       return PosMax;
}

