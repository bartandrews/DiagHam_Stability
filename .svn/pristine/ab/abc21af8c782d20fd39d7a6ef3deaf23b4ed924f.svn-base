////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                     class author : Cecile Repellin                         //
//                                                                            //
//                        class of fermions the CP2                           //
//                           with up to 128 sites                             //
//                                                                            //
//                        last modification : 26/02/2013                      //
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
#include "GeneralTools/Endian.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "FermionOnCP2Long.h"


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

FermionOnCP2Long::FermionOnCP2Long ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// p = number of flux quanta (determines an irreducible representation of SU(3), along with q=0 (LLL))
// totalJz = total value of jz
// totalKz = total value of kz
// memory = amount of memory granted for precalculations

FermionOnCP2Long::FermionOnCP2Long (int nbrFermions, int nbrFluxQuanta, int totalTz, int totalY, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->TotalTz = totalTz;
  this->TotalY = totalY;
  this->TotalR = (this->TotalY + 3*this->TotalTz + 2*this->NbrFermions*this->NbrFluxQuanta)/6;
  this->TotalS = (this->TotalY - 3*this->TotalTz + 2*this->NbrFermions*this->NbrFluxQuanta)/6;
  this->NbrLzValue = (this->NbrFluxQuanta + 1)*(this->NbrFluxQuanta + 2)/2;
  this->LzMax = NbrLzValue - 1;  
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrFluxQuanta, this->NbrFluxQuanta, this->NbrFluxQuanta, 0, 0);
  this->quantumNumberTz = new int [this->NbrLzValue];
  this->quantumNumberY = new int [this->NbrLzValue];
  this->quantumNumberR = new int [this->NbrLzValue];
  this->quantumNumberS = new int [this->NbrLzValue];
  cout << "dim = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->GetQuantumNumbersFromLinearizedIndex(this->quantumNumberTz, this->quantumNumberY, this->quantumNumberR, this->quantumNumberS); 
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
      this->StateLzMax = new int [this->HilbertSpaceDimension];  
      this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrFluxQuanta, this->NbrFluxQuanta, this->NbrFluxQuanta, 0, 0, 0l);
      int TmpLzMax = this->LzMax;
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  while ((this->StateDescription[i] >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  this->StateLzMax[i] = TmpLzMax;
	}
      this->GenerateLookUpTable(memory);      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(ULONGLONG) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
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

FermionOnCP2Long::FermionOnCP2Long(const FermionOnCP2Long& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->TotalTz = fermions.TotalTz;
  this->TotalY = fermions.TotalY;
  this->TotalR = fermions.TotalR;
  this->TotalS = fermions.TotalS;
  this->quantumNumberTz = fermions.quantumNumberTz;
  this->quantumNumberY = fermions.quantumNumberY;  
  this->quantumNumberR = fermions.quantumNumberR;
  this->quantumNumberS = fermions.quantumNumberS; 
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

FermionOnCP2Long::~FermionOnCP2Long ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnCP2Long& FermionOnCP2Long::operator = (const FermionOnCP2Long& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
    }
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalTz = fermions.TotalTz;
  this->TotalY = fermions.TotalY;
  this->TotalR = fermions.TotalR;
  this->TotalS = fermions.TotalS;
  this->quantumNumberTz = fermions.quantumNumberTz;
  this->quantumNumberY = fermions.quantumNumberY;  
  this->quantumNumberR = fermions.quantumNumberR;
  this->quantumNumberS = fermions.quantumNumberS; 
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnCP2Long::Clone()
{
  return new FermionOnCP2Long(*this);
}


// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured
bool FermionOnCP2Long::WriteHilbertSpace (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->HilbertSpaceDimension);
  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);
  WriteLittleEndian(File, this->NbrFermions);
  WriteLittleEndian(File, this->TotalTz);
  WriteLittleEndian(File, this->TotalY);
  if (this->HilbertSpaceDimension != 0)
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateDescription[i]);
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateLzMax[i]);
    }
  else
    {
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateDescription[i]);
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateLzMax[i]);
    }
  File.close();
  return true;
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnCP2Long::PrintState (ostream& Str, int state)
{
  ULONGLONG TmpState = this->StateDescription[state];
  ULONGLONG Tmp;
  Str << "[";
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      Tmp = (TmpState >> i);
      int TmpTz = this->quantumNumberTz[i];
      int TmpY = this->quantumNumberY[i];
//       cout << TmpTz << ", " << TmpY << "; " << (Tmp & 0x1l) << endl;
      if ((Tmp & (ULONGLONG) 0x1l) != (ULONGLONG) 0ul)
	Str << "(" << TmpTz << "," << TmpY << ")";
    }
  Str << "]";
  return Str;
}


// generate all states corresponding to the constraints
// 
// stateDescription = array that gives each state description
// nbrBosons = number of bosons
// currentJ = current value of j for a single particle
// currentJz = current value of jz for a single particle
// currentKz = current value of kz for a single particle
// currentTotalTz = current total value of Jz
// currentTotalY = current total value of Kz
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnCP2Long::GenerateStates(int nbrFermions, int currentTz, int currentTzMax, int currentY, int currentTotalTz, int currentTotalY, long pos)
{
  if (nbrFermions < 0)
    return pos;
  if (currentTotalTz + currentTzMax*nbrFermions < this->TotalTz)
    return pos;
  if (currentTotalY + currentY*nbrFermions < this->TotalY)
    return pos;
  
  if (currentTz < -currentTzMax)
   {
     --currentTzMax;
     currentTz = currentTzMax;
     currentY = currentY - 3;
   }
    
  
  if (nbrFermions == 0)
    {
      if ((currentTotalTz == this->TotalTz) && (currentTotalY == this->TotalY))
	{
	  this->StateDescription[pos] = (ULONGLONG) 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }

  if (currentY < -2*this->NbrFluxQuanta)
    return pos;
  
  long TmpPos = this->GenerateStates(nbrFermions - 1, currentTz - 2, currentTzMax, currentY, currentTotalTz + currentTz, currentTotalY + currentY, pos);
  ULONGLONG Mask = (((ULONGLONG) 0x1ul << 1) - (ULONGLONG) 0x1ul) << this->GetLinearizedIndex(currentTz, currentY, 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, currentTz - 2, currentTzMax, currentY, currentTotalTz, currentTotalY, pos);
};



// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentJz = current value of jz for a single particle
// currentKz = current value of kz for a single particle
// currentTotalTz = current total value of Jz
// currentTotalY = current total value of Kz
// return value = Hilbert space dimension

long FermionOnCP2Long::EvaluateHilbertSpaceDimension(int nbrBosons, int currentTz, int currentTzMax, int currentY, int currentTotalTz, int currentTotalY)
{
  if (nbrBosons < 0)
    return 0l;
  if (currentTotalTz + currentTzMax*nbrBosons < this->TotalTz)
    return 0l;
  if (currentTotalY + currentY*nbrBosons < this->TotalY)
    return 0l;
  
  if (currentTz < -currentTzMax)
   {
     --currentTzMax;
     currentTz = currentTzMax;
     currentY = currentY - 3;
   }
    
  if (nbrBosons == 0)
    {
      if ((currentTotalTz == this->TotalTz) && (currentTotalY == this->TotalY))
      {
	return 1l;
      }
      else	
	return 0l;
    }
    
  if (currentY < -2*this->NbrFluxQuanta)
    return 0l;
  
  long Count = 0;
  
  Count += this->EvaluateHilbertSpaceDimension(nbrBosons - 1, currentTz - 2, currentTzMax, currentY, currentTotalTz + currentTz, currentTotalY + currentY);
  Count += this->EvaluateHilbertSpaceDimension(nbrBosons, currentTz - 2, currentTzMax, currentY, currentTotalTz, currentTotalY);
  return Count;
}

// find state index from an array
//
// stateDescription = array describing the state (stored as tz1,y1,tz2,y2,...)
// return value = corresponding index, -1 if an error occured

int FermionOnCP2Long::FindStateIndexFromArray(int* stateDescription)
{
  ULONGLONG TmpState = (ULONGLONG) 0x0ul;
  for (int i = 0; i < this->NbrFermions; ++i)
    TmpState |= (ULONGLONG) 0x1ul << (this->GetLinearizedIndex(stateDescription[i << 1], stateDescription[1 + (i << 1)], this->NbrFermions));
  int TmpLzMax = this->LzMax;
  while ((TmpState >> TmpLzMax) == (ULONGLONG) 0x0ul)
    --TmpLzMax;
  return this->FindStateIndex(TmpState, TmpLzMax);
}
