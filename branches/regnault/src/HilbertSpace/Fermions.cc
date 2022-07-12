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
//                        last modification : 04/04/2001                      //
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


#include "HilbertSpace/Fermions.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/NumberParticleQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include <iostream.h>

#define MAX_HILBERTSPACE_DIMENSION 32768


// find gcd
//
// m = first number
// n = second number
// return value = gcd of m and n
int FindGCD (int m, int n);

// recurrence for gcd
//
// m = first number
// n = second number
// return value = gcd of m and n
int gcd (int m, int n);

// default constructor
//

Fermions::Fermions () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableSize = 0;
  this->LookUpTableMask = 0;
  this->LookUpPosition = 0;
  this->HilbertSpaceDimension = 0;
  this->FermionDescription = 0;
  this->Parity = 0;
  this->NbrSite = 0;
  this->Sz = 0;
  this->NbrFermion = 0;
  this->FixedQuantumNumberFlag = false;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz nor on
// total number of fermions 
//
// nbrSite = number of site
// memorySize = memory size in bytes allowed for look-up table

Fermions::Fermions (int nbrSite, int memorySize) 
{
  this->Flag.Initialize();
  this->NbrSite = nbrSite;
  this->FixedQuantumNumberFlag = false;
  this->HilbertSpaceDimension = 3;
  for (int i = 1; i < this->NbrSite; i++)
    this->HilbertSpaceDimension *= 3;

  this->LookUpPosition = 0;
  this->LookUpTableSize = 4;
  memorySize >>= 2;
  this->LookUpTableMask = 0xfffffffd;
  while ((this->LookUpPosition <= this->NbrSite) && (memorySize >=  4))
    {
      this->LookUpTableMask <<= 2;
      this->LookUpTableSize <<= 2;
      memorySize >>= 2;
      this->LookUpPosition++;
    }
  this->LookUpTableMask = ~this->LookUpTableMask;
  this->LookUpTable = new int [this->LookUpTableSize];

  this->FermionDescription = new unsigned long [this->HilbertSpaceDimension];
  this->Parity = new unsigned long [this->HilbertSpaceDimension];
  this->GenerateStates (0, 0, 0x00000000);
}

// constructor for complete Hilbert space corresponding to a given total spin projection Sz and a
// given number of fermions
//
// nbrSite = number of site
// n = number of fermions
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table

Fermions::Fermions (int nbrSite, int n, int sz, int memorySize) 
{
  this->Flag.Initialize();
  this->NbrSite = nbrSite;
  this->NbrFermion = n;
  this->Sz = sz;
  if (((this->Sz >= 0) && ((this->NbrFermion % 2) == (this->Sz % 2))) || 
      ((this->Sz < 0) && ((this->NbrFermion % 2) == (-(this->Sz % 2)))))
    this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrSite, 
								      (this->NbrFermion + this->Sz) / 2, 
								      (this->NbrFermion - this->Sz) / 2);
  else
    this->HilbertSpaceDimension = 0;

  cout << this->HilbertSpaceDimension << endl;
  this->FixedQuantumNumberFlag = true;

  this->LookUpPosition = 0;
  this->LookUpTableSize = 4;
  memorySize >>= 1;
  this->LookUpTableMask = 0xfffffffd;
  while ((this->LookUpPosition < this->NbrSite) && (memorySize >=  4))
    {
      this->LookUpTableMask <<= 2;
      this->LookUpTableSize <<= 2;
      memorySize >>= 2;
      this->LookUpPosition++;
    }
  this->LookUpTableMask = ~this->LookUpTableMask;
  this->LookUpTable = new int [this->LookUpTableSize];

  this->FermionDescription = new unsigned long [this->HilbertSpaceDimension];
  this->Parity = new unsigned long [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = this->GenerateStates (0, 0, 0x00000000, 0, 0);
}

// constructor from pre-constructed datas
//
// hilbertSpaceDimension = Hilbert space dimension
// fermionDescription = array describing states
// parity = array describing parity
// nbrSite = number of spin sites
// n = number of fermions
// sz = twice the value of total Sz component
// fixedQuantumNumberFlag = true if hilbert space is restricted to a given quantum number
// lookUpTable = look-up table
// lookUpTableSize = look-Up table size
// lookUpTablePosition = last position described by the look-Up table
// lookUpTableMask = look-Up table mask  

Fermions::Fermions (int hilbertSpaceDimension, unsigned long* fermionDescription, unsigned long* parity,
		    int nbrSite, int n, int sz, bool fixedQuantumNumberFlag, int* lookUpTable, 
		    int lookUpTableSize, int lookUpPosition, unsigned long lookUpTableMask)
{
  this->Flag.Initialize();
  this->LookUpTable = lookUpTable;
  this->LookUpTableMask = lookUpTableMask;
  this->LookUpPosition = lookUpPosition;
  this->LookUpTableSize = lookUpTableSize;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->FermionDescription = fermionDescription;
  this->Parity = parity;
  this->Sz = sz;
  this->NbrSite = nbrSite;
  this->NbrFermion = n;
  this->FixedQuantumNumberFlag = fixedQuantumNumberFlag;
  this->NbrSite = nbrSite;
}
  
// copy constructor (without duplicating datas)
//
// fermions = reference on fermions to copy

Fermions::Fermions (const Fermions& fermions) 
{
  this->Flag = fermions.Flag;
  if (fermions.NbrSite != 0)
    {
      this->NbrSite = fermions.NbrSite;
      this->NbrFermion = fermions.NbrFermion;
      this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
      this->LookUpTable = fermions.LookUpTable;
      this->LookUpTableMask = fermions.LookUpTableMask;
      this->LookUpPosition = fermions.LookUpPosition;
      this->LookUpTableSize = fermions.LookUpTableSize;
      this->FermionDescription = fermions.FermionDescription;
      this->Parity = fermions.Parity;
      this->Sz = fermions.Sz;
      this->FixedQuantumNumberFlag = fermions.FixedQuantumNumberFlag;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableSize = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->HilbertSpaceDimension = 0;
      this->FermionDescription = 0;
      this->Parity = 0;
      this->NbrSite = 0;
      this->NbrFermion = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
    }
}

// destructor
//

Fermions::~Fermions () 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true) && (this->NbrFermion != 0))
    {
      delete[] this->FermionDescription;
      delete[] this->LookUpTable;
      delete[] this->Parity;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on fermions to copy
// return value = reference on current fermions

Fermions& Fermions::operator = (const Fermions& fermions)
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true) && (this->NbrFermion != 0))
    {
      delete[] this->FermionDescription;
      delete[] this->LookUpTable;
      delete[] this->Parity;
    }  
  this->Flag = fermions.Flag;
  if (fermions.NbrFermion != 0)
    {
      this->NbrSite = fermions.NbrSite;
      this->NbrFermion = fermions.NbrFermion;
      this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
      this->LookUpTable = fermions.LookUpTable;
      this->LookUpTableMask = fermions.LookUpTableMask;
      this->LookUpPosition = fermions.LookUpPosition;
      this->LookUpTableSize = fermions.LookUpTableSize;
      this->FermionDescription = fermions.FermionDescription;
      this->Parity = fermions.Parity;
      this->Sz = fermions.Sz;
      this->FixedQuantumNumberFlag = fermions.FixedQuantumNumberFlag;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableSize = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->HilbertSpaceDimension = 0;
      this->FermionDescription = 0;
      this->Parity = 0;
      this->NbrSite = 0;
      this->NbrFermion = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Fermions::Clone()
{
  return new Fermions (*this);
}

// generate all states
//
// statePosition = position for the new states
// sitePosition = site where fermion has to be changed
// currentStateDescription = description of current state
// currentParity = current parity (true if there is an odd number of fermions to 
// the left of the current one)
// return value = number of generated states

int Fermions::GenerateStates(int statePosition,int sitePosition, unsigned long currentStateDescription) 
{
  int NextSitePosition = sitePosition + 1;     
  int NbrGeneratedState = -statePosition;
  sitePosition &= 0x0000000f;
  sitePosition <<= 1;
  if (NextSitePosition != this->NbrSite)
    {
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription);
      
      currentStateDescription |= (0x00000002 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription);
      
      currentStateDescription |= (0x00000001 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription);
    }
  else
    {
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      this->FermionDescription[statePosition] = currentStateDescription;
      this->Parity[statePosition] = this->EvaluateParity(currentStateDescription);
      statePosition++;
      
      currentStateDescription |= (0x00000002 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      this->FermionDescription[statePosition] = (currentStateDescription);
      this->Parity[statePosition] = this->EvaluateParity(currentStateDescription);
      statePosition++;
      
      currentStateDescription |= (0x00000001 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      this->FermionDescription[statePosition] = (currentStateDescription);
      statePosition++;
    }
  NbrGeneratedState += statePosition;
  return NbrGeneratedState;
}

// generate all states corresponding to a given total Sz and a total number of fermions
//
// statePosition = position for the new states
// sitePosition = site where fermion has to be changed
// currentStateDescription = description of current state
// currentN = current number of fermions
// currentSz = total Sz value of current state
// currentParity = current parity (true if there is an odd number of fermions to 
// the left of the current one)
// return value = number of generated states

int Fermions::GenerateStates(int statePosition,int sitePosition, unsigned long currentStateDescription, 
			     int currentN, int currentSz)
{
  int MaxN = (this->NbrSite - sitePosition);
  if ((this->NbrFermion - currentN) > MaxN)
    return 0;    
  int MaxSz = (this->NbrSite - sitePosition);
  if (((currentSz + MaxSz) < this->Sz) || ((currentSz - MaxSz) > this->Sz))
    return 0;
  int NextSitePosition = sitePosition + 1;     
  int NbrGeneratedState = -statePosition;
  sitePosition &= 0x0000000f;
  sitePosition <<= 1;
  if (NextSitePosition != this->NbrSite)
    {
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, 
					    currentStateDescription, currentN, currentSz);
      
      currentStateDescription |= (0x00000002 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, 
					    currentStateDescription, currentN + 1, currentSz - 1);
      
      currentStateDescription |= (0x00000001 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, 
					    currentStateDescription, currentN + 1, currentSz + 1);
    }
  else
    {
      if ((currentSz == this->Sz) && (currentN == this->NbrFermion))
	{
	  if (NextSitePosition == this->LookUpPosition)
	    this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
	  this->FermionDescription[statePosition] = currentStateDescription;
	  this->Parity[statePosition] = this->EvaluateParity(currentStateDescription);
	  statePosition++;
	}
      else
	{
	  currentSz--;
	  currentN += 1;
	  if ((currentSz == this->Sz) && (currentN == this->NbrFermion))
	    {
	      currentStateDescription |= (0x00000002 << sitePosition);
	      if (NextSitePosition == this->LookUpPosition)
		this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
	      this->FermionDescription[statePosition] = currentStateDescription;
	      this->Parity[statePosition] = this->EvaluateParity(currentStateDescription);
	      statePosition++;
	    }
	  else
	    {
	      currentSz += 2;
	      if ((currentSz == this->Sz) && (currentN == this->NbrFermion))
		{
		  currentStateDescription |= (0x00000003 << sitePosition);
		  if (NextSitePosition == this->LookUpPosition)
		    this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
		  this->FermionDescription[statePosition] = currentStateDescription;
		  this->Parity[statePosition] = this->EvaluateParity(currentStateDescription);
		  statePosition++;
		}
	    }
	}
    }
  NbrGeneratedState += statePosition;
  return NbrGeneratedState;
}

// return Hilbert space dimension
//
// return value = Hilbert space dimension

int Fermions::GetHilbertSpaceDimension()
{
  return this->HilbertSpaceDimension;
}

// Evaluate Hilbert space dimension
// 
// nbrSite = number of sites
// nUp = number of fermion with spin up
// nDown = number of fermion with spin down
// return value = Hilbert space dimension

int Fermions::EvaluateHilbertSpaceDimension (int nbrSite, int nUp, int nDown)
{
  int Dim = 1;
  int Div;
  int nUpFac = 1;
  for (int i = 1; i <= nUp; i++)
    nUpFac *= i;
  int nDownFac = 1;
  for (int i = 1; i <= nDown; i++)
    nDownFac *= i;
  for (int i = (nbrSite - nUp - nDown + 1); i <= nbrSite; i++)
    {
      Dim *= i;
      if (nUpFac > 1)
	{
	  Div = FindGCD (Dim, nUpFac);
	  if (Div != 1)
	    {
	      Dim /= Div;
	      nUpFac /= Div;
	    }
	}
      if (nDownFac > 1)
	{
	  Div = FindGCD (Dim, nDownFac);
	  if (Div != 1)
	    {
	      Dim /= Div;
	      nDownFac /= Div;
	    }      
	}
    }
  return Dim;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> Fermions::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->FixedQuantumNumberFlag == true)
    {
      VectorQuantumNumber* TmpQ = new VectorQuantumNumber();
      (*TmpQ) += new NumberParticleQuantumNumber (this->NbrFermion);  
      (*TmpQ) += new SzQuantumNumber (this->Sz);
      TmpList += TmpQ;
    }
  else
    {
      for (int TmpN = 0; TmpN <= this->NbrSite; TmpN++)
	{
	  for (int TmpSz = -TmpN; TmpSz <= TmpN; TmpSz += 2)
	    {
	      VectorQuantumNumber* TmpQ = new VectorQuantumNumber();
	      (*TmpQ) += new NumberParticleQuantumNumber (TmpN);  
	      (*TmpQ) += new SzQuantumNumber (TmpSz);
	      TmpList += TmpQ;
	    }
	}
    }
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* Fermions::GetQuantumNumber (int index)
{ 
  VectorQuantumNumber* TmpQ = new VectorQuantumNumber();
  (*TmpQ) += new NumberParticleQuantumNumber (this->GetNumberParticle(index));  
  (*TmpQ) += new SzQuantumNumber (this->TotalSz(index));
  return TmpQ;
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Fermions::TotalSz (int index)
{
  if (this->FixedQuantumNumberFlag == true)
    return this->Sz;
  unsigned long State = this->FermionDescription[index];
  int TmpSz = 0;
  unsigned long TmpState;
  for (int i = 0; i < this->NbrSite; i++)
    {
      TmpState = State & 0x00000003;
      switch (TmpState)
	{
	case 0x00000003:
	  TmpSz++;
	  break;
	case 0x00000002:
	  TmpSz--;
	  break;
	}
      State >>= 2;
    }
  return TmpSz;
}

// get number of particles for a given state
//
// index = index of the state to test
// return value = number of particles

int Fermions::GetNumberParticle (int index)
{
  if (this->FixedQuantumNumberFlag == true)
    return this->NbrFermion;
  unsigned long State = this->FermionDescription[index];
  int TmpN = 0;
  for (int i = 0; i < this->NbrSite; i++)
    {
      if ((State & 0x00000003) != 0)
	TmpN++;
      State >>= 2;
    }
  return TmpN;
}

//

// evaluate parity of a given state
//
// state = description of the state to evaluate
// return value = parity

unsigned long Fermions::EvaluateParity(unsigned long state)
{
  unsigned long Mask = 0x00000002 << (2 * (this->NbrSite - 1));
  unsigned long ParityMask = 0x00000003 << (2 * (this->NbrSite - 1));
  unsigned long Parity = 0x00000000;
  bool Flag = true;
  for (int i = 0; i < this->NbrSite; i++)
    {
      if ((state & Mask) != 0)
	{
	  Flag = !Flag;
	  if (Flag == true)
	    Parity |= ParityMask;
	}
      Mask >>= 2;
      ParityMask >>= 2;
    } 
  return Parity;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* Fermions::ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter)
{
  if (q.GetQuantumNumberType() != AbstractQuantumNumber::Vector)
    return new Fermions();
  VectorQuantumNumber& TmpQuantumNumber = (VectorQuantumNumber&) q;
  if ((TmpQuantumNumber[0] == 0) || 
      (TmpQuantumNumber[0]->GetQuantumNumberType() != AbstractQuantumNumber::NumberParticle))
    return new Fermions();
  if ((TmpQuantumNumber[1] == 0) || 
      (TmpQuantumNumber[1]->GetQuantumNumberType() != AbstractQuantumNumber::Sz))
    return new Fermions();
  int TmpN = ((NumberParticleQuantumNumber*) TmpQuantumNumber[0])->GetNumberParticle();
  int TmpSz = ((SzQuantumNumber*) (TmpQuantumNumber[1]))->GetSz();
  int HilbertSubspaceDimension = 0;
  if (((TmpSz >= 0) && ((TmpN % 2) == (TmpSz % 2))) || 
      ((TmpSz < 0) && ((TmpN % 2) == (-(TmpSz % 2)))))
    HilbertSubspaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrSite, 
								   (TmpN + TmpSz) / 2, 
								   (TmpN - TmpSz) / 2);
  else
    return new Fermions();
  cout << TmpN << " " << TmpSz << " " << HilbertSubspaceDimension << endl;
  int* TmpConvArray = new int [HilbertSubspaceDimension];
  int Pos = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; i++)
    {
      if ((this->TotalSz(i) == TmpSz) && (this->GetNumberParticle(i) == TmpN))
	{
	  TmpConvArray[Pos++] = i;
	}
    }
  int* ConvArray = new int [HilbertSubspaceDimension];
  unsigned long* SubspaceDescription = new unsigned long [HilbertSubspaceDimension];
  unsigned long* SubspaceParity = new unsigned long [HilbertSubspaceDimension];
  int* SubspaceLookUpTable = new int [this->LookUpTableSize];
  unsigned long TestMask = this->LookUpTableMask;
  for (int i = 0; i < HilbertSubspaceDimension; i++)
    {
      if ((this->FermionDescription[TmpConvArray[i]] & this->LookUpTableMask) != TestMask)
	{
	  TestMask = this->FermionDescription[TmpConvArray[i]] & this->LookUpTableMask;
	  SubspaceLookUpTable[TestMask] = i;
	}
      SubspaceDescription[i] = this->FermionDescription[TmpConvArray[i]];
      SubspaceParity[i] = this->Parity[TmpConvArray[i]];
      ConvArray[i] = TmpConvArray[i];
    }
  converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, HilbertSubspaceDimension, ConvArray);
  return new Fermions (HilbertSubspaceDimension, SubspaceDescription, SubspaceParity, this->NbrSite, TmpN,
		       TmpSz, true, SubspaceLookUpTable, this->LookUpTableSize, 
		       this->LookUpPosition, this->LookUpTableMask);
}

// find state index
//
// state = state description
// return value = corresponding index

int Fermions::FindStateIndex(unsigned long state)
{
  int index = this->LookUpTable[state & this->LookUpTableMask];
  unsigned long* tmpState = &(this->FermionDescription[index]);
  while ((index < this->HilbertSpaceDimension) && (state != *(tmpState++)))
    index++;
  return index;   
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Fermions::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long tmp;
  unsigned long StateDescription = this->FermionDescription[state];  
  Str << this->FindStateIndex(StateDescription) << " : ";
  for (int j = this->NbrSite - 1; j >= 0; j--)
    {
      tmp = (StateDescription >> (j * 2)) & 0x00000003;
      if (tmp == 0x00000002)
	Str << "C(" << j << ",d)";
      else
	if (tmp == 0x00000003)
	  Str << "C(" << j << ",u)";
    }
  Str << "|0>";
  return Str;
}

// find gcd
//
// m = first number
// n = second number
// return value = gcd of m and n

int FindGCD (int m, int n)
{
  if (m < 0)
    m = -m;
  if (n < 0)
    n = -n;
  if (m == n)
    return n;
  if (m < n)
    return gcd(m, n);
  else
    return gcd(n ,m);
}
 
// recurrence for gcd
//
// m = first number
// n = second number
// return value = gcd of m and n

int gcd (int m, int n)
{
  if (m == 0)
    return n;
  else
    return gcd ((n % m), m);
}                                                                                                                         
