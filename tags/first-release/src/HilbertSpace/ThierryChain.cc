////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                        last modification : 26/02/2001                      //
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


#include "HilbertSpace/ThierryChain.h"
#include <iostream>


using std::hex;
using std::dec;


#define M_SQRT3 1.73205080756888
#define M_SQRT6 2.44948974278318
#define MAXHILBERTSPACEDIMENSION 8581300


// default constructor
//

ThierryChain::ThierryChain ()
{
  this->GarbageFlag = 0;
  this->Spin2ChainLength = 0;
  this->Spin3_2ChainLength = 0;
  this->Spin3_2Start = 0;
  this->HilbertSpaceDimension = 0;
  this->Sz = 0;
  this->MaxSpin3_2ContributionToSz = 0;      
  this->ChainDescription = 0;
  this->LookUpTable = 0;
  this->LookUpTableMask = 0;
  this->LookUpPosition = 0;
  this->LookUpPosition2 = 0;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// spin2ChainLength = number of spin 2
// spin3_2ChainLength = number of spin 3/2
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table

ThierryChain::ThierryChain (int spin2ChainLength, int spin3_2ChainLength, int sz, int memorySize) 
{
  this->GarbageFlag = new int;
  *(this->GarbageFlag) = 1;
  this->Spin2ChainLength = spin2ChainLength;
  this->Spin3_2ChainLength = spin3_2ChainLength;
  if ((3 * this->Spin2ChainLength + 2 * this->Spin3_2ChainLength) > 32)
    {
      this->Spin2ChainLength = 1;
      this->Spin3_2ChainLength = 1;
    }
  this->Spin3_2Start = 3 * this->Spin2ChainLength;
  this->Sz = sz;
  this->MaxSpin3_2ContributionToSz = this->Spin3_2ChainLength * 6;

  this->HilbertSpaceDimension = MAXHILBERTSPACEDIMENSION;

  this->LookUpPosition = 0;
  int LookUpTableSize = 1;
  memorySize >>= 2;
  this->LookUpTableMask = 0xffffffff;
  while ((this->LookUpPosition < this->Spin2ChainLength) && (memorySize >=  8))
    {
      this->LookUpTableMask <<= 3;
      LookUpTableSize *= 8;
      memorySize >>= 3;
      this->LookUpPosition++;
    }
  if (memorySize >=  4)
    {
      this->LookUpTableMask = 0xffffffff << this->Spin3_2Start;
      this->LookUpPosition -= this->Spin2ChainLength;
      while ((this->LookUpPosition < this->Spin3_2ChainLength) && (memorySize >=  4))
	{
	  this->LookUpTableMask <<= 2;
	  LookUpTableSize *= 4;
	  memorySize >>= 2;
	  this->LookUpPosition++;
	}
      this->LookUpPosition += this->Spin2ChainLength;
   }
  this->LookUpPosition2 = this->LookUpPosition - this->Spin2ChainLength; 
  this->LookUpTableMask = ~this->LookUpTableMask;
  this->LookUpTable = new int [LookUpTableSize];
  this->ChainDescription = new unsigned long [this->HilbertSpaceDimension];
  this->ChainDescription[0] = 0xffffffff; 
  this->HilbertSpaceDimension = this->GenerateSpin2States (0, 0, 0xffffffff, 
							   (this->Spin2ChainLength * 4) + (this->Spin3_2ChainLength * 3));
}

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

ThierryChain::ThierryChain (const ThierryChain& chain)
{
  if (chain.GarbageFlag != 0)
    {
      this->GarbageFlag = chain.GarbageFlag;
      (*(this->GarbageFlag))++;
      this->Spin2ChainLength = chain.Spin2ChainLength;
      this->Spin3_2ChainLength = chain.Spin3_2ChainLength;
      this->Spin3_2Start = chain.Spin3_2Start;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->Sz = chain.Sz;
      this->MaxSpin3_2ContributionToSz = chain.MaxSpin3_2ContributionToSz;      
      this->ChainDescription = chain.ChainDescription;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpPosition2 = chain.LookUpPosition2;
    }
  else
    {
      this->GarbageFlag = 0;
      this->Spin2ChainLength = 0;
      this->Spin3_2ChainLength = 0;
      this->Spin3_2Start = 0;
      this->HilbertSpaceDimension = 0;
      this->Sz = 0;
      this->MaxSpin3_2ContributionToSz = 0;      
      this->ChainDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->LookUpPosition2 = 0;
    }
}

// destructor
//

ThierryChain::~ThierryChain () 
{
  if (this->GarbageFlag != 0) 
    if ((*(this->GarbageFlag)) == 1)
      {
	delete[] this->ChainDescription;
	delete this->GarbageFlag;
      }
    else
      (*(this->GarbageFlag))--;
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

ThierryChain& ThierryChain::operator = (const ThierryChain& chain)
{
  if (this->GarbageFlag != 0)
    if ((*(this->GarbageFlag)) == 1)
      {
	delete[] this->ChainDescription;
	delete this->GarbageFlag;
      }  
    else
      (*(this->GarbageFlag))--;
  if (chain.GarbageFlag != 0)
    {
      this->GarbageFlag = chain.GarbageFlag;
      (*(this->GarbageFlag))++;
      this->Spin2ChainLength = chain.Spin2ChainLength;
      this->Spin3_2ChainLength = chain.Spin3_2ChainLength;
      this->Spin3_2Start = chain.Spin3_2Start;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->Sz = chain.Sz;
      this->MaxSpin3_2ContributionToSz = chain.MaxSpin3_2ContributionToSz;      
      this->ChainDescription = chain.ChainDescription;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpPosition2 = chain.LookUpPosition2;
    }
  else
    {
      this->GarbageFlag = 0;
      this->Spin2ChainLength = 0;
      this->Spin3_2ChainLength = 0;
      this->Spin3_2Start = 0;
      this->HilbertSpaceDimension = 0;
      this->Sz = 0;
      this->MaxSpin3_2ContributionToSz = 0;      
      this->ChainDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->LookUpPosition2 = 0;
    }
  return *this;
}

// re-initialize chain with another total Sz component
//
// sz = twice the value of total Sz component
// return value = reference on current chain

ThierryChain& ThierryChain::Reinitialize(int sz)
{ 
  this->Sz = sz;
  this->HilbertSpaceDimension = this->GenerateSpin2States (0, 0, 0xffffffff, 
							   (this->Spin2ChainLength * 4) + (this->Spin3_2ChainLength * 3));
  return *this;
}

// generate Spin 2 states for a given total spin projection Sz
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// currentStateDescription = description of current state
// currentSz = total Sz value of current state
// return value = number of generated states

int ThierryChain::GenerateSpin2States(int statePosition, int sitePosition, unsigned currentStateDescription, int currentSz) 
{
  int MaxSz = (this->Spin2ChainLength - sitePosition) * 8 + MaxSpin3_2ContributionToSz;
  int DiffSz = currentSz - this->Sz;
  if (DiffSz > MaxSz)
    return 0;
  int NbrGeneratedState = -statePosition;
  int NextSitePosition = sitePosition + 1;
  sitePosition &= 0x0000000f;
  sitePosition *= 3;
  unsigned long mask;
  if (NextSitePosition != this->Spin2ChainLength)
    {
      mask = ~(0x00000005 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateSpin2States(statePosition, NextSitePosition, currentStateDescription & mask, currentSz);
      
      mask = ~(0x00000006 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateSpin2States(statePosition, NextSitePosition, currentStateDescription & mask, currentSz - 2);
      
      mask = ~(0x00000007 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateSpin2States(statePosition, NextSitePosition, currentStateDescription & mask, currentSz - 4);
      
      mask = ~(0x00000002 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateSpin2States(statePosition, NextSitePosition, currentStateDescription & mask, currentSz - 6);

      mask = ~(0x00000001 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateSpin2States(statePosition, NextSitePosition, currentStateDescription & mask, currentSz - 8);
    }
  else
    {
      mask = ~(0x00000005 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateSpin3_2States(statePosition, 0, currentStateDescription & mask, currentSz);
      
      mask = ~(0x00000006 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateSpin3_2States(statePosition, 0, currentStateDescription & mask, currentSz - 2);
      
      mask = ~(0x00000007 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateSpin3_2States(statePosition, 0, currentStateDescription & mask, currentSz - 4);
      
      mask = ~(0x00000002 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateSpin3_2States(statePosition, 0, currentStateDescription & mask, currentSz - 6);

      mask = ~(0x00000001 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateSpin3_2States(statePosition, 0, currentStateDescription & mask, currentSz - 8);
    }
  NbrGeneratedState += statePosition;
  return NbrGeneratedState;
}

// generate Spin 3/2 states for a given total spin projection Sz
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// currentStateDescription = description of current state
// currentSz = total Sz value of current state
// return value = number of generated states

int ThierryChain::GenerateSpin3_2States(int statePosition, int sitePosition, unsigned currentStateDescription, int currentSz) 
{
  int MaxSz = (this->Spin3_2ChainLength - sitePosition) * 6;
  int DiffSz = currentSz - this->Sz;
  if (DiffSz > MaxSz)
    return 0;
  int NbrGeneratedState = -statePosition;
  int NextSitePosition = sitePosition + 1;
  unsigned long mask;
  sitePosition &= 0x0000000f;
  sitePosition <<= 1;
  sitePosition += this->Spin3_2Start;
  if (NextSitePosition != this->Spin3_2ChainLength)
    {
      mask = ~(0x00000002 << sitePosition);
      if (NextSitePosition == this->LookUpPosition2)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateSpin3_2States(statePosition, NextSitePosition, currentStateDescription & mask, currentSz);
      
      mask = ~(0x00000003 << sitePosition);
      if (NextSitePosition == this->LookUpPosition2)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateSpin3_2States(statePosition, NextSitePosition, currentStateDescription & mask, currentSz - 2);
      
      mask = ~(0x00000001 << sitePosition);
      if (NextSitePosition == this->LookUpPosition2)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateSpin3_2States(statePosition, NextSitePosition, currentStateDescription & mask, currentSz - 4);

      if (NextSitePosition == this->LookUpPosition2)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateSpin3_2States(statePosition, NextSitePosition, currentStateDescription, currentSz - 6);
    }
  else
    {
      if (currentSz == this-> Sz)
	{
	  mask = ~(0x00000002 << sitePosition);
	  if (NextSitePosition == this->LookUpPosition2)
	    {
	      this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
	    }
	  this->ChainDescription[statePosition] = (currentStateDescription & mask);
	  statePosition++;
	}
      else
	{
	  currentSz -= 2;
	  if (currentSz == this-> Sz)
	    {
	      mask = ~(0x00000003 << sitePosition);
	      if (NextSitePosition == this->LookUpPosition2)
		{
		  this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
		}
	      this->ChainDescription[statePosition] = (currentStateDescription & mask);
	      statePosition++;
	    }
	  else
	    {
	      currentSz -= 2;
	      if (currentSz == this-> Sz)
		{
		  mask = ~(0x00000001 << sitePosition);
		  if (NextSitePosition == this->LookUpPosition2)
		    {
		      this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;		      
		    }
		  this->ChainDescription[statePosition] = (currentStateDescription & mask);
		  statePosition++;
		}
	      else
		{
		  currentSz -= 2;
		  if (currentSz == this-> Sz)
		    {
		      if (NextSitePosition == this->LookUpPosition2)
			{
			  this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
			}
		      this->ChainDescription[statePosition] = currentStateDescription;
		      statePosition++;
		    }
		}	  
	    }
	}
    }
  NbrGeneratedState += statePosition;
  return NbrGeneratedState;
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = pointer to double where numerical coefficient has to be stored
// return value = index of resulting state

int ThierryChain::SmiSpj (int i, int j, int state, double* coefficient)
{  
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long tmpMask;
  if (i < this->Spin2ChainLength)
    {
      i *= 3;
      tmpMask = (tmpState >> i) & 0x00000007;
      tmpState &= ~(0x00000007 << i);
      switch (tmpMask)
	{
	case 0x00000002:
	  tmpState |= (0x00000001 << i);
	  (*coefficient) = 2.0;
	  break;
	case 0x00000001:
	  (*coefficient) = M_SQRT6;
	  break;
	case 0x00000000:
	  tmpState |= (0x00000005 << i);
	  (*coefficient) = M_SQRT6;
	  break;
	case 0x00000005:
	  tmpState |= (0x00000006 << i);
	  (*coefficient) = 2.0;
	  break;
	case 0x00000006:
	  return this->HilbertSpaceDimension;
	}
    }
  else
    {
      i -= this->Spin2ChainLength;
      i <<= 1;
      i += this->Spin3_2Start;
      tmpMask = (tmpState >> i) & 0x00000003;
      tmpState &= ~(0x00000003 << i);
      switch (tmpMask)
	{
	case 0x00000001:
	  (*coefficient) = M_SQRT3;
	  break;
	case 0x00000000:
	  tmpState |= (0x00000002 << i);
	  (*coefficient) = 2.0;
	  break;
	case 0x00000002:
	  tmpState |= (0x00000003 << i);
	  (*coefficient) = M_SQRT3;
	  break;
	case 0x00000003:
	  return this->HilbertSpaceDimension;
	}
    }
  if (j < this->Spin2ChainLength)
    {
      j *= 3;
      tmpMask = (tmpState >> j) & 0x00000007;
      tmpState &= ~(0x00000007 << j);
      switch (tmpMask)
	{
	case 0x00000002:
	  return this->HilbertSpaceDimension;
	case 0x00000001:
	  tmpState |= (0x00000002 << j);
	  (*coefficient) *= 2.0;
	  break;
	case 0x00000000:
	  tmpState |= (0x00000001 << j);
	  (*coefficient) *= M_SQRT6;
	  break;
	case 0x00000005:
	  (*coefficient) *= M_SQRT6;
	  break;
	case 0x00000006:
	  tmpState |= (0x00000005 << j);
	  (*coefficient) *= 2.0;
	  break;
	}
    }
  else
    {
      j -= this->Spin2ChainLength;
      j <<= 1;
      j += this->Spin3_2Start;
      tmpMask = (tmpState >> j) & 0x00000003;
      tmpState &= ~(0x00000003 << j);
      switch (tmpMask)
	{
	case 0x00000001:
	  return this->HilbertSpaceDimension;
	case 0x00000000:
	  tmpState |= (0x00000001 << j);
	  (*coefficient) *= M_SQRT3;
	  break;
	case 0x00000002:
	  (*coefficient) *= 2.0;
	  break;
	case 0x00000003:
	  tmpState |= (0x00000002 << j);
	  (*coefficient) *= M_SQRT3;
	  break;
	}
    }
  return this->FindStateIndex(tmpState);
}

// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = position of left Sz operator
// j = position of right Sz operator
// state = index of the state to consider
// return value = corresponding eigenvalue

double ThierryChain::SziSzj (int i, int j, int state)
{
  double x = 1.0;
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long tmpMask;
  if (i < this->Spin2ChainLength)
    {
      i *= 3;
      tmpMask = (tmpState >> i) & 0x00000007;
      switch (tmpMask)
	{
	case 0x00000002:
	  x = 2.0;
	  break;
	case 0x00000000:
	  return 0.0;
	  break;
	case 0x00000005:
	  x = -1.0;
	  break;
	case 0x00000006:
	  x = -2.0;
	  break;
	}
    }
  else
    {
      i -= this->Spin2ChainLength;
      i <<= 1;
      i += this->Spin3_2Start;
      tmpMask = (tmpState >> i) & 0x00000003;
      switch (tmpMask)
	{
	case 0x00000001:
	  x = 1.5;
	  break;
	case 0x00000000:
	  x = 0.5;
	  break;
	case 0x00000002:
	  x = -0.5;
	  break;
	case 0x00000003:
	  x = -1.5;
	  break;
	}
    }
  if (j < this->Spin2ChainLength)
    {
      j *= 3;
      tmpMask = (tmpState >> j) & 0x00000007;
      switch (tmpMask)
	{
	case 0x00000002:
	  x *= 2.0;
	  break;
	case 0x00000000:
	  return 0.0;
	case 0x00000005:
	  x *= -1.0;
	  break;
	case 0x00000006:
	  x *= -2.0;
	  break;
	}
    }
  else
    {
      j -= this->Spin2ChainLength;
      j <<= 1;
      j += this->Spin3_2Start;
      tmpMask = (tmpState >> j) & 0x00000003;
      switch (tmpMask)
	{
	case 0x00000001:
	  x *= 1.5;
	  break;
	case 0x00000000:
	  x *= 0.5;
	  break;
	case 0x00000002:
	  x *= -0.5;
	  break;
	case 0x00000003:
	  x *= -1.5;
	  break;
	}
    }
  return x;
}

// find state index
//
// state = state description
// return value = corresponding index

int ThierryChain::FindStateIndex(unsigned long state)
{
  int index = this->LookUpTable[state & this->LookUpTableMask];
  unsigned long* tmpState = &(this->ChainDescription[index]);
  while ((index < this->HilbertSpaceDimension) && (state != *(tmpState++)))
    index++;
  return index;
    
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& ThierryChain::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long tmp;
//  if (state != this->FindStateIndex(this->ChainDescription[state]))
//    Str << "error ";
   Str << LookUpTable[this->ChainDescription[state] & this->LookUpTableMask] << " " << hex << this->ChainDescription[state] << dec << "   ";
  int max = 3 * this->Spin2ChainLength;
  for (int k = 0; k < max; k += 3)
    {
      tmp = ((this->ChainDescription[state] >> k) & 0x00000007);
      switch (tmp)
	{
	case 0x00000002:
	  Str << "+2 ";
	  break;
	case 0x00000001:
	  Str << "+1 ";
	  break;
	case 0x00000000:
	  Str << " 0 ";
	  break;
	case 0x00000005:
	  Str << "-1 ";
	  break;
	case 0x00000006:
	  Str << "-2 ";
	  break;	  
	}
    }
  max = this->Spin3_2Start + 2 * this->Spin3_2ChainLength;
  for (int k = this->Spin3_2Start; k < max; k += 2)
    {
      tmp = ((this->ChainDescription[state] >> k) & 0x00000003);
      switch (tmp)
	{
	case 0x00000001:
	  Str << "+3/2 ";
	  break;
	case 0x00000000:
	  Str << "+1/2 ";
	  break;
	case 0x00000002:
	  Str << "-1/2 ";
	  break;
	case 0x00000003:
	  Str << "-3/2 ";
	  break;
	}      
    }
  return Str;
}
