////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2008 Gunnar Moeller                   //
//                                                                            //
//                                                                            //
//           class implementing a single slater determinant form              //
//                                                                            //
//                        last modification : 16/01/2008                      //
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


#include "SlaterComponent.h"
#include "SlaterSuperposition.h"
#include "GeneralTools/UnsignedIntegerTools.h"

#include <iostream>
#include <bitset>
#include <cmath>

using std::cout;
using std::endl;
using std::sqrt;
using std::bitset;



// default constructor
SlaterComponent::SlaterComponent()
{
  this->Prefactor=0.0;
  this->StateDescription=0x0l;
  this->NbrParticles=0;
  this->TotalLz=0;
  this->LzMax=0;
  this->LzPositions=0;
}


// standard constructor
SlaterComponent::SlaterComponent(double prefactor, unsigned long stateDescription, int lzMax)
{
  this->Prefactor=prefactor;
  this->StateDescription=stateDescription;
  this->LzMax = lzMax;
  this->NbrParticles=bitcount(StateDescription);
  this->TotalLz = 0;
  this->LzPositions = new int[this->NbrParticles];
  int count=0;
  for (int i=0; i<=lzMax; ++i)
    if (StateDescription & (0x1l<<i))
      {
	this->TotalLz+=(2*i)-LzMax;
	this->LzPositions[count++] = i;
      }
  this->Flag.Initialize();
}

SlaterComponent::SlaterComponent(double prefactor, int nbrParticles, int* positions, int lzMax)
{
  this->Prefactor=prefactor;
  this->LzMax = lzMax;
  this->NbrParticles=nbrParticles;
  this->TotalLz = 0;
  this->LzPositions = positions;
  this->StateDescription=0x0l;
  for (int i=0; i<nbrParticles; ++i)
    {
      this->StateDescription|=0x1l<<(LzPositions[i]);
      this->TotalLz+=(2*LzPositions[i])-LzMax;
    }
  this->Flag.Initialize();
}


// constructor for highest angular momentum state
SlaterComponent::SlaterComponent(int nbrParticles, int lzMax)
{
  this->NbrParticles=nbrParticles;
  this->LzMax=lzMax;
  this->LzPositions = new int[nbrParticles];
  if (nbrParticles > LzMax+1)
    {
      cout << "Too many particles requested!"<<endl;
      exit(1);
    }
  this->Prefactor=1.0;
  this->StateDescription=0x0l;
  this->TotalLz = 0;
  for (int i=LzMax-nbrParticles+1; i<=LzMax; ++i)
    {
      this->StateDescription |= 0x1l << i;
      this->TotalLz += 2*i-LzMax;
      this->LzPositions[i]=i;
    }
  this->Flag.Initialize();
}

// copy constructor 
SlaterComponent::SlaterComponent(const SlaterComponent &toCopy)
{
  this->NbrParticles=toCopy.NbrParticles;
  this->LzMax=toCopy.LzMax;
  this->LzPositions = toCopy.LzPositions;
  this->Prefactor=toCopy.Prefactor;
  this->StateDescription=toCopy.StateDescription;
  this->TotalLz = toCopy.TotalLz;
  this->Flag=toCopy.Flag;
}

SlaterComponent::~SlaterComponent()
{
  if ((this->NbrParticles != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->LzPositions;
    }
}


  // assignment operator
SlaterComponent& SlaterComponent::operator = (const SlaterComponent& toCopy)
{
  if ((this->NbrParticles != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->LzPositions;
    }
  this->NbrParticles=toCopy.NbrParticles;
  this->LzMax=toCopy.LzMax;
  this->LzPositions = toCopy.LzPositions;
  this->Prefactor=toCopy.Prefactor;
  this->StateDescription=toCopy.StateDescription;
  this->TotalLz = toCopy.TotalLz;
  this->Flag=toCopy.Flag;
  return *this;
}

// addition of the prefactors
SlaterComponent& SlaterComponent::operator += (const SlaterComponent& toAdd)
{
  if (*this == toAdd)
    this->Prefactor+=toAdd.Prefactor;
  return *this;
}

// multiplication of a double to the prefactors
SlaterComponent& SlaterComponent::operator *= (double factor)
{
  this->Prefactor *= factor;
  return *this;
}

// comparison of StateDescription
bool operator == (const SlaterComponent &lhs, const SlaterComponent &rhs)
{
  return (lhs.StateDescription == rhs.StateDescription);
}


SlaterSuperposition SlaterComponent::ApplyLMinus()
{
  bitset<10> b;
  b = this->StateDescription;
  cout << "Applying L- to " << b << endl;
  SlaterSuperposition result;
  double LLPlus1=0.25*LzMax*(LzMax+2);
  double MMMinus1;
  unsigned long newStateDescription;
  for (int i=1; i<=this->LzMax; ++i)
    {
      cout << "testing i=" <<i<<endl;
      if ( ((this->StateDescription & (0x3l<<(i-1)))>>(i-1)) == 0x2l )
	{
	  b=(this->StateDescription & (0x2l<<(i-1)))>>(i-1);
	  cout << "masked result: " << b << endl;
	  MMMinus1=0.25*(2*i-LzMax)*(2*i-LzMax-2);
	  newStateDescription = StateDescription ^ (0x3l<<(i-1));
	  b=newStateDescription;
	  cout << "adding term with " << b << endl;
	  SlaterComponent Tmp(sqrt(LLPlus1-MMMinus1)*this->Prefactor, newStateDescription, LzMax);
	  result += Tmp;
	}
    }
  return result;
}


// calculate positions of bits in StateDescription
int* SlaterComponent::GetLzPositions() const
{
  int *result = new int[NbrParticles];
  int found=0, pos=0;
  while ((found < NbrParticles) && (pos <= LzMax))
    {
      if (StateDescription & (0x1l<<pos))
	{
	  result[found]=pos;
	  ++found;
	}
      ++pos;
    }
  if (found != NbrParticles)
    cout << "Problem with Number of Particles in SlaterComponent" << endl;
  return result;
}

// calculate LzValues corresponding to bits in StateDescription
int* SlaterComponent::GetLzValues() const
{
  int *result = this->GetLzPositions();
  for (int i=0; i<NbrParticles; ++i)
    result[i]= 2*result[i]-LzMax;
  return result;
}

// output method
ostream& operator << (ostream & str, const SlaterComponent & s)
{
  if (s.NbrParticles==0)
    str << "0 ";
  else
    {
      int *lzValues = s.GetLzValues();
      if (s.Prefactor!=1.0)
	str << s.Prefactor<<"*";
      str << "Slater[ ";
      if (lzValues[0]&1) str<<lzValues[0]<< "/2";
      else cout <<lzValues[0]/2.0;
      for (int i=1; i<s.NbrParticles; ++i)
	{
	  str << ", ";
	  if (lzValues[i]&1) str<<lzValues[i]<< "/2";
	  else cout <<lzValues[i]/2.0;
	}
      str <<"] ";
    }
  return str;
}
