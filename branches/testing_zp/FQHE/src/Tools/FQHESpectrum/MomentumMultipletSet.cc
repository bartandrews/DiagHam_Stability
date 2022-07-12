////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of set of momentum multiplets                    //
//                                                                            //
//                        last modification : 31/05/2005                      //
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
#include "Polynomial/IntegerPolynomial.h"
#include "Tools/FQHESpectrum/MomentumMultipletSet.h"

#include <iostream>


using std::cout;
using std::endl;


// default constructor
//
// maximumMomentum = twice the maximum momentum

MomentumMultipletSet::MomentumMultipletSet (int maximumMomentum)
{
  this->MaximumMomentum = maximumMomentum;
  this->LValues = new int [this->MaximumMomentum + 1];    
  for (int i = 0; i <= this->MaximumMomentum; ++i)
    this->LValues[i] = 0;
}

// constructor from raw datas
//
// maximumMomentum = twice the maximum momentum
// multiplets = array that contains number of multiplets
// duplicateFlag = true if the array has to be duplicated

MomentumMultipletSet::MomentumMultipletSet (int maximumMomentum, int* multiplets, bool duplicateFlag)
{
  this->MaximumMomentum = maximumMomentum;
  if (duplicateFlag == true)
    {
      this->LValues = new int [this->MaximumMomentum + 1];    
      for (int i = 0; i <= this->MaximumMomentum; ++i)
	this->LValues[i] = multiplets[i];
    }
  else
    this->LValues =  multiplets;
}

// constructor from an integer polynomial, power indicates Lz value up to a translation, coefficient indicates Lz multiplicty (for example, (1+2q+2q^2+2q^3+q^4) * q^6 
// corresponds to two multiplets L=1 and L=2)
//
// polynomial = reference on the integer polynomial

MomentumMultipletSet::MomentumMultipletSet (IntegerPolynomial& polynomial)
{
  int TmpMaximum = polynomial.GetPolynomialDegree();
  int TmpMinimum = 0;
  while ((TmpMinimum <= TmpMaximum) && (polynomial.PolynomialCoefficient(TmpMinimum) == 0l))
    {
      ++TmpMinimum;
    }
  this->MaximumMomentum = TmpMaximum - TmpMinimum;
  this->LValues = new int [this->MaximumMomentum + 1];
  for (int i = 0; i <= this->MaximumMomentum; ++i)
    this->LValues[i] = 0;
  TmpMinimum += TmpMaximum;
  TmpMinimum >>= 1;
  if ((this->MaximumMomentum & 1) == 1)
    ++TmpMinimum;
  this->LValues[this->MaximumMomentum] = (int) polynomial.PolynomialCoefficient(TmpMaximum);
  for (int i = TmpMaximum - 1; i >= TmpMinimum; --i)
    this->LValues[this->MaximumMomentum - ((TmpMaximum - i) << 1)] = (int) (polynomial.PolynomialCoefficient(i) - polynomial.PolynomialCoefficient(i + 1));  
}

// copy constructor
//
// multiplets = set of multiplets to copy 

MomentumMultipletSet::MomentumMultipletSet (const MomentumMultipletSet& multiplets)
{
  this->MaximumMomentum = multiplets.MaximumMomentum;
  this->LValues = new int [this->MaximumMomentum + 1];    
  for (int i = 0; i <= this->MaximumMomentum; ++i)
    this->LValues[i] = multiplets.LValues[i];
}

// destructor 
// 
MomentumMultipletSet::~MomentumMultipletSet()
{
  delete[] this->LValues;
}

// assignment
//
// multiplets = set of multiplets to assign 
// return value = reference on current momentum multiplet set

MomentumMultipletSet& MomentumMultipletSet::operator = (const MomentumMultipletSet& multiplets)
{
  delete[] this->LValues;
  this->MaximumMomentum = multiplets.MaximumMomentum;
  this->LValues = new int [this->MaximumMomentum + 1];    
  for (int i = 0; i <= this->MaximumMomentum; ++i)
    this->LValues[i] = multiplets.LValues[i];
  return *this;
}

// find all multiplets that can obtained when putting a given number of bosons where each of them having the same maximum momentum
//
// nbrBosons = number of bosons
// maximumMomentum = twice the maximum momentum that can have a boson
// return value = reference on corresponding momentum multiplet set

MomentumMultipletSet& MomentumMultipletSet::FindMultipletsForBosons (int nbrBosons, int maximumMomentum)
{
  this->Resize(nbrBosons * maximumMomentum);
  this->LValues[this->MaximumMomentum] = 1;
  int MinimumMomentum = 0;
  if ((this->MaximumMomentum & 1) == 1)
    ++MinimumMomentum;
  for (int i = this->MaximumMomentum - 2; i >= MinimumMomentum; i -= 2)
    this->LValues[i] = this->GetFixedLzBosonHilbertSpaceDimension(nbrBosons, maximumMomentum, (i + (maximumMomentum * nbrBosons)) >> 1);
  for (int i = MinimumMomentum; i < this->MaximumMomentum; i += 2)
    this->LValues[i] -= this->LValues[i + 2];
  return *this;
}
  
// find all multiplets that can obtained when putting a given number of fermions where each of them having the same maximum momentum
//
// nbrFermions = number of fermions
// maximumMomentum = twice the maximum momentum that can have a fermion
// return value = reference on corresponding momentum multiplet set

MomentumMultipletSet& MomentumMultipletSet::FindMultipletsForFermions (int nbrFermions, int maximumMomentum)
{
  this->Resize(nbrFermions * (maximumMomentum - nbrFermions + 1));
  this->LValues[this->MaximumMomentum] = 1;
  int MinimumMomentum = 0;
  if ((this->MaximumMomentum & 1) == 1)
    ++MinimumMomentum;
  for (int i = this->MaximumMomentum - 2; i >= MinimumMomentum; i -= 2)
    this->LValues[i] = this->GetFixedLzFermionHilbertSpaceDimension(nbrFermions, maximumMomentum, (i + (maximumMomentum * nbrFermions)) >> 1);
  for (int i = MinimumMomentum; i < this->MaximumMomentum; i += 2)
    this->LValues[i] -= this->LValues[i + 2];
  return *this;
}
  
// get the total number of states associated to the current set of multiplets
//
// return value = number of states

int MomentumMultipletSet::GetNbrStates()
{
  int NbrStates = 0;
  for (int i = 0; i <= this->MaximumMomentum; ++i)
    NbrStates += (i + 1) * this->LValues[i];
  return NbrStates;
}

// add a set of multiplets to another sets of multiplets (i.e. add the number multiplet for each L value)
// 
// multiplets = set of multiplets to add 
// return value = reference on current momentum multiplet set

MomentumMultipletSet& MomentumMultipletSet::operator += (const MomentumMultipletSet& multiplets)
{
  this->Resize(multiplets.MaximumMomentum);
  for (int i = 0; i < multiplets.MaximumMomentum; ++i)
    this->LValues[i] += multiplets.LValues[i];
  return *this;
}
  
// fuse two sets of multiplets according to momentum addition rules
//
// multiplets1 = first set of multiplets
// multiplets2 = second set of multiplets
// return value = resulting momentum multiplet set

MomentumMultipletSet operator * (MomentumMultipletSet& multiplets1, MomentumMultipletSet& multiplets2)
{
  int TmpMaximumMomentum = multiplets1.MaximumMomentum + multiplets2.MaximumMomentum;
  int* TmpLValues = new int [TmpMaximumMomentum + 1];
  for (int i = 0; i < TmpMaximumMomentum; ++i)
    TmpLValues[i] = 0;
  int Max;
  int Min;
  int Coef;
  for (int i = 0; i <= multiplets1.MaximumMomentum; ++i)
    if (multiplets1.LValues[i] != 0)
      for (int j = 0; j <= multiplets2.MaximumMomentum; ++j)      
	if (multiplets2.LValues[j] != 0)
	  {
	    Max = i + j;
	    Min = i - j;
	    if (Min < 0)
	      Min *= -1;
	    Coef = multiplets1.LValues[i] * multiplets2.LValues[j];
	    while (Min <= Max)
	      {
		TmpLValues[Min] += Coef;
		Min += 2;
	      }
	  }
  return MomentumMultipletSet(TmpMaximumMomentum, TmpLValues);
} 

// fuse two sets of multiplets according to momentum addition rules and add result to another set of multiplets (i.e. add the number multiplet for each L value)
//
// multiplets1 = first set of multiplets
// multiplets2 = second set of multiplets
// return value = reference on current momentum multiplet set

MomentumMultipletSet& MomentumMultipletSet::FuseAndAdd (MomentumMultipletSet& multiplets1, MomentumMultipletSet& multiplets2)
{
  int TmpMaximumMomentum = multiplets1.MaximumMomentum + multiplets2.MaximumMomentum;
  this->Resize(TmpMaximumMomentum);
  int Max;
  int Min;
  int Coef;
  for (int i = 0; i <= multiplets1.MaximumMomentum; ++i)
    if (multiplets1.LValues[i] != 0)
      for (int j = 0; j <= multiplets2.MaximumMomentum; ++j)      
	if (multiplets2.LValues[j] != 0)
	  {
	    Max = i + j;
	    Min = i - j;
	    if (Min < 0)
	      Min *= -1;
	    Coef = multiplets1.LValues[i] * multiplets2.LValues[j];
	    while (Min <= Max)
	      {
		this->LValues[Min] += Coef;
		Min += 2;
	      }
	  }
  return *this;
}

// resize (aka change maximum momentum value, must be greater than the current one)
//
// maximumMomentum = twice the new maximum momentum
// return value = reference on current momentum multiplet set

MomentumMultipletSet& MomentumMultipletSet::Resize(int maximumMomentum)
{  
  if (this->MaximumMomentum >= maximumMomentum)
    return *this;
  int* TmpLValues = new int [maximumMomentum + 1];
  int i = 0;
  for (; i < this->MaximumMomentum; ++i)
    TmpLValues[i] = this->LValues[i];
  for (; i < maximumMomentum; ++i)
    TmpLValues[i] = 0;      
  delete[] this->LValues;
  this->LValues = TmpLValues;
  this->MaximumMomentum = maximumMomentum;
  return *this;
}

// output stream overload (output only half integer momenta if no integer momentum is present, idem for integer momenta)
//
// str = reference on the output stream
// multiplets = set of momentum multiplets to display
// return value = reference on the output stream

ostream& operator << (ostream& str, const MomentumMultipletSet& multiplets)
{
  int i = 0;
  bool Flag = false;
  while ((i <= multiplets.MaximumMomentum) && (Flag == false))
    {
      if (multiplets.LValues[i] != 0)
	Flag = true;
      i += 2;
    }
  if (Flag == false)
    {
      i = 1;
      while (i <= multiplets.MaximumMomentum)
	{
	  str << i << "/2: " << multiplets.LValues[i] << "  ";
	  i += 2;
	}
      return str;
    }
  Flag = false;
  i = 1;
  while ((i <= multiplets.MaximumMomentum) && (Flag == false))
    {
      if (multiplets.LValues[i] != 0)
	Flag = true;
      i += 2;
    }
  if (Flag == false)
    {
      i = 0;
      while (i <= multiplets.MaximumMomentum)
	{
	  str << (i >> 1) << ": " << multiplets.LValues[i] << "  ";
	  i += 2;
	}
      return str;
    }
  i = 0;
  while (i <= multiplets.MaximumMomentum)
    {
      str << (i >> 1) << ": " << multiplets.LValues[i] << "  ";
      ++i;
      if (i <= multiplets.MaximumMomentum)
	{
	  str << i << "/2: " << multiplets.LValues[i] << "  ";
	  ++i;
	}
    }
  return str;
}

// display mulitplet with respect to Lz instead of L (output only half integer momenta if no integer momentum is present, idem for integer momenta)
//
// str = reference on the output stream
// return value = reference on the output stream

ostream& MomentumMultipletSet::PrintLz (ostream& str)
{
  int i = 0;
  bool Flag = false;
  while ((i <= this->MaximumMomentum) && (Flag == false))
    {
      if (this->LValues[i] != 0)
	Flag = true;
      i += 2;
    }
  if (Flag == false)
    {
      i = 1;
      while (i <= this->MaximumMomentum)
	{
	  int Value = 0;
	  for (int j = i; j <= this->MaximumMomentum; j += 2)
	    Value += this->LValues[j];
	  str << i << "/2: " << Value << "  ";
	  i += 2;
	}
      return str;
    }
  Flag = false;
  i = 1;
  while ((i <= this->MaximumMomentum) && (Flag == false))
    {
      if (this->LValues[i] != 0)
	Flag = true;
      i += 2;
    }
  if (Flag == false)
    {
      i = 0;
      while (i <= this->MaximumMomentum)
	{
	  int Value = 0;
	  for (int j = i; j <= this->MaximumMomentum; j += 2)
	    Value += this->LValues[j];
	  str << (i >> 1) << ": " << Value << "  ";
	  i += 2;
	}
      return str;
    }
  i = 0;
  while (i <= this->MaximumMomentum)
    {
      int Value = 0;
      for (int j = i; j <= this->MaximumMomentum; j += 2)
	Value += this->LValues[j];
      str << (i >> 1) << ": " << Value << "  ";
      ++i;
      if (i <= this->MaximumMomentum)
	{
	  Value = 0;
	  for (int j = i; j <= this->MaximumMomentum; j += 2)
	    Value += this->LValues[j];
	  str << i << "/2: " << Value << "  ";
	  ++i;
	}
    }
  return str;
}

// evaluate Hilbert space dimension for bosons with fixed total Lz value and a given L per particle
//
// nbrBosons = number of bosons
// lzMax = two times momentum maximum value for a boson plus one 
// totalLz = momentum total value plus nbrBosons  * (momentum maximum value for a nbrBosons + 1)
// return value = Hilbert space dimension

int MomentumMultipletSet::GetFixedLzBosonHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  if ((nbrBosons == 0) || ((nbrBosons * lzMax) < totalLz))
    return 0;
  if (((nbrBosons * lzMax) == totalLz) || (lzMax == 0) || (totalLz == 0))
    {
      return 1;
    }
  int TmpDim = 0;
  while ((totalLz >= 0) && (nbrBosons > 0))
    {
      TmpDim += this->GetFixedLzBosonHilbertSpaceDimension(nbrBosons, lzMax - 1, totalLz);
      --nbrBosons;
      totalLz -= lzMax;
    }
  return TmpDim; 
}

// evaluate Hilbert space dimension for fermions with fixed total Lz value and a given L per particle
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions  * (momentum maximum value for a fermions + 1)
// return value = Hilbert space dimension

int MomentumMultipletSet::GetFixedLzFermionHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (lzMax < (nbrFermions - 1)))
    return (long) 0;
  int LzTotalMax = ((2 * lzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (LzTotalMax < totalLz)
    return (long) 0;
  if ((nbrFermions == 1) && (lzMax >= totalLz))
    return (long) 1;
  if (LzTotalMax == totalLz)
    return (long) 1;
  return  (this->GetFixedLzFermionHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax)
	   +  this->GetFixedLzFermionHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz));
}

