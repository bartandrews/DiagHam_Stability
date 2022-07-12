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

#ifndef SLATERCOMPONENT_H
#define SLATERCOMPONENT_H

#include "config.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>
using std::ostream;

class SlaterSuperposition;


class SlaterComponent
{
 protected:
  double Prefactor;
  int NbrParticles;
  int* LzPositions;
  unsigned long StateDescription;
  int TotalLz;
  int LzMax;

  GarbageFlag Flag;
  
 public:
  
  // default constructor
  SlaterComponent();

  // standard constructors
  SlaterComponent(double prefactor, unsigned long stateDescription, int lzMax);
  SlaterComponent(double prefactor, int nbrParticles, int* positions, int lzMax);

  // constructor for highest angular momentum state
  SlaterComponent(int nbrParticles, int lzMax);

  // copy constructor 
  SlaterComponent(const SlaterComponent &toCopy);

  // destructor
  ~SlaterComponent();

  // assignment operator
  SlaterComponent& operator = (const SlaterComponent& toCopy);


  // comparison of StateDescription
  friend bool operator == (const SlaterComponent &lhs, const SlaterComponent &rhs);

  // addition of the prefactors
  SlaterComponent& operator += (const SlaterComponent& toAdd);

  // multiplication of a double to the prefactors
  SlaterComponent& operator *= (double factor);

  // angular momentum lowering operator
  SlaterSuperposition ApplyLMinus();

  // accessor methods:
  int GetTotalLz(){ return TotalLz; }
  int GetLzMax(){ return LzMax; }
  double GetPrefactor(){ return Prefactor; }

  int GetNbrParticles(){ return NbrParticles; }
  
  // calculate positions of bits in StateDescription
  // return newly created array
  int* GetLzPositions() const;
  // calculate LzValues corresponding to bits in StateDescription
  // return newly created array
  int* GetLzValues() const;
  // return values of positions as stored (quicker)
  int* GetLzPositionPtr() { return LzPositions; }

  // output method
  friend ostream& operator << (ostream & str, const SlaterComponent & s);
  
  
};

#endif
