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

#ifndef SLATERSUPERPOSITION_H
#define SLATERSUPERPOSITION_H

#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "GeneralTools/List.h"

#include "SlaterComponent.h"



class SlaterSuperposition
{

 protected:

  // list of elements during build-up
  List<SlaterComponent> Elements;

  // number of elements
  int NbrElements;

  // total Angular Momentum
  int LTotal;

  // z-component of total angular momentum
  int LzTotal;
  
  // GarbageFlag
  GarbageFlag Flag;
  
 public:
  
  // default constructor
  SlaterSuperposition();

  // constructor from list of elements
  SlaterSuperposition(int lTotal, int lzTotal, int nbrElements, SlaterComponent *elements);

  // constructor for highest angular momentum eigenstate
  SlaterSuperposition(int nbrParticles, int lzMax);

  // constructor beginning with a single term
  SlaterSuperposition(SlaterComponent &element);

  // copy constructor 
  SlaterSuperposition(const SlaterSuperposition &toCopy);

  // destructor
  ~SlaterSuperposition();

  // assignment operator
  SlaterSuperposition& operator = (const SlaterSuperposition& toCopy);
  
  // addition of an element
  SlaterSuperposition& operator += (const SlaterComponent& toAdd);

  // addition of a superposition
  SlaterSuperposition& operator += (const SlaterSuperposition& toAdd);

  // multiplication with a prefactor
  SlaterSuperposition& operator *= (double factor);

  SlaterSuperposition ApplyLMinus();

  // accessor methods
  void SetLTotal(int lTotal);
  void SetLzTotal(int lzTotal);

  friend ostream& operator << (ostream & str, const SlaterSuperposition & s);
  
};

#endif
