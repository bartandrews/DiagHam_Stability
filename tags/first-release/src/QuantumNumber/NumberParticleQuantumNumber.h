////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of number of particles quantum number               //
//                                                                            //
//                        last modification : 11/05/2001                      //
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


#ifndef NUMBERPATICLEQUANTUMNUMBER_H
#define NUMBERPATICLEQUANTUMNUMBER_H


#include "config.h"
#include "QuantumNumber/AbstractQuantumNumber.h"

#include <iostream>


using std::ostream;


class NumberParticleQuantumNumber : public AbstractQuantumNumber
{

 protected:

  int NumberParticle;

 public:

  // default constructor
  //
  NumberParticleQuantumNumber ();

  // constructor from a numberof particles
  //
  // numberParticle = number of particles
  NumberParticleQuantumNumber (int numberParticle);

  // copy constructor
  //
  // Q = quantum number to copy
  NumberParticleQuantumNumber (const NumberParticleQuantumNumber& Q);

  // destructor
  //
  ~NumberParticleQuantumNumber ();

  // assignement
  //
  // Q = quantum number to copy
  // return value = reference on current quantum number
  NumberParticleQuantumNumber& operator = (const NumberParticleQuantumNumber& Q);

  // clone current quantum number
  //
  // return value = pointer on cloned quantum number  
  AbstractQuantumNumber* Clone ();

  // set number of particles
  //
  // numberParticle = value to assign
  // return value = reference on current quantum number
  NumberParticleQuantumNumber& operator = (int numberParticle);

  // Get number of particles
  //
  // return value = number of particles
  int GetNumberParticle ();

  // add two quantum numbers
  //
  // Q1 = first quantum number
  // Q2 = second quantum number
  // return value = sum of the two quantum numbers
  friend NumberParticleQuantumNumber operator + (const NumberParticleQuantumNumber& Q1, const NumberParticleQuantumNumber& Q2);

  // substract two quantum numbers
  //
  // Q1 = first quantum number
  // Q2 = quantum number to substract
  // return value = sum of the two quantum numbers
  friend NumberParticleQuantumNumber operator - (const NumberParticleQuantumNumber& Q1, const NumberParticleQuantumNumber& Q2);

  // test if two quantum numbers are identical
  //
  // Q1 = first quantum number
  // Q2 = second quantum number
  // return value = true if quantum numbers are identical
  friend bool operator == (const NumberParticleQuantumNumber& Q1, const NumberParticleQuantumNumber& Q2);

  // test if two quantum numbers are different
  //
  // Q1 = first quantum number
  // Q2 = second quantum number
  // return value = true if quantum numbers are different
  friend bool operator != (const NumberParticleQuantumNumber& Q1, const NumberParticleQuantumNumber& Q2);

  // print quantum number
  //
  // Str = reference on current output stream 
  // Q = quantum number to print
  // return value = reference on current output stream 
  friend ostream& operator << (ostream& Str, const NumberParticleQuantumNumber& Q);

};

#endif


