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

  // add a quantum nunber to the current one
  //
  // number = quantum number to add 
  // return value = reference to the current quantum number
  AbstractQuantumNumber& operator += (const AbstractQuantumNumber& number);

  // substract a quantum nunber to the current one
  //
  // number = quantum number to add 
  // return value = reference to the current quantum number
  AbstractQuantumNumber& operator -= (const AbstractQuantumNumber& number);

  // test if two quantum numbers are identical
  //
  // number = quantum number to compare to the current one
  bool IsEqual (const AbstractQuantumNumber& number);

  // test if two quantum numbers are different
  //
  // number = quantum number to compare to the current one
  bool IsDifferent (const AbstractQuantumNumber& number);

  // print quantum number
  //
  // str = reference on current output stream 
  // return value = reference on current output stream 
  ostream& PrintQuantumNumber (ostream& str);

};

#endif


