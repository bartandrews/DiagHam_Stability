////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//              class of hilbert space of one 1d periodic box particle        //
//                                                                            //
//                        last modification : 05/07/2004                      //
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


#ifndef PERIODIC1DONEPARTICLE_H
#define PERIODIC1DONEPARTICLE_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"


class Periodic1DOneParticle : public AbstractHilbertSpace
{

 protected:

  // wave function basis dimension
  int NbrState;

  // lower impulsion
  int LowerImpulsion;

 public:

  // default constructor
  Periodic1DOneParticle();

  // constructor
  //
  // nbrState = wave function basis dimension
  // low = lower impulsion
  Periodic1DOneParticle(int nbrState, int low);

  // copy constructor
  //
  // space = reference on Hilbert space to copy
  Periodic1DOneParticle(const Periodic1DOneParticle& space);

  // destructor
  //
  ~Periodic1DOneParticle();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // assignement
  //
  // space = reference on Hilbert space to assign
  // return value = reference on current Hilbert space
  Periodic1DOneParticle& operator = (const Periodic1DOneParticle& space);

  // get wave function basis dimension
  //
  // return value = wave function basis dimension in the x direction
  virtual int GetNbrState();

  // get lower impulsion
  //
  // return value = lower value in unit of 2 * Pi / L
  virtual int GetLowerImpulsion();

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
					 SubspaceSpaceConverter& converter);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

};

// get wave function basis dimension
//
// return value = wave function basis dimension 

inline int Periodic1DOneParticle::GetNbrState()
{
  return this->NbrState;
}

// get lower impulsion 
//
// return value = lower value in unit of 2 * Pi / L

inline int Periodic1DOneParticle::GetLowerImpulsion()
{
  return this->LowerImpulsion;
}


#endif
