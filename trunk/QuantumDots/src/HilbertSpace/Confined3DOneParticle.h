////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of hilbert space of one 3d box confined particle        //
//                                                                            //
//                        last modification : 26/02/2003                      //
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


#ifndef CONFINED3DONEPARTICLE_H
#define CONFINED3DONEPARTICLE_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"


class Confined3DOneParticle : public AbstractHilbertSpace
{

 protected:

  // wave function basis dimension in the x direction
  int NbrStateX;
  // wave function basis dimension in the y direction
  int NbrStateY;
  // wave function basis dimension in the z direction
  int NbrStateZ;  

 public:

  // constructor
  //
  // nbrStateX = wave function basis dimension in the x direction
  // nbrStateY = wave function basis dimension in the y direction
  // nbrStateZ = wave function basis dimension in the z direction
  Confined3DOneParticle(int nbrStateX, int nbrStateY, int nbrStateZ);

  // copy constructor
  //
  // space = reference on Hilbert space to copy
  Confined3DOneParticle(const Confined3DOneParticle& space);

  // destructor
  //
  ~Confined3DOneParticle();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // assignement
  //
  // space = reference on Hilbert space to assign
  // return value = reference on current Hilbert space
  Confined3DOneParticle& operator = (const Confined3DOneParticle& space);

  // get wave function basis dimension in the x direction
  //
  // return value = wave function basis dimension in the x direction
  int GetNbrStateX();

  // get wave function basis dimension in the y direction
  //
  // return value = wave function basis dimension in the y direction
  int GetNbrStateY();

  // get wave function basis dimension in the z direction
  //
  // return value = wave function basis dimension in the z direction
  int GetNbrStateZ();

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

// get wave function basis dimension in the x direction
//
// return value = wave function basis dimension in the x direction

inline int Confined3DOneParticle::GetNbrStateX()
{
  return this->NbrStateX;
}

// get wave function basis dimension in the y direction
//
// return value = wave function basis dimension in the y direction

inline int Confined3DOneParticle::GetNbrStateY()
{
  return this->NbrStateY;
}

// get wave function basis dimension in the z direction
//
// return value = wave function basis dimension in the z direction

inline int Confined3DOneParticle::GetNbrStateZ()
{
  return this->NbrStateZ;
}

#endif


