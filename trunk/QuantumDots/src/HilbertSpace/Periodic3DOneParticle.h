////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of hilbert space of one 3d periodic box particle        //
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


#ifndef PERIODIC3DONEPARTICLE_H
#define PERIODIC3DONEPARTICLE_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"


class Periodic3DOneParticle : public AbstractHilbertSpace
{

 protected:

  // wave function basis dimension in the x direction
  int NbrStateX;

  // lower impulsion in X
  int LowerImpulsionX;

  // wave function basis dimension in the y direction
  int NbrStateY;

  // lower impulsion in Y
  int LowerImpulsionY;

  // wave function basis dimension in the z direction
  int NbrStateZ;  

  // lower impulsion in Z
  int LowerImpulsionZ;

 public:

  // default constructor
  Periodic3DOneParticle();

  // constructor
  //
  // nbrStateX = wave function basis dimension in the x direction
  // nbrStateY = wave function basis dimension in the y direction
  // nbrStateZ = wave function basis dimension in the z direction
  Periodic3DOneParticle(int nbrStateX, int lowX, int nbrStateY, int lowY, int nbrStateZ, int lowZ);

  // copy constructor
  //
  // space = reference on Hilbert space to copy
  Periodic3DOneParticle(const Periodic3DOneParticle& space);

  // destructor
  //
  ~Periodic3DOneParticle();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // assignement
  //
  // space = reference on Hilbert space to assign
  // return value = reference on current Hilbert space
  Periodic3DOneParticle& operator = (const Periodic3DOneParticle& space);

  // get wave function basis dimension in the x direction
  //
  // return value = wave function basis dimension in the x direction
  virtual int GetNbrStateX();

  // get lower impulsion in X
  //
  // return value = lower value in unit of 2 * Pi / L
  virtual int GetLowerImpulsionX();

  // get wave function basis dimension in the y direction
  //
  // return value = wave function basis dimension in the y direction
  virtual int GetNbrStateY();

  // get lower impulsion in Y
  //
  // return value = lower value in unit of 2 * Pi / L
  virtual int GetLowerImpulsionY();

  // get wave function basis dimension in the z direction
  //
  // return value = wave function basis dimension in the z direction
  virtual int GetNbrStateZ();

  // get lower impulsion in Z
  //
  // return value = lower value in unit of 2 * Pi / L
  virtual int GetLowerImpulsionZ();

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

inline int Periodic3DOneParticle::GetNbrStateX()
{
  return this->NbrStateX;
}

// get lower impulsion in X
//
// return value = lower value in unit of 2 * Pi / L

inline int Periodic3DOneParticle::GetLowerImpulsionX()
{
  return this->LowerImpulsionX;
}

// get wave function basis dimension in the y direction
//
// return value = wave function basis dimension in the y direction

inline int Periodic3DOneParticle::GetNbrStateY()
{
  return this->NbrStateY;
}

// get lower impulsion in Y
//
// return value = lower value in unit of 2 * Pi / L

inline int Periodic3DOneParticle::GetLowerImpulsionY()
{
  return this->LowerImpulsionY;
}

// get wave function basis dimension in the z direction
//
// return value = wave function basis dimension in the z direction

inline int Periodic3DOneParticle::GetNbrStateZ()
{
  return this->NbrStateZ;
}

// get lower impulsion in Z
//
// return value = lower value in unit of 2 * Pi / L

inline int Periodic3DOneParticle::GetLowerImpulsionZ()
{
  return this->LowerImpulsionZ;
}

#endif


