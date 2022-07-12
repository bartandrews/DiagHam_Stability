////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2003-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//                     class of hilbert space of one 3d particle              //
//                                                                            //
//                        last modification : 11/08/2004                      //
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


#ifndef THREEDONEPARTICLE_H
#define THREEDONEPARTICLE_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "HilbertSpace/OneDOneParticle.h"


class ThreeDOneParticle : public AbstractHilbertSpace
{

 protected:

  // wave function basis dimension in the x direction
  OneDOneParticle* StateX;
  // wave function basis dimension in the y direction
  OneDOneParticle* StateY;
  // wave function basis dimension in the z direction
  OneDOneParticle* StateZ;

 public:

  // default constructor
  //
  ThreeDOneParticle ();

  // constructor
  //
  // nbrStateX = wave function basis dimension in the x direction
  // nbrStateY = wave function basis dimension in the y direction
  // nbrStateZ = wave function basis dimension in the z direction
  ThreeDOneParticle (int nbrStateX, int nbrStateY, int nbrStateZ);

  // copy constructor
  //
  // space = reference on Hilbert space to copy
  ThreeDOneParticle (const ThreeDOneParticle& space);

  // destructor
  //
  virtual ~ThreeDOneParticle();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // assignement
  //
  // space = reference on Hilbert space to assign
  // return value = reference on current Hilbert space
  ThreeDOneParticle& operator = (const ThreeDOneParticle& space);

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

  // get the wave function basis in X direction
  //
  // return = pointer to 1D one particle basis
  virtual OneDOneParticle* GetStateX ();

  // get the wave function basis in Y direction
  //
  // return = pointer to 1D one particle basis
  virtual OneDOneParticle* GetStateY ();

  // get the wave function basis in Z direction
  //
  // return = pointer to 1D one particle basis
  virtual OneDOneParticle* GetStateZ ();

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

inline int ThreeDOneParticle::GetNbrStateX ()
{
  return this->StateX->GetNbrState ();
}

// get wave function basis dimension in the y direction
//
// return value = wave function basis dimension in the y direction

inline int ThreeDOneParticle::GetNbrStateY ()
{
  return this->StateY->GetNbrState ();
}

// get wave function basis dimension in the z direction
//
// return value = wave function basis dimension in the z direction

inline int ThreeDOneParticle::GetNbrStateZ ()
{
  return this->StateZ->GetNbrState ();
}

// get the wave function basis in X direction
//
// return = pointer to 1D one particle basis

inline OneDOneParticle* ThreeDOneParticle::GetStateX ()
{
  OneDOneParticle* stateX = (OneDOneParticle*) this->StateX->Clone ();
  return stateX;
}

// get the wave function basis in Y direction
//
// return = pointer to 1D one particle basis

inline OneDOneParticle* ThreeDOneParticle::GetStateY ()
{
  OneDOneParticle* stateY = (OneDOneParticle*) this->StateY->Clone ();
  return stateY;
}

// get the wave function basis in Z direction
//
// return = pointer to 1D one particle basis
inline OneDOneParticle* ThreeDOneParticle::GetStateZ ()
{
  OneDOneParticle* stateZ = (OneDOneParticle*) this->StateZ->Clone ();
  return stateZ;
}


#endif



