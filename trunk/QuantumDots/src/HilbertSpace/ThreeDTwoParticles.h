////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2003-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//                     class of hilbert space of two 3d particles             //
//                                                                            //
//                        last modification : 13/10/2004                      //
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


#ifndef THREEDTWOPARTICLES_H
#define THREEDTWOPARTICLES_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "HilbertSpace/ThreeDOneParticle.h"


class ThreeDTwoParticles : public AbstractHilbertSpace
{

 protected:

  // Hilbert space of the first particle
  ThreeDOneParticle* FirstParticle;

  // Hilbert space of the second particle
  ThreeDOneParticle* SecondParticle;

 public:

  // default constructor
  //
  ThreeDTwoParticles ();

  // constructor from two abstract particles
  //
  // firstParticle = pointer to the first abstract particle
  // secondParticle = pointer to the second abstract particle
  ThreeDTwoParticles (ThreeDOneParticle* firstParticle, ThreeDOneParticle* secondParticle);

  // copy constructor
  //
  // space = reference on Hilbert space to copy
  ThreeDTwoParticles (const ThreeDTwoParticles& space);

  // destructor
  //
  virtual ~ThreeDTwoParticles();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // assignement
  //
  // space = reference on Hilbert space to assign
  // return value = reference on current Hilbert space
  ThreeDTwoParticles& operator = (const ThreeDTwoParticles& space);

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

  // get the Hilbert space description of the first particle
  //
  // return = pointer to the Hilbert space
  virtual ThreeDOneParticle* GetFirstParticleSpace ();

  // get the Hilbert space description of the second particle
  //
  // return = pointer to the Hilbert space
  virtual ThreeDOneParticle* GetSecondParticleSpace ();

};

// get the Hilbert space description of the first particle
//
// return = pointer to the Hilbert space

inline ThreeDOneParticle* ThreeDTwoParticles::GetFirstParticleSpace ()
{
  ThreeDOneParticle* firstParticle = (ThreeDOneParticle*) this->FirstParticle->Clone ();
  return firstParticle;
}

// get the Hilbert space description of the second particle
//
// return = pointer to the Hilbert space

inline ThreeDOneParticle* ThreeDTwoParticles::GetSecondParticleSpace ()
{
  ThreeDOneParticle* secondParticle = (ThreeDOneParticle*) this->SecondParticle->Clone ();
  return secondParticle;
}

#endif
