////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
// class of hilbert space of 3d periodic box particle with X reflexion symmetry//
//                                                                            //
//                        last modification : 05/03/2004                      //
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


#ifndef INPLANEREFLEXIONSYMMETRICPERIODIC3DONEPARTICLE_H
#define INPLANEREFLEXIONSYMMETRICPERIODIC3DONEPARTICLE_H

#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "HilbertSpace/Periodic3DOneParticle.h"


class InPlaneReflexionSymmetricPeriodic3DOneParticle : public Periodic3DOneParticle
{

 protected:

 public:

  // constructor
  //
  // pairX = true if the basis function in x direction is even, false if odd
  // maxX = maximal wave function basis dimension in the x direction
  // nbrStateY = wave function basis dimension in the y direction
  // lowY = lower bound of basis dimension in the y direction
  // nbrStateZ = wave function basis dimension in the z direction
  // lowZ = lower bound of basis dimension in the z direction
  InPlaneReflexionSymmetricPeriodic3DOneParticle (bool pairX, int maxX, int nbrStateY, int lowY, int nbrStateZ, int lowZ);

  // copy constructor
  //
  // space = reference on Hilbert space to copy
  InPlaneReflexionSymmetricPeriodic3DOneParticle(const  InPlaneReflexionSymmetricPeriodic3DOneParticle& space);

  // destructor
  //
  ~InPlaneReflexionSymmetricPeriodic3DOneParticle();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  // AbstractHilbertSpace* Clone();

  // assignement
  //
  // space = reference on Hilbert space to assign
  // return value = reference on current Hilbert space
   InPlaneReflexionSymmetricPeriodic3DOneParticle& operator = (const  InPlaneReflexionSymmetricPeriodic3DOneParticle& space);

  // get sinus wave function basis dimension in the x direction
  //
  // return value = wave function basis dimension in the x direction
  virtual int GetNbrSinusStateX();
  
  // get cosinus wave function basis dimension in the x direction
  //
  // return value = wave function basis dimension in the x direction
  virtual int GetNbrCosinusStateX(); 
  
  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  // List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  // AbstractQuantumNumber* GetQuantumNumber (int index);
  
  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  //AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

};

// get sinus wave function basis dimension in the x direction
//
// return value = wave function basis dimension in the x direction

inline int InPlaneReflexionSymmetricPeriodic3DOneParticle::GetNbrSinusStateX()
{
  return -this->LowerImpulsionX;
}
  
// get cosinus wave function basis dimension in the x direction
//
// return value = wave function basis dimension in the x direction

inline int InPlaneReflexionSymmetricPeriodic3DOneParticle::GetNbrCosinusStateX()
{
   return -this->LowerImpulsionX + 1;
} 

#endif
