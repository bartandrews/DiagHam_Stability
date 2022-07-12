////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2003-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//             class of hilbert space of one periodic 3d particle             //
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


#ifndef PERIODICTHREEDONEPARTICLE_H
#define PERIODICTHREEDONEPARTICLE_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "HilbertSpace/PeriodicOneDOneParticle.h"
#include "HilbertSpace/ThreeDOneParticle.h"


class PeriodicThreeDOneParticle : public ThreeDOneParticle
{

 protected:

 public:

  // default constructor
  //
  PeriodicThreeDOneParticle ();

  // constructor
  //
  // nbrStateX = wave function basis dimension in the x direction
  // lowX = lower impulsion in X direction (in unit of 2 * Pi / Lx)
  // nbrStateY = wave function basis dimension in the y direction
  // lowY = lower impulsion in Y direction (in unit of 2 * Pi / Ly)
  // nbrStateZ = wave function basis dimension in the z direction
  // lowZ = lower impulsion in Z direction (in unit of 2 * Pi / Lz)
  PeriodicThreeDOneParticle (int nbrStateX, int lowX, int nbrStateY, int lowY, int nbrStateZ, int lowZ);

  // copy constructor
  //
  // space = reference on Hilbert space to copy
  PeriodicThreeDOneParticle (const PeriodicThreeDOneParticle& space);

  // destructor
  //
  virtual ~PeriodicThreeDOneParticle();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // assignement
  //
  // space = reference on Hilbert space to assign
  // return value = reference on current Hilbert space
  PeriodicThreeDOneParticle& operator = (const PeriodicThreeDOneParticle& space);

  // get lower impulsion in X
  //
  // return value = lower value in unit of 2 * Pi / L
  virtual int GetLowerImpulsionX ();

  // get lower impulsion in Y
  //
  // return value = lower value in unit of 2 * Pi / L
  virtual int GetLowerImpulsionY ();

  // get lower impulsion in Z
  //
  // return value = lower value in unit of 2 * Pi / L
  virtual int GetLowerImpulsionZ ();

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
};

// get lower impulsion in X
//
// return value = lower value in unit of 2 * Pi / L

inline int PeriodicThreeDOneParticle::GetLowerImpulsionX ()
{
  return ((PeriodicOneDOneParticle*) this->StateX)->GetLowerImpulsion ();
}

// get lower impulsion in Y
//
// return value = lower value in unit of 2 * Pi / L

inline int PeriodicThreeDOneParticle::GetLowerImpulsionY ()
{
  return ((PeriodicOneDOneParticle*) this->StateY)->GetLowerImpulsion ();
}

// get lower impulsion in Z
//
// return value = lower value in unit of 2 * Pi / L

inline int PeriodicThreeDOneParticle::GetLowerImpulsionZ ()
{
  return ((PeriodicOneDOneParticle*) this->StateZ)->GetLowerImpulsion ();
}

// get the wave function basis in X direction
//
// return = pointer to 1D one particle basis

inline OneDOneParticle* PeriodicThreeDOneParticle::GetStateX ()
{
  OneDOneParticle* stateX = (PeriodicOneDOneParticle*) this->StateX->Clone ();
  return stateX;
}

// get the wave function basis in X direction
//
// return = pointer to 1D one particle basis

inline OneDOneParticle* PeriodicThreeDOneParticle::GetStateY ()
{
  OneDOneParticle* stateY = (PeriodicOneDOneParticle*) this->StateY->Clone ();
  return stateY;
}

// get the wave function basis in Z direction
//
// return = pointer to 1D one particle basis
inline OneDOneParticle* PeriodicThreeDOneParticle::GetStateZ ()
{
  OneDOneParticle* stateZ = (PeriodicOneDOneParticle*) this->StateZ->Clone ();
  return stateZ;
}

#endif
