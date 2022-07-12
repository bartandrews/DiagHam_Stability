////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2003-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//class of hilbert space of one 3d periodic box particle with reflexion symmetry of X axis //
//                                                                            //
//                        last modification : 12/08/2004                      //
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


#ifndef PERIODICXREFLEXIONYZPERIODICTHREEDONEPARTICLE_H
#define PERIODICXREFLEXIONYZPERIODICTHREEDONEPARTICLE_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "HilbertSpace/PeriodicOneDOneParticle.h"
#include "HilbertSpace/PeriodicReflexionSymmetryOneDOneParticle.h"
#include "HilbertSpace/ThreeDOneParticle.h"


class PeriodicXReflexionYZPeriodicThreeDOneParticle : public ThreeDOneParticle
{

 protected:

 public:

  // default constructor
  //
  PeriodicXReflexionYZPeriodicThreeDOneParticle ();

  // constructor
  //
  // nbrStateX = wave function basis dimension in the x direction without symmetry redundancy
  // evenX = true if the wave functions in X is even, else odd
  // nbrStateY = wave function basis dimension in the y direction
  // lowY = lower impulsion in Y direction (in unit of 2 * Pi / Lz)
  // nbrStateZ = wave function basis dimension in the z direction
  // lowZ = lower impulsion in Z direction (in unit of 2 * Pi / Lz)
  PeriodicXReflexionYZPeriodicThreeDOneParticle (int nbrStateX, bool evenX, int nbrStateY, int lowY, int nbrStateZ, int lowZ);

  // copy constructor
  //
  // space = reference on Hilbert space to copy
  PeriodicXReflexionYZPeriodicThreeDOneParticle (const PeriodicXReflexionYZPeriodicThreeDOneParticle& space);

  // destructor
  //
  virtual ~PeriodicXReflexionYZPeriodicThreeDOneParticle();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // assignement
  //
  // space = reference on Hilbert space to assign
  // return value = reference on current Hilbert space
  PeriodicXReflexionYZPeriodicThreeDOneParticle& operator = (const PeriodicXReflexionYZPeriodicThreeDOneParticle& space);

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

inline int PeriodicXReflexionYZPeriodicThreeDOneParticle::GetLowerImpulsionX ()
{
  return ((PeriodicReflexionSymmetryOneDOneParticle*) this->StateX)->GetLowerImpulsion ();
}

// get lower impulsion in Y
//
// return value = lower value in unit of 2 * Pi / L

inline int PeriodicXReflexionYZPeriodicThreeDOneParticle::GetLowerImpulsionY ()
{
  return ((PeriodicOneDOneParticle*) this->StateY)->GetLowerImpulsion ();
}

// get lower impulsion in Z
//
// return value = lower value in unit of 2 * Pi / L

inline int PeriodicXReflexionYZPeriodicThreeDOneParticle::GetLowerImpulsionZ ()
{
  return ((PeriodicOneDOneParticle*) this->StateZ)->GetLowerImpulsion ();
}

// get the wave function basis in X direction
//
// return = pointer to 1D one particle basis

inline OneDOneParticle* PeriodicXReflexionYZPeriodicThreeDOneParticle::GetStateX ()
{
  OneDOneParticle* stateX = (PeriodicReflexionSymmetryOneDOneParticle*) this->StateX->Clone ();
  return stateX;
}

// get the wave function basis in X direction
//
// return = pointer to 1D one particle basis

inline OneDOneParticle* PeriodicXReflexionYZPeriodicThreeDOneParticle::GetStateY ()
{
  OneDOneParticle* stateY = (PeriodicOneDOneParticle*) this->StateY->Clone ();
  return stateY;
}

// get the wave function basis in Z direction
//
// return = pointer to 1D one particle basis
inline OneDOneParticle* PeriodicXReflexionYZPeriodicThreeDOneParticle::GetStateZ ()
{
  OneDOneParticle* stateZ = (PeriodicOneDOneParticle*) this->StateZ->Clone ();
  return stateZ;
}

#endif
