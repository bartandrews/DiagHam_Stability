////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//class of hilbert space of one 1d periodic box particle with reflexion symmetry//
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


#ifndef PERIODICREFLEXIONSYMMETRYONEDONEPARTICLE_H
#define PERIODICREFLEXIONSYMMETRYONEDONEPARTICLE_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "HilbertSpace/OneDOneParticle.h"


class PeriodicReflexionSymmetryOneDOneParticle : public OneDOneParticle
{

 protected:

  // Even = true -> even function, else -> odd function
  bool Even;
  // lower impulsion, in unit of 2 * Pi / L
  int LowerImpulsion;

 public:

  // default constructor
  //
  PeriodicReflexionSymmetryOneDOneParticle ();

  // basic constructor
  //
  // nbrState = wave function basis dimension without symmetry redundancy
  // even = true -> even function, else -> odd function
  PeriodicReflexionSymmetryOneDOneParticle (int nbrState, bool even);

  // copy constructor
  //
  // space = reference on Hilbert space to copy
  PeriodicReflexionSymmetryOneDOneParticle (const PeriodicReflexionSymmetryOneDOneParticle& space);

  // destructor
  //
  virtual ~PeriodicReflexionSymmetryOneDOneParticle ();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // assignement
  //
  // space = reference on Hilbert space to assign
  // return value = reference on current Hilbert space
  PeriodicReflexionSymmetryOneDOneParticle& operator = (const PeriodicReflexionSymmetryOneDOneParticle& space);

  // determine whether the function is even or odd
  //
  bool IsEven ();

  // get lower impulsion
  //
  // return value = lower value in unit of 2 * Pi / L
  virtual int GetLowerImpulsion ();

};

// determine whether the function is even or odd
//
inline bool PeriodicReflexionSymmetryOneDOneParticle::IsEven ()
{
  return this->Even;
}

// get lower impulsion 
//
// return value = lower value in unit of 2 * Pi / L

inline int PeriodicReflexionSymmetryOneDOneParticle::GetLowerImpulsion ()
{
  return this->LowerImpulsion;
}

#endif
