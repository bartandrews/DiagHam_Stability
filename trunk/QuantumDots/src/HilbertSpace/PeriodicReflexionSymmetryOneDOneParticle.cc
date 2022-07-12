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


#include "config.h"
#include "HilbertSpace/PeriodicReflexionSymmetryOneDOneParticle.h"


using std::cout;
using std::endl;


// default constructor
//

PeriodicReflexionSymmetryOneDOneParticle::PeriodicReflexionSymmetryOneDOneParticle ()
{
}

// basic constructor
//
// nbrState = wave function basis dimension without symmetry redundancy
// even = true -> even function, else -> odd function

PeriodicReflexionSymmetryOneDOneParticle::PeriodicReflexionSymmetryOneDOneParticle (int nbrState, bool even)
{
 if ((nbrState % 2) == 0)
   cout << "There is one dimension which is not correct: " << nbrState << ". It will be substract by 1" << endl;
   
  this->Even = even;
  if (this->Even)
    {
      this->NbrState = (nbrState / 2) + 1;
      this->LowerImpulsion = 0;
    }
  else
    {
      this->NbrState = nbrState / 2;
      this->LowerImpulsion = 1;
    }
  this->HilbertSpaceDimension = this->NbrState;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}
  
// copy constructor
//
// space = reference on Hilbert space to copy

PeriodicReflexionSymmetryOneDOneParticle::PeriodicReflexionSymmetryOneDOneParticle (const PeriodicReflexionSymmetryOneDOneParticle& space)
{
  this->NbrState = space.NbrState;
  this->LowerImpulsion = space.LowerImpulsion;
  this->Even = space.Even;
}

// destructor
//

PeriodicReflexionSymmetryOneDOneParticle::~PeriodicReflexionSymmetryOneDOneParticle ()
{
}

// assignement
//
// space = reference on Hilbert space to assign
// return value = reference on current Hilbert space

PeriodicReflexionSymmetryOneDOneParticle& PeriodicReflexionSymmetryOneDOneParticle::operator = (const PeriodicReflexionSymmetryOneDOneParticle& space)
{
  this->NbrState = space.NbrState;
  this->LowerImpulsion = space.LowerImpulsion;
  this->Even = space.Even;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* PeriodicReflexionSymmetryOneDOneParticle::Clone()
{
  return new PeriodicReflexionSymmetryOneDOneParticle (*this);
}

