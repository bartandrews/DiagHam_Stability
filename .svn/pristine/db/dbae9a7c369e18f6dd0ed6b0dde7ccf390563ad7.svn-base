////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of bosons on sphere for system size such that            //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                        without a fixed total Lz value                      //
//                                                                            //
//                        last modification : 19/04/2010                      //
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


#ifndef BOSONONSPHEREFULLSHORT_H
#define BOSONONSPHEREFULLSHORT_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/FermionOnSphereFull.h"
#include "Matrix/RealMatrix.h"

#include <iostream>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


class BosonOnSphereFullShort :  public BosonOnSphereShort
{

 protected:


 public:

  // default constructor
  //
  BosonOnSphereFullShort ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // lzMax = maximum Lz value reached by a boson
  BosonOnSphereFullShort (int nbrBosons, int lzMax);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereFullShort(const BosonOnSphereFullShort& bosons);

  // destructor
  //
  virtual ~BosonOnSphereFullShort ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereFullShort& operator = (const BosonOnSphereFullShort& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdA (int index, int m);

  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdA (long index, int m);

  // get Lz component of a component
  //
  // j = index of the component in Hilbert space
  // return value = twice the Lz component
  virtual int GetLzValue(int j = 0);

  // evaluate a density matrix of a subsystem of the whole system described |left><right|. The density matrix is only evaluated for a fixed number of particles
  // 
  // densityMatrix = reference on the temporary storage for the reduced density matrix
  // leftSpace = Hilbert space associated to the left state 
  // leftState = reference on the left state
  // rightSpace = Hilbert space associated to the right state 
  // rightState = reference on the right state
  // return value = reference on the reduced density matrix of the subsytem 
  virtual RealMatrix& EvaluatePartialDensityMatrix (RealMatrix& densityMatrix, BosonOnSphereShort* leftSpace, RealVector& leftState, BosonOnSphereShort* rightSpace, RealVector& rightState);


 protected:

};

// get Lz component of a component
//
// j = index of the component in Hilbert space
// return value = twice the Lz component

inline int BosonOnSphereFullShort::GetLzValue(int j)
{
  return ((FermionOnSphereFull*) this->FermionBasis)->TotalLzValues[j];
}


#endif


