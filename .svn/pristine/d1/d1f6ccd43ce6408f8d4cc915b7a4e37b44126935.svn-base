////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of fermions on sphere truncated to a given P level           //
//      with  LzMax up to 127 (for systems with 128 bit integer support)      //
//           or 63 (on 32 bit systems without 128 bit integer support)        //
//                                                                            //
//                        last modification : 02/10/2012                      //
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


#ifndef FERMIONONSPHEREPTRUNCATEDLONG_H
#define FERMIONONSPHEREPTRUNCATEDLONG_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereLong.h"

#include <iostream>



class FermionOnSpherePTruncatedLong :  public FermionOnSphereLong
{

 protected:

  // root configuration in the occupation basis
  unsigned long ReferenceState;
  // root configuration in the monomial basis
  int* ReferenceStateMonomialBasis;
  
  // truncation level
  int PLevel;

 public:

  // default constuctor
  //
  FermionOnSpherePTruncatedLong();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a fermion
  // pLevel = truncation level
  // referenceState = array that describes the root configuration
  // memory = amount of memory granted for precalculations
  FermionOnSpherePTruncatedLong (int nbrFermions, int& totalLz, int lzMax, int pLevel,
				 int* referenceState, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSpherePTruncatedLong(const FermionOnSpherePTruncatedLong& fermions);

  // destructor
  //
  virtual ~FermionOnSpherePTruncatedLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSpherePTruncatedLong& operator = (const FermionOnSpherePTruncatedLong& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // level = current level for truncation
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int level);

  // evaluate Hilbert space dimension with shifted values for lzMax and totalLz
  //
  // nbrFermions = number of fermions
  // lzMax = two times momentum maximum value for a fermion plus one 
  // totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
  // level = current level for truncation
  // return value = Hilbert space dimension  
  long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int level);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // currentLzMax = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // level = current level
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, int level, long pos);

};


#endif


