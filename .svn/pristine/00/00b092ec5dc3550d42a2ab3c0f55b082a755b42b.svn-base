////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of bosons on sphere truncated to a given P level            //
//   such that LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)   //
//                                                                            //
//                        last modification : 17/10/2012                      //
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


#ifndef BOSONONSPHEREPTRUNCATED_H
#define BOSONONSPHEREPTRUNCATED_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "Matrix/SparseRealMatrix.h"


#include <iostream>


class BosonOnSpherePTruncated :  public BosonOnSphereShort
{

 protected:

  // root configuration in the occupation basis
  unsigned long ReferenceState;
  // root configuration in the monomial basis
  int* ReferenceStateMonomialBasis;
  
  // truncation level
  int PLevel;

  // maximum occupation for a single orbital
  int MaximumOccupation;

 public:

  // basic constructor
  // 
  // nbrBosons = number of fermions
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a boson
  // pLevel = truncation level
  // maximumOccupation = maximum occupation for a single orbital
  // referenceState = array that describes the root configuration
  BosonOnSpherePTruncated (int nbrBosons, int& totalLz, int lzMax, int pLevel,
			   int maximumOccupation, int* referenceState);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSpherePTruncated(const BosonOnSpherePTruncated& bosons);

  // destructor
  //
  virtual ~BosonOnSpherePTruncated ();

  // assignment (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSpherePTruncated& operator = (const BosonOnSpherePTruncated& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // create a state from its MPS description, assuming the resulting state is real
  //
  // bMatrices = array that gives the B matrices 
  // state = reference to vector that will contain the state description
  // mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
  // mPSColumnIndex = column index of the MPS element that has to be evaluated
  // memory = amount of memory that can be use to precompute matrix multiplications  
  // initialIndex = initial index to compute
  // nbrComponents = number of components to compute
  virtual void CreateStateFromMPSDescription (SparseRealMatrix* bMatrices, RealVector& state, int mPSRowIndex, int mPSColumnIndex,
					      long memory = 0l, long initialIndex = 0l, long nbrComponents = 0l);
  
 protected:

};


#endif


