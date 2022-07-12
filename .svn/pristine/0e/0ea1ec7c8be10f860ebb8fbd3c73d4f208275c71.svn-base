////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of fermions on sphere truncated to a given P level           //
//                                                                            //
//                        last modification : 28/09/2012                      //
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


#ifndef FERMIONONSPHEREPTRUNCATED_H
#define FERMIONONSPHEREPTRUNCATED_H


#include "config.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/SparseComplexMatrix.h"


#include <iostream>


class FermionOnSpherePTruncated :  public FermionOnSphere
{

 protected:

  // root configuration in the occupation basis
  unsigned long ReferenceState;
  // root configuration in the monomial basis
  int* ReferenceStateMonomialBasis;
  
  // truncation level
  int PLevel;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a fermion
  // pLevel = truncation level
  // referenceState = array that describes the root configuration
  // memory = amount of memory granted for precalculations
  FermionOnSpherePTruncated (int nbrFermions, int& totalLz, int lzMax, int pLevel,
			     int* referenceState, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSpherePTruncated(const FermionOnSpherePTruncated& fermions);

  // copy constructor, preserving only some specific states 
  //
  // fermions = reference on the hilbert space to copy to copy
  // nbrPreservedStates = number of preserved states
  // preservedStates = array of flags that indicates if the corresponding state has to be preserved 
  //                   (dimension of the array should the one of the original Hilbert space)
  FermionOnSpherePTruncated(const FermionOnSpherePTruncated& fermions, long nbrPreservedStates, bool* preservedStates);

  // destructor
  //
  virtual ~FermionOnSpherePTruncated ();

  // assignment (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSpherePTruncated& operator = (const FermionOnSpherePTruncated& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // create a state from its MPS description
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

  // create a state from its MPS description, inclusing additional quasihole matrices
  //
  // bMatrices = array that gives the B matrices 
  // state = reference to vector that will contain the state description
  // mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
  // mPSColumnIndex = column index of the MPS element that has to be evaluated
  // memory = amount of memory that can be use to precompute matrix multiplications  
  // initialIndex = initial index to compute
  // nbrComponents = number of components to compute
  void CreateStateFromMPSDescription (SparseRealMatrix* bMatrices, SparseComplexMatrix* quasiholeBMatrices, int nbrQuasiholeBMatrices,
				      ComplexVector& state, int mPSRowIndex, int mPSColumnIndex, 
				      long memory = 0l, long initialIndex = 0l, long nbrComponents = 0l);

  // convert a given state from truncated to Haldane basis
  //
  // state = reference on the vector to convert
  // haldaneBasis = reference on the Haldane basis to use
  // return value = converted vector
  RealVector ConvertToHaldaneBasis(RealVector& state, FermionOnSphereHaldaneBasis& haldaneBasis);

  // convert a given state from truncated to Haldane basis
  //
  // state = reference on the vector to convert
  // haldaneBasis = reference on the Haldane basis to use
  // return value = converted vector
  ComplexVector ConvertToHaldaneBasis(ComplexVector& state, FermionOnSphereHaldaneBasis& haldaneBasis);

  // create the Jack polynomial decomposition corresponding to the root partition
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized basis
  virtual RealVector& GenerateJackPolynomial(RealVector& jack, double alpha);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);

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


