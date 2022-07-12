////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of bosons on disk for system size such that             //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 08/07/2008                      //
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


#ifndef BOSONONDISKSHORT_H
#define BOSONONDISKSHORT_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereShort.h"


class BosonOnDiskShort :  public BosonOnSphereShort
{
  friend class BosonOnDiskHaldaneBasisShort;

 public:
  
  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = momentum total value
  // lzMax = maximum angular momentum that a single particle can reach (negative if it has to be deduced from nbrBosons and totalLz)
  BosonOnDiskShort (int nbrBosons, int totalLz, int lzMax = -1);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnDiskShort(const BosonOnDiskShort& bosons);

  // copy constructor, preversing only some specific states 
  // 
  // bosons = reference on the hilbert space to copy to copy
  // nbrPreservedStates = number of preserved states
  // preservedStates = array of flags that indicates if the corresponding state has to be preserved 
  //                   (dimension of the array should the one of the original Hilbert space)
  BosonOnDiskShort (const BosonOnDiskShort& bosons, int nbrPreservedStates, bool* preservedStates);

  // destructor
  //
  ~BosonOnDiskShort ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnDiskShort& operator = (const BosonOnDiskShort& bosons);

  // convert a state such that its components are now expressed in the unnormalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has to be normalized to 1
  // symmetryFactor = if true also remove the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertToUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);

  // convert a state such that its components are now expressed in the normalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has been normalized to 1
  // symmetryFactor = if true also add the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertFromUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);

  // convert a state such that its components are now expressed in the normalized basis, shifting all orbitals
  //
  // state = reference to the state to convert
  // shift = shift to apply to each orbitals
  // reference = set which component has been normalized to 1
  // return value = converted state
  virtual RealVector& ShiftedConvertFromUnnormalizedMonomial(RealVector& state, int shift, long reference = 0);
	
  // evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using real space partition. The entanglement matrix is only evaluated in a given Lz sector.
  // and computed from precalculated particle entanglement matrix
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // radius = radius of the A disk
  // shift = shift to apply to each orbitals
  // entanglementMatrix = reference on the entanglement matrix (will be overwritten)
  // return value = reference on the entanglement matrix
  RealMatrix& EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix (int nbrBosonSector, int lzSector, double radius, int shift, RealMatrix& entanglementMatrix);
  
 protected:

  // core part of the convertion of a state such that its components are now expressed in the normalized basis
  //
  // state = reference to the state to convert
  // sqrtCoefficients = array that contains the normalization coefficients
  // invSqrtCoefficients = array that contains the inverts of the normalization coefficients
  // reference = set which component has been normalized to 1
  // symmetryFactor = if true also add the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertFromUnnormalizedMonomialCore(RealVector& state, double* sqrtCoefficients, double* invSqrtCoefficients, long reference, bool symmetryFactor);

};

#endif


