////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of bosons on disk  using the Haldane basis            //
//                            for system size such that                       //
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


#ifndef BOSONONDISKHALDANEHUGEBASISSHORT_H
#define BOSONONDISKHALDANEHUGEBASISSHORT_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"


class BosonOnDiskHaldaneHugeBasisShort :  public BosonOnSphereHaldaneHugeBasisShort
{


 public:

  // default constructor
  //
  BosonOnDiskHaldaneHugeBasisShort();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a boson
  // maxFileSize = maximum file size (in MBytes)
  // referenceState = array that describes the reference state to start from
  // memory = amount of memory granted for precalculations
  // referenceState = array that describes the reference state to start from
  // symmetricFlag = indicate if a symmetric basis has to be used (only available if totalLz = 0)
  // fullDimension = provide the full (i.e. without squeezing) Hilbert space dimension (0 if it has to be computed)
  BosonOnDiskHaldaneHugeBasisShort (int nbrBosons, int totalLz, int lzMax, unsigned long maxFileSize, int* referenceState, unsigned long memory, bool symmetricFlag, long fullDimension);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memoryHilbert = amount of memory granted to store the Hilbert space (in Mbytes)
  BosonOnDiskHaldaneHugeBasisShort (char* fileName, unsigned long memoryHilbert);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnDiskHaldaneHugeBasisShort(const BosonOnDiskHaldaneHugeBasisShort& bosons);

  // destructor
  //
  ~BosonOnDiskHaldaneHugeBasisShort ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnDiskHaldaneHugeBasisShort& operator = (const BosonOnDiskHaldaneHugeBasisShort& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

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


