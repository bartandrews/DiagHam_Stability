////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of bosons on sphere  using the Haldane basis           //
//                            for system size such that                       //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                          work for huge Hilbert space                       //
//                                                                            //
//                        last modification : 08/02/2009                      //
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


#ifndef BOSONONSPHEREHALDANEHUGEBASISSHORT_H
#define BOSONONSPHEREHALDANEHUGEBASISSHORT_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"

#include <iostream>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


class BosonOnSphereHaldaneHugeBasisShort :  public BosonOnSphereShort
{

 protected:

  // the fermionic huge Hilbert space associated to the bosonic one
  FermionOnSphereHaldaneHugeBasis* FermionHugeBasis;

 public:

  // default constructor
  //
  BosonOnSphereHaldaneHugeBasisShort ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a boson
  // maxFileSize = maximum file size (in MBytes)
  // referenceState = array that describes the reference state to start from
  // memory = amount of memory granted for precalculations
  // symmetricFlag = indicate if a symmetric basis has to be used (only available if totalLz = 0)
  // fullDimension = provide the full (i.e. without squeezing) Hilbert space dimension (0 if it has to be computed)
  BosonOnSphereHaldaneHugeBasisShort (int nbrBosons, int totalLz, int lzMax, unsigned long maxFileSize, int* referenceState, unsigned long memory = 10000000, bool symmetricFlag = false, long fullDimension = 0l);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memoryHilbert = amount of memory granted to store the Hilbert space (in Mbytes)
  BosonOnSphereHaldaneHugeBasisShort (char* fileName, unsigned long memoryHilbert);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereHaldaneHugeBasisShort(const BosonOnSphereHaldaneHugeBasisShort& bosons);

  // destructor
  //
  virtual ~BosonOnSphereHaldaneHugeBasisShort ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereHaldaneHugeBasisShort& operator = (const BosonOnSphereHaldaneHugeBasisShort& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // print a given State using the monomial notation
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintStateMonomial (ostream& Str, int state);

  // create the Jack polynomial decomposition corresponding to the root partition
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized hugebasis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // partialSave = save partial results in a given vector file
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized hugebasis
  RealVector& GenerateJackPolynomial(RealVector& jack, double alpha, long minIndex = 0l, long maxIndex = 0l, char* partialSave = 0);

  // create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized hugebasis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // partialSave = save partial results in a given vector file
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized hugebasis
  RealVector& GenerateSymmetrizedJackPolynomial(RealVector& jack, double alpha, long minIndex = 0l, long maxIndex = 0l, char* partialSave = 0);

  // convert a state such that its components are now expressed in the unnormalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has to be normalized to 1
  // return value = converted state
  virtual RealVector& ConvertToUnnormalizedMonomial(RealVector& state, long reference = 0);

  // convert a state such that its components are now expressed in the normalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has been normalized to 1
  // return value = converted state
  virtual RealVector& ConvertFromUnnormalizedMonomial(RealVector& state, long reference = 0);

  // fuse two states which belong to different Hilbert spaces 
  //
  // outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
  // leftVector = reference on the vector whose Hilbert space will be fuse to the left
  // rightVector = reference on the vector whose Hilbert space will be fuse to the right
  // padding = number of unoccupied one body states that have to be inserted between the fused left and right spaces
  // leftSpace = point to the Hilbert space that will be fuse to the left
  // rightSpace = point to the Hilbert space that will be fuse to the right
  // symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
  // return value = reference on the fused state
  RealVector& FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
			  ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace,
			  bool symmetrizedFlag);

  // use product rule to produce part of the components of a system from a smaller one
  //
  // outputVector = reference on the vector which will contain the product rule state  (without zeroing components which do not occur in the fusion)
  // inputVector = reference on the vector associated to the smaller system
  // inputSpace = pointer to the Hilbert space of the smaller system
  // commonPattern = array describing the shared leftmost pattern between the n-body states in both the smaller and larger system sizes
  // commonPatterSize = number of elements in the commonPattern array
  // addedPattern = array describing the pattern that has to be inserted to go from the smaller system to the larger one
  // addedPatterSize = number of elements in the addedPattern array
  // coefficient = multiplicqtive fqctor to go fron the component of the smaller system to the larger one
  // symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
  // return value = reference on the product rule state
  RealVector& ProductRules (RealVector& outputVector, RealVector& inputVector, ParticleOnSphere* inputSpace, 
			    int* commonPattern, int commonPatterSize, int* addedPattern, int addedPatterSize,
			    double coefficient, bool symmetrizedFlag);

  // compute part of the Jack polynomial square normalization in a given range of indices
  //
  // state = reference on the unnormalized Jack polynomial
  // minIndex = first index to compute 
  // nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
  // return value = quare normalization 
  virtual double JackSqrNormalization (RealVector& outputVector, long minIndex = 0l, long nbrComponents = 0l);

};

#endif
