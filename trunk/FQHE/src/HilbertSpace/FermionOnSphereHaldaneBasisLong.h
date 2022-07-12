////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of fermions on sphere using the Haldane basis            //
//   that allow LzMax up to 127 (for systems with 128 bit integer support)    //
//               or 63 (on 32 bit systems without 128 bit integer support)    //
//                                                                            //
//                        last modification : 22/09/2007                      //
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


#ifndef FERMIONONSPHEREHALDANEBASISLONG_H
#define FERMIONONSPHEREHALDANEBASISLONG_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "Vector/RationalVector.h"
#include "Vector/LongRationalVector.h"

#include <iostream>


class RationalPolynomial;
class LongRationalPolynomial;


class FermionOnSphereHaldaneBasisLong :  public FermionOnSphereLong
{

  friend class FermionOnSphereHaldaneSymmetricBasisLong;

 protected:

  // topmost state 
  ULONGLONG ReferenceState;

  // three temporary arrays used during Hilbert space generation
  ULONGLONG* TmpGeneratedStates;
  int* TmpGeneratedStatesLzMax;
  unsigned long* KeepStateFlag;

 public:

  // default constructor
  //
  FermionOnSphereHaldaneBasisLong();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = reference on twice the momentum total value, totalLz will be recomputed from referenceState and stored in totalLz
  // lzMax = maximum Lz value reached by a fermion
  // memory = amount of memory granted for precalculations
  // referenceState = array that describes the reference state to start from
  FermionOnSphereHaldaneBasisLong (int nbrFermions, int& totalLz, int lzMax, int* referenceState,
				   unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnSphereHaldaneBasisLong (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereHaldaneBasisLong(const FermionOnSphereHaldaneBasisLong& fermions);

  // destructor
  //
  virtual ~FermionOnSphereHaldaneBasisLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereHaldaneBasisLong& operator = (const FermionOnSphereHaldaneBasisLong& fermions);

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // convert a gien state from Haldane basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertToNbodyBasis(RealVector& state, FermionOnSphereLong& nbodyBasis);

  // create the Jack polynomial decomposition corresponding to the root partition
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized basis
  virtual RealVector& GenerateJackPolynomial(RealVector& jack, double alpha);

  // create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized basis
  virtual RealVector& GenerateSymmetrizedJackPolynomial(RealVector& jack, double alpha);

  // check partitions that may lead to singular coefficient in a given Jack polynomial decomposition
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // error = error when comparing two rho values
  // return value = vector with non-zero component being rho factor of possible singular coefficients
  RealVector& CheckPossibleSingularCoefficientsInJackPolynomial(RealVector& jack, double alpha, double error);
  
  // create the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alphaNumerator = numerator of the Jack polynomial alpha coefficient
  // alphaDenominator = numerator of the Jack polynomial alpha coefficient
  // architecture = architecture to use for precalculation
  // symbolicDepth = use symbolic calculation to solve singular values if non zero, if greater than zero it will use symbolic calculation up to a given depth (below that depth it will rely on numerical calculation),
  //                 -1 if the symbolic calculation has to be done up to the point the singular value problem has been solved
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // fileName = optional file name to store temporary calculations
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized basis
  virtual LongRationalVector& GenerateJackPolynomial(LongRationalVector& jack, long alphaNumerator, long alphaDenominator, AbstractArchitecture* architecture, int symbolicDepth = 0, long minIndex = 0l, long maxIndex = 0l, char* fileName = 0);
  
  // create the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur and the resulting state is invariant under the Lz<->-Lz symmetry
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alphaNumerator = numerator of the Jack polynomial alpha coefficient
  // alphaDenominator = numerator of the Jack polynomial alpha coefficient
  // architecture = architecture to use for precalculation
  // symbolicDepth = use symbolic calculation to solve singular values if non zero, if greater than zero it will use symbolic calculation up to a given depth (below that depth it will rely on numerical calculation),
  //                 -1 if the symbolic calculation has to be done up to the point the singular value problem has been solved
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // fileName = optional file name to store temporary calculations
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized basis
  LongRationalVector& GenerateSymmetrizedJackPolynomial(LongRationalVector& jack, long alphaNumerator, long alphaDenominator, AbstractArchitecture* architecture, int symbolicDepth = 0, long minIndex = 0l, long maxIndex = 0l, char* fileName = 0);

  // compute a single coefficient of the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur and using (partial symbolic calculation)
  //
  // index = index of the component to compute
  // numerators = array of polynomials attached to each coefficient numerator
  // denominators = array of polynomials attached to each coefficient denominator
  // tmpMonomial = temporary array for monomial description
  // tmpMonomial2 = temporary array for monomial description
  // rhoRootInvAlphaCoef = coefficient in front of inv alpha in the rho for the root partition
  // rhoRootConstCoef = constant coefficient in the rho for the root partition
  // maxRoot = fermionic expression for the root partition
  // architecture = architecture to use for precalculation
  // currentNbrComponents = current number of components computed for the symbolic of the index-th compoment
  // symmetryFlag = true iof the Lz <-> -Lz symmetry has to be used
  // return value = total number of components computed for the symbolic of the index-th compoment
  long GenerateSingleJackPolynomialCoefficient(long index, LongRationalPolynomial* numerators, LongRationalPolynomial* denominators, long* connectedIndices, long* connectedCoefficients, 
					       unsigned long* tmpMonomial, unsigned long* tmpMonomial2,
					       LongRational& rhoRootInvAlphaCoef, LongRational& rhoRootConstCoef, ULONGLONG maxRoot, AbstractArchitecture* architecture, long currentNbrComponents = 0, bool symmetryFlag = false);

  // compute of many coefficients have to be computed to get a single coefficient of the Jack polynomial decomposition 
  //
  // index = index of the component to compute
  // evaluatedCoeffcients = that indicates which coefficients have already been computed
  // tmpMonomial = temporary array for monomial description
  // tmpMonomial2 = temporary array for monomial description
  // maxRoot = fermionic expression for the root partition
  // return value = true if a fully symbolic calculation has been performed
  long GenerateSingleJackPolynomialCoefficientCountOnly(long index, int* evaluatedCoeffcients, unsigned long* tmpMonomial, unsigned long* tmpMonomial2, ULONGLONG maxRoot);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(ULONGLONG stateDescription, int lzmax);

  // generate all states corresponding to the constraints
  // 
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual int GenerateStates(int lzMax, ULONGLONG referenceState, int pos, long& memory);

  // generate all states (i.e. all possible skew symmetric polynomials with fixed Lz)
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // currentLzMax = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual int RawGenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, int pos);


};

#endif


