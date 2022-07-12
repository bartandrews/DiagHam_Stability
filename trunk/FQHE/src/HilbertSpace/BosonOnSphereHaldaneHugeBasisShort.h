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


class AbstractArchitecture;


class BosonOnSphereHaldaneHugeBasisShort :  public BosonOnSphereShort
{

 protected:

  // the fermionic huge Hilbert space associated to the bosonic one
  FermionOnSphereHaldaneHugeBasis* FermionHugeBasis;

  // temporary array to store monomial representation
  unsigned long* TemporaryMonomial;
  unsigned long* TemporaryMonomial2;

  friend class FQHESphereJackGeneratorOperation;

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

  // check if disk is used to store the Hilbert space
  //
  // return value = true if disk storage is used
  bool CheckDiskStorage();

  // convert a given state from Haldane basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertToNbodyBasis(RealVector& state, BosonOnSphereHaldaneHugeBasisShort& nbodyBasis);

  // convert a given state from the usual n-body basis to the Haldane basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertFromNbodyBasis(RealVector& state, BosonOnSphereHaldaneHugeBasisShort& nbodyBasis);

  // convert a given state from Haldane basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  LongRationalVector ConvertToNbodyBasis(LongRationalVector& state, BosonOnSphereHaldaneHugeBasisShort& nbodyBasis);

  // convert a given state from the usual n-body basis to the Haldane basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  LongRationalVector ConvertFromNbodyBasis(LongRationalVector& state, BosonOnSphereHaldaneHugeBasisShort& nbodyBasis);

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual long AdAdAA (long index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
  //
  // index = index of the state on which the operator has to be applied
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int nbrIndices);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int nbrIndices, double& coefficient);

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

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdA (int index, int m, int n, double& coefficient);

  // print a given State using the monomial notation
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintStateMonomial (ostream& Str, long state);

  // print a given State using the monomial notation, with one column per particle (using space as a seperator)
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintColumnFormattedStateMonomial (ostream& Str, long state);

  // print a given state using the most compact notation
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintCompactState (ostream& Str, long state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState);

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

  // create the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alphaNumerator = numerator of the Jack polynomial alpha coefficient
  // alphaDenominator = numerator of the Jack polynomial alpha coefficient
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // partialSave = save partial results in a given vector file
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized basis
  LongRationalVector& GenerateJackPolynomial(LongRationalVector& jack, long alphaNumerator, long alphaDenominator, long minIndex, long maxIndex, char* partialSave);

  // create the Jack polynomial decomposition corresponding to the root partition and using sparse storage
  //
  // alpha = value of the Jack polynomial alpha coefficient
  // architecture = architecture to use for precalculation
  // partialSave = save partial results in a given vector file
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // memory = amount of memory (in bytes) allowed for temporary vector storage (0 if the whole vector has to be stored in memory)
  // memoryBlock = amount of memory (in bytes) allowed for precomputing state indices
  // resumeFlag = true if the calculation has to be resumed from a previous one (assuming partialSave contains already computed components)
  virtual void GenerateJackPolynomialSparse(double alpha, AbstractArchitecture* architecture, char* partialSave = 0, long minIndex = 0l, long maxIndex = 0l, long memory = 0l, long memoryBlock = 0l, bool resumeFlag = false);

  // create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry and using sparse storage
  //
  // alpha = value of the Jack polynomial alpha coefficient
  // architecture = architecture to use for precalculation
  // partialSave = save partial results in a given vector file
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // memory = amount of memory (in bytes) allowed for temporary vector storage (0 if the whole vector has to be stored in memory)
  // memoryBlock = amount of memory (in bytes) allowed for precomputing state indices
  // resumeFlag = true if the calculation has to be resumed from a previous one (assuming partialSave contains already computed components)
  virtual void GenerateSymmetrizedJackPolynomialSparse(double alpha, AbstractArchitecture* architecture, char* partialSave, long minIndex = 0l, long maxIndex = 0l, long memory = 0l, long memoryBlock = 0l, bool resumeFlag = false);

  // create the Jack polynomial decomposition corresponding to the root partition and using sparse storage
  //
  // alpha = value of the Jack polynomial alpha coefficient
  // architecture = architecture to use for precalculation
  // partialSave = save partial results in a given vector file
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // memory = amount of memory (in bytes) allowed for temporary vector storage (0 if the whole vector has to be stored in memory)
  // memoryBlock = amount of memory (in bytes) allowed for precomputing state indices
  // resumeFlag = true if the calculation has to be resumed from a previous one (assuming partialSave contains already computed components)
  void GenerateJackPolynomialSparse(long alphaNumerator, long alphaDenominator, AbstractArchitecture* architecture, char* partialSave, long minIndex, long maxIndex, long memory, long memoryBlock, bool resumeFlag);

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

  // fuse two states which belong to different Hilbert spaces 
  //
  // outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
  // leftVector = reference on the vector whose Hilbert space will be fuse to the left
  // rightVector = reference on the vector whose Hilbert space will be fuse to the right
  // padding = number of unoccupied one body states that have to be inserted between the fused left and right spaces
  // leftSpace = point to the Hilbert space that will be fuse to the left
  // rightSpace = point to the Hilbert space that will be fuse to the right
  // symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
  // coefficient = optional multiplicative factor to apply to the fused state 
  // return value = reference on the fused state
  RealVector& FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
			  ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace,
			  bool symmetrizedFlag = false, double coefficient = 1.0);

  // fuse multiple states which belong to different Hilbert spaces 
  //
  // outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
  // nbrInputVectors = number of input vectors
  // inputVectors = input vectors whose Hilbert space will be fuse from  left to right
  // paddings = number of unoccupied one body states that have to be inserted between two consecutive fused spaces
  // inputSpaces = point to the Hilbert space that will be fuse to the left
  // symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
  // coefficient = optional multiplicative factor to apply to the fused state 
  // return value = reference on the fused state
  virtual RealVector& FuseMultipleStates (RealVector& outputVector, int nbrInputVectors, RealVector* inputVectors, int* paddings, 
					  ParticleOnSphere** inputSpaces, bool symmetrizedFlag, double coefficient);

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

  // compute part of the Jack polynomial scalar product in a given range of indices
  //
  // state1 = reference on the first unnormalized Jack polynomial
  // state2 = reference on the second unnormalized Jack polynomial
  // minIndex = first index to compute 
  // nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
  // return value = quare normalization 
  virtual double JackScalarProduct (RealVector& state1, RealVector& state2, long minIndex, long nbrComponents);

  // compute part of the Jack polynomial square normalization in a given range of indices
  //
  // state1 = reference on the first unnormalized Jack polynomial
  // state2 = reference on the second unnormalized Jack polynomial
  // minIndex = first index to compute 
  // nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
  // return value = quare normalization 
  LongRational JackScalarProduct (LongRationalVector& state1, LongRationalVector& state2, long minIndex, long nbrComponents);

 protected :
  
  // core part of the Jack generator using the Lz<->-Lz symmetry and the factorized algorithm
  //
  // invAlpha = inverse of the Jack polynomial alpha coefficient
  // maxRoot = root partition (in fermionic binary representation)
  // partialSave = save partial results in a given vector file
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // indexArray = array where state indices are stored
  // stateArray = array use to store computed state description
  // componentArray = array where computed component numerical factors are stored
  // nbrComputedComponentArray = number of connected components associated to each state through the Jack generator
  // rhoArray = rho factor associated to each state
  virtual void GenerateSymmetrizedJackPolynomialFactorizedCore(double invAlpha, unsigned long maxRoot, long minIndex, long maxIndex, unsigned long** stateArray, double** componentArray, long** indexArray, int* nbrComputedComponents, double* rhoArray);

  // core part of the Jack generator the factorized algorithm
  //
  // invAlpha = inverse of the Jack polynomial alpha coefficient
  // maxRoot = root partition (in fermionic binary representation)
  // partialSave = save partial results in a given vector file
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // indexArray = array where state indices are stored
  // stateArray = array use to store computed state description
  // componentArray = array where computed component numerical factors are stored
  // nbrComputedComponentArray = number of connected components associated to each state through the Jack generator
  // rhoArray = rho factor associated to each state
  virtual void GenerateJackPolynomialFactorizedCore(double invAlpha, unsigned long maxRoot, long minIndex, long maxIndex, unsigned long** stateArray, double** componentArray, long** indexArray, int* nbrComputedComponents, double* rhoArray);

  // core part of the Jack generator using the factorized algorithm
  //
  // invAlpha = inverse of the Jack polynomial alpha coefficient
  // maxRoot = root partition (in fermionic binary representation)
  // partialSave = save partial results in a given vector file
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // indexArray = array where state indices are stored
  // stateArray = array use to store computed state description
  // componentArray = array where computed component numerical factors are stored
  // nbrComputedComponentArray = number of connected components associated to each state through the Jack generator
  // rhoArray = rho factor associated to each state
  void GenerateJackPolynomialFactorizedCore(LongRational invAlpha, unsigned long maxRoot, long minIndex, long maxIndex, unsigned long** stateArray, LongRational** componentArray, long** indexArray, int* nbrComputedComponents, LongRational* rhoArray);

  // core part of multiple state fuse 
  //
  // outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
  // nbrInputVectors = number of input vectors
  // inputVectors = input vectors whose Hilbert space will be fuse from  left to right
  // paddings = number of unoccupied one body states that have to be inserted between two consecutive fused spaces
  // inputSpaces = point to the Hilbert space that will be fuse to the left
  // currentPosition = index of the current space to fuse
  // currentState = current fermionic state obtained by fusing previous states
  // currentCoefficient = current multiplicative coefficient
  // symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
  virtual void CoreFuseMultipleStates (RealVector& outputVector, int nbrInputVectors, RealVector* inputVectors, int* paddings, 
				       BosonOnSphereShort** inputSpaces, int currentPosition, unsigned long currentState, int currentPadding, 
				       double currentCoefficient, bool symmetrizedFlag);

};

// check if disk is used to store the Hilbert space
//
// return value = true if disk storage is used

inline  bool BosonOnSphereHaldaneHugeBasisShort::CheckDiskStorage()
{
  return this->FermionHugeBasis->CheckDiskStorage();
}

// print a given state using the most compact notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

inline ostream& BosonOnSphereHaldaneHugeBasisShort::PrintCompactState (ostream& Str, long state)
{
  return this->FermionHugeBasis->PrintCompactState(Str, state);
}

#endif
