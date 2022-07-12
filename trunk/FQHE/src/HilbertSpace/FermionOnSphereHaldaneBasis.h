////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of fermions on sphere using the Haldane basis            //
//                                                                            //
//                        last modification : 06/07/2006                      //
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


#ifndef FERMIONONSPHEREHALDANEBASIS_H
#define FERMIONONSPHEREHALDANEBASIS_H


#include "config.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "Vector/RationalVector.h"
#include "Vector/LongRationalVector.h"

#include <iostream>


class FermionOnSphere;
class RationalPolynomial;
class LongRationalPolynomial;


class FermionOnSphereHaldaneBasis :  public FermionOnSphere
{

  friend class FermionOnSphereHaldaneSymmetricBasis;
  friend class BosonOnSphereHaldaneBasisShort;
  friend class FermionOnSpherePTruncated;

 protected:

  // topmost state 
  unsigned long ReferenceState;
  // flag to indicate if the reference state is Lz<->-Lz invariant
  bool SymmetricReferenceState;

  // three temporary arrays used during Hilbert space generation
  unsigned long* TmpGeneratedStates;
  int* TmpGeneratedStatesLzMax;
  unsigned long* KeepStateFlag;

 public:

  // default constructor
  //
  FermionOnSphereHaldaneBasis();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = reference on twice the momentum total value, totalLz will be recomputed from referenceState and stored in totalLz
  // lzMax = maximum Lz value reached by a fermion
  // memory = amount of memory granted for precalculations
  // referenceState = array that describes the reference state to start from
  FermionOnSphereHaldaneBasis (int nbrFermions, int& totalLz, int lzMax, int* referenceState,
			       unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnSphereHaldaneBasis (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereHaldaneBasis(const FermionOnSphereHaldaneBasis& fermions);

  // destructor
  //
  virtual ~FermionOnSphereHaldaneBasis ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereHaldaneBasis& operator = (const FermionOnSphereHaldaneBasis& fermions);

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  virtual List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  virtual AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  virtual AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
						 SubspaceSpaceConverter& converter);

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnSphere* targetSpace);

  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  virtual int GetTargetHilbertSpaceDimension();

  // convert a given state from Haldane basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertToNbodyBasis(RealVector& state, FermionOnSphere& nbodyBasis);

  // convert a given state from the usual n-body basis to the Haldane basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertFromNbodyBasis(RealVector& state, FermionOnSphere& nbodyBasis);

  // convert a given state from Haldane basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  LongRationalVector ConvertToNbodyBasis(LongRationalVector& state, FermionOnSphere& nbodyBasis);

  // convert a given state from the usual n-body basis to the Haldane basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  LongRationalVector ConvertFromNbodyBasis(LongRationalVector& state, FermionOnSphere& nbodyBasis);

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

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdA (int index, int m, int n, double& coefficient);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

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
					       LongRational& rhoRootInvAlphaCoef, LongRational& rhoRootConstCoef, unsigned long maxRoot, AbstractArchitecture* architecture, long currentNbrComponents = 0, bool symmetryFlag = false);

  // compute of many coefficients have to be computed to get a single coefficient of the Jack polynomial decomposition 
  //
  // index = index of the component to compute
  // evaluatedCoeffcients = that indicates which coefficients have already been computed
  // tmpMonomial = temporary array for monomial description
  // tmpMonomial2 = temporary array for monomial description
  // maxRoot = fermionic expression for the root partition
  // return value = true if a fully symbolic calculation has been performed
  long GenerateSingleJackPolynomialCoefficientCountOnly(long index, int* evaluatedCoeffcients, unsigned long* tmpMonomial, unsigned long* tmpMonomial2, unsigned long maxRoot);

  // check partitions that may lead to singular coefficient in a given Jack polynomial decomposition, assuming only rational numbers occur
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alphaNumerator = numerator of the Jack polynomial alpha coefficient
  // alphaDenominator = numerator of the Jack polynomial alpha coefficient
  // checkConnectivity = if true, compute how many componets are involved in the calculation of a given singular coefficients
  // return value = vector with non-zero component being rho factor of possible singular coefficients
  RationalVector& CheckPossibleSingularCoefficientsInJackPolynomial(RationalVector& jack, long alphaNumerator, long alphaDenominator, bool checkConnectivity);

  // check partitions that may lead to singular coefficient in a given Jack polynomial decomposition
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // error = error when comparing two rho values
  // return value = vector with non-zero component being rho factor of possible singular coefficients
  RealVector& CheckPossibleSingularCoefficientsInJackPolynomial(RealVector& jack, double alpha, double error);
  
 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);

  // evaluate upper bound for the Haldane basis
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

  // evaluate upper bound for the Haldane basis with shifted values for lzMax and totalLz
  //
  // nbrFermions = number of fermions
  // lzMax = two times momentum maximum value for a fermion plus one 
  // totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
  // return value = Hilbert space dimension  
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int lzMax, unsigned long referenceState, long pos, long& memory);

  // generate all states (i.e. all possible skew symmetric polynomials with fixed Lz)
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // currentLzMax = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, long pos);

};


#endif


