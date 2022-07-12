////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of fermions on a torus                       //
//                                                                            //
//                        last modification : 18/07/2002                      //
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


#ifndef FERMIONONTORUS_H
#define FERMIONONTORUS_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorus.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/SparseComplexMatrix.h"


class FermionOnTorus :  public ParticleOnTorus
{

  friend class FermionOnTorusWithMagneticTranslations;

 protected:

  // number of fermions
  int NbrFermions;
  // number of fermions plus 1
  int IncNbrFermions;
  // maximum momentum value reached by a fermion
  int KyMax;
  // number of Lz values in a state
  int NbrLzValue;

  // index of the momentum orbit
  int TotalKy;
  // index of the momentum orbit
  bool TotalKyFlag;

  // GCD of KyMax and NbrFermions (momemta are defined modulo MomentumModulo)
  int MomentumModulo;
  // value that has to be substracted to the momentum for each translation of the canonical form research
  int MomentumIncrement;
  // shift that has to be done on a state for each translation of the canonical form research
  int StateShift;
  // complementary shift (with respect to MaxMomentum) to StateShift
  int ComplementaryStateShift;
  // mask corresponding to StateShift
  unsigned long MomentumMask;

  // array describing each state
  unsigned long* StateDescription;
  // array giving maximum Lz value reached for a fermion in a given state
  int* StateKyMax;

  // temporary state used when applying ProdA operator
  unsigned long ProdATemporaryState;
  // Lz maximum value associated to temporary state used when applying ProdA operator
  int ProdALzMax;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given maxMomentum sector
  int LookUpTableMemorySize;
  // shift used in each maxMomentum sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used maxMomentum value of the state an the second 
  int** LookUpTable;

  // a table containing ranging from 0 to 2^MaximumSignLookUp - 1
  double* SignLookUpTable;
  // a table containing the mask on the bits to keep for each shift that is requested by sign evaluation
  unsigned long* SignLookUpTableMask;
  // number to evalute size of SignLookUpTable
  int MaximumSignLookUp;

  // pointer to the target space when an index is require after applying basic operation
  FermionOnTorus* TargetSpace;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  FermionOnTorus (int nbrFermions, int maxMomentum);

  // constructor with a constraint of the total momentum of states
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // momentumConstraint = index of the momentum orbit
  FermionOnTorus (int nbrFermions, int maxMomentum, int momentumConstraint);

  // constructor from full datas (with no constraint on the total momentum)
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // hilbertSpaceDimension = Hilbert space dimension
  // stateDescription = array describing each state
  // stateKyMax = array giving maximum Lz value reached for a fermion in a given state
  FermionOnTorus (int nbrFermions, int maxMomentum, int hilbertSpaceDimension, 
		  unsigned long* stateDescription, int* stateKyMax);

  // constructor from full datas
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // momentumConstraint = index of the momentum orbit
  // hilbertSpaceDimension = Hilbert space dimension
  // stateDescription = array describing each state
  // stateKyMax = array giving maximum Lz value reached for a fermion in a given state
  FermionOnTorus (int nbrFermions, int maxMomentum, int momentumConstraint, int hilbertSpaceDimension, 
		  unsigned long* stateDescription, int* stateKyMax);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnTorus(const FermionOnTorus& fermions);

  // destructor
  //
  ~FermionOnTorus ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnTorus& operator = (const FermionOnTorus& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnTorus* targetSpace);
  
  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  virtual int GetTargetHilbertSpaceDimension();

  // return number of particles of the target space
  //
  // return value = number of particles of the target space
  virtual int GetTargetNbrParticles();

  // get the number of orbitals
  //
  // return value = number of orbitals
  virtual int GetNbrOrbitals();

  // get the number of particles
  //
  // return value = number of particles
  virtual int GetNbrParticles();

  // get the particle statistic 
  //
  // return value = particle statistic
  int GetParticleStatistic();

  // get momemtum value of a given state
  //
  // index = state index
  // return value = state momentum
  int GetMomentumValue(int index);

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
					 SubspaceSpaceConverter& converter);

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2[KyMax])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AA (int index, int n1, int n2);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int nbrIndices);

  // apply a^+_m1 a^+_m2 operator to the state produced using AAA method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient);

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
  // m = index for creation operator
  // return value =  resulting multiplicative factor 
  double AdA (int index, int m);

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdA (int index, int m, int n, double& coefficient);

  // return matrix representation of the annihilation operator a_i
  //
  // i = operator index
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  Matrix& A (int i, Matrix& M);

  // return matrix representation ofthw creation operator a^+_i
  //
  // i = operator index
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  Matrix& Ad (int i, Matrix& M);

  // apply an annihilation operator a_i and return the index in the target space
  //
  // i = state index
  // m = index of annihilation operator
  // return value = index in the target space
  virtual int A (int i, int m, double &coefficient); 

  // apply a creation operator a_i and return the index in the target space
  //
  // i = state index
  // m = index of annihilation operator
  // return value = index in the target space
  virtual int Ad (int i, int m, double &coefficient);

  // get a pointer to the target space
  // return value = target space
  //ParticleOnTorus* GetTargetSpace()
  //{ return (ParticleOnTorus*) this->TargetSpace;}

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrFermionSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int kySector, RealVector& groundState);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Ky sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem
  // nbrFermionSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int kySector, ComplexVector& groundState);

  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given Ky sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem
  // nbrFermionSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // return value = entanglement matrix of the subsytem
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrFermionSector, int kySector, RealVector& groundState);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given Ky sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem
  // nbrFermionSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // return value = entanglement matrix of the subsytem
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrFermionSector, int kySector, ComplexVector& groundState);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. 
  // The entanglement matrix is only evaluated in a given Ky sector and fixed number of particles for the part A 
  // but without the Ky constraint for the part B
  // 
  // subsytemSize = number of states that belong to the subsytem
  // nbrFermionSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // return value = entanglement matrix of the subsytem
  virtual RealMatrix EvaluatePartialEntanglementMatrixFullKyPartB (int subsytemSize, int nbrFermionSector, int kySector, RealVector& groundState);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. 
  // The entanglement matrix is only evaluated in a given Ky sector and fixed number of particles for the part A 
  // but without the Ky constraint for the part B
  // 
  // subsytemSize = number of states that belong to the subsytem
  // nbrFermionSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // return value = entanglement matrix of the subsytem
  virtual ComplexMatrix EvaluatePartialEntanglementMatrixFullKyPartB (int subsytemSize, int nbrFermionSector, int kySector, ComplexVector& groundState);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Ky sector.
  // 
  // nbrFermionSector = number of particles that belong to the subsytem 
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrFermionSector, int kySector, RealVector& groundState);

  // create a state from its MPS description
  //
  // bMatrices = array that gives the B matrices 
  // twistMatrix = reference on the twist matrix to insert in the trace
  // state = reference to vector that will contain the state description
  // mPSSumIndices = diagonal indices that have to be kept in the trace
  // nbrMPSSumIndices = number of diagonal indices that have to be kept in the trace
  // memory = amount of memory that can be use to precompute matrix multiplications  
  // initialIndex = initial index to compute
  // nbrComponents = number of components to compute
  virtual void CreateStateFromMPSDescription (SparseRealMatrix* bMatrices, SparseRealMatrix& twistMatrix, RealVector& state, 
					      int* mPSSumIndices, int nbrMPSSumIndices,
					      long memory, long initialIndex, long nbrComponents);
  
  // create a state from its MPS description
  //
  // bMatrices = array that gives the B matrices 
  // twistMatrix = reference on the twist matrix to insert in the trace
  // state = reference to vector that will contain the state description
  // mPSSumIndices = diagonal indices that have to be kept in the trace
  // nbrMPSSumIndices = number of diagonal indices that have to be kept in the trace
  // memory = amount of memory that can be use to precompute matrix multiplications  
  // initialIndex = initial index to compute
  // nbrComponents = number of components to compute
  virtual void CreateStateFromMPSDescription (SparseComplexMatrix* bMatrices, SparseRealMatrix& twistMatrix, ComplexVector& state, 
					      int* mPSSumIndices, int nbrMPSSumIndices,
					      long memory, long initialIndex, long nbrComponents);
  
  // request whether state with given index satisfies a general Pauli exclusion principle
  //
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  // return value = true if teh state satisfies the general Pauli exclusion principle
  bool HasPauliExclusions(int index, int pauliK, int pauliR);

  // apply a magnetic translation along x to a given state
  //
  // index = state index 
  // return value = translated state index
  virtual int ApplyXMagneticTranslation(int index);

  // apply a magnetic translation along x to a given state
  //
  // index = state index 
  // sign = additional sign due to the particle statistics
  // return value = translated state index
  virtual int ApplyXMagneticTranslation(int index, double& sign);

  // transform a state expressed on a torus with a given angle to a state expressed on the same trous but a different angle
  //
  // inputVector = reference on the input vector
  // inputAngle = angle (in radian) between the two vectors that span the torus on which the input state is defined
  // inputAspectRatio = length ratio of the two vectors that span the torus on which the input state is defined
  // outputAngle = angle (in radian) between the two vectors that span the torus on which the output state is defined
  // outputAspectRatio = length ratio of the two vectors that span the torus on which the output state is defined
  // firstComponent = first component of the input vector that has to be symmetrized
  // nbrComponents = number of components of the input vector that have to be symmetrized
  // return value = transformed state
  virtual ComplexVector ChangeTorusAngle (ComplexVector& inputVector, double inputAngle, double inputAspectRatio, double outputAngle, double outputAspectRatio,
					  unsigned long firstComponent, unsigned long nbrComponents);


 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  int FindStateIndex(unsigned long stateDescription, int maxMomentum);

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // return value = Hilbert space dimension
  int EvaluateHilbertSpaceDimension(int nbrFermions, int maxMomentum);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(int memory);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // currentKyMax = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrFermions, int maxMomentum, int currentKyMax, int pos);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion in the state
  // currentKyMax = momentum maximum value for fermions that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentMomentum = current value of the momentum
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrFermions, int maxMomentum, int currentKyMax, int pos, int currentMomentum);

  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in complementary Hilbert space
  // nbrIndex = number of indices to consider in complementary Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnTorus* complementaryHilbertSpace,  ParticleOnTorus* destinationHilbertSpace,
							  RealVector& groundState,  RealSymmetricMatrix* densityMatrix);

  // core part of the C4 rotation
  //
  // inputState = reference on the state that has to be rotated
  // inputSpace = Hilbert space associated to the input state
  // outputState = reference on the rotated state
  // minIndex = minimum index that has to be computed
  // nbrIndices = number of indices that have to be computed
  // clockwise = the rotation is done clockwise
  // return value = reference on the rotated state
  virtual ComplexVector& CoreC4Rotation (ComplexVector& inputState, ParticleOnTorus* inputSpace, ComplexVector& outputState, int minIndex, int nbrIndices, bool clockwise);
  
  // convert a fermionic state to its monomial representation
  //
  // initialState = initial fermionic state in its fermionic representation
  // finalState = reference on the array where the monomial representation has to be stored
  virtual void ConvertToMonomial(unsigned long initialState, unsigned long*& finalState);

  // convert a fermionic state from its monomial representation
  //
  // initialState = array where the monomial representation is stored
  // return value = fermionic state in its fermionic representation
  virtual unsigned long ConvertFromMonomial(unsigned long* initialState);

  // symmetrized a product of two uncoupled states 
  //
  // outputVector = reference on the vector which will contain the symmetrized state
  // leftVector = reference on the vector associated to the first color
  // rightVector = reference on the vector associated to the second color
  // leftSpace = pointer to the Hilbert space of the first color
  // rightSpace = pointer to the Hilbert space of the second color
  // unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
  // return value = symmetrized state
  virtual void SymmetrizeU1U1StateCore (ComplexVector& symmetrizedVector, ComplexVector& leftVector, ComplexVector& rightVector, ParticleOnTorus* leftSpace, ParticleOnTorus* rightSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents);
  
  // symmetrized a product of several uncoupled states 
  //
  // outputState = reference on the output state
  // inputStates = states which will be symmetrized
  // inputSpaces = Hilbert spaces attached to each states
  // nbrStates = number of states to symmetrize
  // firstComponent = first component to symmetrize within the first Hilbert space of inputSpaces
  // nbrComponents = number of components to symmetrize within the first Hilbert space of inputSpaces
  virtual void SymmetrizeU1U1StateCore (ComplexVector& outputState, ComplexVector* inputStates, ParticleOnTorus** inputSpaces, int nbrStates, unsigned long firstComponent, unsigned long nbrComponents);
  
  // symmetrize a vector by grouping distant and equally separated orbitals, core part
  //
  // inputVector = reference on the vector to symmetrize
  // nbrOrbitals = number of orbitals to group together
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
  // first component = index of the first vector component 
  // last component = index of the last component
  virtual void SymmetrizeSingleStateGroupingDistantOrbitalsCore (ComplexVector& inputVector, ComplexVector* symmetrizedVectors, int nbrOrbitals, unsigned long firstComponent, unsigned long nbrComponents, bool twistedTorus = false);

  // symmetrize a vector by keeping only a subset of equally separated orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // firstOrbitalIndex = index of the first orbital to keep
  // symmetrizedVectors = array on the symmetrize states ranging from the smallest Ky to the largest Ky
  // periodicity = momentum periodicity (should be a multiple of the number of orbitals)
  // firstComponent = first component of the input vector that has to be symmetrized
  // nbrComponents = number of components of the input vector that have to be symmetrized
  // return value = symmetrized state
  virtual void SymmetrizeSingleStatePeriodicSubsetOrbitalCore (ComplexVector& inputVector, ComplexVector** symmetrizedVectors, int firstOrbitalIndex, int periodicity, 
							       unsigned long firstComponent, unsigned long nbrComponents);
  
  // symmetrize a vector by keeping only a subset of equally separated orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // firstOrbitalIndex = index of the first orbital to keep
  // symmetrizedVectors = array on the symmetrize states ranging from the smallest Ky to the largest Ky
  // periodicity = momentum periodicity (should be a multiple of the number of orbitals)
  // phase = an optional phase (in pi units) that can be added for each kept and occupied orbital
  // firstComponent = first component of the input vector that has to be symmetrized
  // nbrComponents = number of components of the input vector that have to be symmetrized
  // return value = symmetrized state
  void SymmetrizeSingleStatePeriodicSubsetOrbitalCore (ComplexVector& inputVector, ComplexVector** symmetrizedVectors, int firstOrbitalIndex, int periodicity, double phase, 
						       unsigned long firstComponent, unsigned long nbrComponents);

};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnTorus::GetParticleStatistic()
{
  return ParticleOnTorus::FermionicStatistic;
}

// convert a fermionic state to its monomial representation
//
// initialState = initial fermionic state in its fermionic representation
// finalState = reference on the array where the monomial representation has to be stored

inline void FermionOnTorus::ConvertToMonomial(unsigned long initialState, unsigned long*& finalState)
{
  int Index = 0;
  for (long j = this->KyMax; j >= 0l; --j)
    if (((initialState >> j) & 1ul) != 0ul)
      finalState[Index++] = (unsigned long) j;
}


// convert a fermionic state from its monomial representation
//
// initialState = array where the monomial representation is stored
// return value = fermionic state in its fermionic representation

inline unsigned long FermionOnTorus::ConvertFromMonomial(unsigned long* initialState)
{
  unsigned long TmpState = 0x0ul;  
  for (int j = 0; j < this->NbrFermions; ++j)
    TmpState |= 0x1ul << initialState[j];
  return TmpState;
}

// get the number of orbitals
//
// return value = number of orbitals

inline int FermionOnTorus::GetNbrOrbitals()
{
  return this->KyMax;
}

// get the number of particles
//
// return value = number of particles

inline int FermionOnTorus::GetNbrParticles()
{
  return this->NbrFermions;
}

#endif


