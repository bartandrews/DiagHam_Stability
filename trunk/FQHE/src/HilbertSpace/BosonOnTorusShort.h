////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of bosons on torusfor system size such that            //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 08/12/2009                      //
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


#ifndef BOSONONTORUSSHORT_H
#define BOSONONTORUSSHORT_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorus.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>
#include <limits>


class BosonOnTorusShort :  public ParticleOnTorus
{

  friend class BosonOnSquareLatticeWannierSpace;
  friend class BosonOnTorusWithSpin;
  friend class BosonOnTorusWithMagneticTranslationsShort;

 protected:

  // number of bosons
  int NbrBosons;
  // number of bosons plus 1
  int IncNbrBosons;
  // maximum bosonic occupation per orbital
  int MaxOccupation;

  // maximum momentum value reached by a boson
  int KyMax;
  // number of Ky values in a state
  int NbrKyValue;

  // index of the momentum orbit
  int TotalKy;
  // index of the momentum orbit
  bool TotalKyFlag;
  //  GCD of nbrBosons and maxMomentum
  int MomentumModulo;
  //  modulo to apply on th Ky momentum
  int KyMomentumModulo;

  // shift that has to be done on a state for each translation of the canonical form research
  int StateShift;
  // mask that corresponds to last bit that can be set to one
  unsigned long LastMomentumMask;

  // array describing each state
  unsigned long* StateDescription;
  // array giving maximum Lz value reached for a boson in a given state
  int* StateKyMax;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given lzmax sector
  unsigned long LookUpTableMemorySize;
  // shift used in each lzmax sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used lzmax value of the state an the second 
  int** LookUpTable;
  // look-up table for Hilbert spaces larger than 2^31
  long** LargeLookUpTable;

  // temporary state used when applying operators
  unsigned long* TemporaryState;
  int TemporaryStateKyMax;
  // temporary state used when applying ProdA operator
  unsigned long* ProdATemporaryState;
   int ProdATemporaryStateKyMax;

  // pointer to the target space when an index is require after applying basic operation
  BosonOnTorusShort* TargetSpace;

 
 public:

  // default constructor
  // 
  BosonOnTorusShort ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  BosonOnTorusShort (int nbrBosons, int maxMomentum);

  // constructor with a constraint of the total momentum of states
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  // momentumConstraint = index of the momentum orbit
  BosonOnTorusShort (int nbrBosons, int maxMomentum, int momentumConstraint);
  
  // constructor with a constraint of the total momentum of states and the maximum occupation per orbital
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  // momentumConstraint = index of the momentum orbit
  // maxOccupation = maximum bosonic occupation per orbital
  BosonOnTorusShort (int nbrBosons, int maxMomentum, int momentumConstraint, int maxOccupation);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnTorusShort(const BosonOnTorusShort& bosons);

  // destructor
  //1
  ~BosonOnTorusShort ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnTorusShort& operator = (const BosonOnTorusShort& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  void SetTargetSpace(ParticleOnTorus* targetSpace);

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

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient);

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
  // coefficient = will be multiplied by the prefactor of the bosonic ladder operator
  // return value = index in the target space
  virtual int A (int i, int m, double &coefficient);

  // apply a creation operator a_i and return the index in the target space
  //
  // i = state index
  // m = index of annihilation operator
  // coefficient = will be multiplied by the prefactor of the bosonic ladder operator
  // return value = index in the target space
  virtual int Ad (int i, int m, double &coefficient);

  // get a pointer to the target space
  // return value = target space
  //ParticleOnTorus* GetTargetSpace()
  //{ return (ParticleOnTorus*) this->TargetSpace;}

  // print a given State using the monomial notation
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintStateMonomial (ostream& Str, long state);
  
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
  
  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state   
  virtual int AdA (int index, int m, int n, double& coefficient);
  
  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m  
  virtual double AdA (int index, int m);
  
  // apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor   
  virtual double AA (int index, int n1, int n2);
  
  // apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Ky sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int , RealVector& groundState);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Ky sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int , ComplexVector& groundState);

  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given Ky sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // return value = entanglement matrix of the subsytem
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrBosonSector, int kySector, RealVector& groundState);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given Ky sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // return value = entanglement matrix of the subsytem
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrBosonSector, int kySector, ComplexVector& groundState);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. 
  // The entanglement matrix is only evaluated in a given Ky sector and fixed number of particles for the part A 
  // but without the Ky constraint for the part B
  // 
  // subsytemSize = number of states that belong to the subsytem
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // return value = entanglement matrix of the subsytem
  virtual RealMatrix EvaluatePartialEntanglementMatrixFullKyPartB (int subsytemSize, int nbrBosonSector, int kySector, RealVector& groundState);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. 
  // The entanglement matrix is only evaluated in a given Ky sector and fixed number of particles for the part A 
  // but without the Ky constraint for the part B
  // 
  // subsytemSize = number of states that belong to the subsytem
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // return value = entanglement matrix of the subsytem
  virtual ComplexMatrix EvaluatePartialEntanglementMatrixFullKyPartB (int subsytemSize, int nbrBosonSector, int kySector, ComplexVector& groundState);

  // apply a magnetic translation along x to a given state
  //
  // index = state index 
  // return value = translated state index
  virtual int ApplyXMagneticTranslation(int index);

  // remove part of each Fock state, discarding component if the Fock state does not a given pattern
  //
  // inputVector = state to truncate
  // reducedSpace = Hilbert space where the truncated state will lie
  // pattern = array describing the pattern 
  // patternSize = pattern size
  // patternShift = indicate where the pattern has to be applied
  // return value = trucated state
  virtual RealVector TruncateStateWithPatternConstraint(RealVector& inputVector, ParticleOnTorus* reducedSpace, int* pattern, int patternSize, int patternShift = 0);
  
  // symmetrize a vector with even number of orbitals 
  //
  // outputVector = reference on the vector which will contain the symmetrozed state
  // leftVector = reference on the vector to be symmetrized
  // leftSpace = pointer to the Hilbert space
  // unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
  // return value = symmetrized state
  virtual RealVector SymmetrizeU1U1SingleState (RealVector& leftVector, ParticleOnTorus* leftSpace, bool oneInTwoFlag, bool unnormalizedBasisFlag = false, AbstractArchitecture* architecture = 0);
  
  // symmetrize a vector with even number of orbitals 
  //
  // outputVector = reference on the vector which will contain the symmetrozed state
  // leftVector = reference on the vector to be symmetrized
  // leftSpace = pointer to the Hilbert space
  // groupNeighbouringOrbitals = group neighbouring orbitals instead of grouping orbitals separated by Nphi
  // return value = symmetrized state
  virtual ComplexVector SymmetrizeU1U1SingleState (ComplexVector& leftVector, ParticleOnTorus* leftSpace, bool groupNeighbouringOrbitals,  AbstractArchitecture* architecture = 0);
  
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

  // request whether state with given index satisfies a general Pauli exclusion principle
  //
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  // return value = true if teh state satisfies the general Pauli exclusion principle
  virtual bool HasPauliExclusions(int index, int pauliK, int pauliR);
  
 protected:

  // convert a bosonic state into its fermionic counterpart
  //
  // initialState = reference on the array where initialbosonic  state is stored
  // initialStateLzMax = reference on the initial bosonic state maximum Lz value
  // return value = corresponding fermionic state
  unsigned long BosonToFermion(unsigned long*& initialState, int& initialStateLzMax);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialState = initial fermionic state
  // initialStateLzMax = initial fermionic state maximum Lz value
  // finalState = reference on the array where the bosonic state has to be stored
  // finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored
  void FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState, int& finalStateLzMax);

  // find state index
  //
  // stateDescription = array describing the state
  // lzmax = maximum Lz value reached by a boson in the state
  // return value = corresponding index
  int FindStateIndex(unsigned long stateDescription, int lzmax);

  // find state index when the Hilbert space is larger than 2^31
  //
  // stateDescription = array describing the state
  // lzmax = maximum Lz value reached by a boson in the state
  // return value = corresponding index
  long FindStateLargeIndex(unsigned long stateDescription, int lzmax);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrBosons, int maxMomentum);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson in the state
  // currentKyMax = momentum maximum value for bosons that are still to be placed
  // currentMomentum = current value of the momentum
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrBosons, int maxMomentum, int currentKyMax, int currentMomentum);

  // evaluate Hilbert space dimension, including a truncation on the occupation numbers
  //
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson in the state
  // currentKyMax = momentum maximum value for bosons that are still to be placed
  // currentMomentum = current value of the momentum
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimensionTruncatedOccupation(int nbrBosons, int maxMomentum, int currentKyMax, int currentMomentum);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(int memory);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson in the state
  // currentKyMax = momentum maximum value for bosons that are still to be placed
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int nbrBosons, int maxMomentum, int currentKyMax,long  pos);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson in the state
  // currentKyMax = momentum maximum value for bosons that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentMomentum = current value of the momentum
  // return value = position from which new states have to be stored
  long GenerateStates(int nbrBosons, int maxMomentum, int currentKyMax, long pos, int currentMomentum);

  // generate all states corresponding to the constraints, including a truncation on the occupation numbers
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson in the state
  // currentKyMax = momentum maximum value for bosons that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentMomentum = current value of the momentum
  // return value = position from which new states have to be stored
  long GenerateStatesTruncatedOccupation(int nbrBosons, int maxMomentum, int currentKyMax, long pos, int currentMomentum);

  // convert a bosonic state to its monomial representation
  //
  // initialState = initial  bosonic state
  // initialStateLzMax = initial bosonic state maximum Ky value
  // finalState = reference on the array where the monomial representation has to be stored
  void ConvertToMonomial(unsigned long* initialState, int initialStateKyMax, unsigned long*& finalState);

  // convert a bosonic state to its monomial representation
  //
  // initialState = initial bosonic state in its fermionic representation
  // initialStateLzMax = initial bosonic state maximum Lz value
  // finalState = reference on the array where the monomial representation has to be stored
  void ConvertToMonomial(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState);

  // convert a bosonic state from its monomial representation
  //
  // initialState = array where the monomial representation is stored
  // finalState = bosonic state
  // return value = maximum Ky value reached by a particle
  int ConvertFromMonomial(unsigned long* initialState, unsigned long*& finalState);

  // convert a bosonic state from its monomial representation
  //
  // initialState = array where the monomial representation is stored
  // return value = bosonic state in its fermionic representation
  unsigned long ConvertFromMonomial(unsigned long* initialState);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Ky sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int kySector, RealVector& groundState);
  HermitianMatrix  EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int kySector, ComplexVector& groundState);
  
  // symmetrized a product of two uncoupled states 
  //
  // outputVector = reference on the vector which will contain the symmetrized state
  // leftVector = reference on the vector associated to the first color
  // rightVector = reference on the vector associated to the second color
  // leftSpace = pointer to the Hilbert space of the first color
  // rightSpace = pointer to the Hilbert space of the second color
  // unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
  // return value = symmetrized state
  virtual void SymmetrizeU1U1StateCore (RealVector& symmetrizedVector, RealVector& leftVector, RealVector& rightVector, ParticleOnTorus* leftSpace, ParticleOnTorus* rightSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents);
  
  // symmetrized a product of two uncoupled states 
  //
  // outputVector = reference on the vector which will contain the symmetrozed state
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
  
  // symmetrize a vector with even number of orbitals
  //
  // outputVector = reference on the vector which will contain the symmetrozed state
  // leftVector = reference on the vector associated to the first color
  // leftSpace = pointer to the Hilbert space of the first color
  // return value = symmetrized state
  virtual void SymmetrizeU1U1SingleStateCore (RealVector& symmetrizedVector, RealVector& leftVector, ParticleOnTorus* leftSpace, unsigned long firstComponent, unsigned long nbrComponents);
  
  // symmetrize a vector by grouping neighbouring orbitals, core part
  //
  // inputVector = reference on the vector to symmetrize
  // nbrOrbitals = number of orbitals to group together
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
  // first component = index of the first vector component 
  // last component = index of the last component
  virtual void SymmetrizeSingleStateGroupingNeighbouringOrbitalsCore (ComplexVector& inputVector, ComplexVector* symmetrizedVectors, int nbrOrbitals, 
								unsigned long firstComponent, unsigned long nbrComponents);

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
  

};

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnTorusShort::GetParticleStatistic()
{
  return ParticleOnTorus::BosonicStatistic;
}

// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initialbosonic  state is stored
// initialStateLzMax = reference on the initial bosonic state maximum Lz value
// return value = corresponding fermionic state

inline unsigned long BosonOnTorusShort::BosonToFermion(unsigned long*& initialState, int& initialStateLzMax)
{
  unsigned long TmpState = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= initialStateLzMax; ++i)
    {
      TmpState |= ((1ul << initialState[i]) - 1ul) << Shift;
      Shift += initialState[i];
      ++Shift;
    }
  return TmpState;
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// initialStateLzMax = initial fermionic state maximum Lz value
// finalState = reference on the array where the bosonic state has to be stored
// finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored

inline void BosonOnTorusShort::FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState, int& finalStateLzMax)
{
  finalStateLzMax = 0;
  while (initialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialState - 1ul) ^ (~initialState);
      TmpState &= ~(TmpState >> 1);
#ifdef  __64_BITS__
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0f0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ffff0000ul) != 0) << 4;      
      TmpPower |= ((TmpState & 0xffffffff00000000ul) != 0) << 5;      
#else
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ul) != 0) << 4;      
#endif
      finalState[finalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState >>= TmpPower;
      ++finalStateLzMax;
      initialStateLzMax -= TmpPower;
    }
  --finalStateLzMax;
}

// convert a bosonic state to its monomial representation
//
// initialState = initial  bosonic state
// initialStateKyMax = initial bosonic state maximum Ky value
// finalState = reference on the array where the monomial representation has to be stored

inline void BosonOnTorusShort::ConvertToMonomial(unsigned long* initialState, int initialStateKyMax, unsigned long*& finalState)
{
  int Index = 0;
  for (int i = initialStateKyMax; i >= 0; --i)
    for (unsigned long j = 0ul; j < initialState[i]; ++j)
      finalState[Index++] = i;
}

// convert a bosonic state to its monomial representation
//
// initialState = initial bosonic state in its fermionic representation
// initialStateKyMax = initial bosonic state maximum Ky value
// finalState = reference on the array where the monomial representation has to be stored

inline void BosonOnTorusShort::ConvertToMonomial(unsigned long initialState, int initialStateKyMax, unsigned long*& finalState)
{
  int Index = 0;
  int TmpKy = initialStateKyMax - this->NbrBosons + 1;
  while (initialStateKyMax >= 0)
    {
      while ((initialStateKyMax >= 0) && (((initialState >> initialStateKyMax) & 0x1ul) != 0x0ul))
	{
	  finalState[Index++] = TmpKy;
	  --initialStateKyMax;
	}
      while ((initialStateKyMax >= 0) && (((initialState >> initialStateKyMax) & 0x1ul) == 0x0ul))
	{
	  --TmpKy;
	  --initialStateKyMax;
	}
    }
}

// convert a bosonic state from its monomial representation
//
// initialState = array where the monomial representation is stored
// finalState = bosonic state
// return value = maximum Ky value reached by a particle

inline int BosonOnTorusShort::ConvertFromMonomial(unsigned long* initialState, unsigned long*& finalState)
{
  int TmpKyMax = initialState[0]; 
  for (int i = 0; i < TmpKyMax; ++i)
    finalState[i] = 0ul;
  for (int i = 0; i < this->NbrBosons; ++i)
    finalState[initialState[i]]++;
  return TmpKyMax;
}

// convert a bosonic state from its monomial representation
//
// initialState = array where the monomial representation is stored
// return value = bosonic state in its fermionic representation

inline unsigned long BosonOnTorusShort::ConvertFromMonomial(unsigned long* initialState)
{
  unsigned long Tmp = 0x0ul;
  for (int i = 0; i < this->NbrBosons; ++i)
    Tmp |= 0x1ul << (initialState[i] + ((unsigned long) (this->NbrBosons - i)) - 1ul);
  return Tmp;
}

// get the number of orbitals
//
// return value = number of orbitals

inline int BosonOnTorusShort::GetNbrOrbitals()
{
  return this->KyMax;
}

// get the number of particles
//
// return value = number of particles

inline int BosonOnTorusShort::GetNbrParticles()
{
  return this->NbrBosons;
}

#endif


