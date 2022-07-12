////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of fermions on sphere                       //
//                                                                            //
//                        last modification : 24/06/2002                      //
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


#ifndef FERMIONONSPHEREEFFICIENT_H
#define FERMIONONSPHEREEFFICIENT_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "MathTools/BinomialCoefficients.h"


#include <iostream>

// precalculation table used to invert a state
//

static  unsigned long FermionOnSphereEfficientInvertTable[] = {0x0ul, 0x80ul, 0x40ul, 0xc0ul, 0x20ul, 0xa0ul, 0x60ul, 0xe0ul, 0x10ul, 0x90ul, 0x50ul, 
						      0xd0ul, 0x30ul, 0xb0ul, 0x70ul, 0xf0ul, 0x8ul, 0x88ul, 0x48ul, 0xc8ul, 0x28ul, 0xa8ul, 
						      0x68ul, 0xe8ul, 0x18ul, 0x98ul, 0x58ul, 0xd8ul, 0x38ul, 0xb8ul, 0x78ul, 0xf8ul, 0x4ul, 
						      0x84ul, 0x44ul, 0xc4ul, 0x24ul, 0xa4ul, 0x64ul, 0xe4ul, 0x14ul, 0x94ul, 0x54ul, 0xd4ul, 
						      0x34ul, 0xb4ul, 0x74ul, 0xf4ul, 0xcul, 0x8cul, 0x4cul, 0xccul, 0x2cul, 0xacul, 0x6cul, 
						      0xecul, 0x1cul, 0x9cul, 0x5cul, 0xdcul, 0x3cul, 0xbcul, 0x7cul, 0xfcul, 0x2ul, 0x82ul, 
						      0x42ul, 0xc2ul, 0x22ul, 0xa2ul, 0x62ul, 0xe2ul, 0x12ul, 0x92ul, 0x52ul, 0xd2ul, 0x32ul, 
						      0xb2ul, 0x72ul, 0xf2ul, 0xaul, 0x8aul, 0x4aul, 0xcaul, 0x2aul, 0xaaul, 0x6aul, 0xeaul, 
						      0x1aul, 0x9aul, 0x5aul, 0xdaul, 0x3aul, 0xbaul, 0x7aul, 0xfaul, 0x6ul, 0x86ul, 0x46ul, 
						      0xc6ul, 0x26ul, 0xa6ul, 0x66ul, 0xe6ul, 0x16ul, 0x96ul, 0x56ul, 0xd6ul, 0x36ul, 0xb6ul, 
						      0x76ul, 0xf6ul, 0xeul, 0x8eul, 0x4eul, 0xceul, 0x2eul, 0xaeul, 0x6eul, 0xeeul, 0x1eul, 
						      0x9eul, 0x5eul, 0xdeul, 0x3eul, 0xbeul, 0x7eul, 0xfeul, 0x1ul, 0x81ul, 0x41ul, 0xc1ul, 
						      0x21ul, 0xa1ul, 0x61ul, 0xe1ul, 0x11ul, 0x91ul, 0x51ul, 0xd1ul, 0x31ul, 0xb1ul, 0x71ul, 
						      0xf1ul, 0x9ul, 0x89ul, 0x49ul, 0xc9ul, 0x29ul, 0xa9ul, 0x69ul, 0xe9ul, 0x19ul, 0x99ul, 
						      0x59ul, 0xd9ul, 0x39ul, 0xb9ul, 0x79ul, 0xf9ul, 0x5ul, 0x85ul, 0x45ul, 0xc5ul, 0x25ul, 
						      0xa5ul, 0x65ul, 0xe5ul, 0x15ul, 0x95ul, 0x55ul, 0xd5ul, 0x35ul, 0xb5ul, 0x75ul, 0xf5ul, 
						      0xdul, 0x8dul, 0x4dul, 0xcdul, 0x2dul, 0xadul, 0x6dul, 0xedul, 0x1dul, 0x9dul, 0x5dul, 
						      0xddul, 0x3dul, 0xbdul, 0x7dul, 0xfdul, 0x3ul, 0x83ul, 0x43ul, 0xc3ul, 0x23ul, 0xa3ul, 
						      0x63ul, 0xe3ul, 0x13ul, 0x93ul, 0x53ul, 0xd3ul, 0x33ul, 0xb3ul, 0x73ul, 0xf3ul, 0xbul, 
						      0x8bul, 0x4bul, 0xcbul, 0x2bul, 0xabul, 0x6bul, 0xebul, 0x1bul, 0x9bul, 0x5bul, 0xdbul, 
						      0x3bul, 0xbbul, 0x7bul, 0xfbul, 0x7ul, 0x87ul, 0x47ul, 0xc7ul, 0x27ul, 0xa7ul, 0x67ul, 
						      0xe7ul, 0x17ul, 0x97ul, 0x57ul, 0xd7ul, 0x37ul, 0xb7ul, 0x77ul, 0xf7ul, 0xful, 0x8ful, 
						      0x4ful, 0xcful, 0x2ful, 0xaful, 0x6ful, 0xeful, 0x1ful, 0x9ful, 0x5ful, 0xdful, 0x3ful, 
						      0xbful, 0x7ful, 0xfful};


class FermionOnSphereEfficient :  public ParticleOnSphere
{

  friend class FermionOnSphereUnnormalizedBasis;
  friend class FermionOnSphereHaldaneBasis;
  friend class FermionOnSphereHaldaneLargeBasis;
  friend class FermionOnSphereSymmetricBasis;
  friend class FermionOnSphereHaldaneSymmetricBasis;

  friend class FermionOnSphereWithSpin;
  friend class FermionOnSphereWithSpinAllSz;
  friend class FermionOnSphereWithSpinLong;
  friend class FermionOnSphereWithSU3Spin;
  friend class FermionOnSphereWithSU4Spin;  
  friend class FermionOnSphereWithSpinSqueezedBasis;
  friend class FermionOnSphereTwoLandauLevels;

  friend class BosonOnSphereShort;
  friend class BosonOnSphereHaldaneBasisShort;
  friend class BosonOnSphereHaldaneHugeBasisShort;
  friend class BosonOnSphereBernevigBasisShort;

  friend class BosonOnDiskHaldaneBasisShort;

  friend class BosonOnDiskShort;

  friend class PairedCFOnSphereThreeEightsState;

  friend class ParticleOnSphereProjectorHamiltonian;

 protected:

  // number of fermions
  int NbrFermions;
  // number of fermions plus 1
  int IncNbrFermions;
  // momentum total value
  int TotalLz;
  // maximum Lz value reached by a fermion
  int LzMax;
  // number of Lz values in a state
  int NbrLzValue;

  // array describing each state
  unsigned long* StateDescription;

  // shift for higher bits
  int HighBitLookUpShift;
  // mask for low bits
  unsigned LowBitMask;
  // look-up table for higher bits
  unsigned* HighBitLookUpTable;
  // look-up table with two entries : the first one used lzmax value of the state an the second 
  unsigned* LowBitLookUpTable;

  // a table containing ranging from 0 to 2^MaximumSignLookUp - 1
  double* SignLookUpTable;
  // a table containing the mask on the bits to keep for each shift that is requested by sign evaluation
  unsigned long* SignLookUpTableMask;
  // number to evalute size of SignLookUpTable
  int MaximumSignLookUp;

  // temporary state used when applying ProdA operator
  unsigned long ProdATemporaryState;

  // pointer to the target space when an index is require after applying basic operation
  FermionOnSphereEfficient* TargetSpace;

  // shift to apply to a state before inverting its expression
  int InvertShift;
  // shift to apply to a state after inverting its expression
  int InvertUnshift;

 public:

  // default constuctor
  //
  FermionOnSphereEfficient();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a fermion
  // memory = amount of memory granted for precalculations
  FermionOnSphereEfficient (int nbrFermions, int totalLz, int lzMax, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereEfficient(const FermionOnSphereEfficient& fermions);

  // destructor
  //
  virtual ~FermionOnSphereEfficient ();

  // assignment (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereEfficient& operator = (const FermionOnSphereEfficient& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic();

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

  // apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
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

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual long AdA (long index, int m, int n, double& coefficient);

  // apply creation operator to a word, using the conventions
  // for state-coding and quantum numbers of this space
  // state = word to be acted upon
  // m = Lz value of particle to be added
  // coefficient = reference on the double where the multiplicative factor has to be stored
  virtual unsigned long Ad (unsigned long state, int m, double& coefficient);


  // check whether HilbertSpace implements ordering of operators
  //
  virtual bool HaveOrder ();
  
  // check whether a given operator \prod c^\dagger_m \prod c_n increases or decreases the index of a state
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // return value = 1, if created state is of higher value, 0 if equal, and -1 if lesser value
  virtual int CheckOrder (int* m, int* n, int nbrIndices);

  // get the variance of the state
  // index = index of state to consider
  virtual int StateVariance (int index);

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // print a given State using the monomial notation
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintStateMonomial (ostream& Str, long state);

  // evaluate wave function in real space using a given basis and only for agiven range of components
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
					int firstComponent, int nbrComponent);                                
  
  // evaluate wave functions in real space using a given basis and only for agiven range of components
  //
  // states = array of vector corresponding to the state in the Fock basis
  // nbrStates = number of states in the states array
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // waveFuntions = array where the  wave function values at the given location will be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  virtual void EvaluateWaveFunctions (RealVector* states, int nbrStates, RealVector& position, AbstractFunctionBasis& basis,
				      Complex* waveFuntions, int firstComponent, int nbrComponent);                                
  
  // initialize evaluation of wave function in real space using a given basis and only for a given range of components and
  //
  // timeCoherence = true if time coherence has to be used
  virtual void InitializeWaveFunctionEvaluation (bool timeCoherence = false);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrFermionSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int lzSector, RealVector& groundState);

  // evaluate a density matrix of a shited subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // subsystemShift = shift (in number of states) of the subsytem with repect to the leftmost state (i.e -LzMax)
  // nbrFermionSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluateShiftedPartialDensityMatrix (int subsytemSize, int subsystemShift, int nbrFermionSector, int lzSector, RealVector& groundState);

  // compute part of the Schmidt decomposition, allowing cut in the reduced denisty matrix eigenvalue space
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrFermionSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // eigenvalueCut = discard all contribution from the reduced density matrix whose eigenvalues is lower than eigenvalueCut
  // rebuiltSchmidtGroundState = reference on the state to whose current sector contribution to the Schmidt decomposition will be added 
  // diagonalizedDensityMatrix = reference on the diagonalized reduced density matrix associated to the current sector (with down ordered diagonal elements)
  // transformationMatrix =  reference on the transformation matric that diagonalizes the reduced density matrix
  // return value = reference on rebuiltSchmidtGroundState
  virtual RealVector& EvaluatePartialSchmidtDecomposition(int subsytemSize, int nbrFermionSector, int lzSector, double eigenvalueCut,
							  RealVector& groundState, RealVector& rebuiltSchmidtGroundState,
							  RealDiagonalMatrix& diagonalizedDensityMatrix, RealMatrix& transformationMatrix);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrFermionSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrFermionSector, int lzSector, RealVector& groundState, AbstractArchitecture* architecture = 0);

  // compute particule-hole symmetric state from a given state
  //
  // state = vector corresponding to the state to symmetrize
  // holeBasis = n-body basis on which the symmetrized state has to be expressed
  virtual RealVector ParticleHoleSymmetrize (RealVector& state, FermionOnSphereEfficient& holeBasis);

  // find state index from a string
  //
  // stateDescription = string describing the state
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndex(char* stateDescription);

  // carefully test whether state is in Hilbert-space and find corresponding state index
  //
  // stateDescription = unsigned integer describing the state
  // highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
  // return value = corresponding index, or dimension of space, if not found
  int CarefulFindStateIndex(unsigned long stateDescription, int lzMax);

  // convert a fermionic state to its monomial representation
  //
  // index = index of the fermionic state
  // finalState = reference on the array where the monomial representation has to be stored
  virtual void GetMonomial(long index, unsigned long*& finalState);

  // convert a state such that its components are now expressed in the unnormalized basis
  //
  // state = reference to the state to convert
  // reference = set which component as to be normalized to 1
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
  virtual RealVector& FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
				  ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace, bool symmetrizedFlag = false, double coefficient = 1.0);

  // compute part of the Jack polynomial square normalization in a given range of indices
  //
  // state = reference on the unnormalized Jack polynomial
  // minIndex = first index to compute 
  // nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
  // return value = quare normalization 
  virtual double JackSqrNormalization (RealVector& outputVector, long minIndex = 0l, long nbrComponents = 0l);

  // remove part of each Fock state, discarding component if the Fock state does not a given pattern
  //
  // inputVector = state to truncate
  // reducedSpace = Hilbert space where the truncated state will lie
  // pattern = array describing the pattern 
  // patternSize = pattern size
  // patternShift = indicate where the pattern has to be applied
  // return value = trucated state
  virtual RealVector TruncateStateWithPatternConstraint(RealVector& inputVector, ParticleOnSphere* reducedSpace, int* pattern, int patternSize, int patternShift = 0);

  // get Lz component of a component
  //
  // j = index of the component in Hilbert space
  // return value = twice the Lz component
  virtual int GetLzValue(int j=0);


 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription);

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state (ignored in this class)
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzMax);


  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

  // evaluate Hilbert space dimension with shifted values for lzMax and totalLz
  //
  // nbrFermions = number of fermions
  // lzMax = two times momentum maximum value for a fermion plus one 
  // totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
  // return value = Hilbert space dimension  
  long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

  // generate look-up table for sign calculation
  // 
  virtual void GenerateSignLookUpTable();

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // currentLzMax = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, int pos);

  // get Lz<->-Lz symmetric state of a given state 
  //
  // initialState = reference on the state whose symmetric counterpart has to be computed
  virtual unsigned long GetSymmetricState (unsigned long initialState);

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

};


// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnSphereEfficient::GetParticleStatistic()
{
  return ParticleOnSphere::FermionicStatistic;
}


// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state (ignored in this class)
// return value = corresponding index
inline int FermionOnSphereEfficient::FindStateIndex(unsigned long stateDescription)
{
  int Result = HighBitLookUpTable[stateDescription>>this->HighBitLookUpShift];
  Result += LowBitLookUpTable[stateDescription&this->LowBitMask];
  return Result;
}


// get Lz<->-Lz symmetric state of a given state 
//
// initialState = reference on the state whose symmetric counterpart has to be computed

inline unsigned long FermionOnSphereEfficient::GetSymmetricState (unsigned long initialState)
{
  initialState <<= this->InvertShift;
#ifdef __64_BITS__
  unsigned long TmpState = FermionOnSphereEfficientInvertTable[initialState & 0xff] << 56;
  TmpState |= FermionOnSphereEfficientInvertTable[(initialState >> 8) & 0xff] << 48;
  TmpState |= FermionOnSphereEfficientInvertTable[(initialState >> 16) & 0xff] << 40;
  TmpState |= FermionOnSphereEfficientInvertTable[(initialState >> 24) & 0xff] << 32;
  TmpState |= FermionOnSphereEfficientInvertTable[(initialState >> 32) & 0xff] << 24;
  TmpState |= FermionOnSphereEfficientInvertTable[(initialState >> 40) & 0xff] << 16;
  TmpState |= FermionOnSphereEfficientInvertTable[(initialState >> 48) & 0xff] << 8;
  TmpState |= FermionOnSphereEfficientInvertTable[initialState >> 56]; 
#else
  unsigned long TmpState = FermionOnSphereEfficientInvertTable[initialState & 0xff] << 24;
  TmpState |= FermionOnSphereEfficientInvertTable[(initialState >> 8) & 0xff] << 16;
  TmpState |= FermionOnSphereEfficientInvertTable[(initialState >> 16) & 0xff] << 8;
  TmpState |= FermionOnSphereEfficientInvertTable[initialState >> 24];
#endif	
  TmpState >>= this->InvertUnshift;
  return TmpState;
}

// convert a fermionic state to its monomial representation
//
// index = index of the fermionic state
// finalState = reference on the array where the monomial representation has to be stored

inline void FermionOnSphereEfficient::GetMonomial(long index, unsigned long*& finalState)
{
  this->ConvertToMonomial(this->StateDescription[index], finalState);
}


// convert a fermionic state to its monomial representation
//
// initialState = initial fermionic state in its fermionic representation
// finalState = reference on the array where the monomial representation has to be stored

inline void FermionOnSphereEfficient::ConvertToMonomial(unsigned long initialState, unsigned long*& finalState)
{
  int Index = 0;
  for (long j = this->LzMax; j >= 0l; --j)
    if (((initialState >> j) & 1ul) != 0ul)
      finalState[Index++] = (unsigned long) j;
}


// convert a fermionic state from its monomial representation
//
// initialState = array where the monomial representation is stored
// return value = fermionic state in its fermionic representation

inline unsigned long FermionOnSphereEfficient::ConvertFromMonomial(unsigned long* initialState)
{
  unsigned long TmpState = 0x0ul;  
  for (int j = 0; j < this->NbrFermions; ++j)
    TmpState |= 0x1ul << initialState[j];
  return TmpState;
 }

#endif


