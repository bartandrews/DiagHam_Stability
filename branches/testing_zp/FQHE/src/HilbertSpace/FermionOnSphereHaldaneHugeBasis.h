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


#ifndef FERMIONONSPHEREHALDANEHUGEBASIS_H
#define FERMIONONSPHEREHALDANEHUGEBASIS_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "GeneralTools/List.h"

#include <iostream>


class FermionOnSphereHaldaneHugeBasis :  public ParticleOnSphere
{

  friend class BosonOnSphereHaldaneHugeBasisShort;

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
  // number of states that are used to described the highest Lz part of the full Hilbert space  
  long PartialHilbertSpaceDimension;

  // shift to get the global index from the local index in a given lowest Lz part of the full Hilbert space

   long* StateDescriptionIndexShift;
  // index of the lowest Lz part file for a given highest Lz part
  int* StateDescriptionFileIndex;
  // convertion array to translate from the highest Lz part to its corresponding index in StateDescription
  int* StateHighestLzToIndex;
  // right shift to apply to get the highest Lz part
  int StateHighestLzShift;

  // mask to apply to get the highest Lz part of a state
  unsigned long HighestLzStateMask;

  // temporary buffers used to store the lowest Lz part of the full Hilbert space  
  unsigned long** StateDescriptionBuffers;
  // temporary buffers used to store positon of eqch mqximum Lz sector in the lowest Lz part of the full Hilbert space  
  unsigned long** StateDescriptionLzSectorBuffers;
  // number of temporary buffers
  int NbrBuffers;
  // file index of each buffer
  int* BufferIndices;
  // age of each buffer
  int* BufferAges;
  // convertion array to get buffer id from file id
  int* FileToBuffer;

  // name of the files that contains the lowest Lz part of the full Hilbert space 
  char** StateDescriptionFileNames;
  // size (in sizeof(unsigned long) units) of the files that contains the lowest Lz part of the full Hilbert space  
  unsigned long* StateDescriptionFileSizes;
  // number of files used to store the lowest Lz part of the full Hilbert space  
  int NbrFiles;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given lzmax sector
  unsigned long LookUpTableMemorySize;
  // shift used in each lzmax sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used lzmax value of the state an the second 
  long** LookUpTable;

  // a table containing ranging from 0 to 2^MaximumSignLookUp - 1
  double* SignLookUpTable;
  // a table containing the mask on the bits to keep for each shift that is requested by sign evaluation
  unsigned long* SignLookUpTableMask;
  // number to evalute size of SignLookUpTable
  int MaximumSignLookUp;

  // temporary state used when applying ProdA operator
  unsigned long ProdATemporaryState;
  // Lz maximum value associated to temporary state used when applying ProdA operator
  int ProdALzMax;

  // pointer to the target space when an index is require after applying basic operation
  FermionOnSphereHaldaneHugeBasis* TargetSpace;

  // topmost state 
  unsigned long ReferenceState;

  // three temporary arrays used during Hilbert space generation
  unsigned long* TmpGeneratedStates;
  int* TmpGeneratedStatesLzMax;
  unsigned long* KeepStateFlag;

  // temporary list used during first Hilbert space size evaluation
  List<unsigned long> FileSizes;

  // maximum size for a file
  long SizeLimit;

  // shift to apply to a state before inverting its expression (useful only for symmetric case calculation)
  int InvertShift;
  // shift to apply to a state after inverting its expression (useful only for symmetric case calculation)
  int InvertUnshift;

  // size if the full Hilbert space (i.e. without squeezing)
  long FullLargeHilbertSpaceDimension;

  // glag to indicate disk storage is useless, keep everything in memory
  bool NoDiskFlag;
  // compute the symmetrized basis
  bool SymmetricFlag;

  // number of bytes dedicated to the haeder in the Hilbert space description file
  long FileHeaderSize;
  // name of the file where the Hilbert space is stored
  char* HilbertSpaceFileName;

  // a subset of the complete Hilbert space
  unsigned long* SparseHilbertSpaceDescription;
  // size of the Hilbert space subset
  long SparseHilbertSpaceDimension;
  // size of the Hilbert space chunck associated to one element of the sparse Hilbert space
  long SparseHilbertSpaceChunckSize;
  // size of the Hilbert space chunck associated to the last element of the sparse Hilbert space
  long SparseHilbertSpaceRemainderChunckSize;
  // array which associates a part of the Hilbert space to a buffer
  int* SparseHilbertSpaceBufferIndices;
  // temporary buffers used to store part of the Hilbert space
  unsigned long** SparseBuffers;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a fermion
  // maxFileSize = maximum file size (in MBytes)
  // memory = amount of memory granted for precalculations
  // referenceState = array that describes the reference state to start from
  // symmetricFlag = indicate if a symmetric basis has to be used (only available if totalLz = 0)
  // fullDimension = provide the full (i.e. without squeezing) Hilbert space dimension (0 if it has to be computed)
  FermionOnSphereHaldaneHugeBasis (int nbrFermions, int totalLz, int lzMax, unsigned long maxFileSize, int* referenceState, unsigned long memory = 10000000, bool symmetricFlag = false, long fullDimension = 0l);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memoryHilbert = amount of memory granted to store the Hilbert space (in Mbytes)
  // memory = amount of memory allowed for precalculations
  FermionOnSphereHaldaneHugeBasis(char* fileName, long memoryHilbert, long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereHaldaneHugeBasis(const FermionOnSphereHaldaneHugeBasis& fermions);

  // destructor
  //
  virtual ~FermionOnSphereHaldaneHugeBasis ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereHaldaneHugeBasis& operator = (const FermionOnSphereHaldaneHugeBasis& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  bool WriteHilbertSpace (char* fileName);

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

  // print a given State using the monomial notation
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintStateMonomial (ostream& Str, int state);

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
  
  // initialize evaluation of wave function in real space using a given basis and only for a given range of components and
  //
  // timeCoherence = true if time coherence has to be used
  virtual void InitializeWaveFunctionEvaluation (bool timeCoherence = false);
  
  // create the Jack polynomial decomposition corresponding to the root partition
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // partialSave = save partial results in a given vector file
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized basis
  virtual RealVector& GenerateJackPolynomial(RealVector& jack, double alpha, long minIndex = 0l, long maxIndex = 0l, char* partialSave = 0);

  // create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // partialSave = save partial results in a given vector file
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized basis
  virtual RealVector& GenerateSymmetrizedJackPolynomial(RealVector& jack, double alpha, long minIndex = 0l, long maxIndex = 0l, char* partialSave = 0);

  // create the Jack polynomial decomposition corresponding to the root partition, using an optimized version of the code
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // partialSave = save partial results in a given vector file
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized basis
  virtual RealVector& OptimizedGenerateJackPolynomial(RealVector& jack, double alpha, long minIndex = 0l, long maxIndex = 0l, char* partialSave = 0);

  // create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry, using an optimized version of the code
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // partialSave = save partial results in a given vector file
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized basis
  virtual RealVector& OptimizedGenerateSymmetrizedJackPolynomial(RealVector& jack, double alpha, long minIndex = 0l, long maxIndex = 0l, char* partialSave = 0);

  // compute part of the Jack polynomial square normalization in a given range of indices
  //
  // state = reference on the unnormalized Jack polynomial
  // minIndex = first index to compute 
  // nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
  // return value = quare normalization 
  virtual double JackSqrNormalization (RealVector& outputVector, long minIndex = 0l, long nbrComponents = 0l);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // return value = corresponding index
  long FindStateIndex(unsigned long stateDescription);

  // find state index assuming the whole Hilbert space is stored in memory
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  long FindStateIndexMemory(unsigned long stateDescription, int lzmax);

  // find state index when hilbert space storage is based on sparse algorithm
  //
  // stateDescription = unsigned integer describing the state
  // return value = corresponding index
  long FindStateIndexSparse(unsigned long stateDescription);

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
  long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

  // evaluate upper bound for the Haldane basis with shifted values for lzMax and totalLz
  //
  // nbrFermions = number of fermions
  // lzMax = two times momentum maximum value for a fermion plus one 
  // totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
  // return value = Hilbert space dimension  
  long ShiftedEvaluateHilbertSpaceDimension2(int nbrFermions, int lzMax, int totalLz);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

  // generate look-up table associated to current Hilbert space assuming a huge basis
  // 
  // fileName = name of the binary file
  // memory = memory size that can be allocated for the look-up table
  void GenerateLookUpTableHugeBasis(char* fileName, unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long  GenerateStates(int lzMax, unsigned long referenceState, long pos, long& memory);

  // generate all states (i.e. all possible skew symmetric polynomials with fixed Lz)
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // currentLzMax = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long  RawGenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, long pos);

  // load a part of the Hilbert space associated to one element of the sparse Hilbert basis
  //
  // sparseIndex = index associated to the element of the sparse Hilbert space 
  // buffer = reference on the pointer to the part of the total space 
  // bufferSize = reference on the size of the buffer
  void LoadSparseBuffer(long sparseIndex, unsigned long*& buffer, long& bufferSize);

  // load a lowest Lz part file into memory
  //
  // fileIndex = index of the file to read
  // return file = index of the buffer that has been loaded
  int LoadLowestLzBuffer(int fileIndex);

  // find position of each sector with a given maximum Lz
  // 
  // stateDescription = array that contains state description
  // dimension = dimension of the stateDescription array
  // lzSectors = array that where the position of each sector with a given maximum Lz will be stored
  // lzMax = maximum momentum value for a fermion
  void FindLzMaxSectors(unsigned long* stateDescription, long dimension, unsigned long* lzSectors, int lzMax);

  // test if a given state corresponds to its canonical expression
  //
  // initialState = state that has to be tested as a canonical state
  // return value = true if it is a canonical state
  bool IsCanonicalState (unsigned long initialState);

  // get Lz<->-Lz symmetric state of a given state 
  //
  // initialState = reference on the state whose symmetric counterpart has to be computed
  virtual unsigned long GetSymmetricState (unsigned long initialState);

  // convert a bosonic state to its monomial representation
  //
  // initialState = initial bosonic state in its fermionic representation
  // initialStateLzMax = initial bosonic state maximum Lz value
  // finalState = reference on the array where the monomial representation has to be stored
  virtual void ConvertToMonomial(unsigned long initialState, int*& finalState);

  // convert a bosonic state from its monomial representation
  //
  // initialState = array where the monomial representation is stored
  // return value = bosonic state in its fermionic representation
  virtual unsigned long ConvertFromMonomial(int* initialState);

  // find squeezed partitions that are connected throught the Jack calculation algorithm
  //
  // nbrPartitions = number of partitions whose connection has to be computed
  // partitionIndices = indices of the partitions whose connection has to be computed
  // nbrConnectedPartitions = array where the number of partitions connected to a given one will be stored
  // connectedPartitionIndices = array where the index of the partitions connected to a given one will be stored
  // factors = numerical factor that relates two connected partitions
  // rootPartition = Jack root partition
  virtual void GetConnectedSqueezedPartitions(long nbrPartitions, long* partitionIndices, int* nbrConnectedPartitions, 
					      long** connectedPartitionIndices, double** factors, unsigned long rootPartition);

};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnSphereHaldaneHugeBasis::GetParticleStatistic()
{
  return ParticleOnSphere::FermionicStatistic;
}

// test if a given state corresponds to its canonical expression
//
// initialState = state that has to be tested as a canonical state
// return value = true if it is a canonical state

inline bool FermionOnSphereHaldaneHugeBasis::IsCanonicalState (unsigned long initialState)
{
  initialState <<= this->InvertShift;
#ifdef __64_BITS__
  unsigned long TmpState = FermionOnSphereInvertTable[initialState & 0xff] << 56;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 8) & 0xff] << 48;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 16) & 0xff] << 40;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 24) & 0xff] << 32;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 32) & 0xff] << 24;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 40) & 0xff] << 16;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 48) & 0xff] << 8;
  TmpState |= FermionOnSphereInvertTable[initialState >> 56]; 
#else
  unsigned long TmpState = FermionOnSphereInvertTable[initialState & 0xff] << 24;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 8) & 0xff] << 16;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 16) & 0xff] << 8;
  TmpState |= FermionOnSphereInvertTable[initialState >> 24];
#endif	
  initialState >>= this->InvertShift;
  TmpState >>= this->InvertUnshift;
  if (TmpState < initialState)
    return false;
  else
    return true;
}

// get Lz<->-Lz symmetric state of a given state 
//
// initialState = reference on the state whose symmetric counterpart has to be computed

inline unsigned long FermionOnSphereHaldaneHugeBasis::GetSymmetricState (unsigned long initialState)
{
  initialState <<= this->InvertShift;
#ifdef __64_BITS__
  unsigned long TmpState = FermionOnSphereInvertTable[initialState & 0xff] << 56;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 8) & 0xff] << 48;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 16) & 0xff] << 40;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 24) & 0xff] << 32;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 32) & 0xff] << 24;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 40) & 0xff] << 16;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 48) & 0xff] << 8;
  TmpState |= FermionOnSphereInvertTable[initialState >> 56]; 
#else
  unsigned long TmpState = FermionOnSphereInvertTable[initialState & 0xff] << 24;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 8) & 0xff] << 16;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 16) & 0xff] << 8;
  TmpState |= FermionOnSphereInvertTable[initialState >> 24];
#endif	
  TmpState >>= this->InvertUnshift;
  return TmpState;
}

// convert a bosonic state to its monomial representation
//
// initialState = initial bosonic state in its fermionic representation
// initialStateLzMax = initial bosonic state maximum Lz value
// finalState = reference on the array where the monomial representation has to be stored

inline void FermionOnSphereHaldaneHugeBasis::ConvertToMonomial(unsigned long initialState, int*& finalState)
{
  int Index = 0;
  for (int j = this->LzMax; j >= 0; --j)
    if (((initialState >> j) & 1ul) != 0ul)
      finalState[Index++] = j;
}

// convert a bosonic state from its monomial representation
//
// initialState = array where the monomial representation is stored
// return value = bosonic state in its fermionic representation

inline unsigned long FermionOnSphereHaldaneHugeBasis::ConvertFromMonomial(int* initialState)
{
  unsigned long TmpState = 0x0ul;  
  for (int j = 0; j < this->NbrFermions; ++j)
    TmpState |= 0x1ul << initialState[j];
  return TmpState;
 }

#endif


