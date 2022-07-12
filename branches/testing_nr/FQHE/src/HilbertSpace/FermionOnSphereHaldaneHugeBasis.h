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

 protected:

  // Hilbert space dimension
  unsigned long HugeHilbertSpaceDimension;

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
  unsigned long PartialHilbertSpaceDimension;
  // shift to get the global index from the local index in a given lowest Lz part of the full Hilbert space 
  unsigned long* StateDescriptionIndexShift;
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
  unsigned long SizeLimit;

  // shift to apply to a state before inverting its expression (useful only for symmetric case calculation)
  int InvertShift;
  // shift to apply to a state after inverting its expression (useful only for symmetric case calculation)
  int InvertUnshift;

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
  FermionOnSphereHaldaneHugeBasis (int nbrFermions, int totalLz, int lzMax, unsigned long maxFileSize, int* referenceState,
				   unsigned long memory = 10000000, bool symmetricFlag = true);

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
  
 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // return value = corresponding index
  unsigned long FindStateIndex(unsigned long stateDescription);

  // evaluate upper bound for the Haldane basis
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual unsigned long EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

  // evaluate upper bound for the Haldane basis with shifted values for lzMax and totalLz
  //
  // nbrFermions = number of fermions
  // lzMax = two times momentum maximum value for a fermion plus one 
  // totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
  // return value = Hilbert space dimension  
  unsigned long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

  // evaluate upper bound for the Haldane basis with shifted values for lzMax and totalLz
  //
  // nbrFermions = number of fermions
  // lzMax = two times momentum maximum value for a fermion plus one 
  // totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
  // return value = Hilbert space dimension  
  unsigned long ShiftedEvaluateHilbertSpaceDimension2(int nbrFermions, int lzMax, int totalLz);

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
  unsigned long  GenerateStates(int lzMax, unsigned long referenceState, unsigned long pos, long& memory);

  // generate all states (i.e. all possible skew symmetric polynomials with fixed Lz)
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // currentLzMax = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  unsigned long  RawGenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, unsigned long pos);

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
  void FindLzMaxSectors(unsigned long* stateDescription, unsigned long dimension, unsigned long* lzSectors, int lzMax);

  // test if a given state corresponds to its canonical expression
  //
  // initialState = state that has to be tested as a canonical state
  // return value = true if it is a canonical state
  bool IsCanonicalState (unsigned long initialState);

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
  unsigned long TmpState = InvertTable[initialState & 0xff] << 56;
  TmpState |= InvertTable[(initialState >> 8) & 0xff] << 48;
  TmpState |= InvertTable[(initialState >> 16) & 0xff] << 40;
  TmpState |= InvertTable[(initialState >> 24) & 0xff] << 32;
  TmpState |= InvertTable[(initialState >> 32) & 0xff] << 24;
  TmpState |= InvertTable[(initialState >> 40) & 0xff] << 16;
  TmpState |= InvertTable[(initialState >> 48) & 0xff] << 8;
  TmpState |= InvertTable[initialState >> 56]; 
#else
  unsigned long TmpState = InvertTable[initialState & 0xff] << 24;
  TmpState |= InvertTable[(initialState >> 8) & 0xff] << 16;
  TmpState |= InvertTable[(initialState >> 16) & 0xff] << 8;
  TmpState |= InvertTable[initialState >> 24];
#endif	
  initialState >>= this->InvertShift;
  TmpState >>= this->InvertUnshift;
  if (TmpState < initialState)
    return false;
  else
    return true;
}

#endif


