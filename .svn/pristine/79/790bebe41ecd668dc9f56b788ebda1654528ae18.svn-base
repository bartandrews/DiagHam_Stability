////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of fermions on sphere using                     //
//                            the Lz <-> -Lz symmetry                         //
//                                                                            //
//                        last modification : 10/03/2007                      //
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


#ifndef FERMIONONSPHERESYMMETRICBASIS_H
#define FERMIONONSPHERESYMMETRICBASIS_H


#include "config.h"
#include "HilbertSpace/FermionOnSphere.h"

#include <iostream>


#ifdef __64_BITS__
#define FERMION_SPHERE_SYMMETRIC_BIT  0x8000000000000000ul
#define FERMION_SPHERE_SYMMETRIC_MASK 0x7ffffffffffffffful
#else
#define FERMION_SPHERE_SYMMETRIC_BIT  0x80000000ul
#define FERMION_SPHERE_SYMMETRIC_MASK 0x7ffffffful
#endif



// precalculation table used to invert a state
//

static  unsigned long InvertTable[] = {0x0ul, 0x80ul, 0x40ul, 0xc0ul, 0x20ul, 0xa0ul, 0x60ul, 0xe0ul, 0x10ul, 0x90ul, 0x50ul, 
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


class FermionOnSphereSymmetricBasis :  public FermionOnSphere
{

  friend class FermionOnSphereHaldaneSymmetricBasis;

  friend class BosonOnSphereSymmetricBasisShort;

 protected:

  // shift to apply to a state before inverting its expression
  int InvertShift;
  // shift to apply to a state after inverting its expression
  int InvertUnshift;


  // signature associated to temporary state used when applying ProdA operator
  unsigned long ProdASignature;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // lzMax = maximum Lz value reached by a fermion
  // memory = amount of memory granted for precalculations
  FermionOnSphereSymmetricBasis (int nbrFermions, int lzMax, unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnSphereSymmetricBasis (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereSymmetricBasis(const FermionOnSphereSymmetricBasis& fermions);

  // destructor
  //
  virtual ~FermionOnSphereSymmetricBasis ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereSymmetricBasis& operator = (const FermionOnSphereSymmetricBasis& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

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

  // get information about any additional symmetry of the Hilbert space
  //
  // return value = symmetry id
  virtual int GetHilbertSpaceAdditionalSymmetry();

  // convert a given state from Lz-symmetric basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertToNbodyBasis(RealVector& state, FermionOnSphere& nbodyBasis);

  // convert a given state from the usual n-body basis to the Lz-symmetric basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphere& nbodyBasis);

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

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AdA (int index, int m, int n, double& coefficient);

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
  
  // compute particule-hole symmetric state from a given state
  //
  // state = vector corresponding to the state to symmetrize
  // holeBasis = n-body basis on which the symmetrized state has to be expressed
  virtual RealVector ParticleHoleSymmetrize (RealVector& state, FermionOnSphere& holeBasis);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);

  // get canonical expression of a given state
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state
  unsigned long GetCanonicalState (unsigned long initialState);

  // get symmetry of a given state 
  //
  // initialState = referennce state whose symmetry has to be computed
  void GetStateSymmetry (unsigned long& initialState);

  // get canonical expression of a given state and its symmetry
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state (with symmetry bit)
  unsigned long GetSignedCanonicalState (unsigned long initialState);

};

// get canonical expression of a given state
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state

inline unsigned long FermionOnSphereSymmetricBasis::GetCanonicalState (unsigned long initialState)
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
    return TmpState;
  else
    return initialState;
}

// get symmetry of a given state 
//
// initialState = reference on the state whose symmetry has to be computed

inline void FermionOnSphereSymmetricBasis::GetStateSymmetry (unsigned long& initialState)
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
  if (TmpState != initialState)    
    initialState |= FERMION_SPHERE_SYMMETRIC_BIT;
}

// get canonical expression of a given state and its symmetry
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state (with symmetry bit)

inline unsigned long FermionOnSphereSymmetricBasis::GetSignedCanonicalState (unsigned long initialState)
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
    {
      return (TmpState | FERMION_SPHERE_SYMMETRIC_BIT);
    }
  else
    if (TmpState != initialState)
      return (initialState | FERMION_SPHERE_SYMMETRIC_BIT);
    else
      return initialState;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

inline int FermionOnSphereSymmetricBasis::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  stateDescription &= FERMION_SPHERE_SYMMETRIC_MASK;
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = (this->StateDescription[PosMid] & FERMION_SPHERE_SYMMETRIC_MASK);
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = (this->StateDescription[PosMid] & FERMION_SPHERE_SYMMETRIC_MASK);
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    return PosMin;
}


#endif


