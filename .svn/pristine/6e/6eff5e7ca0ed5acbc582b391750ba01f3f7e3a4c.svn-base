////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of fermions on sphere with spin with             //
//                           Lz<->-Lz and Sz<->-Sz symmetries                 //
//    that allow LzMax up to 60 (for systems with 128 bit integer support)    //
//               or 27 (on 32 bit systems without 128 bit integer support)    //
//                                                                            //
//                        last modification : 26/09/2008                      //
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


#ifndef FERMIONONSPHEREWITHSPINLZSZSYMMETRYLONG_H
#define FERMIONONSPHEREWITHSPINLZSZSYMMETRYLONG_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"

#include <iostream>


#ifdef __128_BIT_LONGLONG__
#define FERMION_SPHERE_SU2_SYMMETRIC_BIT_LONG       (((ULONGLONG) 0xf000000000000000ul) << 64)
#define FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT_LONG    (((ULONGLONG) 0x8000000000000000ul) << 64)
#define FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT_LONG    (((ULONGLONG) 0x4000000000000000ul) << 64)
#define FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT_LONG (((ULONGLONG) 0xc000000000000000ul) << 64)
#define FERMION_SPHERE_SU2_LZSZ_SYMMETRIC_BIT_LONG  (((ULONGLONG) 0x2000000000000000ul) << 64)
#define FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG (((ULONGLONG) 0xe000000000000000ul) << 64)
#define FERMION_SPHERE_SU2_LZ_SYMMETRIC_TEST_LONG   (((ULONGLONG) 0xa000000000000000ul) << 64)
#define FERMION_SPHERE_SU2_SZ_SYMMETRIC_TEST_LONG   (((ULONGLONG) 0x6000000000000000ul) << 64)
#define FERMION_SPHERE_SU2_SINGLETPARITY_BIT_LONG   (((ULONGLONG) 0x1000000000000000ul) << 64)
#define FERMION_SPHERE_SU2_SZ_TRANSFORMATION_BIT_LONG (((ULONGLONG) 0x0800000000000000ul) << 64)
#define FERMION_SPHERE_SU2_LZ_TRANSFORMATION_BIT_LONG (((ULONGLONG) 0x0400000000000000ul) << 64)
#define FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT_LONG 123
#define FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT_LONG 122
#define FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG ((((ULONGLONG) 0x03fffffffffffffful) << 64) | ((ULONGLONG) 0xfffffffffffffffful))
#define FERMION_SPHERE_SU2_SZ_MASK_LONG ((((ULONGLONG) 0x0155555555555555ul) << 64) | ((ULONGLONG) 0x5555555555555555ul))
#define FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG 124
#else
#define FERMION_SPHERE_SU2_SYMMETRIC_BIT_LONG       (((ULONGLONG) 0xf0000000ul) << 32)
#define FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT_LONG    (((ULONGLONG) 0x80000000ul) << 32)
#define FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT_LONG    (((ULONGLONG) 0x40000000ul) << 32)
#define FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT_LONG (((ULONGLONG) 0xc0000000ul) << 32)
#define FERMION_SPHERE_SU2_LZSZ_SYMMETRIC_BIT_LONG  (((ULONGLONG) 0x20000000ul) << 32)
#define FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG (((ULONGLONG) 0xe0000000ul) << 32)
#define FERMION_SPHERE_SU2_LZ_SYMMETRIC_TEST_LONG   (((ULONGLONG) 0xa0000000ul) << 32)
#define FERMION_SPHERE_SU2_SZ_SYMMETRIC_TEST_LONG   (((ULONGLONG) 0x60000000ul) << 32)
#define FERMION_SPHERE_SU2_SINGLETPARITY_BIT_LONG   (((ULONGLONG) 0x10000000ul) << 32)
#define FERMION_SPHERE_SU2_SZ_TRANSFORMATION_BIT_LONG (((ULONGLONG) 0x08000000ul) << 32)
#define FERMION_SPHERE_SU2_LZ_TRANSFORMATION_BIT_LONG (((ULONGLONG) 0x04000000ul) << 32)
#define FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT_LONG 59
#define FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT_LONG 58
#define FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG ((((ULONGLONG) 0x03fffffful) << 32) | ((ULONGLONG) 0xfffffffful))
#define FERMION_SPHERE_SU2_SZ_MASK_LONG ((((ULONGLONG) 0x01555555ul) << 32) | ((ULONGLONG) 0x55555555ul))
#define FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG 60
#endif

static  ULONGLONG FermionOnSphereWithSpinLzInvertTableLong[] = {0x0ul, 0x40ul, 0x80ul, 0xc0ul, 0x10ul, 0x50ul, 0x90ul, 0xd0ul, 0x20ul, 0x60ul, 0xa0ul, 0xe0ul, 0x30ul, 0x70ul, 0xb0ul, 0xf0ul,
								0x4ul, 0x44ul, 0x84ul, 0xc4ul, 0x14ul, 0x54ul, 0x94ul, 0xd4ul, 0x24ul, 0x64ul, 0xa4ul, 0xe4ul, 0x34ul, 0x74ul, 0xb4ul, 0xf4ul,
								0x8ul, 0x48ul, 0x88ul, 0xc8ul, 0x18ul, 0x58ul, 0x98ul, 0xd8ul, 0x28ul, 0x68ul, 0xa8ul, 0xe8ul, 0x38ul, 0x78ul, 0xb8ul, 0xf8ul,
								0xcul, 0x4cul, 0x8cul, 0xccul, 0x1cul, 0x5cul, 0x9cul, 0xdcul, 0x2cul, 0x6cul, 0xacul, 0xecul, 0x3cul, 0x7cul, 0xbcul, 0xfcul,
								0x1ul, 0x41ul, 0x81ul, 0xc1ul, 0x11ul, 0x51ul, 0x91ul, 0xd1ul, 0x21ul, 0x61ul, 0xa1ul, 0xe1ul, 0x31ul, 0x71ul, 0xb1ul, 0xf1ul,
								0x5ul, 0x45ul, 0x85ul, 0xc5ul, 0x15ul, 0x55ul, 0x95ul, 0xd5ul, 0x25ul, 0x65ul, 0xa5ul, 0xe5ul, 0x35ul, 0x75ul, 0xb5ul, 0xf5ul,
								0x9ul, 0x49ul, 0x89ul, 0xc9ul, 0x19ul, 0x59ul, 0x99ul, 0xd9ul, 0x29ul, 0x69ul, 0xa9ul, 0xe9ul, 0x39ul, 0x79ul, 0xb9ul, 0xf9ul,
								0xdul, 0x4dul, 0x8dul, 0xcdul, 0x1dul, 0x5dul, 0x9dul, 0xddul, 0x2dul, 0x6dul, 0xadul, 0xedul, 0x3dul, 0x7dul, 0xbdul, 0xfdul,
								0x2ul, 0x42ul, 0x82ul, 0xc2ul, 0x12ul, 0x52ul, 0x92ul, 0xd2ul, 0x22ul, 0x62ul, 0xa2ul, 0xe2ul, 0x32ul, 0x72ul, 0xb2ul, 0xf2ul,
								0x6ul, 0x46ul, 0x86ul, 0xc6ul, 0x16ul, 0x56ul, 0x96ul, 0xd6ul, 0x26ul, 0x66ul, 0xa6ul, 0xe6ul, 0x36ul, 0x76ul, 0xb6ul, 0xf6ul,
								0xaul, 0x4aul, 0x8aul, 0xcaul, 0x1aul, 0x5aul, 0x9aul, 0xdaul, 0x2aul, 0x6aul, 0xaaul, 0xeaul, 0x3aul, 0x7aul, 0xbaul, 0xfaul,
								0xeul, 0x4eul, 0x8eul, 0xceul, 0x1eul, 0x5eul, 0x9eul, 0xdeul, 0x2eul, 0x6eul, 0xaeul, 0xeeul, 0x3eul, 0x7eul, 0xbeul, 0xfeul,
								0x3ul, 0x43ul, 0x83ul, 0xc3ul, 0x13ul, 0x53ul, 0x93ul, 0xd3ul, 0x23ul, 0x63ul, 0xa3ul, 0xe3ul, 0x33ul, 0x73ul, 0xb3ul, 0xf3ul,
								0x7ul, 0x47ul, 0x87ul, 0xc7ul, 0x17ul, 0x57ul, 0x97ul, 0xd7ul, 0x27ul, 0x67ul, 0xa7ul, 0xe7ul, 0x37ul, 0x77ul, 0xb7ul, 0xf7ul,
								0xbul, 0x4bul, 0x8bul, 0xcbul, 0x1bul, 0x5bul, 0x9bul, 0xdbul, 0x2bul, 0x6bul, 0xabul, 0xebul, 0x3bul, 0x7bul, 0xbbul, 0xfbul,
								0xful, 0x4ful, 0x8ful, 0xcful, 0x1ful, 0x5ful, 0x9ful, 0xdful, 0x2ful, 0x6ful, 0xaful, 0xeful, 0x3ful, 0x7ful, 0xbful, 0xfful};


class FermionOnSphereWithSpinLzSzSymmetryLong :  public FermionOnSphereWithSpinLong
{

 private:

  // indicate that both the Lz and Sz symmetry have the same parity
  bool LzSzSameParityFlag;

 protected:

  // shift to apply to a state before inverting its expression
  int InvertShift;
  // shift to apply to a state after inverting its expression
  int InvertUnshift;

  // signature associated to temporary state used when applying ProdA operator
  ULONGLONG ProdASignature;

  // additional sign due to the parity sector for the Lz<->-Lz symmetry
  double LzParitySign;

  // additional sign due to the parity sector for the Sz<->-Sz symmetry
  double SzParitySign;

 public:

  // default constructor 
  //
  FermionOnSphereWithSpinLzSzSymmetryLong ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // lzMax = twice the maximum Lz value reached by a fermion
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // minusLzParity = select the  Lz <-> -Lz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinLzSzSymmetryLong (int nbrFermions, int lzMax, bool minusSzParity, bool minusLzParity, unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinLzSzSymmetryLong (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpinLzSzSymmetryLong(const FermionOnSphereWithSpinLzSzSymmetryLong& fermions);

  // destructor
  //
  ~FermionOnSphereWithSpinLzSzSymmetryLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpinLzSzSymmetryLong& operator = (const FermionOnSphereWithSpinLzSzSymmetryLong& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // convert a given state from symmetric basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector  
  virtual RealVector ConvertToNbodyBasis(RealVector& state, FermionOnSphereWithSpinLong& nbodyBasis);

  // convert a given state from the usual n-body basis to the symmetric basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  virtual RealVector ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereWithSpinLong& nbodyBasis);

  // apply a^+_m1_d a^+_m2_d a_n1_d a_n2_d operator to a given state (with m1+m2=n1+n2, only spin down)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // n1 = first index for annihilation operator (spin down)
  // n2 = second index for annihilation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_m1_u a^+_m2_u a_n1_u a_n2_u operator to a given state (with m1+m2=n1+n2, only spin up)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_d_m1 a^+_u_m2 a_d_n1 a_u_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_m_d a_m_d operator to a given state (only spin down)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AddAd (int index, int m);

  // apply a^+_m_u a_m_u operator to a given state  (only spin up)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AduAu (int index, int m);

  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin up)
  // return value =  multiplicative factor 
  virtual double AuAu (int index, int n1, int n2);

  // apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin down)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double AdAd (int index, int n1, int n2);

  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double AuAd (int index, int n1, int n2);

  // apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdu (int m1, int m2, double& coefficient);

  // apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAdd (int m1, int m2, double& coefficient);

  // apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdd (int m1, int m2, double& coefficient);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int* spinIndices, int nbrIndices);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // spinIndices = integer that gives the spin indices associated to each annihilation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int spinIndices, int nbrIndices);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // spinIndices = integer that gives the spin indices associated to each creation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int spinIndices, int nbrIndices, double& coefficient);

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
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(ULONGLONG stateDescription, int lzmax);

  // get canonical expression of a given state
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state
  virtual ULONGLONG GetCanonicalState (ULONGLONG initialState);

  // get symmetry of a given state 
  //
  // initialState = referennce state whose symmetry has to be computed
  virtual void GetStateSymmetry (ULONGLONG& initialState);

  // get canonical expression of a given state and its symmetry
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state (with symmetry bit)
  virtual ULONGLONG GetSignedCanonicalState (ULONGLONG initialState);

  // compute the parity of the number of spin singlet 
  //
  // initialState = reference on the state whose parity has to be evaluated
  void GetStateSingletParity(ULONGLONG& initialState);

  // read Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description is stored
  // return value = true if no error occured
  virtual bool ReadHilbertSpace (char* fileName);

  // factorized code that is used to symmetrize result of the AdxAdy operations
  //
  // state = reference on the state that has been produced with the AdxAdy operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state
  virtual int SymmetrizeAdAdResult(ULONGLONG& state, double& coefficient);

};

// get canonical expression of a given state
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state

inline ULONGLONG FermionOnSphereWithSpinLzSzSymmetryLong::GetCanonicalState (ULONGLONG initialState)
{
  initialState <<= this->InvertShift;
#ifdef __128_BIT_LONGLONG__
  ULONGLONG TmpState = FermionOnSphereWithSpinLzInvertTableLong[initialState & ((ULONGLONG) 0xfful)] << 120;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 8) & ((ULONGLONG) 0xfful)] << 112;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 16) & ((ULONGLONG) 0xfful)] << 104;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 24) & ((ULONGLONG) 0xfful)] << 96;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 32) & ((ULONGLONG) 0xfful)] << 88;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 40) & ((ULONGLONG) 0xfful)] << 80;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 48) & ((ULONGLONG) 0xfful)] << 72;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 56) & ((ULONGLONG) 0xfful)] << 64; 
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 64) & ((ULONGLONG) 0xfful)] << 56;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 72) & ((ULONGLONG) 0xfful)] << 48;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 80) & ((ULONGLONG) 0xfful)] << 40;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 88) & ((ULONGLONG) 0xfful)] << 32;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 96) & ((ULONGLONG) 0xfful)] << 24;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 104) & ((ULONGLONG) 0xfful)] << 16;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 112) & ((ULONGLONG) 0xfful)] << 8;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[initialState >> 120]; 
#else
  ULONGLONG TmpState = FermionOnSphereWithSpinLzInvertTableLong[initialState & ((ULONGLONG) 0xfful)] << 56;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 8) & ((ULONGLONG) 0xfful)] << 48;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 16) & ((ULONGLONG) 0xfful)] << 40;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 24) & ((ULONGLONG) 0xfful)] << 32;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 32) & ((ULONGLONG) 0xfful)] << 24;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 40) & ((ULONGLONG) 0xfful)] << 16;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 48) & ((ULONGLONG) 0xfful)] << 8;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[initialState >> 56];  
#endif
  initialState >>= this->InvertShift;
  TmpState >>= this->InvertUnshift;
  ULONGLONG TmpState2 = ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU2_SZ_MASK_LONG;
  TmpState2 |= TmpState2 << 1;
  TmpState2 ^= initialState;
  ULONGLONG TmpState3 = ((TmpState >> 1) ^ TmpState) & FERMION_SPHERE_SU2_SZ_MASK_LONG;
  TmpState3 |= TmpState3 << 1;
  TmpState3 ^= TmpState;
  if (TmpState2 < initialState)
    initialState = TmpState2;
  if (TmpState3 < TmpState)
    TmpState = TmpState3;
  if (initialState < TmpState)
    return initialState;
  else
    return TmpState;
}

// get symmetry of a given state 
//
// initialState = reference on the state whose symmetry has to be computed

inline void FermionOnSphereWithSpinLzSzSymmetryLong::GetStateSymmetry (ULONGLONG& initialState)
{
  ULONGLONG TmpSymmetryMask = ((ULONGLONG) 0x0ul);
  initialState <<= this->InvertShift;
#ifdef __128_BIT_LONGLONG__
  ULONGLONG TmpState = FermionOnSphereWithSpinLzInvertTableLong[initialState & ((ULONGLONG) 0xfful)] << 120;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 8) & ((ULONGLONG) 0xfful)] << 112;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 16) & ((ULONGLONG) 0xfful)] << 104;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 24) & ((ULONGLONG) 0xfful)] << 96;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 32) & ((ULONGLONG) 0xfful)] << 88;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 40) & ((ULONGLONG) 0xfful)] << 80;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 48) & ((ULONGLONG) 0xfful)] << 72;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 56) & ((ULONGLONG) 0xfful)] << 64; 
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 64) & ((ULONGLONG) 0xfful)] << 56;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 72) & ((ULONGLONG) 0xfful)] << 48;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 80) & ((ULONGLONG) 0xfful)] << 40;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 88) & ((ULONGLONG) 0xfful)] << 32;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 96) & ((ULONGLONG) 0xfful)] << 24;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 104) & ((ULONGLONG) 0xfful)] << 16;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 112) & ((ULONGLONG) 0xfful)] << 8;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[initialState >> 120]; 
#else
  ULONGLONG TmpState = FermionOnSphereWithSpinLzInvertTableLong[initialState & ((ULONGLONG) 0xfful)] << 56;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 8) & ((ULONGLONG) 0xfful)] << 48;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 16) & ((ULONGLONG) 0xfful)] << 40;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 24) & ((ULONGLONG) 0xfful)] << 32;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 32) & ((ULONGLONG) 0xfful)] << 24;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 40) & ((ULONGLONG) 0xfful)] << 16;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 48) & ((ULONGLONG) 0xfful)] << 8;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[initialState >> 56];  
#endif
  initialState >>= this->InvertShift;
  TmpState >>= this->InvertUnshift;
  if (TmpState != initialState)    
    TmpSymmetryMask |= FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT_LONG;
  ULONGLONG TmpMask =  ((TmpState >> 1) ^ TmpState) & FERMION_SPHERE_SU2_SZ_MASK_LONG;
  TmpMask |= TmpMask << 1;
  if ((TmpState ^ TmpMask) != initialState)
    TmpSymmetryMask |= FERMION_SPHERE_SU2_LZSZ_SYMMETRIC_BIT_LONG;
  TmpState = initialState;
  TmpMask =  ((TmpState >> 1) ^ TmpState) & FERMION_SPHERE_SU2_SZ_MASK_LONG;
  TmpMask |= TmpMask << 1;
  TmpState ^= TmpMask; 
  if (TmpState != initialState)    
    TmpSymmetryMask |= FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT_LONG;
  initialState |= TmpSymmetryMask;
}

// get canonical expression of a given state and its symmetry and which transformations have been done to get the canonical expression
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state (with symmetry bit)

inline ULONGLONG FermionOnSphereWithSpinLzSzSymmetryLong::GetSignedCanonicalState (ULONGLONG initialState)
{
  initialState <<= this->InvertShift;
#ifdef __128_BIT_LONGLONG__
  ULONGLONG TmpState = FermionOnSphereWithSpinLzInvertTableLong[initialState & ((ULONGLONG) 0xfful)] << 120;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 8) & ((ULONGLONG) 0xfful)] << 112;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 16) & ((ULONGLONG) 0xfful)] << 104;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 24) & ((ULONGLONG) 0xfful)] << 96;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 32) & ((ULONGLONG) 0xfful)] << 88;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 40) & ((ULONGLONG) 0xfful)] << 80;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 48) & ((ULONGLONG) 0xfful)] << 72;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 56) & ((ULONGLONG) 0xfful)] << 64; 
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 64) & ((ULONGLONG) 0xfful)] << 56;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 72) & ((ULONGLONG) 0xfful)] << 48;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 80) & ((ULONGLONG) 0xfful)] << 40;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 88) & ((ULONGLONG) 0xfful)] << 32;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 96) & ((ULONGLONG) 0xfful)] << 24;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 104) & ((ULONGLONG) 0xfful)] << 16;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 112) & ((ULONGLONG) 0xfful)] << 8;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[initialState >> 120]; 
#else
  ULONGLONG TmpState = FermionOnSphereWithSpinLzInvertTableLong[initialState & ((ULONGLONG) 0xfful)] << 56;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 8) & ((ULONGLONG) 0xfful)] << 48;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 16) & ((ULONGLONG) 0xfful)] << 40;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 24) & ((ULONGLONG) 0xfful)] << 32;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 32) & ((ULONGLONG) 0xfful)] << 24;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 40) & ((ULONGLONG) 0xfful)] << 16;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[(initialState >> 48) & ((ULONGLONG) 0xfful)] << 8;
  TmpState |= FermionOnSphereWithSpinLzInvertTableLong[initialState >> 56];  
#endif
  initialState >>= this->InvertShift;
  TmpState >>= this->InvertUnshift;
  ULONGLONG TmpState2 = ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU2_SZ_MASK_LONG;
  TmpState2 |= TmpState2 << 1;
  TmpState2 ^= initialState;
  ULONGLONG TmpState3 = ((TmpState >> 1) ^ TmpState) & FERMION_SPHERE_SU2_SZ_MASK_LONG;
  TmpState3 |= TmpState3 << 1;
  TmpState3 ^= TmpState;
  ULONGLONG TmpSymmetryMask = ((ULONGLONG) 0x0ul);
  if (initialState != TmpState)
    TmpSymmetryMask |= FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT_LONG;
  if (TmpState2 != initialState)
    TmpSymmetryMask |= FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT_LONG;
  if (TmpState3 != initialState)
    TmpSymmetryMask |= FERMION_SPHERE_SU2_LZSZ_SYMMETRIC_BIT_LONG;
  if (TmpState2 < initialState)
    initialState = TmpState2 | FERMION_SPHERE_SU2_SZ_TRANSFORMATION_BIT_LONG;
  if (TmpState3 < TmpState)
    TmpState = TmpState3 | FERMION_SPHERE_SU2_SZ_TRANSFORMATION_BIT_LONG;
  if ((TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG) <  (initialState & FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG))
    initialState = TmpState | FERMION_SPHERE_SU2_LZ_TRANSFORMATION_BIT_LONG;
  return (initialState | TmpSymmetryMask);
}

// compute the parity of the number of spin singlet 
//
// initialState = reference on the state whose parity has to be evaluated

inline void FermionOnSphereWithSpinLzSzSymmetryLong::GetStateSingletParity(ULONGLONG& initialState)
{
  ULONGLONG TmpState = initialState;
  TmpState &= (TmpState >> 1);
  TmpState &= FERMION_SPHERE_SU2_SZ_MASK_LONG;
#ifdef __128_BIT_LONGLONG__
  TmpState ^= (TmpState >> 64);
#endif
  TmpState ^= (TmpState >> 32);
  TmpState ^= (TmpState >> 16);
  TmpState ^= (TmpState >> 8);
  TmpState ^= (TmpState >> 4);
  TmpState ^= (TmpState >> 2);
  initialState |= (TmpState & ((ULONGLONG) 0x1ul)) << FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG;
}

// factorized code that is used to symmetrize result of the AdxAdy operations
//
// state = reference on the state that has been produced with the AdxAdy operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int FermionOnSphereWithSpinLzSzSymmetryLong::SymmetrizeAdAdResult(ULONGLONG& state, double& coefficient)
{
  ULONGLONG TmpState2 = state;
  state = this->GetSignedCanonicalState(state);
  switch (state & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG)
    {
    case 0x0ul:
      {
	if (this->LzSzSameParityFlag == false)
	  return this->HilbertSpaceDimension;	
	this->GetStateSingletParity(TmpState2);
	if (((1.0 - 2.0 * ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul)))) * this->SzParitySign) < 0.0)
	  return this->HilbertSpaceDimension;
	if ((this->ProdASignature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG) != ((ULONGLONG) 0x0ul))
	  {
	    if ((this->ProdASignature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG) == FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG)
	      coefficient *= 2;
	    else
	      coefficient *= M_SQRT2;
	  }
	int NewLzMax = 1 + (this->LzMax << 1);
	while (((state & FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG) >> NewLzMax) == ((ULONGLONG) 0x0ul))
	  --NewLzMax;
	return this->FindStateIndex(state, NewLzMax);
      }
      break;
    case FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT_LONG :
      {
	if (this->LzSzSameParityFlag == false)
	  return this->HilbertSpaceDimension;	
	this->GetStateSingletParity(TmpState2);
	double TmpSign = (1.0 - 2.0 * ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul))));
	coefficient *= 1.0 + ((double) (((state >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT_LONG) 
					 | (state >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT_LONG)) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->SzParitySign) - 1.0);
	if ((this->ProdASignature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG) == 0x0ul)
	  coefficient *= M_SQRT1_2;
	else
	  if ((this->ProdASignature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG) == FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG)
	    coefficient *= M_SQRT2;
	int NewLzMax = 1 + (this->LzMax << 1);
	while (((state & FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG) >> NewLzMax) == ((ULONGLONG) 0x0ul))
	  --NewLzMax;
	return this->FindStateIndex(state, NewLzMax);
        }
      break;
    case  FERMION_SPHERE_SU2_SZ_SYMMETRIC_TEST_LONG :
      {
	this->GetStateSingletParity(TmpState2);
	double TmpSign = (1.0 - 2.0 * ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul))));
	if ((TmpSign * this->LzParitySign) < 0.0)
	  return this->HilbertSpaceDimension;
	coefficient *= 1.0 + ((double) ((state >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT_LONG) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->SzParitySign) - 1.0);
	if ((this->ProdASignature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG) == ((ULONGLONG) 0x0ul))
	  coefficient *= M_SQRT1_2;
	else
	  if ((this->ProdASignature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG) == FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG)
	    coefficient *= M_SQRT2;
	int NewLzMax = 1 + (this->LzMax << 1);
	while (((state & FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG) >> NewLzMax) == ((ULONGLONG) 0x0ul))
	  --NewLzMax;
	return this->FindStateIndex(state, NewLzMax);
      }
      break;
    case FERMION_SPHERE_SU2_LZ_SYMMETRIC_TEST_LONG :
      {
	this->GetStateSingletParity(TmpState2);
	double TmpSign = (1.0 - 2.0 * ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul))));
	if ((TmpSign * this->SzParitySign) < 0.0)
	  return this->HilbertSpaceDimension;
	coefficient *= 1.0 + ((double) ((state >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT_LONG) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->LzParitySign) - 1.0);
	if ((this->ProdASignature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG) == ((ULONGLONG) 0x0ul))
	  coefficient *= M_SQRT1_2;
	else
	  if ((this->ProdASignature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG) == FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG)
	    coefficient *= M_SQRT2;
	int NewLzMax = 1 + (this->LzMax << 1);
	while (((state & FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG) >> NewLzMax) == ((ULONGLONG) 0x0ul))
	  --NewLzMax;
	return this->FindStateIndex(state, NewLzMax);
      }
      break;      
    default:
      {
	if ((TmpState2 & FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG) != (state & FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG))
	  {
	    this->GetStateSingletParity(TmpState2);
	    double TmpSign = (1.0 - 2.0 * ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & 0x1ul)));
	    coefficient *= 1.0 + ((double) ((state >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT_LONG) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->LzParitySign) - 1.0);
	    coefficient *= 1.0 + ((double) ((state >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT_LONG) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->SzParitySign) - 1.0);
	  }
	if ((this->ProdASignature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG) != FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG)
	  {
	    if ((this->ProdASignature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG) == ((ULONGLONG) 0x0ul))
	      coefficient *= 0.5;
	    else
	      coefficient *= M_SQRT1_2;
	  }
	int NewLzMax = 1 + (this->LzMax << 1);
	while (((state & FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG) >> NewLzMax) == ((ULONGLONG) 0x0ul))
	  --NewLzMax;
	return this->FindStateIndex(state, NewLzMax);      
      }
    }
  return this->HilbertSpaceDimension;
}

#endif


