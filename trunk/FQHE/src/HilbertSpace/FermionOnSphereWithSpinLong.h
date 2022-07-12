////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of fermions on sphere with spin that allow LzMax up to         //
//                  63 (for systems with 128 bit integer support)             //
//               or 31 (on 32 bit systems without 128 bit integer support)    //
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


#ifndef FERMIONONSPHEREWITHSPINLONG_H
#define FERMIONONSPHEREWITHSPINLONG_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"

#include <iostream>


class FermionOnSphere;


class FermionOnSphereWithSpinLong :  public ParticleOnSphereWithSpin
{

  friend class FermionOnSphereWithSpinLzSzSymmetryLong;
  friend class FermionOnSphereWithSpinSzSymmetryLong;
  friend class FermionOnSphereWithSpinLzSymmetryLong;
  friend class FermionOnSphereWithSpinHaldaneBasisLong;

 protected:

  // number of fermions
  int NbrFermions;
  // number of fermions plus 1
  int IncNbrFermions;
  // momentum total value
  int TotalLz;
  // maximum Lz value reached by a fermion
  int LzMax;
  // number of fermions with spin up / down
  int NbrFermionsUp;
  int NbrFermionsDown;
  // number of Lz values in a stat
  int NbrLzValue;
  // twice the total spin value
  int TotalSpin;

  // array describing each state
  ULONGLONG* StateDescription;
  // array giving maximum Lz value reached for a fermion in a given state
  int* StateHighestBit;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given lzmax sector
  unsigned long LookUpTableMemorySize;
  // shift used in each lzmax sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used lzmax value of the state an the second 
  int** LookUpTable;

  // a table containing ranging from 0 to 2^MaximumSignLookUp - 1
  double* SignLookUpTable;
  // a table containing the mask on the bits to keep for each shift that is requested by sign evaluation
  ULONGLONG* SignLookUpTableMask;
  // number to evalute size of SignLookUpTable
  int MaximumSignLookUp;

  // temporary state used when applying ProdA operator
  ULONGLONG ProdATemporaryState;
  // Lz maximum value associated to temporary state used when applying ProdA operator
  int ProdALzMax;

  // target space for operations leaving the Hilbert-space
  FermionOnSphereWithSpinLong* TargetSpace;

 public:

  // default constructor
  //
  FermionOnSphereWithSpinLong();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalSpin = twce the total spin value
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinLong (int nbrFermions, int totalLz, int lzMax, int totalSpin, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpinLong(const FermionOnSphereWithSpinLong& fermions);

  // destructor
  //
  ~FermionOnSphereWithSpinLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpinLong& operator = (const FermionOnSphereWithSpinLong& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  int GetParticleStatistic();

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  void SetTargetSpace(ParticleOnSphereWithSpin* targetSpace);
  
  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  int GetTargetHilbertSpaceDimension();

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

  // apply creation operator to a word, using the conventions
  // for state-coding and quantum numbers of this space
  // state = word to be acted upon
  // m = Lz value of particle to be added
  // s = spin index of particle to be added (0=down, 1=up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  virtual ULONGLONG Ad (ULONGLONG state, int m, int s, double& coefficient);

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

  // apply a^+_m_s a_m_s operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // sigma = internal degree of freedom label of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdsigmaAsigma (int index, int m, int sigma);

  // apply a^+_m_s a_m_s operator to a given state)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // sigma = internal degree of freedom label of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdsigmaAsigma (long index, int m, int sigma);
  
  // apply a^+_m1_s1 a_m2_s2 operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m1 = index of the creation operator
  // sigma1 = internal degree of freedom label of the creation operator
  // m2 = index of the annihilation operator
  // sigma2 = internal degree of freedom label of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdsigmaAsigma (int index, int m1, int sigma1, int m2, int sigma2, double& coefficient);

  // apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // sigma1 = SU(2) index for the first annihilation operator
  // sigma2 = SU(2) index for the second annihilation operator
  // return value =  multiplicative factor 
  virtual double AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2);

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

  // apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is is 0 for up and 1 for down
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // sigma1 = SU(2) index for the first creation operator
  // sigma2 = SU(2) index for the second creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient);

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

  // apply a^+_n1_d a_n2_u operator to a given state. 
  //
  // index = index of the state on which the operator has to be applied
  // n = first index for annihilation operator (spin up)
  // m = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value =  index of the destination state 
  virtual int AddAu (int index, int m, int n, double& coefficient);

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

  // create an SU(2) state from two U(1) state
  //
  // upState = vector describing the up spin part of the output state
  // upStateSpace = reference on the Hilbert space associated to the up spin part
  // downState = vector describing the down spin part of the output state
  // downStateSpace = reference on the Hilbert space associated to the down spin part  
  // return value = resluting SU(2) state
  virtual RealVector ForgeSU2FromU1(RealVector& upState, FermionOnSphere& upStateSpace, RealVector& downState, FermionOnSphere& downStateSpace);

  // create a U(1) state from an SU(2) state
  //
  // state = vector describing the SU(2) state
  // u1Space = reference on the Hilbert space associated to the U(1) state
  // return value = resulting U(1) state
  virtual RealVector ForgeU1FromSU2(RealVector& state, FermionOnSphere& u1Space);

  // Evaluate the Density Matrix of the spin up fermions in a sector with a fixed lzUp 
  //
  // lzUp = twice total momentum of up fermions.
  // groundstate = reference on the total system groundstate
  // return value = density matrix of the subsystem of spins up fermions.
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrixSpinSeparation (int lzUp, RealVector & groundstate);

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

 protected:

  // factorized code for any a^+_m_x a_n_y operator 
  //
  // index = index of the state on which the operator has to be applied
  // m = global index of the creation operator
  // n = global index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int GenericAdA(int index, int m, int n, double& coefficient);

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(ULONGLONG stateDescription, int lzmax);


  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // totalSpin = number of particles with spin up
  // return value = Hilbert space dimension      
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin);

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // currentLzMax = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // totalSz = spin total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  int OldGenerateStates(int nbrFermions, int lzMax, int totalLz, int totalSz);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion in the state
  // totalLz = momentum total value
  // totalSpin = number of particles with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int lzMax, int totalLz, int totalSpin, long pos);

};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnSphereWithSpinLong::GetParticleStatistic()
{
  return ParticleOnSphereWithSpin::FermionicStatistic;
}

// apply a^+_m_s a_m_s operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// sigma = internal degree of freedom label of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

inline double FermionOnSphereWithSpinLong::AdsigmaAsigma (int index, int m, int sigma)
{
  return ((double) ((this->StateDescription[index] >> ((m << 1) + sigma)) & ((ULONGLONG) 0x1ul)));
}

// apply a^+_m_s a_m_s operator to a given state)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// sigma = internal degree of freedom label of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

inline double FermionOnSphereWithSpinLong::AdsigmaAsigma (long index, int m, int sigma)
{
  return ((double) ((this->StateDescription[index] >> ((m << 1) + sigma)) & ((ULONGLONG) 0x1ul)));
}

// apply a^+_m1_s1 a_m2_s2 operator to a given state
//
// index = index of the state on which the operator has to be applied
// m1 = index of the creation operator
// sigma1 = internal degree of freedom label of the creation operator
// m2 = index of the annihilation operator
// sigma2 = internal degree of freedom label of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int FermionOnSphereWithSpinLong::AdsigmaAsigma (int index, int m1, int sigma1, int m2, int sigma2, double& coefficient)
{
  return this->GenericAdA(index, (m1 << 1) + sigma1, (m2 << 1) + sigma2, coefficient);
}
// factorized code for any a^+_m_x a_n_y operator 
//
// index = index of the state on which the operator has to be applied
// m = global index of the creation operator
// n = global index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int FermionOnSphereWithSpinLong::GenericAdA(int index, int m, int n, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  if ((n > StateHighestBit) || ((State & (((ULONGLONG) 0x1ul) << n)) == ((ULONGLONG) 0x0ul)) )
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = -this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(State >> (n + 64)) & this->SignLookUpTableMask[n + 64]];
  coefficient *= this->SignLookUpTable[(State >> (n + 80))  & this->SignLookUpTableMask[n + 80]];
  coefficient *= this->SignLookUpTable[(State >> (n + 96)) & this->SignLookUpTableMask[n + 96]];
  coefficient *= this->SignLookUpTable[(State >> (n + 112)) & this->SignLookUpTableMask[n + 112]];
#endif
  State &= ~(((ULONGLONG) 0x1ul) << n);
  if (State != ((ULONGLONG) 0x0ul))
    while ((State >> NewLargestBit) == ((ULONGLONG) 0x0ul))
      --NewLargestBit;

  if ((State & (((ULONGLONG) 0x1ul) << m)) != ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
  coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
  coefficient *= this->SignLookUpTable[(State >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
  coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
  coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(State >> (m + 64)) & this->SignLookUpTableMask[m + 64]];
  coefficient *= this->SignLookUpTable[(State >> (m + 80))  & this->SignLookUpTableMask[m + 80]];
  coefficient *= this->SignLookUpTable[(State >> (m + 96)) & this->SignLookUpTableMask[m + 96]];
  coefficient *= this->SignLookUpTable[(State >> (m + 112)) & this->SignLookUpTableMask[m + 112]];
#endif
    }
  State |= ((ULONGLONG) 0x1ul) << m;
  return this->FindStateIndex(State, NewLargestBit);
}

// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(2) index for the first annihilation operator
// sigma2 = SU(2) index for the second annihilation operator
// return value =  multiplicative factor 

inline double FermionOnSphereWithSpinLong::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 1;
  n1 += 1 - sigma1;
  n2 <<= 1;
  n2 += 1 - sigma2;
 if (((this->ProdATemporaryState & (((ULONGLONG) 0x1ul) << n1)) == 0) || ((this->ProdATemporaryState & (((ULONGLONG) 0x1ul) << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n1);
  if (this->ProdATemporaryState != ((ULONGLONG) ((ULONGLONG) 0x0ul)))
    {
      while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
	--this->ProdALzMax;
    }
  else
    this->ProdALzMax = 0;
  return Coefficient;
}

// apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is 0 for up and 1
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// sigma1 = SU(2) index for the first creation operator
// sigma2 = SU(2) index for the second creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int FermionOnSphereWithSpinLong::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  m1 += 1 - sigma1;
  m2 <<= 1;
  m2 += 1 - sigma2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80))  & this->SignLookUpTableMask[m2 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80))  & this->SignLookUpTableMask[m1 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

#endif


