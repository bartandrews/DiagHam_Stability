////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//       in real space with translation invariance in two directions          //
//                          supporting up to 64 sites                         //
//                                                                            //
//                                                                            //
//                        last modification : 16/06/2016                      //
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


#ifndef FERMIONONLATTICEWITHSPINREALSPACEAND2DTRANSLATIONLONG_H
#define FERMIONONLATTICEWITHSPINREALSPACEAND2DTRANSLATIONLONG_H

#include "config.h"
#include "HilbertSpace/FermionOnTorusWithSpinAndMagneticTranslationsLong.h"

#include <iostream>


#ifdef __128_BIT_LONGLONG__
#define FERMION_LATTICE_REALSPACE_SU2_SZ_MASK_LONG ((((ULONGLONG) (0x5555555555555555ul)) << 64) | ((ULONGLONG) (0x5555555555555555ul)))
#else
#define FERMION_LATTICE_REALSPACE_SU2_SZ_MASK_LONG ((((ULONGLONG) (0x55555555ul)) << 32) | ((ULONGLONG) (0x55555555ul)))
#endif


class FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong : public FermionOnTorusWithSpinAndMagneticTranslationsLong
{

  friend class FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong;

 protected:

  // total number of sites
  int NbrSite;
  // number of sites per unit cell
  int NbrSitePerUnitCell;
  
  // flag to indicate that the Hilbert space should preserve Sz
  bool SzFlag;

  // number of momentum sectors in the x direction 
  int MaxXMomentum;
  // bit shift that has to applied to perform a translation in the x direction 
  int StateXShift;
  // binary mask for the StateXShift first bits 
  ULONGLONG XMomentumMask;
  // bit shift to apply to move the first StateXShift bits at the end of a state description
  int ComplementaryStateXShift;

  // number of momentum sectors in the y direction 
  int MaxYMomentum;
  // bit shift that has to applied to perform a translation in the y direction 
  int StateYShift;
  // binary mask for the StateYShift first bits 
  ULONGLONG YMomentumMask;
  // binary mask for the StateYShift first bits of each group
  ULONGLONG YMomentumFullMask;
  // binary mask for the ~YMomentumFullMask
  ULONGLONG ComplementaryYMomentumFullMask;
  // bit shift to apply to move the first StateYShift bits at the end of a state description
  int ComplementaryStateYShift;
  // number of bits that are related by a translation along the y direction 
  int YMomentumBlockSize;
  // binary mask corresponding to YMomentumBlockSize
  ULONGLONG YMomentumBlockMask;
  // number of independant blockse related by translations in the y direction 
  int NbrYMomentumBlocks;

  // parity of the number of fermions, 0x1ul if even, 0x0ul if odd
  ULONGLONG NbrFermionsParity;

 protected:
    
  // target space for operations leaving the Hilbert space
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong* TargetSpace;

 public:

  // default constructor
  // 
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong (int nbrFermions, int nbrSite, int xMomentum, int maxXMomentum,
							 int yMomentum, int maxYMomentum, unsigned long memory = 10000000);
  
  // basic constructor when Sz is preserved
  // 
  // nbrFermions = number of fermions
  // totalSpin = twice the value of Sz
  // nbrSite = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong (int nbrFermions, int totalSpin, int nbrSite, int xMomentum, int maxXMomentum,
							 int yMomentum, int maxYMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong(const FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong& fermions);

  // destructor
  //
  ~FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong& operator = (const FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnSphereWithSpin* targetSpace);
  
  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  virtual int GetTargetHilbertSpaceDimension();

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

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

  // apply a^+_m_u a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AduAu (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // apply a^+_m_d a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AddAd (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // apply a^+_m_u a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation/annihilation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AduAd (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // apply a^+_m_d a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AddAu (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

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

  // apply a_n1_u a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
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
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AduAdu (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AddAdd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AduAdd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // apply a^+_m_u operator to the state produced using Au method (without destroying it)
  //
  // m = first index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int Adu (int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // apply a^+_m_d operator to the state produced using Au method (without destroying it)
  //
  // m = first index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int Add (int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // find state index from a string
  //
  // stateDescription = string describing the state
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndex(char* stateDescription);
  
  // convert a state defined in the real space basis into a state in the (Kx,Ky) basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertToKxKyBasis(ComplexVector& state, ParticleOnSphere* space);

  // convert a state defined in the (Kx,Ky) basis into a state in the real space basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertFromKxKyBasis(ComplexVector& state, ParticleOnSphere* space);

  // convert a given state from a given  n-body basis basis to another one
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis where state is defined
  // return value = converted vector
  virtual ComplexVector ConvertToNbodyBasis(ComplexVector& state, FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong* nbodyBasis);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = subsystem momentum along the x direction
  // kySector = subsystem momentum along the x direction
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in given momentum and Sz sectors.
  //
  // nbrParticleSector = number of particles that belong to the subsytem
  // kxSector = subsystem momentum along the x direction
  // kySector = subsystem momentum along the x direction
  // szSector  = twice the total Sz value of the subsytem
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int szSector, int kxSector, int kySector, 
									 ComplexVector& groundState, AbstractArchitecture* architecture = 0);
  
  // get the linearized index (e.g. used for the creation/annihilation operators) from the lattice position, sanitizing the input data first
  // 
  // xPosition = x coordinate of the unit cell
  // yPosition = y coordinate of the unit cell
  // orbitalIndex = index of the orbital within the unit cell
  // return value = linearized index
  virtual int GetLinearizedIndexSafe(int xPosition, int yPosition, int orbitalIndex);

  // get the linearized index (e.g. used for the creation/annihilation operators) from the lattice position
  // 
  // xPosition = x coordinate of the unit cell
  // yPosition = y coordinate of the unit cell
  // orbitalIndex = index of the orbital within the unit cell
  // return value = linearized index
  virtual int GetLinearizedIndex(int xPosition, int yPosition, int orbitalIndex);

  // get the lattice position from the linearized index (e.g. used for the creation/annihilation operators)
  // 
  // index = linearized index
  // xPosition = reference on the x coordinate of the unit cell
  // yPosition =reference on the  y coordinate of the unit cell
  // orbitalIndex = reference on the index of the orbital within the unit cell
  virtual void GetLinearizedIndex(int index, int& xPosition, int& yPosition, int& orbitalIndex);

  // get the momentum along the x axis
  // 
  // return avlue = momentum along the x axis
  virtual int GetKxMomentum();

  // get the momentum along the y axis
  // 
  // return avlue = momentum along the y axis
  virtual int GetKyMomentum();

  // get the maximum momentum along the x axis (i.e. the number of momentum sectors)
  // 
  // return avlue = maximum momentum along the x axis
  virtual int GetMaxXMomentum();
  
  // get the maximum momentum along the y axis (i.e. the number of momentum sectors)
  // 
  // return avlue = maximum momentum along the y axis
  virtual int GetMaxYMomentum();

  // Apply the Sz operator to flip all the spins and compute the fermonic sign if any
  //
  // index = index of the initial state  on state description
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state  
  virtual int ApplySzSymmetry (int index, double& coefficient, int nbrTranslationX, int nbrTranslationY);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(ULONGLONG stateDescription, int maxMomentum);

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions);
  
  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // nbrSpinUp = number of fermions with spin up
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int nbrSpinUp);

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();

  // generate all states corresponding to the constraints (core part of the method)
  //
  // return value = Hilbert space dimension
  virtual long CoreGenerateStates();

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int currentSite, long pos);
  
  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // nbrSpinUp = number of fermions with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, long pos);

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

  // apply a^+_m_sigma a_n_sigma operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator including the orbital and the spin index
  // n = index of the annihilation operator including the orbital and the spin index
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AdsigmaAsigma (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // apply a_n1_sigma a_n2_sigma operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdsigmaAdsigma call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin up)
  // return value =  multiplicative factor 
  virtual double AsigmaAsigma (int index, int n1, int n2);

  // apply a^+_m1_sigma a^+_m2_sigma operator to the state produced using AsigmaAsigma method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AdsigmaAdsigma (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // apply a^+_m1_sigma a^+_m2_sigma operator to the state, assuming a different target space
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int AdsigmaAdsigma (int index, int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // apply a^+_m1_u a^+_m2_u operator to the state, assuming a different target space
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int AduAdu (int index, int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // apply a^+_m1_d a^+_m2_d operator to the state, assuming a different target space 
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int AddAdd (int index, int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // apply a^+_m1_u a^+_m2_d operator to the state, assuming a different target space 
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int AduAdd (int index, int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // apply a^+_m1_d a^+_m2_u operator to the state, assuming a different target space 
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int AddAdu (int index, int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = reference on the state that has been produced with the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state  
  virtual int SymmetrizeAdAdResult(ULONGLONG& state, double& coefficient, 
				   int& nbrTranslationX, int& nbrTranslationY);

  // factorized code that is used to symmetrize the result of any operator action when the target space is different from the current space
  //
  // state = reference on the state that has been produced with the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state  
  int SymmetrizeAdAdResultTarget(ULONGLONG& state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual ULONGLONG FindCanonicalForm(ULONGLONG stateDescription, int& nbrTranslationX, int& nbrTranslationY);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(ULONGLONG stateDescription);

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  virtual int FindOrbitSize(ULONGLONG stateDescription);

  // apply a single translation in the x direction for a state description
  //
  // stateDescription = reference on the state description
  virtual void ApplySingleXTranslation(ULONGLONG& stateDescription);

  // apply a single translation in the y direction for a state description
  //
  // stateDescription = reference on the state description  
  virtual void ApplySingleYTranslation(ULONGLONG& stateDescription);

  // get the fermonic sign when performing a single translation in the x direction on a state description, and apply the single translation
  //
  // stateDescription = reference on state description
  // return value = 0 if the sign is +1, 1 if the sign is -1
  virtual ULONGLONG GetSignAndApplySingleXTranslation(ULONGLONG& stateDescription);

  // get the fermonic sign when performing a single translation in the y direction on a state description, and apply the single translation
  //
  // stateDescription = reference on state description
  // return value = 0 if the sign is +1, 1 if the sign is -1
  virtual ULONGLONG GetSignAndApplySingleYTranslation(ULONGLONG& stateDescription);

  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in source Hilbert space
  // nbrIndex = number of indices to consider in source Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnTorusWithSpinAndMagneticTranslations* complementaryHilbertSpace,  
								  ParticleOnTorusWithSpinAndMagneticTranslations* destinationHilbertSpace,
								  ComplexVector& groundState, HermitianMatrix* densityMatrix);

};


// apply a^+_m_u a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AduAu (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->AdsigmaAsigma(index, (m << 1) + 1, (n << 1) + 1, coefficient, nbrTranslationX, nbrTranslationY);
}
  
// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AddAd (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->AdsigmaAsigma(index, (m << 1), (n << 1), coefficient, nbrTranslationX, nbrTranslationY);
}
  
// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AduAd (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->AdsigmaAsigma(index, (m << 1) + 1, (n << 1), coefficient, nbrTranslationX, nbrTranslationY);
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AddAu (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->AdsigmaAsigma(index, (m << 1), (n << 1) + 1, coefficient, nbrTranslationX, nbrTranslationY);
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 

inline double FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AuAu (int index, int n1, int n2)
{
  return this->AsigmaAsigma(index, (n1 << 1) + 1, (n2 << 1) + 1);
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

inline double FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AdAd (int index, int n1, int n2)
{
  return this->AsigmaAsigma(index, (n1 << 1), (n2 << 1));
}

// apply a_n1_u a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

inline double FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AuAd (int index, int n1, int n2)
{
  return this->AsigmaAsigma(index, (n1 << 1) + 1, (n2 << 1));
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AduAdu (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->AdsigmaAdsigma((m1 << 1) + 1, (m2 << 1) + 1, coefficient, nbrTranslationX, nbrTranslationY);
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AddAdd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->AdsigmaAdsigma((m1 << 1), (m2 << 1), coefficient, nbrTranslationX, nbrTranslationY);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AduAdd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->AdsigmaAdsigma((m1 << 1) + 1, (m2 << 1), coefficient, nbrTranslationX, nbrTranslationY);
}

// apply a^+_m1_u a^+_m2_u operator to the state, assuming a different target space
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AduAdu (int index, int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->AdsigmaAdsigma(index, (m1 << 1) + 1, (m2 << 1) + 1, coefficient, nbrTranslationX, nbrTranslationY);
}

// apply a^+_m1_d a^+_m2_d operator to the state, assuming a different target space 
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AddAdd (int index, int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->AdsigmaAdsigma(index, (m1 << 1), (m2 << 1), coefficient, nbrTranslationX, nbrTranslationY);
}

// apply a^+_m1_u a^+_m2_d operator to the state, assuming a different target space 
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AduAdd (int index, int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->AdsigmaAdsigma(index, (m1 << 1) + 1, (m2 << 1), coefficient, nbrTranslationX, nbrTranslationY);
}
  
// apply a^+_m1_d a^+_m2_u operator to the state, assuming a different target space 
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AddAdu (int index, int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->AdsigmaAdsigma(index, (m1 << 1), (m2 << 1) + 1, coefficient, nbrTranslationX, nbrTranslationY);
}
  
// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = index of the destination state  

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::SymmetrizeAdAdResult(ULONGLONG& state, double& coefficient, 
										   int& nbrTranslationX, int& nbrTranslationY)
{
  state = this->FindCanonicalForm(state, nbrTranslationX, nbrTranslationY);
  int TmpMaxMomentum = 2 * this->NbrSite - 1;
  while ((state >> TmpMaxMomentum) == ((ULONGLONG) 0x0ul))
    {
      --TmpMaxMomentum;
    }
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      nbrTranslationX = (this->MaxXMomentum - nbrTranslationX) % this->MaxXMomentum;
      nbrTranslationY = (this->MaxYMomentum - nbrTranslationY) % this->MaxYMomentum;
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> ((nbrTranslationY * this->MaxXMomentum) + nbrTranslationX)) & 0x1ul))); 
    }
  return TmpIndex;
}

// factorized code that is used to symmetrize the result of any operator action when the target space is different from the current space
//
// state = reference on the state that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = index of the destination state  

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::SymmetrizeAdAdResultTarget(ULONGLONG& state, double& coefficient, 
											 int& nbrTranslationX, int& nbrTranslationY)
{
  state = this->TargetSpace->FindCanonicalForm(state, nbrTranslationX, nbrTranslationY);
  int TmpMaxMomentum = 2 * this->NbrSite - 1;
  while ((state >> TmpMaxMomentum) == ((ULONGLONG) 0x0ul))
    {
      --TmpMaxMomentum;
    }
  int TmpIndex = this->TargetSpace->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->TargetSpace->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][this->TargetSpace->NbrStateInOrbit[TmpIndex]];
      nbrTranslationX = (this->MaxXMomentum - nbrTranslationX) % this->MaxXMomentum;
      nbrTranslationY = (this->MaxYMomentum - nbrTranslationY) % this->MaxYMomentum;
      coefficient *= 1.0 - (2.0 * ((double) ((this->TargetSpace->ReorderingSign[TmpIndex] >> ((nbrTranslationY * this->MaxXMomentum) + nbrTranslationX)) & 0x1ul))); 
    }
  return TmpIndex;
}


// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline ULONGLONG FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::FindCanonicalForm(ULONGLONG stateDescription, int& nbrTranslationX, int& nbrTranslationY)
{
  ULONGLONG CanonicalState = stateDescription;
  ULONGLONG stateDescriptionReference = stateDescription;  
  ULONGLONG TmpStateDescription;  
  nbrTranslationX = 0;
  nbrTranslationY = 0;
  TmpStateDescription = stateDescription;
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);      
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslationX = n;	      
	  nbrTranslationY = 0;	      
	}
    }
  for (int m = 1; m < this->MaxYMomentum; ++m)
    {
      this->ApplySingleYTranslation(stateDescription);      
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslationX = 0;	      
	  nbrTranslationY = m;	      
	}
      TmpStateDescription = stateDescription;
      for (int n = 1; n < this->MaxXMomentum; ++n)
	{
	  this->ApplySingleXTranslation(TmpStateDescription);      
	  if (TmpStateDescription < CanonicalState)
	    {
	      CanonicalState = TmpStateDescription;
	      nbrTranslationX = n;	      
	      nbrTranslationY = m;	      
	    }
	}
    }
  return CanonicalState;
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::TestMomentumConstraint(ULONGLONG stateDescription)
{
  ULONGLONG TmpStateDescription = stateDescription;
  ULONGLONG TmpStateDescription2 = stateDescription;
  int XSize = 1;
  ULONGLONG TmpSign = this->GetSignAndApplySingleXTranslation(TmpStateDescription);   
  ULONGLONG TmpSign2 = 0x0ul;
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      TmpSign ^= this->GetSignAndApplySingleXTranslation(TmpStateDescription);      
    }
  if ((((this->XMomentum * XSize) + ((((int) TmpSign) * this->MaxXMomentum) >> 1)) % this->MaxXMomentum) != 0)
    return false;
  int YSize = this->MaxYMomentum;
  int TmpXSize = 0;
  TmpSign = 0x0ul;
  TmpStateDescription2 = stateDescription;
  for (int m = 1; m < YSize; ++m)
    {
      TmpSign ^= this->GetSignAndApplySingleYTranslation(TmpStateDescription2); 
      TmpSign2 = TmpSign;
      TmpStateDescription = TmpStateDescription2;
      TmpXSize = 0;
      while ((TmpXSize < XSize) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  TmpSign2 ^= this->GetSignAndApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize)
	{
	  YSize = m;
	}
      else
	{
	  TmpXSize = 0;
	}
    } 
  if (YSize == this->MaxYMomentum)
    {
      TmpSign ^= this->GetSignAndApplySingleYTranslation(TmpStateDescription2); 
      TmpSign2 = TmpSign;
    }
  if ((((2 * this->YMomentum * YSize * this->MaxXMomentum)
	+ (2 * this->XMomentum * TmpXSize * this->MaxYMomentum)
	+ (((int) TmpSign2) * this->MaxXMomentum * this->MaxYMomentum)) % (2 * this->MaxXMomentum * this->MaxYMomentum)) != 0)
    return false;
  return true;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::FindOrbitSize(ULONGLONG stateDescription)
{
  ULONGLONG TmpStateDescription = stateDescription;
  ULONGLONG TmpStateDescription2 = stateDescription;
  int XSize = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }
  int YSize = this->MaxYMomentum;
  for (int m = 1; m < YSize; ++m)
    {
      this->ApplySingleYTranslation(stateDescription); 
      TmpStateDescription = TmpStateDescription2;
      int TmpXSize = 0;
      while ((TmpXSize < XSize) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  this->ApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize)
	{
	  YSize = m;
	}
    }
  return (XSize * YSize);
}

// apply a single translation in the x direction for a state description
//
// stateDescription = reference on the state description

inline void FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::ApplySingleXTranslation(ULONGLONG& stateDescription)
{
  stateDescription = (stateDescription >> this->StateXShift) | ((stateDescription & this->XMomentumMask) << this->ComplementaryStateXShift);
}

// apply a single translation in the y direction for a state description
//
// stateDescription = reference on the state description

inline void FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::ApplySingleYTranslation(ULONGLONG& stateDescription)
{
  stateDescription = (((stateDescription & this->ComplementaryYMomentumFullMask) >> this->StateYShift) | 
		      ((stateDescription & this->YMomentumFullMask) << this->ComplementaryStateYShift));
}

// get the fermonic sign when performing a single translation in the x direction on a state description, and apply the single translation
//
// stateDescription = reference on state description
// return value = 0 if the sign is +1, 1 if the sign is -1

inline ULONGLONG FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::GetSignAndApplySingleXTranslation(ULONGLONG& stateDescription)
{
  ULONGLONG TmpSign =  stateDescription >> this->StateXShift;
  stateDescription = TmpSign | ((stateDescription & this->XMomentumMask) << this->ComplementaryStateXShift);
#ifdef __128_BIT_LONGLONG__
  TmpSign ^= (TmpSign >> 64);
#endif
  TmpSign ^= (TmpSign >> 32);
  TmpSign ^= (TmpSign >> 16);
  TmpSign ^= (TmpSign >> 8);
  TmpSign ^= (TmpSign >> 4);
  TmpSign ^= (TmpSign >> 2);
  TmpSign ^= (TmpSign >> 1);
  TmpSign &= this->NbrFermionsParity;
  return TmpSign;
}

// get the fermonic sign when performing a single translation in the y direction on a state description, and apply the single translation
//
// stateDescription = reference on state description
// return value = 0 if the sign is +1, 1 if the sign is -1

inline ULONGLONG FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::GetSignAndApplySingleYTranslation(ULONGLONG& stateDescription)
{
  ULONGLONG TmpState = ((ULONGLONG) 0x0ul);
  ULONGLONG TmpSign =  ((ULONGLONG) 0x0ul);
  ULONGLONG TmpSign2;
  ULONGLONG TmpSign3;
  for (int i = 0; i < this->NbrYMomentumBlocks; ++i)
    {
      TmpSign2 = (stateDescription & this->YMomentumBlockMask) >> this->StateYShift;
      TmpSign3 = (stateDescription & this->YMomentumMask) << this->ComplementaryStateYShift;
      TmpState |= (TmpSign2 | TmpSign3) << (this->YMomentumBlockSize * i);
#ifdef __128_BIT_LONGLONG__
      TmpSign2 ^= (TmpSign2 >> 64);
#endif
      TmpSign2 ^= (TmpSign2 >> 32);
      TmpSign2 ^= (TmpSign2 >> 16);
      TmpSign2 ^= (TmpSign2 >> 8);
      TmpSign2 ^= (TmpSign2 >> 4);
      TmpSign2 ^= (TmpSign2 >> 2);
      TmpSign2 ^= (TmpSign2 >> 1);
#ifdef __128_BIT_LONGLONG__
      TmpSign3 ^= (TmpSign3 >> 64);
#endif
      TmpSign3 ^= (TmpSign3 >> 32);
      TmpSign3 ^= (TmpSign3 >> 16);
      TmpSign3 ^= (TmpSign3 >> 8);
      TmpSign3 ^= (TmpSign3 >> 4);
      TmpSign3 ^= (TmpSign3 >> 2);
      TmpSign3 ^= (TmpSign3 >> 1);
      TmpSign2 *= TmpSign3;
      TmpSign2 &= ((ULONGLONG) 0x1ul);
      TmpSign ^= TmpSign2;
      stateDescription >>= this->YMomentumBlockSize;
    }
  stateDescription = TmpState;
  return TmpSign;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->GetHilbertSpaceDimension();
}

// get the linearized index (e.g. used for the creation/annihilation operators) from the lattice position, sanitizing the input data first
// 
// xPosition = x coordinate of the unit cell
// yPosition = y coordinate of the unit cell
// orbitalIndex = index of the orbital within the unit cell
// return value = linearized index

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::GetLinearizedIndexSafe(int xPosition, int yPosition, int orbitalIndex)
{
  orbitalIndex %= this->NbrSitePerUnitCell;
  if (orbitalIndex < 0)
    orbitalIndex +=  this->NbrSitePerUnitCell;
  xPosition %= this->MaxXMomentum;
  if (xPosition < 0)
    xPosition +=  this->MaxXMomentum;
  yPosition %= this->MaxYMomentum;
  if (yPosition < 0)
    yPosition +=  this->MaxYMomentum;
  return this->GetLinearizedIndex(xPosition, yPosition, orbitalIndex); 
}

// get the linearized index (e.g. used for the creation/annihilation operators) from the lattice position
// 
// xPosition = x coordinate of the unit cell
// yPosition = y coordinate of the unit cell
// orbitalIndex = index of the orbital within the unit cell
// return value = linearized index

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::GetLinearizedIndex(int xPosition, int yPosition, int orbitalIndex)
{
  return (((xPosition * this->MaxYMomentum) + yPosition) * this->NbrSitePerUnitCell) + orbitalIndex;
}

// get the lattice position from the linearized index (e.g. used for the creation/annihilation operators)
// 
// index = linearized index
// xPosition = reference on the x coordinate of the unit cell
// yPosition =reference on the  y coordinate of the unit cell
// orbitalIndex = reference on the index of the orbital within the unit cell

inline void FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::GetLinearizedIndex(int index, int& xPosition, int& yPosition, int& orbitalIndex)
{
  orbitalIndex = index % this->NbrSitePerUnitCell;
  index /= this->NbrSitePerUnitCell;
  xPosition = index / this->MaxYMomentum;
  yPosition = index % this->MaxYMomentum;
}

// get the momentum along the x axis
// 
// return avlue = momentum along the x axis

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::GetKxMomentum()
{
  return this->XMomentum;
}

// get the momentum along the y axis
// 
// return avlue = momentum along the y axis

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::GetKyMomentum()
{
  return this->YMomentum;
}

// get the maximum momentum along the x axis (i.e. the number of momentum sectors)
// 
// return avlue = maximum momentum along the x axis

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::GetMaxXMomentum()
{
  return this->MaxXMomentum;
}
  
// get the maximum momentum along the y axis (i.e. the number of momentum sectors)
// 
// return avlue = maximum momentum along the y axis

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::GetMaxYMomentum()
{
  return this->MaxYMomentum;
}
  
// apply a^+_m_d a_m_d operator to a given state (only spin down)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

inline double FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AddAd (int index, int m)
{
  return this->FermionOnTorusWithSpinAndMagneticTranslationsLong::AddAd (index, m);
}

// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

inline double FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::AduAu (int index, int m)
{
  return this->FermionOnTorusWithSpinAndMagneticTranslationsLong::AduAu (index, m);
}

// Apply the Sz operator to flip all the spins and compute the fermonic sign if any
//
// index = index of the initial state  on state description
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = index of the destination state  

inline int FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::ApplySzSymmetry (int index, double& coefficient, int nbrTranslationX, int nbrTranslationY)
{
  
  ULONGLONG TmpState = this->StateDescription[index];
  ULONGLONG TmpState2 = ((TmpState >> 1) ^ TmpState) & FERMION_LATTICE_REALSPACE_SU2_SZ_MASK_LONG;
  TmpState2 |= TmpState2 << 1;
  TmpState2 ^= TmpState; 

  TmpState &= (TmpState >> 1);
  TmpState &= FERMION_LATTICE_REALSPACE_SU2_SZ_MASK_LONG;
#ifdef __128_BIT_LONGLONG__
  TmpState ^= (TmpState >> 64);
#endif
  TmpState ^= (TmpState >> 32);
  TmpState ^= (TmpState >> 16);
  TmpState ^= (TmpState >> 8);
  TmpState ^= (TmpState >> 4);
  TmpState ^= (TmpState >> 2);
  coefficient = 1.0;
  if ((TmpState & ((ULONGLONG) 0x1ul)) != ((ULONGLONG) 0x0ul))
    coefficient = -1.0;
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  return this->SymmetrizeAdAdResult(TmpState2, coefficient, nbrTranslationX, nbrTranslationY);
}

#endif


