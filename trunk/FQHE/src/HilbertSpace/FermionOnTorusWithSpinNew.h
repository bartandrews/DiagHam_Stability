////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of fermion with spin on a torus                  //
//                                                                            //
//                        last modification : 26/11/2007                      //
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


#ifndef FERMIONONTORUSWITHSPINNEW_H
#define FERMIONONTORUSWITHSPINNEW_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSpin.h"

using std::cout;
using std::endl;
using std::dec;
using std::hex;


class FermionOnTorusWithSpinNew :  public ParticleOnTorusWithSpin
{

  // number of fermions
  int NbrFermions;

  // number of fermions with spin up
  int NbrFermionsUp;

  // number of fermions with spin down
  int NbrFermionsDown;
  
  // total number of fermions plus 1
  int IncNbrFermions;

  // total value of Spin
  int TotalSpin;

  // maximum momentum value reached by a fermion
  int MaxMomentum;
  // number of momentum values in a state (= MaxMomentum +1)
  int NbrMomentum;
  // number of different electron orbitals in a state (= 2*(MaxMomentum +1))
  int NbrFermionStates;
  // GCD of MaxMomentum and NbrFermions (momemta are defined modulo MomentumModulo)
  int MomentumModulo;
  // momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)
  int YMomentum;

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
  // array giving maximum momentum value reached for a fermion in a given state
  int* StateHighestBit;
  
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
  // a table containing parity of the sum of 1 bits for all integer ranging from 0 to 2^MaximumSignLookUp - 1 (1 if odd)
  int* NbrParticleLookUpTable;

  // temporary state used when applying ProdA operator
  unsigned long ProdATemporaryState;
  // Highest Bit value associated to temporary state used when applying ProdA operator
  int ProdAHighestBit;
  // Index of the initial state when applying ProdA operator
  int ProdAIndex;


 public:

  
  // basic constructor
  // 
  // nbrFermions= number of fermions
  // totalSpin = twice the total spin value
  // maxMomentum = momentum maximum value for a fermion
  // xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
  // yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)
  
  FermionOnTorusWithSpinNew (int nbrFermions, int totalSpin, int maxMomentum, int yMomentum);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnTorusWithSpinNew(const FermionOnTorusWithSpinNew& fermions);

  // destructor
  //
  ~FermionOnTorusWithSpinNew();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnTorusWithSpinNew& operator = (const FermionOnTorusWithSpinNew& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  int GetParticleStatistic();

  // get momemtum value in the y direction of a given state
  //
  // index = state index
  // return value = state momentum in the y direction
  int GetYMomentumValue(int index);

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

    // apply sum_m au^+_m au_m operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int SumAduAu (int index, double& coefficient);

  // apply sum_m ad^+_m ad_m operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int SumAddAd (int index, double& coefficient);

  // apply au^+_m au_m operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for density operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAu (int index, int m, double& coefficient);

  // apply ad^+_m ad_m operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for density operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAd (int index, int m, double& coefficient);

  // apply au^+_m1 au^+_m2 au_n1 au_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply ad^+_m1 ad^+_m2 ad_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply ad^+_m1 ad^+_m2 au_n1 au_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AddAddAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply au^+_m1 au^+_m2 ad_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AduAduAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply au^+_m1 au^+_m2 au_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AduAduAuAd (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply ad^+_m1 au^+_m2 au_n1 au_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AddAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply ad^+_m1 ad^+_m2 ad_n1 au_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AddAddAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply au^+_m1 ad^+_m2 ad_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AduAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply au^+_m1 ad^+_m2 au_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AduAddAuAd (int index, int m1, int m2, int n1, int n2, double& coefficient);


  // verbose versions
  int AduAduAuAuV (int index, int m1, int m2, int n1, int n2, double& coefficient);


  // apply a^+_m_d a_m_d operator to a given state (only spin down)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  double AddAd (int index, int m);

  // apply a^+_m_u a_m_u operator to a given state  (only spin up)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  double AduAu (int index, int m);

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
  double AuAu (int index, int n1, int n2);
  double AuAuV (int index, int n1, int n2);

  // apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin down)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  double AdAd (int index, int n1, int n2);

  // apply a_n1_u a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  double AuAd (int index, int n1, int n2);

  // apply a_n1_d a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  double AdAu (int index, int n1, int n2);

  
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
  int AduAdu (int m1, int m2, double& coefficient);
  int AduAduV (int m1, int m2, double& coefficient);

  // apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AddAdd (int m1, int m2, double& coefficient);

  // apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AduAdd (int m1, int m2, double& coefficient);  

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

 private:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  int FindStateIndex(unsigned long stateDescription, int maxMomentum);

  
  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // totalSpin = twice z-projection of total spin
  // maxMomentum = momentum maximum value for a fermion
  // return value = Hilbert space dimension
  int EvaluateHilbertSpaceDimension(int nbrFermions, int totalSpin, int maxMomentum);
  
  
  // evaluate Hilbert space dimension using recursive algorithm
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalMomentum = momentum total value
  // totalSpin = number of particles with spin up
  // return value = Hilbert space dimension
  
  long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalMomentum, int totalSpin);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(int memory);

  // generate look-up table associated to sign calculations
  // 
  void GenerateSignLookUpTable();

  // generate all states corresponding to the constraints
  // 
  // return value = hilbert space dimension
  int GenerateStates();
 
  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalMomentum = momentum total value
  // totalSpinUp = number of particles with spin up
  // pos = position in StateDescription array where to store states
  // return value = Hilbert space dimension
  
  long RawGenerateStates(int nbrFermions, int lzMax, int totalMomentum, int totalSpinUp, long pos);
};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnTorusWithSpinNew::GetParticleStatistic()
{
  return ParticleOnTorusWithSpin::FermionicStatistic;
}


// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(2) index for the first annihilation operator
// sigma2 = SU(2) index for the second annihilation operator
// return value =  multiplicative factor 

inline double FermionOnTorusWithSpinNew::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2)
{
  cout << "warning : AsigmaAsigma not defined in FermionOnTorusWithSpinNew" << endl;
  return 0.0;
}

// apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is is 0 for up and 1 for down
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// sigma1 = SU(3) index for the first creation operator
// sigma2 = SU(3) index for the second creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int FermionOnTorusWithSpinNew::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient)
{
  cout << "warning : AdsigmaAdsigma not defined in FermionOnTorusWithSpinNew" << endl;
  return this->HilbertSpaceDimension;
}

#endif // FERMIONONTORUSWITHSPINANDMAGNETICTRANSLATIONS_H


