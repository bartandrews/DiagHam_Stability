////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of fermions on a torus with spin                  //
//                                                                            //
//                        last modification : 10/09/2002                      //
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
#include "HilbertSpace/ParticleOnTorusWithSpin.h"

#include <iostream>


class FermionOnTorusWithSpin :  public ParticleOnTorusWithSpin
{

  // number of fermions
  int NbrFermions;
  // number of fermions plus 1
  int IncNbrFermions;
  // maximum momentum value reached by a fermion
  int MaxMomentum;
  // number of Lz values in a state
  int NbrLzValue;

  // twice the total spin momentum to be used as constraint
  int TotalSpinMomentum;
  // indicate if there is a constraint on the total spin momentum
  bool TotalSpinMomentumFlag; 

  // index of the momentum orbit
  int MomentumConstraint;
  // index of the momentum orbit
  bool MomentumConstraintFlag;

  // array describing each state
  unsigned int* StateDescription;
  // array giving maximum Lz value reached for a fermion in a given state
  int* StateMaxMomentum;

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
  // number to evalute size of SignLookUpTable
  int MaximumSignLookUp;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  FermionOnTorusWithSpin (int nbrFermions, int maxMomentum);

  // constructor with a constraint on total spin momentum
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // totalSpinMomentum = twice the total spin momentum to be used as constraint
  FermionOnTorusWithSpin (int nbrFermions, int maxMomentum, int totalSpinMomentum);

  // constructor with a constraint on total spin momentum and total momentum
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // totalSpinMomentum = twice the total spin momentum to be used as constraint
  // momentumConstraint = index of the momentum orbit
  FermionOnTorusWithSpin (int nbrFermions, int maxMomentum, int totalSpinMomentum, int momentumConstaint);

  // constructor from full datas
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // totalSpinMomentum = twice the total spin momentum to be used as constraint
  // hilbertSpaceDimension = Hilbert space dimension
  // stateDescription = array describing each state
  // stateMaxMomentum = array giving maximum Lz value reached for a fermion in a given state
  // momentumConstraintFlag = flag for momementum constraint
  // momentumConstraint = index of the momentum orbit
 FermionOnTorusWithSpin (int nbrFermions, int maxMomentum, int totalSpinMomentum, int hilbertSpaceDimension, unsigned int* stateDescription, 
			  int* StateMaxMomentum, bool momentumConstraintFlag = false, int momentumConstraint = 0);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnTorusWithSpin(const FermionOnTorusWithSpin& fermions);

  // destructor
  //
  ~FermionOnTorusWithSpin ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnTorusWithSpin& operator = (const FermionOnTorusWithSpin& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  int GetParticleStatistic();

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  AbstractQuantumNumber* GetQuantumNumber (int index);

  // get momemtum value of a given state
  //
  // index = state index
  // return value = state momentum
  int GetMomentumValue(int index);

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
  int SumAudAu (int index, double& coefficient);

  // apply sum_m ad^+_m ad_m operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int SumAddAd (int index, double& coefficient);

  // apply au^+_m au_m operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for density operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AudAu (int index, int m, double& coefficient);

  // apply ad^+_m ad_m operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for density operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AddAd (int index, int m, double& coefficient);

  // apply au^+_m1 au^+_m2 au_n1 au_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AudAudAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

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
  int AudAudAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply au^+_m1 au^+_m2 au_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AudAudAuAd (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply ad^+_m1 au^+_m2 au_n1 au_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AddAudAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

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
  int AudAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply au^+_m1 ad^+_m2 au_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AudAddAuAd (int index, int m1, int m2, int n1, int n2, double& coefficient);

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
  int FindStateIndex(unsigned int stateDescription, int maxMomentum);

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // return value = Hilbert space dimension
  int EvaluateHilbertSpaceDimension(int nbrFermions, int maxMomentum);

  // evaluate Hilbert space dimension for a given total spin momentum
  //
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // momemtum = twice the total spin momentum
  // return value = Hilbert space dimension
  int EvaluateHilbertSpaceDimension(int nbrFermions, int maxMomentum, int momemtum);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(int memory);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // currentMaxMomentum = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, int pos);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion in the state
  // currentMaxMomentum = momentum maximum value for fermions that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentTotalSpinMomentum = total spin momemtum of the fermions that have already been placed
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, int pos, int currentTotalSpinMomentum);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion in the state
  // currentMaxMomentum = momentum maximum value for fermions that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentTotalSpinMomentum = total spin momemtum of the fermions that have already been placed
  // currentMomentum = current value of the momentum
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, int pos, int currentTotalSpinMomentum, int currentMomentum);

};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnTorusWithSpin::GetParticleStatistic()
{
  return ParticleOnTorusWithSpin::FermionicStatistic;
}

#endif


