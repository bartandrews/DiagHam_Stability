////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of bosons on sphere for system size such that            //
//            LzMax < 127 (for systems with 128 bit integer support)          //
//           or 63 (on 32 bit systems without 128 bit integer support)        //
//                     and forbidding orbital multiple occupancy              //
//                                                                            //
//                        last modification : 14/02/2017                      //
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


#ifndef BOSONONSPHERELONGHARDCORE_H
#define BOSONONSPHERELONGHARDCORE_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereLong.h"


#include <iostream>



class BosonOnSphereLongHardcore :  public FermionOnSphereLong
{

  friend class BosonOnS2xS2;
  
 protected:

  // pointer to the target space when an index is require after applying basic operation
  BosonOnSphereLongHardcore* TargetSpace;

 public:

  // default constuctor
  //
  BosonOnSphereLongHardcore();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a boson
  // memory = amount of memory granted for precalculations
  BosonOnSphereLongHardcore (int nbrBosons, int totalLz, int lzMax, unsigned long memory = 10000000);

  // constructor using an external array for state description
  // 
  // nbrBosons = number of bosons
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a boson
  // stateDescription = array that gives each state description (data are not duplicated)
  // hilbertSpaceDimension = Hilbert space dimension
  // memory = amount of memory granted for precalculations
  BosonOnSphereLongHardcore (int nbrBosons, int totalLz, int lzMax, ULONGLONG* stateDescription, long hilbertSpaceDimension, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereLongHardcore(const BosonOnSphereLongHardcore& bosons);

  // copy constructor, preserving only some specific states 
  //
  // bosons = reference on the hilbert space to copy to copy
  // nbrPreservedStates = number of preserved states
  // preservedStates = array of flags that indicates if the corresponding state has to be preserved 
  //                   (dimension of the array should the one of the original Hilbert space)
  BosonOnSphereLongHardcore(const BosonOnSphereLongHardcore& bosons, long nbrPreservedStates, bool* preservedStates);

  // destructor
  //
  virtual ~BosonOnSphereLongHardcore ();

  // assignment (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereLongHardcore& operator = (const BosonOnSphereLongHardcore& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic();

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnSphere* targetSpace);

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

  // apply a_n1 a_n2 operator to a given state without keeping it in cache
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AA (int index, int n1, int n2, double& coefficient);

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

  // apply a^+_m1 a^+_m2 operator to the state 
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAd (int index, int m1, int m2, double& coefficient);

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
  virtual ULONGLONG Ad (ULONGLONG state, int m, double& coefficient);
  
  // apply a_n  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
  //
  // index = index of the state on which the operator has to be applied
  // n = index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A (int index, int n);

  // apply a_n  operator to a given state. 
  //
  // index = index of the state on which the operator has to be applied
  // n = index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value =  index of the resulting state 
  virtual int A (int index, int n, double& coefficient);

  // apply a^+_n  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
  //
  // index = index of the state on which the operator has to be applied
  // n = index for annihilation operator
  // return value =  multiplicative factor 
  virtual double Ad (int index, int n);

  // apply a_n operator to the state produced using the A or Ad method (without destroying it)
  //
  // n = first index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int A (int n, double& coefficient);

  // apply a^+_m operator to the state produced using the A or Ad method (without destroying it)
  //
  // m = first index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad (int m, double& coefficient);


 protected:

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory = 10000000);

  // generate look-up table for sign calculation
  // 
  virtual void GenerateSignLookUpTable();

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a boson in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);

};


// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnSphereLongHardcore::GetParticleStatistic()
{
  return ParticleOnSphere::BosonicStatistic;
}

#endif
