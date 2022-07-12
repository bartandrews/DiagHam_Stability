////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of fermions on sphere without fixed total Lz             //
//                                                                            //
//                        last modification : 19/04/2010                      //
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


#ifndef FERMIONONSPHEREFULL_H
#define FERMIONONSPHEREFULL_H


#include "config.h"
#include "HilbertSpace/FermionOnSphere.h"

#include <iostream>



class FermionOnSphereFull :  public FermionOnSphere
{

  friend class BosonOnSphereFullShort;
  friend class FermionOnTorus;

 protected:

  // array that contains the total Lz of each state
  int* TotalLzValues;

 public:

  // default constructor
  //
  FermionOnSphereFull();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // lzMax = maximum Lz value reached by a fermion
  // memory = amount of memory granted for precalculations
  // referenceState = array that describes the reference state to start from
  FermionOnSphereFull (int nbrFermions, int lzMax, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereFull(const FermionOnSphereFull& fermions);

  // destructor
  //
  virtual ~FermionOnSphereFull ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereFull& operator = (const FermionOnSphereFull& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get Lz component of a component
  //
  // j = index of the component in Hilbert space
  // return value = twice the Lz component
  virtual int GetLzValue(int j = 0);

  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  double AdA (int index, int m);

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 

  int AdA (int index, int m, int n, double& coefficient);

  // convert the vector with a given Lz to the full space (all Lz components)
  // inputState = input vector
  // inputSpace = input Hilbert space with given Lz
  // return value = vector in the full Hilbert space

  void ConvertToAllLz (ComplexVector& inputState, ParticleOnSphere* inputSpace, ComplexVector& outputState);


 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax);

  // generate all states (i.e. all possible skew symmetric polynomials)
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // currentLzMax = momentum maximum value for fermions that are still to be placed
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int lzMax, int currentLzMax, long pos);

  // carefully test whether state is in Hilbert-space and find corresponding state index
  //
  // stateDescription = unsigned integer describing the state
  // highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
  // return value = corresponding index, or dimension of space, if not found
  virtual int FindStateIndex(unsigned long stateDescription, int lzMax);

};

// get Lz component of a component
//
// j = index of the component in Hilbert space
// return value = twice the Lz component

inline int FermionOnSphereFull::GetLzValue(int j)
{
  return this->TotalLzValues[j];
}


#endif


