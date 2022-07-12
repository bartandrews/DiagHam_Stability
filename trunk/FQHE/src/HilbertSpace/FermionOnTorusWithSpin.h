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


#ifndef FERMIONONTORUSWITHSPIN_H
#define FERMIONONTORUSWITHSPIN_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"

#include <iostream>


class FermionOnTorusWithSpin :  public FermionOnSphereWithSpin
{

  friend class FermionOnTorusWithSpinAndMagneticTranslations;
  friend class FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations;

 public:

  // constructor with a constraint on total spin momentum and total momentum
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // totalSpinMomentum = twice the total spin momentum to be used as constraint
  // momentumConstraint = index of the momentum orbit
  FermionOnTorusWithSpin (int nbrFermions, int maxMomentum, int totalSpinMomentum, int momentumConstaint);

  // constructor without  constraint on total spin momentum
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // momentumConstraint = index of the momentum orbit
  FermionOnTorusWithSpin (int nbrFermions, int maxMomentum, int momentumConstaint);

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


 private:

  // evaluate Hilbert space dimension for a given total spin momentum
  //
  // nbrFermions = number of fermions
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // nbrSpinUp = number of particles with spin up
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKy, int currentTotalKy, int nbrSpinUp);

  // evaluate Hilbert space dimension for a given total momentum
  //
  // nbrFermions = number of fermions
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKy, int currentTotalKy);

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

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion in the state
  // currentMaxMomentum = momentum maximum value for fermions that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentMomentum = current value of the momentum
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, int pos, int currentMomentum);

};


#endif


