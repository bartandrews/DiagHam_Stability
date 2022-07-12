////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                  class of fermions on the 4D manifold S2 x S2              //
//               with an exclusion principle along within each sphere         //
//                                                                            //
//                        last modification : 08/12/2016                      //
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


#ifndef FERMIONONS2XS2WITHEXCLUSIONPRINCIPLE_H
#define FERMIONONS2XS2WITHEXCLUSIONPRINCIPLE_H

#include "config.h"
#include "HilbertSpace/FermionOnS2xS2.h"

#include <iostream>



class FermionOnS2xS2WithExclusionPrinciple : public FermionOnS2xS2
{

 protected:

  unsigned long* TemporaryCanonicalArray;
  unsigned long* TemporaryCanonicalArray3x3Block;

 public:

  // default constructor
  // 
  FermionOnS2xS2WithExclusionPrinciple ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrFluxQuanta1 = number of flux quanta for the first sphere
  // nbrFluxQuanta2 = number of flux quanta for the second sphere
  // totalLz1 = total angular momentum for the first sphere
  // totalLz2 = total angular momentum for the second sphere
  // memory = amount of memory granted for precalculations
  FermionOnS2xS2WithExclusionPrinciple (int nbrFermions, int nbrFluxQuanta1, int nbrFluxQuanta2, int totalLz1, int totalLz2, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnS2xS2WithExclusionPrinciple(const FermionOnS2xS2WithExclusionPrinciple& fermions);

  // destructor
  //
  ~FermionOnS2xS2WithExclusionPrinciple ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnS2xS2WithExclusionPrinciple& operator = (const FermionOnS2xS2WithExclusionPrinciple& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // request whether state with given index satisfies a general Pauli exclusion principle
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  virtual bool HasPauliExclusions(int index, int pauliK, int pauliR);

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos);

  // request whether state with given index satisfies a general Pauli exclusion principle
  //
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  // 
  virtual unsigned long FindCanonical(unsigned long state, int xPosition, int yPosition);

  // request whether state with given index satisfies a general Pauli exclusion principle
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  virtual int FindCanonical(unsigned long state, int currentPosition);

  // request whether state with given index satisfies a general Pauli exclusion principle
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  virtual int FindCanonical3x3Block(unsigned long state, int currentPosition);

  // test if a configuration satisfies the core exclusion principle 
  //
  // state = configuration to test
  // return value = true if the configuration satisfies the core exclusion principle
  virtual bool CheckCoreExclusion(unsigned long state);

};


#endif


