////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2014 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of bosons on the 4D manifold S2 x S2,              //
// forbidding mutliple orbital occupations and nearest orbital occupations    //
//                                                                            //
//                        last modification : 12/10/2016                      //
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


#ifndef BOSONONS2XS2HARDCORENONEARESTNEIGHBORS_H
#define BOSONONS2XS2HARDCORENONEARESTNEIGHBORS_H


#include "config.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpaceHardcore.h"

#include <iostream>



class BosonOnS2xS2HardcoreNoNearestNeighbors : public BosonOnSquareLatticeMomentumSpaceHardcore
{
  
 protected:

  // number of flux quanta for the first sphere
  int NbrFluxQuanta1;
  // number of flux quanta for the second sphere
  int NbrFluxQuanta2;

  // total angular momentum along z for the first sphere
  int TotalLz1;
  // total angular momentum along z for the second sphere
  int TotalLz2;
  
  
 public:

  // default constructor
  // 
  BosonOnS2xS2HardcoreNoNearestNeighbors ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrFluxQuanta1 = number of flux quanta for the first sphere
  // nbrFluxQuanta2 = number of flux quanta for the second sphere
  // totalLz1 = total angular momentum for the first sphere
  // totalLz2 = total angular momentum for the second sphere
  // memory = amount of memory granted for precalculations
  BosonOnS2xS2HardcoreNoNearestNeighbors (int nbrBosons, int nbrFluxQuanta1, int nbrFluxQuanta2, int totalLz1, int totalLz2, unsigned long memory = 10000000);

  // copy constructor (without duplicating data)
  //
  // space = reference on the hilbert space to copy
  BosonOnS2xS2HardcoreNoNearestNeighbors(const BosonOnS2xS2HardcoreNoNearestNeighbors& space);

  // destructor
  //
  ~BosonOnS2xS2HardcoreNoNearestNeighbors ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnS2xS2HardcoreNoNearestNeighbors & operator = (const BosonOnS2xS2HardcoreNoNearestNeighbors & bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // convert a state such that its components are now expressed in the unnormalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has to be normalized to 1
  // symmetryFactor = if true also remove the symmetry factors
  // return value = converted state
  //  virtual RealVector& ConvertToUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);

  // print a given State using the monomial notation
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  //  virtual ostream& PrintStateMonomial (ostream& Str, long state);

  // request whether state with given index satisfies a general Pauli exclusion principle
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  virtual bool HasPauliExclusions(int index, int pauliK, int pauliR);

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy);

  // generate all states corresponding to the constraints
  // 
  // stateDescription = array that gives each state description
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentFermionicPosition = current fermionic position within the state description
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int currentFermionicPosition, long pos);

  // get a linearized index from the two angular momenta
  //
  // lz = z projection of the angular momentum for the first sphere
  // kz = z projection of the angular momentum for the second sphere
  // return value = linearized index 
  virtual int GetLinearizedIndex(int lz, int kz);

  // get the two angular momenta associated to a given linearized index
  //
  // index = linearized index 
  // lz = reference on the z projection of the angular momentum for the first sphere
  // kz = reference on the  z projection of the angular momentum for the second sphere
  virtual void GetLinearizedIndex(int index, int& lz, int& kz);

  // request whether state with given index satisfies a general Pauli exclusion principle
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  virtual unsigned long FindCanonical(unsigned long state, int xPosition, int yPosition);

};

// get a linearized index from the two angular momenta
//
// lz = z projection of the angular momentum for the first sphere
// kz = z projection of the angular momentum for the second sphere
// return value = linearized index 

inline int BosonOnS2xS2HardcoreNoNearestNeighbors::GetLinearizedIndex(int lz, int kz)
{
  return ((lz * this->NbrSiteY) + kz);
}

// get the two angular momenta associated to a given linearized index
//
// index = linearized index 
// lz = reference on the z projection of the angular momentum for the first sphere
// kz = reference on the  z projection of the angular momentum for the second sphere

inline void BosonOnS2xS2HardcoreNoNearestNeighbors::GetLinearizedIndex(int index, int& lz, int& kz)
{
  lz = index / this->NbrSiteY;
  kz = index % this->NbrSiteY;
}


#endif

