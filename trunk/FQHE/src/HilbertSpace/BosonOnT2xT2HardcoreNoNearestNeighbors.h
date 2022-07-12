////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            Diagham  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2014 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of bosons on the 4D manifold T2 x T2,              //
// forbidding mutliple orbital occupations and nearest orbital occupations    //
//                                                                            //
//                        last modification : 16/02/2016                      //
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


#ifndef BOSONONT2XT2HARDCORENONEARESTNEIGHBORS_H
#define BOSONONT2XT2HARDCORENONEARESTNEIGHBORS_H


#include "config.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpaceHardcore.h"

#include <iostream>



class BosonOnT2xT2HardcoreNoNearestNeighbors : public BosonOnSquareLatticeMomentumSpaceHardcore
{
  
 protected:

  
 public:

  // default constructor
  // 
  BosonOnT2xT2HardcoreNoNearestNeighbors ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrFluxQuanta1 = number of flux quanta for the first torus
  // nbrFluxQuanta2 = number of flux quanta for the second torus
  // totalKy1 = total momentum along y for the first torus
  // totalKy2 = total momentum along y for the second torus
  // memory = amount of memory granted for precalculations
  BosonOnT2xT2HardcoreNoNearestNeighbors (int nbrBosons, int nbrFluxQuanta1, int nbrFluxQuanta2, int totalKy1, int totalKy2, unsigned long memory = 10000000);

  // copy constructor (without duplicating data)
  //
  // space = reference on the hilbert space to copy
  BosonOnT2xT2HardcoreNoNearestNeighbors(const BosonOnT2xT2HardcoreNoNearestNeighbors& space);

  // destructor
  //
  ~BosonOnT2xT2HardcoreNoNearestNeighbors ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnT2xT2HardcoreNoNearestNeighbors & operator = (const BosonOnT2xT2HardcoreNoNearestNeighbors & bosons);

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

  // get the occupation of a given orbital in a state, not assuming the coordinates are valid, assume periodic bounadary condtions along x
  //
  // index = state index
  // xPosition = orbital x coordinates
  // yPosition = orbital y coordinates
  // return value = orbital occupation
  virtual unsigned long GetSafeOccupationWithPBC(int index, int xPosition, int yPosition);

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // occupationLastOrbital = occupation of the last orbital at the current kx value
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int occupationLastOrbital);

  // generate all states corresponding to the constraints
  // 
  // stateDescription = array that gives each state description
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentFermionicPosition = current fermionic position within the state description
  // occupationLastOrbital = occupation of the last orbital at the current kx value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int currentFermionicPosition, int occupationLastOrbital, long pos);

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

inline int BosonOnT2xT2HardcoreNoNearestNeighbors::GetLinearizedIndex(int lz, int kz)
{
  return ((lz * this->NbrSiteY) + kz);
}

// get the two angular momenta associated to a given linearized index
//
// index = linearized index 
// lz = reference on the z projection of the angular momentum for the first sphere
// kz = reference on the  z projection of the angular momentum for the second sphere

inline void BosonOnT2xT2HardcoreNoNearestNeighbors::GetLinearizedIndex(int index, int& lz, int& kz)
{
  lz = index / this->NbrSiteY;
  kz = index % this->NbrSiteY;
}

// get the occupation of a given orbital in a state, not assuming the coordinates are valid
//
// index = state index
// xPosition = orbital x coordinates
// yPosition = orbital y coordinates
// return value = orbital occupation

inline unsigned long BosonOnT2xT2HardcoreNoNearestNeighbors::GetSafeOccupationWithPBC(int index, int xPosition, int yPosition)
{
  while (xPosition < 0)
    {
      xPosition += this->NbrSiteX;
    }
  while (xPosition >= this->NbrSiteX)
    {
      xPosition -= this->NbrSiteX;
    }
  while (yPosition < 0)
    {
      yPosition += this->NbrSiteY;
    }
  while (yPosition >= this->NbrSiteY)
    {
      yPosition -= this->NbrSiteY;
    }
  return ((this->StateDescription[index] >> ((xPosition * this->NbrSiteY) + yPosition)) & 0x1ul);
}


#endif

