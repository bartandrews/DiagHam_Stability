////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of fermions on lattice in real space              //
//             where only the parity of the particle number is fixed          //
//                                                                            //
//                        last modification : 10/08/2015                      //
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


#ifndef FERMIONONLATTICEREALSPACEFIXEDPARITY_H
#define FERMIONONLATTICEREALSPACEFIXEDPARITY_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeRealSpace.h"

#include <iostream>



class FermionOnLatticeRealSpaceFixedParity : public FermionOnLatticeRealSpace
{


 friend class FermionOnLatticeRealSpaceAnd2DTranslation;

 protected:

  // parity of the particle number (0 for even or 1 for odd)
  int Parity;

  // total number of sites
  int NbrSite;
  
 public:

  // default constructor
  // 
  FermionOnLatticeRealSpaceFixedParity ();

  // basic constructor
  // 
  // nbrSite = number of sites
  // parity = parity of the particle number (0 for even or 1 for odd)
  FermionOnLatticeRealSpaceFixedParity (int nbrSite, int parity);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeRealSpaceFixedParity(const FermionOnLatticeRealSpaceFixedParity& fermions);

  // destructor
  //
  ~FermionOnLatticeRealSpaceFixedParity ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeRealSpaceFixedParity& operator = (const FermionOnLatticeRealSpaceFixedParity& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // find state index from a string
  //
  // stateDescription = string describing the state
  // return value = corresponding index, -1 if an error occured
  // virtual int FindStateIndex(char* stateDescription);

  // evaluate the tensor product of three states and apply a Gutzwiller projection
  //
  // state1 = reference on the first state 
  // space1 = Hilbert space associated to the first space
  // state2 = reference on the second state 
  // space2 = Hilbert space associated to the second space
  // state3 = reference on the third state 
  // space3 = Hilbert space associated to the third space
  virtual void TripleTensorProductAndGutzwillerProjection (ComplexVector& state1, FermionOnLatticeRealSpaceFixedParity* space1,
							   ComplexVector& state2, FermionOnLatticeRealSpaceFixedParity* space2,
							   ComplexVector& state3, FermionOnLatticeRealSpaceFixedParity* space3);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);
  
};

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

inline int FermionOnLatticeRealSpaceFixedParity::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  return ((int) (stateDescription >> 1));    
}
  

#endif


