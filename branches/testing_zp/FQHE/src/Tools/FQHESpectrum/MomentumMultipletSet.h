////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of set of momentum multiplets                    //
//                                                                            //
//                        last modification : 31/05/2005                      //
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


#ifndef MOMENTUMMULTIPLETSET_H
#define MOMENTUMMULTIPLETSET_H


#include "config.h"

#include <iostream>


using std::ostream;


class IntegerPolynomial;


class MomentumMultipletSet
{

 protected:

  // twice the maximum momentum
  int MaximumMomentum;

  // array that gives degeneracy associated to each L value (i.e. number of time each L value appeared in the spectrum)
  int* LValues;

 public:

  // default constructor
  //
  // maximumMomentum = twice the maximum momentum
  MomentumMultipletSet (int maximumMomentum = 0);

  // constructor from raw datas
  //
  // maximumMomentum = twice the maximum momentum
  // multiplets = array that contains number of multiplets
  // duplicateFlag = true if the array has to be duplicated
  MomentumMultipletSet (int maximumMomentum, int* multiplets, bool duplicateFlag = false);

  // constructor from an integer polynomial, power indicates Lz value up to a translation, coefficient indicates Lz multiplicty (for example, (1+2q+2q^2+2q^3+q^4) * q^6 
  // corresponds to two multiplets L=1 and L=2)
  //
  // polynomial = reference on the integer polynomial
  MomentumMultipletSet (IntegerPolynomial& polynomial);

  // copy constructor
  //
  // multiplets = set of multiplets to copy 
  MomentumMultipletSet (const MomentumMultipletSet& multiplets);

  // destructor 
  // 
  ~MomentumMultipletSet();

  // assignment
  //
  // multiplets = set of multiplets to assign 
  // return value = reference on current momentum multiplet set
  MomentumMultipletSet& operator = (const MomentumMultipletSet& multiplets);

  // get twice the maximum momentum that is reached in the current set of multiplets
  //
  // return value = maximum momentum  
  int GetMaximumMomentum();

  // get reference on the number of multiplets for a given L value
  //
  // lvalue = twice the L value
  // return value = reference on the corresponding number of multiplets
  int& operator [] (int lvalue);

  // find all multiplets that can obtained when putting a given number of bosons where each of them having the same maximum momentum
  //
  // nbrBosons = number of bosons
  // maximumMomentum = twice the maximum momentum that can have a boson
  // return value = reference on corresponding momentum multiplet set
  MomentumMultipletSet& FindMultipletsForBosons (int nbrBosons, int maximumMomentum);
  
  // find all multiplets that can obtained when putting a given number of fermions where each of them having the same maximum momentum
  //
  // nbrFermions = number of fermions
  // maximumMomentum = twice the maximum momentum that can have a fermion
  // return value = reference on corresponding momentum multiplet set
  MomentumMultipletSet& FindMultipletsForFermions (int nbrFermions, int maximumMomentum);
  
  // get the total number of states associated to the current set of multiplets
  //
  // return value = number of states
  int GetNbrStates();
    
  // add a set of multiplets to another sets of multiplets (i.e. add the number multiplet for each L value)
  // 
  // multiplets = set of multiplets to add 
  // return value = reference on current momentum multiplet set
  MomentumMultipletSet& operator += (const MomentumMultipletSet& multiplets);
  
  // fuse two sets of multiplets according to momentum addition rules
  //
  // multiplets1 = set of first multiplets
  // multiplets2 = set of second multiplets
  // return value = resulting momentum multiplet set
  friend MomentumMultipletSet operator * (MomentumMultipletSet& multiplets1, MomentumMultipletSet& multiplets2); 

  // fuse two sets of multiplets according to momentum addition rules and add result to another set of multiplets (i.e. add the number multiplet for each L value)
  //
  // multiplets1 = first set of multiplets
  // multiplets2 = second set of multiplets
  // return value = reference on current momentum multiplet set
  MomentumMultipletSet& FuseAndAdd (MomentumMultipletSet& multiplets1, MomentumMultipletSet& multiplets2); 

  // resize (aka change maximum momentum value, must be greater than the current one)
  //
  // maximumMomentum = twice the new maximum momentum
  // return value = reference on current momentum multiplet set
  MomentumMultipletSet& Resize(int maximumMomentum);

  // output stream overload (output only half integer momenta if no integer momentum is present, idem for integer momenta)
  //
  // str = reference on the output stream
  // multiplets = set of momentum multiplets to display
  // return value = reference on the output stream
  friend ostream& operator << (ostream& str, const MomentumMultipletSet& multiplets);

  // display mulitplet with respect to Lz instead of L (output only half integer momenta if no integer momentum is present, idem for integer momenta)
  //
  // str = reference on the output stream
  // return value = reference on the output stream
  ostream& PrintLz (ostream& str);

 private:

  // evaluate Hilbert space dimension for bosons with fixed total Lz value and a given L per particle
  //
  // nbrBosons = number of bosons
  // lzMax = two times momentum maximum value for a boson plus one 
  // totalLz = momentum total value plus nbrBosons  * (momentum maximum value for a plusnbrBosons + 1)
  // return value = Hilbert space dimension
  int GetFixedLzBosonHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

  // evaluate Hilbert space dimension for fermions with fixed total Lz value and a given L per particle
  //
  // nbrFermions = number of fermions
  // lzMax = two times momentum maximum value for a fermion plus one 
  // totalLz = momentum total value plus nbrFermions  * (momentum maximum value for a plus nbrFermions + 1)
  // return value = Hilbert space dimension
  int GetFixedLzFermionHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

};

// get twice the maximum momentum that is reached in the current set of multiplets
//
// return value = maximum momentum  

inline int MomentumMultipletSet::GetMaximumMomentum()
{
  return this->MaximumMomentum;
}

// get reference on the number of multiplets for a given L value
//
// lvalue = twice the L value
// return value = reference on the corresponding number of multiplets

inline int& MomentumMultipletSet::operator [] (int lvalue)
{
  return this->LValues[lvalue];
}

#endif
