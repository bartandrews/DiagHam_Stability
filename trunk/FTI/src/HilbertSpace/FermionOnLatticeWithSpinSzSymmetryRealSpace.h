////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//                  in real space with the Sz<->-Sz symmetry                  //
//                                                                            //
//                        last modification : 18/08/2014                      //
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


#ifndef FERMIONONLATTICEWITHSPINSZSYMMETRYREALSPACE_H
#define FERMIONONLATTICEWITHSPINSZSYMMETRYREALSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"

#include <iostream>



class FermionOnLatticeWithSpinSzSymmetryRealSpace : public FermionOnSphereWithSpinSzSymmetry
{

 protected:

  // total number of sites
  int NbrSite;
  
  // flag to indicate that the Hilbert space should preserve Sz
  bool SzFlag;

  // indices of the orbitals that are kept when performing an orbital cut
  int* KeptOrbitals;

 public:

  // default constructor
  // 
  FermionOnLatticeWithSpinSzSymmetryRealSpace ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinSzSymmetryRealSpace (int nbrFermions, int nbrSite, bool minusSzParity, unsigned long memory = 10000000);
  
  // basic constructor when Sz is preserved
  // 
  // nbrFermions = number of fermions
  // totalSpin = twice the Sz projection
  // nbrSite = number of sites
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinSzSymmetryRealSpace (int nbrFermions, int totalSpin, int nbrSite, bool minusSzParity, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeWithSpinSzSymmetryRealSpace(const FermionOnLatticeWithSpinSzSymmetryRealSpace& fermions);

  // destructor
  //
  ~FermionOnLatticeWithSpinSzSymmetryRealSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeWithSpinSzSymmetryRealSpace& operator = (const FermionOnLatticeWithSpinSzSymmetryRealSpace& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // find state index from a string
  //
  // stateDescription = string describing the state
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndex(char* stateDescription);
  

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);
  
  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions);
  
  // evaluate Hilbert space dimension with a fixed number of fermions with spin up
  //
  // nbrFermions = number of fermions
  // nbrSpinUp = number of fermions with spin up
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int nbrSpinUp);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentSite, long pos);
  
  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real space
  // nbrSpinUp = number of fermions with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, long pos);

  // generate all states corresponding to the constraints
  // 
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // return value = number of generated states
  virtual long GenerateStatesWithSzSymmetry(bool minusParity);
  

};


#endif


