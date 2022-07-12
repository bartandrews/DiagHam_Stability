////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//       in real space with translation invariance in one direction           //
//                          and the Sz<->-Sz symmetry                         //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                        last modification : 30/07/2014                      //
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


#ifndef FERMIONONLATTICEWITHSPINSZSYMMETRYREALSPACEAND1DTRANSLATION_H
#define FERMIONONLATTICEWITHSPINSZSYMMETRYREALSPACEAND1DTRANSLATION_H

#include "config.h"
#include "HilbertSpace/FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations.h"

#include <iostream>



class FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation : public FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations
{

  friend class FermionOnSquareLatticeWithSU4SpinMomentumSpace;

 protected:

  // total number of sites
  int NbrSite;
  
  // flag to indicate that the Hilbert space should preserve Sz
  bool SzFlag;

 public:

  // default constructor
  // 
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // momnetum = momentum sector
  // periodicity = periodicity with respect to site numbering 
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation (int nbrFermions, int nbrSite, int momentum, int periodicity, bool minusSzParity, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation(const FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation& fermions);

  // destructor
  //
  ~FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation& operator = (const FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation& fermions);

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
  // maxMomentum = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int maxMomentum);

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int currentSite, long pos);

};


#endif


