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
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                        last modification : 11/07/2014                      //
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


#ifndef FERMIONONLATTICEWITHSPINREALSPACEAND1DTRANSLATION_H
#define FERMIONONLATTICEWITHSPINREALSPACEAND1DTRANSLATION_H

#include "config.h"
#include "HilbertSpace/FermionOnTorusWithSpinAndMagneticTranslations.h"

#include <iostream>



class FermionOnLatticeWithSpinRealSpaceAnd1DTranslation : public FermionOnTorusWithSpinAndMagneticTranslations
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
  FermionOnLatticeWithSpinRealSpaceAnd1DTranslation ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // momnetum = momentum sector
  // periodicity = periodicity with respect to site numbering 
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinRealSpaceAnd1DTranslation (int nbrFermions, int nbrSite, int momentum, int periodicity, unsigned long memory = 10000000);

  // basic constructor when Sz is preserved
  // 
  // nbrFermions = number of fermions
  // totalSpin = twice the total spin value
  // nbrSite = number of sites in the x direction
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinRealSpaceAnd1DTranslation (int nbrFermions, int totalSpin, int nbrSite, int momentum, int periodicity, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeWithSpinRealSpaceAnd1DTranslation(const FermionOnLatticeWithSpinRealSpaceAnd1DTranslation& fermions);

  // destructor
  //
  ~FermionOnLatticeWithSpinRealSpaceAnd1DTranslation ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeWithSpinRealSpaceAnd1DTranslation& operator = (const FermionOnLatticeWithSpinRealSpaceAnd1DTranslation& fermions);

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
  virtual long RawGenerateStates(int nbrFermions, int currentSite, long pos);
  
 
};


#endif


