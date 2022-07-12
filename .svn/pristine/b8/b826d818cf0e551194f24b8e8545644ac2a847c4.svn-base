////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//            class of fermions on lattice with spin  and Gutzwiller          //
//   projection in real space with translation invariance in one direction    //
//                                                                            //
//                       class author: Nicolas Regnault                       //
//                                                                            //
//                        last modification : 13/07/2014                      //
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


#ifndef FERMIONONLATTICEWITHSPINANDGUTZWILLERPROJECTIONREALSPACEAND1DTRANSLATION_H
#define FERMIONONLATTICEWITHSPINANDGUTZWILLERPROJECTIONREALSPACEAND1DTRANSLATION_H


#include "config.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd1DTranslation.h"

#include <iostream>



class FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation : public FermionOnLatticeWithSpinRealSpaceAnd1DTranslation
{

  friend class FermionOnSquareLatticeWithSU4SpinMomentumSpace;

 protected:


 public:

  // default constructor
  // 
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // momentum = momentum sector
  // periodicity = periodicity with respect to site numbering 
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation (int nbrFermions, int nbrSite, int momentum, int periodicity, unsigned long memory = 10000000);
  
  // basic constructor when Sz is conserved
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // momentum = momentum sector
  // periodicity = periodicity with respect to site numbering 
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation (int nbrFermions, int totalSpin, int nbrSite, int momentum, int periodicity, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation(const FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation& fermions);

  // destructor
  //
  ~FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation& operator = (const FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions);
  
  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // nbrSpinUp = number of spin up
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int nbrSpinUp);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int currentSite, long pos);

  // generate all states corresponding to the constraints, knowing the number of holes
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // nbrHoles = number of unoccupied sites
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStatesWithHoleCounting(int nbrFermions, int currentSite, int nbrHoles, long pos);
  
  
  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index
  // pos = position in StateDescription array where to store states	
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, long pos);
  
  // generate all states corresponding to the constraints with Sz conserved
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index
  // nbrHoles = number of unoccupied sites
  // nbrSpinUp = number of fermions with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored

  virtual long RawGenerateStatesWithHoleCounting(int nbrFermions, int currentSite, int nbrHoles, int nbrSpinUp, long pos);

};


#endif


