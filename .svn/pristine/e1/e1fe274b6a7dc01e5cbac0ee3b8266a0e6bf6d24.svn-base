////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//       in real space with translation invariance in two directions          //
//        with a constraint on the mininum number of on-site singlet          //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                        last modification : 14/07/2016                      //
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


#ifndef FERMIONONLATTICEWITHSPINREALSPACEAND2DTRANSLATIONMNINNBRSINGLETSLONG_H
#define FERMIONONLATTICEWITHSPINREALSPACEAND2DTRANSLATIONMNINNBRSINGLETSLONG_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong.h"

#include <iostream>



class FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong : public FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong
{


 protected:

  // minimum number of on-site singlets
  int MinNbrSinglets;

 public:

  // default constructor
  // 
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // minNbrSinglets = minimum number of on-site singlets
  // nbrSite = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong (int nbrFermions, int minNbrSinglets, int nbrSite, int xMomentum, int maxXMomentum,
								   int yMomentum, int maxYMomentum, unsigned long memory = 10000000);
  
  // basic constructor when Sz is preserved
  // 
  // nbrFermions = number of fermions
  // minNbrSinglets = minimum number of on-site singlets
  // totalSpin = twice the value of Sz
  // nbrSite = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong (int nbrFermions, int minNbrSinglets, int totalSpin, int nbrSite, int xMomentum, int maxXMomentum,
								   int yMomentum, int maxYMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong(const FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong& fermions);

  // destructor
  //
  ~FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong& operator = (const FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentSite = current site index
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentSite, int nbrSinglets);
  
  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentSite = current site index
  // nbrSpinUp = number of fermions with spin up
  // nbrSinglets = number of on-site singlets
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentSite, int nbrSpinUp, int nbrSinglets);

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // nbrSinglets = number of on-site singlets
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int currentSite, int nbrSinglets, long pos);
  
  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // nbrSpinUp = number of fermions with spin up
  // nbrSinglets = number of on-site singlets
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, int nbrSinglets, long pos);


};

#endif


