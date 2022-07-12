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
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//                        last modification : 02/04/2018                      //
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


#ifndef FERMIONONHONEYCOMBLATTICEWITHSPINSZSYMMETRYREALSPACEPLAQUETTEEXCLUSIONAND2DTRANSLATION_H
#define FERMIONONHONEYCOMBLATTICEWITHSPINSZSYMMETRYREALSPACEPLAQUETTEEXCLUSIONAND2DTRANSLATION_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation.h"

#include <iostream>


class FermionOnHoneycombLatticeWithSpinSzSymmetryRealSpacePlaquetteExclusionAnd2DTranslation : public FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation
{

 protected:
   
   // array listing all the site indices for each hexagon
 int** ListIndicesPerPlaquette;
 // array listing the indices larger than a reference index in the same plaquette
 int*** LargerIndicesInPlaquette;
 // array listing the number of indices in a plaquette larger than the reference index
 int** NbrLargerIndicesInPlaquette;


 protected:
    
  // target space for operations leaving the Hilbert space
  FermionOnHoneycombLatticeWithSpinSzSymmetryRealSpacePlaquetteExclusionAnd2DTranslation* TargetSpace;

 public:

  // default constructor
  // 
  FermionOnHoneycombLatticeWithSpinSzSymmetryRealSpacePlaquetteExclusionAnd2DTranslation ();
  
  // basic constructor when Sz is preserved
  // 
  // nbrFermions = number of fermions
  // totalSpin = twice the value of Sz
  // nbrSite = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  FermionOnHoneycombLatticeWithSpinSzSymmetryRealSpacePlaquetteExclusionAnd2DTranslation (int nbrFermions, int totalSpin, int nbrSite, bool minusSzParity, int xMomentum, int maxXMomentum,
						     int yMomentum, int maxYMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnHoneycombLatticeWithSpinSzSymmetryRealSpacePlaquetteExclusionAnd2DTranslation(const FermionOnHoneycombLatticeWithSpinSzSymmetryRealSpacePlaquetteExclusionAnd2DTranslation& fermions);

  // destructor
  //
  ~FermionOnHoneycombLatticeWithSpinSzSymmetryRealSpacePlaquetteExclusionAnd2DTranslation ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnHoneycombLatticeWithSpinSzSymmetryRealSpacePlaquetteExclusionAnd2DTranslation& operator = (const FermionOnHoneycombLatticeWithSpinSzSymmetryRealSpacePlaquetteExclusionAnd2DTranslation& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();


 protected:

  // evaluate Hilbert space dimension with a fixed number of fermions with spin up
  //
  // nbrFermions = number of fermions
  // nbrSpinUp = number of fermions with spin up
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int nbrSpinUp, int currentSite, unsigned long currentConfiguration);

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();
  
  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // nbrSpinUp = number of fermions with spin up
  // currentSite = current site linearized index
  // currentConfiguration = current configuraton
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int nbrSpinUp, int currentSite, unsigned long currentConfiguration, long pos);
 
  // compute the charge on all three hexagons surrounding one site and find the largest one
  //
  // siteIndex = site index
  // configuraton =hilbert space configuraton
  // return value = maximum charge
  virtual int FindMaximumChargeSurroundingPlaquettes (int siteIndex, unsigned long configuraton);
  
  // initialize all of the arrays that will be used to implement hexagon exclusion conditions
  virtual void InitializeHexagonArrays();

};

#endif


