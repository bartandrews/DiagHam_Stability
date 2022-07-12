////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                     class author: Cecile Repellin                          //
//                                                                            //
//                                                                            //
//              class of fermions on a honeycomb lattice in real space        //
//                        with plaquette cluster exclusion                    //
//                                                                            //
//                        last modification : 01/04/2018                      //
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


#ifndef FERMIONONHONEYCOMBLATTICEWITHSPINREALSPACEPLAQUETTEEXCLUSION_H
#define FERMIONONHONEYCOMBLATTICEWITHSPINREALSPACEPLAQUETTEEXCLUSION_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"

#include <iostream>



class FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion : public FermionOnLatticeWithSpinRealSpace
{ 
 friend class FermionOnLatticeRealSpaceAnd2DTranslation;
 friend class FermionOnLatticeWithSpinRealSpaceAnd2DTranslation;
 friend class BosonOnLatticeRealSpaceAnd2DTranslation;

 protected:

 // number of sites in the x direction
 int NbrSitesX;
 // number of sites in the y direction
 int NbrSitesY;  
 
 // array listing all the site indices for each hexagon
 int** ListIndicesPerPlaquette;
 // array listing the indices larger than a reference index in the same plaquette
 int*** LargerIndicesInPlaquette;
 // array listing the number of indices in a plaquette larger than the reference index
 int** NbrLargerIndicesInPlaquette;

 public:

  // default constructor
  // 
  FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSitesX = number of sites in the x direction
  // nbrSitesY = number of sites in the y direction
  // memory = amount of memory granted for precalculations
  FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion (int nbrFermions, int nbrSitesX, int nbrSitesY, unsigned long memory = 10000000);
  
  // constructor for conserved Sz
  // 
  // nbrFermions = number of fermions
  // totalSpin = twice the total spin value
  // nbrSitesX = number of sites in the x direction
  // nbrSitesY = number of sites in the y direction
  // memory = amount of memory granted for precalculations
  FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion (int nbrFermions, int totalSpin, int nbrSitesX, int nbrSitesY, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion(const FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion& fermions);

  // destructor
  //
  ~FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion& operator = (const FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

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
  // currentSite = current site linearized index
  // currentSiteY = y coordinate of the current site
  // currentConfiguration = configuraton at the current coordinate
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentSite, unsigned long currentConfiguration);
  
  // evaluate Hilbert space dimension with a fixed number of fermions with spin up
  //
  // nbrFermions = number of fermions
  // nbrSpinUp = number of fermions with spin up
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int nbrSpinUp, int currentSite, unsigned long currentConfiguration);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site linearized index
  // currentConfiguration = configuraton at the current state
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentSite, unsigned long currentConfiguration, long pos);
  
  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // nbrSpinUp = number of fermions with spin up
  // currentSite = current site linearized index
  // currentConfiguration = current configuraton
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int nbrSpinUp, int currentSite, unsigned long currentConfiguration, long pos);

  // find linearized index associated with lattice coordinates
  //
  // siteAlpha = A or B site of honeycomb lattice
  // siteX = coordinate along X axis
  // siteY = coordinate along Y axis
  // return value = linearized index
  int FindSiteIndex(int siteAlpha, int siteX, int siteY);
  
  // find lattice coordinates associated with linearized index
  //
  // index = linearized index of given site
  // siteAlpha = reference to A or B site of honeycomb lattice
  // siteX = reference to coordinate along X axis
  // siteY = reference to coordinate along Y axis
  void FindSiteCoordinates(int index, int& siteAlpha, int& siteX, int& siteY);
  
  // compute the charge on all three hexagons surrounding one site and find the largest one
  //
  // siteIndex = site index
  // configuraton =hilbert space configuraton
  // return value = maximum charge
  virtual int FindMaximumChargeSurroundingPlaquettes (int siteIndex, unsigned long configuraton);
  
  // initialize all of the arrays that will be used to implement hexagon exclusion conditions
  virtual void InitializeHexagonArrays();

};


// find linearized index associated with lattice coordinates
//
// siteAlpha = A or B site of honeycomb lattice
// siteX = coordinate along X axis
// siteY = coordinate along Y axis
// return value = linearized index
//
inline int FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::FindSiteIndex(int siteAlpha, int siteX, int siteY)
{
  if (siteX < 0)
    siteX += this->NbrSitesX;
  if (siteX >= this->NbrSitesX)
    siteX -= this->NbrSitesX;
  
  if (siteY < 0)
    siteY += this->NbrSitesY;
  if (siteY >= this->NbrSitesY)
    siteY -= this->NbrSitesY;
  
  int TmpIndex = siteAlpha + 2 * (this->NbrSitesY * siteX + siteY);
  
  return TmpIndex;
}


// find lattice coordinates associated with linearized index
//
// index = linearized index of given site
// siteAlpha = reference to A or B site of honeycomb lattice
// siteX = reference to coordinate along X axis
// siteY = reference to coordinate along Y axis
inline void FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::FindSiteCoordinates(int index, int& siteAlpha, int& siteX, int& siteY)
{
  siteAlpha = index & 0x1ul;
  index = index >> 0x1ul;
  siteY = index % (this->NbrSitesY);
  siteX = (index - siteY) / this->NbrSitesY;
}


#endif


