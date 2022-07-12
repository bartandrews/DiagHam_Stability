////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice in real space           //
//                with C4 symmetry and nearest neighbor exclusion             //
//                                                                            //
//                        last modification : 25/02/2018                      //
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


#ifndef FERMIONONSQUARELATTICEREALSPACEANDC4SYMMETRYWITHEXCLUSION_H
#define FERMIONONSQUARELATTICEREALSPACEANDC4SYMMETRYWITHEXCLUSION_H

#include "config.h"
#include "HilbertSpace/FermionOnSquareLatticeRealSpaceAndC4Symmetry.h"

#include <iostream>



class FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion : public FermionOnSquareLatticeRealSpaceAndC4Symmetry
{

 protected:

  // masks used to detected excluded sites around a given position
  unsigned long* ExcludedSiteMasks;
  // masks used to indicate the site around which a given exclusion rule is defined
  unsigned long* ExcludedSiteCenters;
  // number of masks in ExcludedSiteMasks
  int NbrExcludedSiteMasks;

 public:

  // default constructor
  // 
  FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSitesX = number of sites in the x (or y) direction
  // symmetrySector = C4 symmetry sector (either 0, 1, 2 or 3)
  // excludedSites = array of the linearized coordinates of the excluded sites around each site
  // nbrExcludedSites = number of excluded sites around each site
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion (int nbrFermions, int nbrSitesX, int symmetrySector, int** excludedSites,
							     int* nbrExcludedSites, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion(const FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion& operator = (const FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

 protected:


  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();


};


#endif


