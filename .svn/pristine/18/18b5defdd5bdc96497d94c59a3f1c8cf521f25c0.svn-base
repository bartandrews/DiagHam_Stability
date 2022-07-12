////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice in real space           //
//                        with nearest neighbor exclusion                     //
//                                                                            //
//                        last modification : 11/02/2018                      //
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


#ifndef FERMIONONSQUARELATTICEREALSPACENNEXCLUSION_H
#define FERMIONONSQUARELATTICEREALSPACENNEXCLUSION_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeRealSpace.h"

#include <iostream>



class FermionOnSquareLatticeRealSpaceNNExclusion : public FermionOnLatticeRealSpace
{

 friend class FermionOnLatticeRealSpaceAnd2DTranslation;
 friend class FermionOnLatticeWithSpinRealSpaceAnd2DTranslation;
 friend class BosonOnLatticeRealSpaceAnd2DTranslation;

 protected:

 // number of sites in the x direction
 int NbrSitesX;
 // number of sites in the y direction
 int NbrSitesY;  

 public:

  // default constructor
  // 
  FermionOnSquareLatticeRealSpaceNNExclusion ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSitesX = number of sites in the x direction
  // nbrSitesY = number of sites in the y direction
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeRealSpaceNNExclusion (int nbrFermions, int nbrSitesX, int nbrSitesY, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeRealSpaceNNExclusion(const FermionOnSquareLatticeRealSpaceNNExclusion& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeRealSpaceNNExclusion ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeRealSpaceNNExclusion& operator = (const FermionOnSquareLatticeRealSpaceNNExclusion& fermions);

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

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentSite = current site linearized index
  // currentSiteY = y coordinate of the current site
  // previousXConfiguration = configuraton at the previous x coordinate
  // currentXConfiguration = configuraton at the current x coordinate
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentSite, int currentSiteY, unsigned long previousXConfiguration, unsigned long currentXConfiguration);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site linearized index
  // currentSiteY = y coordinate of the current site
  // previousXConfiguration = configuraton at the previous x coordinate
  // currentXConfiguration = configuraton at the current x coordinate
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentSite, int currentSiteY, unsigned long previousXConfiguration, unsigned long currentXConfiguration, long pos);

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  int FindStateIndex(unsigned long stateDescription, int lzmax);

};


#endif


