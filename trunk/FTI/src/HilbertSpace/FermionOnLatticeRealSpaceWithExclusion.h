////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of fermions on lattice in real space              //
//                   and an exclusion rule for neighboring sites              //
//                                                                            //
//                        last modification : 17/02/2018                      //
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


#ifndef FERMIONONLATTICEREALSPACEWITHEXCLUSION_H
#define FERMIONONLATTICEREALSPACEWITHEXCLUSION_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeRealSpace.h"

#include <iostream>



class FermionOnLatticeRealSpaceWithExclusion : public FermionOnLatticeRealSpace
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
  FermionOnLatticeRealSpaceWithExclusion ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // excludedSites = array of the linearized coordinates of the excluded sites around each site
  // nbrExcludedSites = number of excluded sites around each site
  // memory = amount of memory granted for precalculations
  FermionOnLatticeRealSpaceWithExclusion (int nbrFermions, int nbrSite, int** excludedSites,
					  int* nbrExcludedSites, unsigned long memory = 10000000);

  // constructor from prebuilt exclusion rules
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // excludedSiteMasks = masks used to detected excluded sites around a given position
  // excludedSiteCenters = masks used to indicate the site around which a given exclusion rule is defined
  // nbrExcludedSiteMasks = number of masks in ExcludedSiteMasks
  // memory = amount of memory granted for precalculations
  FermionOnLatticeRealSpaceWithExclusion (int nbrFermions, int nbrSite, unsigned long* excludedSiteMasks,
					  unsigned long* excludedSiteCenters, int nbrExcludedSiteMasks, unsigned long memory = 10000000);

  // copy constructor (without duplicating data)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeRealSpaceWithExclusion(const FermionOnLatticeRealSpaceWithExclusion& fermions);

  // destructor
  //
  ~FermionOnLatticeRealSpaceWithExclusion ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeRealSpaceWithExclusion& operator = (const FermionOnLatticeRealSpaceWithExclusion& fermions);

  // clone Hilbert space (without duplicating data)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

 protected:

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();

};


#endif


