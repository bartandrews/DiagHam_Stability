////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                        class of fermions on lattice                        //
//       in real space with translation invariance in two directions          //
//                   and an exclusion rule for neighboring sites              //
//                                                                            //
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


#ifndef FERMIONONLATTICEREALSPACEAND2DTRANSLATIONWITHEXCLUSION_H
#define FERMIONONLATTICEREALSPACEAND2DTRANSLATIONWITHEXCLUSION_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceWithExclusion.h"

#include <iostream>
#include <bitset>


class  FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion : public FermionOnLatticeRealSpaceAnd2DTranslation
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
  FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // excludedSites = array of the linearized coordinates of the excluded sites (for the unit cell at (0,0), first index being the orbital index)
  // nbrExcludedSites = number of excluded sites per orbital within the unit cell at (0,0)
  // memory = amount of memory granted for precalculations
  FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion (int nbrFermions, int nbrSite, int xMomentum, int maxXMomentum,
							  int yMomentum, int maxYMomentum, int** excludedSites,
							  int* nbrExcludedSites, unsigned long memory = 10000000);

  // constructor from prebuilt exclusion rules
  // 
  // nbrFermions = number of fermions
  // nbrSite = total number of sites 
  // xMomentum = momentum sector in the x direction
  // xTranslation = translation that has to be applied on the site index to connect two sites with a translation in the x direction
  // yMomentum = momentum sector in the y direction
  // yPeriodicity = periodicity in the y direction with respect to site numbering 
  // excludedSiteMasks = masks used to detected excluded sites around a given position
  // excludedSiteCenters = masks used to indicate the site around which a given exclusion rule is defined
  // nbrExcludedSiteMasks = number of masks in ExcludedSiteMasks
  // memory = amount of memory granted for precalculations
  FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion (int nbrFermions, int nbrSite, int xMomentum,  int maxXMomentum,
							  int yMomentum,  int maxYMomentum, unsigned long* excludedSiteMasks,
							  unsigned long* excludedSiteCenters, int nbrExcludedSiteMasks,
							  unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion(const FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion& fermions);

  // destructor
  //
  ~FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion& operator = (const FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = subsystem momentum along the x direction
  // kySector = subsystem momentum along the x direction
  // groundState = reference on the total system ground state  
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // get an Hilbert space without the translation symmetries, but preserving any other properties
  //
  // return value = Hilbert space without translations
  virtual FermionOnLatticeRealSpace* GetHilbertSpaceWithoutTranslations();

 protected:

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();

};

// get an Hilbert space without the translation symmetries, but preserving any other properties
//
// return value = Hilbert space without translations

inline FermionOnLatticeRealSpace* FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion::GetHilbertSpaceWithoutTranslations()
{
  return new FermionOnLatticeRealSpaceWithExclusion(this->NbrFermions, this->NbrSite, 
						    this->ExcludedSiteMasks, this->ExcludedSiteCenters, this->NbrExcludedSiteMasks);
}


#endif


