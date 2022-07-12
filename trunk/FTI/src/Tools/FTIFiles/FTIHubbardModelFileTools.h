////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     set of functions used to managed files related to Hubbard models       //
//                                                                            //
//                        last modification : 19/06/2014                      //
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


#ifndef FTIHUBBARDMODELFILETOOLS_H
#define FTIHUBBARDMODELFILETOOLS_H

#include "config.h"


// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference to the number sites
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// return value = true if no error occured
bool FTIHubbardModelFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& nbrSites, bool& statistics);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference to the number sites
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference to flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, bool& statistics, bool& gutzwiller);

// try to guess system information from file name
//
// filename = vector file name
// nbrSites = reference to the number sites
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// parity = reference to the particle number parity
// return value = true if no error occured
bool FTIHubbardModelFindSystemInfoFromVectorFileName(char* filename, int& nbrSites, bool& statistics, int& parity);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference to the number sites
// szValue = reference on the value of the total spin
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference to flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, bool& statistics, bool& gutzwiller);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// xMomentum = reference on the momentum sector in the x direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWith1DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& xMomentum, int& xPeriodicity,
								      bool& statistics, bool& gutzwiller);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// yPeriodicity = reference on the periodicity in the y direction with respect to site numbering
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& xMomentum, int& yMomentum, 
								      int& xPeriodicity, int& yPeriodicity, bool& statistics, bool& gutzwiller);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// szValue = reference on the value of the total spin
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, int& szSymmetry, bool& statistics, bool& gutzwiller);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// szValue = reference on the value of the total spin
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, bool& statistics, bool& gutzwiller);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// szValue = reference on the value of the total spin
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// yPeriodicity = reference on the periodicity in the y direction with respect to site numbering
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, int& xMomentum, int& yMomentum, 
								      int& xPeriodicity, int& yPeriodicity, bool& statistics, bool& gutzwiller);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// szValue = reference on the value of the total spin
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// yPeriodicity = reference on the periodicity in the y direction with respect to site numbering
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, int& szSymmetry, 
								      int& xMomentum, int& yMomentum, int& xPeriodicity, int& yPeriodicity, bool& statistics, bool& gutzwiller);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// szValue = reference on the value of the total spin
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// minNbrSinglets = minimum number of on-site singlets
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// yPeriodicity = reference on the periodicity in the y direction with respect to site numbering
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, int& szSymmetry, int& minNbrSinglets, 
								      int& xMomentum, int& yMomentum, int& xPeriodicity, int& yPeriodicity, bool& statistics, bool& gutzwiller);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrUnitCellX = reference on the number of unit cell  in the x direction
// nbrUnitCellY = reference on the number of unit cell  in the y direction
// nbrSiteInUnitCellX  = reference on the number of site in each unit cell  in the x direction
// nbrSiteInUnitCellY  = reference on the number of site in each unit cell  in the y direction
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHofstadterdModelFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrUnitCellX, int& nbrUnitCellY, int& nbrSiteInUnitCellX,int& nbrSiteInUnitCellY, bool& statistics, bool& gutzwiller);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// nbrUnitCellX = reference on the number of unit cell  in the x direction
// nbrUnitCellY = reference on the number of unit cell  in the y direction
// nbrSiteInUnitCellX  = reference on the number of site in each unit cell  in the x direction
// nbrSiteInUnitCellY  = reference on the number of site in each unit cell  in the y direction
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHofstadterdModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& xMomentum, int& yMomentum,  int& nbrUnitCellX, int& nbrUnitCellY, int& nbrSiteInUnitCellX, int& nbrSiteInUnitCellY, bool& statistics, bool& gutzwiller);


// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// fluxInUnitCell, = reference on the number of flux in the unit cell
// nbrUnitCellX = reference on the number of unit cell  in the x direction
// nbrUnitCellY = reference on the number of unit cell  in the y direction
// nbrSiteInUnitCellX  = reference on the number of site in each unit cell  in the x direction
// nbrSiteInUnitCellY  = reference on the number of site in each unit cell  in the y direction
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// translationFlag = bool indicating if translations are used
// return value = true if no error occured

bool FTIHofstadterdModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& xMomentum, int& yMomentum, int & fluxInUnitCell, int& nbrUnitCellX, int& nbrUnitCellY, int& nbrSiteInUnitCellX,int& nbrSiteInUnitCellY, bool& statistics, bool& gutzwiller,  bool& translationFlag);


// try to guess system information from file name
//
// filename = vector file name
// szValue = reference on the value of the total spin
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// minNbrSinglets = minimum number of on-site singlets
// return value = true if no error occured
bool FTIHofstadterModelWithSzFindSystemInfoFromVectorFileName(char* filename, int& szValue, int& szSymmetry, int& minNbrSinglets);

// try to guess system information from file name
//
// filename = vector file name
// szValue = reference on the value of the total spin
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// minNbrSinglets = minimum number of on-site singlets
// usingSzSymmetryFlag = bool indicating if the SzSymmetry is used
// usingNbrSingletConstraintFlag = bool indicating if the constraint on the number of singlets is used
// return value = true if no error occured

bool FTIHofstadterModelWithSzFindSystemInfoFromVectorFileName(char* filename,int & szValue, int& szSymmetry, int& minNbrSinglets, bool & usingSzSymmetryFlag, bool & usingNbrSingletConstraintFlag);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// szValue = reference on the value of the total spin
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// yPeriodicity = reference on the periodicity in the y direction with respect to site numbering
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, int& szSymmetry, int& xMomentum, int& yMomentum, int& xPeriodicity, int& yPeriodicity, bool& statistics, bool& gutzwiller, bool& clusterExclusion);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// szValue = reference on the value of the total spin
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// yPeriodicity = reference on the periodicity in the y direction with respect to site numbering
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, int& xMomentum, int& yMomentum, int& xPeriodicity, int& yPeriodicity, bool& statistics, bool& gutzwiller, bool& clusterExclusion);

#endif
