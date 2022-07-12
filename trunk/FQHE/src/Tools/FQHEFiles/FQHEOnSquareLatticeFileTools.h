////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             set of functions used to managed files related to FQHE         //
//                          on square lattice disk                            //
//                                                                            //
//                        last modification : 01/03/2011                      //
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


#ifndef FQHEONSQUARELATTICEFILETOOLS_H
#define FQHEONSQUARELATTICEFILETOOLS_H

#include "config.h"


// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FQHEOnSquareLatticeFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, bool& statistics);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumX = reference to the momentum along the y direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& momentumX, int& momentumY, bool& statistics);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumX = reference to the momentum along the y direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FQHEOnSquareLatticeWannierFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& momentumX, int& momentumY, bool& statistics);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumX = reference to the momentum along the y direction
// totalSz = reference to the Sz value
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& momentumX, int& momentumY, int& totalSz, bool& statistics);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumX = reference to the momentum along the y direction
// totalSz = reference to the Sz value
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& totalSz, bool& statistics);

// try to guess system information from file name 
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumX = reference to the momentum along the y direction
// mass = reference to the mass term
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& momentumX, int& momentumY, double& mass, bool& statistics);

// try to guess system information from file name for a cubic lattice
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// nbrSiteZ = reference to the number sites along the z direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FQHEOnCubicLatticeFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& nbrSiteZ, bool& statistics);

// try to guess system information from file name for a cubic lattice
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// nbrSiteZ = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumY = reference to the momentum along the y direction
// momentumZ = reference to the momentum along the z direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FQHEOnCubicLatticeFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& nbrSiteZ, int& momentumX, int& momentumY, int& momentumZ, bool& statistics);

// try to guess system information from file name for a Hofstadter lattice model
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// NbrCellX = number of magnetic unit cells along X
// NbrCellY = number of magnetic unit cells along Y
// Interaction = onsite interaction
// FluxPerCell = number of flux quanta per unit cell
// NbrState = number of the eigenstate
// Statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// Hardcore = flag indicating hard-core bosons
// EmbeddingFlag = flag indicating whether embedding is used
// Axis = character indicating axis of Landau gauge
// GammaX = periodic boundary conditions along X
// GammaY = periodic boundary conditions along Y
// MomentumX = momentum along x-direction
// MomentumY = momentum along y-direction
// UnitCellX = size of magnetic unit cell along x
// UnitCellY = size of magnetic unit cell along y
// return value = true if no error occured
// 
bool FQHEOnSquareLatticeFindSystemInfoFromVectorFileName_Hofstadter(char* filename, int& NbrParticles, int& NbrCellX, int& NbrCellY, double& Interaction, int& FluxPerCell, int& NbrState, bool& Statistics, bool& Hardcore, bool& EmbeddingFlag, char& Axis, double& GammaX, double& GammaY, int& MomentumX, int& MomentumY, int& UnitCellX, int& UnitCellY);


// try to guess system information from file name for a Hofstadter lattice model
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// NbrCellX = number of magnetic unit cells along X
// NbrCellY = number of magnetic unit cells along Y
// Interaction = onsite interaction
// FluxPerCell = number of flux quanta per unit cell
// NbrState = number of the eigenstate
// Statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// Hardcore = flag indicating hard-core bosons
// EmbeddingFlag = flag indicating whether embedding is used
// Axis = character indicating axis of Landau gauge
// GammaX = periodic boundary conditions along X
// GammaY = periodic boundary conditions along Y
// MomentumX = momentum along x-direction
// MomentumY = momentum along y-direction
// UnitCellX = size of magnetic unit cell along x
// UnitCellY = size of magnetic unit cell along y
// enlargeCell = true if unit cell contains 2 * FluxPerCell magnetic fluxes
// muS = amplitude of the symmetry-breaking on-site potential
// nbrBands = number of bands that have to be considered
// return value = true if no error occured
// 
bool FQHEOnSquareLatticeFindSystemInfoFromVectorFileName_Hofstadter(char* filename, int& NbrParticles, int& NbrCellX, int& NbrCellY, double& Interaction, int& FluxPerCell, int& NbrState, bool& Statistics, bool& Hardcore, bool& EmbeddingFlag, char& Axis, double& GammaX, double& GammaY, int& MomentumX, int& MomentumY, int& UnitCellX, int& UnitCellY, bool& enlargeCell, double& muS, int& nbrBands);

#endif
