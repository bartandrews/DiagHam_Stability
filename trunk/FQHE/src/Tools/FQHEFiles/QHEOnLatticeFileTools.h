////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    set of functions used to managed files related to QHE on lattice         //
//                                                                            //
//                        last modification : 06/06/2006                      //
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


#ifndef QHEONLATTICEFILETOOLS_H
#define QHEONLATTICEFILETOOLS_H

#include "config.h"

// see if we have a generalised lattice, or standard square lattice
// filename = filename to be tested
// return = true if general lattice
bool FQHEOnLatticeHaveGeneralLattice(char* filename);


bool FQHEOnLatticeFindSystemInfoFromFileName1(char* filename, int& nbrParticles, int& lx, int& ly, int &flux, bool& statistics, bool& hardcore);

// try to guess system information from file name
//
// filename = file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lx, ly = indication of lattice geometry
// interaction = strength of interaction parameter
// flux = number of flux quanta
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// hardcore = returns true if hardcore bosons encountered
// return value = true if no error occured
bool FQHEOnLatticeFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& lx, int& ly, double &interaction, int &flux, bool& statistics, bool &hardcore);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lx, ly = indication of lattice geometry
// interaction = strength of interaction parameter
// flux = number of flux quanta
// nbrState = number of eigenstate
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// hardcore = returns true if hardcore bosons encountered
// return value = true if no error occured
bool FQHEOnLatticeFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lx, int& ly, double &interaction, int &flux, int &nbrState, bool& statistics, bool &hardcore);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lx, ly = indication of lattice geometry
// interaction = strength of interaction parameter
// kyMomentum = momentum in y-direction when using symmetries
// flux = number of flux quanta
// nbrState = number of eigenstate
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// hardcore = returns true if hardcore bosons encountered
// return value = true if no error occured
bool FQHEOnLatticeFindSystemInfoWithKyFromVectorFileName(char* filename, int& nbrParticles, int& lx, int& ly, int &kyMomentum, double &interaction, int &flux, int &nbrState, bool& statistics, bool &hardcore);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// interaction = strength of interaction parameter
// flux = number of flux quanta
// nbrState = number of eigenstate
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// hardcore = returns true if hardcore bosons encountered
// return value = true if no error occured
bool FQHEOnLatticeFindSystemInfoFromGeneralVectorFileName(char* filename, int& nbrParticles, double &interaction, int &flux, int &nbrState, bool& statistics, bool &hardcore, bool &haveContFlux, double &contFlux);



#endif
