////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     set of functions used to managed files related to QHE on torus         //
//                                                                            //
//                        last modification : 12/04/2010                      //
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


#ifndef FQHEONTORUSFILETOOLS_H
#define FQHEONTORUSFILETOOLS_H

#include "config.h"


// try to guess system information from file name
//
// filename = file name
// nbrParticles = reference to the number of particles
// kyMax = reference to the maximum momentum for a single particle
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnTorusFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& kyMax, bool& statistics);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// kyMax = reference to the maximum momentum for a single particle
// ky = reference to the y projection of the angular momentum
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnTorusFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& ky, bool& statistics);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// kyMax = reference to the momentum for a single particle
// kx = reference to the x momentum
// ky = reference to the y momentum
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnTorusFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& kx, int& ky, bool& statistics);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// kyMax = reference to the momentum for a single particle
// kx = reference to the x momentum
// ky = reference to the y momentum
// ratio = reference on the aspect ratio
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnTorusFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& kx, int& ky, double& ratio, bool& statistics);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// kyMax = reference to the momentum for a single particle
// ky = reference to the y momentum
// ratio = reference on the aspect ratio
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnTorusFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& ky, double& ratio, bool& statistics);

// try to guess system information from file name
//
// filename = vector file name
// NbrParticles = reference to the number of particles 
// NbrFluxQuanta (kyMax) = reference to the momentum for a single particle
// Momentum (ky) = reference to the y momentum
// Ratio = reference on the aspect ratio
// Statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnTorusFindSystemInfoFromVectorFileName_SpectralResponse(char* filename, int& NbrParticles, int& NbrFluxQuanta, int& Momentum, double& Ratio, bool& Statistics);

// try to guess system information from file name for system suth an SU(2) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles
// kyMax = reference to the maximum momentum for a single particle
// ky = reference to the y projection of the angular momentum
// sz = reference to twice the z projection of the total spin
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnTorusWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& ky, int& sz, bool& statistics);

// try to guess system information from file name for system suth an SU(2) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles
// kyMax = reference to the maximum momentum for a single particle
// kx = reference to the x momentum
// ky = reference to the y projection of the angular momentum
// sz = reference to twice the z projection of the total spin
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnTorusWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& kx, int& ky, int& sz, bool& statistics);

// try to guess system information from file name for system suth an SU(3) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles
// kyMax = reference to the maximum momentum for a single particle
// ky = reference to the y projection of the angular momentum
// totalTz = reference to twice the total Tz value
// totalY = reference to three time the total Y value
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnTorusWithSU3SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& ky, int& totalTz, int& totalY, bool& statistics);

// try to guess system information from file name for system suth an SU(7) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles
// kyMax = reference to the maximum momentum for a single particle
// ky = reference to the y projection of the angular momentum
// totalSz = reference to twice the total spin value
// totalIz = reference to the total isospin value
// totalPz = reference to the total entanglement value
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnTorusWithSU4SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& ky, int& totalSz, int& totalIz, int& totalPz, bool& statistics);

#endif
