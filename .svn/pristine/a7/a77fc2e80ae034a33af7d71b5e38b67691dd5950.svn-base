////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  set of functions used to manage file names                //
//                   related to the FQHE on cylinder geometry                 //
//                                                                            //
//                        last modification : 08/06/2016                      //
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


#ifndef FQHEONCYLINDERFILETOOLS_H
#define FQHEONCYLINDERFILETOOLS_H

#include "config.h"


// try to guess system information from file name
//
// filename = file name
// nbrParticles = reference to the number of particles 
// lzMax = reference to twice the maximum momentum for a single particle
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons)
// ratio = reference to the cylinder aspect ratio
// perimeter = reference to the cylinder perimeter
// return value = true if no error occured
bool FQHEOnCylinderFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& lzMax, bool& statistics,
					      double& ratio, double& perimeter);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles
// lzMax = reference to twice the maximum momentum for a single particle
// ky = reference to the momentum along the cylinder perimeter
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// ratio = reference to the cylinder aspect ratio
// perimeter = reference to the cylinder perimeter
// return value = true if no error occured
bool FQHEOnCylinderFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& ky, bool& statistics,
						    double& ratio, double& perimeter);


// try to guess system information from file name for system with an SU(2) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles
// lzMax = reference to twice the maximum momentum for a single particle
// ky = reference to the momentum along the cylinder perimeter
// sz = reference to twice the z projection of the total spin
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// ratio = reference to the cylinder aspect ratio
// perimeter = reference to the cylinder perimeter
// return value = true if no error occured
bool FQHEOnCylinderWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& ky, int& sz, bool& statistics,
							    double& ratio, double& perimeter);

#endif
