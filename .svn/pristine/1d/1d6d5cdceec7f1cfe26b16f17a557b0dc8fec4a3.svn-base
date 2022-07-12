////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             set of functions used to squeezed (aka Haldane) basis          //
//                                                                            //
//                        last modification : 24/02/2009                      //
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


#ifndef FQHESQUEEZEDBASISTOOLS_H
#define FQHESQUEEZEDBASISTOOLS_H

#include "config.h"


// get the root parition from a file
// 
// rootFileName = name of the file that contains the root description
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the maximum Lz value
// referenceState = array where the root partition description will be stored
// return value = true if no error occured
bool FQHEGetRootPartition (char* rootFileName, int& nbrParticles, int& lzMax, int*& referenceState);

// get the root partition from a file in the SU2 case
// 
// rootFileName = name of the file that contains the root description
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the maximum Lz value
// referenceStates = array where the root partition descriptions will be stored
// nbrReferenceStates = number of root partitions that have been extracted
// return value = true if no error occured
bool FQHEGetRootPartitionSU2 (char* rootFileName, int& nbrParticles, int& lzMax, 
			      int**& referenceStates, int& nbrReferenceStates);
			      
// get the root partition from a file in the SU2 case
// 
// rootFileName = name of the file that contains the root description
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the maximum Lz value
// referenceStates = array where the root partition descriptions will be stored
// nbrReferenceStates = number of root partitions that have been extracted
// texturelessFlag = flag to indicate whether or not to consider spin texture when performing squeezing
// return value = true if no error occured
bool FQHEGetRootPartitionSU2 (char* rootFileName, int& nbrParticles, int& lzMax, 
			      int**& referenceStates, int& nbrReferenceStates, bool &texturelessFlag);			      

#endif

