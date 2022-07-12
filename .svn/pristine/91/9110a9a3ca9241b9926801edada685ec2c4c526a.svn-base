////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         set of functions used to handle sphere pseudopotential files       //
//                                                                            //
//                        last modification : 26/02/2009                      //
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


#ifndef FQHESPHEREPSEUDOPOTENTIALTOOLS_H
#define FQHESPHEREPSEUDOPOTENTIALTOOLS_H

#include "config.h"



// get pseudopototentials for spinless particles on sphere from file
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMax = reference on twice the maximum Lz value
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
// onebodyPotential =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state), null pointer if none
// return value = true if no error occured
bool FQHESphereGetPseudopotentials (char* fileName, int lzMax, double* pseudoPotentials, double*& oneBodyPotential);

// get the one-body potential for spinless particles on sphere from file
// 
// fileName = name of the file that contains the one-body potential description
// lzMax = reference on twice the maximum Lz value
// onebodyPotential =  reference on the one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state)
// return value = true if no error occured
bool FQHESphereGetOneBodyPotentials (char* fileName, int& lzMax, double*& oneBodyPotential);

// get pseudopototentials for particles on sphere with SU(2) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMax = reference on twice the maximum Lz value
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// return value = true if no error occured
bool FQHESphereSU2GetPseudopotentials (char* fileName, int lzMax, double** pseudoPotentials,
				       double*& oneBodyPseudopotentialUpUp, double*& oneBodyPseudopotentialDownDown);

// get pseudopototentials for particles on sphere with SU(2) spin from file, including a tunneling term
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMax = reference on twice the maximum Lz value
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyPotentialUpDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// return value = true if no error occured
bool FQHESphereSU2GetPseudopotentials (char* fileName, int lzMax, double** pseudoPotentials,
				       double*& oneBodyPseudopotentialUpUp, double*& oneBodyPseudopotentialDownDown, double*& oneBodyPseudopotentialUpDown);

// get pseudopototentials for particles on sphere with SU(2) spin from file, including a tunneling term and pairing
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMax = reference on twice the maximum Lz value
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyPotentialUpDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyPotentialUpDown = reference on the one body potential for the pairing term (sorted from component on the lowest Lz state to component on the highest Lz state), null pointer if none
// return value = true if no error occured
bool FQHESphereSU2GetPseudopotentialsWithPairing (char* fileName, int lzMax, double** pseudoPotentials,
						  double*& oneBodyPseudopotentialUpUp, double*& oneBodyPseudopotentialDownDown, double*& oneBodyPseudopotentialUpDown,
						  double*& oneBodyPseudopotentialPairing);

// get pseudopototentials for particles on sphere with SU(2) spin from file including all possible interaction terms
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMaxUp = reference on twice the maximum Lz value for particle with spin up
// lzMaxDown = reference on twice the maximum Lz value for particle with spin down
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as upup-upup, updown-upup, downdown-upup, 
//                                                                     upup-updown, updown-updown, downdown-updown,
//                                                                     upup-downdown, updown-downdown, downdown-downdown)
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialUpDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// return value = true if no error occured
//bool FQHESphereFullSU2GetPseudopotentials (char* fileName, int lzMaxUp, lzMaxDown, double** pseudoPotentials,
//					   double* oneBodyPseudopotentialUpUp, double* oneBodyPseudopotentialDownUp,
//					   double* oneBodyPseudopotentialUpDown, double* oneBodyPseudopotentialDownDown);


// get pseudopototentials for particles on sphere with two landau levels
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMax = reference on twice the maximum Lz value of the LLL
// pseudoPotentials = array or arrays of pseudo-potentials. 9 in total which go in ascending order of index p = 0-8 where label p = l*3 + r where l and r label the ll indices on the left and right of the 
//                   interaction respectively and take values 0:up-up, 1: down-down, 2: up-down.
// return value = true if no error occured
bool FQHESphereTwoLandauLevelGetPseudopotentials (char* fileName, int lzMax, double** pseudoPotentials);

// get pseudopototentials for particles on sphere with SU(4) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMax = twice the maximum Lz value of the LLL
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                    first index refered to the spin sector (sorted as 1-1, 1-2, 1-3, 2-2, 2-3, 3-3)
// return value = true if no error occured
bool FQHESphereSU3GetPseudopotentials (char* fileName, int lzMax, double** pseudoPotentials);

// get pseudopototentials for particles on sphere with SU(4) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMax = twice the maximum Lz value of the LLL
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as pplus-upplus, upplus-upminus, upplus-downplus, upplus-downminus, 
//                                                                     upminus-upminus, upminus-downplus, upminus-downminus, downplus-downplus, downplus-downminus, downminus-downminus)
// return value = true if no error occured
bool FQHESphereSU4GetPseudopotentials (char* fileName, int lzMax, double** pseudoPotentials);

#endif

