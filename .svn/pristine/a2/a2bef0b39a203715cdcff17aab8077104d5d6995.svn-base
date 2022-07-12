////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          set of functions used to handle torus pseudopotential files       //
//                                                                            //
//                        last modification : 17/04/2012                      //
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

#ifndef FQHETORUSPSEUDOPOTENTIALTOOLS_H
#define FQHETORUSPSEUDOPOTENTIALTOOLS_H

#include "config.h"



// get pseudopototentials for particles on torus from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrPseudoPotentials = reference on the number of pseudopotentials
// pseudoPotentials = reference on the array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
// return value = true if no error occured
bool FQHETorusGetPseudopotentials (char* fileName, int& nbrPseudoPotentials, double*& pseudoPotentials);

// get pseudopototentials for particles on torus from file with an interaction name
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrPseudoPotentials = reference on the number of pseudopotentials
// pseudoPotentials = reference on the array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
// interactionName = naming convention read from definition, or generated from LL index
// return value = true if no error occured
bool FQHETorusGetPseudopotentials (char* fileName, int& nbrPseudoPotentials, double*& pseudoPotentials, char*& interactionName);


// get pseudopototentials for particles on torus from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrFluxQuanta = number of flux quanta
// nbrPseudoPotentials = reference on the number of pseudopotentials
// pseudoPotentials = reference on the array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
// oneBodyPseudoPotentials  = array with the one-body pseudo-potentials
// return value = true if no error occured
bool FQHETorusGetPseudopotentials (char* fileName, int nbrFluxQuanta, int& nbrPseudoPotentials, double*& pseudoPotentials, double*& oneBodyPotentials);


// get pseudopototentials for particles on torus from file
// 
// fileName = name of the file that contains the pseudopotantial description
// haveCoulomb = flag indicating if Coulomb terms are present
// landauLevel = index of Coulomb Landau-level
// nbrPseudoPotentials = reference on the number of pseudopotentials
// pseudoPotentials = reference on the array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
// interactionName = naming convention read from definition, or generated from LL index
// return value = true if no error occured

bool FQHETorusGetPseudopotentials (char* fileName, bool haveCoulomb, int &landauLevel, int& nbrPseudoPotentials, double*& pseudoPotentials, char*& interactionName);


// get pseudopototentials for particles on torus with SU(2) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// return value = true if no error occured
bool FQHETorusSU2GetPseudopotentials (char* fileName, int* nbrPseudoPotentials, double** pseudoPotentials);

// get pseudopototentials for particles on torus with SU(2) spin from file, including the one-body pseudopotentials if any
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrFluxQuanta = number of flux quanta
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// oneBodyPseudoPotentials  = array with the one body pseudo-potentials (sorted from the ky=0 to the ky=nphi-1)
//                            first index refered to the spin sector (sorted as up-up, down-down, up-down)
// return value = true if no error occured
bool FQHETorusSU2GetPseudopotentials (char* fileName, int nbrFluxQuanta, int* nbrPseudoPotentials, double** pseudoPotentials, double** oneBodyPseudoPotentials);

// get the one-body pseudopototentials for particles on torus with SU(2) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrFluxQuanta = number of flux quanta
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// return value = true if no error occured
bool FQHETorusSU2GetOneBodyPseudopotentials (char* fileName, int nbrFluxQuanta, double*& oneBodyPotentialUpUp, double*& oneBodyPotentialDownDown, double*& oneBodyPotentialUpDown);

// get pseudopototentials for particles on torus with SU(3) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as 11, 12, 13, 22, 23, 33)
// return value = true if no error occured
bool FQHETorusSU3GetPseudopotentials (char* fileName, int* nbrPseudoPotentials, double** pseudoPotentials);

// get pseudopototentials for particles on torus with SU(3) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrFluxQuanta = number of flux quanta
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as 11, 12, 13, 22, 23, 33)
// oneBodyPseudoPotentials  = array with the one body pseudo-potentials (sorted from the ky=0 to the ky=nphi-1)
//                            first index refered to the spin sector (sorted as 11, 12, 13, 22, 23, 33)
// return value = true if no error occured
bool FQHETorusSU3GetPseudopotentials (char* fileName, int nbrFluxQuanta, int* nbrPseudoPotentials, double** pseudoPotentials, double** oneBodyPseudoPotentials);

// get pseudopototentials for particles on torus with SU(3) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrFluxQuanta = number of flux quanta
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// oneBodyPseudoPotentials  = array with the one body pseudo-potentials (sorted from the ky=0 to the ky=nphi-1)
//                            first index refered to the spin sector (sorted as 11, 12, 13, 22, 23, 33)
// return value = true if no error occured
bool FQHETorusSU3GetOneBodyPseudopotentials (char* fileName, int nbrFluxQuanta, double** oneBodyPseudoPotentials);

// get pseudopototentials for particles on torus with SU(4) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as pplus-upplus, upplus-upminus, upplus-downplus, upplus-downminus, 
//                                                                     upminus-upminus, upminus-downplus, upminus-downminus, downplus-downplus, downplus-downminus, downminus-downminus)
// return value = true if no error occured
bool FQHETorusSU4GetPseudopotentials (char* fileName, int* nbrPseudoPotentials, double** pseudoPotentials);


#endif

