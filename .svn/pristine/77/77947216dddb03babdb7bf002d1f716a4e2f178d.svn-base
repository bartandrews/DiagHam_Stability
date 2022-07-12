////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of Landau level spectrum on sphere                 //
//                                                                            //
//                        last modification : 14/01/2005                      //
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


#ifndef LANDAUSPECTRUMONSPHERE_H
#define LANDAUSPECTRUMONSPHERE_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class LandauSpectrumOnSphere
{

 protected:

  // number of particles
  int NbrParticles;

  // number of Landau levels filled with composite fermions
  int NbrLandauLevels;

  // twice the value of the momentum in the lowest pseudo-Landau level
  int TwiceS;
  
  // array describing the occupation of the Landau levels
  int** LevelOccupation;

  // garbage flag to avoid duplication of occupation array
  GarbageFlag Flag;

 public:

  // default constructor
  //
  LandauSpectrumOnSphere();

  // constructor from datas
  //
  // nbrLandauLevels = number of Landau levels
  // momentum = twice the value of the momentum in the lowest Landau level
  // description = pointer to the string containing the description
  LandauSpectrumOnSphere(int nbrLandauLevels, int momentum, char* description);

  // copy constructor
  //
  // spectrum = reference on the spectrum to copy
  LandauSpectrumOnSphere(const LandauSpectrumOnSphere& spectrum);

  // destructor
  //
  ~LandauSpectrumOnSphere();

  // assignement
  //
  // spectrum = reference on the spectrum to assign
  // return value = reference on current spectrum
  LandauSpectrumOnSphere& operator = (const LandauSpectrumOnSphere& spectrum);

  // get occupation of a given state in the spectrum
  //
  // level = index of the Landau level
  // momentum = index of the state in the givne Landau level) (0 being the the one with the lowest Lz projection)
  // return value = reference on the inetger containing state occupation
  int& operator () (int level, int momentum);

  // get number of particles that lie in the spectrum
  //
  // return value = number of particles
  int GetNbrParticles();

  // print Landau level spectrum
  //
  // str = reference on the output stream
  // return value = reference on the output stream
  ostream& PrintSpectrum(ostream& str);

  // evaluate number of changes needed to go from one spectrum to another
  //
  // spectrum = reference on the spectrum to compare with
  // return value = number of changes (-1 if an error occurs)
  int EvaluateDistance(LandauSpectrumOnSphere& spectrum);

 protected:

  
  // get occupation information from a formatted string
  //
  // descriptionString = pointer to the string containing the description
  // descriptionArray = reference on the array where description has to be stored
  // return value = number of particles (0 if an error occured)
  int ParseOccupationDescription (char* descriptionString, int**& descriptionArray);


};

// get occupation of a given state in the spectrum
//
// level = index of the Landau level
// momentum = index of the state in the givne Landau level) (0 being the the one with the lowest Lz projection)
// return value = reference on the inetger containing state occupation

inline int& LandauSpectrumOnSphere::operator () (int level, int momentum)
{
  return this->LevelOccupation[level][momentum];
}

// get number of particles that lie in the spectrum
//
// return value = number of particles

inline int LandauSpectrumOnSphere::GetNbrParticles()
{
  return this->NbrParticles;
}

#endif
