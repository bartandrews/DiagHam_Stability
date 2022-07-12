////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2002-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//            class of hilbert space of one particle in magnetic field        //
//              (res
//                               and a quatum well       //
//                                                                            //
//                        last modification : 10/13/2005                      //
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


#ifndef PARTICLEINWELLMAGNETICFIELDTWOSUBBANDS_H
#define PARTICLEINWELLMAGNETICFIELDTWOSUBBANDS_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"


class ParticleInWellMagneticFieldTwoSubbands : public AbstractHilbertSpace
{

 protected:

  // Landau index of the first subband
  int LandauSubbandIndex1;
  // index of the second subband
  int LandauSubbandIndex2;

  // z confinement index of the first subband
  int ZConfinementSubbandIndex1;
  // z confinement index of the second subband
  int ZConfinementSubbandIndex2;
  
  // number of states per Landau level
  int LandauDegeneracy;

  // quantum well length (in Angstroem)
  double QuantumWellLength;
  // magnetic length (in Angstroem)
  
  

 public:

  // constructor
  //
  // landauDegeneracy = number of states per Landau level
  ParticleInWellMagneticFieldTwoSubbands(int landauDegeneracy);

  // copy constructor
  //
  // space = reference on Hilbert space to copy
  ParticleInWellMagneticFieldTwoSubbands(const ParticleInWellMagneticFieldTwoSubbands& space);

  // destructor
  //
  ~ParticleInWellMagneticFieldTwoSubbands();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // assignement
  //
  // space = reference on Hilbert space to assign
  // return value = reference on current Hilbert space
  ParticleInWellMagneticFieldTwoSubbands& operator = (const ParticleInWellMagneticFieldTwoSubbands& space);

  // get the subband index corresponding to a given base state 
  //
  // index = index of the base state
  // return value = subband index (0 or 1)
  int GetLandauSubbandIndex(int index);

  // get the y momentum index corresponding to a given base state 
  //
  // index = index of the base state
  // return value = y momentum index (in 2pi/Ly unit)
  int GetYMomentum (int index);

  // return a list of all possible quantum numbers
  //
  // return value = pointer to corresponding quantum number
  List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q,
					 SubspaceSpaceConverter& converter);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream
  ostream& PrintState (ostream& Str, int state);

};


// get the subband index corresponding to a given base state 
//
// index = index of the base state
// return value = subband index (0 or 1)

inline int ParticleInWellMagneticFieldTwoSubbands::GetLandauSubbandIndex(int index)
{
  return (index & 1);
}

// get the y momentum index corresponding to a given base state 
//
// index = index of the base state
// return value = y momentum index (in 2pi/Ly unit)

inline int ParticleInWellMagneticFieldTwoSubbands::GetYMomentum (int index)
{
  return (index >> 1);
}

#endif
