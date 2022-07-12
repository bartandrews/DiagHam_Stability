////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2002-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//            class of hilbert space of one particle in magnetic field        //
//                                                                            //
//                        last modification : 04/22/2004                      //
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


#ifndef VERTICALPERIODICPARTICLEINMAGNETICFIELD_H
#define VERTICALPERIODICPARTICLEINMAGNETICFIELD_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"


class VerticalPeriodicParticleInMagneticField : public AbstractHilbertSpace
{

 protected:

  // quantum number of kinetic momentum in Z direction
  int NumberM;

  // number of Landau states in plane
  int NbrStateR;

  // wave function basis dimension in the z direction
  int NbrStateZ;  
  
  // lower limit of Z impulsion
  int LowerImpulsionZ;

 public:

  // constructor
  //
  // nbrStateR = number of Landau states in plane 
  // nbrStateZ = wave function basis dimension in the z direction
  // lowerImpulsionZ = LowerImpulsionZ
  // numberM = quantum number of kinetic momentum in Z direction
  VerticalPeriodicParticleInMagneticField(int numberM, int nbrStateR, int nbrStateZ, int lowerImpulsionZ);

  // copy constructor
  //
  // space = reference on Hilbert space to copy
  VerticalPeriodicParticleInMagneticField(const VerticalPeriodicParticleInMagneticField& space);

  // destructor
  //
  ~VerticalPeriodicParticleInMagneticField();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // assignement
  //
  // space = reference on Hilbert space to assign
  // return value = reference on current Hilbert space
  VerticalPeriodicParticleInMagneticField& operator = (const VerticalPeriodicParticleInMagneticField& space);

  // get wave function basis dimension in plane
  //
  // return value = wave function basis dimension in the x direction
  int GetNbrStateR();

  // get wave function basis dimension in the z direction
  //
  // return value = wave function basis dimension in the z direction
  int GetNbrStateZ();

  // get lower limit of Z impulsion
  // 
  // return value = lower limit of Z impulsion
  int GetLowerImpulsionZ();

  // get quantum number of kinetic momentum
  //
  // return = the kinetic momentum quantum number value
  int GetQuantumNumberM();

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

// get wave function basis dimension in plane
//
// return value = wave function basis dimension in the x direction

inline int VerticalPeriodicParticleInMagneticField::GetNbrStateR()
{
  return this->NbrStateR;
}


// get wave function basis dimension in the z direction
//
// return value = wave function basis dimension in the z direction

inline int VerticalPeriodicParticleInMagneticField::GetNbrStateZ()
{
  return this->NbrStateZ;
}

// get lower limit of Z impulsion
// 
// return value = lower limit of Z impulsion

inline int VerticalPeriodicParticleInMagneticField::GetLowerImpulsionZ()
{
  return this->LowerImpulsionZ;
}

// get quantum number of kinetic momentum
//
//return = the kinetic momentum quantum number value

inline int VerticalPeriodicParticleInMagneticField::GetQuantumNumberM()
{
  return this->NumberM;
}

#endif
