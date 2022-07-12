////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of bosons on the 4D space  torus x sphere              //
//            with magnetic translations and for system size such that        //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 22/02/2017                      //
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


#ifndef BOSONONT2XS2WITHMAGNETICTRANSLATIONSSHORT_H
#define BOSONONT2XS2WITHMAGNETICTRANSLATIONSSHORT_H


#include "config.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"
#include "Matrix/HermitianMatrix.h"

#include <iostream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;


class BosonOnT2xS2WithMagneticTranslationsShort :  public BosonOnTorusWithMagneticTranslationsShort
{

 protected:

  // number of flux quanta piercing the torus
  int NbrFluxQuantumTorus;
  // number of flux quanta piercing the sphere
  int NbrFluxQuantumSphere;
  // projection of the total angular momentum along the z axis for the sphere
  int TotalLz;
  // projection of the total angular momentum along the z axis, using the disk convention
  int ShiftedTotalLz;
  // number of orbitals on the sphere
  int NbrLzValues;

 public:

  // default constructor
  // 
  BosonOnT2xS2WithMagneticTranslationsShort ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrFluxQuantumTorus = number of flux quanta piercing the torus
  // kxMomentum = momentum in the x direction for the torus
  // kyMomentum = momentum in the y direction for the torus
  // nbrFluxQuantumSphere = number of flux quanta piercing the sphere
  // totalLz = projection of the total angular momentum along the z axis for the sphere
  BosonOnT2xS2WithMagneticTranslationsShort (int nbrBosons, int nbrFluxQuantumTorus, int kxMomentum, int kyMomentum,
					     int nbrFluxQuantumSphere, int totalLz);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnT2xS2WithMagneticTranslationsShort(const BosonOnT2xS2WithMagneticTranslationsShort& bosons);

  // destructor
  //
  ~BosonOnT2xS2WithMagneticTranslationsShort ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnT2xS2WithMagneticTranslationsShort& operator = (const BosonOnT2xS2WithMagneticTranslationsShort& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // print a given State using the monomial notation
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintStateMonomial (ostream& Str, long state);
  
  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentKyMax = torus momentum maximum value for bosons that are still to be placed
  // currentLz = current maximum angular momentum for bosons that are still to be placed
  // currentKy = current value of the momentum along y for the torus 
  // currentTotalLz = current total angular momentum along the z axis for the sphere  
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKyMax, int currentLz, int currentKy, int currentTotalLz);

  // generate all states with both the kx and ky constraint
  // 
  // return value = new dimension of the Hilbert space
  virtual long GenerateStates();

  // generate all states corresponding to the ky/Lz constraints without taking care of the kx constraint
  // 
  // nbrBosons = number of bosons
  // currentKyMax = torus momentum maximum value for bosons that are still to be placed
  // currentLz = current maximum angular momentum for bosons that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentKy = current value of the momentum along y for the torus 
  // currentTotalLz = current total angular momentum along the z axis for the sphere  
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrBosons, int currentKyMax, int currentLz, long pos, int currentKy, int currentTotalLz);

  // get a linearized index from the two momenta
  //
  // ky = momentum along the y direction for the torus
  // lz = z projection of the angular momentum for the sphere
  // return value = linearized index 
  virtual int GetLinearizedIndex(int ky, int lz);

  // get the two momenta associated to a given linearized index
  //
  // index = linearized index 
  // ky = reference on the momentum along the y direction for the torus
  // lz = reference on the z projection of the angular momentum for the sphere
  virtual void GetLinearizedIndex(int index, int& ky, int& lz);

};

// get a linearized index from the two momenta
//
// ky = momentum along the y direction for the torus
// lz = z projection of the angular momentum for the sphere
// return value = linearized index 

inline int BosonOnT2xS2WithMagneticTranslationsShort::GetLinearizedIndex(int ky, int lz)
{
  return ((ky * this->NbrLzValues) + lz);
}

// get the two momenta associated to a given linearized index
//
// index = linearized index 
// ky = reference on the momentum along the y direction for the torus
// lz = reference on the z projection of the angular momentum for the sphere

inline void BosonOnT2xS2WithMagneticTranslationsShort::GetLinearizedIndex(int index, int& ky, int& lz)
{
  lz = index % this->NbrLzValues;
  ky = index / this->NbrLzValues;
}


#endif


