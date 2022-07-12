////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//        class of fermions on sphere with spin and two Landau levels         //
//                                                                            //
//                        last modification : 09/12/2016                      //
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


#ifndef FERMIONONSPHEREWITHSPINTWOLANDAULEVELS_H
#define FERMIONONSPHEREWITHSPINTWOLANDAULEVELS_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"

#include <iostream>

#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::dec;
using std::hex;
using std::ios;
using std::ofstream;

class FermionOnSphere;
class FermionOnSphereWithSpin;
class ParticleOnSphereWithSpin;


class FermionOnSphereWithSpinTwoLandauLevels :  public FermionOnSphereWithSU4Spin
{


  friend class BosonOnSphereWithSU2Spin;


 protected:

  // number of particles with spin up
  int NbrFermionsUp;
  // number of particles with spin down
  int NbrFermionsDown;

  // number of orbitals minus one for the lowest Landau level
  int LzMaxLLL;
  // number of orbitals minus one for the second Landau level
  int LzMax2LL;
  // angular momentum shift for the lowest Landau level
  int LzShiftLLL;
  // angular momentum shift for the second Landau level
  int LzShift2LL;

 public:

  // default constructor
  //
  FermionOnSphereWithSpinTwoLandauLevels();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion in the lowest Landau level
  // totalSpin = twice the total spin value
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinTwoLandauLevels (int nbrFermions, int totalLz, int lzMax, int totalSpin, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpinTwoLandauLevels(const FermionOnSphereWithSpinTwoLandauLevels& fermions);

  // destructor
  //
  ~FermionOnSphereWithSpinTwoLandauLevels ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpinTwoLandauLevels& operator = (const FermionOnSphereWithSpinTwoLandauLevels& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // totalSpin = twice the total spin value
  // return value = Hilbert space dimension
  virtual int EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin);

  // evaluate Hilbert space dimension with shifted values for lzMax and totalLz
  //
  // nbrFermions = number of fermions
  // lzMax = two times momentum maximum value for a fermion plus one 
  // totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
  // totalSpin = number of particles with spin up
  // return value = Hilbert space dimension  
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // totalSpin = number of particles with spin up
  // totalIsospin = number of particles with isospin plus
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int lzMax, int totalLz, int totalSpin, long pos);

  // convert a state to its monomial representation
  //
  // initialState = initial bosonic state in its fermionic representation
  // finalStateLLLUp = reference on the array where the monomial spin up in the lowest Landau level representation has to be stored
  // finalState2LLUp = reference on the array where the monomial spin up in the second Landau level representation has to be stored
  // finalStateLLLDown = reference on the array where the monomial spin down in the lowest Landau level representation has to be stored
  // finalState2LLDown = reference on the array where the monomial spin down in the second Landau level representation has to be stored
  // nbrParticlesLLLUp = reference on the  number of spin up particles in the lowest Landau level
  // nbrParticlesLLLDown = reference on the  number of spin down particles in the lowest Landau level
  virtual void ConvertToMonomial(unsigned long initialState, unsigned long*& finalStateLLLUp, unsigned long*& finalState2LLUp, 
				 unsigned long*& finalStateLLLDown, unsigned long*& finalState2LLDown, 
				 int& nbrParticlesLLLUp,  int& nbrParticlesLLLDown);

  // convert a state from its monomial representation
  //
  // initialStateLLLUp = array where the monomial spin up, lowest Landau level representation is stored
  // initialState2LLUp = array where the monomial spin up, second Landau level representation is stored
  // initialStateLLLDown = array where the monomial spin down, lowest Landau level representation is stored
  // initialState2LLDown = array where the monomial spin down, second Landau level representation is stored
  // nbrParticlesLLLUp = number of spin up particles in the lowest Landau level
  // nbrParticlesLLLDown = number of spin down particles in the lowest Landau level
  // return value = state in its fermionic representation
  virtual unsigned long ConvertFromMonomial(unsigned long* initialStateLLLUp, unsigned long* initialState2LLUp,
					    unsigned long* finalStateLLLDown, unsigned long* finalState2LLDown, 
					    int nbrParticlesLLLUp,  int nbrParticlesLLLDown);

};

// convert a state to its monomial representation
//
// initialState = initial bosonic state in its fermionic representation
// finalStateLLLUp = reference on the array where the monomial spin up in the lowest Landau level representation has to be stored
// finalState2LLUp = reference on the array where the monomial spin up in the second Landau level representation has to be stored
// finalStateLLLDown = reference on the array where the monomial spin down in the lowest Landau level representation has to be stored
// finalState2LLDown = reference on the array where the monomial spin down in the second Landau level representation has to be stored
// nbrParticlesLLLUp = reference on the  number of spin up particles in the lowest Landau level
// nbrParticlesLLLDown = reference on the  number of spin down particles in the lowest Landau level

inline void FermionOnSphereWithSpinTwoLandauLevels::ConvertToMonomial(unsigned long initialState, unsigned long*& finalStateLLLUp, unsigned long*& finalState2LLUp, 
								      unsigned long*& finalStateLLLDown, unsigned long*& finalState2LLDown, 
								      int& nbrParticlesLLLUp,  int& nbrParticlesLLLDown)
{
  nbrParticlesLLLUp = 0;
  nbrParticlesLLLDown = 0;
  int Index2LLUp = 0;
  int Index2LLDown = 0;
  unsigned long Tmp;
  for (int j = this->LzMax; j >= 0; --j)
    {
      Tmp = (initialState >> (j << 2)) & 0xful;
      if ((Tmp & 0x1ul) != 0x0ul)
	{
	  finalStateLLLDown[nbrParticlesLLLDown++] = (unsigned long) j;  
	}
      if ((Tmp & 0x2ul) != 0x0ul)
	{
	  finalState2LLDown[Index2LLDown++] = (unsigned long) j;  
	}
      if ((Tmp & 0x4ul) != 0x0ul)
	{
	  finalStateLLLUp[nbrParticlesLLLUp++] = (unsigned long) j;
	}
      if ((Tmp & 0x8ul) != 0x0ul)
	{
	  finalState2LLUp[Index2LLUp++] = (unsigned long) j;

	}
    }
}

// convert a state from its monomial representation
//
// initialStateLLLUp = array where the monomial spin up, lowest Landau level representation is stored
// initialState2LLUp = array where the monomial spin up, second Landau level representation is stored
// initialStateLLLDown = array where the monomial spin down, lowest Landau level representation is stored
// initialState2LLDown = array where the monomial spin down, second Landau level representation is stored
// nbrParticlesLLLUp = number of spin up particles in the lowest Landau level
// nbrParticlesLLLDown = number of spin down particles in the lowest Landau level
// return value = state in its fermionic representation

inline unsigned long FermionOnSphereWithSpinTwoLandauLevels::ConvertFromMonomial(unsigned long* initialStateLLLUp, unsigned long* initialState2LLUp,
										 unsigned long* initialStateLLLDown, unsigned long* initialState2LLDown, 
										 int nbrParticlesLLLUp,  int nbrParticlesLLLDown)
{
  unsigned long TmpState = 0x0ul;  
  for (int j = 0; j < nbrParticlesLLLUp; ++j)
    TmpState |= 0x4ul << (initialStateLLLUp[j] << 4);
  for (int j = nbrParticlesLLLUp; j < this->NbrFermionsUp; ++j)
    TmpState |= 0x8ul << (initialState2LLUp[j] << 4);
  for (int j = 0; j < nbrParticlesLLLDown; ++j)
    TmpState |= 0x1ul << (initialStateLLLDown[j] << 4);
  for (int j = nbrParticlesLLLDown; j < this->NbrFermionsDown; ++j)
    TmpState |= 0x2ul << (initialState2LLDown[j] << 4);
  return TmpState;
}

#endif


