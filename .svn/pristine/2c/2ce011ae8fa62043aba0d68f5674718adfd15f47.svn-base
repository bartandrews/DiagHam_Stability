////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2005 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with spin without            //
//                            sign precalculation table                       //
//        supporting system partial polarization (i.e. an infinite Zeeman     //
//                            energy on some orbitals)                        //
//                                                                            //
//                        last modification : 03/10/2018                      //
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


#ifndef FERMIONONSPHEREWITHSPINPARTIALPOLARIZATION_H
#define FERMIONONSPHEREWITHSPINPARTIALPOLARIZATION_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"

#include <iostream>
#include <map>

using std::map;

using std::cout;
using std::endl;
using std::ostream;

class FermionOnSphere;


class FermionOnSphereWithSpinPartialPolarization :  public FermionOnSphereWithSpin
{

  friend class FermionOnSphereWithSpinSzProjection;
  friend class FermionOnSphereWithSpinLzSzSymmetry;
  friend class FermionOnSphereWithSpinSzSymmetry;
  friend class FermionOnSphereWithSpinLzSymmetry;
  friend class FermionOnSphereWithSU4Spin;
  friend class FermionOnSphereWithSU4SpinLong;
  friend class FermionOnSphereWithSpinAllSz;
  friend class FermionOnSphereWithSpinHaldaneBasis;
  friend class FermionOnSphereWithSpinAllSzLzSymmetry;
  friend class FermionOnSphereWithSpinHaldaneLzSzSymmetry;

  friend class BosonOnSphereTwoLandauLevels;
  friend class BosonOnSphereShort;
  friend class BosonOnSphereWithSU2Spin;
  
 protected:

  // number of orbitals (from the lowest momenta) that are fully spin polarized (i.e. where only spin up are allowed)
  int NbrSpinPolarizedOrbitals;

 public:

  // default constructor
  //
  FermionOnSphereWithSpinPartialPolarization();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalSpin = twice the total spin value
  // nbrSpinPolarizedOrbitals = number of orbitals (from the lowest momenta) that are fully spin polarized (i.e. where only spin up are allowed)
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinPartialPolarization (int nbrFermions, int totalLz, int lzMax, int totalSpin, int nbrSpinPolarizedOrbitals, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpinPartialPolarization(const FermionOnSphereWithSpinPartialPolarization& fermions);

  // destructor
  //
  ~FermionOnSphereWithSpinPartialPolarization ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpinPartialPolarization& operator = (const FermionOnSphereWithSpinPartialPolarization& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  void SetTargetSpace(ParticleOnSphereWithSpin* targetSpace);
  
  
 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // nbrNUp = number of particles with quantum number up
  // nbrNDown = number of particles with quantum number down
  // return value = Hilbert space dimension      
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int nbrNUp, int nbrNDown);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion in the state
  // totalLz = momentum total value
  // nbrNUp = number of particles with quantum number up
  // nbrNDown = number of particles with quantum number down
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int lzMax, int totalLz, int nbrNUp, int nbrNDown, long pos);

};

#endif


