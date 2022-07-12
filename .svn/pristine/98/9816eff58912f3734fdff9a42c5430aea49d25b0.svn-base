////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//      altenate version of the class for bosons on sphere with SU(2) spin    //
//        supporting system partial polarization (i.e. an infinite Zeeman     //
//                            energy on some orbitals                         //
//                                                                            //
//                        last modification : 21/03/2018                      //
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


#ifndef BOSONONSPHEREWITHSU2SPINPARTIALPOLARIZATION_H
#define BOSONONSPHEREWITHSU2SPINPARTIALPOLARIZATION_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"

#include <iostream>

using std::cout;
using std::endl;


class BosonOnSphereWithSU2SpinPartialPolarization :  public BosonOnSphereWithSU2Spin
{

  friend class BosonOnSphereWithSU4Spin;
  friend class BosonOnSphereWithSU2SpinSzSymmetry;
  friend class BosonOnSphereWithSU2SpinLzSymmetry;
  friend class BosonOnSphereWithSU2SpinLzSzSymmetry;
  friend class FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation;

 protected:

  // number of orbitals (from the lowest momenta) that are fully spin polarized (i.e. where only spin up are allowed)
  int NbrSpinPolarizedOrbitals;

 public:

  // default constructor
  // 
  BosonOnSphereWithSU2SpinPartialPolarization ();

  // basic constructor without any constraint on Sz
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // nbrSpinPolarizedOrbitals = number of orbitals (from the lowest momenta) that are fully spin polarized (i.e. where only spin up are allowed)
  BosonOnSphereWithSU2SpinPartialPolarization (int nbrBosons, int totalLz, int lzMax, int nbrSpinPolarizedOrbitals);

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // totalSpin = twice the total spin value
  // nbrSpinPolarizedOrbitals = number of orbitals (from the lowest momenta) that are fully spin polarized (i.e. where only spin up are allowed)
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSU2SpinPartialPolarization (int nbrBosons, int totalLz, int lzMax, int totalSpin, int nbrSpinPolarizedOrbitals, 
					       unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  // 
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSU2SpinPartialPolarization (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereWithSU2SpinPartialPolarization(const BosonOnSphereWithSU2SpinPartialPolarization& bosons);

  // destructor
  //
  virtual ~BosonOnSphereWithSU2SpinPartialPolarization ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereWithSU2SpinPartialPolarization& operator = (const BosonOnSphereWithSU2SpinPartialPolarization& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

 protected:
  
  // read Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description is stored
  // return value = true if no error occured
  virtual bool ReadHilbertSpace (char* fileName);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // nbrNUp = number of particles with quantum number up
  // nbrNDown = number of particles with quantum number down
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int nbrNUp, int nbrNDown);

  // evaluate Hilbert space dimension without the Sz constraint
  //
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // lzMaxUp = momentum maximum value for a boson in the state with up
  // lzMaxDown = momentum maximum value for a boson in the state with down
  // totalLz = momentum total value
  // nbrNUp = number of particles with up
  // nbrNDown = number of particles with down
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int lzMaxUp, int lzMaxDown, int totalLz, 
			      int nbrNUp, int nbrNDown, long pos);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // currentFermionicPositionUp = current fermionic position within the state description for the spin up
  // currentFermionicPositionDown = current fermionic position within the state description for the spin down
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int lzMax, int totalLz, 
			      int currentFermionicPositionUp, int currentFermionicPositionDown, long pos);

};

#endif


