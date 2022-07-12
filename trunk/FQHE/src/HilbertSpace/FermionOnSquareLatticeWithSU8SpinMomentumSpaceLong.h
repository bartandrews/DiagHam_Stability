////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(8) spin         //
//                    in momentum space for more than 8 orbitals              //
//                                                                            //
//                        last modification : 15/05/2020                      //
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


#ifndef FERMIONONSQUARELATTICEWITHSU8SPINMOMENTUMSPACELONG_H
#define FERMIONONSQUARELATTICEWITHSU8SPINMOMENTUMSPACELONG_H

#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSU8SpinLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"


#include <iostream>



class FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong : public FermionOnSphereWithSU8SpinLong
{

 protected:

  // number of sites in the x direction
  int NbrSiteX;
  // number of sites in the y direction
  int NbrSiteY;

  // momentum along the x direction
  int KxMomentum;
  // momentum along the y direction
  int KyMomentum;

  // flag to indicate that the Hilbert space should preserve Sz (i.e. N_{1+2+3+4} and N_{5+6+7+8})
  bool SzFlag;
  // flag to indicate that the Hilbert space should preserve the SU(4) sector (i.e. N_{1+2}, N_{3+4}, N_{5+6} and N_{7+8})
  bool SU4Flag;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory = 10000000);
  
  // constructor when preserving only a SU(4) degree of freedom 
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // nbrParticles12 = number of particles with sigma=1 or sigma=2
  // nbrParticles34 = number of particles with sigma=3 or sigma=4
  // nbrParticles56 = number of particles with sigma=5 or sigma=6
  // nbrParticles78 = number of particles with sigma=7 or sigma=8
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum,
						  int nbrParticles12, int nbrParticles34, int nbrParticles56, int nbrParticles78,
						  unsigned long memory = 10000000);

  // constructor when conserving all 8 Cartan related quantum numbers
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // nbrParticleSigma = array that provides the number of particles with a given internal degree of freedom
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY,
						  int kxMomentum, int kyMomentum,
						  int* nbrParticleSigma, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong(const FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong& operator = (const FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

 protected:
 
  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy);
  
  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // nbrParticles1 = number of particles with sigma=1
  // nbrParticles2 = number of particles with sigma=2
  // nbrParticles3 = number of particles with sigma=3
  // nbrParticles4 = number of particles with sigma=4
  // nbrParticles5 = number of particles with sigma=5
  // nbrParticles6 = number of particles with sigma=6
  // nbrParticles7 = number of particles with sigma=7
  // nbrParticles8 = number of particles with sigma=8
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy,
					     int nbrParticles1, int nbrParticles2, int nbrParticles3, int nbrParticles4,
					     int nbrParticles5, int nbrParticles6, int nbrParticles7, int nbrParticles8);

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // nbrParticles12 = number of particles with sigma=1 or sigma=2
  // nbrParticles34 = number of particles with sigma=3 or sigma=4
  // nbrParticles56 = number of particles with sigma=5 or sigma=6
  // nbrParticles78 = number of particles with sigma=7 or sigma=8
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy,
					     int nbrParticles12, int nbrParticles34, int nbrParticles56, int nbrParticles78);
  
  // evaluate the Hilbert space dimension for the 8 component case from the one component Hilbert space dimensions
  //
  // nbrParticles1 = number of particles with sigma=1
  // nbrParticles2 = number of particles with sigma=2
  // nbrParticles3 = number of particles with sigma=3
  // nbrParticles4 = number of particles with sigma=4
  // nbrParticles5 = number of particles with sigma=5
  // nbrParticles6 = number of particles with sigma=6
  // nbrParticles7 = number of particles with sigma=7
  // nbrParticles8 = number of particles with sigma=8
  virtual long EvaluateHilbertSpaceDimension(int nbrParticles1, int nbrParticles2, int nbrParticles3, int nbrParticles4,
					     int nbrParticles5, int nbrParticles6, int nbrParticles7, int nbrParticles8);

  // evaluate Hilbert space dimension for fermions for a single band
  //
  // nbrParticles = number of nbrParticles
  // kxMomentum = total momentum along x
  // kyMomentum = total momentum along y
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension  
  long EvaluateHilbertSpaceDimensionOneBand(int nbrParticles, int kxMomentum, int kyMomentum,
					    int currentKx, int currentKy, int currentTotalKx, int currentTotalKy);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos);
  
  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // nbrParticles1 = number of particles with sigma=1
  // nbrParticles2 = number of particles with sigma=2
  // nbrParticles3 = number of particles with sigma=3
  // nbrParticles4 = number of particles with sigma=4
  // nbrParticles5 = number of particles with sigma=5
  // nbrParticles6 = number of particles with sigma=6
  // nbrParticles7 = number of particles with sigma=7
  // nbrParticles8 = number of particles with sigma=8
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy,
					     int nbrParticles1, int nbrParticles2, int nbrParticles3, int nbrParticles4,
					     int nbrParticles5, int nbrParticles6, int nbrParticles7, int nbrParticles8, long pos);

  // generate all states corresponding to the constraints from the one component Hilbert spaces
  // 
  // nbrParticles1 = number of particles with sigma=1
  // nbrParticles2 = number of particles with sigma=2
  // nbrParticles3 = number of particles with sigma=3
  // nbrParticles4 = number of particles with sigma=4
  // nbrParticles5 = number of particles with sigma=5
  // nbrParticles6 = number of particles with sigma=6
  // nbrParticles7 = number of particles with sigma=7
  // nbrParticles8 = number of particles with sigma=8
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrParticles1, int nbrParticles2, int nbrParticles3, int nbrParticles4,
			      int nbrParticles5, int nbrParticles6, int nbrParticles7, int nbrParticles8);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // nbrParticles12 = number of particles with sigma=1 or sigma=2
  // nbrParticles34 = number of particles with sigma=3 or sigma=4
  // nbrParticles56 = number of particles with sigma=5 or sigma=6
  // nbrParticles78 = number of particles with sigma=7 or sigma=8
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy,
			      int nbrParticles12, int nbrParticles34, int nbrParticles56, int nbrParticles78, long pos);

  // generate all states corresponding to the constraints for a single band
  // 
  // stateDescriptions = array where the many-body basis configurations will be stored
  // nbrFermions = number of fermions
  // kxMomentum = total momentum along x
  // kyMomentum = total momentum along y
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStatesOneBand(ULONGLONG* stateDescriptions, int nbrFermions, int kxMomentum, int kyMomentum,
				     int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos);
  
};


#endif


