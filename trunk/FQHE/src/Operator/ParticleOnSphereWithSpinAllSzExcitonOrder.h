////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of particle on sphere density-density operator           //
//                                                                            //
//                        last modification : 10/12/2002                      //
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


#ifndef PARTICLEONSPHEREWITHSPINALLSZEXCITONORDER_H
#define PARTICLEONSPHEREWITHSPINALLSZEXCITONORDER_H

#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"


class ParticleOnSphereWithSpinAllSzExcitonOrder : public AbstractOperator
{

 protected:

  // hilbert space associated to the particles
  FermionOnSphereWithSpinAllSz* Particle;
  
  // indices attached to the a+_i a+_j a_k a_l
  // index of the leftmost creation operator
  int CreationIndex1;
  // index of the rightmost creation operator
  int CreationIndex2;
  // index of the leftmost annihilation operator
  int AnnihilationIndex1;
  // index of the rightmost annihilation operator
  int AnnihilationIndex2;
  
 public:
  
  // constructor from default datas
  //
  // particle = hilbert space associated to the particles
  // creationIndex1 = index of the leftmost creation operator
  // creationIndex2 = index of the rightmost creation operator
  // annihilationIndex1 = index of the leftmost annihilation operator
  // annihilationIndex2 = index of the rightmost annihilation operator
  ParticleOnSphereWithSpinAllSzExcitonOrder(FermionOnSphereWithSpinAllSz* particle, int creationIndex1, int creationIndex2, int annihilationIndex1, int annihilationIndex2);

  // destructor
  //
  ~ParticleOnSphereWithSpinAllSzExcitonOrder();
  
  // clone operator without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractOperator* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which operator acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where operator acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();

  // Get exciton order parameter <c^*_up,n c^*_down,0 c_up,o c_down,n>
  //
  // vSource = groundstate vector
  // firstIndex = index of the first exciton orbital
  // secondIndex = index of the second exciton orbital
  double GetExcitonOrderLevel2(RealVector& vSource);
  
};

#endif
