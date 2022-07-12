////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class density-density operator for particle with spin          //
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


#ifndef PARTICLEONSPHEREWITHSPINALLSZDENSITYODDCHANNEL_H
#define PARTICLEONSPHEREWITHSPINALLSZDENSITYODDCHANNEL_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
//#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"

class ParticleOnSphereWithSpinAllSzDensityOddChannel : public AbstractOperator
{

 protected:

  // hilbert space associated to the particles
  FermionOnSphereWithSpinAllSz* Particle;
  
  // momentum index of the leftmost creation operator
  int CreationMomentumIndex1;
  // momentum index of the leftmost annihilation operator
  int AnnihilationMomentumIndex1;
  
 public:
  
  // constructor from default datas
  //
  // particle = hilbert space associated to the particles
  // creationMomentumIndex1 = momentum index of the leftmost creation operator (from 0 to 2S)
  // creationMomentumIndex2 = momentum index of the leftmost creation operator (from 0 to 2S)
  // annihilationMomentumIndex1 = momentum index of the leftmost annihilation operator (from 0 to 2S)
  // annihilationMomentumIndex2 = momentum index of the rightmost annihilation operator(from 0 to 2S)
  ParticleOnSphereWithSpinAllSzDensityOddChannel(FermionOnSphereWithSpinAllSz* particle,  int creationMomentumIndex1,  int annihilationMomentumIndex1);

  // destructor
  //
  ~ParticleOnSphereWithSpinAllSzDensityOddChannel();
  
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
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);
   
  // multiply a vector by the current operator for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
				  int firstComponent, int nbrComponent);
  

};

#endif
