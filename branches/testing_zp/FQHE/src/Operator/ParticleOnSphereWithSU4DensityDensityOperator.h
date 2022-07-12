////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of particle on sphere density-density operator           //
//                     with SU(4) internal degree of freedom                  //
//                                                                            //
//                        last modification : 22/10/2007                      //
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


#ifndef PARTICLEONSPHEREWITHSU4DENSITYDENSITYOPERATOR_H
#define PARTICLEONSPHEREWITHSU4DENSITYDENSITYOPERATOR_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/ParticleOnSphereWithSU4Spin.h"


class ParticleOnSphereWithSU4DensityDensityOperator : public AbstractOperator
{

 protected:

  // hilbert space associated to the particles
  ParticleOnSphereWithSU4Spin* Particle;
  
  // indices attached to the a+_{sigma_1,i_1} a+_{sigma_2,i_2} a_{sigma_3,i_3} a_{sigma_4,i_4}
  // momentum index of the leftmost creation operator
  int CreationMomentumIndex1;
  // symmetry index of the leftmost annihilation operator
  int CreationSymmetryIndex1;
  // momentum index of the rightmost creation operator
  int CreationMomentumIndex2;
  // symmetry index of the rightmost annihilation operator
  int CreationSymmetryIndex2;
  // momentum index of the leftmost annihilation operator
  int AnnihilationMomentumIndex1;
  // symmetry index of the leftmost annihilation operator
  int AnnihilationSymmetryIndex1;
  // momentum index of the rightmost annihilation operator
  int AnnihilationMomentumIndex2;
  // symmetry index of the rightmost annihilation operator
  int AnnihilationSymmetryIndex2;

  // additional sign that can arise due to operator reordering
  double SignFactor;
  
 public:
  
  // constructor from default datas
  //
  // particle = hilbert space associated to the particles
  // creationMomentumIndex1 = momentum index of the leftmost creation operator (from 0 to 2S)
  // creationSymmetryIndex1 = symmetry index of the leftmost creation operator (0 for (up,plus), 1 for (up,minus), 2 for (down,plus), 3 for (down,minus))
  // creationMomentumIndex2 = momentum index of the leftmost creation operator (from 0 to 2S)
  // creationSymmetryIndex2 = symmetry index of the rightmost creation operator (0 for (up,plus), 1 for (up,minus), 2 for (down,plus), 3 for (down,minus))
  // annihilationMomentumIndex1 = momentum index of the leftmost annihilation operator (from 0 to 2S)
  // annihilationSymmetryIndex1 = symmetry index of the leftmost annihilation operator (0 for (up,plus), 1 for (up,minus), 2 for (down,plus), 3 for (down,minus))
  // annihilationMomentumIndex2 = momentum index of the rightmost annihilation operator(from 0 to 2S)
  // annihilationSymmetryIndex2 = symmetry index of the rightmost annihilation operator (0 for (up,plus), 1 for (up,minus), 2 for (down,plus), 3 for (down,minus))
  ParticleOnSphereWithSU4DensityDensityOperator(ParticleOnSphereWithSU4Spin* particle,  int creationMomentumIndex1, int creationSymmetryIndex1,
						int creationMomentumIndex2, int creationSymmetryIndex2, int annihilationMomentumIndex1, int annihilationSymmetryIndex1,
						int annihilationMomentumIndex2, int annihilationSymmetryIndex2);

  // destructor
  //
  ~ParticleOnSphereWithSU4DensityDensityOperator();
  
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
  RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
			       int firstComponent, int nbrComponent);
  
};

#endif
