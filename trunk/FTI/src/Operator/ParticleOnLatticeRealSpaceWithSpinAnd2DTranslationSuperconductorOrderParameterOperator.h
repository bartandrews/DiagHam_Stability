////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class superconductor order parameter operator              // 
//             for particle with spin on a lattice using translations         //
//                                                                            //
//                        last modification : 09/05/2015                      //
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


#ifndef PARTICLEONLATTICEREALSPACEWITHSPINAND2DTRANSLATIONSUPERCONDUCTORORDERPARAMETEROPERATOR_H
#define PARTICLEONLATTICEREALSPACEWITHSPINAND2DTRANSLATIONSUPERCONDUCTORORDERPARAMETEROPERATOR_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/ParticleOnSphereWithSpinSuperconductorOrderParameterOperator.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"


class ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator : public ParticleOnSphereWithSpinSuperconductorOrderParameterOperator
{

 protected:


 public:
  
  // constructor from default datas
  //
  // particle = hilbert space associated to the particles with the small number of particles
  // creationMomentumIndex1 = momentum index of the leftmost creation operator (from 0 to 2S)
  // creationSymmetryIndex1 = symmetry index of the leftmost creation operator (0 for down, 1 for up)
  // creationMomentumIndex2 = momentum index of the rightmost creation operator (from 0 to 2S)
  // creationSymmetryIndex2 = symmetry index of the rightmost creation operator (0 for down, 1 for up)
  ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator(FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* particle,  int creationMomentumIndex1, int creationSymmetryIndex1,
											 int creationMomentumIndex2, int creationSymmetryIndex2);
  
  // constructor for operator such as a+_{sigma_1,i_1} a+_{sigma_2,i_2} +/- a+_{sigma_3,i_1} a+_{sigma_4,i_2}
  //
  // particle = hilbert space associated to the particles with the small number of particles
  // creationMomentumIndex1 = momentum index of the leftmost creation operator (from 0 to 2S)
  // creationSymmetryIndex1 = symmetry index of the leftmost creation operator (0 for down, 1 for up)
  // creationMomentumIndex2 = momentum index of the rightmost creation operator (from 0 to 2S)
  // creationSymmetryIndex2 = symmetry index of the rightmost creation operator (0 for down, 1 for up)
  // creationSymmetryIndex1SecondTerm = symmetry index of the rightmost creation operator for the second term
  // creationSymmetryIndex2SecondTerm = symmetry index of the leftmost creation operator for the second term
  // sign = sign in front of the second term
  ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator(FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* particle,  int creationMomentumIndex1, int creationSymmetryIndex1,
											 int creationMomentumIndex2, int creationSymmetryIndex2,
											 int creationSymmetryIndex1SecondTerm, int creationSymmetryIndex2SecondTerm, 
											 double sign);

  // copy constructor
  //
  ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator(ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator& oper);

  // destructor
  //
  ~ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator();
  
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
  
  // evaluate part of the matrix element, within a given of indices
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = corresponding matrix element
  Complex PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent);

  // multiply a vector by the current operator for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
				     int firstComponent, int nbrComponent);
  

};

#endif
