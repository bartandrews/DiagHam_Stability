////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                        class author: Cecile Repellin                       //
//                                                                            //
//                                                                            //
//                  class of particle on lattice 1-body operator              //
//                                                                            //
//                        last modification : 10/04/2018                      //
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

#ifndef PARTICLEONLATTICEREALSPACEAND2DTRANSLATIONONEBODYOPERATOR_H
#define PARTICLEONLATTICEREALSPACEAND2DTRANSLATIONONEBODYOPERATOR_H



#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/ParticleOnSphereDensityOperator.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"

class ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator : public ParticleOnSphereDensityOperator
{

 protected:
  
 public:
  
  // constructor from default datas
  //
  // particle = hilbert space associated to the particles
  // creationIndex = index of the creation operator
  // annihilationIndex = index of the annihilation operator
  ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator(ParticleOnSphere* particle, int creationIndex=0,int annihilationIndex=0);

  // copy constructor
  //
  // oper = reference on the operator to copy
  ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator(const ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator & oper);

  // destructor
  //
  ~ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator();
  
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

  // change indices of creation / annihilation operators
  // creationIndex = index of the creation operator
  // annihilationIndex = index of the annihilation operator
  void SetCreationAnnihilationIndex (int creationIndex, int annihilationIndex);
  
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
