////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of particle on lattice 1-body operator              //
//                                                                            //
//                        last modification : 09/04/2008                      //
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

#ifndef PARTICLEONLATTICETRANSLATIONOPERATOR_H
#define PARTICLEONLATTICETRANSLATIONOPERATOR_H

#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/ParticleOnLattice.h"

class ParticleOnLatticeTranslationOperator : public AbstractOperator
{

 protected:

  // hilbert space associated to the particles
  ParticleOnLattice* Particle;

  // x-component of the translation
  int Rx;
  // y-component of the translation
  int Ry;
  
  // indices of the creation operators
  int* CreationIndices;
  // index of the annihilation operators
  int* AnnihilationIndices;
  // complex phase of the operator pairs from magnetic translational invariance
  Complex* TranslationPhases;
  
 public:
  
  // constructor from default datas
  //
  // particle = hilbert space associated to the particles
  // rx = x-component of the desired translation
  // ry = y-component of the desired translation
  ParticleOnLatticeTranslationOperator(ParticleOnLattice* particle, int rx=0, int ry=1);

  // copy constructor
  //
  // oper = reference on the operator to copy
  ParticleOnLatticeTranslationOperator(const ParticleOnLatticeTranslationOperator& oper);

  // destructor
  //
  ~ParticleOnLatticeTranslationOperator();
  
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

  // set components of translation vector
  //
  // rx = x-component of the desired translation
  // ry = y-component of the desired translation
  void SetTranslationComponents(int rx, int ry);

  // evaluate part of the matrix element, within a given of indices
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = corresponding matrix element
  Complex PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent);

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
