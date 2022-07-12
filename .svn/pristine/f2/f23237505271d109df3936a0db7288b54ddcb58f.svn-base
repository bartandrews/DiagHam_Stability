////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of S^2 operator for aa 1D system with translation         //
//                                                                            //
//                        last modification : 23/03/2016                      //
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


#ifndef SPINWITHTRANSLATIONS2OPERATOR_H
#define SPINWITHTRANSLATIONS2OPERATOR_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/AbstractSpinChainWithTranslations.h"

#include <iostream>


class MathematicaOutput;


class SpinWith1DTranslationS2Operator : public AbstractOperator
{

 protected:

  // pointer to the Hilbert space
  AbstractSpinChainWithTranslations* Chain;
  // number of spins
  int NbrSpin;

 public:
  
  // constructor from default data
  //
  // chain = pointer to the Hilbert space
  // nbrSpin = number of spins
  SpinWith1DTranslationS2Operator(AbstractSpinChainWithTranslations* chain, int nbrSpin);

  // destructor
  //
  ~SpinWith1DTranslationS2Operator();
  
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
  Complex PartialMatrixElement (RealVector& V1, RealVector& V2, 
				long firstComponent, long nbrComponent);

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
  RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
			       int firstComponent, int nbrComponent);

  // multiply a vector by the current operator for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
				  int firstComponent, int nbrComponent);

};

#endif
