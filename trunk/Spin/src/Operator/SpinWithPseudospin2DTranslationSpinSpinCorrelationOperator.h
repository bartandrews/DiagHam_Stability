////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                      Class author Cecile Repellin                          //
//                                                                            //
//                                                                            //
//              class of spin-spin operator for spin-spin model with 2D       //
//                             translations                                   //
//                                                                            //
//                        last modification : 21/03/2017                      //
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


#ifndef SPINWITHPSEUDOSPIN2DTRANSLATIONSPINSPINCORRELATIONOPERATOR_H
#define SPINWITHPSEUDOSPIN2DTRANSLATIONSPINSPINCORRELATIONOPERATOR_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation.h"

#include <iostream>




class SpinWithPseudospin2DTranslationSpinSpinCorrelationOperator : public AbstractOperator
{

 protected:

  // pointer to the Hilbert space
  Spin1_2ChainWithPseudospinAnd2DTranslation* Chain;
  // number of spins
  int NbrSpin;
  
 
  int SiteIndex1;
  
  int SiteIndex2;
  
  
  // number of spin chain along the x direction
  int NbrSpinX;
  // number of spin chain along the y direction
  int NbrSpinY;
  // momentum along the x direction
  int XMomentum;
  // momentum along the y direction
  int YMomentum;
  
  Complex** ExponentialFactors;

 public:
  
  // constructor from default datas
  //
  // chain = pointer to the Hilbert space
  // nbrSpin = number of spins
  SpinWithPseudospin2DTranslationSpinSpinCorrelationOperator(Spin1_2ChainWithPseudospinAnd2DTranslation* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, int siteIndex1, int siteIndex2);

  // destructor
  //
  ~SpinWithPseudospin2DTranslationSpinSpinCorrelationOperator();
  
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
  ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
			       int firstComponent, int nbrComponent);
  
protected:
  
  // evaluate all exponential factors
  //   
  virtual void EvaluateExponentialFactors();

};

#endif
