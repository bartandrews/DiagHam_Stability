////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of two dimension Ising hamiltonian with transverse field     //
//                             and 2d translations                            //
//                                                                            //
//                        last modification : 18/12/2015                      //
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


#ifndef TWODIMENSIONALTRANSVERSEFIELDAND2DTRANSLATIONISINGHAMILTONIAN_H
#define TWODIMENSIONALTRANSVERSEFIELDAND2DTRANSLATIONISINGHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Hamiltonian/AbstractHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class TwoDimensionalTransverseFieldAnd2DTranslationIsingHamiltonian : public AbstractHamiltonian
{

 protected:
  
  //pointer to Hilbert space of the associated system
  AbstractSpinChain* Chain;

  // total number of spins
  int NbrSpin;
  // number of spin chain along the x direction
  int NbrSpinX;
  // number of spin chain along the y direction
  int NbrSpinY;

  // amplitude of the Ising term
  double JFactor;

  // half amplitude of the Zeeman term along x
  double HxFactor;
  // amplitude of the Zeeman term along z
  double HzFactor;

  // array to store the diagonal contribution of the Hamiltonian
  double* SzSzContributions;

  // momentum along the x direction
  int XMomentum;
  // momentum along the y direction
  int YMomentum;

  //array containing all the phase factors that are needed when computing matrix elements
  Complex** ExponentialFactors;

 public:

  // constructor from default data
  //
  // chain = pointer to Hilbert space of the associated system
  // xMomentum = momentum along the x direction
  // nbrSpinX = number of spin along the x direction
  // yMomentum = momentum along the y direction
  // nbrSpinY = number of spin along the y direction
  // jFactor = amplitude of the Ising term
  // hxFactor = amplitude of the Zeeman term along x
  // hzFactor = amplitude of the Zeeman term along z
  TwoDimensionalTransverseFieldAnd2DTranslationIsingHamiltonian(AbstractSpinChain* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, double jFactor,
								double hxFactor, double hzFactor);

  // destructor
  //
  ~TwoDimensionalTransverseFieldAnd2DTranslationIsingHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					     int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						     int firstComponent, int nbrComponent);

 protected:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

  // evaluate all exponential factors
  //   
  virtual void EvaluateExponentialFactors();

  // get a linearized position index from the 2d coordinates
  //
  // xPosition = position along the x direction
  // yPosition = position along the y direction
  // return value = linearized index
  virtual int GetLinearizedIndex(int xPosition, int yPosition);

};

// get a linearized position index from the 2d coordinates
//
// xPosition = position along the x direction
// yPosition = position along the y direction
// return value = linearized index

inline int TwoDimensionalTransverseFieldAnd2DTranslationIsingHamiltonian::GetLinearizedIndex(int xPosition, int yPosition)
{
  return ((xPosition * this->NbrSpinY) + yPosition);
}

#endif
