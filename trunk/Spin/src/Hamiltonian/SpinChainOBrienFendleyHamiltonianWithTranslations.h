////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of spin chain hamiltonian implementing the              //
//                   O'Brien-Fendley model with translations                  //
//                                                                            //
//                        last modification : 02/10/2019                      //
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


#ifndef SPINCHAINOBRIENFENDLEYHAMILTONIANWITHTRANSLATIONS_H
#define SPINCHAINOBRIENFENDLEYHAMILTONIANWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChainWithTranslations.h"
#include "Hamiltonian/SpinChainAKLTHamiltonianWithTranslations.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class SpinChainOBrienFendleyHamiltonianWithTranslations : public SpinChainAKLTHamiltonianWithTranslations
{

 protected:
  
 public:

  // default constructor
  //
  SpinChainOBrienFendleyHamiltonianWithTranslations();

  // constructor from default data
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpin = number of spin
  SpinChainOBrienFendleyHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin);

  // destructor
  //
  ~SpinChainOBrienFendleyHamiltonianWithTranslations();

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
  // and store result in another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						  int firstComponent, int nbrComponent);

 
 protected:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

};

#endif
