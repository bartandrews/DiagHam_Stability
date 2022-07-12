////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class of spin chain generalized AKLT hamiltonian with translations     //
//             using local projectors on the spin 3 and spin 4 sectors        //
//                                                                            //
//                        last modification : 11/11/2016                      //
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


#ifndef SPINCHAINAKLTP3P4HAMILTONIANWITHTRANSLATIONS_H
#define SPINCHAINAKLTP3P4HAMILTONIANWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChainWithTranslations.h"
#include "Hamiltonian/SpinChainAKLTHamiltonianWithTranslations.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class SpinChainAKLTP3P4HamiltonianWithTranslations : public SpinChainAKLTHamiltonianWithTranslations 
{

 protected:
  
  // numerical factor in front of the spin 3 projector
  double P3Factor;
  // numerical factor in front of the spin 4 projector
  double P4Factor;

 public:

  // default constructor
  //
  SpinChainAKLTP3P4HamiltonianWithTranslations();

  // constructor from default datas
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpin = number of spin
  // p3Factor = numerical factor in front of the spin 3 projector
  // p4Factor = numerical factor in front of the spin 4 projector
  SpinChainAKLTP3P4HamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin, double p3Factor = 1.0, double p4Factor = 1.0);

  // destructor
  //
  ~SpinChainAKLTP3P4HamiltonianWithTranslations();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set chain
  // 
  // chain = pointer on Hilbert space of the associated system
  // return value = reference on current Hamiltonian
  SpinChainAKLTP3P4HamiltonianWithTranslations& SetChain(AbstractSpinChainWithTranslations* chain);

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

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
