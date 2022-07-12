////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class of pair-hopping hanmiltonian written in the spin-1 language      //
//                        with translations at k=0 or pi                      //
//                                                                            //
//                        last modification : 13/03/2019                      //
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


#ifndef PAIRHOPPINGREALHAMILTONIANWITHTRANSLATIONS_H
#define PAIRHOPPINGREALHAMILTONIANWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChainWithTranslations.h"
#include "Hamiltonian/SpinChainRealHamiltonianWithTranslations.h"

#include <iostream>


using std::ostream;


class PairHoppingRealHamiltonianWithTranslations : public SpinChainRealHamiltonianWithTranslations 
{

 protected:
  
  // value that defines the filling factor p/(2p+1)
  int PValue;
  
  // optional global energy shift
  double GlobalEnergyShift;

 public:

  // default constructor
  //
  PairHoppingRealHamiltonianWithTranslations();

  // constructor from default data
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpin = number of spin
  // pValue = value that defines the filling factor p/(2p+1)
  PairHoppingRealHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin, int pValue);

  // destructor
  //
  ~PairHoppingRealHamiltonianWithTranslations();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

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
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
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
  virtual RealVector* LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
					       int firstComponent, int nbrComponent);

 
 protected:
 
   // evaluate all cosinus/sinus that are needed when computing matrix elements
  //
  void EvaluateCosinusTable();

};

#endif
