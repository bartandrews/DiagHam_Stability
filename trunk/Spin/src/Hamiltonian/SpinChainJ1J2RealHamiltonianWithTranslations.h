////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of spin chain J1-J2hamiltonian with translations           //
//                         at inversion symmetric points                      //
//                                                                            //
//                        last modification : 13/07/2016                      //
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


#ifndef SPINCHAINJ1J2REALHAMILTONIANWITHTRANSLATIONS_H
#define SPINCHAINJ1J2REALHAMILTONIANWITHTRANSLATIONS_H


#include "config.h"
#include "Hamiltonian/SpinChainRealHamiltonianWithTranslations.h"

#include <iostream>


using std::ostream;



class SpinChainJ1J2RealHamiltonianWithTranslations : public SpinChainRealHamiltonianWithTranslations
{

 protected:

  // coupling constant between spin in the xx and yy direction for nearest neighbors
  double J1;
  // coupling constant between spin in the zz direction  for nearest neighbors
  double J1z;
  // half coupling constant between spin in the xx and yy direction  for nearest neighbors
  double HalfJ1;
  // coupling constant between spin in the xx and yy direction for second nearest neighbors
  double J2;
  // coupling constant between spin in the zz direction for second nearest neighbors 
  double J2z;
  // half coupling constant between spin in the xx and yy direction for second nearest neighbors 
  double HalfJ2;
  

 public:

  // default constructor
  //
  SpinChainJ1J2RealHamiltonianWithTranslations();

  // constructor from default data
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpin = number of spin
  // j1 = coupling constants between nearest neighbor spins
  // j2 = coupling constants between second nearest neighbor spins
  SpinChainJ1J2RealHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin, double j1, double j2);

  // destructor
  //
  ~SpinChainJ1J2RealHamiltonianWithTranslations();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

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
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

};

#endif
