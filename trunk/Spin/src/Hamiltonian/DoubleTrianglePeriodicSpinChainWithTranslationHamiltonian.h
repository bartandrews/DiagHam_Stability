////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of double triangle spin chain  hamiltonian           //
//                with periodic boundary conditions and translations          //
//                                                                            //
//                        last modification : 25/06/2010                      //
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


#ifndef DOUBLETRIANGLEPERIODICSPINCHAINWITHTRANSLATIONHAMILTONIAN_H
#define DOUBLETRIANGLEPERIODICSPINCHAINWITHTRANSLATIONHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/SpinChainHamiltonianWithTranslations.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian : public SpinChainHamiltonianWithTranslations
{

 protected:
  
  // nearest neighbour coupling constant and half its value
  double J1;
  double HalfJ1;
  // second nearest neighbour coupling constant and half its value
  double J2;
  double HalfJ2;

  // nearest neighbour constant between spins along z
  double Jz1;
  // second nearest neighbour constant between spins along z
  double Jz2;

  // number of triangles
  int NbrTriangles;

 public:

  // constructor from default datas
  //
  // chain = reference on Hilbert space of the associated system
  // nbrSpin = number of spins
  // j1 = nearest neighbour coupling constant
  // j2 = second nearest neighbour coupling constant
  // djz1 = constant added to the nearest neighbour constant between spins along z
  // djz2 = constant added to the second nearest neighbour constant between spins along z
  DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian(AbstractSpinChainWithTranslations* chain, int nbrSpin, double j1, double j2, double djz1, double djz2);

  // destructor
  //
  ~DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian();

  // set chain
  // 
  // chain = pointer on Hilbert space of the associated system
  // return value = reference on current Hamiltonian
  DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian& SetChain(AbstractSpinChainWithTranslations* chain);

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


 protected :
 
  // evaluate all matrix elements
  //   
  virtual void EvaluateDiagonalMatrixElements();

  // evaluate all cosinus/sinus that are needed when computing matrix elements
  //
  virtual void EvaluateCosinusTable();

};

#endif
