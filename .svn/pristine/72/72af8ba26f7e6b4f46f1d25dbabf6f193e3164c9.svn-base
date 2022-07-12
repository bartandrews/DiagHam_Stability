////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of Potts 3 chain hamiltonian with translations             //
//                     for the dual O'Brien-Fendley model                     //
//                                                                            //
//                        last modification : 02/12/2019                      //
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


#ifndef POTTS3CHAINDUALOBRIENFENDLEYHAMILTONIANWITHTRANSLATIONS_H
#define POTTS3CHAINDUALOBRIENFENDLEYHAMILTONIANWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/Potts3ChainWithTranslations.h"
#include "Hamiltonian/SpinChainHamiltonianWithTranslations.h"
#include "MathTools/Complex.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class Potts3ChainDualOBrienFendleyHamiltonianWithTranslations : public SpinChainHamiltonianWithTranslations
{

 protected:
  
  // number of sites minus one
  int ReducedNbrSpin;

 public:

  // constructor from default datas
  //
  // chain = reference on Hilbert space of the associated system
  // nbrSpin = number of spin
  Potts3ChainDualOBrienFendleyHamiltonianWithTranslations(Potts3ChainWithTranslations* chain, int nbrSpin);

  // destructor
  //
  ~Potts3ChainDualOBrienFendleyHamiltonianWithTranslations();

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

 
 protected:
 
  // evaluate all cosinus/sinus that are needed when computing matrix elements
  //
  virtual void EvaluateCosinusTable();

  // evaluate all matrix elements
  //   
  virtual void EvaluateDiagonalMatrixElements();

};

#endif
