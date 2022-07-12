////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of Potts 3 chain hamiltonian                     //
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


#ifndef POTTS3CHAINDUALOBRIENFENDLEYHAMILTONIAN_H
#define POTTS3CHAINDUALOBRIENFENDLEYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/Potts3Chain.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "MathTools/Complex.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class Potts3ChainDualOBrienFendleyHamiltonian : public AbstractHamiltonian
{

 protected:
  
  Potts3Chain* Chain;

  
  // true if the chain is periodic
  bool PeriodicFlag;

  // number of sites
  int NbrSpin;
  // number of sites minus one
  int ReducedNbrSpin;

  double* SzSzContributions;

  // flag for fast multiplication algorithm
  bool FastMultiplicationFlag;
  // number of non-null term in the hamiltonian for each state
  int NbrInteractionPerComponent;
  // index of the state obtained for each term of the hamiltonian when applying on a given state
  int* InteractionPerComponentIndex;

 public:

  // default constructor
  //
  Potts3ChainDualOBrienFendleyHamiltonian();

  // constructor from default datas
  //
  // chain = reference on Hilbert space of the associated system
  // nbrSpin = number of spin
  // periodicFlag = true if the chain is periodic
  // memory = amount of memory that can be used from precalculations (in bytes)
  Potts3ChainDualOBrienFendleyHamiltonian(Potts3Chain* chain, int nbrSpin, bool periodicFlag, long memory);

  // destructor
  //
  ~Potts3ChainDualOBrienFendleyHamiltonian();

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
  
 private:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

};

#endif
