////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of Potts 3 chain hamiltonian                     //
//                                                                            //
//                        last modification : 04/06/2012                      //
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


#ifndef POTTS3CHAINHAMILTONIAN_H
#define POTTS3CHAINHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/Potts3Chain.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "MathTools/Complex.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class Potts3ChainHamiltonian : public AbstractHamiltonian
{

 protected:
  
  Potts3Chain* Chain;

  // magnitude of the coupling term
  double JFactor;
  // phase of the coupling term (in PI units)
  double PhiJ; 
  // complex version of the J coupling
  Complex JFullFactor;
  // magnitude of the Zeeman term
  double FFactor;
  // phase of the Zeeman term (in PI units)
  double PhiF;  
  
  // true if the chain is periodic
  bool PeriodicFlag;
  // type of boundary conditions if the chain is periodic (0 for 1, 1 for exp(i 2 \pi / 3), -1 1 for exp(i 2 \pi / 3)) 
  double BoundaryCondition;
  // constant factor when applying periodic boundary conditions 
  Complex BoundaryFactor;

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
  Potts3ChainHamiltonian();

  // constructor from default datas
  //
  // chain = reference on Hilbert space of the associated system
  // nbrSpin = number of spin
  // jFactor = magnitude of the coupling term
  // phiJ = phase of the coupling term (in PI units)
  // fFactor = magnitude of the Zeeman term
  // phiF = phase of the Zeeman term (in PI units)
  // periodicFlag = true if the chain is periodic
  // boundaryCondition = type of boundary conditions if the chain is periodic (0 for 1, 1 for exp(i 2 \pi / 3), -1 1 for exp(i 2 \pi / 3)) 
  // memory = amount of memory that can be used from precalculations (in bytes)
  Potts3ChainHamiltonian(Potts3Chain* chain, int nbrSpin, double jFactor, double phiJ, double fFactor, double phiF, bool periodicFlag, double boundaryCondition, long memory);

  // destructor
  //
  ~Potts3ChainHamiltonian();

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
  
 private:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

};

#endif
