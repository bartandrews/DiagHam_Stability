////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of Potts 3 chain hamiltonian with a natural boundary term      //
//                                                                            //
//                        last modification : 15/12/2014                      //
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


#ifndef POTTS3CHAINNATURALBOUNDARYTERMHAMILTONIAN_H
#define POTTS3CHAINNATURALBOUNDARYTERMHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/Potts3ChainHamiltonian.h"
#include "MathTools/Complex.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class Potts3ChainNaturalBoundaryTermHamiltonian : public Potts3ChainHamiltonian
{

 protected:
  
  // perturbation order for the edge mode development
  int PerturbationOrder;  

  // first factor coming from the filter function when using the first order correction
  double FilterFunctionComponent0;
  // second factor coming from the filter function when using the first order correction
  double FilterFunctionComponent1;
  // third factor coming from the filter function when using the first order correction
  double FilterFunctionComponent2;

 public:

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
  // perturbationOrder = perturbation order for the edge mode development
  // filterFunctionComponent0 = first factor coming from the filter function when using the first order correction
  // filterFunctionComponent1 = second factor coming from the filter function when using the first order correction
  // filterFunctionComponent2 = third factor coming from the filter function when using the first order correction
  // memory = amount of memory that can be used from precalculations (in bytes)
  Potts3ChainNaturalBoundaryTermHamiltonian(Potts3Chain* chain, int nbrSpin, double jFactor, double phiJ, double fFactor, double phiF, 
					    double boundaryCondition, int perturbationOrder, 
					    double filterFunctionComponent0, double filterFunctionComponent1, double filterFunctionComponent2, long memory);

  // destructor
  //
  ~Potts3ChainNaturalBoundaryTermHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

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
  
 private:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

};

#endif
