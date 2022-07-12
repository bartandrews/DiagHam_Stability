////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of XYZ chain with a natural boundary term             //
//                                                                            //
//                        last modification : 27/10/2014                      //
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


#ifndef SPINCHAINXYZNATURALBOUNDARYTERMHAMILTONIAN_H
#define SPINCHAINXYZNATURALBOUNDARYTERMHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/Spin1_2Chain.h"
#include "Hamiltonian/SpinChainXYZHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class SpinChainXYZNaturalBoundaryTermHamiltonian : public SpinChainXYZHamiltonian
{

 protected:

  // perturbation order for the edge mode development
  int PerturbationOrder;  
  // true if the boundary term is Jy dominated instead of Jx dominated
  bool JyDominatedBoundaryTerm;

  
 public:

  // constructor
  //
  // chain = pointer to the Hilbert space of the system
  // nbrSpin = number of spins
  // jxFactor = coupling along the x direction
  // jyFactor = coupling along the y direction
  // jzFactor = coupling along the z direction
  // hFactor = Zeeman term 
  // boundaryCondition = boundary condition to apply (0 for open chain, 1 for periodic, -1 for antiperiodic)
  // perturbationOrder = perturbation order for the edge mode development
  // fixedParityFlag = true if the parity if fixed for the Hilbert space
  // fixedParity = value of the parity if fixed for the Hilbert space
  // jyDominatedBoundaryTerm = true if the boundary term is Jy dominated instead of Jx dominated
  SpinChainXYZNaturalBoundaryTermHamiltonian(Spin1_2Chain* chain, int nbrSpin, 
					     double jxFactor, double jyFactor, double jzFactor, double hFactor,
					     double boundaryCondition = 0.0, int perturbationOrder = 0, 
					     bool fixedParityFlag = false, int fixedParity = 0, bool jyDominatedBoundaryTerm = false);
  
  // destructor
  //
  ~SpinChainXYZNaturalBoundaryTermHamiltonian();

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
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						  int firstComponent, int nbrComponent);


 protected:

  // evaluate diagonal matrix elements
  // 
  virtual void EvaluateDiagonalMatrixElements();

 

};

#endif
