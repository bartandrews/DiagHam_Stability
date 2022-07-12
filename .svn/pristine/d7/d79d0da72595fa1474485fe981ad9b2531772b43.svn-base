////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of spin chain hamiltonian with fully anisotropic          //
//                  coupling constants and a magnetic field                   //
//                         along the x and z directions                       //
//                                                                            //
//                        last modification : 26/08/2015                      //
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


#ifndef SPINCHAINREALFULLHAMILTONIAN_H
#define SPINCHAINREALFULLHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Hamiltonian/AbstractHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class SpinChainRealFullHamiltonian : public AbstractHamiltonian
{

 protected:
  
  //pointer to Hilbert space of the associated system
  AbstractSpinChain* Chain;

  // length of the spin chain
  int NbrSpin;

  // array containing coupling constants between spins along x + y
  double* JPlusFactor;
  // array containing coupling constants between spins along x - y 
  double* JMinusFactor;
  // array containing coupling constants between spins along z
  double* JzFactor;

  // amplitude of the Zeeman term along x
  double* HxFactor;
  // amplitude of the Zeeman term along z
  double* HzFactor;

  // flag to indicate if  periodic boundary conditions should be used
  bool PeriodicBoundaryConditions;

  // array to store the diagonal contribution of the Hamiltonian
  double* SzSzContributions;

 public:

  // constructor from default data
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpin = number of spin
  // jx = array containing coupling constants between spins along x
  // jy = array containing coupling constants between spins along x
  // jz = array containing coupling constants between spins along x
  // hxFactor = array containing the amplitude of the Zeeman term along x
  // hzFactor = array containing the amplitude of the Zeeman term along z
  // periodicBoundaryConditions = true if periodic boundary conditions have to be used
  SpinChainRealFullHamiltonian(AbstractSpinChain* chain, int nbrSpin, double* jx, double* jy, double* jz, 
			       double* hx, double* hz, bool periodicBoundaryConditions = false);

  // destructor
  //
  ~SpinChainRealFullHamiltonian();

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
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

};

#endif
