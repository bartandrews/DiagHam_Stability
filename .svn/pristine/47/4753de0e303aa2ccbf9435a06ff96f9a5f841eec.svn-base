////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of spin chain hamiltonian with a triple product term       //
//                                                                            //
//                        last modification : 21/06/2015                      //
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


#ifndef SPINCHAINTRIPLEPRODUCTHAMILTONIAN_H
#define SPINCHAINTRIPLEPRODUCTHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Hamiltonian/AbstractHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class SpinChainTripleProductHamiltonian : public AbstractHamiltonian
{

 protected:
  
  AbstractSpinChain* Chain;

  // array containing coupling constants between spins along x and z
  double* J;
  double* HalfJ;

  // array containing coupling constants between spins along z
  double* Jz;

  // array containing the coupling constants of the triple product
  double* Chi;
  double* HalfChi;

  // amplitude of the Zeeman term along z
  double* Hz;

  // flag to indicate if  periodic boundary conditions should be used
  bool PeriodicBoundaryConditions;

  // length of the spin chain
  int NbrSpin;

  // array to store the diagonal contribution of the Hamiltonian
  double* SzSzContributions;

 public:

  // constructor from default datas
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpin = number of spin
  // j = array containing coupling constants between spins
  // chi = array containing the coupling constants of the triple product
  // periodicBoundaryConditions = true if periodic boundary conditions have to be used
  SpinChainTripleProductHamiltonian(AbstractSpinChain* chain, int nbrSpin, double* j, double* chi, bool periodicBoundaryConditions = false);

  // constructor from default datas
  //
  // chain = reference on Hilbert space of the associated system
  // nbrSpin = number of spin
  // j = array containing coupling constants between spins along x and z
  // jz = array containing coupling constants between spins along z
  // chi = array containing the coupling constants of the triple product
  // periodicBoundaryConditions = true if periodic boundary conditions have to be used
  SpinChainTripleProductHamiltonian(AbstractSpinChain* chain, int nbrSpin, double* j, double* jz, double* chi, bool periodicBoundaryConditions = false);

  // constructor from default datas
  //
  // chain = reference on Hilbert space of the associated system
  // nbrSpin = number of spin
  // j = array containing coupling constants between spins along x and z
  // jz = array containing coupling constants between spins along z
  // chi = array containing the coupling constants of the triple product
  // hz = array containing the amplitude of the Zeeman term along z
  // periodicBoundaryConditions = true if periodic boundary conditions have to be used
  SpinChainTripleProductHamiltonian(AbstractSpinChain* chain, int nbrSpin, double* j, double* jz, double* chi, double* hz, bool periodicBoundaryConditions = false);

  // destructor
  //
  ~SpinChainTripleProductHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set chain
  // 
  // chain = pointer on Hilbert space of the associated system
  // return value = reference on current Hamiltonian
  SpinChainTripleProductHamiltonian& SetChain(AbstractSpinChain* chain);

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
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
				     int firstComponent, int nbrComponent);

 protected:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

};

#endif
