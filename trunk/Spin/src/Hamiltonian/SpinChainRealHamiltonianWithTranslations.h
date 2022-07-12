////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of spin chain hamiltonian with translations             //
//                         at inversion symmetric points                      //
//                                                                            //
//                        last modification : 23/05/2016                      //
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


#ifndef SPINCHAINREALHAMILTONIANWITHTRANSLATIONS_H
#define SPINCHAINREALHAMILTONIANWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChainWithTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"

#include <iostream>


using std::ostream;



class SpinChainRealHamiltonianWithTranslations : public AbstractHamiltonian
{

 protected:
  
  AbstractSpinChainWithTranslations* Chain;

  // coupling constant between spin in the xx and yy direction 
  double J;
  // coupling constant between spin in the zz direction 
  double Jz;
  // half coupling constant between spin in the xx and yy direction 
  double HalfJ;
  // coupling constant between next nearest neighbour spins in the zz direction 
  double NNNCoupling;

  
  // number of spin 
  int NbrSpin;

  // array conating all matrix diagonal elements
  double* SzSzContributions;

  //array containing all the cosinus that are needed when computing matrix elements
  double* CosinusTable;
  //array containing all the sinus that are needed when computing matrix elements
  double* SinusTable;
  //array containing all the complex phase that are needed when computing matrix elements
  double* ExponentialTable;

 public:

  // default constructor
  //
  SpinChainRealHamiltonianWithTranslations();

  // constructor from default datas
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpin = number of spin
  // j = coupling constants between spins
  // nnCoupling = term to add to ZZ nearest-neighbour interaction
  // nnnCoupling = nearest-neighbour interaction in Z direction
  SpinChainRealHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin, double j, double nnCoupling, double nnnCoupling);

  // destructor
  //
  ~SpinChainRealHamiltonianWithTranslations();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set chain
  // 
  // chain = pointer on Hilbert space of the associated system
  // return value = reference on current Hamiltonian
  SpinChainRealHamiltonianWithTranslations& SetChain(AbstractSpinChainWithTranslations* chain);

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
 
 protected:
 
  // evaluate all cosinus/sinus that are needed when computing matrix elements
  //
  void EvaluateCosinusTable();

  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

};

#endif
