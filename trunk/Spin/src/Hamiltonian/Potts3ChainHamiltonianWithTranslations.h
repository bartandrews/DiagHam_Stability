////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of Potts 3 chain hamiltonian with translations            //
//                                                                            //
//                        last modification : 01/01/2014                      //
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


#ifndef POTTS3CHAINHAMILTONIANWITHTRANSLATIONS_H
#define POTTS3CHAINHAMILTONIANWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/Potts3ChainWithTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "MathTools/Complex.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class Potts3ChainHamiltonianWithTranslations : public AbstractHamiltonian
{

 protected:
  
  Potts3ChainWithTranslations* Chain;

  // magnitude of the coupling term
  double JFactor;
  // phase of the coupling term (in PI units)
  double PhiJ; 
  // complex version of the J coupling
  Complex JFullFactor;
  // complex version of the J coupling, for each possible number of translations
  Complex* JFullFactors;
  // magnitude of the Zeeman term
  double FFactor;
  // phase of the Zeeman term (in PI units)
  double PhiF;  


  // true if the chain is periodic
  bool PeriodicFlag;
  // type of boundary conditions if the chain is periodic (0 for 1, 1 for exp(i 2 \pi / 3), -1 1 for exp(i 2 \pi / 3)) 
  double BoundaryCondition;

  // number of sites
  int NbrSpin;
  // number of sites minus one
  int ReducedNbrSpin;
  // momentum sector
  int Momentum;

  double* SzSzContributions;
  // constant factors for each state when applying periodic boundary conditions 
  Complex* BoundaryFactors;


 public:

  // constructor from default datas
  //
  // chain = reference on Hilbert space of the associated system
  // nbrSpin = number of spin
  // momentum = momentum sector
  // jFactor = magnitude of the coupling term
  // phiJ = phase of the coupling term (in PI units)
  // fFactor = magnitude of the Zeeman term
  // phiF = phase of the Zeeman term (in PI units)
  // boundaryCondition = type of boundary conditions if the chain is periodic (0 for 1, 1 for exp(i 2 \pi / 3), -1 1 for exp(i 2 \pi / 3)) 
  Potts3ChainHamiltonianWithTranslations(Potts3ChainWithTranslations* chain, int nbrSpin, int momentum, double jFactor, double phiJ, double fFactor, double phiF, double boundaryCondition);

  // destructor
  //
  ~Potts3ChainHamiltonianWithTranslations();

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

 
 private:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

};

#endif
