////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of Z2 interacting chain with translations               //
//                                                                            //
//                        last modification : 06/06/2014                      //
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


#ifndef SPINCHAINZ2INTERACTINGHAMILTONIANWITHTRANSLATIONS_H
#define SPINCHAINZ2INTERACTINGHAMILTONIANWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class SpinChainZ2InteractingHamiltonianWithTranslations : public AbstractHamiltonian
{

 protected:
  
  // pointer to the Hilbert space of the system
  Spin1_2ChainWithTranslations* Chain;

  // coupling along the z direction
  double JFactor;
  // Zeeman term
  double FFactor;
  // coupling along the x direction
  double InteractionStrength;

  // boundary conditions (0 for open chain, 1 for periodic, -1 for antiperiodic)
  double BoundaryCondition;

  // number of spins
  int NbrSpin;

  // precalculation array where the diagonal elements are stored
  double* SzSzContributions;
  // precalculation array where the parity of each state is stored
  double* Parities;

  //array containing all the complex phase that are needed when computing matrix elements
  Complex* ExponentialTable;

 public:

  // default constructor
  //
  SpinChainZ2InteractingHamiltonianWithTranslations();

  // constructor
  //
  // chain = pointer to the Hilbert space of the system
  // nbrSpin = number of spins
  // jFactor = coupling along the z direction
  // fFactor = Zeeman term
  // interactionStrength = coupling along the x direction
  // boundaryCondition = boundary condition to apply (1 for periodic, -1 for antiperiodic)
  SpinChainZ2InteractingHamiltonianWithTranslations(Spin1_2ChainWithTranslations* chain, int nbrSpin, 
						    double jFactor, double fFactor, double interactionStrength, 
						    double boundaryCondition = 1.0);

  // destructor
  //
  ~SpinChainZ2InteractingHamiltonianWithTranslations();

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


 protected:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

  // evaluate all cosine/sine that are needed when computing matrix elements
  //
  void EvaluateCosineTable();

};

#endif
