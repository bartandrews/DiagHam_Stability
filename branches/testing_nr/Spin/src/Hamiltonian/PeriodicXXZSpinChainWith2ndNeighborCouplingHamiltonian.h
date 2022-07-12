////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of periodic XXZ spin chain hamiltonian with second            //
//    nearest neighbor coupling in z direction and alternative coupling to    //
//                              the magnetic field                            //
//                                                                            //
//                        last modification : 04/02/2003                      //
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


#ifndef PERIODICXXZSPINCHAINWITHSECONDNEARESTNEIGHBORCOUPLINGHAMILTONIAN_H
#define PERIODICXXZSPINCHAINWITHSECONDNEARESTNEIGHBORCOUPLINGHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Hamiltonian/AbstractHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian : public AbstractHamiltonian
{

 protected:
  
  // hilbert space describing the spin chain
  AbstractSpinChain* Chain;

  // coupling constant between nearest neighbor Sz components (in jx coupling constant unit)
  double Jz;
  // coupling constant between second nearest neighbor Sz components (in jx coupling constant unit)
  double Jz2;

  // magnetic field times the coupling constant (in jx coupling constant unit)
  double BField;
  
  // number of spins
  int NbrSpin;

  // arary containing the hamiltonian diagonal elements
  double* SzSzContributions;

 public:

  // constructor from default datas
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpin = number of spin (rounded to the lower even number)
  // jz = coupling constant between nearest neighbor Sz components (in jx coupling constant unit)
  // jz2 = coupling constant between second nearest neighbor Sz components (in jx coupling constant unit)
  // bField = magnetic field times the coupling constant (in jx coupling constant unit)
  PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian(AbstractSpinChain* chain, int nbrSpin, double jz, double jz2,
								   double bField);

  // destructor
  //
  ~PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian();

  // set chain
  // 
  // chain = pointer on Hilbert space of the associated system
  // return value = reference on current Hamiltonian
  PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian& SetChain(AbstractSpinChain* chain);

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

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of idinces 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
			       int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
				  int firstComponent, int nbrComponent);

  // return a list of left interaction operators
  //
  // return value = list of left interaction operators
  List<Matrix*> LeftInteractionOperators();  

  // return a list of right interaction operators 
  //
  // return value = list of right interaction operators
  List<Matrix*> RightInteractionOperators();  

  // Output Stream overload
  //
  // Str = reference on output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian& H);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian& H);

 private:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

};

#endif
