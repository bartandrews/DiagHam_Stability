////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of two dimension Shastry Sutherland model              //
//                                                                            //
//                        last modification : 07/06/2017                      //
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


#ifndef SHASTRYSUTHERLANDHAMILTONIAN_H
#define SHASTRYSUTHERLANDHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Hamiltonian/AbstractHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class ShastrySutherlandHamiltonian : public AbstractHamiltonian
{

 protected:
  
  //pointer to Hilbert space of the associated system
  AbstractSpinChain* Chain;

  // total number of spins
  int NbrSpin;
  // number of spin chain along the x direction
  int NbrSpinX;
  // number of spin chain along the y direction
  int NbrSpinY;

  // amplitude of the nearest neighbour Heisenberg term
  double JFactor;
  // amplitude of the dimer Heisenberg term
  double JpFactor;

  // array to store the diagonal contribution of the Hamiltonian
  double* SzSzContributions;

 public:
   
   // default constructor
  //
  ShastrySutherlandHamiltonian();

  // constructor from default data
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpinX = number of spin along the x direction
  // nbrSpinY = number of spin along the y direction
  // jFactor = amplitude of the nearest neighbour Heisenberg term
  // jpFactor = amplitude of the dimer Heisenberg term
  ShastrySutherlandHamiltonian(AbstractSpinChain* chain, int nbrSpinX, int nbrSpinY, double jFactor, double jpFactor);

  // destructor
  //
  ~ShastrySutherlandHamiltonian();

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

  // get a linearized position index from the 2d coordinates
  //
  // xPosition = position along the x direction
  // yPosition = position along the y direction
  // return value = linearized index
  virtual int GetLinearizedIndex(int xPosition, int yPosition);

  // get a linearized position index from the 2d coordinates, wihtout assuming that the input parametes are lower than the maximum ones
  //
  // xPosition = position along the x direction
  // yPosition = position along the y direction
  // return value = linearized index
  virtual int GetSafeLinearizedIndex(int xPosition, int yPosition);

};

// get a linearized position index from the 2d coordinates
//
// xPosition = position along the x direction
// yPosition = position along the y direction
// return value = linearized index

inline int ShastrySutherlandHamiltonian::GetLinearizedIndex(int xPosition, int yPosition)
{
  return ((xPosition * this->NbrSpinY) + yPosition);
}

// get a linearized position index from the 2d coordinates, wihtout assuming that the input parametes are lower than the maximum ones
//
// xPosition = position along the x direction
// yPosition = position along the y direction
// return value = linearized index

inline int ShastrySutherlandHamiltonian::GetSafeLinearizedIndex(int xPosition, int yPosition)
{
  xPosition %= this->NbrSpinX;
  yPosition %= this->NbrSpinY;
  return ((xPosition * this->NbrSpinY) + yPosition);
}

#endif
