////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                         Class author: Cecile Repellin                      //
//                                                                            //
//                                                                            //
//         class of two dimensional spin model on the kagome lattice          //
//                                                                            //
//                        last modification : 18/06/2016                      //
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


#ifndef TWODIMENSIONALKAGOMELATTICEHAMILTONIAN_H
#define TWODIMENSIONALKAGOMELATTICEHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Hamiltonian/TwoDimensionalTransverseFieldIsingHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class TwoDimensionalKagomeLatticeHamiltonian : public TwoDimensionalTransverseFieldIsingHamiltonian
{

 protected:
  
  // amplitude of the SzSz term for downward triangles
  double JDownFactor;
  // amplitude of the SxSx + SySy term
  double JEasyPlaneFactor;
  // amplitude of the SxSx + SySy term
  double JDownEasyPlaneFactor;
  //amplitude of the D3 breaking term
  double JBreakD3Factor;
  
   // flag for implementation of hermitian symmetry
  bool HermitianSymmetryFlag;
  //offset for tilted lattice
  int Offset;
  
 public:

   // default constructor
   //
   TwoDimensionalKagomeLatticeHamiltonian();
   
  // constructor from default data
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpinX = number of spin along the x direction
  // nbrSpinY = number of spin along the y direction
  // jFactor = amplitude of the Ising term
  // hxFactor = amplitudes of the Zeeman term along x
  // hzFactor = amplitudes of the Zeeman term along z
  // periodicBoundaryConditions = true if periodic boundary conditions have to be used
  TwoDimensionalKagomeLatticeHamiltonian(AbstractSpinChain* chain, int nbrSpinX, int nbrSpinY, double jFactor, double jDownFactor, double jEasyPlaneFactor, double jDownEasyPlaneFactor, bool periodicBoundaryConditions = true, int offset = 0);

  // destructor
  //
  ~TwoDimensionalKagomeLatticeHamiltonian();

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
  
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
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
  
  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* HermitianLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						  int firstComponent, int nbrComponent);

     // ask if Hamiltonian implements hermitian symmetry operations
  //
  virtual bool IsHermitian();
  
 protected:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

  // get a linearized position index from the 2d coordinates
  //
  // xPosition = position along the x direction
  // yPosition = position along the y direction
  // atomicIndex = position in the unit cell
  // return value = linearized index
  virtual int GetLinearizedIndex(int xPosition, int yPosition, int atomicIndex);

};

// get a linearized position index from the 2d coordinates
//
// xPosition = position along the x direction
// yPosition = position along the y direction
// return value = linearized index

inline int TwoDimensionalKagomeLatticeHamiltonian::GetLinearizedIndex(int xPosition, int yPosition, int atomicIndex)
{
  if (xPosition < 0)
    xPosition += this->NbrSpinX;
  if (xPosition >= this->NbrSpinX)
    xPosition -= this->NbrSpinX;
  if (yPosition < 0)
    yPosition += this->NbrSpinY;
  return (3 * ((xPosition * this->NbrSpinY) + yPosition) + atomicIndex);
}

#endif
