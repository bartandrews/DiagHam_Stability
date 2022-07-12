////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                         Class author: Cecile Repellin                      //
//                                                                            //
//                                                                            //
//       class of two dimensional spin model on the triangular lattice        //
//                   with pseudospin and translations                         //
//                                                                            //
//                     last modification : 08/12/2016                         //
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


#ifndef TWODIMENSIONALTRIANGULARLATTICEWITHPSEUDOSPINAND2DTRANSLATIONHAMILTONIAN_H
#define TWODIMENSIONALTRIANGULARLATTICEWITHPSEUDOSPINAND2DTRANSLATIONHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/Spin1_2ChainWithPseudospinAnd2DTranslation.h"
#include "Hamiltonian/TwoDimensionalTriangularLatticeWithPseudospinHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian : public TwoDimensionalTriangularLatticeWithPseudospinHamiltonian
{

 protected:
   
  // momentum along the x direction
  int XMomentum;
  // momentum along the y direction
  int YMomentum;

  //array containing all the phase factors that are needed when computing matrix elements
  Complex** ExponentialFactors;
  
//   // shift to apply to go from precalculation index to the corresponding index in the HilbertSpace
//   int PrecalculationShift;
//   // amount of memory (in bytes) that can be used to store precalculated matrix elements
//   long Memory;
//   // flag for fast multiplication algorithm
//   bool FastMultiplicationFlag;
//   // step between each precalculated index (main part: start at 0, FastMultiplicationStep, 2*FastMultiplicationStep, ...)
//   int FastMultiplicationStep;
//   // number of non-null term in the hamiltonian for each state
//   int* NbrInteractionPerComponent;
//   // index of the state obtained for each term of the hamiltonian when applying on a given state
//   int** InteractionPerComponentIndex;
//   // multiplicative coefficient obtained for each term of the hamiltonian when applying on a given state and with a given destination state
//   Complex** InteractionPerComponentCoefficient;
  
 public:

  // constructor
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpinX = number of spin along the x direction
  // nbrSpinY = number of spin along the y direction
  // jFactor = amplitude of the Ising term
  // hxFactor = amplitudes of the Zeeman term along x
  // hzFactor = amplitudes of the Zeeman term along z
  // periodicBoundaryConditions = true if periodic boundary conditions have to be used
  TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian (Spin1_2ChainWithPseudospin* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, double jFactor, bool periodicBoundaryConditions, int offset, long memory = -1);
  
  // constructor
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpinX = number of spin along the x direction
  // nbrSpinY = number of spin along the y direction
  // jFactor = amplitude of the Ising term
  // hxFactor = amplitudes of the Zeeman term along x
  // hzFactor = amplitudes of the Zeeman term along z
  // periodicBoundaryConditions = true if periodic boundary conditions have to be used
  TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian (Spin1_2ChainWithPseudospin* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, double jFactor, double jBreakD3Factor, bool periodicBoundaryConditions, int offset, long memory = -1);

  // destructor
  //
  ~TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian();

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
  
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
  
  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
//   virtual ComplexVector* HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
// 						  int firstComponent, int nbrComponent);
  
//    // multiply a vector by the current hamiltonian for a given range of indices 
//   // and add result to another vector, low level function (no architecture optimization)
//   // using partial fast multiply option
//   //
//   // vSource = vector to be multiplied
//   // vDestination = vector at which result has to be added
//   // firstComponent = index of the first component to evaluate
//   // nbrComponent = number of components to evaluate
//   // return value = reference on vector where result has been stored
//   virtual ComplexVector& HermitianLowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
// 							     int firstComponent, int nbrComponent);
//   
//   // multiply a et of vectors by the current hamiltonian for a given range of indices 
//   // and add result to another et of vectors, low level function (no architecture optimization)
//   // using partial fast multiply option
//   //
//   // vSources = array of vectors to be multiplied
//   // vDestinations = array of vectors at which result has to be added
//   // nbrVectors = number of vectors that have to be evaluated together
//   // firstComponent = index of the first component to evaluate
//   // nbrComponent = number of components to evaluate
//   // return value = pointer to the array of vectors where result has been stored
//   virtual ComplexVector* HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
//                                                                      int firstComponent, int nbrComponent);
//   
//   // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//   //
//   // firstComponent = index of the first component that has to be precalcualted
//   // nbrComponent  = number of components that has to be precalcualted
//   // return value = number of non-zero matrix element
//   virtual long PartialFastMultiplicationMemory(int firstComponent, int nbrComponent);
// 
//   // firstComponent = index of the first component that has to be precalcualted
//   // nbrComponent  = number of components that has to be precalcualted
//   virtual void PartialEnableFastMultiplication(int firstComponent, int nbrComponent);

  
  // test the amount of memory needed for fast multiplication algorithm
  //
  // allowedMemory = amount of memory that cam be allocated for fast multiplication
  // return value = amount of memory needed
//   long FastMultiplicationMemory(long allowedMemory);
  
  
  // enable fast multiplication algorithm
  //
//   void EnableFastMultiplication();
  
  
  // ask if Hamiltonian implements hermitian symmetry operations
  //
  virtual bool IsHermitian();
  
 protected:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();
  
  // evaluate all exponential factors
  //   
  virtual void EvaluateExponentialFactors();

  


};

#endif
