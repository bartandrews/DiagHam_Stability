////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of two dimension spin model that could host             //
//                a Read-Rezayi Z3 phase with 2d translations                 //
//                                                                            //
//                        last modification : 27/07/2018                      //
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


#ifndef TWODIMENSIONALRRAND2DTRANSLATIONISINGHAMILTONIAN_H
#define TWODIMENSIONALRRAND2DTRANSLATIONISINGHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Hamiltonian/TwoDimensionalHeisenbergAnd2DTranslationHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class TwoDimensionalRRAnd2DTranslationHamiltonian : public TwoDimensionalHeisenbergAnd2DTranslationHamiltonian
{

 protected:
  
  // amplitude of the Heisenberg coupling between nearest neighbors
  double J1Factor;
  // amplitude of the (S_i S_j)^2 nearest neighbor coupling
  double J2Factor;
  // amplitude of the (S_i S_j)^3 nearest neighbor coupling
  double J3Factor;
  // amplitude of the chiral term
  double JcFactor;
  // half the amplitude of the chiral term
  double HalfJcFactor;

 public:

  // default constructor
  //
  TwoDimensionalRRAnd2DTranslationHamiltonian();

  // constructor from default data
  //
  // chain = pointer to Hilbert space of the associated system
  // xMomentum = momentum along the x direction
  // nbrSpinX = number of spin along the x direction
  // yMomentum = momentum along the y direction
  // nbrSpinY = number of spin along the y direction
  // j1Factor = amplitude of the Heisenberg coupling between nearest neighbors
  // j2Factor = amplitude of the (S_i S_j)^2 nearest neighbor coupling
  // j3Factor = amplitude of the (S_i S_j)^3 nearest neighbor coupling
  // jcFactor = amplitude of the chiral term
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  TwoDimensionalRRAnd2DTranslationHamiltonian(AbstractSpinChain* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, 
					      double j1Factor, double j2Factor, double j3Factor, double jcFactor, 
					      AbstractArchitecture* architecture, long memory = -1l);

  // destructor
  //
  ~TwoDimensionalRRAnd2DTranslationHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

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

  // core part of the AddMultiply method
  // 
  // chain = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  virtual void HermitianEvaluateAddMultiplyComponent(AbstractSpinChain* chain, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method for a set of vectors
  // 
  // chain = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void HermitianEvaluateAddMultiplyComponent(AbstractSpinChain* chain, int index, ComplexVector* vSources, 
						     ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // evaluate the off-diagonal contribution for one type of Hamiltonian terms ( all (S_i S_j)^n )
  //
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // index = index of the many-body state to act on
  // dimension = total Hilbert space dimension
  // vDestination = vector at which result has to be added
  // coefficient = global multiplicative coefficient
  virtual void EvaluateOffDiagonalPowerHeisenbergContribution(int i, int j, int index, int dimension, ComplexVector& vDestination, Complex& coefficient);

  // evaluate the off-diagonal contribution for one type of Hamiltonian terms ( all (S_i S_j)^n )
  //
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // index = index of the many-body state to act on
  // dimension = total Hilbert space dimension
  // vDestinations = vectors to which results have to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // coefficients = global multiplicative coefficients
  virtual void EvaluateOffDiagonalPowerHeisenbergContribution(int i, int j, int index, int dimension, ComplexVector* vDestinations, int nbrVectors, Complex* coefficients);

  // evaluate the off-diagonal chiral contribution for a single term ( S_i (S_j ^ S_k) )
  //
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // k = linearized position of the second spin
  // index = index of the many-body state to act on
  // dimension = total Hilbert space dimension
  // vDestination = vector at which result has to be added
  // coefficient = global multiplicative coefficient
  virtual void EvaluateOffDiagonalChiralContribution(int i, int j, int k, int index, int dimension, ComplexVector& vDestination, Complex& coefficient);

  // evaluate the off-diagonal chiral contribution for a single term ( S_i (S_j ^ S_k) )
  //
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // k = linearized position of the second spin
  // index = index of the many-body state to act on
  // dimension = total Hilbert space dimension
  // vDestinations = vectors to which results have to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // coefficients = global multiplicative coefficients
  virtual void EvaluateOffDiagonalChiralContribution(int i, int j, int k, int index, int dimension, ComplexVector* vDestinations, int nbrVectors, Complex* coefficients);

  // evaluate the off-diagonal contribution for one type of Hamiltonian terms ( all (S_i S_j)^n )
  //
  // chain = pointer to the Hilbert space
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // index = index of the many-body state to act on
  // vSource = vector  to be multiplied
  // vDestination = vector to which result has to be added
  virtual void HermitianEvaluateOffDiagonalPowerHeisenbergContribution(AbstractSpinChain* chain, int i, int j, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // evaluate the off-diagonal chiral contribution for a single term ( S_i (S_j ^ S_k) )
  //
  // chain = pointer to the Hilbert space
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // k = linearized position of the second spin
  // index = index of the many-body state to act on
  // vSource = vector  to be multiplied
  // vDestination = vector to which result has to be added
  virtual void HermitianEvaluateOffDiagonalChiralContribution(AbstractSpinChain* chain, int i, int j, int k, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the FastMultiplication method
  // 
  // chain = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray  
  virtual void EvaluateFastMultiplicationComponent(AbstractSpinChain* chain, int index, 
						   int* indexArray, Complex* coefficientArray, long& position);

  // core part of the PartialFastMultiplication for all Hamiltonian terms (S_i S_j)^n 
  //
  // chain = pointer to the Hilbert space
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // index = index of the many-body state to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray  
  virtual void CoreEvaluateFastMultiplicationComponentPowerHeisenberg(AbstractSpinChain* chain, int i, int j, int index, 
									    int* indexArray, Complex* coefficientArray, long& position);

  //  core part of the PartialFastMultiplication for a single term ( S_i (S_j ^ S_k) )
  //
  // chain = pointer to the Hilbert space
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // k = linearized position of the second spin
  // index = index of the many-body state to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray  
  virtual void CoreEvaluateFastMultiplicationComponentChiral(AbstractSpinChain* chain, int i, int j, int k, int index, 
								   int* indexArray, Complex* coefficientArray, long& position);

  // core part of the PartialFastMultiplicationMemory
  // 
  // chain = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations  
  virtual void EvaluateFastMultiplicationMemoryComponent(AbstractSpinChain* chain, int firstComponent, int lastComponent, long& memory);

  // core part of the PartialFastMultiplicationMemory for all Hamiltonian terms (S_i S_j)^n 
  //
  // chain = pointer to the Hilbert space
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // index = index of the many-body state to act on
  // nbrInteractionPerComponent = array that contains the number of interaction per component
  // memory = reference on the amount of memory required for precalculations  
  virtual void CoreEvaluateFastMultiplicationMemoryComponentPowerHeisenberg(AbstractSpinChain* chain, int i, int j, int index, int* nbrInteractionPerComponent, long& memory);             

  //  core part of the PartialFastMultiplicationMemory for a single term ( S_i (S_j ^ S_k) )
  //
  // chain = pointer to the Hilbert space
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // k = linearized position of the second spin
  // index = index of the many-body state to act on
  // nbrInteractionPerComponent = array that contains the number of interaction per component
  // memory = reference on the amount of memory required for precalculations  
  virtual void CoreEvaluateFastMultiplicationMemoryComponentChiral(AbstractSpinChain* chain, int i, int j, int k, int index, int* nbrInteractionPerComponent, long& memory);

};

#endif
