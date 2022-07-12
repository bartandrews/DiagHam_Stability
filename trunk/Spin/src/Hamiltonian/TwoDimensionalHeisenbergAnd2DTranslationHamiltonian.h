////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of two dimension Heisenberg model                  //
//                             and 2d translations                            //
//                                                                            //
//                        last modification : 13/05/2018                      //
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


#ifndef TWODIMENSIONALHEISENBERGAND2DTRANSLATIONISINGHAMILTONIAN_H
#define TWODIMENSIONALHEISENBERGAND2DTRANSLATIONISINGHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Architecture/AbstractArchitecture.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class TwoDimensionalHeisenbergAnd2DTranslationHamiltonian : public AbstractHamiltonian
{

 protected:
  
  // architecture used for precalculation
  AbstractArchitecture* Architecture;

  //pointer to Hilbert space of the associated system
  AbstractSpinChain* Chain;

  // total number of spins
  int NbrSpin;
  // number of spin chain along the x direction
  int NbrSpinX;
  // number of spin chain along the y direction
  int NbrSpinY;

  // amplitude of the Heisenberg XX coupling between nearest neighbors
  double JFactor;
  // amplitude of the Heisenberg Z coupling between nearest neighbors
  double JzFactor;
  // half of the Heisenberg XX coupling half of the between nearest neighbors
  double HalfJFactor;

  // array to store the diagonal contribution of the Hamiltonian
  double* SzSzContributions;

  // momentum along the x direction
  int XMomentum;
  // momentum along the y direction
  int YMomentum;

  //array containing all the phase factors that are needed when computing matrix elements
  Complex** ExponentialFactors;

  // shift to apply to go from precalculation index to the corresponding index in the HilbertSpace
  int PrecalculationShift;

  // flag for fast multiplication algorithm
  bool FastMultiplicationFlag;
  // step between each precalculated index
  int FastMultiplicationStep;

  // stored interactions per component
  int *NbrInteractionPerComponent;

  // number of tasks for load balancing
  int NbrBalancedTasks;
  // load balancing array for parallelisation, indicating starting indices
  long* LoadBalancingArray;

  // indices of matrix elements per component
  int **InteractionPerComponentIndex;
  // coefficients of matrix elements per component
  Complex** InteractionPerComponentCoefficient;

  // flag for implementation of hermitian symmetry
  bool HermitianSymmetryFlag;
  
 public:

  // default constructor
  //
  TwoDimensionalHeisenbergAnd2DTranslationHamiltonian();

  // constructor from default data
  //
  // chain = pointer to Hilbert space of the associated system
  // xMomentum = momentum along the x direction
  // nbrSpinX = number of spin along the x direction
  // yMomentum = momentum along the y direction
  // nbrSpinY = number of spin along the y direction
  // jFactor = Heisenberg XX coupling constant between nearest neighbors
  // jzFactor = Heisenberg Z coupling constant between nearest neighbors
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  TwoDimensionalHeisenbergAnd2DTranslationHamiltonian(AbstractSpinChain* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, 
						      double jFactor, double jzFactor, AbstractArchitecture* architecture, long memory = -1l);

  // destructor
  //
  ~TwoDimensionalHeisenbergAnd2DTranslationHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // ask if Hamiltonian implements hermitian symmetry operations
  //
  virtual bool IsHermitian();

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
 
  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							      int firstComponent, int nbrComponent);

 protected:
 
  // evaluate all matrix elements
  //   
  virtual void EvaluateDiagonalMatrixElements();

  // evaluate all exponential factors
  //   
  virtual void EvaluateExponentialFactors();

  // get a linearized position index from the 2d coordinates
  //
  // xPosition = unit cell position along the x direction
  // yPosition = unit cell position along the y direction
  // return value = linearized index
  virtual int GetLinearizedIndex(int xPosition, int yPosition);

  // get a linearized position index from the 2d coordinates, wihtout assuming that the input parametes are lower than the maximum ones
  //
  // xPosition = unit cell position along the x direction
  // yPosition = unit cell position along the y direction
  // return value = linearized index
  virtual int GetSafeLinearizedIndex(int xPosition, int yPosition);

  // evaluate all interaction factors
  //   
  //  virtual void EvaluateInteractionFactors();

  // test the amount of memory needed for fast multiplication algorithm
  //
  // allowedMemory = amount of memory that cam be allocated for fast multiplication
  // return value = amount of memory needed
  virtual long FastMultiplicationMemory(long allowedMemory);

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // return value = number of non-zero matrix element
  virtual long PartialFastMultiplicationMemory(int firstComponent, int lastComponent);

  // enable fast multiplication algorithm
  //
  virtual void EnableFastMultiplication();

  // enable fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // nbrComponent  = index of the last component that has to be precalcualted
  virtual void PartialEnableFastMultiplication(int firstComponent, int nbrComponent);

  // core part of the FastMultiplication method
  // 
  // chain = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray  
  virtual void EvaluateFastMultiplicationComponent(AbstractSpinChain* chain, int index, 
						   int* indexArray, Complex* coefficientArray, long& position);

  // core part of the PartialFastMultiplicationMemory
  // 
  // chain = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations  
  virtual void EvaluateFastMultiplicationMemoryComponent(AbstractSpinChain* chain, int firstComponent, int lastComponent, long& memory);

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

};

// get a linearized position index from the 2d coordinates
//
// xPosition = unit cell position along the x direction
// yPosition = unit cell position along the y direction
// index = site index within the unit cell
// return value = linearized index

inline int TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::GetLinearizedIndex(int xPosition, int yPosition)
{
  return ((xPosition * this->NbrSpinY) + yPosition);
}

// get a linearized position index from the 2d coordinates, wihtout assuming that the input parametes are lower than the maximum ones
//
// xPosition = unit cell position along the x direction
// yPosition = unit cell position along the y direction
// return value = linearized index

inline int TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::GetSafeLinearizedIndex(int xPosition, int yPosition)
{
  xPosition %= this->NbrSpinX;
  yPosition %= this->NbrSpinY;
  return ((xPosition * this->NbrSpinY) + yPosition);
}

#endif
