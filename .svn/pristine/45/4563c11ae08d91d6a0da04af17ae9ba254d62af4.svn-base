////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of spin chain hamiltonian with translations             //
//                                                                            //
//                        last modification : 04/03/2002                      //
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


#ifndef SPINCHAINHAMILTONIANWITHTRANSLATIONS_H
#define SPINCHAINHAMILTONIANWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChainWithTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class SpinChainHamiltonianWithTranslations : public AbstractHamiltonian
{

 protected:
  
  // architecture used for precalculation
  AbstractArchitecture* Architecture;

  // pointer to the Hilbert space
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
  Complex* ExponentialTable;

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
  SpinChainHamiltonianWithTranslations();

  // constructor from default datas
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpin = number of spin
  // j = coupling constants between spins
  // nnCoupling = term to add to ZZ nearest-neighbour interaction
  // nnnCoupling = nearest-neighbour interaction in Z direction
  
  SpinChainHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin, double j, double nnCoupling, double nnnCoupling);

  // destructor
  //
  ~SpinChainHamiltonianWithTranslations();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // ask if Hamiltonian implements hermitian symmetry operations
  //
  virtual bool IsHermitian();

  // set chain
  // 
  // chain = pointer on Hilbert space of the associated system
  // return value = reference on current Hamiltonian
  SpinChainHamiltonianWithTranslations& SetChain(AbstractSpinChainWithTranslations* chain);

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

  // save precalculations in a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be stored
  // return value = true if no error occurs
  bool SavePrecalculation (char* fileName);

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  virtual Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  virtual Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

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
  // and store result in another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						  int firstComponent, int nbrComponent);

  // return a list of left interaction operators
  //
  // return value = list of left interaction operators
  List<Matrix*> LeftInteractionOperators();  

  // return a list of right interaction operators 
  //
  // return value = list of right interaction operators
  List<Matrix*> RightInteractionOperators();  

  // get the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
  //
  // nbrThreads = number of threads requested
  // segmentIndices = array returning the reference to an array of the first index of each of the segments
  // return value = true if no error occured
  virtual bool GetLoadBalancing(int nbrTasks, long* &segmentIndices);
    
 protected:
 
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

  // evaluate all cosinus/sinus that are needed when computing matrix elements
  //
  virtual void EvaluateCosinusTable();

  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

};

// ask if Hamiltonian implements hermitian symmetry operations
//

inline bool SpinChainHamiltonianWithTranslations::IsHermitian()
{
  return this->HermitianSymmetryFlag;
}

#endif
