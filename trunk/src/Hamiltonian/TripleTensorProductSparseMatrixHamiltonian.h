////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of hamiltonian defined as a linear combination            //
//                           of triple tensor products                        //
//                                                                            //
//                        last modification : 31/07/2016                      //
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


#ifndef TRIPLETENSORPRODUCTSPARSEMATRIXHAMILTONIAN_H
#define TRIPLETENSORPRODUCTSPARSEMATRIXHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Hamiltonian/TensorProductSparseMatrixHamiltonian.h"
#include "Matrix/SparseRealMatrix.h"


using std::ostream;


class TripleTensorProductSparseMatrixHamiltonian : public TensorProductSparseMatrixHamiltonian
{

  friend class VectorTensorMultiplicationCoreOperation;
  friend class VectorSparseTensorMultiplyOperation;

 protected:

  // middle matrices of each tensor product
  SparseRealMatrix* MiddleMatrices;  

  // number of row for the middle matrix
  int MiddleMatrixNbrRow;
  
 public:

  // default contructor 
  //
  TripleTensorProductSparseMatrixHamiltonian();

  // contructor 
  //
  // nbrTensorProducts = number of tensor products whose linear combination defined the Hamiltonian 
  // leftMatrices = left matrices of each tensor product
  // middleMatrices = middle matrices of each tensor product
  // rightMatrices = right matrices of each tensor product
  // coefficients = coefficients of the ensor product linear combination
  // architecture = architecture to use for precalculation
  TripleTensorProductSparseMatrixHamiltonian(int nbrTensorProducts, SparseRealMatrix* leftMatrices, SparseRealMatrix* middleMatrices, SparseRealMatrix* rightMatrices, 
					     double* coefficients, AbstractArchitecture* architecture);

  // destructor
  //
  ~TripleTensorProductSparseMatrixHamiltonian();

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
  

  // get the linearized index corresponding to a set of indices
  //
  // leftIndex = left index 
  // middleIndex = middle index 
  // rightIndex = right index 
  // return value = linearized index
  virtual int GetLinearizedIndex (int leftIndex, int middleIndex, int rightIndex);

  // get a set of indices from their linearized version
  //
  // linearizedIndex = linearized index
  // leftIndex = reference ont the left index 
  // middleIndex = reference ont the middle index 
  // rightIndex = reference ont the right index 
  virtual void GetIndicesFromLinearizedIndex (int linearizedIndex, int& leftIndex, int& middleIndex, int& rightIndex);

 protected:
  
  // initialize the temporary arrays
  //
  void InitializeTemporaryArrays();

  // core part of the tensor-multiplication
  //
  // tensorIndex = index of tensore to consider
  // localTemporaryArray = temporary array used to store the partial multiplication
  // vSource = vector to be multiplied
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  virtual void LowLevelAddMultiplyTensorCore(int tensorIndex, double** localTemporaryArray,
					     RealVector& vSource, int firstComponent, int nbrComponent);

  // core part of the tensor-multiplication
  //
  // tensorIndex = index of tensore to consider
  // localTemporaryArray = temporary array used to store the partial multiplication
  // vSource = vector to be multiplied
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  virtual void LowLevelAddMultiplyTensorCore(int tensorIndex, Complex** localTemporaryArray,
					     ComplexVector& vSource, int firstComponent, int nbrComponent);

  // core part of the tensor-multiplication (second part computing the final result for one tensor product)
  //
  // tensorIndex = index of tensore to consider
  // localTemporaryArray = temporary array used to store the partial multiplication
  // vDestination = vector where the result will be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  virtual void LowLevelAddMultiplyTensorCoreDestination(int tensorIndex, double** localTemporaryArray, RealVector& vDestination, 
							int firstComponent, int nbrComponent);

  // core part of the tensor-multiplication (second part computing the final result for one tensor product)
  //
  // tensorIndex = index of tensore to consider
  // localTemporaryArray = temporary array used to store the partial multiplication
  // vDestination = vector where the result will be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  virtual void LowLevelAddMultiplyTensorCoreDestination(int tensorIndex, Complex** localTemporaryArray, ComplexVector& vDestination, 
							int firstComponent, int nbrComponent);

};

// get the linearized index corresponding to a set of indices
//
// leftIndex = left index 
// middleIndex = middle index 
// rightIndex = right index 
// return value = linearized index

inline int TripleTensorProductSparseMatrixHamiltonian::GetLinearizedIndex (int leftIndex, int middleIndex, int rightIndex)
{
  return ((((leftIndex * this->MiddleMatrixNbrRow) + middleIndex) * this->RightMatrixNbrRow) + rightIndex);
}

// get a set of indices from their linearized version
//
// linearizedIndex = linearized index
// leftIndex = reference ont the left index 
// middleIndex = reference ont the middle index 
// rightIndex = reference ont the right index 

inline void TripleTensorProductSparseMatrixHamiltonian::GetIndicesFromLinearizedIndex (int linearizedIndex, int& leftIndex, int& middleIndex, int& rightIndex)
{
  rightIndex = linearizedIndex % this->RightMatrixNbrRow;
  linearizedIndex /= this->RightMatrixNbrRow;;
  leftIndex = linearizedIndex / this->MiddleMatrixNbrRow;
  middleIndex =  linearizedIndex % this->MiddleMatrixNbrRow;
}

#endif
