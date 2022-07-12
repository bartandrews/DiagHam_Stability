////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//  class of hamiltonian defined as a linear combination of tensor products   //
//               focusing on a single block of  tensor product                //
//                                                                            //
//                        last modification : 07/01/2013                      //
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


#ifndef TENSORPRODUCTSPARSEMATRIXSELECTEDBLOCKHAMILTONIAN_H
#define TENSORPRODUCTSPARSEMATRIXSELECTEDBLOCKHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Hamiltonian/TensorProductSparseMatrixHamiltonian.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/RealMatrix.h"


using std::ostream;

class MathematicaOutput;
class Matrix;
class AbstractArchitecture;
class AbstractArchitecture;


class TensorProductSparseMatrixSelectedBlockHamiltonian : public TensorProductSparseMatrixHamiltonian
{

  friend class VectorSparseTensorMultiplyOperation;

 protected:

  // linearized indices that define the selected block 
  long* BlockIndices;

  int* InvertBlockIndices;

  // selected block size
  int BlockSize;

  // table that contains all the linearized indices attached to one left matrix index
  long** BlockIndexProductTable;
  // number of  linearized indices attached to one left matrix index
  int* BlockIndexProductTableNbrElements;
  // first index in the hilbert space where a given left matrix index occurs
  int* BlockIndexProductTableShift;
  // true if the BlockIndexProductTable* arrays have been provided to the hamiltonian and not computed locally
  bool ExternalBlockIndexProductTable;

  // flag for fast multiplication algorithm
  bool FastMultiplicationFlag;
  
  //pointer to the architecture
  AbstractArchitecture* Architecture;

  // temporary arrays used to store the full tensor product sparse matrix
  long* TemporaryRowPointers;
  long* TemporaryRowLastPointers;
  double* TemporaryMatrixElements;
  int* TemporaryMatrixColumnIndices;
  // effective dimension of TemporaryRowPointers and TemporaryRowLastPointers
  int EffectiveHilbertSpaceDimension;

  // shift to apply to go from precalculation index to the corresponding index in the HilbertSpace
  int PrecalculationShift;

 public:

  // contructor 
  //
  // nbrTensorProducts = number of tensor products whose linear combination defined the Hamiltonian 
  // leftMatrices = left matrices of each tensor product
  // rightMatrices = right matrices of each tensor product
  // coefficients = coefficients of the ensor product linear combination
  // blockSize = number of indices in the selected block
  // blockIndices = pairs of indices (for resp. the left and right matrix) that define the selected block 
  // architecture = architecture to use for precalculation
  // memory = amount of memory that can be used for precalculations (in bytes)
  TensorProductSparseMatrixSelectedBlockHamiltonian(int nbrTensorProducts, SparseRealMatrix* leftMatrices,  SparseRealMatrix* rightMatrices, double* coefficients,
						    int blockSize, long* blockIndices, AbstractArchitecture* architecture, long memory);

  // contructor providing an efficient block index scheme
  //
  // nbrTensorProducts = number of tensor products whose linear combination defined the Hamiltonian 
  // leftMatrices = left matrices of each tensor product
  // rightMatrices = right matrices of each tensor product
  // coefficients = coefficients of the ensor product linear combination
  // blockSize = number of indices in the selected block
  // blockIndices = pairs of indices (for resp. the left and right matrix) that define the selected block 
  // blockIndexProductTable = table that contains all the linearized indices attached to one left matrix index
  // blockIndexProductTableNbrElements = number of  linearized indices attached to one left matrix index
  // blockIndexProductTableShift = first index in the hilbert space where a given left matrix index occurs
  // architecture = architecture to use for precalculation
  // memory = amount of memory that can be used for precalculations (in bytes)  
  TensorProductSparseMatrixSelectedBlockHamiltonian(int nbrTensorProducts, SparseRealMatrix* leftMatrices,  
						    SparseRealMatrix* rightMatrices, double* coefficients,
						    int blockSize, long* blockIndices, 
						    long** blockIndexProductTable, int* blockIndexProductTableNbrElements,
						    int* blockIndexProductTableShift,
						    AbstractArchitecture* architecture, long memory);
  // destructor
  //
  ~TensorProductSparseMatrixSelectedBlockHamiltonian();

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

 protected:

  // test the amount of memory needed for fast multiplication algorithm
  //
  // return value = amount of memory needed
  virtual long FastMultiplicationMemory();

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // nbrComponent  = number of components that has to be precalcualted
  // return value = number of non-zero matrix elements that have to be stored
  virtual long PartialFastMultiplicationMemory(int firstComponent, int nbrComponent);

  // enable fast multiplication algorithm
  //
  virtual void EnableFastMultiplication();

  // enable fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // nbrComponent  = number of components that has to be precalcualted
  virtual void PartialEnableFastMultiplication(int firstComponent, int nbrComponent);
  
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

#endif
