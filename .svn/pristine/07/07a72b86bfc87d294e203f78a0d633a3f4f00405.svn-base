////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of hamiltonian defined as a linear combination of         //
//                           multiple tensor products                         //
//                                                                            //
//                        last modification : 30/07/2016                      //
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


#ifndef MULTIPLETENSORPRODUCTSPARSEMATRIXHAMILTONIAN_H
#define MULTIPLETENSORPRODUCTSPARSEMATRIXHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/SparseRealMatrix.h"


using std::ostream;
class MathematicaOutput;
class Matrix;


class MultipleTensorProductSparseMatrixHamiltonian : public AbstractHamiltonian
{

  friend class VectorTensorMultiplicationCoreOperation;
  friend class VectorSparseTensorMultiplyOperation;

 protected:

  // Hilbert space assocaited to the Hamiltonian 
  AbstractHilbertSpace* HilbertSpace;
  
  // number of tensor products whose linear combination defined the Hamiltonian 
  int NbrTensorProducts;
  // number of matrices within each multiple tensor products
  int NbrMatrixPerTensorProducts;
  // left matrices of each tensor product
  SparseRealMatrix** Matrices;  
  // coefficients of the ensor product linear combination
  double* Coefficients;

  // global shift to apply to the diagonal matrix elements
  double HamiltonianShift;

  // a temporary array used to perform the tensor-vector multiplication (thread safe)
  double** TemporaryArray;
  // a temporary array used to perform the tensor-vector multiplication (thread safe, complex version)
  Complex** ComplexTemporaryArray;
  
  // pointer to the architecture
  AbstractArchitecture* Architecture;

 public:

  // default contructor 
  //
  MultipleTensorProductSparseMatrixHamiltonian();

  // contructor 
  //
  // nbrTensorProducts = number of multiple tensor products whose linear combination defined the Hamiltonian 
  // nbrMatrixPerTensorProducts = number of matrices within each multiple tensor products
  // matrices = matrices of each multiple tensor product (first entry) and whithin each multiple tensor product (second entry)
  // coefficients = coefficients of the ensor product linear combination
  // architecture = architecture to use for precalculation
  MultipleTensorProductSparseMatrixHamiltonian(int nbrTensorProducts, int nbrMatrixPerTensorProducts, SparseRealMatrix** matrices, 
					       double* coefficients, AbstractArchitecture* architecture);

  // destructor
  //
  ~MultipleTensorProductSparseMatrixHamiltonian();

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
  virtual RealVector& ConjugateLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
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
  virtual RealVector* ConjugateLowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& ConjugateLowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
  virtual ComplexVector* ConjugateLowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							   int firstComponent, int nbrComponent);

  // ask if Hamiltonian implements methods applying the conjugate of the Hamiltonian
  //
  virtual bool IsConjugate();

  // test if the hamiltonian is compatible with the hamiltonian-vector multiplication operations
  //
  // return value = true if compatible (otherwise, any parallelization will be disable for these operations)
  virtual bool IsHamiltonianVectorOperationCompatible();
  
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

// ask if Hamiltonian implements methods applying the conjugate of the Hamiltonian
//

inline bool MultipleTensorProductSparseMatrixHamiltonian::IsConjugate()
{
  return true;
}

// test if the hamiltonian is compatible with the hamiltonian-vector multiplication operations
//
// return value = true if compatible (otherwise, any parallelization will be disable for these operations)
  
inline bool MultipleTensorProductSparseMatrixHamiltonian::IsHamiltonianVectorOperationCompatible()
{
  return false;
}
  
#endif
