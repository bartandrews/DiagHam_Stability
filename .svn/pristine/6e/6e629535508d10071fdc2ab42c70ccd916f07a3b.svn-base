////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of abstract hamiltonian                      //
//                                                                            //
//                        last modification : 28/02/2001                      //
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


#ifndef ABSTRACTHAMILTONIAN_H
#define ABSTRACTHAMILTONIAN_H


#include "config.h"
#include "GeneralTools/List.h"
#include "GeneralTools/GarbageFlag.h"


class ComplexVector;
class RealVector;
class Vector;
class RealTriDiagonalSymmetricMatrix;
class RealSymmetricMatrix;
class HermitianMatrix;
class RealMatrix;
class ComplexMatrix;
class Complex;
class Matrix;
class AbstractHilbertSpace;
class AbstractBitmapPicture;
class SparseRealMatrix;
class SparseComplexMatrix;
class IntegerMatrix;
class LongIntegerMatrix;


class AbstractHamiltonian
{

  friend class GenericHamiltonianPrecalculationOperation;

 protected:
  
  GarbageFlag Flag;

  // flag to indicate if hamiltonian-vector multiplication is done on the left hand side
  bool LeftHamiltonianVectorMultiplicationFlag;

 public:

  // default constructor
  //
  AbstractHamiltonian();

  // destructor
  //
  virtual ~AbstractHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  //  virtual AbstractHamiltonian* Clone () = 0;

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  virtual void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace) = 0;

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  virtual AbstractHilbertSpace* GetHilbertSpace () = 0;

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  virtual int GetHilbertSpaceDimension () = 0;
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  virtual void ShiftHamiltonian (double shift) = 0;

  // save precalculations in a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be stored
  // return value = true if no error occurs
  virtual bool SavePrecalculation (char* fileName);
  
  // store Hamiltonian into an hermitian matrix
  //
  // M = reference on matrix where Hamiltonian has to be stored
  // return value = reference on  corresponding hermitian matrix
  virtual HermitianMatrix& GetHamiltonian (HermitianMatrix& M);
  
  // store Hamiltonian into a complex matrix
  //
  // M = reference on matrix where Hamiltonian has to be stored
  // return value = reference on  corresponding complex matrix
  virtual ComplexMatrix& GetHamiltonian (ComplexMatrix& M);
  
  // store Hamiltonian into a sparse complex sparse matrix
  //
  // M = reference on matrix where Hamiltonian has to be stored
  // return value = reference on  corresponding complex matrix
  virtual SparseComplexMatrix& GetHamiltonian (SparseComplexMatrix& M);
  
  // store real part of Hamiltonian into a real symmetric matrix
  //
  // M = reference on matrix where Hamiltonian has to be stored
  // return value = reference on  corresponding real symmetric matrix 
  virtual RealSymmetricMatrix& GetHamiltonian (RealSymmetricMatrix& M);
  
  // store real part of Hamiltonian into a real matrix
  //
  // M = reference on matrix where Hamiltonian has to be stored
  // return value = reference on  corresponding real matrix 
  virtual RealMatrix& GetHamiltonian (RealMatrix& M);

  // store real part of Hamiltonian into a real sparse matrix
  //
  // M = reference on matrix where Hamiltonian has to be stored
  // return value = reference on  corresponding real matrix 
  virtual SparseRealMatrix& GetHamiltonian (SparseRealMatrix& M);

  // store the real part of Hamiltonian into an integer matrix
  //
  // M = reference on matrix where Hamiltonian has to be stored
  // scalingFactor = use an additional scaling factor before converting coefficients into integers
  // return value = reference on corresponding matrix 
  virtual IntegerMatrix& GetHamiltonian (IntegerMatrix& M, double scalingFactor = 1.0);
  
  // store the real part of Hamiltonian into a long integer matrix
  //
  // M = reference on matrix where Hamiltonian has to be stored
  // scalingFactor = use an additional scaling factor before converting coefficients into integers
  // return value = reference on corresponding matrix 
  virtual LongIntegerMatrix& GetHamiltonian (LongIntegerMatrix& M, double scalingFactor = 1.0);
  
  // store real part of Hamiltonian into a matrix
  //
  // M = reference on matrix where Hamiltonian has to be stored
  // return value = reference on  corresponding matrix 
  virtual Matrix& GetHamiltonian (Matrix& M);
  
  // return matrix representation of current Hamiltonian
  //
  // return value = reference to representation
  virtual Matrix* GetHamiltonian ();
  
  // store Hamiltonian into a picture (drawing non zero element in black)
  //
  // error = absolute minimum value to be considered as non zero element
  // return value = pointer to the picture associated to the matrix
  virtual AbstractBitmapPicture* GetHamiltonianPicture (double error);

  // store Hamiltonian into a picture (drawing non zero element with a color scale)
  //
  // error = absolute minimum value to be considered as non zero element
  // return value = pointer to the picture associated to the matrix
  virtual AbstractBitmapPicture* GetHamiltonianColorPicture (double error);
    
  // return a list of left interaction operators
  //
  // return value = list of left interaction operators
  virtual List<Matrix*> LeftInteractionOperators();  

  // return a list of right interaction operators
  //
  // return value = list of right interaction operators
  virtual List<Matrix*> RightInteractionOperators();

  // get the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
  // nbrThreads = number of threads requested
  // segmentIndices = array returning the reference to an array of the first index of each of the segments
  //
  virtual bool GetLoadBalancing(int nbrTasks, long* &segmentIndices);

  // set the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
  // nbrThreads = number of threads requested
  // segmentIndices = array returning the first index of each of the segments
  //
  virtual bool SetLoadBalancing(int nbrTasks, long* segmentIndices);

  // ask if Hamiltonian implements methods using hermitian symmetry 
  //
  virtual bool IsHermitian();

  // ask if Hamiltonian implements methods applying the conjugate of the Hamiltonian
  //
  virtual bool IsConjugate();

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

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
				       int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
				       int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
				       int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination);

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
  // and add result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
					  int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
					  int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a set of vectors by the current hamiltonian and store result in another set of vectors
  // low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and store result in another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
					       int firstComponent, int nbrComponent);

  // multiply a set of vector by the current hamiltonian for a given range of indices 
  // and store result in another set of vector, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
					       int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
					       int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vector sat which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors);

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
  // and add result in another set of vectors, low level function (no architecture optimization)
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors,
						  int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						  int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  virtual ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					  int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					  int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
					  int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination);

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
  // vDestination = vector where result has to be stored
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
					     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a set of vectors by the current hamiltonian and store result in another set of vectors
  // low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors);

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

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and store result in another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						  int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						  int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);
  
  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors);

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
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  virtual RealVector& ConjugateLowLevelMultiply(RealVector& vSource, RealVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& ConjugateLowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
				       int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = reference on vector where result has been stored
  virtual RealVector& ConjugateLowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  virtual RealVector& ConjugateLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination);

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

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = reference on vector where result has been stored
  virtual RealVector& ConjugateLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						   int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						   int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a set of vectors by the current hamiltonian and store result in another set of vectors
  // low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* ConjugateLowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors);

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

  // multiply a set of vector by the current hamiltonian for a given range of indices 
  // and store result in another set of vector, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* ConjugateLowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vector sat which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* ConjugateLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* ConjugateLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							   int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result in another set of vectors, low level function (no architecture optimization)
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* ConjugateLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors,
							   int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							   int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  virtual ComplexVector& ConjugateLowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination);

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

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = reference on vector where result has been stored
  virtual ComplexVector& ConjugateLowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						   int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						   int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  virtual ComplexVector& ConjugateLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& ConjugateLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					     int firstComponent, int nbrComponent);
 

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = reference on vector where result has been stored
  virtual ComplexVector& ConjugateLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						      int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						      int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a set of vectors by the current hamiltonian and store result in another set of vectors
  // low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* ConjugateLowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors);

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

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and store result in another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* ConjugateLowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							   int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							   int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);
  
  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* ConjugateLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* ConjugateLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							      int firstComponent, int nbrComponent);
 

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* ConjugateLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  virtual RealVector& HermitianLowLevelMultiply(RealVector& vSource, RealVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& HermitianLowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = reference on vector where result has been stored
  virtual RealVector& HermitianLowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  virtual RealVector& HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination);

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

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = reference on vector where result has been stored
  virtual RealVector& HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						   int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						   int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a set of vectors by the current hamiltonian and store result in another set of vectors
  // low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* HermitianLowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and store result in another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* HermitianLowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							int firstComponent, int nbrComponent);

  // multiply a set of vector by the current hamiltonian for a given range of indices 
  // and store result in another set of vector, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* HermitianLowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vector sat which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* HermitianLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors);

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

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result in another set of vectors, low level function (no architecture optimization)
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* HermitianLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors,
							   int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							   int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  virtual ComplexVector& HermitianLowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& HermitianLowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					  int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = reference on vector where result has been stored
  virtual ComplexVector& HermitianLowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					  int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
					  int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  virtual ComplexVector& HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination);

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
 

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = reference on vector where result has been stored
  virtual ComplexVector& HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
					     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);

  // multiply a set of vectors by the current hamiltonian and store result in another set of vectors
  // low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* HermitianLowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and store result in another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* HermitianLowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						  int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and store result in another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* HermitianLowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						  int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						  int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);
  
  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							      int firstComponent, int nbrComponent);
 

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // sourceStart = source vector first index
  // sourceStep = step to add to go to the following source vector index
  // sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
  // sourceNbrComponent = number of component to take into account in the source vector
  // destinationStart = destination vector first index
  // destinationStep = step to add to go to the following destination vector index
  // destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
  // destinationNbrComponent = number of component to take into account in the destination vector
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							      int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							      int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent);


  // multiply a vector by the current hamiltonian and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vector where result has been stored
  virtual Vector& Multiply(Vector& vSource, Vector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual Vector& Multiply(Vector& vSource, Vector& vDestination, 
			   int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vector where result has been stored
  virtual Vector& AddMultiply(Vector& vSource, Vector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual Vector& AddMultiply(Vector& vSource, Vector& vDestination, 
			      int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors where result has to be stored
  // return value = pointer to the array of vectors where result has been stored
  virtual Vector* MultipleMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual Vector* MultipleMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors, 
				   int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual Vector* MultipleAddMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual Vector* MultipleAddMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors,
				      int firstComponent, int nbrComponent);


  // multiply a vector by the current hamiltonian and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vector where result has been stored
  virtual Vector& ConjugateMultiply(Vector& vSource, Vector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual Vector& ConjugateMultiply(Vector& vSource, Vector& vDestination, 
			   int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vector where result has been stored
  virtual Vector& ConjugateAddMultiply(Vector& vSource, Vector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual Vector& ConjugateAddMultiply(Vector& vSource, Vector& vDestination, 
			      int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors where result has to be stored
  // return value = pointer to the array of vectors where result has been stored
  virtual Vector* ConjugateMultipleMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual Vector* ConjugateMultipleMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors, 
					    int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual Vector* ConjugateMultipleAddMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual Vector* ConjugateMultipleAddMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors,
					       int firstComponent, int nbrComponent);


  // multiply a vector by the current hamiltonian and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vector where result has been stored
  virtual Vector& HermitianMultiply(Vector& vSource, Vector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual Vector& HermitianMultiply(Vector& vSource, Vector& vDestination, 
				    int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vector where result has been stored
  virtual Vector& HermitianAddMultiply(Vector& vSource, Vector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual Vector& HermitianAddMultiply(Vector& vSource, Vector& vDestination, 
				       int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors where result has to be stored
  // return value = pointer to the array of vectors where result has been stored
  virtual Vector* HermitianMultipleMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual Vector* HermitianMultipleMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors, 
					    int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSource = array of vectors to be multiplied
  // vDestination = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // return value = pointer to the array of vectors where result has been stored
  virtual Vector* HermitianMultipleAddMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual Vector* HermitianMultipleAddMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors,
					       int firstComponent, int nbrComponent);

  // test if the hamiltonian is compatible with the hamiltonian-vector multiplication operations
  //
  // return value = true if compatible (otherwise, any parallelization will be disable for these operations)
  virtual bool IsHamiltonianVectorOperationCompatible();
  
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
  

  
};

// test if the hamiltonian is compatible with the hamiltonian-vector multiplication operations
//
// return value = true if compatible (otherwise, any parallelization will be disable for these operations)
  
inline bool AbstractHamiltonian::IsHamiltonianVectorOperationCompatible()
{
  return true;
}
  
#endif
