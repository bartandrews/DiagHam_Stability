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
class Complex;
class Matrix;
class AbstractHilbertSpace;
class AbstractBitmapPicture;


class AbstractHamiltonian
{

 protected:

  GarbageFlag Flag;

 public:

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
  
  // store Hamiltonian into an hermitian matrix
  //
  // M = reference on matrix where Hamiltonian has to be stored
  // return value = reference on  corresponding hermitian matrix
  virtual HermitianMatrix& GetHamiltonian (HermitianMatrix& M);
  
  // store real part of Hamiltonian into a real symmetric matrix
  //
  // M = reference on matrix where Hamiltonian has to be stored
  // return value = reference on  corresponding real symmetric matrix 
  virtual RealSymmetricMatrix& GetHamiltonian (RealSymmetricMatrix& M);
  
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
  AbstractBitmapPicture* GetHamiltonianPicture (double error);

  // store Hamiltonian into a picture (drawing non zero element with a color scale)
  //
  // error = absolute minimum value to be considered as non zero element
  // return value = pointer to the picture associated to the matrix
  AbstractBitmapPicture* GetHamiltonianColorPicture (double error);
    
  // return a list of left interaction operators
  //
  // return value = list of left interaction operators
  virtual List<Matrix*> LeftInteractionOperators();  

  // return a list of right interaction operators
  //
  // return value = list of right interaction operators
  virtual List<Matrix*> RightInteractionOperators();  

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  virtual Complex MatrixElement (RealVector& V1, RealVector& V2) = 0;
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  virtual Complex MatrixElement (ComplexVector& V1, ComplexVector& V2) = 0;

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

};

#endif
