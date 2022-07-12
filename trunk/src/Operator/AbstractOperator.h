////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of abstract operator                       //
//                                                                            //
//                        last modification : 22/03/2002                      //
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


#ifndef ABSTRACTOPERATOR_H
#define ABSTRACTOPERATOR_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include <iostream>


class Complex;
class Matrix;
class RealSymmetricMatrix;
class HermitianMatrix;
class AbstractHilbertSpace;
class IntegerMatrix;
class RealMatrix;
class ComplexMatrix;
class LongIntegerMatrix;
#ifdef USE_OUTPUT
class MathematicaOutput;
#endif


using std::ostream;


class AbstractOperator
{

 protected:

  GarbageFlag Flag;

 public:

  // destructor
  //
  virtual ~AbstractOperator();

  // clone operator without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  virtual AbstractOperator* Clone () = 0;

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  virtual void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace) = 0;

  // get Hilbert space on which operator acts
  //
  // return value = pointer to used Hilbert space
  virtual AbstractHilbertSpace* GetHilbertSpace () = 0;

  // return dimension of Hilbert space where operator acts
  //
  // return value = corresponding matrix elementdimension
  virtual int GetHilbertSpaceDimension () = 0;
  
  // store operator into an hermitian matrix
  //
  // M = reference on matrix where operator has to be stored
  // return value = reference on  corresponding hermitian matrix
  virtual HermitianMatrix& GetOperator (HermitianMatrix& M);
  
  // store the real part of an operator into an integer matrix
  //
  // M = reference on matrix where the operator has to be stored
  // return value = reference on  corresponding real symmetric matrix 
  virtual RealSymmetricMatrix& GetOperator (RealSymmetricMatrix& M);
  
  // store the real part of an operator into a real matrix
  //
  // M = reference on matrix where the operator has to be stored
  // return value = reference on corresponding real matrix 
  virtual RealMatrix& GetOperator (RealMatrix& M);
  
  // store the real part of an operator into a complex matrix
  //
  // M = reference on matrix where the operator has to be stored
  // return value = reference on corresponding complex matrix 
  virtual ComplexMatrix& GetOperator (ComplexMatrix& M);
  
  // store the real part of an operator into a long integer matrix
  //
  // M = reference on matrix where the operator has to be stored
  // scalingFactor = use an additional scaling factor before converting coefficients into integers
  // return value = reference on corresponding matrix 
  virtual IntegerMatrix& GetOperator (IntegerMatrix& M, double scalingFactor = 1.0);
  
  // store the real part of Hamiltonian into a long integer matrix
  //
  // M = reference on matrix where Hamiltonian has to be stored
  // scalingFactor = use an additional scaling factor before converting coefficients into integers
  // return value = reference on corresponding matrix 
  virtual LongIntegerMatrix& GetOperator (LongIntegerMatrix& M, double scalingFactor = 1.0);
  
  // store real part of operator into a matrix
  //
  // M = reference on matrix where operator has to be stored
  // return value = reference on  corresponding matrix 
  virtual Matrix& GetOperator (Matrix& M);
  
  // return matrix representation of current operator
  //
  // return value = reference to representation
  virtual Matrix* GetOperator ();
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  virtual Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate part of the matrix element, within a given of indices
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = corresponding matrix element
  virtual Complex PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent);

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  virtual Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // evaluate part of the matrix element, within a given of indices
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = corresponding matrix element
  virtual Complex PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent);

  // multiply a vector by the current operator for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored  
  virtual Vector& Multiply(Vector& vSource, Vector& vDestination, 
			   int firstComponent, int nbrComponent);


  // multiply a vector by the current operator and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  virtual RealVector& Multiply(RealVector& vSource, RealVector& vDestination);

  // multiply a vector by the current operator for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
				       int firstComponent, int nbrComponent);

  // multiply a vector by the current operator for a given range of indices 
  // and add result to another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
					  int firstComponent, int nbrComponent);
 
  // multiply a vector by the current operator and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vector where result has been stored
  virtual ComplexVector& Multiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current operator for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					  int firstComponent, int nbrComponent);

  // multiply a vector by the current operator for a given range of indices 
  // and add result to another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					     int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current operator for a given range of indices 
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

  // multiply a set of vectors by the current operator for a given range of indices 
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

  // multiply a set of vectors by the current operator for a given range of indices 
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

  // multiply a set of vectors by the current operator for a given range of indices 
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

  // Output Stream overload
  //
  // Str = reference on output stream
  // O = Operator to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, AbstractOperator& O);

#ifdef USE_OUTPUT
  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // O = Operator to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, AbstractOperator& O);
#endif

};

#endif
