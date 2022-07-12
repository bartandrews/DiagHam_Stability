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
  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination) = 0;

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
				       int firstComponent, int nbrComponent) = 0;

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination) = 0;

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
					  int firstComponent, int nbrComponent) = 0;

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  virtual ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) = 0;

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					  int firstComponent, int nbrComponent) = 0;

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination) = 0;

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					     int firstComponent, int nbrComponent) = 0;
 

  // multiply a vector by the current hamiltonian and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
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
  // return value = reference on vectorwhere result has been stored
  virtual Vector& AddMultiply(Vector& vSource, Vector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual Vector& AddMultiply(Vector& vSource, Vector& vDestination, 
			      int firstComponent, int nbrComponent);

#ifdef __SMP__

  // multiply a vector by the current hamiltonian using threads
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // nbrProcess = number of process to run
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, int nbrProcess);

  // multiply a vector by the current hamiltonian using threads
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // nbrProcess = number of process to run
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, int nbrProcess);

#endif

  // Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step
  //
  // dimension = maximum iteration number
  // M = reference on complex tridiagonal hermitian matrix where result has to be stored
  // V1 = reference on complex vector used as first vector (will contain last produced vector at the end)
  // return value = reference on complex tridiagonal hermitian matrix
  virtual RealTriDiagonalSymmetricMatrix& Lanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, ComplexVector& V1);

  // Tridiagonalize hamiltonian using Lanczos algorithm with partial base re-orthogonalization
  //
  // dimension = maximum iteration number
  // M = reference on complex tridiagonal hermitian matrix where result has to be stored
  // V1 = reference on real vector used as first vector (will contain last produced vector at the end)
  // step = number of iterations before re-orthogonalizing whole base
  // return value = reference on complex tridiagonal hermitian matrix
  virtual RealTriDiagonalSymmetricMatrix& ReorthogonalizedLanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, 
								   RealVector& V1, int step = 1);

  // Tridiagonalize hamiltonian using Lanczos algorithm with base re-orthogonalization
  //
  // dimension = maximum iteration number
  // M = reference on real tridiagonal hermitian matrix where result has to be stored
  // V1 = reference on real vector used as first vector (will contain last produced vector at the end)
  // return value = reference on complex tridiagonal hermitian matrix
  virtual RealTriDiagonalSymmetricMatrix& FullReorthogonalizedLanczos (int dimension, RealTriDiagonalSymmetricMatrix& M,  RealVector& V1);

  // Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step
  //
  // dimension = maximum iteration number
  // M = reference on complex tridiagonal hermitian matrix where result has to be stored
  // V1 = reference on complex vector used as first vector (will contain last produced vector at the end)
  // return value = reference on complex tridiagonal hermitian matrix
  virtual RealTriDiagonalSymmetricMatrix& Lanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, RealVector& V1);

  // Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step
  //
  // dimension = maximum iteration number
  // M = reference on real tridiagonal symmetric matrix where result has to be stored
  // V1 = reference on real vector used as first vector (will contain last produced vector at the end)
  // V2 = reference on real vector used as second vector (will contain next to last produced vector at the end)
  // useV2 = true if V2 already contains result of second Lanczos iteration, in that case M is supposed to give
  //         results of previous Lanczos iteration
  // return value = reference on real tridiagonal symmetric matrix

  RealTriDiagonalSymmetricMatrix& Lanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, 
					   RealVector& V1, RealVector& V2, bool useV2 = false);

};

#endif
