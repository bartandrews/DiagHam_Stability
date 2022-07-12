////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                    class for n dimensional real vector                     //
//                                                                            //
//                        last modification : 04/01/2001                      //
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


#ifndef REALVECTOR_H
#define REALVECTOR_H


#include "config.h"
#include "Vector/Vector.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>
#include <fstream>


using std::ostream;
using std::ifstream;
using std::ofstream;


class Complex;
class ComplexVector;
class RealPtrVector;
class MathematicaOutput;
class BlockDiagonalMatrix;
class Matrix;
class RealMatrix;
class RealAntisymmetricMatrix;
class RealDiagonalMatrix;
class RealSymmetricMatrix;
class RealTriDiagonalSymmetricMatrix;
class ComplexVector;


class RealVector : public Vector
{

  friend class DelocalizedRealVector;
  friend class ComplexVector;
  friend class RealPtrVector;
  friend class AbstractHamiltonian;
  friend class RealMatrix;
  friend class ComplexMatrix;
  friend class HermitianMatrix;
  friend class RealSymmetricMatrix;
  friend class RealUpperTriangularMatrix;
  friend class RealLowerTriangularMatrix;
  friend class RealDiagonalMatrix;
  friend class RealAntisymmetricMatrix;
  friend class RealTriDiagonalSymmetricMatrix;
  friend class SparseRealMatrix;


  friend Complex operator * (const ComplexVector& V1, const RealVector& V2);
  friend Complex operator * (const RealVector& V1, const ComplexVector& V2);
  friend double operator * (RealPtrVector& V1, RealVector& V2);
  friend double operator * (RealVector& V1, RealPtrVector& V2);
  friend ComplexVector operator + (const RealVector& V1, const ComplexVector& V2);
  friend ComplexVector operator + (const ComplexVector& V1, const RealVector& V2);
  friend ComplexVector operator - (const RealVector& V1, const ComplexVector& V2);
  friend ComplexVector operator - (const ComplexVector& V1, const RealVector& V2);
  friend RealMatrix operator + (const RealMatrix& M1, const RealMatrix& M2);
  friend RealMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const RealMatrix& M2);
  friend RealMatrix operator + (const RealMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);
  friend RealMatrix operator - (const RealMatrix& M1, const RealMatrix& M2);
  friend RealMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, const RealMatrix& M2);
  friend RealMatrix operator - (const RealMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);
  friend RealMatrix operator * (const RealMatrix& M1, const RealMatrix& M2);
  friend RealMatrix operator * (const RealDiagonalMatrix& M1, const RealMatrix& M2);
  friend RealMatrix operator * (const RealMatrix& M1, const RealDiagonalMatrix& M2);
  friend RealMatrix operator * (const RealMatrix& M, double x);
  friend RealMatrix operator * (double x, const RealMatrix& M);
  friend RealMatrix operator / (const RealMatrix& M, double x);
  friend RealMatrix operator / (const RealMatrix& M1, const RealDiagonalMatrix& M2);

  friend RealVector Cross(RealVector& V);
  friend RealVector Cross(RealVector& V1, RealVector& V2);

  friend ostream& operator << (ostream& Str, const RealMatrix& P);
#ifdef USE_OUTPUT
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const RealMatrix& P);
#endif

 protected:
  
  double* Components;
  GarbageFlag Flag;

 public:

  // default constructor
  //
  RealVector();

  // constructor for an empty real vector 
  //
  // size = Vector Dimension 
  // zeroFlag = true if all coordinates have to be set to zero
  RealVector(int size, bool zeroFlag = false);

  // constructor for an empty real vector bigger than 2^31
  //
  // size = Vector Dimension 
  // zeroFlag = true if all coordinates have to be set to zero
  RealVector(long size, bool zeroFlag = false);

  // constructor from an array of doubles
  //
  // array = array of doubles to become Components of vector
  // size = Vector Dimension  
  RealVector(double* array, int size);

  // constructor from an array of doubles
  //
  // array = array of doubles to become Components of vector
  // size = Vector Dimension
  RealVector(double* array, long size);

  // copy constructor
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  RealVector(const RealVector& vector, bool duplicateFlag = false);

  // copy constructor from a complex vector (keep only real part and datas are duplicated)
  //
  // vector = vector to copy
  RealVector(const ComplexVector& vector);

  // copy constructor from a vector (duplicate datas if necessary)
  //
  // vector = vector to copy
  RealVector(const Vector& vector);

#ifdef __MPI__
  // constructor from informations sent using MPI
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts or sends the vector
  // broadcast = true if the vector is broadcasted
  RealVector(MPI::Intracomm& communicator, int id, bool broadcast = true);
#endif

  // destructor
  //
  virtual ~RealVector ();

  // assignement
  //
  // vector = vector to assign
  // return value = reference on current vector
  RealVector& operator = (const RealVector& vector);

  // assignement from a complex vector (keep only real part and datas are duplicated)
  //
  // vector = vector to assign
  // return value = reference on current vector
  RealVector& operator = (const ComplexVector& vector);

  // assignement from a vector (duplicate datas if necessary)
  //
  // vector = vector to assign
  // return value = reference on current vector
  RealVector& operator = (const Vector& vector);

  // Resize vector
  //
  // dimension = new dimension
  virtual void Resize (int dimension);

  // Resize vector
  //
  // dimension = new dimension
  virtual void Resize (long dimension);

  // Resize vector and set to zero all components that have been added
  //
  // dimension = new dimension
  virtual void ResizeAndClean (int dimension);

  // copy a vector into another
  //
  // vector = vector to copy
  // return value = reference on current vector
  RealVector& Copy (RealVector& vector);

  // copy a vector into another
  //
  // vector = vector to copy
  // coefficient = optional coefficient which multiply source to copy
  // return value = reference on current vector
  RealVector& Copy (RealVector& vector, double coefficient);

  // create a new vector with same size and same type but without duplicating datas
  //
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  virtual Vector* EmptyClone(bool zeroFlag = false);

  // create an array of new vectors with same size and same type but non-initialized components
  //
  // nbrVectors = number of vectors to sreate
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to the array of new vectors
  virtual Vector* EmptyCloneArray(int nbrVectors, bool zeroFlag = false);

  // put all vector components to zero
  //
  // return value = reference on current vector
  virtual Vector& ClearVector ();

  // put select vector components to zero
  // start = start index
  // nbrComponent = number of components to set to zero
  // return value = reference on current vector
  virtual Vector& ClearVectorSegment (long start, long nbrComponent);

  // change sign of a vector
  //
  // return value = reference on current vector
  RealVector& operator - ();

  // return a new vector with opposite sign form a given source vector
  //
  // V1 = source vector
  // return value = new vector
  friend RealVector operator - (RealVector& V1);

  // scalar product between two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = result of scalar product
  friend double operator * (const RealVector& V1, const RealVector& V2);

  // do part of the scalar product between two vectors in a given range of indices
  //
  // vRight = right vector of the scalar product
  // firstComponent = index of the first component to consider
  // nbrComponent = number of components to consider
  // step = increment between to consecutive indices to consider
  // return value = result of the partial scalar product
  double PartialScalarProduct (RealVector& vRight, int firstComponent, int nbrComponent, int step = 1);

  // sum two vectors
  //
  // V1 = vector to add
  // return value = reference on current vector
  RealVector& operator += (RealVector& V1);

  // sum two vectors
  //
  // vector = vector to add
  // return value = reference on current vector
  Vector& operator += (Vector& vector);

  // substract two vectors
  //
  // V1 = first vector
  // return value = reference on current vector
  RealVector& operator -= (RealVector& V1);

  // sum two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = resulting vector
  friend RealVector operator + (RealVector& V1, RealVector& V2);

  // substract two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = resulting vector
  friend RealVector operator - (RealVector& V1, RealVector& V2);

  // add a linear combination to a given vector
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  RealVector& AddLinearCombination (const double& x, RealVector& V);

  // add a linear combination to a given vector, for a given range of indices
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  RealVector& AddLinearCombination (double x, RealVector& V, int firstComponent, int nbrComponent);

  // add a linear combination of two vectors to a given vector
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // return value = reference on current vector
  RealVector& AddLinearCombination (double x1, RealVector& v1, double x2, RealVector& v2);

  // add a linear combination of two vectors to a given vector, for a given range of indices
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on current vector
  RealVector& AddLinearCombination (double x1, RealVector& v1, double x2, 
				    RealVector& v2, int firstComponent, int nbrComponent);

  // multiply a vector with a real number on the right hand side
  //
  // V1 = vector to multiply
  // d = real to use
  // return value = resulting vector
  friend RealVector operator * (RealVector& V1, double d);

  // multiply a vector with a real number on the left hand side
  //
  // V1 = vector to multiply
  // d = real to use
  // return value = resulting vector
  friend RealVector operator * (double d, RealVector& V1);

  // multiply a vector with a real number on the right hand side
  //
  // d = real to use
  // return value = reference on current vector
  RealVector& operator *= (double d);

  // divide a vector by a real number on the right hand side
  //
  // d = real to use
  // return value = reference on current vector
  RealVector& operator /= (double d);

  // left multiply a vector with a real symmetric matrix (without using temporary vector)
  //
  // M = matrix to use
  // return value = reference on current vector
  RealVector& operator *= (const RealSymmetricMatrix&  M);

  // left multiply a vector with a real tridiagonal symmetric matrix (without using temporary vector)
  //
  // M = matrix to use
  // return value = reference on current vector
  RealVector& operator *= (const RealTriDiagonalSymmetricMatrix&  M);

  // left multiply a vector with a symmetric matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply  
  // return value = reference on current vector
  RealVector& Multiply (const RealSymmetricMatrix&  M, RealVector& V);

  // do a partial left multication of a vector with a real symmetric matrix and store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  RealVector& Multiply (const RealSymmetricMatrix&  M, RealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with a symmetric matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& Multiply (const RealSymmetricMatrix&  M, RealVector& V, int sourceStart, int sourceStep, int destStart, int destStep);

  // left multiply a vector with a symmetric matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& Multiply (const RealSymmetricMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
			int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with a real matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  RealVector& AddMultiply (const RealSymmetricMatrix&  M, RealVector& V);

  // do a partial left multication of a vector with a real symmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  RealVector& AddMultiply (const RealSymmetricMatrix&  M, RealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with a antisymmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply  
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& AddMultiply (const RealSymmetricMatrix&  M, RealVector& V, int sourceStart, 
			   int sourceStep, int destStart, int destStep);

  // left multiply a vector with a antisymmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply  
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& AddMultiply (const RealSymmetricMatrix&  M, RealVector& V, int sourceStart, 
			   int sourceNbrComponent, int sourceStep, int destStart, int destStep);

  // left multiply a vector with an antisymmetric matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply  
  // return value = reference on current vector
  RealVector& Multiply (const RealAntisymmetricMatrix&  M, RealVector& V);

  // do a partial left multication of a vector with a real antisymmetric matrix and store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  RealVector& Multiply (const RealAntisymmetricMatrix&  M, RealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with an antisymmetric matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& Multiply (const RealAntisymmetricMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
			int destStart, int destStep);

  // left multiply a vector with an antisymmetric matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& Multiply (const RealAntisymmetricMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
			int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with an antisymmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  RealVector& AddMultiply (const RealAntisymmetricMatrix&  M, RealVector& V);

  // do a partial left multication of a vector with a real antisymmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  RealVector& AddMultiply (const RealAntisymmetricMatrix&  M, RealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with an antisymmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply  
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& AddMultiply (const RealAntisymmetricMatrix&  M, RealVector& V, int sourceStart, int sourceStep, int destStart, int destStep);

  // left multiply a vector with an antisymmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply  
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& AddMultiply (const RealAntisymmetricMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
			   int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with a real diagonal matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply  
  // return value = reference on current vector
  RealVector& Multiply (const RealDiagonalMatrix&  M, RealVector& V);

  // do a partial left multication of a vector with a real matrix and store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  RealVector& Multiply (const RealDiagonalMatrix&  M, RealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with a real diagonal matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& Multiply (const RealDiagonalMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
			int destStart, int destStep);

  // left multiply a vector with a real diagonal matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& Multiply (const RealDiagonalMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
			int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with a real diagonal matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  RealVector& AddMultiply (const RealDiagonalMatrix&  M, RealVector& V);
  
  // do a partial left multication of a vector with a real matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  RealVector& AddMultiply (const RealDiagonalMatrix&  M, RealVector& V, int sourceStart, int sourceNbrComponent);
    
  // left multiply a vector with a real diagonal matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply  
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& AddMultiply (const RealDiagonalMatrix&  M, RealVector& V, int sourceStart, int sourceStep, int destStart, int destStep);

  // left multiply a vector with a real diagonal matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply  
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& AddMultiply (const RealDiagonalMatrix&  M, RealVector& V, int sourceStart, 
			   int sourceNbrComponent, int sourceStep, int destStart, int destStep);

  // left multiply a vector with a real matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  RealVector& Multiply (const RealMatrix&  M, RealVector& V);

  // left multiply a vector with the transpose of a real matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  RealVector& TransposeMultiply (const RealMatrix&  M, RealVector& V);

  // do a partial left multication of a vector with a real matrix and store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  RealVector& Multiply (const RealMatrix&  M, RealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with a real matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector  
  RealVector& Multiply (const RealMatrix&  M, RealVector& V, int sourceStart, int sourceStep, int destStart, int destStep);

  // left multiply a vector with a real matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector  
  RealVector& Multiply (const RealMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
			int sourceNbrComponent, int destStart, int destStep);

  // do a partial left multication of a vector with a real matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  RealVector& AddMultiply (const RealMatrix&  M, RealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with a real matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& AddMultiply (const RealMatrix&  M, RealVector& V, int sourceStart, int sourceStep, int destStart, int destStep);

  // left multiply a vector with a real matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& AddMultiply (const RealMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
			   int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with a block-diagonal matrix and use to store result 
  // in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  RealVector& Multiply (const BlockDiagonalMatrix&  M, RealVector& V);

  // left multiply a vector with a block-diagonal matrix and add result to the current vector
  // in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  RealVector& AddMultiply (const BlockDiagonalMatrix&  M, RealVector& V);

  // left multiply a vector with a block-diagonal matrix and use to store result 
  // in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector  
  RealVector& Multiply (const BlockDiagonalMatrix&  M, RealVector& V, int sourceStart, 
			int sourceStep, int destStart, int destStep);

  // left multiply a vector with a block-diagonal matrix and use to store result 
  // in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector  
  RealVector& Multiply (const BlockDiagonalMatrix&  M, RealVector& V, int sourceStart, 
			int sourceStep, int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with a block-diagonal matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& AddMultiply (const BlockDiagonalMatrix&  M, RealVector& V, int sourceStart, 
			   int sourceStep, int destStart, int destStep);

  // left multiply a vector with a block-diagonal matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& AddMultiply (const BlockDiagonalMatrix&  M, RealVector& V, int sourceStart, 
			   int sourceStep, int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with a matrix and use to store result in current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  RealVector& Multiply (const Matrix&  M, RealVector& V);

  // left multiply a vector with a matrix and add result to current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  RealVector& AddMultiply (const Matrix&  M, RealVector& V);

  // do a partial left multication of a vector with a real matrix and store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& Multiply (const Matrix&  M, RealVector& V, int sourceStart, int sourceNbrComponent);

  // do a partial left multication of a vector with a real matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& AddMultiply (const Matrix&  M, RealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with a matrix and use to store result in current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& Multiply (const Matrix&  M, RealVector& V, int sourceStart, int sourceStep, 
			int destStart, int destStep);

  // left multiply a vector with a matrix and add current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& AddMultiply (const Matrix&  M, RealVector& V, int sourceStart, int sourceStep, 
			   int destStart, int destStep);

  // left multiply a vector with a matrix and use to store result in current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& Multiply (const Matrix&  M, RealVector& V, int sourceStart, int sourceStep, 
			int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with a matrix and add current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  RealVector& AddMultiply (const Matrix&  M, RealVector& V, int sourceStart, int sourceStep, 
			   int sourceNbrComponent, int destStart, int destStep);

  // return vector i-th coordinate (without testing if position is valid)
  //
  // i = coordinate position
  virtual double& operator [] (int i);

  // return vector i-th coordinate (without testing if position is valid)
  //
  // i = coordinate position
  virtual double& operator [] (long i);

  // get vector norm
  //
  // return value = vector norm
  double Norm();
  
  // get square of vector norm
  //
  // return value = square of vector norm
  double SqrNorm ();
  
  // normalize vector
  //
  // return value = reference on current vector
  RealVector& Normalize();

  // orthonormalized a vector with respect to a set of orthonormalized vectors
  //
  // vectors = vector array corresponding to the set
  // nbrVectors = number of vectors in the set
  // return value = resulting vector norm (can be used to see if vector is can be decomposed on vector set)
  double Orthonormalized (RealVector* vectors, int nbrVectors);

  // Extract a subvector from a given vector
  //
  // firstCoordinate = Coordinate where extraction has to begin
  // lastCoordinate = Coordinate where extraction has to stop (extract also include this last coordinate)
  // step = distance to the next coordinate (1 means to take the folowing)
  // return value = return corresponding subvector
  RealVector Extract(int firstCoordinate, int lastCoordinate, int step = 1);

  // get a serie of components
  //
  // indices = array where component indices are stored
  // nbrIndices = number of components to retrieve
  // components = array where retrieved components will be stored
  virtual void GetMultipleComponents(long* indices, long nbrIndices, double* components);

  // Merge a subvector into a given vector
  //
  // V = vector to merge
  // firstCoordinate = Coordinate where merge has to begin
  // step = distance to the next coordinate in the destination vector (1 means to take the following)
  // return value = reference to the current Vector
  RealVector& Merge(RealVector& V, int firstCoordinate, int step = 1);

  // reverse elements of the current vector (i.e. exchanging i <-> N - i)
  //
  // return value = reference to the current Vector
  RealVector& ReverseVector();

  // write vector in a file 
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool WriteVector (const char* fileName);
  bool ByteWriteVector (const char* fileName);

  // write vector in a file in ascii mode
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool WriteAsciiVector (const char* fileName);

  // read vector from a file 
  //
  // fileName = name of the file where the vector has to be read
  // return value = true if no error occurs
  bool ReadVector (const char* fileName);
  bool ByteReadVector (const char* fileName);

  // read vector from a file, only within a given range of indices
  //
  // fileName = name of the file where the vector has to be read
  // minIndex = index of the first component to read
  // maxIndex = index of the last component to read
  // return value = true if no error occurs
  bool ReadVector (const char* fileName, long minIndex, long maxIndex);

  // read vector dimension from a file, without loading the full vector 
  //
  // fileName = name of the file where the vector has to be read
  // return value = vector dimension
  long ReadVectorDimension (const char* fileName);

  // test if a vector can be read from a file (matching the right type), without loading the full vector 
  //
  // fileName = name of the file where the vector has to be read
  // return value = true if the vector can be read
  bool ReadVectorTest (const char* fileName);
  
  // Output Stream overload
  //
  // str = reference on output stream
  // v = vector to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, RealVector& v);

  // print a single component
  //
  // str = reference on output stream
  // index = index of the component to print
  // return value  = reference on output stream
  virtual ostream& PrintComponent(ostream& str, long index);

  // output the vector in a sparse display
  //
  // str = reference on output stream
  // error = numerical accuracy below which a vector component is considered to be equal to zero
  // return value = reference on output stream  
  virtual ostream& PrintNonZero(ostream& str, double error = MACHINE_PRECISION);

  // output the vector in a sparse display, using labels for the component indices
  //
  // str = reference on output stream
  // componentLabels = array of labels for the component indices
  // error = numerical accuracy below which a vector component is considered to be equal to zero
  // return value = reference on output stream  
  virtual ostream& PrintNonZero(ostream& str, char** componentLabels, double error = MACHINE_PRECISION);

  // output file stream overload
  //
  // file = reference on output file stream
  // vector = reference on vector to save
  // return value = reference on output file stream
//  friend ofstream& operator << (ofstream& file, const RealVector& vector);

  // input file stream overload
  //
  // file = reference on output file stream
  // vector = reference on vector to load
  // return value = reference on output file stream
  friend ifstream& operator >> (ifstream& file, RealVector& vector);

#ifdef __MPI__

  // send a vector to a given MPI process
  // 
  // communicator = reference on the communicator to use
  // id = id of the destination MPI process
  // return value = reference on the current vector
  virtual Vector& SendVector(MPI::Intracomm& communicator, int id);

  // broadcast a vector to all MPI processes associated to the same communicator
  // 
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // return value = reference on the current vector
  virtual Vector& BroadcastVector(MPI::Intracomm& communicator,  int id);

  // broadcast part of vector to all MPI processes associated to the same communicator
  // 
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // firstComponent = index of the first component (useless if the method is not called by the MPI process which broadcasts the vector)
  // nbrComponent = number of component (useless if the method is not called by the MPI process which broadcasts the vector)
  // return value = reference on the current vector
  virtual Vector& BroadcastPartialVector(MPI::Intracomm& communicator, int id, int firstComponent = 0, int nbrComponent = 0);

  // receive a vector from a MPI process
  // 
  // communicator = reference on the communicator to use 
  // id = id of the source MPI process
  // return value = reference on the current vector
  virtual Vector& ReceiveVector(MPI::Intracomm& communicator, int id);

  // add current vector to the current vector of a given MPI process
  // 
  // communicator = reference on the communicator to use 
  // id = id of the destination MPI process
  // return value = reference on the current vector
  virtual Vector& SumVector(MPI::Intracomm& communicator, int id);

  // reassemble vector from a scattered one
  // 
  // communicator = reference on the communicator to use 
  // id = id of the destination MPI process
  // return value = reference on the current vector
  Vector& ReassembleVector(MPI::Intracomm& communicator, int id);

  // gather vector from a scattered one
  // 
  // communicator = reference on the communicator to use 
  // id = id of the destination MPI process
  // return value = reference on the current vector
  Vector& GatherVector(MPI::Intracomm& communicator, int id);

  // create a new vector on each MPI node which is an exact clone of the broadcasted one
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  virtual Vector* BroadcastClone(MPI::Intracomm& communicator, int id);

  // create a new vector on given MPI node which is an exact clone of the sent one but with only part of the data
  // 
  // communicator = reference on the communicator to use
  // id = id of the destination MPI process
  // firstComponent = index of the first component 
  // nbrComponent = number of component to send
  // return value = reference on the current vector
  virtual Vector& SendPartialClone(MPI::Intracomm& communicator, int id, int firstComponent, int nbrComponent);

  // scatter this vector across all MPI nodes with the given load balancing information
  // 
  // communicator = reference on the communicator to use
  // minimumIndices = lowest index for each thread
  // maximumIndices = largest index for each thread
  // id = id of the process to send the vector
  // return value = reference on the current vector
  Vector& ScatterPartialClones(MPI::Intracomm& communicator, long *minimumIndices, long *maximumIndices, int id);

  // create a new vector on each MPI node with same size and same type but non-initialized components
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  virtual Vector* BroadcastEmptyClone(MPI::Intracomm& communicator, int id, bool zeroFlag = false);

#endif

  // query whether blas is being used:
  // verbose = flag indicating whether to print a message on screen
  bool HaveBlas(bool verbose = false);

};

// return vector i-th coordinate (without testing if position is valid)
//
// i = coordinate position

inline double& RealVector::operator [] (int i)
{
  return this->Components[i];
}
 
// return vector i-th coordinate (without testing if position is valid)
//
// i = coordinate position

inline double& RealVector::operator [] (long i)
{
  return this->Components[i];
}
 
// print a single component
//
// str = reference on output stream
// index = index of the component to print
// return value  = reference on output stream

inline ostream& RealVector::PrintComponent(ostream& str, long index)
{
  str << this->Components[index];
  return str;
}


#endif

