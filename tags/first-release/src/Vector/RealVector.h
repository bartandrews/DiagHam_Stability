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
class MathematicaOutput;
class BlockDiagonalMatrix;
class Matrix;
class ComplexVector;


class RealVector : public Vector
{

  friend class ComplexVector;
  friend class AbstractHamiltonian;
  friend class RealMatrix;
  friend class HermitianMatrix;
  friend class RealSymmetricMatrix;
  friend class RealUpperTriangularMatrix;
  friend class RealLowerTriangularMatrix;
  friend class RealDiagonalMatrix;
  friend class RealAntisymmetricMatrix;
  friend class RealTriDiagonalSymmetricMatrix;
  friend Complex operator * (const ComplexVector& V1, const RealVector& V2);
  friend Complex operator * (const RealVector& V1, const ComplexVector& V2);
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
  friend RealMatrix operator * (const RealMatrix& M, double x);
  friend RealMatrix operator * (double x, const RealMatrix& M);
  friend RealMatrix operator / (const RealMatrix& M, double x);
  friend ostream& operator << (ostream& Str, const RealMatrix& P);
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const RealMatrix& P);


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

  // constructor from an array of doubles
  //
  // array = array of doubles with real in even position and imaginary part in odd position
  // size = Vector Dimension  
  RealVector(double* array, int size);

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

  // destructor
  //
  ~RealVector ();

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
  void Resize (int dimension);

  // Resize vector and set to zero all components that have been added
  //
  // dimension = new dimension
  void ResizeAndClean (int dimension);

  // copy a vector into another
  //
  // vector = vector to copy
  // coefficient = optional coefficient which multiply source to copy
  // return value = reference on current vector
  RealVector& Copy (RealVector& vector, double coefficient = 1.0);

  // create a new vector with same size and same type but without duplicating datas
  //
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  Vector* EmptyClone(bool zeroFlag = false);

  // put all vector components to zero
  //
  // return value = reference on current vector
  Vector& ClearVector ();

  // change sign of a vector
  //
  // return value = reference on current vector
  RealVector& operator - ();

  // return a new vector with opposite sign form a given source vector
  //
  // V1 = source vector
  // return value = new vector
  friend RealVector operator - (const RealVector& V1);

  // scalar product between two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = result of scalar product
  friend double operator * (const RealVector& V1, const RealVector& V2);

  // sum two vectors
  //
  // V1 = vector to add
  // return value = reference on current vector
  RealVector& operator += (const RealVector& V1);

  // sum two vectors
  //
  // vector = vector to add
  // return value = reference on current vector
  Vector& operator += (const Vector& vector);

  // substract two vectors
  //
  // V1 = first vector
  // return value = reference on current vector
  RealVector& operator -= (const RealVector& V1);

  // sum two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = resulting vector
  friend RealVector operator + (const RealVector& V1, const RealVector& V2);

  // substract two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = resulting vector
  friend RealVector operator - (const RealVector& V1, const RealVector& V2);

  // add a linear combination to a given vector
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  RealVector& AddLinearCombination (const double& x, const RealVector& V);

  // add a linear combination to a given vector, for a given range of indices
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  RealVector& AddLinearCombination (double x, const RealVector& V, int firstComponent, int nbrComponent);

  // add a linear combination of two vectors to a given vector
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // return value = reference on current vector
  RealVector& AddLinearCombination (double x1, const RealVector& v1, double x2, const RealVector& v2);

  // add a linear combination of two vectors to a given vector, for a given range of indices
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on current vector
  RealVector& AddLinearCombination (double x1, const RealVector& v1, double x2, 
				    const RealVector& v2, int firstComponent, int nbrComponent);

  // multiply a vector with a real number on the right hand side
  //
  // V1 = vector to multiply
  // d = real to use
  // return value = resulting vector
  friend RealVector operator * (const RealVector& V1, double d);

  // multiply a vector with a real number on the left hand side
  //
  // V1 = vector to multiply
  // d = real to use
  // return value = resulting vector
  friend RealVector operator * (double d, const RealVector& V1);

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

  // left multiply a vector with an antissymmetric matrix and use to store result in current vector (without creating temporary vector)
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

  // left multiply a vector with an antissymmetric matrix and use to store result in current vector (without creating temporary vector)
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
  double& operator [] (int i);

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
  
  // Merge a subvector into a given vector
  //
  // V = vector to merge
  // firstCoordinate = Coordinate where merge has to begin
  // step = distance to the next coordinate in the destination vector (1 means to take the following)
  // return value = reference to the current Vector
  RealVector& Merge(const RealVector& V, int firstCoordinate, int step = 1);
  
  // write vector in a file 
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool WriteVector (char* fileName);

  // write vector in a file in ascii mode
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool WriteAsciiVector (char* fileName);

  // read vector from a file 
  //
  // fileName = name of the file where the vector has to be read
  // return value = true if no error occurs
  bool ReadVector (char* fileName);

  // Output Stream overload
  //
  // str = reference on output stream
  // v = vector to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, const RealVector& v);

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

};

// return vector i-th coordinate (without testing if position is valid)
//
// i = coordinate position

inline double& RealVector::operator [] (int i)
{
  return this->Components[i];
}
 

#endif

