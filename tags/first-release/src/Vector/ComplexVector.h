////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                   class for n dimensional complex vector                   //
//                                                                            //
//                        last modification : 17/01/2001                      //
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


#ifndef COMPLEXVECTOR_H
#define COMPLEXVECTOR_H


#include "config.h"
#include "Vector/Vector.h"
#include "Vector/RealVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexTriDiagonalHermitianMatrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Complex.h"
#include "GeneralTools/GarbageFlag.h"


#include <iostream>


using std::ostream;


class MathematicaOutput;
class RealMatrix;
class ComplexMatrix;


class ComplexVector : public Vector
{

  friend class RealVector;
  friend class AbstractHamiltonian;
  friend class ComplexMatrix;
  friend class HermitianMatrix;
  friend class RealTriDiagonalSymmetricMatrix;
  friend class ComplexTriDiagonalHermitianMatrix;
  friend ComplexMatrix operator + (const ComplexMatrix& M1, const ComplexMatrix& M2);
  friend ComplexMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const ComplexMatrix& M2);
  friend ComplexMatrix operator + (const ComplexMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);
  friend ComplexMatrix operator - (const ComplexMatrix& M1, const ComplexMatrix& M2);
  friend ComplexMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, const ComplexMatrix& M2);
  friend ComplexMatrix operator - (const ComplexMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);
  friend ComplexMatrix operator * (const ComplexMatrix& M1, const ComplexMatrix& M2);
  friend ComplexMatrix operator * (const ComplexMatrix& M, double x);
  friend ComplexMatrix operator * (double x, const ComplexMatrix& M);
  friend ComplexMatrix operator / (const ComplexMatrix& M, double x);
  friend ostream& operator << (ostream& Str, const ComplexMatrix& P);
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexMatrix& P);

 protected:
  
  double* RealComponents;
  double* ImaginaryComponents;
  GarbageFlag Flag;

 public:

  // default constructor
  //
  ComplexVector();

  // constructor for an empty real vector
  //
  // size = Vector Dimension 
  // zeroFlag = true if all coordinates have to be set to zero
  ComplexVector(int size, bool zeroFlag = false);
  
  // constructor from arrays of doubles
  //
  // real = array of doubles corresponding to real part
  // imaginary = array of doubles corresponding to imaginary part
  // size = Vector Dimension 
  ComplexVector(double* real, double* imaginary, int size) ;

  // copy constructor
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  ComplexVector(const ComplexVector& vector, bool duplicateFlag = false);

  // copy constructor from a real vector
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  ComplexVector(const RealVector& vector, bool duplicateFlag = false);

  // destructor
  //
  ~ComplexVector ();

  // assignement
  //
  // vector = vector to assign
  // return value = reference on current vector
  ComplexVector& operator = (const ComplexVector& vector);

  // assignement from a real vector
  //
  // vector = vector to assign
  // return value = reference on current vector
  ComplexVector& operator = (const RealVector& vector);

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
  ComplexVector& Copy (ComplexVector& vector, double coefficient = 1.0);

  // copy a vector into another
  //
  // vector = vector to copy
  // coefficient = optional coefficient which multiply source to copy
  // return value = reference on current vector
  ComplexVector& Copy (ComplexVector& vector, const Complex& coefficient);

  // create a new vector with same size and same type but non-initialized components
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
  ComplexVector& operator - ();

  // return a new vector with opposite sign form a given source vector
  //
  // V1 = source vector
  // return value = new vector
  friend ComplexVector operator - (const ComplexVector& V1);

  // scalar product between two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = result of scalar product
  friend Complex operator * (const ComplexVector& V1, const ComplexVector& V2);

  // scalar product between two vectors
  //
  // V1 = first vector
  // V2 = second vector (real vector)
  // return value = result of scalar product
  friend Complex operator * (const ComplexVector& V1, const RealVector& V2);

  // scalar product between two vectors
  //
  // V1 = first vector (real vector)
  // V2 = second vector
  // return value = result of scalar product
  friend Complex operator * (const RealVector& V1, const ComplexVector& V2);

  // sum two vectors
  //
  // V1 = vector to add
  // return value = reference on current vector
  ComplexVector& operator += (const ComplexVector& V1);

  // sum two vectors
  //
  // V1 = real vector to add
  // return value = reference on current vector
  ComplexVector& operator += (const RealVector& V1);

  // sum two vectors
  //
  // vector = vector to add
  // return value = reference on current vector
  Vector& operator += (const Vector& vector);

  // add a linear combination to a given vector
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  ComplexVector& AddLinearCombination (double x, const ComplexVector& V);

  // add a linear combination to a given vector
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector  
  ComplexVector& AddLinearCombination (const Complex& x, const ComplexVector& V);

  // add a linear combination of two vectors to a given vector
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // return value = reference on current vector
  ComplexVector& AddLinearCombination (double x1, const ComplexVector& v1, double x2, 
				       const ComplexVector& v2);

  // add a linear combination of two vectors to a given vector
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // return value = reference on current vector
  ComplexVector& AddLinearCombination (const Complex& x1, const ComplexVector& v1, const Complex& x2, 
				       const ComplexVector& v2);
  // substract two vectors
  //
  // V1 = first vector
  // return value = reference on current vector
  ComplexVector& operator -= (const ComplexVector& V1);

  // substract two vectors
  //
  // V1 = first real vector
  // return value = reference on current vector
  ComplexVector& operator -= (const RealVector& V1);

  // sum two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = resulting vector
  friend ComplexVector operator + (const ComplexVector& V1, const ComplexVector& V2);

  // sum two vectors with left one real
  //
  // V1 = first vector (real)
  // V2 = second vector
  // return value = resulting vector
  friend ComplexVector operator + (const RealVector& V1, const ComplexVector& V2);

  // sum two vectors with right one real
  //
  // V1 = first vector
  // V2 = second vector (real)
  // return value = resulting vector
  friend ComplexVector operator + (const ComplexVector& V1, const RealVector& V2);

  // substract two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = resulting vector
  friend ComplexVector operator - (const ComplexVector& V1, const ComplexVector& V2);

  // substract two vectors with left one real
  //
  // V1 = first vector (real)
  // V2 = second vector
  // return value = resulting vector
  friend ComplexVector operator - (const RealVector& V1, const ComplexVector& V2);

  // substract two vectors with rightt one real
  //
  // V1 = first vector 
  // V2 = second vector (real)
  // return value = resulting vector
  friend ComplexVector operator - (const ComplexVector& V1, const RealVector& V2);

  // multiply a vector with a real number on the right hand side
  //
  // V1 = vector to multiply
  // d = real to use
  // return value = resulting vector
  friend ComplexVector operator * (const ComplexVector& V1, double d);

  // multiply a vector with a complex number on the right hand side
  //
  // V1 = vector to multiply
  // d = complex to use
  // return value = resulting vector
  friend ComplexVector operator * (const ComplexVector& V1, const Complex& d);

  // multiply a vector with a real number on the left hand side
  //
  // V1 = vector to multiply
  // d = real to use
  // return value = resulting vector
  friend ComplexVector operator * (double d, const ComplexVector& V1);

  // multiply a vector with a complex number on the left hand side
  //
  // V1 = vector to multiply
  // d = complex to use
  // return value = resulting vector
  friend ComplexVector operator * (const Complex& d, const ComplexVector& V1);

  // multiply a vector with a real number on the right hand side
  //
  // d = real to use
  // return value = reference on current vector
  ComplexVector& operator *= (double d);

  // divide a vector with a real number on the right hand side
  //
  // d = real to use
  // return value = reference on current vector  
  ComplexVector& operator /= (double d);

  // multiply a vector with a complex number on the right hand side
  //
  // d = complex to use
  // return value = reference on current vector
  ComplexVector& operator *= (const Complex& d);

  // divide a vector with a complex number on the right hand side
  //
  // d = complex to use
  // return value = reference on current vector  
  ComplexVector& operator /= (const Complex& d);

  // left multiply a vector with a real matrix (using temporary vector)
  //
  // M = matrix to use
  // return value = reference on current vector
  ComplexVector& operator *= (const RealMatrix&  M);

  // left multiply a vector with a complex matrix (using temporary vector)
  //
  // M = matrix to use
  // return value = reference on current vector
  ComplexVector& operator *= (const ComplexMatrix&  M);

  // left multiply a vector with an hermtian conjugated complex matrix (using temporary vector)
  //
  // M = matrix to use
  // return value = reference on current vector
  ComplexVector& operator &= (const ComplexMatrix&  M);

  // left multiply a vector with an hermitian matrix
  //
  // M = matrix to use
  // return value = reference on current vector
  ComplexVector& operator *= (const HermitianMatrix&  M);

  // left multiply a vector with a complex tridiagonal hermitian matrix (without using temporary vector)
  //
  // M = matrix to use
  // return value = reference on current vector
  ComplexVector& operator *= (const ComplexTriDiagonalHermitianMatrix&  M);

  // left multiply a vector with a real tridiagonal symmetric matrix (without using temporary vector)
  //
  // M = matrix to use
  // return value = reference on current vector
  ComplexVector& operator *= (const RealTriDiagonalSymmetricMatrix&  M);

  // left multiply a vector with an hermitian matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  ComplexVector& Multiply (const HermitianMatrix&  M, ComplexVector& V);

  // left multiply a vector with a real matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  ComplexVector& Multiply (const RealMatrix&  M, ComplexVector& V);

  // left multiply a vector with a complex matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  ComplexVector& Multiply (const ComplexMatrix&  M, ComplexVector& V);

  // left multiply a vector with a matrix and use to store result in current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  ComplexVector& Multiply (const Matrix&  M, ComplexVector& V);

  // left multiply a vector with a matrix and add result to current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  ComplexVector& AddMultiply (const Matrix&  M, ComplexVector& V);

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
  ComplexVector& Multiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceStep, 
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
  ComplexVector& AddMultiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceStep, 
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
  ComplexVector& Multiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceStep, 
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
  ComplexVector& AddMultiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceStep, 
			      int sourceNbrComponent, int destStart, int destStep);

  // return reference on vector i-th coordinate
  //
  // i = coordinate position
  Complex operator [] (int i);
  
  // return reference on real part of vector i-th coordinate
  //
  // i = coordinate position
  double& Re (int i);
  
  // return reference on imaginary part of vector i-th coordinate
  //
  // i = coordinate position
  double& Im (int i);
  
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
  ComplexVector& Normalize();
  
  // Extract a subvector from a given vector
  //
  // firstCoordinate = Coordinate where extraction has to begin
  // lastCoordinate = Coordinate where extraction has to stop (extract also include this last coordinate)
  // step = distance to the next coordinate (1 means to take the folowing)
  // return value = return corresponding subvector
  ComplexVector Extract(int firstCoordinate, int lastCoordinate, int step = 1);
  
  // Merge a subvector into a given vector
  //
  // V = vector to merge
  // firstCoordinate = Coordinate where merge has to begin
  // step = distance to the next coordinate in the destination vector (1 means to take the following)
  // return value = reference to the current Vector
  ComplexVector& Merge(const ComplexVector& V, int firstCoordinate, int step = 1);
  
  // Merge a real subvector into a complex given vector
  //
  // V = real vector to merge
  // firstCoordinate = Coordinate where merge has to begin
  // step = distance to the next coordinate in the destination vector (1 means to take the following)
  // return value = reference to the current Vector
  ComplexVector& Merge(const RealVector& V, int firstCoordinate, int step = 1);
  
  // Output Stream overload
  //
  friend ostream& operator << (ostream& Str, const ComplexVector& P);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // v = vector to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexVector& v);

};

// return reference on vector i-th coordinate
//
// i = coordinate position

inline Complex ComplexVector::operator [] (int i) 
{
  return Complex(this->RealComponents[i],  this->ImaginaryComponents[i]);
}
  
// return reference on real part of vector i-th coordinate
//
// i = coordinate position

inline double& ComplexVector::Re (int i)
{
  return this->RealComponents[i];
}

// return reference on imaginary part of vector i-th coordinate
//
// i = coordinate position

inline double& ComplexVector::Im (int i)
{
  return this->ImaginaryComponents[i];
}
  

#endif

