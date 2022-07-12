////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                     class for tensor product real vector                   //
//                                                                            //
//                        last modification : 23/03/2001                      //
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


#ifndef TENSORPRODUCTREALVECTOR_H
#define TENSORPRODUCTREALVECTOR_H


#include "config.h"
#include "TensorProduct/TensorProductVector.h"
#include "TensorProduct/TensorProductIndex.h"
#include "Vector/RealVector.h"


#include <iostream>


class MathematicaOutput;
class RealSymmetricMatrix;
class RealTriDiagonalSymmetricMatrix;
class RealMatrix;
class OneSpaceTensor;
class TwoSpaceTensor;

class TensorProductRealVector : public TensorProductVector
{

 protected:
  
  RealVector GlobalVector;

 public:

  // default constructor
  //
  TensorProductRealVector();

  // constructor for an empty tensor product real vector 
  //
  // struture = reference on tensor product structure
  // zeroFlag = true if all coordinates have to be set to zero
  TensorProductRealVector(AbstractTensorProductStructure* structure, bool zeroFlag = false);

  // constructor from an array of doubles
  //
  // struture = reference on tensor product structure
  // array = array of doublesn
  TensorProductRealVector(AbstractTensorProductStructure* structure, double* array);

  // constructor from a real vector
  //
  // struture = reference on tensor product structure
  // v = reference on vector containing datas
  TensorProductRealVector(AbstractTensorProductStructure* structure, RealVector& v);

  // copy constructor
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  TensorProductRealVector(const TensorProductRealVector& vector, bool duplicateFlag = false);

  // destructor
  //
  ~TensorProductRealVector ();

  // assignement
  //
  // vector = vector to assign
  // return value = reference on current vector
  TensorProductRealVector& operator = (const TensorProductRealVector& vector);

  // Resize vector
  //
  // dimension = new dimension
  void Resize (int dimension);

  // Resize vector and set to zero all components that have been added
  //
  // dimension = new dimension
  void ResizeAndClean (int dimension);

  // Resize vector
  //
  // structure = new product tensor structure
  void Resize (AbstractTensorProductStructure* structure);

  // Resize vector and set to zero all components that have been added
  //
  // structure = new product tensor structure
  void ResizeAndClean (AbstractTensorProductStructure* structure);

  // change sign of a vector
  //
  // return value = reference on current vector
  TensorProductRealVector& operator - ();

  // return a new vector with opposite sign form a given source vector
  //
  // V1 = source vector
  // return value = new vector
  friend TensorProductRealVector operator - (TensorProductRealVector& V1);

  // scalar product between two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = result of scalar product
  friend double operator * (TensorProductRealVector& V1, TensorProductRealVector& V2);

  // sum two vectors
  //
  // V1 = vector to add
  // return value = reference on current vector
  TensorProductRealVector& operator += (TensorProductRealVector& V1);

  // substract two vectors
  //
  // V1 = first vector
  // return value = reference on current vector
  TensorProductRealVector& operator -= (TensorProductRealVector& V1);

  // sum two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = resulting vector
  friend TensorProductRealVector operator + (TensorProductRealVector& V1, TensorProductRealVector& V2);

  // substract two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = resulting vector
  friend TensorProductRealVector operator - (TensorProductRealVector& V1, TensorProductRealVector& V2);

  // add a linear combination to a given vector
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  TensorProductRealVector& AddLinearCombination (double x, TensorProductRealVector& V);

  // add a linear combination of two vectors to a given vector
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // return value = reference on current vector
  TensorProductRealVector& AddLinearCombination (double x1, TensorProductRealVector& v1, double x2, 
						 TensorProductRealVector& v2);

  // multiply a vector with a real number on the right hand side
  //
  // V1 = vector to multiply
  // d = real to use
  // return value = resulting vector
  friend TensorProductRealVector operator * (TensorProductRealVector& V1, double d);

  // multiply a vector with a real number on the left hand side
  //
  // V1 = vector to multiply
  // d = real to use
  // return value = resulting vector
  friend TensorProductRealVector operator * (double d, TensorProductRealVector& V1);

  // multiply a vector with a real number on the right hand side
  //
  // d = real to use
  // return value = reference on current vector
  TensorProductRealVector& operator *= (double d);

  // divide a vector by a real number on the right hand side
  //
  // d = real to use
  // return value = reference on current vector
  TensorProductRealVector& operator /= (double d);

  // left multiply a vector with a real symmetric matrix (without using temporary vector)
  //
  // M = matrix to use
  // return value = reference on current vector
  TensorProductRealVector& operator *= (const RealSymmetricMatrix&  M);

  // left multiply a vector with a real tridiagonal symmetric matrix (without using temporary vector)
  //
  // M = matrix to use
  // return value = reference on current vector
  TensorProductRealVector& operator *= (const RealTriDiagonalSymmetricMatrix&  M);

  // left multiply a vector with a symmetric matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply  
  // return value = reference on current vector
  TensorProductRealVector& Multiply (const RealSymmetricMatrix&  M, TensorProductRealVector& V);

  // left multiply a vector with a symmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply  
  // start = vector first coordinate to modify
  // step = step to add to go to the following vector coordinate
  // return value = reference on current vector
  TensorProductRealVector& AddMultiply (const RealSymmetricMatrix&  M, TensorProductRealVector& V, int start, int step);

  // left multiply a vector with a real matrix and use to store result 
  // in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  TensorProductRealVector& Multiply (const RealMatrix&  M, TensorProductRealVector& V);

  // left multiply a vector with a one space tensor and use to store result 
  // in current vector (without creating temporary vector)
  //
  // T = tensor to use
  // V = vector to multiply  
  // return value = reference on current vector
  TensorProductRealVector& Multiply (const OneSpaceTensor& T, TensorProductRealVector& V);

  // left multiply a vector with a one space tensor for a given range of indices 
  // and use to store result in current vector (without creating temporary vector)
  //
  // T = tensor to use
  // V = vector to multiply  
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on current vector
  TensorProductRealVector& Multiply (const OneSpaceTensor& T, TensorProductRealVector& V,
				     int firstComponent, int nbrComponent);

  // left multiply a vector with a one space tensor and 
  // add result to current vector (without creating temporary vector)
  //
  // T = tensor to use
  // V = vector to multiply  
  // return value = reference on current vector
  TensorProductRealVector& AddMultiply (const OneSpaceTensor& T, TensorProductRealVector& V);

  // left multiply a vector with a one space tensor for a given range of indices 
  // and add result to current vector (without creating temporary vector)
  //
  // T = tensor to use
  // V = vector to multiply  
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on current vector
  TensorProductRealVector& AddMultiply (const OneSpaceTensor& T, TensorProductRealVector& V,
					int firstComponent, int nbrComponent);

  // left multiply a vector with a two space tensor and 
  // use to store result in current vector (without creating temporary vector)
  //
  // T = tensor to use
  // V = vector to multiply  
  // return value = reference on current vector
  TensorProductRealVector& Multiply (const TwoSpaceTensor& T, TensorProductRealVector& V);

  // return vector i-th global coordinate (without testing if position is valid)
  //
  // i = global coordinate position
  // return value = reference on indexed coordinate
  double& operator [] (int i);

  // return vector coordinate corresponding to a tensor product index
  //
  // index = tensor product index
  // return value = reference on indexed coordinate
  double& operator [] (TensorProductIndex& index);

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
  TensorProductRealVector& Normalize();

  
};


#endif

