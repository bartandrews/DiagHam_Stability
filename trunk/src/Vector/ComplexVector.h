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
#include "MathTools/Complex.h"
#include "GeneralTools/GarbageFlag.h"


#include <iostream>


using std::ostream;


#ifdef USE_OUTPUT
class MathematicaOutput;
#endif
class RealMatrix;
class ComplexMatrix;
class ComplexDiagonalMatrix;

class ComplexVector : public Vector
{

  friend class RealVector;
  friend class DelocalizedRealVector;
  friend class ComplexMatrix;
  friend class HermitianMatrix;
  friend class ComplexSkewSymmetricMatrix;
  friend class ComplexUpperTriangularMatrix;
  friend class ComplexLowerTriangularMatrix;
  friend class RealTriDiagonalSymmetricMatrix;
  friend class ComplexTriDiagonalHermitianMatrix;
  friend class ComplexDiagonalMatrix;
  friend ComplexMatrix operator + (const ComplexMatrix& M1, const ComplexMatrix& M2);
  friend ComplexMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const ComplexMatrix& M2);
  friend ComplexMatrix operator + (const ComplexMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);
  friend ComplexMatrix operator - (const ComplexMatrix& M1, const ComplexMatrix& M2);
  friend ComplexMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, const ComplexMatrix& M2);
  friend ComplexMatrix operator - (const ComplexMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);
  friend ComplexMatrix operator * (const ComplexMatrix& M1, const ComplexMatrix& M2);
  friend ComplexMatrix operator * (const ComplexMatrix& M1, const RealDiagonalMatrix& M2);
  friend ComplexMatrix operator * (const RealDiagonalMatrix& M1, const ComplexMatrix& M2);
  friend ComplexMatrix operator * (const ComplexMatrix& M1, const ComplexDiagonalMatrix& M2);
  friend ComplexMatrix operator * (const ComplexDiagonalMatrix& M1, const ComplexMatrix& M2);

  friend ComplexMatrix operator * (const ComplexMatrix& M, double x);
  friend ComplexMatrix operator * (double x, const ComplexMatrix& M);
  friend ComplexMatrix operator / (const ComplexMatrix& M, double x);
  friend ComplexMatrix operator / (const ComplexMatrix& M1, const RealDiagonalMatrix& M2);
  friend ostream& operator << (ostream& Str, const ComplexMatrix& P);
#ifdef USE_OUTPUT
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexMatrix& P);
#endif
 protected:
  
  Complex* Components;
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

  
  // constructor for an empty real vector bigger than 2^31
  //
  // size = Vector Dimension 
  // zeroFlag = true if all coordinates have to be set to zero
  ComplexVector(long size, bool zeroFlag);

  
  // constructor from arrays of doubles
  //
  // real = array of doubles corresponding to real part
  // imaginary = array of doubles corresponding to imaginary part
  // size = Vector Dimension 
  ComplexVector(double* real, double* imaginary, int size);

  // constructor from large arrays of doubles
  //
  // real = array of doubles corresponding to real part
  // imaginary = array of doubles corresponding to imaginary part
  // size = Vector Dimension 
  ComplexVector(double* real, double* imaginary, long size);

  // constructor from Complex array
  //
  // components = array of Complex values
  // size = Vector Dimension 
  ComplexVector(Complex *components, int size);

  // constructor from large Complex array
  //
  // components = array of Complex values
  // size = Vector Dimension 
  ComplexVector(Complex *components, long size);

  // copy constructor
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  ComplexVector(const ComplexVector& vector, bool duplicateFlag = false);

  // copy constructor from a real vector
  //
  // vector = vector to copy
  ComplexVector(const RealVector& vector);

#ifdef __MPI__
  // constructor from informations sent using MPI
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts or sends the vector
  // broadcast = true if the vector is broadcasted
  ComplexVector(MPI::Intracomm& communicator, int id, bool broadcast = true);
#endif

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

  // Resize long vector
  //
  // dimension = new dimension
  
  void Resize (long dimension);

  // Resize vector and set to zero all components that have been added
  //
  // dimension = new dimension
  void ResizeAndClean (int dimension);

  // Resize long vector and set to zero all components that have been added
  //
  // dimension = new dimension
  void ResizeAndClean (long dimension);

  // copy a vector into another
  //
  // vector = vector to copy
  // return value = reference on current vector
  ComplexVector& Copy (ComplexVector& vector);

  // copy a vector into another
  //
  // vector = vector to copy
  // coefficient = optional coefficient which multiply source to copy
  // return value = reference on current vector
  ComplexVector& Copy (ComplexVector& vector, double coefficient);

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

  // create an array of new vectors with same size and same type but non-initialized components
  //
  // nbrVectors = number of vectors to sreate
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to the array of new vectors
  Vector* EmptyCloneArray(int nbrVectors, bool zeroFlag = false);

  // put all vector components to zero
  //
  // return value = reference on current vector
  Vector& ClearVector ();

  // put select vector components to zero
  //
  // start = start index
  // nbrComponent = number of components to set to zero
  // return value = reference on current vector
  Vector& ClearVectorSegment (long start, long nbrComponent);

  // apply standard phase conventions for this vector
  //
  // maxIndex = component with maximum amplitude that was set to real
  ComplexVector& SetStandardPhase(long &maxIndex);

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

  // do part of the scalar product between two vectors in a given range of indices
  //
  // vRight = right vector of the scalar product
  // firstComponent = index of the first component to consider
  // nbrComponent = number of components to consider
  // step = increment between to consecutive indices to consider
  // return value = result of the partial scalar product
  Complex PartialScalarProduct (const ComplexVector& vRight, int firstComponent, int nbrComponent, int step = 1);

  // Euclidian scalar product between two vectors (i.e. \sum_i V1[i] * V2[i])
  //
  // V1 = first vector
  // V2 = second vector
  // return value = result of scalar product
  friend Complex EuclidianScalarProduct (const ComplexVector& V1, const ComplexVector& V2);

  // assuming the vector is real up to a global phase, compute this global phase
  //
  // return value = global phase factor (i.e. exp(i phi))
  Complex GlobalPhase();

  // sum two vectors
  //
  // V1 = vector to add
  // return value = reference on current vector
  ComplexVector& operator += (ComplexVector& V1);

  // sum two vectors
  //
  // V1 = real vector to add
  // return value = reference on current vector
  ComplexVector& operator += (RealVector& V1);

  // sum two vectors
  //
  // vector = vector to add
  // return value = reference on current vector
  Vector& operator += (Vector& vector);

  // add a linear combination to a given vector
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  ComplexVector& AddLinearCombination (double x, const ComplexVector& V);

  // add a linear combination to a given vector, for a given range of indices
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  ComplexVector& AddLinearCombination (double x, const ComplexVector& V, int firstComponent, int nbrComponent);

  // add a linear combination to a given vector
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector  
  ComplexVector& AddLinearCombination (const Complex& x, const ComplexVector& V);

  // add a linear combination to a given vector, for a given range of indices
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  ComplexVector& AddLinearCombination (const Complex& x, const ComplexVector& V, int firstComponent, int nbrComponent);

  // add a linear combination to a given vector, for a given range of indices
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  ComplexVector& AddLinearCombination (const Complex& x, const RealVector& V, int firstComponent, int nbrComponent);

  // add a linear combination of two vectors to a given vector
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // return value = reference on current vector
  ComplexVector& AddLinearCombination (double x1, const ComplexVector& v1, double x2, 
				       const ComplexVector& v2);

  // add a linear combination of two vectors to a given vector, for a given range of indices
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on current vector
  ComplexVector& AddLinearCombination (double x1, const ComplexVector& v1, double x2, 
				       const ComplexVector& v2, int firstComponent, int nbrComponent);

  // add a linear combination of two vectors to a given vector
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // return value = reference on current vector
  ComplexVector& AddLinearCombination (const Complex& x1, const ComplexVector& v1, const Complex& x2, 
				       const ComplexVector& v2);

  // add a linear combination of two vectors to a given vector, for a given range of indices
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on current vector
  ComplexVector& AddLinearCombination (const Complex& x1, const ComplexVector& v1, const Complex& x2, 
				       const ComplexVector& v2, int firstComponent, int nbrComponent);

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
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  ComplexVector& Multiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with a matrix and add current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  ComplexVector& AddMultiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceNbrComponent);

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

  // left multiply a vector with an hermitian matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  ComplexVector& AddMultiply (const HermitianMatrix&  M, ComplexVector& V);

  // do a partial left multication of a vector with an hermitian matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  ComplexVector&Multiply (const HermitianMatrix&  M, ComplexVector& V, int sourceStart, int sourceNbrComponent);

  // do a partial left multication of a vector with an hermitian matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  ComplexVector& AddMultiply (const HermitianMatrix&  M, ComplexVector& V, int sourceStart, int sourceNbrComponent);

  // return reference on vector i-th coordinate
  //
  // i = coordinate position
  virtual Complex& operator [] (int i);
  
  // return reference on vector i-th coordinate
  //
  // i = coordinate position
  virtual Complex& operator [] (long i);
  
  // return reference on real part of vector i-th coordinate
  //
  // i = coordinate position
  virtual double& Re (int i);
  
  // return reference on real part of vector i-th coordinate
  //
  // i = coordinate position
  virtual double& Re (long i);
  
  // return reference on imaginary part of vector i-th coordinate
  //
  // i = coordinate position
  virtual double& Im (int i);
  
  // return reference on imaginary part of vector i-th coordinate
  //
  // i = coordinate position
  virtual double& Im (long i);
  
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
  
  // conjugate the vector
  //
  // return value = reference on current vector
  ComplexVector& Conjugate();
  
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

  // reverse elements of the current vector (i.e. exchanging i <-> N - i)
  //
  // return value = reference to the current Vector
  ComplexVector& ReverseVector();

  
  // write vector in a file 
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool ByteWriteVector (const char* fileName);

  // write vector in a file 
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool WriteVector (const char* fileName);


  // write vector in a file in ascii mode
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool WriteAsciiVector (const char* fileName);

  // read vector from a file 
  //
  // fileName = name of the file where the vector has to be read
  // return value = true if no error occurs
  bool ByteReadVector (const char* fileName);

  // blocked version of read
  bool ReadVector (const char* fileName);
    

  // read vector from a file, only within a given range of indices
  //
  // fileName = name of the file where the vector has to be read
  // minIndex = index of the first component to read (if negative, start from the end of vector)
  // maxIndex = index of the last component to read (negative or zero is considered as the last component)
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
  friend ostream& operator << (ostream& Str, const ComplexVector& P);

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

  // compare if any entries differ between the current and a reference vector
  //
  // V = reference vector
  // threshold = threshold for reporting of differences
  // str = stream to log output
  // return value = reference on output stream
  ostream& CompareVector(ComplexVector& V, double threshold=1e-13, ostream &str = std::cout);

#ifdef USE_OUTPUT
  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // v = vector to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexVector& v);
#endif

#ifdef __MPI__

  // send a vector to a given MPI process
  // 
  // communicator = reference on the communicator to use
  // id = id of the destination MPI process
  // return value = reference on the current vector
  Vector& SendVector(MPI::Intracomm& communicator, int id);

  // broadcast a vector to all MPI processes associated to the same communicator
  // 
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // return value = reference on the current vector
  Vector& BroadcastVector(MPI::Intracomm& communicator,  int id);

  // broadcast part of vector to all MPI processes associated to the same communicator
  // 
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // firstComponent = index of the first component (useless if the method is not called by the MPI process which broadcasts the vector)
  // nbrComponent = number of component (useless if the method is not called by the MPI process which broadcasts the vector)
  // return value = reference on the current vector
  Vector& BroadcastPartialVector(MPI::Intracomm& communicator, int id, int firstComponent = 0, int nbrComponent = 0);

  // receive a vector from a MPI process
  // 
  // communicator = reference on the communicator to use 
  // id = id of the source MPI process
  // return value = reference on the current vector
  Vector& ReceiveVector(MPI::Intracomm& communicator, int id);

  // add current vector to the current vector of a given MPI process
  // 
  // communicator = reference on the communicator to use 
  // id = id of the destination MPI process
  // return value = reference on the current vector
  Vector& SumVector(MPI::Intracomm& communicator, int id);

  // reassemble vector from a scattered one
  // 
  // communicator = reference on the communicator to use 
  // id = id of the destination MPI process
  // return value = reference on the current vector
  Vector& ReassembleVector(MPI::Intracomm& communicator, int id);

  // create a new vector on each MPI node which is an exact clone of the broadcasted one
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  Vector* BroadcastClone(MPI::Intracomm& communicator, int id);

  // create a new vector on given MPI node which is an exact clone of the sent one but with only part of the data
  // 
  // communicator = reference on the communicator to use
  // id = id of the destination MPI process
  // firstComponent = index of the first component 
  // nbrComponent = number of component to send
  // return value = reference on the current vector
  Vector& SendPartialClone(MPI::Intracomm& communicator, int id, int firstComponent, int nbrComponent);

  // scatter this vector across all MPI nodes with the given load balancing information
  // 
  // communicator = reference on the communicator to use
  // mininumIndices = lowest index for each thread
  // maximumIndices = largest index for each thread
  // id = id of the process to send the vector
  // return value = reference on the current vector
  Vector& ScatterPartialClones(MPI::Intracomm& communicator, long *mininumIndices, long *maximumIndices, int id);

  // create a new vector on each MPI node with same size and same type but non-initialized components
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  Vector* BroadcastEmptyClone(MPI::Intracomm& communicator, int id, bool zeroFlag = false);

#endif

};

// return reference on vector i-th coordinate
//
// i = coordinate position

inline Complex& ComplexVector::operator [] (int i) 
{
  return this->Components[i];
}
  
// return reference on vector i-th coordinate
//
// i = coordinate position

inline Complex& ComplexVector::operator [] (long i) 
{
  return this->Components[i];
}
  
// return reference on real part of vector i-th coordinate
//
// i = coordinate position

inline double& ComplexVector::Re (int i)
{
  return (this->Components[i].Re);
}

// return reference on real part of vector i-th coordinate
//
// i = coordinate position

inline double& ComplexVector::Re (long i)
{
  return (this->Components[i].Re);
}

// return reference on imaginary part of vector i-th coordinate
//
// i = coordinate position

inline double& ComplexVector::Im (int i)
{
  return (this->Components[i].Im);
}
  
// return reference on imaginary part of vector i-th coordinate
//
// i = coordinate position

inline double& ComplexVector::Im (long i)
{
  return (this->Components[i].Im);
}
  
// print a single component
//
// str = reference on output stream
// index = index of the component to print
// return value  = reference on output stream

inline ostream& ComplexVector::PrintComponent(ostream& str, long index)
{
  str << this->Components[index];
  return str;
}


#endif

