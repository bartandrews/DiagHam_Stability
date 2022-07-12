////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of tensor acting on two consecutive spaces            //
//                                                                            //
//                        last modification : 03/04/2001                      //
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


#ifndef TWOSPACETENSOR_H
#define TWOSPACETENSOR_H


#include "Matrix/Matrix.h" 
#include "Tensor/Tensor.h"


class OneSpaceTensor;
class RealDiagonalMatrix;
class RealSymmetricMatrix;
class RealAntisymmetricMatrix;
class RealTriDiagonalSymmetricMatrix;
class SpaceDecomposition;


class TwoSpaceTensor : public Tensor
{

  friend class TensorProductRealVector;

 protected:

  int FirstTargetSpace;

 public:

  Matrix* ElementaryMatrix;

  // default constructor
  //
  TwoSpaceTensor();

  // constructor from standard datas
  //
  // struture = reference on tensor product structure
  // elementaryMatrix = matrix to use
  // firstTargetSpace = first space on which matrix acts
  TwoSpaceTensor(AbstractTensorProductStructure* structure, Matrix& elementaryMatrix, int firstTargetSpace);

  // constructor from the tensor product of two matrices
  //
  // struture = reference on tensor product structure
  // m1 = pointer to the matrix to use for first target space
  // m2 = pointer to the matrix to use for second target space
  // firstTargetSpace = first space on which matrix acts
  // coef = global multiplicative factor
  TwoSpaceTensor(AbstractTensorProductStructure* structure, Matrix* m1, Matrix* m2, int firstTargetSpace, 
		 double coef = 1.0);

  // constructor from the tensor product of two matrices with a given space decomposition
  //
  // struture = reference on tensor product structure
  // m1 = reference on matrix to use for first target space
  // m2 = reference on matrix to use for second target space
  // decomposition1 = first space decomposition
  // decomposition2 = second space decomposition
  // firstTargetSpace = first space on which matrix acts
  // coef = global multiplicative factor
  TwoSpaceTensor(AbstractTensorProductStructure* structure, RealSymmetricMatrix& m1, 
		 RealSymmetricMatrix& m2, SpaceDecomposition& decomposition1,
		 SpaceDecomposition& decomposition2, int firstTargetSpace, double coef = 1.0);

  // constructor from the tensor product of two matrices with a given space decomposition
  //
  // struture = reference on tensor product structure
  // m1 = reference on matrix to use for first target space
  // m2 = reference on matrix to use for second target space
  // decomposition1 = first space decomposition
  // decomposition2 = second space decomposition
  // firstTargetSpace = first space on which matrix acts
  // coef = global multiplicative factor
  TwoSpaceTensor(AbstractTensorProductStructure* structure, Matrix& m1, 
		 Matrix& m2, SpaceDecomposition& decomposition1,
		 SpaceDecomposition& decomposition2, int firstTargetSpace, double coef = 1.0);

  // constructor from a given matrix assuming the other is identity (with a space index permutation)
  //
  // struture = reference on tensor product structure
  // m = reference on matrix to use
  // permutation = two dimensional array describing permutation
  // firstTargetSpace = first space on which matrix acts
  // coef = global multiplicative factor
  TwoSpaceTensor(AbstractTensorProductStructure* structure, Matrix& m, 
		 int** permutation, int firstTargetSpace, double coef = 1.0);

  // constructor from a given matrix asociated to a given space 
  // and array containing positions where elements are to be put (assuming identity on other spaces)
  //
  // struture = reference on tensor product structure
  // m = reference on matrix to use
  // totalSpaceDimension = dimension of the total space
  // arraySize = size of used array
  // subspaceSize = indicate size of each subspace associated to a fixed coordinate 
  // for spaces where identity acts
  // spaceIndex = array containing index in space where matrix acts for each fixed coordinate
  // globalSpaceIndex = array containing global index in space where matrix acts for each fixed coordinate
  // targetSpace = target space index
  // coef = global multiplicative factor
  TwoSpaceTensor(AbstractTensorProductStructure* structure, Matrix& m, 
		 int totalSpaceDimension, int arraySize, int* subspaceSize,
		 int** spaceIndex, int** globalSpaceIndex, int targetSpace, double coef = 1.0);

  // constructor from the product of two one space tensors
  //
  // t1 = first one space tensor
  // t2 = second one space tensor
  // coef = global multiplicative factor
  TwoSpaceTensor(const OneSpaceTensor& t1, const OneSpaceTensor& t2, double coef);
  
  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  TwoSpaceTensor(const TwoSpaceTensor& M);

  // destructor
  //
  ~TwoSpaceTensor();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  TwoSpaceTensor& operator = (const TwoSpaceTensor& M);

  // return pointer on a clone matrix (without duplicating datas)
  //
  // retrun value = pointer on new matrix 
  Matrix* Clone ();

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  void SetMatrixElement(int i, int j, double x);

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  void SetMatrixElement(int i, int j, const Complex& x);

  // add a value to a matrix element
  //
  // i = line position
  // j = column position
  // x = value to add to matrix element
  void AddToMatrixElement(int i, int j, double x);

  // add a value  a matrix element
  //
  // i = line position
  // j = column position
  // x = value to add to matrix element
  void AddToMatrixElement(int i, int j, const Complex& x);

  // Resize matrix
  //
  // nbrRow = new number of rows
  // nbrColumn = new number of columns
  void Resize (int nbrRow, int nbrColumn);

  // Resize matrix and set to zero all elements that have been added
  //
  // nbrRow = new number of rows
  // nbrColumn = new number of columns
  void ResizeAndClean (int nbrRow, int nbrColumn);

  // Resize tensor
  //
  // structure = new product tensor structure
  void Resize (const TensorProductStructure& structure);

  // Resize tensor and set to zero all components that have been added
  //
  // structure = new product tensor structure
  void ResizeAndClean (const TensorProductStructure& structure);

  // Add product of two one space tensors to the current tensor
  //
  // t1 = first one space tensor
  // t2 = second one space tensor
  // return value = reference on current tensor
  TwoSpaceTensor& AddProduct (const OneSpaceTensor& t1, const OneSpaceTensor& t2);

  // evaluate matrix trace
  //
  // return value = matrix trace 
  double Tr ();

  // evaluate matrix determinant
  //
  // return value = matrix determinant 
  double Det ();

  // add a one space tensor to the current two space tensor
  //
  // tensor = one space tensor to add 
  // return value = reference to current tensor
  TwoSpaceTensor& operator += (const OneSpaceTensor& tensor);
 
  // add tensor product of two matrices to the current tensor (assuming identical space tensor pattern)
  //
  // m1 = pointer to the matrix associated to first target space
  // m2 = pointer to the matrix associated to first target space
  // coef = global multiplicative factor
  // return value = reference to current tensor
  TwoSpaceTensor& AddTensorProductMatrices (Matrix* m1, Matrix* m2, double coef = 1.0);

  // add tensor product of two matrices to the current tensor with a given space decomposition
  //
  // m1 = reference on the matrix associated to first target space
  // m2 = reference on  the matrix associated to second  target space
  // decomposition1 = first space decomposition
  // decomposition2 = second space decomposition
  // coef = global multiplicative factor
  // return value = reference to current tensor
  TwoSpaceTensor& AddTensorProductMatrices(Matrix& m1, Matrix& m2, 
					   SpaceDecomposition& decomposition1, 
					   SpaceDecomposition& decomposition2, double coef = 1.0);

  // add tensor product of two matrices to the current tensor with a given space decomposition
  //
  // m1 = reference on the matrix associated to first target space
  // m2 = reference on  the matrix associated to second  target space
  // decomposition1 = first space decomposition
  // decomposition2 = second space decomposition
  // coef = global multiplicative factor
  // return value = reference to current tensor
  TwoSpaceTensor& AddTensorProductMatrices(RealSymmetricMatrix& m1, RealSymmetricMatrix& m2, 
					   SpaceDecomposition& decomposition1, 
					   SpaceDecomposition& decomposition2, double coef = 1.0);

  // add tensor product of two matrices to the current tensor with a given space decomposition
  //
  // m1 = reference on the matrix associated to first target space
  // m2 = reference on  the matrix associated to second  target space
  // decomposition1 = first space decomposition
  // decomposition2 = second space decomposition
  // coef = global multiplicative factor
  // return value = reference to current tensor
  TwoSpaceTensor& AddTensorProductMatrices(RealAntisymmetricMatrix& m1, RealAntisymmetricMatrix& m2, 
					   SpaceDecomposition& decomposition1, 
					   SpaceDecomposition& decomposition2, double coef = 1.0);

  // substract tensor product of two matrices to the cuurent tensor 
  // (assuming identical space tensor pattern)
  //
  // m1 = pointer to the matrix associated to first target space
  // m2 = pointer to the matrix associated to first target space
  // return value = reference to current tensor
  TwoSpaceTensor& SubTensorProductMatrices (Matrix* m1, Matrix* m2);

  // add a real diagonal matrix to matrix of first target space
  //
  // matrix = matrix to add
  void AddToFirstSpace(RealDiagonalMatrix& matrix);
  
  // add a matrix to matrix of first target space
  //
  // matrix = matrix to add
  void AddToFirstSpace(Matrix& matrix);
  
  // add a real symmetric matrix to matrix of first target space
  //
  // matrix = matrix to add
  void AddToFirstSpace(RealSymmetricMatrix& matrix);
  
  void AddToFirstSpace(RealDiagonalMatrix& matrix, SpaceDecomposition& decomposition1, 
		       SpaceDecomposition& decomposition2);

  // add a real diagonal matrix to matrix of second target space
  //
  // matrix = matrix to add
  void AddToSecondSpace(RealDiagonalMatrix& matrix);

  // add a matrix to matrix of second target space
  //
  // matrix = matrix to add
  void AddToSecondSpace(Matrix& matrix);
  
  // add a real symmetric matrix to matrix of second target space
  //
  // matrix = matrix to add
  void AddToSecondSpace(RealSymmetricMatrix& matrix);
  
  void AddToSecondSpace(RealDiagonalMatrix& matrix, SpaceDecomposition& decomposition1, 
			SpaceDecomposition& decomposition2);

  // add a real tridiagonal symmetric matrix to matrix of first target space
  //
  // matrix = matrix to add
  void AddToFirstSpace(RealTriDiagonalSymmetricMatrix& matrix);

  // add a real tridiagonal symmetric matrix to matrix of second target space
  //
  // matrix = matrix to add
  void AddToSecondSpace(RealTriDiagonalSymmetricMatrix& matrix);

  // Output Stream overload
  //
  // str = reference on output stream
  // tensor = tensor to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, const TwoSpaceTensor& tensor);

};

#endif
