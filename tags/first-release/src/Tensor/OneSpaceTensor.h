////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of tensor acting on one space                   //
//                                                                            //
//                        last modification : 30/03/2001                      //
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


#ifndef ONESPACETENSOR_H
#define ONESPACETENSOR_H


#include "Matrix/Matrix.h" 
#include "Tensor/Tensor.h"


class OneSpaceTensor;


class OneSpaceTensor : public Tensor
{

  friend class TensorProductRealVector;
  friend class TwoSpaceTensor;

 protected:

  Matrix* ElementaryMatrix;
  int TargetSpace;

 public:

  // default constructor
  //
  OneSpaceTensor();

  // constructor from standard datas
  //
  // struture = reference on tensor product structure
  // elementaryMatrix = reference on matrix to use
  // targetSpace = space on which matrix acts
  OneSpaceTensor(AbstractTensorProductStructure* structure, Matrix& elementaryMatrix, int targetSpace);

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  OneSpaceTensor(const OneSpaceTensor& M);

  // destructor
  //
  ~OneSpaceTensor();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  OneSpaceTensor& operator = (const OneSpaceTensor& M);

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

  // evaluate matrix trace
  //
  // return value = matrix trace 
  double Tr ();

  // evaluate matrix determinant
  //
  // return value = matrix determinant 
  double Det ();

};

#endif
