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


#include "Tensor/OneSpaceTensor.h"


// default constructor
//

OneSpaceTensor::OneSpaceTensor() 
{
  this->Structure = 0;
  this->ElementaryMatrix = 0;
  this->TargetSpace = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = 0;
  this->TrueNbrColumn = 0;
  this->MatrixType = 0;
}

// constructor from standard datas
//
// struture = reference on tensor product structure
// elementaryMatrix = reference on matrix to use
// targetSpace = space on which matrix acts

OneSpaceTensor::OneSpaceTensor(AbstractTensorProductStructure* structure, Matrix& elementaryMatrix, 
			       int targetSpace) 
{
  this->Structure = structure;
  this->TargetSpace = targetSpace;
  this->ElementaryMatrix = elementaryMatrix.Clone();
  this->NbrRow = this->Structure->GetDimension(this->TargetSpace);
  this->NbrColumn = this->NbrRow;
  this->TrueNbrRow = this->Structure->GetTotalDimension();
  this->TrueNbrColumn = this->TrueNbrRow;
  this->MatrixType = this->ElementaryMatrix->GetMatrixType();
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

OneSpaceTensor::OneSpaceTensor(const OneSpaceTensor& T) 
{
  this->Structure = T.Structure;
  this->TargetSpace = T.TargetSpace;
  this->ElementaryMatrix = T.ElementaryMatrix->Clone();
  this->NbrRow = T.NbrRow;
  this->NbrColumn = T.NbrColumn;
  this->TrueNbrRow = T.TrueNbrRow;
  this->TrueNbrColumn = T.TrueNbrColumn;
  this->MatrixType = T.MatrixType;
}

// destructor
//

OneSpaceTensor::~OneSpaceTensor() 
{
  if (this->ElementaryMatrix != 0)
    delete this->ElementaryMatrix;
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

OneSpaceTensor& OneSpaceTensor::operator = (const OneSpaceTensor& T) 
{
  this->Structure = T.Structure;
  this->TargetSpace = T.TargetSpace;
  this->ElementaryMatrix = T.ElementaryMatrix->Clone();
  this->NbrRow = T.NbrRow;
  this->NbrColumn = T.NbrColumn;
  this->TrueNbrRow = T.TrueNbrRow;
  this->TrueNbrColumn = T.TrueNbrColumn;
  this->MatrixType = T.MatrixType;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* OneSpaceTensor::Clone ()
{
  return (Matrix*) new OneSpaceTensor (*this);
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void OneSpaceTensor::SetMatrixElement(int i, int j, double x) 
{
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void OneSpaceTensor::SetMatrixElement(int i, int j, const Complex& x) 
{
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void OneSpaceTensor::AddToMatrixElement(int i, int j, double x) 
{
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void OneSpaceTensor::AddToMatrixElement(int i, int j, const Complex& x) 
{
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void OneSpaceTensor::Resize (int nbrRow, int nbrColumn) 
{
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void OneSpaceTensor::ResizeAndClean (int nbrRow, int nbrColumn) 
{
}

// Resize tensor
//
// structure = new product tensor structure

void OneSpaceTensor::Resize (const TensorProductStructure& structure) 
{
}

// Resize tensor and set to zero all components that have been added
//
// structure = new product tensor structure

void OneSpaceTensor::ResizeAndClean (const TensorProductStructure& structure) 
{
}

// evaluate matrix trace
//
// return value = matrix trace 

double OneSpaceTensor::Tr () 
{
  return ((this->ElementaryMatrix->Tr () * (double) this->Structure->GetTotalDimension()) / 
	  (double) this->Structure->GetDimension(this->TargetSpace));
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double OneSpaceTensor::Det () 
{
  return 1.0;
}
