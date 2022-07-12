////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                               class of tensor                              //
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


#ifndef TENSOR_H
#define TENSOR_H


#include "Matrix/Matrix.h" 
#include "TensorProduct/TensorProductStructure.h"


class Tensor : public Matrix
{

 protected:

  AbstractTensorProductStructure* Structure;

 public:

  // virtual destructor
  //
  virtual ~Tensor ();

  // get number of row for a given space
  //
  // return value = number of row
  virtual int GetNbrRow (int space);

  // get number of column for a given space
  //
  // return value = number of column
  virtual int GetNbrColumn (int space);

  // get tensor product structure associated to the current tensor
  //
  // return value = reference on tensor product structure
  virtual AbstractTensorProductStructure* GetTensorProductStructure ();
  
  // Resize tensor
  //
  // structure = new product tensor structure
  virtual void Resize (const TensorProductStructure& structure) = 0;

  // Resize tensor and set to zero all components that have been added
  //
  // structure = new product tensor structure
  virtual void ResizeAndClean (const TensorProductStructure& structure) = 0;

};

#endif
