////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                     base class for tensor product vector                   //
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


#ifndef TENSORPRODUCTVECTOR_H
#define TENSORPRODUCTVECTOR_H


#include "config.h"
#include "Vector/Vector.h"
#include "TensorProduct/AbstractTensorProductStructure.h"


class TensorProductVector : public Vector
{
  
 protected:

  AbstractTensorProductStructure* Structure;
  
 public:

  // virtual destructor
  //
  virtual ~TensorProductVector ();

  // Get Vector dimension for a given space
  //
  // return value = dimension
  int GetVectorDimension(int space);
  
  // Resize vector
  //
  // structure = new product tensor structure
  virtual void Resize (AbstractTensorProductStructure* structure) = 0;

  // Resize vector and set to zero all components that have been added
  //
  // structure = new product tensor structure
  virtual void ResizeAndClean (AbstractTensorProductStructure* structure) = 0;

};


#endif

