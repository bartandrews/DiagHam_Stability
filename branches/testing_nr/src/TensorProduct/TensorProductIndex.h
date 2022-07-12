////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                       class for tensor product index                       //
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


#ifndef TENSORPRODUCTINDEX_H
#define TENSORPRODUCTINDEX_H


#include "TensorProduct/TensorProductStructure.h"


class TensorProductIndex
{

 private:

  TensorProductStructure Structure;
  int* Indices;

 public:
  
  // default constructor
  //
  TensorProductIndex();
  
  // constructor tensor product structure
  //
  // structure = reference on tensor product structure
  TensorProductIndex(const TensorProductStructure& structure);
  
  // copy constructor
  //
  // structure = reference on structure to copy
  TensorProductIndex(const TensorProductIndex& index);
  
  // destructor
  //
  ~TensorProductIndex();
  
  // assignement
  //
  // structure = reference on structure to assign
  // return value = reference on current structure
  TensorProductIndex& operator = (const TensorProductIndex& index);

  // return index of a given space
  //
  // space = space index
  // return value = refrence on space dimension
  int& operator [] (int space);

  // return global index corresponding to current tensor product index
  //
  // return value = global index
  int GetGlobalIndex();

};


// return index of a given space
//
// space = space index
// return value = space dimension

inline int& TensorProductIndex::operator [] (int space) 
{
  return this->Indices[space];
}

// return global index corresponding to current tensor product index
//
// return value = global index

inline int TensorProductIndex::GetGlobalIndex()
{
  int pos = this->Indices[0];
  for (int i = 1; i < this->Structure.GetNbrSpace(); i++)
    {
      pos += this->Structure.GetIncrement(i) * this->Indices[i];
    }
  return pos;
}


#endif
