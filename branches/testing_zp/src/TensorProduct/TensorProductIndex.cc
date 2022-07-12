////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                     class for tensor product structure                     //
//                                                                            //
//                        last modification : 22/03/2001                      //
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


#include "TensorProduct/TensorProductIndex.h"


// default constructor
//

TensorProductIndex::TensorProductIndex() 
{
  this->Structure = TensorProductStructure();
  this->Indices = 0;
}
  
// constructor tensor product structure
//
// structure = reference on tensor product structure

TensorProductIndex::TensorProductIndex(const TensorProductStructure& structure) 
{
  this->Structure = structure;
  this->Indices = new int [this->Structure.GetNbrSpace()];
}
  
// copy constructor
//
// structure = reference on structure to copy

TensorProductIndex::TensorProductIndex(const TensorProductIndex& index) 
{
  this->Structure = index.Structure;
  if (index.Indices == 0)
    {
      this->Indices = 0;
    }
  else
    {
      this->Indices = new int [this->Structure.GetNbrSpace()];
      for (int i = 0; i < this->Structure.GetNbrSpace(); i++)
	this->Indices[i] = index.Indices[i];
    }
}
  
// destructor
//

TensorProductIndex::~TensorProductIndex() 
{
  if (this->Indices != 0)
    {
      delete[] this->Indices;
    }
}
  
// assignement
//
// structure = reference on structure to assign
// return value = reference on current structure

TensorProductIndex& TensorProductIndex::operator = (const TensorProductIndex& index) 
{
  if (this->Indices != 0)
    {
      delete[] this->Indices;
    }
  if (index.Indices == 0)
    {
      this->Indices = 0;
    }
  else
    {
      this->Indices = new int [this->Structure.GetNbrSpace()];
      for (int i = 0; i < this->Structure.GetNbrSpace(); i++)
	this->Indices[i] = index.Indices[i];
    }
  return *this;
}
