////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                 class for abstract tensor product structure                //
//                                                                            //
//                        last modification : 08/06/2001                      //
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


#ifndef ABSTRACTTENSORPRODUCTSTRUCTURE_H
#define ABSTRACTTENSORPRODUCTSTRUCTURE_H


#include "config.h"

#include <iostream>


using std::ostream;


class AbstractTensorProductStructure
{

 protected:

  int NbrSpace;
  int StructureID;

 public:
  
  enum StructureTypes
  {
    Simple = 0x0001,
    Composite = 0x0002,
    Full = 0x0004
  };

  // virtual destructor
  //
  virtual ~AbstractTensorProductStructure();
  
  // get tensor product structure type
  //
  // return value = tensor product structure type
  virtual int GetTensorProductStructureType();

  // return dimension of a given space
  //
  // space = space index
  // return value = space dimension
  virtual int GetDimension(int space) = 0;

  // return total space dimension
  //
  // return value = dimension
  virtual int GetTotalDimension() = 0;

  // return number of spaces
  //
  // return value = number of spaces
  virtual int GetNbrSpace();

  // return tensor product structure of a space obtained by tensor product of two tensor spaces
  //
  // structure1 = first structure
  // structure2 = second structure
  // return value = structure corresponding to the new tensor space
  friend AbstractTensorProductStructure* operator * (const AbstractTensorProductStructure& structure1, 
						     const AbstractTensorProductStructure& structure2);

  // test if two tensor product structures are equivalent
  //
  // structure1 = first structure
  // structure2 = second structure
  // return value = true if structures are equivalent
  friend bool operator == (const AbstractTensorProductStructure& structure1, 
			   const AbstractTensorProductStructure& structure2);

  // test if two tensor product structures are different
  //
  // structure1 = first structure
  // structure2 = second structure
  // return value = true if structures are different
  friend bool operator != (const AbstractTensorProductStructure& structure1, 
			   const AbstractTensorProductStructure& structure2);

  // print information on current tensor product structure
  //
  // str = output stream
  // structure = structure to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, const AbstractTensorProductStructure& structure);

};

// get tensor product structure type
//
// return value = tensor product structure type

inline int AbstractTensorProductStructure::GetTensorProductStructureType()
{
  return this->StructureID;
}

// return number of spaces
//
// return value = number of spaces

inline int AbstractTensorProductStructure::GetNbrSpace()
{
  return this->NbrSpace;
}

#endif
