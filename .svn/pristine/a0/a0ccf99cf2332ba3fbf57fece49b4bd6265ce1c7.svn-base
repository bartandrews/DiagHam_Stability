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


#ifndef TENSORPRODUCTSTRUCTURE_H
#define TENSORPRODUCTSTRUCTURE_H


#include "TensorProduct/AbstractTensorProductStructure.h"


#include <iostream>


using std::ostream;


class TensorProductStructure : public AbstractTensorProductStructure
{

 private:

  int* SpaceDimension;
  int* Increment;

  int* GarbageFlag;

 public:
  
  // default constructor
  //
  TensorProductStructure();
  
  // constructor from number of spaces
  //
  // nbrSpace = number of spaces
  TensorProductStructure(int nbrSpace);
  
  // copy constructor
  //
  // structure = reference on structure to copy
  TensorProductStructure(const TensorProductStructure& structure);
  
  // destructor
  //
  ~TensorProductStructure();
  
  // assignement
  //
  // structure = reference on structure to assign
  // return value = reference on current structure
  TensorProductStructure& operator = (const TensorProductStructure& structure);

  // return dimension of a given space
  //
  // space = space index
  // return value = space dimension
  int GetDimension(int space);

  // set dimension of a given space
  //
  // space = space index
  // dimension = new dimension
  void SetDimension(int space, int dimension);

  // return increment to go to the following component for a given space
  //
  // space = space index
  // return value = increment
  int GetIncrement(int space);

  // return total space dimension
  //
  // return value = dimension
  int GetTotalDimension();

  // return number of spaces
  //
  // return value = number of spaces
  int GetNbrSpace();

  // return tensor product structure of a space obtained by tensor product of two tensor spaces
  //
  // structure1 = first structure
  // structure2 = second structure
  // return value = structure corresponding to the new tensor space
  friend TensorProductStructure operator * (const TensorProductStructure& structure1, 
						const TensorProductStructure& structure2);

  // test if two tensor product structures are equivalent
  //
  // structure1 = first structure
  // structure2 = second structure
  // return value = true if structures are equivalent
  friend bool operator == (const TensorProductStructure& structure1, 
			   const TensorProductStructure& structure2);

  // test if two tensor product structures are different
  //
  // structure1 = first structure
  // structure2 = second structure
  // return value = true if structures are different
  friend bool operator != (const TensorProductStructure& structure1, 
			   const TensorProductStructure& structure2);

  // print information on current tensor product structure
  //
  // str = output stream
  // structure = structure to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, const TensorProductStructure& structure);

 private:

  // evaluate all increments
  //
  void EvaluateIncrement();

};


// return dimension of a given space
//
// space = space index
// return value = space dimension

inline int TensorProductStructure::GetDimension(int space)
{
  return this->SpaceDimension[space];
}

// return increment to go to the following component for a given space
//
// space = space index
// return value = increment

inline int TensorProductStructure::GetIncrement(int space)
{
  return this->Increment[space];
}

// return total space dimension
//
// return value = dimension

inline int TensorProductStructure::GetTotalDimension()
{
  return this->Increment[this->NbrSpace];
}

// return number of spaces
//
// return value = number of spaces

inline int TensorProductStructure::GetNbrSpace()
{
  return this->NbrSpace;
}

#endif
