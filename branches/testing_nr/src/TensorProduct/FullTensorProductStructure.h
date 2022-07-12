////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                   class for full tensor product structure                  //
//         (with index description of each state in original spaces)          //
//                                                                            //
//                        last modification : 14/06/2001                      //
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


#ifndef FULLTENSORPRODUCTSTRUCTURE_H
#define FULLTENSORPRODUCTSTRUCTURE_H


#include "TensorProduct/AbstractTensorProductStructure.h"


#include <iostream>


using std::ostream;


class FullTensorProductStructure : public AbstractTensorProductStructure
{

 private:

  int TotalSpaceDimension;
  int* SpaceDimension;
  int* Increment;
  
  int** Indices;

  int* GarbageFlag;

 public:
  
  // default constructor
  //
  FullTensorProductStructure();
  
  // constructor from number of spaces
  //
  // nbrSpace = number of spaces
  // spaceDimension = Total space dimension
  FullTensorProductStructure(int nbrSpace, int spaceDimension);
  
  // copy constructor
  //
  // structure = reference on structure to copy
  FullTensorProductStructure(const FullTensorProductStructure& structure);
  
  // destructor
  //
  ~FullTensorProductStructure();
  
  // assignement
  //
  // structure = reference on structure to assign
  // return value = reference on current structure
  FullTensorProductStructure& operator = (const FullTensorProductStructure& structure);

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

  // get index of a given state on a given space
  //
  // state = state index
  // space = space index
  // return value = reference on corresponding index  
  int& operator () (int state, int space);

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
  friend FullTensorProductStructure operator * (const FullTensorProductStructure& structure1, 
						const FullTensorProductStructure& structure2);

  // test if two tensor product structures are equivalent
  //
  // structure1 = first structure
  // structure2 = second structure
  // return value = true if structures are equivalent
  friend bool operator == (const FullTensorProductStructure& structure1, 
			   const FullTensorProductStructure& structure2);

  // test if two tensor product structures are different
  //
  // structure1 = first structure
  // structure2 = second structure
  // return value = true if structures are different
  friend bool operator != (const FullTensorProductStructure& structure1, 
			   const FullTensorProductStructure& structure2);

  // print information on current tensor product structure
  //
  // str = output stream
  // structure = structure to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, const FullTensorProductStructure& structure);

 private:

  // evaluate all increments
  //
  void EvaluateIncrement();

};


// return dimension of a given space
//
// space = space index
// return value = space dimension

inline int FullTensorProductStructure::GetDimension(int space)
{
  return this->SpaceDimension[space];
}

// return increment to go to the following component for a given space
//
// space = space index
// return value = increment

inline int FullTensorProductStructure::GetIncrement(int space)
{
  return this->Increment[space];
}

// return total space dimension
//
// return value = dimension

inline int FullTensorProductStructure::GetTotalDimension()
{
  return this->TotalSpaceDimension;
}

// return number of spaces
//
// return value = number of spaces

inline int FullTensorProductStructure::GetNbrSpace()
{
  return this->NbrSpace;
}

// get index of a given state on a given space
//
// state = state index
// space = space index
// return value = reference on corresponding index  

inline int& FullTensorProductStructure::operator () (int state, int space)
{
  return this->Indices[space][state];
}

#endif
