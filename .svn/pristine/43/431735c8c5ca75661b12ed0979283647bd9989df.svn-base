////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                class for composite tensor product structure                //
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


#ifndef COMPOSITETENSORPRODUCTSTRUCTURE_H
#define COMPOSITETENSORPRODUCTSTRUCTURE_H


#include "TensorProduct/AbstractTensorProductStructure.h"
#include "TensorProduct/TensorProductStructure.h"
#include "HilbertSpace/SpaceDecomposition.h"
#include "GeneralTools/List.h"


#include <iostream>


using std::ostream;


class CompositeTensorProductStructure : public AbstractTensorProductStructure
{

 protected:

  List<TensorProductStructure> SpaceStructure;
  SpaceDecomposition SpaceDescription;
  
 public:
  
  // default constructor
  //
  CompositeTensorProductStructure();
  
  // constructor from datas
  //
  // spaceStructure = list of all tensor space structure describing each subspace
  CompositeTensorProductStructure(List<TensorProductStructure> spaceStructure);
  
  // copy constructor
  //
  // structure = reference on structure to copy
  CompositeTensorProductStructure(const CompositeTensorProductStructure& structure);
  
  // destructor
  //
  ~CompositeTensorProductStructure();
  
  // assignement
  //
  // structure = reference on structure to assign
  // return value = reference on current structure
  CompositeTensorProductStructure& operator = (const CompositeTensorProductStructure& structure);

  // get number of subspaces
  //
  // return value = number of subspaces
  int GetNbrSubspace ();

  // return dimension of a given space
  //
  // space = space index
  // return value = space dimension
  int GetDimension(int space);

  // return dimension of a given space contained in a given subspace
  //
  // space = space index
  // subspace = subspace index
  // return value = space dimension
  int GetDimension(int space, int subspace);
  
  // return dimension of a given subspace
  //
  // subspace = subspace index
  // return value = subspace dimension
  int GetSubspaceDimension(int subspace);
  
  // return tensor product structure of given subspace
  //
  // subspace = subspace index
  // return value = subspace tensor product structure
  TensorProductStructure& GetSubspaceTensorProductStructure(int subspace);
  
  // return increment to go to a given subspace
  //
  // subspace = subspace index
  // return value = increment
  int GetSubspaceIncrement(int subspace);
  
  // return total space dimension
  //
  // return value = dimension
  int GetTotalDimension();
  
  // get list of tensor product structures that composed the current tensor space
  //
  // return value = reference on list of tensor product structures
  List<TensorProductStructure>& GetTensorProductStructures();
  
  // return tensor product structure of a space obtained by tensor product of two tensor spaces
  //
  // structure1 = first structure
  // structure2 = second structure
  // return value = structure corresponding to the new tensor space
  friend CompositeTensorProductStructure operator * (CompositeTensorProductStructure& structure1, 
						     CompositeTensorProductStructure& structure2);
  
  // test if two tensor product structures are equivalent
  //
  // structure1 = first structure
  // structure2 = second structure
  // return value = true if structures are equivalent
  friend bool operator == (CompositeTensorProductStructure& structure1, 
			   CompositeTensorProductStructure& structure2);
  
  // test if two tensor product structures are different
  //
  // structure1 = first structure
  // structure2 = second structure
  // return value = true if structures are different
  friend bool operator != (CompositeTensorProductStructure& structure1, 
			   CompositeTensorProductStructure& structure2);
  
  // print information on current tensor product structure
  //
  // str = output stream
  // structure = structure to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, CompositeTensorProductStructure& structure);
  
 private:
  
  // evaluate subspace description
  //
  SpaceDecomposition EvaluateSubspaceDescription();
  
};

// return total space dimension
//
// return value = dimension

inline int CompositeTensorProductStructure::GetTotalDimension() 
{
  return this->SpaceDescription.GetSpaceDimension();
}

// get number of subspaces
//
// return value = nimber of subspaces

inline int CompositeTensorProductStructure::GetNbrSubspace () 
{
  return this->SpaceStructure.GetNbrElement();
}

// return dimension of a given subspace
//
// subspace = subspace index
// return value = subspace dimension

inline int CompositeTensorProductStructure::GetSubspaceDimension(int subspace) 
{
  return this->SpaceDescription.GetSubspaceDimension(subspace);
}

// return increment to go to a given subspace
//
// subspace = subspace index
// return value = increment

inline int CompositeTensorProductStructure::GetSubspaceIncrement(int subspace) 
{
  return this->SpaceDescription.GetSubspacePosition(subspace);
}

// return dimension of a given space
//
// space = space index
// return value = space dimension

inline int CompositeTensorProductStructure::GetDimension(int space)
{
  return 0;
}

// get list of tensor product structures that composed the current tensor space
//
// return value = reference on list of tensor product structures

inline List<TensorProductStructure>& CompositeTensorProductStructure::GetTensorProductStructures()
{
  return this->SpaceStructure;
}
  
#endif
