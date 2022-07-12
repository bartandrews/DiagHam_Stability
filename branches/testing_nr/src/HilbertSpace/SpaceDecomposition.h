////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of space decomposition                       //
//                                                                            //
//                        last modification : 29/05/2001                      //
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


#ifndef SPACEDECOMPOSITION_H
#define SPACEDECOMPOSITION_H


#include "config.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"


class SpaceDecomposition
{

protected:

  int SpaceDimension;
  int NbrSubspace;
  
  int* GarbageFlag;
  int* SubspacePosition;
  SubspaceSpaceConverter* SubspaceDescription;

 public:

  // default constructor
  //
  SpaceDecomposition ();

  // constructor from datas
  //
  // spaceDimension = dimension of the total space
  // nbrSubspace = number of subspaces
  // subspacePosition = array containing index of each subspace first component (without duplicating array)
  SpaceDecomposition (int spaceDimension, int nbrSubspace, int* subspacePosition);

  // copy constructor (without duplicating datas)
  //
  // decomposition = reference on space description to copy
  SpaceDecomposition (const SpaceDecomposition& decomposition);

  // destructor
  //
  ~SpaceDecomposition ();

  // assignement (without duplicating datas)
  //
  // decomposition = reference on space decomposition to copy
  // return value = reference on current decomposition
  SpaceDecomposition& operator = (const SpaceDecomposition& decomposition);

  // return total space dimension
  //
  // return value = space dimension
  int GetSpaceDimension();

  // get number of subspaces
  //
  // return value = number of subspaces
  int GetNbrSubspace();

  // return subspace dimension of a given subspace
  //
  // index = subspace index
  // return value = subspace dimension
  int GetSubspaceDimension(int index);

  // return index of the first component of a given subspace
  //
  // index = subspace index
  // return value = subspace position
  int GetSubspacePosition(int index);
  
  // return description of a given subspace
  //
  // index = subspace index
  // return value = reference on subspace description
  SubspaceSpaceConverter& GetSubspaceDescription(int index);

  // find in which subspace lies a given component of total space
  //
  // index = index of the component 
  // return value = subspace index
  int FindComponent(int index);

  // print decomposition descripriton
  //
  // str = reference on output stream
  // subspace = reference on subspace to describe
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, const SpaceDecomposition& subspace);

};

// return total space dimension
//
// return value = space dimension

inline int SpaceDecomposition::GetSpaceDimension() 
{
  return this->SpaceDimension;
}

// get number of subspaces
//
// return value = number of subspaces

inline int SpaceDecomposition::GetNbrSubspace()
{
  return this->NbrSubspace;
}

// return subspace dimension of a given subspace
//
// index = subspace index
// return value = subspace dimension

inline int SpaceDecomposition::GetSubspaceDimension(int index) 
{
  return this->SubspaceDescription[index].GetSubspaceDimension();
}

// return index of the first component of a given subspace
//
// index = subspace index
// return value = subspace position

inline int SpaceDecomposition::GetSubspacePosition(int index) 
{
  return this->SubspacePosition[index];
}
  
// return description of a given subspace
//
// index = subspace index
// return value = reference on subspace description

inline SubspaceSpaceConverter& SpaceDecomposition::GetSubspaceDescription(int index) 
{
  return this->SubspaceDescription[index];
}

#endif


