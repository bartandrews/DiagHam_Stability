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


#include "HilbertSpace/SpaceDecomposition.h"

#include <iostream>


using std::endl;


// default constructor
//

SpaceDecomposition::SpaceDecomposition () 
{
  this->SpaceDimension = 0;
  this->NbrSubspace = 0;
  this->GarbageFlag = 0;
  this->SubspacePosition = 0;
  this->SubspaceDescription = 0;
}

// constructor from datas
//
// spaceDimension = dimension of the total space
// nbrSubspace = number of subspaces
// subspacePosition = array containing index of each subspace first component (without duplicating array)

SpaceDecomposition::SpaceDecomposition (int spaceDimension, int nbrSubspace, int* subspacePosition) 
{
  this->SpaceDimension = spaceDimension;
  this->NbrSubspace = nbrSubspace;
  this->GarbageFlag = new int;
  (*(this->GarbageFlag)) = 1;
  this->SubspacePosition = subspacePosition;
  int SubspaceDim = 0;
  int* SubspaceConverterArray;
  this->SubspaceDescription = new SubspaceSpaceConverter [this->NbrSubspace];
  for (int i = 0; i < (this->NbrSubspace - 1); i++)
    {
      SubspaceDim =  this->SubspacePosition[i + 1] - this->SubspacePosition[i];
      SubspaceConverterArray = new int [SubspaceDim];
      for (int k = 0; k < SubspaceDim; k++)
	SubspaceConverterArray[k] = subspacePosition[i] + k;
      this->SubspaceDescription[i] = SubspaceSpaceConverter (this->SpaceDimension, 
							     SubspaceDim, SubspaceConverterArray);
    }
  this->SubspacePosition[this->NbrSubspace - 1] = subspacePosition[this->NbrSubspace - 1];
  SubspaceDim = this->SpaceDimension - this->SubspacePosition[this->NbrSubspace - 1];
  SubspaceConverterArray = new int [SubspaceDim];
  for (int k = 0; k < SubspaceDim; k++)
    SubspaceConverterArray[k] = this->SubspacePosition[this->NbrSubspace - 1] + k;
  this->SubspaceDescription[this->NbrSubspace - 1] = SubspaceSpaceConverter (this->SpaceDimension, 
									     SubspaceDim, 
									     SubspaceConverterArray);
}

// copy constructor (without duplicating datas)
//
// decomposition = reference on space decomposition to copy

SpaceDecomposition::SpaceDecomposition (const SpaceDecomposition& decomposition) 
{
  if (decomposition.GarbageFlag == 0)
    {
      this->SpaceDimension = 0;
      this->NbrSubspace = 0;
      this->GarbageFlag = 0;
      this->SubspacePosition = 0;
      this->SubspaceDescription = 0;
    }
  else
    {
      this->SpaceDimension = decomposition.SpaceDimension;
      this->NbrSubspace = decomposition.NbrSubspace;
      this->GarbageFlag = decomposition.GarbageFlag;
      (*(this->GarbageFlag))++;
      this->SubspacePosition = decomposition.SubspacePosition;
      this->SubspaceDescription = decomposition.SubspaceDescription;
    }
}

// destructor
//

SpaceDecomposition::~SpaceDecomposition () 
{
  if (this->GarbageFlag != 0)
    {
      if ((*(this->GarbageFlag)) == 1)
	{
	  delete[] this->SubspacePosition;
	  delete[] this->SubspaceDescription;
	  delete this->GarbageFlag;
	}
      else
	(*(this->GarbageFlag))--;
    }
}

// assignement (without duplicating datas)
//
// converter = reference on chain to copy
// return value = reference on current chain

SpaceDecomposition& SpaceDecomposition::operator = (const SpaceDecomposition& decomposition) 
{
  if (this->GarbageFlag != 0)
    {
      if ((*(this->GarbageFlag)) == 1)
	{
	  delete[] this->SubspacePosition;
	  delete[] this->SubspaceDescription;
	  delete this->GarbageFlag;
	}
      else
	(*(this->GarbageFlag))--;
    }
  if (decomposition.GarbageFlag == 0)
    {
      this->SpaceDimension = 0;
      this->NbrSubspace = 0;
      this->GarbageFlag = 0;
      this->SubspacePosition = 0;
      this->SubspaceDescription = 0;
    }
  else
    {
      this->SpaceDimension = decomposition.SpaceDimension;
      this->NbrSubspace = decomposition.NbrSubspace;
      this->GarbageFlag = decomposition.GarbageFlag;
      (*(this->GarbageFlag))++;
      this->SubspacePosition = decomposition.SubspacePosition;
      this->SubspaceDescription = decomposition.SubspaceDescription;
    }
  return *this;
}

// find in which subspace lies a given component of total space
//
// index = index of the component 
// return value = subspace index

int SpaceDecomposition::FindComponent(int index) 
{
  int i = 1;
  while ((i < this->NbrSubspace) && (index < this->SubspacePosition[i++]));
  return --i; 
}

// print decomposition descripriton
//
// str = reference on output stream
// subspace = reference on subspace to describe
// return value = reference on output stream

ostream& operator << (ostream& str, const SpaceDecomposition& subspace) 
{
  if (subspace.GarbageFlag == 0)
    return str;
  str << "number of subsapces = " << subspace.NbrSubspace << endl;
  for (int i = 0; i < subspace.NbrSubspace; i++)
    str << "subspace " << i << ":" << endl << subspace.SubspaceDescription[i] << endl;
  return str;
}
