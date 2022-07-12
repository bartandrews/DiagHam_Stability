////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of subspace-space converter                    //
//                                                                            //
//                        last modification : 18/04/2001                      //
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


#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "GeneralTools/ListIterator.h"
#include <iostream.h>


// default constructor
//

SubspaceSpaceConverter::SubspaceSpaceConverter () 
{
  this->SubspaceSpaceConverterArray = 0;
  this->SpaceDimension = 0;
  this->SubspaceDimension = 0;
  this->GarbageFlag = 0;
}

// constructor from datas
//
// spaceDimension = dimension of the total space
// subspaceDimension = subspace dimension
// converterArray = array containing index of each subspace base vector in total space base (without duplicating array)

SubspaceSpaceConverter::SubspaceSpaceConverter (int spaceDimension, int subspaceDimension, int* converterArray) 
{
  this->SpaceDimension = spaceDimension;
  this->SubspaceDimension = subspaceDimension;
  this->SubspaceSpaceConverterArray = converterArray;
  this->GarbageFlag = new int;
  (*(this->GarbageFlag)) = 1;
}

// constructor from datas, assuming subspace has consecutive indices in total space
//
// spaceDimension = dimension of the total space
// subspaceDimension = subspace dimension
// firstIndex = first index in total space

SubspaceSpaceConverter::SubspaceSpaceConverter (int spaceDimension, int subspaceDimension, int firstIndex)
{
  this->SpaceDimension = spaceDimension;
  this->SubspaceDimension = subspaceDimension;
  this->SubspaceSpaceConverterArray = new int [this->SubspaceDimension];
  for (int i = 0; i < this->SubspaceDimension; i++)
    this->SubspaceSpaceConverterArray[i] = firstIndex++;
  this->GarbageFlag = new int;
  (*(this->GarbageFlag)) = 1;
}

// constructor from a list describing distinct subspace of given space
//
// subspaces = list of subspaces

SubspaceSpaceConverter::SubspaceSpaceConverter (List<SubspaceSpaceConverter> subspaces) 
{
  this->GarbageFlag = new int;
  (*(this->GarbageFlag)) = 1;
  this->SpaceDimension = subspaces[0].SpaceDimension;
  this->SubspaceDimension = 0;
  ListIterator<SubspaceSpaceConverter> IterSubspaces (subspaces);
  SubspaceSpaceConverter* TmpSubspace;
  while ((TmpSubspace = IterSubspaces()))
    {
      this->SubspaceDimension += TmpSubspace->SubspaceDimension;
    } 
  IterSubspaces.DefineList(subspaces);
  this->SubspaceSpaceConverterArray = new int [SubspaceDimension];
  int Pos = 0;
  while ((TmpSubspace = IterSubspaces()))
    for (int i = 0; i < TmpSubspace->SubspaceDimension; i++)
      this->SubspaceSpaceConverterArray[Pos++] = TmpSubspace->SubspaceSpaceConverterArray[i];
}

// copy constructor (without duplicating datas)
//
// converter = reference on chain to copy

SubspaceSpaceConverter::SubspaceSpaceConverter (const SubspaceSpaceConverter& converter) 
{
  if (converter.GarbageFlag == 0)
    {
      this->SubspaceSpaceConverterArray = 0;
      this->SpaceDimension = 0;
      this->SubspaceDimension = 0;
      this->GarbageFlag = 0;
    }
  else
    {
      this->SubspaceSpaceConverterArray = converter.SubspaceSpaceConverterArray;
      this->SpaceDimension = converter.SpaceDimension;
      this->SubspaceDimension = converter.SubspaceDimension;
      this->GarbageFlag = converter.GarbageFlag;
      (*(this->GarbageFlag))++;
    }
}

// destructor
//

SubspaceSpaceConverter::~SubspaceSpaceConverter () 
{
  if (this->GarbageFlag != 0)
    {
      if (*(this->GarbageFlag) > 1)
	(*(this->GarbageFlag))--;
      else
	{
	  delete[] this->SubspaceSpaceConverterArray;
	  delete this->GarbageFlag;
	}
    }
}

// assignement (without duplicating datas)
//
// converter = reference on chain to copy
// return value = reference on current chain

SubspaceSpaceConverter& SubspaceSpaceConverter::operator = (const SubspaceSpaceConverter& converter) 
{
  if (this->GarbageFlag != 0)
    {
      if (*(this->GarbageFlag) > 1)
	(*(this->GarbageFlag))--;
      else
	{
	  delete[] this->SubspaceSpaceConverterArray;
	  delete this->GarbageFlag;
	  this->GarbageFlag = 0;
	}
    }
  if (converter.GarbageFlag == 0)
    {
      this->SubspaceSpaceConverterArray = 0;
      this->SpaceDimension = 0;
      this->SubspaceDimension = 0;
      this->GarbageFlag = 0;
    }
  else
    {
      this->SpaceDimension = converter.SpaceDimension;
      this->SubspaceDimension = converter.SubspaceDimension;
      this->GarbageFlag = converter.GarbageFlag;
      (*(this->GarbageFlag))++;
      this->SubspaceSpaceConverterArray = converter.SubspaceSpaceConverterArray;
    }
  return *this;
}

// return total space dimension
//
// return value = space dimension

int SubspaceSpaceConverter::GetSpaceDimension() 
{
  return this->SpaceDimension;
}

// return subspace dimension
//
// return value = subspace dimension

int SubspaceSpaceConverter::GetSubspaceDimension() 
{
  return this->SubspaceDimension;
}

// project a vector from total space to subspace
//
// source = vector that has to be projected
// destination = vector where result has to be stored
// return value = reference on resulting vector

RealVector& SubspaceSpaceConverter::SpaceToSubspace (RealVector& source, RealVector& destination) 
{
  for (int i = 0; i < this->SubspaceDimension; i++)
    destination[i] = source[this->SubspaceSpaceConverterArray[i]];
  return destination;
}

// project a vector from total space to subspace
//
// source = vector that has to be projected
// destination = vector where result has to be stored
// return value = reference on resulting vector

ComplexVector& SubspaceSpaceConverter::SpaceToSubspace (ComplexVector& source, ComplexVector& destination) 
{
  return destination;
}

// evaluate projected vector in total space
//
// source = projected vector
// destination = vector where result has to be stored
// return value = reference on resulting vector

RealVector& SubspaceSpaceConverter::SubspaceToSpace (RealVector& source, RealVector& destination) 
{
  for (int i = 0; i < this->SpaceDimension; i++)
    destination[i] = 0.0;  
  for (int i = 0; i < this->SubspaceDimension; i++)
    destination[this->SubspaceSpaceConverterArray[i]] = source[i];
  return destination;
}

// evaluate projected vector in total space
//
// source = projected vector
// destination = vector where result has to be stored
// return value = reference on resulting vector

ComplexVector& SubspaceSpaceConverter::SubspaceToSpace (ComplexVector& source, ComplexVector& destination) 
{
  return destination;
}

// print subspace descripriton
//
// str = reference on output stream
// subspace = reference on subspace to describe
// return value = reference on output stream

ostream& operator << (ostream& str, const SubspaceSpaceConverter& subspace)
{
  if (subspace.GarbageFlag != 0)
    {
      str << "space dimension = " << subspace.SpaceDimension << endl;
      str << "subspace dimension = " << subspace.SubspaceDimension << endl;
      str << "kept indices = ";
      for (int i = 0; i < subspace.SubspaceDimension; i++)
	str << subspace.SubspaceSpaceConverterArray[i] << " ";
      str << endl;
    }
  return str;
}
