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


#ifndef SUBSPACESPACECONVERTER_H
#define SUBSPACESPACECONVERTER_H


#include "config.h"
#include "GeneralTools/List.h"

#include <iostream>


using std::ostream;


class RealVector;
class ComplexVector;


class SubspaceSpaceConverter
{

  friend class RealSymmetricMatrix;
  friend class RealAntisymmetricMatrix;
  friend class RealDiagonalMatrix;
  friend class ComplexDiagonalMatrix;

protected:

  int SpaceDimension;
  int SubspaceDimension;
  
  int* GarbageFlag;
  int* SubspaceSpaceConverterArray;

 public:

  // default constructor
  //
  SubspaceSpaceConverter ();

  // constructor from datas
  //
  // spaceDimension = dimension of the total space
  // subspaceDimension = subspace dimension
  // converterArray = array containing index of each subspace base vector in total space base (without duplicating array)
  SubspaceSpaceConverter (int spaceDimension, int subspaceDimension, int* converterArray);

  // constructor from datas, assuming subspace has consecutive indices in total space
  //
  // spaceDimension = dimension of the total space
  // subspaceDimension = subspace dimension
  // firstIndex = first index in total space
  SubspaceSpaceConverter (int spaceDimension, int subspaceDimension, int firstIndex);

  // constructor from a list describing distinct subspace of given space
  //
  // subspaces = list of subspaces
  SubspaceSpaceConverter (List<SubspaceSpaceConverter> subspaces);

  // copy constructor (without duplicating datas)
  //
  // converter = reference on chain to copy
  SubspaceSpaceConverter (const SubspaceSpaceConverter& chain);

  // destructor
  //
  ~SubspaceSpaceConverter ();

  // assignement (without duplicating datas)
  //
  // converter = reference on chain to copy
  // return value = reference on current chain
  SubspaceSpaceConverter& operator = (const SubspaceSpaceConverter& converter);

  // return total space dimension
  //
  // return value = space dimension
  int GetSpaceDimension();

  // return subspace dimension
  //
  // return value = subspace dimension
  int GetSubspaceDimension();

  // return space index corresponding to a subspace index
  //
  // index = subspace index
  // return value = sapce index
  int GetSpaceIndex(int index);

  // project a vector from total space to subspace
  //
  // source = vector that has to be projected
  // destination = vector where result has to be stored
  // return value = reference on resulting vector
  RealVector& SpaceToSubspace (RealVector& source, RealVector& destination);

  // project a vector from total space to subspace
  //
  // source = vector that has to be projected
  // destination = vector where result has to be stored
  // return value = reference on resulting vector
  ComplexVector& SpaceToSubspace (ComplexVector& source, ComplexVector& destination);

  // evaluate projected vector in total space
  //
  // source = projected vector
  // destination = vector where result has to be stored
  // return value = reference on resulting vector
  RealVector& SubspaceToSpace (RealVector& source, RealVector& destination);

  // evaluate projected vector in total space
  //
  // source = projected vector
  // destination = vector where result has to be stored
  // return value = reference on resulting vector
  ComplexVector& SubspaceToSpace (ComplexVector& source, ComplexVector& destination);

  // print subspace descripriton
  //
  // str = reference on output stream
  // subspace = reference on subspace to describe
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, const SubspaceSpaceConverter& subspace);

};

// return space index corresponding to a subspace index
//
// index = subspace index
// return value = sapce index

inline int SubspaceSpaceConverter::GetSpaceIndex (int index)
{
  return this->SubspaceSpaceConverterArray[index];
}

#endif


