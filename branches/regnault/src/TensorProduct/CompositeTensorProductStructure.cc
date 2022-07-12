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


#include "TensorProduct/CompositeTensorProductStructure.h"
#include "GeneralTools/ListIterator.h"

#include <iostream>


using std::endl;
using std::ostream;


// default constructor
//

CompositeTensorProductStructure::CompositeTensorProductStructure() 
{  
  this->NbrSpace = 0;
  this->StructureID = AbstractTensorProductStructure::Composite;
}
  
// constructor from datas
//
// spaceStructure = list of all tensor space structure describing each subspace

CompositeTensorProductStructure::CompositeTensorProductStructure(List<TensorProductStructure> spaceStructure) 
{
  this->SpaceStructure = spaceStructure;
  this->NbrSpace = spaceStructure[0].GetNbrSpace();
  this->StructureID = AbstractTensorProductStructure::Composite;
  this->SpaceDescription = this->EvaluateSubspaceDescription();
}
  
// copy constructor
//
// structure = reference on structure to copy

CompositeTensorProductStructure::CompositeTensorProductStructure(const CompositeTensorProductStructure& structure) 
{
  this->SpaceStructure = structure.SpaceStructure;
  this->NbrSpace = structure.NbrSpace;
  this->StructureID = structure.StructureID;
  this->SpaceDescription = structure.SpaceDescription;
}
  
// destructor
//

CompositeTensorProductStructure::~CompositeTensorProductStructure() 
{
}
  
// assignement
//
// structure = reference on structure to assign
// return value = reference on current structure

CompositeTensorProductStructure& CompositeTensorProductStructure::operator = (const CompositeTensorProductStructure& structure) 
{
  this->SpaceStructure = structure.SpaceStructure;
  this->NbrSpace = structure.NbrSpace;
  this->StructureID = structure.StructureID;
  this->SpaceDescription = structure.SpaceDescription;
  return *this;
}

// return dimension of a given space contained in a given subspace
//
// space = space index
// subspace = subspace index
// return value = space dimension

int CompositeTensorProductStructure::GetDimension(int space, int subspace) 
{
  return this->SpaceStructure[subspace].GetDimension(space);
}

// return tensor product structure of given subspace
//
// subspace = subspace index
// return value = subspace tensor product structure

TensorProductStructure& CompositeTensorProductStructure::GetSubspaceTensorProductStructure(int subspace) 
{
  return this->SpaceStructure[subspace];
}

// return tensor product structure of a space obtained by tensor product of two tensor spaces
//
// structure1 = first structure
// structure2 = second structure
// return value = structure corresponding to the new tensor space

CompositeTensorProductStructure operator * (CompositeTensorProductStructure& structure1, 
					    CompositeTensorProductStructure& structure2)
{
}

// test if two tensor product structures are equivalent
//
// structure1 = first structure
// structure2 = second structure
// return value = true if structures are equivalent

bool operator == (CompositeTensorProductStructure& structure1, 
		  CompositeTensorProductStructure& structure2)
{
  if (structure1.SpaceStructure.GetNbrElement() != structure2.SpaceStructure.GetNbrElement())
    return false;
  TensorProductStructure* TmpStructure1;
  ListIterator<TensorProductStructure>  IterStructure1 (structure1.SpaceStructure);
  TensorProductStructure* TmpStructure2;
  ListIterator<TensorProductStructure>  IterStructure2 (structure2.SpaceStructure);
  while (((TmpStructure1 = IterStructure1())) && ((TmpStructure2 = IterStructure2())))
    if ((*TmpStructure1) != (*TmpStructure2))
      return false;
  return true;  
}

// test if two tensor product structures are different
//
// structure1 = first structure
// structure2 = second structure
// return value = true if structures are different

bool operator != (CompositeTensorProductStructure& structure1, 
		  CompositeTensorProductStructure& structure2)
{
  if (structure1.SpaceStructure.GetNbrElement() != structure2.SpaceStructure.GetNbrElement())
    return true;
  TensorProductStructure* TmpStructure1;
  ListIterator<TensorProductStructure>  IterStructure1 (structure1.SpaceStructure);
  TensorProductStructure* TmpStructure2;
  ListIterator<TensorProductStructure>  IterStructure2 (structure2.SpaceStructure);
  while (((TmpStructure1 = IterStructure1())) && ((TmpStructure2 = IterStructure2())))
    if ((*TmpStructure1) != (*TmpStructure2))
      return true;
  return false;  
}

// print information on current tensor product structure
//
// str = output stream
// structure = structure to print
// return value = reference on output stream

ostream& operator << (ostream& str, CompositeTensorProductStructure& structure) 
{
  str << "Total Space Dimension = " << structure.GetTotalDimension() << "  Number of subspaces = " 
      << structure.SpaceStructure.GetNbrElement() << endl;
  TensorProductStructure* TmpStructure;
  ListIterator<TensorProductStructure>  IterStructure (structure.SpaceStructure);
  while ((TmpStructure = IterStructure()))
    str << *TmpStructure; 
  return str;
}

// evaluate subspace description
//

SpaceDecomposition CompositeTensorProductStructure::EvaluateSubspaceDescription() 
{
  TensorProductStructure* TmpStructure;
  ListIterator<TensorProductStructure>  IterStructure (this->SpaceStructure);
  int* SubspacePosition = new int [this->SpaceStructure.GetNbrElement()];
  int Pos = 0;
  int i = 0;
  while ((TmpStructure = IterStructure()))
    {
      SubspacePosition[i++] = Pos;
      Pos += TmpStructure->GetTotalDimension();
    }
  return SpaceDecomposition(Pos, this->SpaceStructure.GetNbrElement(), SubspacePosition);
}
