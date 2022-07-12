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


#include "TensorProduct/TensorProductStructure.h"

#include <iostream>


using std::endl;
using std::ostream;


// default constructor
//

TensorProductStructure::TensorProductStructure() 
{
  this->NbrSpace = 0;
  this->SpaceDimension = 0;
  this->Increment = 0;
  this->GarbageFlag = 0;
  this->StructureID = AbstractTensorProductStructure::Simple;
}
  
// constructor from number of spaces
//
// nbrSpace = number of spaces

TensorProductStructure::TensorProductStructure(int nbrSpace) 
{
  this->GarbageFlag = new int;
  (*(this->GarbageFlag)) = 1;
  this->NbrSpace = nbrSpace;
  this->SpaceDimension = new int [this->NbrSpace];
  this->Increment = new int [this->NbrSpace + 1];
  for (int i = 0; i < this->NbrSpace; i++)
    {
      this->SpaceDimension[i] = 1;
      this->Increment[i] = 1;
    }
  this->Increment[this->NbrSpace] = 1;
  this->StructureID = AbstractTensorProductStructure::Simple;
}
  
// copy constructor
//
// structure = reference on structure to copy

TensorProductStructure::TensorProductStructure(const TensorProductStructure& structure) 
{
  if (structure.NbrSpace == 0)
    {
      this->NbrSpace = 0;
      this->SpaceDimension = 0;
      this->Increment = 0;
      this->GarbageFlag = 0;
    }
  else
    {
      this->GarbageFlag = structure.GarbageFlag;
      (*(this->GarbageFlag))++;
      this->NbrSpace = structure.NbrSpace;
      this->SpaceDimension = structure.SpaceDimension;
      this->Increment = structure.Increment;
    }
  this->StructureID = AbstractTensorProductStructure::Simple;
}
  
// destructor
//

TensorProductStructure::~TensorProductStructure() 
{
  if (this->NbrSpace != 0)
    {
      if ((*(this->GarbageFlag)) == 1)
	{
	  delete[] this->Increment;
	  delete[] this->SpaceDimension;  
	  delete this->GarbageFlag;
	}
      else
	{
	  (*(this->GarbageFlag))--;
	}
    }
}
  
// assignement (without duplicating datas)
//
// structure = reference on structure to assign
// return value = reference on current structure

TensorProductStructure& TensorProductStructure::operator = (const TensorProductStructure& structure) 
{
  if (this->NbrSpace != 0)
    {
      if ((*(this->GarbageFlag)) == 1)
	{
	  delete[] this->Increment;
	  delete[] this->SpaceDimension;  
	  delete this->GarbageFlag;
	}
      else
	{
	  (*(this->GarbageFlag))--;
	}
    }
  if (structure.NbrSpace == 0)
    {
      this->NbrSpace = 0;
      this->SpaceDimension = 0;
      this->Increment = 0;
      this->GarbageFlag = 0;
    }
  else
    {
      this->GarbageFlag = structure.GarbageFlag;
      (*(this->GarbageFlag))++;
      this->NbrSpace = structure.NbrSpace;
      this->SpaceDimension = structure.SpaceDimension;
      this->Increment = structure.Increment;
    }
  this->StructureID = AbstractTensorProductStructure::Simple;
  return *this;
}

// set dimension of a given space
//
// space = space index
// dimension = new dimension

void TensorProductStructure::SetDimension(int space, int dimension) 
{
  this->SpaceDimension[space] = dimension;
  this->EvaluateIncrement();
}

// evaluate all increments
//

void TensorProductStructure::EvaluateIncrement() 
{
  this->Increment[0] = 1;
  for (int i = 0; i < NbrSpace; i++)
    this->Increment[i + 1] = this->Increment[i] * this->SpaceDimension[i];
}

// return tensor product structure of a space obtained by tensor product of two tensor spaces
//
// structure1 = first structure
// structure2 = second structure
// return value = structure corresponding to the new tensor space

TensorProductStructure operator * (const TensorProductStructure& structure1, const TensorProductStructure& structure2)
{
  TensorProductStructure Tmp (structure1.NbrSpace + structure2.NbrSpace);
  int Space = 0;
  for (int i = 0; i < structure1.NbrSpace; i++)
    Tmp.SetDimension(Space++, structure1.SpaceDimension[i]);
  for (int i = 0; i < structure2.NbrSpace; i++)
    Tmp.SetDimension(Space++, structure2.SpaceDimension[i]);
  return Tmp;
}

// test if two tensor product structures are equivalent
//
// structure1 = first structure
// structure2 = second structure
// return value = true if structures are equivalent

bool operator == (const TensorProductStructure& structure1, const TensorProductStructure& structure2)
{
  if (structure1.GarbageFlag == structure2.GarbageFlag)
    return true;
  if (structure1.NbrSpace != structure2.NbrSpace)
    return false;
   bool Flag = true;
   for (int i = 0; (i < structure1.NbrSpace) && (Flag == true); i++)
     if (structure1.SpaceDimension[i] != structure2.SpaceDimension[i])
       Flag = false;
   return Flag;
}

// test if two tensor product structures are different
//
// structure1 = first structure
// structure2 = second structure
// return value = true if structures are different

bool operator != (const TensorProductStructure& structure1, const TensorProductStructure& structure2)
{
  if (structure1.GarbageFlag == structure2.GarbageFlag)
    return false;
   if (structure1.NbrSpace != structure2.NbrSpace)
    return true;
   bool Flag = false;
   for (int i = 0; (i < structure1.NbrSpace) && (Flag == false); i++)
     if (structure1.SpaceDimension[i] != structure2.SpaceDimension[i])
       Flag = false;
   return Flag;
}

// print information on current tensor product structure
//
// str = output stream
// structure = structure to print
// return value = reference on output stream

ostream& operator << (ostream& str, const TensorProductStructure& structure)
{
  str << "number of spaces = " << structure.NbrSpace << "  total dimension = " << structure.Increment[structure.NbrSpace] << endl;
  for (int i = 0; i < structure.NbrSpace; i++)
    str << "space " << i << "  dimension = " << structure.SpaceDimension[i] << endl;
  return str;
}
