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


#include "TensorProduct/AbstractTensorProductStructure.h"
#include "TensorProduct/TensorProductStructure.h"

#include <iostream>


using std::ostream;


// virtual destructor
//

AbstractTensorProductStructure::~AbstractTensorProductStructure()
{
}

// return tensor product structure of a space obtained by tensor product of two tensor spaces
//
// structure1 = first structure
// structure2 = second structure
// return value = structure corresponding to the new tensor space

AbstractTensorProductStructure* operator * (const AbstractTensorProductStructure& structure1, 
					    const AbstractTensorProductStructure& structure2)
{
  if ((structure1.StructureID == AbstractTensorProductStructure::Simple) && 
      (structure2.StructureID == AbstractTensorProductStructure::Simple))
    return new TensorProductStructure(((TensorProductStructure&) structure1) * ((TensorProductStructure&) structure2));
  return 0;
}

// test if two tensor product structures are equivalent
//
// structure1 = first structure
// structure2 = second structure
// return value = true if structures are equivalent

bool operator == (const AbstractTensorProductStructure& structure1, 
		  const AbstractTensorProductStructure& structure2)
{
  if (structure1.StructureID != structure2.StructureID)
    return false;
  if (structure1.StructureID == AbstractTensorProductStructure::Simple)
    return (((TensorProductStructure&) structure1) == ((TensorProductStructure&) structure2));
  return false; 
}

// test if two tensor product structures are different
//
// structure1 = first structure
// structure2 = second structure
// return value = true if structures are different

bool operator != (const AbstractTensorProductStructure& structure1, 
		  const AbstractTensorProductStructure& structure2)
{
  if (structure1.StructureID != structure2.StructureID)
    return true;
  if (structure1.StructureID == AbstractTensorProductStructure::Simple)
    return (((TensorProductStructure&) structure1) != ((TensorProductStructure&) structure2));
  return false; 
}

// print information on current tensor product structure
//
// str = output stream
// structure = structure to print
// return value = reference on output stream

ostream& operator << (ostream& str, const AbstractTensorProductStructure& structure)
{
  switch (structure.StructureID)
    {
      case AbstractTensorProductStructure::Simple:
	str << ((TensorProductStructure&) structure);
      break;
    }
  return str;
}

  
