////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of abstract interaction                      //
//                                                                            //
//                        last modification : 05/04/2001                      //
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


#include "Interaction/AbstractInteraction.h"


// destructor
//

AbstractInteraction::~AbstractInteraction() 
{
}

// set index of space where left interaction acts
//
// index = space index

void AbstractInteraction::SetLeftSpaceIndex (int index)
{
  this->LeftSpaceIndex = index;
}

// set index of space where right interaction acts
//
// index = space index

void AbstractInteraction::SetRightSpaceIndex (int index)
{
  this->RightSpaceIndex = index;
}

// set structure of tensor space where interactions act
//
// structure = tensor space structure

void AbstractInteraction::SetTensorProductStructure (AbstractTensorProductStructure* structure)
{
  this->Structure = structure;
}
