////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//                          class of abstract potential                       //
//                                                                            //
//                        last modification : 02/11/2004                      //
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


#ifndef ABSTRACTPOTENTIAL_H
#define ABSTRACTPOTENTIAL_H


#include "config.h"


class AbstractPotential
{

 protected:

 public:

  // destructor
  //
  virtual ~AbstractPotential();

  // shift the potential with a given quantity
  //
  // delta = shift value
  virtual void ShiftPotential(double delta) = 0;

  // save the diagram of atoms in a file
  //
  // fileName = name of the file to stock the diagram
  virtual void SaveDiagram(char* fileName) = 0;

  // load the diagram of atoms from a file
  //
  // fileName = name of the file in which the diagram is stocked
  virtual void LoadDiagram(char* fileName) = 0;

  // save the potential in a file
  //
  // fileName = name of the file to stock the potential
  virtual void SavePotential(char* fileName) = 0;

  // load the potential from a file
  //
  // fileName = name of the file in which the potential is stocked
  virtual void LoadPotential(char* fileName) = 0;

  // save the whole diagram presentation in a bitmap file
  //
  // fileName = name of the file to stock the diagram presentation  
  virtual void SaveBmpPicture(char* fileName) = 0;

};

#endif
