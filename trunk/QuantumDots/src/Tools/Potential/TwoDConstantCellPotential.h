////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//           class of potential in two directions with constant cells         //
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


#ifndef TWODCONSTANTCELLPOTENTIAL_H
#define TWODCONSTANTCELLPOTENTIAL_H


#include "config.h"
#include "Tools/Potential/AbstractPotential.h"


class TwoDConstantCellPotential : public AbstractPotential
{

 protected:

  // number of cells in X direction
  int NumberX;

  // number of cells in Y direction
  int NumberY;

 public:

  // destructor
  //
  virtual ~TwoDConstantCellPotential();

  // get the number of cells in X direction
  //
  int GetNumberCellX();

  // get the number of cells in Y direction
  //
  int GetNumberCellY();

  // shift the potential with a given quantity
  //
  // delta = shift value
  virtual void ShiftPotential(double delta);

  // save the diagram of atoms in a file
  //
  // fileName = name of the file to stock the diagram
  virtual void SaveDiagram(char* fileName);

  // load the diagram of atoms from a file
  //
  // fileName = name of the file in which the diagram is stocked
  virtual void LoadDiagram(char* fileName);

  // save the potential in a file
  //
  // fileName = name of the file to stock the potential
  virtual void SavePotential(char* fileName);

  // load the potential from a file
  //
  // fileName = name of the file in which the potential is stocked
  virtual void LoadPotential(char* fileName);

  // assign the potential a value at a given position 
  //
  // i = x coordinate of the considered cell
  // j = y coordinate of the considered cell 
  // value = value of potential
  virtual void SetPotential(int i, int j, double& value) = 0;

  // get the potential at a given position
  //
  // i = x coordinate of the considered cell
  // j = y coordinate of the considered cell
  // return value = the potential in the cell
  virtual double GetPotential(int i, int j) = 0;

  // save the whole diagram presentation in a bitmap file
  //
  // fileName = name of the file to stock the diagram presentation
  virtual void SaveBmpPicture(char* fileName) = 0;

};

// get the number of cells in X direction
//

inline int TwoDConstantCellPotential::GetNumberCellX()
{
  return this->NumberX;
}

// get the number of cells in Y direction
//

inline int TwoDConstantCellPotential::GetNumberCellY()
{
  return this->NumberY;
}

#endif
