////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//            class of potential in one direction with constant cells         //
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


#ifndef ONEDCONSTANTCELLPOTENTIAL_H
#define ONEDCONSTANTCELLPOTENTIAL_H


#include "config.h"
#include "Tools/Potential/AbstractPotential.h"


class OneDConstantCellPotential : public AbstractPotential
{

 protected:

  // number of cells
  int Number;

  // size of potential (in Angstrom unit)
  double Size;

  // the width of cells
  double* CellWidth;

  // potential in each cell
  double* PotentialValue;

 public:

  // constructor
  //
  // number = number of cells
  // cellWidth = array containing the width of cells
  // potentialValue = array containing the potential of cells
  OneDConstantCellPotential(int number, double* &cellWidth, double* &potentialValue);

  // destructor
  //
  virtual ~OneDConstantCellPotential();

  // get the number of cells
  //
  int GetNumberCell();

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
  // position = position of the considered cell
  // value = value of potential
  virtual void SetPotential(int position, double value);

  // get the potential at a given position
  //
  // position = position of the considered cell
  // return value = the potential in the cell
  virtual double GetPotential(int position);

  // get the size of the potential sample
  //
  // return = size value (in Angstrom unit)
  double GetSize();

  // get the width of a cell
  //
  // position = position of the cell
  // return = the width of the given cell
  double GetCellWidth(int position);

  // save the whole diagram presentation in a bitmap file
  //
  // fileName = name of the file to stock the diagram presentation
  virtual void SaveBmpPicture(char* fileName);

};

// get the number of cells
//

inline int OneDConstantCellPotential::GetNumberCell()
{
  return this->Number;
}

// assign the potential a value at a given position 
//
// position = position of the considered cell
// value = value of potential

inline void OneDConstantCellPotential::SetPotential(int position, double value)
{
  this->PotentialValue[position] = value;
}

// get the potential at a given position
//
// position = position of the considered cell
// return value = the potential in the cell

inline double OneDConstantCellPotential::GetPotential(int position)
{
  return this->PotentialValue[position];
}

// get the size of the potential sample
//
// return = size value (in Angstrom unit)

inline double OneDConstantCellPotential::GetSize()
{
  return this->Size;
}

// get the width of a cell
//
// position = position of the cell
// return = the width of the given cell

inline double OneDConstantCellPotential::GetCellWidth(int position)
{
  return this->CellWidth[position];
}

#endif
