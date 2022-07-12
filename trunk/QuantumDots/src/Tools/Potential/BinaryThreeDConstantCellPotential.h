////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//           class of potential in three directions for binary atoms          //
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


#ifndef BINARYTHREEDCONSTANTCELLPOTENTIAL_H
#define BINARYTHREEDCONSTANTCELLPOTENTIAL_H


#include "config.h"
#include "Tools/Potential/ThreeDConstantCellPotential.h"

class PicRGB;
class AbstractBitmapPicture;

class BinaryThreeDConstantCellPotential : public ThreeDConstantCellPotential
{

 protected:

  // type of atoms, true: alloy atom, false: host atom, order: zyx
  bool*** Alloy;

  // the potential in 3D, the first index is for z, the second for y and the third for x
  double*** PotentialValue;

 public:

  // default constructor
  //
  BinaryThreeDConstantCellPotential ();

  // constructor
  //
  // numberX, numberY, numberZ = number of cells in X, Y and Z directions respectively
  BinaryThreeDConstantCellPotential(int numberX, int numberY, int numberZ);

  // destructor
  //
  virtual ~BinaryThreeDConstantCellPotential();

  // save the diagram of atoms in a file
  //
  // fileName = name of the file to stock the diagram
  virtual void SaveDiagram(char* fileName);

  // load the diagram of atoms from a file
  //
  // fileName = name of the file in which the diagram is stocked
  virtual void LoadDiagram(char* fileName);

  // save the potential description of atoms in a file (binary mode)
  //
  // fileName = name of the file to stock the potential description
  virtual void SaveBinaryPotential(char* fileName);

  // load the potential description of atoms from a file (binary mode)
  //
  // fileName = name of the file in which the potential description is stocked
  virtual void LoadBinaryPotential(char* fileName);

  // assign the potential a value at a given position
  //
  // i = x coordinate of the considered cell
  // j = y coordinate of the considered cell
  // k = z coordinate of the considered cell
  // value = value of potential
  virtual void SetPotential(int i, int j, int k, double value);

  // get the potential at a given position
  //
  // i = x coordinate of the considered cell
  // j = y coordinate of the considered cell
  // k = z coordinate of the considered cell
  // return value = the potential in the cell
  virtual double GetPotential(int i, int j, int k);

  // save the whole diagram presentation in a bitmap file
  //
  // fileName = name of the file to stock the diagram presentation
  virtual void SaveBmpPicture(char* fileName);

  // save the whole diagram presentation in a bitmap file
  //
  // sizeX = the size in X direction of each cell in pixel
  // sizeY = the size in Y direction of each cell in pixel
  // fileName = name of the file to stock the diagram presentation
  void SaveBmpPicture(char* fileName, int sizeX, int sizeY);

  // save the diagram to a bitmap file
  //
  // under = the number of monolayers omitted under the considered structure
  // above = the number of monolayers omitted above the considered structure
  // startX = point to start in diagram in X direction
  // endX = point to end in diagram in X direction
  // startY = point to start in diagram in Y direction
  // endY = point to end in diagram in Y direction
  // choice = choice of orientation: 1 for XY, 2 for XZ, 3 for YZ
  // sizeX = the size in X direction of each cell in pixel
  // sizeY = the size in Y direction of each cell in pixel
  // InN = color (in RGB definition) of InN cell
  // GaN = color (in RGB definition) of GaN cell
  // background = color (in RGB definition) of background
  // NbrX = number of cell displayed in X direction
  // fileName = name of the file to store the picture
  // return = true if succeeded, false otherwise
  bool SaveBmpPicture(int under, int above, int startX, int endX, int startY, int endY, int choice, int sizeX, int sizeY, PicRGB& InN, PicRGB& GaN, PicRGB& background, int NbrX, char* fileName);
 
};

// fill a cell with a given color
//
// startX: X coordination to start
// sizeX:  size in X direction
// startY: Y coordination to start
// sizeY:  size in Y direction
// Col: color in RGB definition
// picture: pointer to picture which is filled
void CellFill(int startX, int sizeX, int startY, int sizeY, PicRGB& Col, AbstractBitmapPicture* picture);

// assign the potential a value at a given position 
//
// i = x coordinate of the considered cell
// j = y coordinate of the considered cell 
// k = z coordinate of the considered cell 
// value = value of potential

inline void BinaryThreeDConstantCellPotential::SetPotential(int i, int j, int k, double value)
{
  this->PotentialValue[k][j][i] = value;
}

// get the potential at a given position
//
// i = x coordinate of the considered cell
// j = y coordinate of the considered cell
// k = z coordinate of the considered cell   
// return value = the potential in the cell

inline double BinaryThreeDConstantCellPotential::GetPotential(int i, int j, int k)
{
  return this->PotentialValue[k][j][i];
}

#endif
