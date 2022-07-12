////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//             class of potential in two directions for binary atoms          //
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


#ifndef BINARYTWODCONSTANTCELLPOTENTIAL_H
#define BINARYTWODCONSTANTCELLPOTENTIAL_H


#include "config.h"
#include "Tools/Potential/TwoDConstantCellPotential.h"


class ThreeDConstantCellPotential;
class RealVector;

class BinaryTwoDConstantCellPotential : public TwoDConstantCellPotential
{

 protected:

  // type of atoms, true: alloy atom, false: host atom, order: yx
  bool** Alloy;

  // the potential in 2D, the first index is for y, the second for x
  double** PotentialValue;

 public:

  // constructor
  //
  // numberX, numberY = number of cells in X and Y directions respectively
  BinaryTwoDConstantCellPotential(int numberX, int numberY);

  // destructor
  //
  virtual ~BinaryTwoDConstantCellPotential();

  // construct a potential by averaging a 3D potential in the Z direction
  // 
  // potential = 3D potential with constant value in each cell
  // coefficient = the coefficient for each monolayer
  void ConstructEffectivePotential(ThreeDConstantCellPotential* potential, RealVector& coefficient);

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

  // assign the potential a value at a given position
  //
  // i = x coordinate of the considered cell
  // j = y coordinate of the considered cell
  // value = value of potential
  virtual void SetPotential(int i, int j, double& value);

  // get the potential at a given position
  //
  // i = x coordinate of the considered cell
  // j = y coordinate of the considered cell
  // return value = the potential in the cell
  virtual double GetPotential(int i, int j);

  // save the whole diagram presentation in a bitmap file
  //
  // fileName = name of the file to stock the diagram presentation
  virtual void SaveBmpPicture(char* fileName);

};

// assign the potential a value at a given position 
//
// i = x coordinate of the considered cell
// j = y coordinate of the considered cell
// value = value of potential

inline void BinaryTwoDConstantCellPotential::SetPotential(int i, int j, double& value)
{
  this->PotentialValue[j][i] = value;
}

// get the potential at a given position
//
// i = x coordinate of the considered cell
// j = y coordinate of the considered cell
// return value = the potential in the cell

inline double BinaryTwoDConstantCellPotential::GetPotential(int i, int j)
{
  return this->PotentialValue[j][i];
}

#endif
