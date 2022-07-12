////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//                  class of potential of a dot embedded in a well            //
//                                                                            //
//                        last modification : 02/17/2004                      //
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


#ifndef DOTEMBEDDEDWELLTHREEDCONSTANTCELLPOTENTIAL_H
#define DOTEMBEDDEDWELLTHREEDCONSTANTCELLPOTENTIAL_H


#include "config.h"
#include "Tools/Potential/ThreeDConstantCellPotential.h"


class DotEmbeddedWellThreeDConstantCellPotential : public ThreeDConstantCellPotential
{
 protected:

    // type of atoms, 0 for well barrier, 1 for bulk and 2 for dot atoms
  short*** Alloy;

  // the potential in 3D, the first index is for z, the second for y and the third for x
  double*** PotentialValue;

  // height of well barrier below the dot
  int UnderBarrier;

  // height barrier just below the wetting layer
  int BelowWettingLayer;

  // width of the wetting layer
  int WettingWidth;

  // base radius of the dot
  int BaseRadius;

  // height of the dot
  int DotHeight;

  // top radius of the dot
  int TopRadius;

 public:

  // constructor
  //
  // numberX, numberY, numberZ = number of cells in X, Y and Z directions respectively
  // wettingWidth= width of wetting layers in cell unit
  // underBarrier = height of the layer serving as well barrier
  // belowWettingLayer = height of the layer just under the wetting layer
  // baseRadius = base radius of the truncated cone
  // dotHeight = height of the dot
  // topRadius = base radius of the truncated cone
  DotEmbeddedWellThreeDConstantCellPotential(int numberX, int numberY, int numberZ, int underBarrier, int belowWettingLayer, int wettingWidth, int baseRadius, int dotHeight, int topRadius);

  // destructor
  //
  ~DotEmbeddedWellThreeDConstantCellPotential();

  // construct potential from physical parameters 
  //
  // wellPotential = potential in the well barrier (0 reference: potential in the bulk, outside the dot)
  // dotPotential = potential in the dot (0 reference: potential in the bulk, outside the dot)
  // anisotropy = anisotropy of the dot, equal to Ry / Rx
  void ConstructPotential(double wellPotential, double dotPotential, double anisotropy = 0.0);

  // add two barriers of the same potential just below and above the dot + WL (for Vienna experimenters)
  //
  // belowBarrier = number of barrier layers just below the WL
  // aboveBarrier = number of barrier layers just above  the dot
  // potential = potential to add in these barriers (in eV unit)
  
  void AddBarrierPotential (int belowBarrier, int aboveBarrier, double potential);

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

 private:

  // determine if a cell is in an elliptical dot or wetting layer
  //
  // x = x coordinate of the cell
  // y = y coordinate of the cell
  // z = z coordinate of the cell
  // anisotropy = anisotropy of the dot, equal to Ry / Rx
  // return = true if the cell is in the dot, false otherwise  
  bool InTheEllipticalDot(int x, int y, int z, double anisotropy = 0.0);

};

// assign the potential a value at a given position 
//
// i = x coordinate of the considered cell
// j = y coordinate of the considered cell
// k = z coordinate of the considered cell
// value = value of potential

inline void DotEmbeddedWellThreeDConstantCellPotential::SetPotential(int i, int j, int k, double value)
{
  this->PotentialValue[k][j][i] = value;
}

// get the potential at a given position
//
// i = x coordinate of the considered cell
// j = y coordinate of the considered cell
// k = z coordinate of the considered cell
// return value = the potential in the cell

inline double DotEmbeddedWellThreeDConstantCellPotential::GetPotential(int i, int j, int k)
{
  return this->PotentialValue[k][j][i];
}

#endif
