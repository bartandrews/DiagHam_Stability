////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//                     class of potential of an elliptical dot                //
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


#ifndef ELLIPTICALDOTTHREEDCONSTANTCELLPOTENTIAL_H
#define ELLIPTICALDOTTHREEDCONSTANTCELLPOTENTIAL_H


#include "config.h"
#include "Tools/Potential/ThreeDConstantCellPotential.h"


class EllipticalDotThreeDConstantCellPotential : public ThreeDConstantCellPotential
{
 protected:

    // type of atoms, 0 for well barrier, 1 for bulk and 2 for dot atoms
  short*** Alloy;

  // the potential in 3D, the first index is for z, the second for y and the third for x
  double*** PotentialValue;

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

  // anisotropy of the dot, epsilon defined by Rx / Ry = 1 + epsilon
  double Anisotropy;
  
  // cell size in Z direction (in Angstrom unit)
  double CellSizeZ;

 public:

  // constructor
  //
  // numberX, numberY, numberZ = number of cells in X, Y and Z directions respectively
  // wettingWidth= width of wetting layers in cell unit
  // belowWettingLayer = height of the layer just under the wetting layer
  // baseRadius = base radius of the truncated cone
  // dotHeight = height of the dot
  // topRadius = base radius of the truncated cone
  // anisotropy = anisotropy factor
  EllipticalDotThreeDConstantCellPotential(int numberX, int numberY, int numberZ, int belowWettingLayer, int wettingWidth, int baseRadius, int dotHeight, int topRadius, double anisotropy);

  // destructor
  //
  ~EllipticalDotThreeDConstantCellPotential();

  // construct potential from physical parameters 
  //
  // dotPotential = potential in the dot (0 reference: potential in the bulk, outside the dot)
  void ConstructPotential(double dotPotential);

  // construct potential from physical parameters 
  //
  // dotPotential = potential in the dot (0 reference: potential in the bulk, outside the dot)
  // maxPotential = potential max caused by strain
  void ConstructPotential(double dotPotential, double maxPotential);

  // shift the potential with a given quantity
  //
  // delta = shift value
  virtual void ShiftPotential(double delta);
 
  // determine if a cell is in the dot or wetting layer
  //
  // x = x coordinate of the cell
  // y = y coordinate of the cell
  // z = z coordinate of the cell
  // return = true if the cell is in the dot, false otherwise
  bool InTheDot(int x, int y, int z);

  // set the size of each cell
  //
  void SetCellSizeZ (double lz);

  // get the number of monolayers under the wetting layer
  //
  int GetUnder();

  // get the height of one monolayer under the wetting layer
  //
  // z = the indice of the monolayer
  double GetUnderSize(int z);

  // get the potential of one monolayer under the wetting layer
  //
  // z = the indice of the monolayer
  double GetUnderPotentialValue(int z);

  // get the number of monolayers above the dot
  //
  int GetAbove();

  // get the height of one monolayer above the dot
  //
  // z = the indice of the monolayer, which is counted from 0
  double GetAboveSize(int z);

  // get the potential of one monolayer above the dot
  //
  // z = the indice of the monolayer, which is counted from 0
  double GetAbovePotentialValue(int z);

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

};

// set the size of each cell
//

inline void EllipticalDotThreeDConstantCellPotential::SetCellSizeZ (double lz)
{
  this->CellSizeZ = lz;
}

// get the number of monolayers under the wetting layer
//

inline int EllipticalDotThreeDConstantCellPotential::GetUnder()
{
  return this->BelowWettingLayer;
}

// get the height of one monolayer under the wetting layer
//
// z = the indice of the monolayer

inline double EllipticalDotThreeDConstantCellPotential::GetUnderSize(int z)
{
  return this->CellSizeZ;
}

// get the potential of one monolayer under the wetting layer
//
// z = the indice of the monolayer

inline double EllipticalDotThreeDConstantCellPotential::GetUnderPotentialValue(int z)
{
  return this->PotentialValue[0][0][z];
}

// get the number of monolayers above the dot
//

inline int EllipticalDotThreeDConstantCellPotential::GetAbove()
{
  return (this->NumberZ - this->BelowWettingLayer - this->WettingWidth - this->DotHeight);
}

// get the height of one monolayer above the dot
//
// z = the indice of the monolayer, which is counted from 0

inline double EllipticalDotThreeDConstantCellPotential::GetAboveSize(int z)
{
  return this->CellSizeZ;
}

// get the potential of one monolayer above the dot
//
// z = the indice of the monolayer

inline double EllipticalDotThreeDConstantCellPotential::GetAbovePotentialValue(int z)
{
  return this->PotentialValue[0][0][z + this->BelowWettingLayer + this->WettingWidth + this->DotHeight];
}

// assign the potential a value at a given position 
//
// i = x coordinate of the considered cell
// j = y coordinate of the considered cell
// k = z coordinate of the considered cell
// value = value of potential

inline void EllipticalDotThreeDConstantCellPotential::SetPotential(int i, int j, int k, double value)
{
  this->PotentialValue[k][j][i] = value;
}

// get the potential at a given position
//
// i = x coordinate of the considered cell
// j = y coordinate of the considered cell
// k = z coordinate of the considered cell
// return value = the potential in the cell

inline double EllipticalDotThreeDConstantCellPotential::GetPotential(int i, int j, int k)
{
  return this->PotentialValue[k][j][i];
}

#endif
