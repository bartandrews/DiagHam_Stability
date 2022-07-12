////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//           class of potential in three directions for quantum dots          //
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


#ifndef HARDBOXPYRAMIDQUANTUMDOTTHREEDCONSTANTCELLPOTENTIAL_H
#define HARDBOXPYRAMIDQUANTUMDOTTHREEDCONSTANTCELLPOTENTIAL_H


#include "config.h"
#include "Tools/Potential/ThreeDConstantCellPotential.h"


class HardBoxPyramidQuantumDotThreeDConstantCellPotential : public ThreeDConstantCellPotential
{

 protected:

  // type of atoms, true: alloy atom, false: host atom, order: zyx
  bool*** Alloy;

  // the potential in 3D, the first index is for z, the second for y and the third for x
  double*** PotentialValue;

  // number of mono-layers under the wetting layer
  int Under;

  // size of each layer under the wetting layer
  double* UnderSize;

  // potential along growth direction, in each mono-layer under the wetting layer
  double* UnderPotentialValue;

  // number of mono-layers above the dot
  int Above;

  // sizeof each layer above the dot
  double* AboveSize;

  // potential along growth direction, in each mono-layer above the dot
  double* AbovePotentialValue;

  // width of wetting layer
  int WettingWidth;

  // base radius of the truncated pyramid
  int BaseRadius;

  // top radius of the truncated pyramid
  int TopRadius;

  // concentration of the alloy in the dot and wetting layer
  double Concentration;

 public:
  
  // constructor from geometric parameters
  //
  // numberX, numberY, numberZ = number of cells in X, Y and Z directions respectively
  // under = number of monolayers under the wetting layer
  // above = number of monolayers above the dot
  // wettingWidth= width of wetting layers in cell unit
  // baseRadius = base radius of the truncated pyramid
  // topRadius = base radius of the truncated pyramid
  HardBoxPyramidQuantumDotThreeDConstantCellPotential(int numberX, int numberY, int numberZ, int under, int above, int wettingWidth, int baseRadius, int topRadius);

  // destructor
  //
  virtual ~HardBoxPyramidQuantumDotThreeDConstantCellPotential();

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

  // contruct the potential from electric fields
  //
  // noInNProbability = probability for having InN without InN neighborhood
  // withInNProbability = probability for having InN with InN neighborhood
  // downField, wettingField, dotField, upField = the electric field in respective parts
  // offset = offset bulk-alloy material
  // cellSizeZ = height of a monolayer
  // scratch = true if constructed from nothing, else from existing Diagram
  // fileName = file to store parameters
  void ConstructPotential(double noInNProbability, double withInNProbability, double downField, double wettingField, double dotField, double upField, double offset, double cellSizeZ, bool scratch, char* fileName);

  // shift the potential with a given quantity
  //
  // delta = shift value
  virtual void ShiftPotential(double delta);

  // determine if there is any In in the first neigbor region (6 possibilities)
  //
  // m, n, h = three coordinations of the considered cell
  // return = true if there is any InN neighborhood, false otherwise
  bool Neighborhood(int m, int n, int h);

  // determine if a cell is in the dot or wetting layer
  //
  // x = x coordinate of the cell
  // y = y coordinate of the cell
  // z = z coordinate of the cell
  // return = true if the cell is in the dot, false otherwise
  bool InTheDot(int x, int y, int z); 
 
  // save the diagram of atoms in a file
  //
  // fileName = name of the file to stock the diagram
  virtual void SaveDiagram(char* fileName);

  // load the diagram of atoms from a file
  //
  // fileName = name of the file in which the diagram is stocked
  virtual void LoadDiagram(char* fileName);

  // save the potential in a file as a reduced form
  //
  // fileName = name of the file to stock the potential
  virtual void SavePotentialWithConstantField(char* fileName);

  // load the potential from a file as a reduced form
  //
  // fileName = name of the potential file
  virtual void LoadPotentialWithConstantField(char* fileName);

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

  // get the concentration
  //
  // return = concentration value, in ration of 1.0
  double GetConcentration();

  // save the whole diagram presentation in a bitmap file
  //
  // fileName = name of the file to stock the diagram presentation
  virtual void SaveBmpPicture(char* fileName);

};

// get the number of monolayers under the wetting layer
//

inline int HardBoxPyramidQuantumDotThreeDConstantCellPotential::GetUnder()
{
  return this->Under;
}

// get the height of one monolayer under the wetting layer
//
// z = the indice of the monolayer

inline double HardBoxPyramidQuantumDotThreeDConstantCellPotential::GetUnderSize(int z)
{
  return this->UnderSize[z];
}

// get the potential of one monolayer under the wetting layer
//
// z = the indice of the monolayer

inline double HardBoxPyramidQuantumDotThreeDConstantCellPotential::GetUnderPotentialValue(int z)
{
  return this->UnderPotentialValue[z];
}

// get the number of monolayers above the dot
//

inline int HardBoxPyramidQuantumDotThreeDConstantCellPotential::GetAbove()
{
  return this->Above;
}

// get the height of one monolayer above the dot
//
// z = the indice of the monolayer, which is counted from 0

inline double HardBoxPyramidQuantumDotThreeDConstantCellPotential::GetAboveSize(int z)
{
  return this->AboveSize[z];
}

// get the potential of one monolayer above the dot
//
// z = the indice of the monolayer

inline double HardBoxPyramidQuantumDotThreeDConstantCellPotential::GetAbovePotentialValue(int z)
{
  return this->AbovePotentialValue[z];
}

// assign the potential a value at a given position 
//
// i = x coordinate of the considered cell
// j = y coordinate of the considered cell 
// k = z coordinate of the considered cell 
// value = value of potential

inline void HardBoxPyramidQuantumDotThreeDConstantCellPotential::SetPotential(int i, int j, int k, double value)
{
  if (k < this->Under)
    this->UnderPotentialValue[k] = value;
  else
    if (k >= (this->NumberZ - this->Above))
      this->AbovePotentialValue[k] = value;
    else
      this->PotentialValue[k - this->Under][j][i] = value;
}

// get the potential at a given position
//
// i = x coordinate of the considered cell
// j = y coordinate of the considered cell
// k = z coordinate of the considered cell   
// return value = the potential in the cell

inline double HardBoxPyramidQuantumDotThreeDConstantCellPotential::GetPotential(int i, int j, int k)
{
  if (k < this->Under)
    return this->UnderPotentialValue[k];
  else
    if (k >= (this->NumberZ - this->Above))
      return this->AbovePotentialValue[k];
    else
      return this->PotentialValue[k - this->Under][j][i];
}

// get the concentration
//
// return = concentration value, in ration of 1.0

inline double HardBoxPyramidQuantumDotThreeDConstantCellPotential::GetConcentration()
{
  return this->Concentration;
}

#endif
