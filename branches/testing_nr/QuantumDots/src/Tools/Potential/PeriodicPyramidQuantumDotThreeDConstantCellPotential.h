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


#ifndef PERIODICPYRAMIDQUANTUMDOTTHREEDCONSTANTCELLPOTENTIAL_H
#define PERIODICPYRAMIDQUANTUMDOTTHREEDCONSTANTCELLPOTENTIAL_H


#include "config.h"
#include "Tools/Potential/ThreeDConstantCellPotential.h"


class PeriodicPyramidQuantumDotThreeDConstantCellPotential : public ThreeDConstantCellPotential
{
  
 protected:

  // type of atoms, true: alloy atom, false: host atom, order: zyx
  bool*** Alloy;

  // the potential in 3D, the first index is for z, the second for y and the third for x
  double*** PotentialValue;

  // number of mono-layers under the wetting layer
  int Under;

  // number of mono-layers above the dot
  int Above;

  // width of wetting layer
  int WettingWidth;

  // base radius of the truncated pyramid
  int BaseRadius;

  // top radius of the truncated pyramid
  int TopRadius;

  // concentration of the alloy in the dot and wetting layer
  double Concentration;

  // integer of launched ConstructPotential processes, used for the random process to prevent the same draw
  int LaunchNumber;

 public:

  // constructor from geometric parameters
  //
  // numberX, numberY, numberZ = number of cells in X, Y and Z directions respectively
  // under = number of monolayers under the wetting layer
  // above = number of monolayers above the dot
  // wettingWidth= width of wetting layers in cell unit
  // baseRadius = base radius of the truncated pyramid
  // topRadius = base radius of the truncated pyramid
  PeriodicPyramidQuantumDotThreeDConstantCellPotential(int numberX, int numberY, int numberZ, int under, int above, int wettingWidth, int baseRadius, int topRadius);

  // destructor
  //
  virtual ~PeriodicPyramidQuantumDotThreeDConstantCellPotential();

  // constructor the potential from piezoelectric parameter
  //
  // noInNProbability = probability for having InN without InN neighborhood
  // withInNProbability = probability for having InN with InN neighborhood
  // piezoField = difference of piezoelectric field
  // cellSizeZ = height of a monolayer
  // offset = offset bulk-alloy material
  // scratch = true if constructed from nothing, else from existing Diagram
  // fileName = file to store parameters
  void ConstructPotential(double noInNProbability, double withInNProbability, double piezoField, double cellSizeZ, double offset, bool scratch, char* fileName);

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

// assign the potential a value at a given position 
//
// i = x coordinate of the considered cell
// j = y coordinate of the considered cell
// k = z coordinate of the considered cell
// value = value of potential

inline void PeriodicPyramidQuantumDotThreeDConstantCellPotential::SetPotential(int i, int j, int k, double value)
{
  this->PotentialValue[k][j][i] = value;
}

// get the potential at a given position
//
// i = x coordinate of the considered cell
// j = y coordinate of the considered cell
// k = z coordinate of the considered cell
// return value = the potential in the cell

inline double PeriodicPyramidQuantumDotThreeDConstantCellPotential::GetPotential(int i, int j, int k)
{
  return this->PotentialValue[k][j][i];
}

// get the concentration
//
// return = concentration value, in ration of 1.0

inline double PeriodicPyramidQuantumDotThreeDConstantCellPotential::GetConcentration()
{
  return this->Concentration;
}

#endif
