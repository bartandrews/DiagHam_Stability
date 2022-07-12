////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//                     class of potential of a tetrapod dot                   //
//                                                                            //
//                        last modification : 11/17/2004                      //
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


#ifndef TETRAPODTHREEDCONSTANTCELLPOTENTIAL_H
#define TETRAPODTHREEDCONSTANTCELLPOTENTIAL_H


#include "config.h"
#include "Tools/Potential/BinaryThreeDConstantCellPotential.h"
#include "Vector/RealVector.h"


class TetrapodThreeDConstantCellPotential : public BinaryThreeDConstantCellPotential
{

 protected:
  /*
    // type of atoms, false for bulk and true for tetrapod atoms
  bool*** Alloy;

  // the potential in 3D, the first index is for z, the second for y and the third for x
  double*** PotentialValue;
  */
  // height of the barrier just below the tetrapod (in cell unit)
  int BelowTetrapod;

  // radius of the spherical dot in Angstrom unit
  double DotRadius;

  // length of the four arms in Angstrom unit
  double ArmLength;

  // radius of the arm in Angstrom unit
  double ArmRadius;

  // cell size, supposed as cubic cell, in Angstrom unit
  double CellSize;

 private:

  // coordinates of the center of the dot
  RealVector Center;
  
  // coordinates of the other verticex
  RealVector UpperVertice;
  RealVector LowerVertice1;
  RealVector LowerVertice2;
  RealVector LowerVertice3;

 public:

  // constructor
  //
  // numberX, numberY, numberZ = number of cells in X, Y and Z directions respectively
  // belowTetrapod = height of the barrier just below the tetrapod (in cell unit)
  // dotRadius = radius of the spherical dot in Angstrom unit
  // armLength = length of the four arms in Angstrom unit
  // armRadius = radius of the arm in Angstrom unit
  // cellSize = cell size, supposed as cubic cell, in Angstrom unit
  TetrapodThreeDConstantCellPotential (int numberX, int numberY, int numberZ, int belowTetrapod, double dotRadius, double armLength, double armRadius, double cellSize);

  // destructor
  //
  ~TetrapodThreeDConstantCellPotential();

  // construct potential from physical parameters 
  //
  // dotPotential = potential in the dot (0 reference: potential in the bulk, outside the dot)
  void ConstructPotential(double dotPotential);

  // shift the potential with a given quantity
  //
  // delta = shift value
  virtual void ShiftPotential(double delta);
 
  // determine if a cell is in the dot
  //
  // x = x coordinate of the cell
  // y = y coordinate of the cell
  // z = z coordinate of the cell
  // return = true if the cell is in the dot, false otherwise
  bool InTheDot(int x, int y, int z);

  // determine if a cell is in the arms
  //
  // x = x coordinate of the cell
  // y = y coordinate of the cell
  // z = z coordinate of the cell
  // return = true if the cell is in the dot, false otherwise
  bool InTheArms(int x, int y, int z);

  // determine if a cell is in one arm
  //
  // x = x coordinate of the cell
  // y = y coordinate of the cell
  // z = z coordinate of the cell
  // a, b = vector of coordinates of the two ends
  // r = radius of the arm
  // return = true if the cell is in the dot, false otherwise
  bool InOneArm(double x, double y, double z, RealVector& a, RealVector& b, double r);
  /*
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
  */
};
/*
// assign the potential a value at a given position 
//
// i = x coordinate of the considered cell
// j = y coordinate of the considered cell
// k = z coordinate of the considered cell
// value = value of potential

inline void TetrapodThreeDConstantCellPotential::SetPotential(int i, int j, int k, double value)
{
  this->PotentialValue[k][j][i] = value;
}

// get the potential at a given position
//
// i = x coordinate of the considered cell
// j = y coordinate of the considered cell
// k = z coordinate of the considered cell
// return value = the potential in the cell

inline double TetrapodThreeDConstantCellPotential::GetPotential(int i, int j, int k)
{
  return this->PotentialValue[k][j][i];
}
*/
#endif
