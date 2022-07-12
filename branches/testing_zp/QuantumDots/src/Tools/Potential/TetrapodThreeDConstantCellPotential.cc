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


#include "Tools/Potential/TetrapodThreeDConstantCellPotential.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using std::ifstream;
using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


// numberX, numberY, numberZ = number of cells in X, Y and Z directions respectively
// belowTetrapod = height of the barrier just below the tetrapod (in cell unit)
// dotRadius = radius of the spherical dot in Angstrom unit
// armLength = length of the four arms in Angstrom unit
// armRadius = radius of the arm in Angstrom unit
// cellSize = cell size, supposed as cubic cell, in Angstrom unit

TetrapodThreeDConstantCellPotential::TetrapodThreeDConstantCellPotential (int numberX, int numberY, int numberZ, int belowTetrapod, double dotRadius, double armLength, double armRadius, double cellSize)
{
  this->NumberX = numberX;
  this->NumberY = numberY;
  this->NumberZ = numberZ;
  this->BelowTetrapod = belowTetrapod;
  this->DotRadius = dotRadius;
  this->ArmLength = armLength;
  this->ArmRadius = armRadius;
  this->CellSize = cellSize;

  this->Alloy = new bool** [this->NumberZ];
  this->PotentialValue = new double** [this->NumberZ];
  for (int k = 0; k < this->NumberZ; ++k)
    {
      this->Alloy[k] = new bool* [this->NumberY];
      this->PotentialValue[k] = new double* [this->NumberY];
      for (int j = 0; j < this->NumberY; ++j)
	{
	  this->Alloy[k][j] = new bool [this->NumberX];
	  this->PotentialValue[k][j] = new double [this->NumberX];
	}
    }

  double R = (this->DotRadius + this->ArmLength);
  double a = R * 4.0 / sqrt(6.0);
  double x = a / sqrt(3.0);
  double d = a * sqrt(3.0) / 6.0;
  double h = a * sqrt(6.0) / 3.0;
  double r = a * sqrt(6.0) / 12.0;
  double delta = ((double) this->BelowTetrapod - 0.1) * this->CellSize;

  this->Center = RealVector (3);
  this->Center[0] = ((double) this->NumberX) * this->CellSize / 2.0;
  this->Center[1] = ((double) this->NumberY) * this->CellSize / 2.0;
  this->Center[2] = delta + r;
  
  this->LowerVertice1 = RealVector (3);
  this->LowerVertice1[0] = x + this->Center[0];
  this->LowerVertice1[1] = this->Center[1];
  this->LowerVertice1[2] = delta;

  this->LowerVertice2 = RealVector (3);
  this->LowerVertice2[0] = -d + this->Center[0];
  this->LowerVertice2[1] = a / 2.0 + this->Center[1];
  this->LowerVertice2[2] = delta;

  this->LowerVertice3 = RealVector (3);
  this->LowerVertice3[0] = -d + this->Center[0];
  this->LowerVertice3[1] = -a / 2.0 + this->Center[1];
  this->LowerVertice3[2] = delta;

  this->UpperVertice = RealVector (3);
  this->UpperVertice[0] = this->Center[0];
  this->UpperVertice[1] = this->Center[1];
  this->UpperVertice[2] = delta + h;
}

// destructor
//
TetrapodThreeDConstantCellPotential::~TetrapodThreeDConstantCellPotential()
{
  delete[] this->Alloy;
  delete[] this->PotentialValue;
}

// construct potential from physical parameters 
//
// dotPotential = potential in the dot (0 reference: potential in the bulk, outside the dot)

void TetrapodThreeDConstantCellPotential::ConstructPotential(double dotPotential)
{
  for (int k = 0; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	if (this->InTheDot(i, j, k) || this->InTheArms(i, j, k))
	  {
	    this->Alloy[k][j][i] = true;
	    this->PotentialValue[k][j][i] = dotPotential;  
	  }
	else
	  {
	    this->Alloy[k][j][i] = false;
	    this->PotentialValue[k][j][i] = 0.0;
	  } 
}

// shift the potential with a given quantity
//
// delta = shift value

void TetrapodThreeDConstantCellPotential::ShiftPotential(double delta)
{
  for (int k = 0; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	this->PotentialValue[k][j][i] += delta;
}
  
// determine if a cell is in the dot
//
// x = x coordinate of the cell
// y = y coordinate of the cell
// z = z coordinate of the cell
// return = true if the cell is in the dot, false otherwise

bool TetrapodThreeDConstantCellPotential::InTheDot(int x, int y, int z)
{
  double X = ((double) x) * this->CellSize - this->Center[0];
  double Y = ((double) y) * this->CellSize - this->Center[1];
  double Z = ((double) z) * this->CellSize - this->Center[2];

  double r = sqrt (X * X + Y * Y + Z * Z);
  if (r < this->DotRadius)
    return true;
  else
    return false;
}

// determine if a cell is in the arms
//
// x = x coordinate of the cell
// y = y coordinate of the cell
// z = z coordinate of the cell
// return = true if the cell is in the dot, false otherwise

bool TetrapodThreeDConstantCellPotential::InTheArms(int x, int y, int z)
{
  double X = ((double) x) * this->CellSize;
  double Y = ((double) y) * this->CellSize;
  double Z = ((double) z) * this->CellSize;

  if (this->InOneArm(X, Y, Z, this->Center, this->UpperVertice, this->ArmRadius))
    return true;
  if (this->InOneArm(X, Y, Z, this->Center, this->LowerVertice1, this->ArmRadius))
    return true;
  if (this->InOneArm(X, Y, Z, this->Center, this->LowerVertice2, this->ArmRadius))
    return true;
  if (this->InOneArm(X, Y, Z, this->Center, this->LowerVertice3, this->ArmRadius))
    return true;
  
  return false;  
}

// determine if a cell is in one arm
//
// x = x coordinate of the cell
// y = y coordinate of the cell
// z = z coordinate of the cell
// a, b = vector of coordinates of the two ends
// r = radius of the arm
// return = true if the cell is in the dot, false otherwise

bool TetrapodThreeDConstantCellPotential::InOneArm(double x, double y, double z, RealVector& a, RealVector& b, double r)
{
  if (r <= 0)
    return false;

  // test if the point is in good direction
  double tmp = (x - a[0]) * (b[0] - a[0]) + (y - a[1]) * (b[1] - a[1]) + (z - a[2]) * (b[2] - a[2]);
  double l = (b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]) + (b[2] - a[2]) * (b[2] - a[2]);
  double tmp1 = tmp / l;
  if ((tmp1 < 0.0) || (tmp1 > 1.0))
    return false;

  // test if the point is not too far from the line
  double tmp2 = (x - a[0]) * (x - a[0]) + (y - a[1]) * (y - a[1]) + (z - a[2]) * (z - a[2]);
  double tmp3 = tmp2 - tmp * tmp / l;
  if (tmp3 < (r * r))
    return true;

  return false;
}
/*
// save the diagram of atoms in a file
//
// fileName = name of the file to stock the diagram

void TetrapodThreeDConstantCellPotential::SaveDiagram(char* fileName)
{
  ofstream file(fileName);
  if (!file.is_open())
    cout << "Error when open the file: " << fileName << " to write diagram." << endl;
  for (int k = 0; k < this->NumberZ; ++k)
    {
      for (int j = 0; j < this->NumberY; ++j)
	{
	  for (int i = 0; i < this->NumberX; ++i)
	    file << this->Alloy[k][j][i] << " ";
	  file << '\n';
	}
      file << '\n';
    }
  file.close();
}

// load the diagram of atoms from a file
//
// fileName = name of the file in which the diagram is stocked

void TetrapodThreeDConstantCellPotential::LoadDiagram(char* fileName)
{
  ifstream file(fileName);
  if (! file.is_open())
    {
      cout << "Error when open the diagram file: " << fileName << " . Exit now" << endl;
      exit(1);
    }
  for (int k = 0; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	file >> this->Alloy[k][j][i];
  file.close();
}

// save the whole diagram presentation in a bitmap file
//
// fileName = name of the file to stock the diagram presentation

void TetrapodThreeDConstantCellPotential::SaveBmpPicture(char* fileName)
{
}
*/
