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


#include "Tools/Potential/EllipticalDotThreeDConstantCellPotential.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using std::ifstream;
using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


// constructor
//
// numberX, numberY, numberZ = number of cells in X, Y and Z directions respectively
// wettingWidth= width of wetting layers in cell unit
// belowWettingLayer = height of the layer just under the wetting layer
// baseRadius = base radius of the truncated cone
// dotHeight = height of the dot
// topRadius = base radius of the truncated cone
// anisotropy = anisotropy factor

EllipticalDotThreeDConstantCellPotential::EllipticalDotThreeDConstantCellPotential(int numberX, int numberY, int numberZ, int belowWettingLayer, int wettingWidth, int baseRadius, int dotHeight, int topRadius, double anisotropy)
{
  this->NumberX = numberX;
  this->NumberY = numberY;
  this->NumberZ = numberZ;
  this->BelowWettingLayer = belowWettingLayer;
  this->WettingWidth = wettingWidth;
  this->DotHeight = dotHeight;
  this->BaseRadius = baseRadius;
  this->TopRadius = topRadius;
  this->Anisotropy = anisotropy;
  this->Alloy = new short** [this->NumberZ];
  this->PotentialValue = new double** [this->NumberZ];
  for (int k = 0; k < this->NumberZ; ++k)
    {
      this->Alloy[k] = new short* [this->NumberY];
      this->PotentialValue[k] = new double* [this->NumberY];
      for (int j = 0; j < this->NumberY; ++j)
	{
	  this->Alloy[k][j] = new short [this->NumberX];
	  this->PotentialValue[k][j] = new double [this->NumberX];
	}
    }
}

// destructor
//
EllipticalDotThreeDConstantCellPotential::~EllipticalDotThreeDConstantCellPotential()
{
  delete[] this->Alloy;
  delete[] this->PotentialValue;
}

// construct potential from physical parameters 
//
// dotPotential = potential in the dot (0 reference: potential in the bulk, outside the dot)

void EllipticalDotThreeDConstantCellPotential::ConstructPotential(double dotPotential)
{
  // below the wetting layer
  for (int k = 0; k < this->BelowWettingLayer; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	{
	  this->Alloy[k][j][i] = 1;
	  this->PotentialValue[k][j][i] = 0.0;
	}
  // wetting layer
  for (int k = this->BelowWettingLayer; k < (this->BelowWettingLayer + this->WettingWidth); ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	{
	  this->Alloy[k][j][i] = 2;
	  this->PotentialValue[k][j][i] = dotPotential;
	}
  // dot
  for (int k = this->BelowWettingLayer + this->WettingWidth; k < (this->BelowWettingLayer + this->WettingWidth + this->DotHeight); ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	if (this->InTheDot(i, j, k))
	  {
	    this->Alloy[k][j][i] = 2;
	    this->PotentialValue[k][j][i] = dotPotential;  
	  }
	else
	  {
	    this->Alloy[k][j][i] = 1;
	    this->PotentialValue[k][j][i] = 0.0;
	  }
  // above the dot  
  for (int k = this->BelowWettingLayer + this->WettingWidth + this->DotHeight; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	{
	  this->Alloy[k][j][i] = 1;
	  this->PotentialValue[k][j][i] = 0.0;
	}
}

// construct potential from physical parameters 
//
// dotPotential = potential in the dot (0 reference: potential in the bulk, outside the dot)
// maxPotential = potential max caused by strain

void EllipticalDotThreeDConstantCellPotential::ConstructPotential(double dotPotential, double maxPotential)
{
  // below the wetting layer
  for (int k = 0; k < this->BelowWettingLayer; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	{
	  this->Alloy[k][j][i] = 1;
	  this->PotentialValue[k][j][i] = 0.0;
	}
  // wetting layer
  for (int k = this->BelowWettingLayer; k < (this->BelowWettingLayer + this->WettingWidth); ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	{
	  this->Alloy[k][j][i] = 2;
	  this->PotentialValue[k][j][i] = dotPotential;
	}
  // dot
  for (int k = this->BelowWettingLayer + this->WettingWidth; k < (this->BelowWettingLayer + this->WettingWidth + this->DotHeight); ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	if (this->InTheDot(i, j, k))
	  {
	    double RkX = double(this->BaseRadius) - double(k - (this->BelowWettingLayer + this->WettingWidth)) * double(this->BaseRadius - this->TopRadius) / double(this->DotHeight);
	    double RkY = RkX * (1.0 + this->Anisotropy);
	    double Tmp = double((i - this->NumberX / 2) * (i - this->NumberX / 2)) / (RkX * RkX) - double((j - this->NumberY / 2) * (j - this->NumberY / 2)) / (RkY * RkY);
	    this->Alloy[k][j][i] = 2;
	    this->PotentialValue[k][j][i] = dotPotential + maxPotential * Tmp;  
	  }
	else
	  {
	    this->Alloy[k][j][i] = 1;
	    this->PotentialValue[k][j][i] = 0.0;
	  }
  // above the dot  
  for (int k = this->BelowWettingLayer + this->WettingWidth + this->DotHeight; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	{
	  this->Alloy[k][j][i] = 1;
	  this->PotentialValue[k][j][i] = 0.0;
	}
}

// shift the potential with a given quantity
//
// delta = shift value

void EllipticalDotThreeDConstantCellPotential::ShiftPotential(double delta)
{
  for (int k = 0; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	this->PotentialValue[k][j][i] += delta;
}
  
// determine if a cell is in the dot or wetting layer
//
// x = x coordinate of the cell
// y = y coordinate of the cell
// z = z coordinate of the cell
// return = true if the cell is in the dot, false otherwise

bool EllipticalDotThreeDConstantCellPotential::InTheDot(int x, int y, int z)
{
  if (z < this->BelowWettingLayer)
    return false;
  if ((z >= (this->BelowWettingLayer)) && (z < (this->BelowWettingLayer + this->WettingWidth)))
    return true;
  if (z >= (this->BelowWettingLayer + this->WettingWidth + this->DotHeight))
    return false;
  double RkX = double(this->BaseRadius) - double(z - (this->BelowWettingLayer + this->WettingWidth)) * double(this->BaseRadius - this->TopRadius) / double(this->DotHeight);
  double RkY = RkX * (1.0 + this->Anisotropy);
  double TmpX = double((x - this->NumberX / 2) * (x - this->NumberX / 2)) / (RkX * RkX);
  double TmpY = double((y - this->NumberY / 2) * (y - this->NumberY / 2)) / (RkY * RkY);
  if ((TmpX + TmpY) < 1.0)
    return true;
  else
    return false;  
}

// save the diagram of atoms in a file
//
// fileName = name of the file to stock the diagram

void EllipticalDotThreeDConstantCellPotential::SaveDiagram(char* fileName)
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

void EllipticalDotThreeDConstantCellPotential::LoadDiagram(char* fileName)
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

void EllipticalDotThreeDConstantCellPotential::SaveBmpPicture(char* fileName)
{
}
