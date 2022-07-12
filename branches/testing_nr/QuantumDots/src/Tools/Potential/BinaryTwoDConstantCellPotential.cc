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


#include "Tools/Potential/BinaryTwoDConstantCellPotential.h"
#include "BitmapTools/BitmapPicture/BmpFormat.h"
#include "BitmapTools/Color/PicRGB.h"
#include "Tools/Potential/ThreeDConstantCellPotential.h"
#include "Vector/RealVector.h"

#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <math.h>

using std::ifstream;
using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


// constructor
//
// numberX, numberY = number of cells in X and Y directions respectively

BinaryTwoDConstantCellPotential::BinaryTwoDConstantCellPotential(int numberX, int numberY)
{
  this->NumberX = numberX;
  this->NumberY = numberY;
  this->Alloy = new bool* [this->NumberY];
  this->PotentialValue = new double* [this->NumberY];
  for (int j = 0; j < this->NumberY; ++j)
    {
      this->Alloy[j] = new bool [this->NumberX];
      this->PotentialValue[j] = new double [this->NumberX];
    }    
}

// destructor
//

BinaryTwoDConstantCellPotential::~BinaryTwoDConstantCellPotential()
{
  delete[] this->Alloy;
  delete[] this->PotentialValue;
}

// construct a potential by averaging a 3D potential in the Z direction
// 
// potential = 3D potential with constant value in each cell
// coefficient = the coefficient for each monolayer

void BinaryTwoDConstantCellPotential::ConstructEffectivePotential(ThreeDConstantCellPotential* potential, RealVector& coefficient)
{
  if ((potential->GetNbrCellX() != this->NumberX) || (potential->GetNbrCellY() != this->NumberY))
    cout << "The dimensions in X or Y direction are not the same in 3D and 2D potentials" << endl;
  if (coefficient.GetVectorDimension() != potential->GetNbrCellZ())
    cout << "The dimension in Z of coeffients and 3D potential are not the same" << endl;
  int Height = potential->GetNbrCellZ();
  double tmp = 0.0;
  for (int j = 0; j < this->NumberY; ++j)
    for (int i = 0; i < this->NumberX; ++i)
      {
	tmp = 0.0;
	for (int k = 0; k < Height; ++k)
	  tmp += (potential->GetPotential(i, j, k) * coefficient[k]);
	this->PotentialValue[j][i] = tmp;
      }  
}

// save the diagram of atoms in a file
//
// fileName = name of the file to stock the diagram

// shift the potential with a given quantity
//
// delta = shift value

void BinaryTwoDConstantCellPotential::ShiftPotential(double delta)
{
  for (int j = 0; j < this->NumberY; ++j)
    for (int i = 0; i < this->NumberX; ++i)
      this->PotentialValue[j][i] += delta;
} 

void BinaryTwoDConstantCellPotential::SaveDiagram(char* fileName)
{
  ofstream file(fileName);
  if (!file.is_open())
    cout << "Error when open the file: " << fileName << " to write diagram." << endl;
  for (int j = 0; j < this->NumberY; ++j)
    {
      for (int i = 0; i < this->NumberX; ++i)
	file << this->Alloy[j][i] << " ";
      file << '\n';
    }
  file.close();
}

// load the diagram of atoms from a file
//
// fileName = name of the file in which the diagram is stocked

void BinaryTwoDConstantCellPotential::LoadDiagram(char* fileName)
{
  ifstream file(fileName);
  if (! file.is_open())
    {
      cout << "Error when open the diagram file: " << fileName << " . Exit now" << endl;
      exit(1);
    }
  for (int j = 0; j < this->NumberY; ++j)
    for (int i = 0; i < this->NumberX; ++i)
      file >> this->Alloy[j][i];
  file.close();
}

// save the whole diagram presentation in a bitmap file
//
// fileName = name of the file to stock the diagram presentation

void BinaryTwoDConstantCellPotential::SaveBmpPicture(char* fileName)
{
}
