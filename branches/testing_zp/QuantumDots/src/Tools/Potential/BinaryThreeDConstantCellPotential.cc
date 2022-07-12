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


#include "Tools/Potential/BinaryThreeDConstantCellPotential.h"
#include "BitmapTools/BitmapPicture/BmpFormat.h"
#include "BitmapTools/Color/PicRGB.h"
#include "GeneralTools/Endian.h"

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


// default constructor
//

BinaryThreeDConstantCellPotential::BinaryThreeDConstantCellPotential ()
{
}

// constructor
//
// numberX, numberY, numberZ = number of cells in X, Y and Z directions respectively

BinaryThreeDConstantCellPotential::BinaryThreeDConstantCellPotential(int numberX, int numberY, int numberZ)
{
  this->NumberX = numberX;
  this->NumberY = numberY;
  this->NumberZ = numberZ;
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
}

// destructor
//

BinaryThreeDConstantCellPotential::~BinaryThreeDConstantCellPotential()
{
  delete[] this->Alloy;
  delete[] this->PotentialValue;
}

// save the diagram of atoms in a file
//
// fileName = name of the file to stock the diagram

void BinaryThreeDConstantCellPotential::SaveDiagram(char* fileName)
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

void BinaryThreeDConstantCellPotential::LoadDiagram(char* fileName)
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

// save the potential description of atoms in a file (binary mode)
//
// fileName = name of the file to stock the potential description

void BinaryThreeDConstantCellPotential::SaveBinaryPotential(char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (! File.is_open())
    {
      cout << "Error when open the potential description  file: " << fileName << " . Exit now" << endl;
      exit(1);
    }
  for (int k = 0; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	WriteLittleEndian(File, this->PotentialValue[k][j][i]);
  File.close();
}

// load the potential description of atoms from a file (binary mode)
//
// fileName = name of the file in which the potential description is stocked

void BinaryThreeDConstantCellPotential::LoadBinaryPotential(char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (! File.is_open())
    {
      cout << "Error when open the potential description  file: " << fileName << " . Exit now" << endl;
      exit(1);
    }
  for (int k = 0; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	ReadLittleEndian(File, this->PotentialValue[k][j][i]);
  File.close();
}

// save the whole diagram presentation in a bitmap file
//
// fileName = name of the file to stock the diagram presentation

void BinaryThreeDConstantCellPotential::SaveBmpPicture(char* fileName)
{
  PicRGB background (255, 255, 255);
  PicRGB InN (0, 0, 0);
  PicRGB GaN (231, 231, 231);  
  if (!this->SaveBmpPicture(0, 0, 0, this->NumberX, 0, this->NumberY, 1, 5, 5, InN, GaN, background, 5, fileName));
    cout << "Error when saving the bitmap file: " << fileName << endl;
}

// save the whole diagram presentation in a bitmap file
//
// sizeX = the size in X direction of each cell in pixel
// sizeY = the size in Y direction of each cell in pixel
// fileName = name of the file to stock the diagram presentation

void BinaryThreeDConstantCellPotential::SaveBmpPicture(char* fileName, int sizeX, int sizeY)
{
  PicRGB background (255, 255, 255);
  PicRGB InN (0, 0, 0);
  PicRGB GaN (231, 231, 231);  
  if (!this->SaveBmpPicture(0, 0, 0, this->NumberX, 0, this->NumberY, 1, sizeX, sizeY, InN, GaN, background, 5, fileName));
    cout << "Error when saving the bitmap file: " << fileName << endl;
}

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

bool BinaryThreeDConstantCellPotential::SaveBmpPicture(int under, int above, int startX, int endX, int startY, int endY, int choice, int sizeX, int sizeY, PicRGB& InN, PicRGB& GaN, PicRGB& background, int NbrX, char* fileName)
{
  int Layers = 0;
  switch(choice)
    {
    case 1:
      Layers = this->NumberZ - under - above;
      break;
    case 2:
      Layers = this->NumberY - under - above;
      break;
    default:
      Layers = this->NumberX - under - above;
    }
  int NbrY = (Layers - 1) / NbrX + 1;
  int Border = 5;
  int tmpLx = (endX - startX) * sizeX + Border;
  int tmpLy = (endY - startY) * sizeY + Border;
  int Lx = tmpLx * NbrX - Border;
  int Ly = tmpLy * NbrY - Border;

  BmpFormat* picture = new BmpFormat(Lx, Ly, background);

  int tmpX, tmpY, tmpZ, originX, originY, tmpOriginX;
  bool dopage;
  for (int k = 0; k < Layers; ++k)
    {
      tmpY = k / NbrX;
      tmpX = k - tmpY * NbrX;
      tmpZ = k + under;
      originX = tmpX * tmpLx;
      originY = tmpY * tmpLy;

      for (int j = startY; j < endY; ++j)
	{
	  tmpOriginX = originX;
	  for (int i = startX; i < endX; ++i)
	    {
	      switch(choice)
		{
		case 1:
		  dopage = this->Alloy[tmpZ][j][i];
		  break;
		case 2:
		  dopage = this->Alloy[j][tmpZ][i];
		default:
		  dopage = this->Alloy[j][i][tmpZ];
		}
	      if (dopage)
		CellFill(tmpOriginX, sizeX, originY, sizeY, InN, picture);
	      else
		CellFill(tmpOriginX, sizeX, originY, sizeY, GaN, picture);

	      tmpOriginX += sizeX;
	    }
	  originY += sizeY;
	}
    }
  
  picture->SavePicture(fileName);

  return true;
}

// fill a cell with a given color
//
// startX: X coordination to start
// sizeX:  size in X direction
// startY: Y coordination to start
// sizeY:  size in Y direction
// Col: color in RGB definition
// picture: pointer to picture which is filled

void CellFill(int startX, int sizeX, int startY, int sizeY, PicRGB& Col, AbstractBitmapPicture* picture)
{
  for (int i = startX; i < startX + sizeX; ++i)
    for (int j = startY; j < startY + sizeY; ++j)
      picture->SetPixel(i, j, Col);
}
