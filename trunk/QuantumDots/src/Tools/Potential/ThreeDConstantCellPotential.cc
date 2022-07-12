////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//         class of potential in three directions with constant cells         //
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


#include "Tools/Potential/ThreeDConstantCellPotential.h"
#include "GeneralTools/Endian.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;


// destructor
//

ThreeDConstantCellPotential::~ThreeDConstantCellPotential()
{
}

// shift the potential with a given quantity
//
// delta = shift value

void ThreeDConstantCellPotential::ShiftPotential(double delta)
{
}

// save the diagram of atoms in a file
//
// fileName = name of the file to stock the diagram

void ThreeDConstantCellPotential::SaveDiagram(char* fileName)
{
}

// load the diagram of atoms from a file
//
// fileName = name of the file in which the diagram is stocked

void ThreeDConstantCellPotential::LoadDiagram(char* fileName)
{
}

// save the potential in a file
//
// fileName = name of the file to stock the potential
// the file contains NumberZ blocks, each block has NumberY lines and NumberX columns

void ThreeDConstantCellPotential::SavePotential(char* fileName)
{
  ofstream file(fileName);
  if (!file.is_open())
    cout << "Error when open the file: " << fileName << " to write potential." << endl;
  double tmp = 0.0;
  for (int k = 0; k < this->NumberZ; ++k)
    {
      for (int j = 0; j < this->NumberY; ++j)
	{
	  for (int i = 0; i < this->NumberX; ++i)
	    {
	      tmp = this->GetPotential(i, j, k);
	      file << tmp << '\t';	  
	    }
	  file << '\n';
	}
      file << '\n';
    }
  file.close();
}

// load the potential from a file
//
// fileName = name of the file in which the potential is stocked

void ThreeDConstantCellPotential::LoadPotential(char* fileName)
{
  ifstream file(fileName);
  if (!file.is_open())
    cout << "Error when open the file: " << fileName << " to load potential." << endl;
  double tmp = 0.0;
  for (int k = 0; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)   
      for (int i = 0; i < this->NumberX; ++i)
	{
	  file >> tmp;
	  this->SetPotential(i, j, k, tmp);	  
	}  
  file.close();
}

// save the potential description of atoms in a file (binary mode)
//
// fileName = name of the file to stock the potential description

void ThreeDConstantCellPotential::SaveBinaryPotential(char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (! File.is_open())
    {
      cout << "Error when open the potential description  file: " << fileName << " . Exit now" << endl;
      exit(1);
    }
  double Tmp;
  for (int k = 0; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	{
	  Tmp = this->GetPotential(i, j, k);
	  WriteLittleEndian(File, Tmp);
	}
  File.close();
}

// load the potential description of atoms from a file (binary mode)
//
// fileName = name of the file in which the potential description is stocked

void ThreeDConstantCellPotential::LoadBinaryPotential(char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (! File.is_open())
    {
      cout << "Error when open the potential description  file: " << fileName << " . Exit now" << endl;
      exit(1);
    }
  double Tmp;
  for (int k = 0; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	{
	  ReadLittleEndian(File, Tmp);
	  this->SetPotential(i, j, k, Tmp);
	}
  File.close();
}

// save the whole diagram presentation in a bitmap file
//
// fileName = name of the file to stock the diagram presentation

void ThreeDConstantCellPotential::SaveBmpPicture(char* fileName)
{
}
