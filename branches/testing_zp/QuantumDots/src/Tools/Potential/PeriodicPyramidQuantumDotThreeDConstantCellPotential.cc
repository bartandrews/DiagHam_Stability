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


#include "Tools/Potential/PeriodicPyramidQuantumDotThreeDConstantCellPotential.h"

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


// constructor from geometric parameters
//
// numberX, numberY, numberZ = number of cells in X, Y and Z directions respectively
// under = number of monolayers under the wetting layer
// above = number of monolayers above the dot
// wettingWidth= width of wetting layers in cell unit
// baseRadius = base radius of the truncated pyramid
// topRadius = base radius of the truncated pyramid

PeriodicPyramidQuantumDotThreeDConstantCellPotential::PeriodicPyramidQuantumDotThreeDConstantCellPotential(int numberX, int numberY, int numberZ, int under, int above, int wettingWidth, int baseRadius, int topRadius)
{
  this->NumberX = numberX;
  this->NumberY = numberY;
  this->NumberZ = numberZ;
  this->Under = under;
  this->Above = above;
  this->WettingWidth = wettingWidth;
  this->BaseRadius = baseRadius;
  this->TopRadius = topRadius;
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
  this->LaunchNumber = 0;
}

//destructor
//

PeriodicPyramidQuantumDotThreeDConstantCellPotential::~PeriodicPyramidQuantumDotThreeDConstantCellPotential()
{
  delete[] this->Alloy;
  delete[] this->PotentialValue;
}

// constructor the potential from piezoelectric parameter
//
// noInNProbability = probability for having InN without InN neighborhood
// withInNProbability = probability for having InN with InN neighborhood
// piezoField = difference of piezoelectric field
// cellSizeZ = height of a monolayer
// offset = offset bulk-alloy material
// scratch = true if constructed from nothing, else from existing Diagram
// fileName = file to store parameters

void PeriodicPyramidQuantumDotThreeDConstantCellPotential::ConstructPotential(double noInNProbability, double withInNProbability, double piezoField, double cellSizeZ, double offset, bool scratch, char* fileName)
{
  srand(time(NULL) + this->LaunchNumber);
  this->LaunchNumber += 2;

  if (scratch)
    for (int k = 0; k < this->NumberZ; ++k)      
      for (int j = 0; j < this->NumberY; ++j)
	for (int i = 0; i < this->NumberX; ++i)
	  this->Alloy[k][j][i] = false;

  int** Height = new int* [this->NumberX];
  int temp;
  for (int j = 0; j < this->NumberY; ++j)
    {
      Height[j] = new int [this->NumberX];
      for (int i = 0; i < this->NumberX; ++i)
	{
	  temp = 0;
	  for (int k = this->Under + this->WettingWidth; this->InTheDot(i, j, k); ++k)
	    ++temp;
	  Height[j][i] = temp + this->WettingWidth;
	}
    }

  double** Barrier = new double* [this->NumberY];
  double** Dot = new double* [this->NumberY];

  double*** ReferencePotential = new double** [this->NumberZ];
  for (int k = 0; k < this->NumberZ; ++k)
    {
      ReferencePotential[k] = new double* [this->NumberY];
      for (int j = 0; j < this->NumberY; ++j)
	ReferencePotential[k][j] = new double[this->NumberX];
    }

  for (int j = 0; j < this->NumberY; ++j)
    {
      Barrier[j] = new double [this->NumberX];
      Dot[j] = new double [this->NumberX];
      for (int i = 0; i < this->NumberX; ++i)
	{
	  Barrier[j][i] = piezoField * Height[j][i] / double(this->NumberZ);
	  Dot[j][i] = Barrier[j][i] - piezoField;
	  Barrier[j][i] *= cellSizeZ;
	  Dot[j][i] *= cellSizeZ;
	}
    }

  for (int j = 0; j < this->NumberY; ++j)
    for (int i = 0; i < this->NumberX; ++i)
      {
	if (scratch)
	  this->Alloy[0][j][i] = false;
	ReferencePotential[0][j][i] = 0.0;
	this->PotentialValue[0][j][i] = 0.0;
      }

  for (int k = 1; k < this->Under; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	{ 
	  if (scratch)
	    this->Alloy[k][j][i] = false;
	  ReferencePotential[k][j][i] = ReferencePotential[k - 1][j][i] + Barrier[j][i];
	  this->PotentialValue[k][j][i] = ReferencePotential[k - 1][j][i] + Barrier[j][i];
	}

  double TmpPotential;
  double inNCounter = 0.0;
  double dotCounter = 0.0;
  for (int k = this->Under; k < (this->Under + this->WettingWidth); ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	  {
	    ReferencePotential[k][j][i] = ReferencePotential[k - 1][j][i] + Dot[j][i];  
	    TmpPotential = ReferencePotential[k - 1][j][i] + Dot[j][i];
	    if (scratch)
	      if (this->Neighborhood(i, j, k))
		if ((double(rand())/RAND_MAX) < withInNProbability)		    
		  Alloy[k][j][i] = true;		  		    
		else
		  Alloy[k][j][i] = false;
	      else
		if ((double(rand())/RAND_MAX) < noInNProbability)	      		    
		  Alloy[k][j][i] = true;		  
		else
		  Alloy[k][j][i] = false;
	    dotCounter += 1.0;
	    if (this->Alloy[k][j][i])
	      {
		this->PotentialValue[k][j][i] = TmpPotential - offset;
		inNCounter += 1.0;
	      }
	    else
	      this->PotentialValue[k][j][i] = TmpPotential;
	  }

  for (int k = this->Under + this->WettingWidth; k < (this->NumberZ - this->Above); ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	{
	  if (InTheDot(i, j, k))
	    {
	      ReferencePotential[k][j][i] = ReferencePotential[k - 1][j][i] + Dot[j][i];
	      TmpPotential = ReferencePotential[k - 1][j][i] + Dot[j][i];
	      if (scratch)
		if (this->Neighborhood(i, j, k))
		  if ((double(rand())/RAND_MAX) < withInNProbability)		    
		    Alloy[k][j][i] = true;		  		    
		  else
		    Alloy[k][j][i] = false;
		else
		  if ((double(rand())/RAND_MAX) < noInNProbability)	      		    
		    Alloy[k][j][i] = true;		  
		  else
		    Alloy[k][j][i] = false;
	      dotCounter += 1.0;
	      if (this->Alloy[k][j][i])
		{
		  this->PotentialValue[k][j][i] = TmpPotential - offset;
		  inNCounter += 1.0;
		}
	      else
		this->PotentialValue[k][j][i] = TmpPotential;
	    }
	  else
	    {
	      ReferencePotential[k][j][i] = ReferencePotential[k - 1][j][i] + Barrier[j][i];
	      TmpPotential = ReferencePotential[k - 1][j][i] + Barrier[j][i];
	      if (scratch)
		this->Alloy[k][j][i] = false;
	      this->PotentialValue[k][j][i] = TmpPotential;

	    }
	}

  for (int k = this->NumberZ - this->Above; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	{
	  if (scratch)
	    this->Alloy[k][j][i] = false;
	  ReferencePotential[k][j][i] = ReferencePotential[k - 1][j][i] + Barrier[j][i];
	  TmpPotential = ReferencePotential[k - 1][j][i] + Barrier[j][i];
	  this->PotentialValue[k][j][i] = TmpPotential;
	}
  this->Concentration = inNCounter / dotCounter;
  delete[] ReferencePotential; delete[] Barrier; delete[] Dot; delete[] Height; 
}

// shift the potential with a given quantity
//
// delta = shift value

void PeriodicPyramidQuantumDotThreeDConstantCellPotential::ShiftPotential(double delta)
{
  for (int k = 0; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	this->PotentialValue[k][j][i] += delta;
}

// determine if there is any In in the first neigbor region (6 possibilities)
//
// m, n, h = three coordinations of the considered cell
// return = true if there is any InN neighborhood, false otherwise

bool PeriodicPyramidQuantumDotThreeDConstantCellPotential::Neighborhood(int m, int n, int h)
{
  int m1 = m - 1;
  int m2 = m + 1;
  int n1 = n - 1;
  int n2 = n + 1;
  int h1 = h - 1;
  int h2 = h + 1;
  if (m1 < 0) m1 = 0;
  if (m2 >= this->NumberX) m2 = this->NumberX - 1;
  if (n1 < 0) n1 = 0;
  if (n2 >= this->NumberY) n2 = this->NumberY - 1;
  if (h1 < 0) h1 = 0;
  if (h2 >= this->NumberZ) h2 = this->NumberZ - 1;

  for (int i = m1; i <= m2; ++i)
    if (this->Alloy[h][n][i])
      return true;
  for (int j = n1; j <= n2; ++j)
    if (this->Alloy[h][j][m])
      return true;
  for (int k = h1; k <= h2; ++k)
    if (this->Alloy[k][n][m])
      return true;

  return false;
}

// determine if a cell is in the dot or wetting layer
//
// x = x coordinate of the cell
// y = y coordinate of the cell
// z = z coordinate of the cell
// return = true if the cell is in the dot, false otherwise

bool PeriodicPyramidQuantumDotThreeDConstantCellPotential::InTheDot(int x, int y, int z)
{
  //position of base center
  int Cx = this->NumberX/2;
  int Cy = this->NumberY/2;

  if ((z >= this->Under) && (z < (this->Under + this->WettingWidth)))
    return true; // in the wetting layer

  if (z >= (this->Under + this->WettingWidth) && (z < (this->NumberZ - this->Above)))
    {
      double t1 = fabs((double)(x - Cx));
      double t2 = fabs((double)(y - Cy));
      double Rz = double(this->BaseRadius) - double(z - this->Under - this->WettingWidth + 1) * double(this->BaseRadius - this->TopRadius) / double(this->NumberZ - this->WettingWidth - this->Under - this->Above);
      double in1 = t1 - Rz + t2 / sqrt(3.0);
      double in2 = t2 - Rz * sqrt(3.0) / 2.0;
      // in the dot
      if ((in1 < 0) && (in2 < 0))
	return true;
      else
	return false;
    }
  return false;
}

// save the diagram of atoms in a file
//
// fileName = name of the file to stock the diagram

void PeriodicPyramidQuantumDotThreeDConstantCellPotential::SaveDiagram(char* fileName)
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

void PeriodicPyramidQuantumDotThreeDConstantCellPotential::LoadDiagram(char* fileName)
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

void PeriodicPyramidQuantumDotThreeDConstantCellPotential::SaveBmpPicture(char* fileName)
{
}
