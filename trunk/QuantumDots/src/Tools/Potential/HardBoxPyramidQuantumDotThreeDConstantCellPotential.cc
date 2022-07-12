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


#include "Tools/Potential/HardBoxPyramidQuantumDotThreeDConstantCellPotential.h"

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

HardBoxPyramidQuantumDotThreeDConstantCellPotential::HardBoxPyramidQuantumDotThreeDConstantCellPotential(int numberX, int numberY, int numberZ, int under, int above, int wettingWidth, int baseRadius, int topRadius)
{
  this->NumberX = numberX;
  this->NumberY = numberY;
  this->NumberZ = numberZ;
  this->Under = under;
  this->Above = above;
  this->WettingWidth = wettingWidth;
  this->BaseRadius = baseRadius;
  this->TopRadius = topRadius;
  this->UnderSize = new double[this->Under];
  this->UnderPotentialValue = new double[this->Under];
  this->AboveSize = new double[this->Above];
  this->AbovePotentialValue = new double[this->Above];
  this->Alloy = new bool** [this->NumberZ];
  for (int k = 0; k < this->NumberZ; ++k)
    {
      this->Alloy[k] = new bool* [this->NumberY];
      for (int j = 0; j < this->NumberY; ++j)
	this->Alloy[k][j] = new bool [this->NumberX];
    }
  this->PotentialValue = new double** [this->NumberZ - this->Under - this->Above];
  for (int k = 0; k < (this->NumberZ - this->Under - this->Above); ++k)
    {
      this->PotentialValue[k] = new double* [this->NumberY];
	for (int j = 0; j < this->NumberY; ++j)
	  this->PotentialValue[k][j] = new double [this->NumberX];
    }
}

//destructor
//

HardBoxPyramidQuantumDotThreeDConstantCellPotential::~HardBoxPyramidQuantumDotThreeDConstantCellPotential()
{
  delete[] this->Alloy;
  delete[] this->PotentialValue;
  delete[] this->UnderSize;
  delete[] this->UnderPotentialValue;
  delete[] this->AboveSize;
  delete[] this->AbovePotentialValue;
}

// contruct the potential from electric fields
//
// noInNProbability = probability for having InN without InN neighborhood
// withInNProbability = probability for having InN with InN neighborhood
// downField, wettingField, dotField, upField = the electric field in respective parts
// offset = offset bulk-alloy material
// cellSizeZ = height of a monolayer
// scratch = true if constructed from nothing, else from existing Diagram
// fileName = file to store parameters

void HardBoxPyramidQuantumDotThreeDConstantCellPotential::ConstructPotential(double noInNProbability, double withInNProbability, double downField, double wettingField, double dotField, double upField, double offset, double cellSizeZ, bool scratch, char* fileName)
{
  srand(time(NULL));
  //position of base center
  double tmp1 = 0.0, tmp2 = 0.0;

  // under the wetting layer
  for (int k = 0; k < this->Under; ++k)
    {
      tmp1 = downField * (k + 1 - this->Under) * cellSizeZ;
      if (scratch)
        for (int j = 0; j < NumberY; ++j)
	  for (int i = 0; i < NumberX; ++i)
	    this->Alloy[k][j][i] = false;

      this->UnderPotentialValue[k] = tmp1;
      this->UnderSize[k] = cellSizeZ;
    }
  double inNCounter = 0.0;
  double dotCounter = 0.0;
  // wetting layer
  for (int k = this->Under; k < (this->Under + this->WettingWidth); ++k)
    {
      tmp1 = (k + 1 - this->Under) * wettingField * cellSizeZ;
      tmp2 = tmp1 - offset;
      for (int j = 0; j < NumberY; ++j)
	for (int i = 0; i < NumberX; ++i)
	  {
	    if (scratch)
	      {
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
	      }
	    dotCounter += 1.0;
	    if (Alloy[k][j][i])
	      {
		this->PotentialValue[k - this->Under][j][i] = tmp2;
		inNCounter += 1.0;
	      }
	    else
	      this->PotentialValue[k - this->Under][j][i] = tmp1;
	  }
    }

  // quantum dot
  bool** base = new bool* [NumberY];
  bool** interface = new bool* [NumberY];
  for (int i = 0; i < NumberY; ++i)
    {
      base[i] = new bool [NumberX];
      interface[i] = new bool [NumberX];
    }

  for (int j = 0; j < NumberY; ++j)
    for (int i = 0; i < NumberX; ++i)
      {
	interface[j][i] = false;
	base[j][i] = false;
      }

  // the first mono-layer of the quantum dot
  tmp1 = (wettingField * this->WettingWidth + dotField) * cellSizeZ; tmp2 = tmp1 - offset; double tmp3 = wettingField * this->WettingWidth * cellSizeZ;
  for (int j = 0; j < NumberY; ++j)
    for (int i = 0; i < NumberX; ++i)
      {
	// in the dot
	if (this->InTheDot(i, j, this->Under + this->WettingWidth))
	  {
	    dotCounter += 1.0;
	    if (scratch)
	      {
		if (this->Neighborhood(i, j, this->Under + this->WettingWidth))
		  if ((double(rand())/RAND_MAX) < withInNProbability)		    
		    Alloy[this->Under + this->WettingWidth][j][i] = true;		     		    
		  else
		    Alloy[this->Under + this->WettingWidth][j][i] = false;
		else
		  if ((double(rand())/RAND_MAX) < noInNProbability)		    
		    Alloy[this->Under + this->WettingWidth][j][i] = true;				    
		  else
		    Alloy[this->Under + this->WettingWidth][j][i] = false;
	      }

	    if (Alloy[this->Under + this->WettingWidth][j][i])
	      {
		this->PotentialValue[this->WettingWidth][j][i] = tmp2;
		inNCounter += 1.0;
	      }
	    else
	      this->PotentialValue[this->WettingWidth][j][i] = tmp1;
	    base[j][i] = true;
	    interface[j][i] = true;
	  }
	// outside the dot
	else
	  {
	    this->PotentialValue[this->WettingWidth][j][i] = tmp3;
	    if (scratch)
	      Alloy[this->Under + this->WettingWidth][j][i] = false;
	  }
      }

  // the following mono-layers in the dot
  for (int k = this->Under + this->WettingWidth + 1; k < NumberZ - this->Above; ++k)
    {
      tmp1 = (wettingField * this->WettingWidth + dotField * (k + 1 - this->Under - this->WettingWidth)) * cellSizeZ; tmp2 = tmp1 - offset;
      for (int j = 0; j < NumberY; ++j)
	for (int i = 0; i < NumberX; ++i)
	  {
	    // outside the base
	    if (!base[j][i])
	      {
		if (scratch)
		  Alloy[k][j][i] = false;
		this->PotentialValue[k - this->Under][j][i] = tmp3;
	      }
	    // inside the base
	    else
	      {
		// in the dot
		if (this->InTheDot(i, j, k))
		  {
		    dotCounter += 1.0;
		    if (scratch)
		      {
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
		      }
		    if (Alloy[k][j][i])
		      {
			this->PotentialValue[k - this->Under][j][i] = tmp2;
			inNCounter += 1.0;
		      }
		    else
		      this->PotentialValue[k - this->Under][j][i] = tmp1;
		  }
		// between the dot and the base
		else
		  {
		    if (interface[j][i])
		      {
			dotCounter += 1.0;
			interface[j][i] = false;
			if (scratch)
			  {
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
			  }
			
			if (Alloy[k][j][i])
			  {
			    this->PotentialValue[k - this->Under][j][i] = tmp2;
			    inNCounter += 1.0;
			  }
			else
			  this->PotentialValue[k - this->Under][j][i] = tmp1;
		      }
		    else
		      {
			this->PotentialValue[k - this->Under][j][i] = this->PotentialValue[k - this->Under - 1][j][i];
			if (scratch)
			  Alloy[k][j][i] = false;
		      }
		  }
	      }
	  }
    }

  int tmpIndice = this->NumberZ - this->Above;
  // above the dot
  for (int k = 0; k < this->Above; ++k)
    {
      tmp1 = (wettingField * this->WettingWidth + dotField * (NumberZ - this->WettingWidth - this->Under - this->Above) + upField * (k + 1)) * cellSizeZ;
      if (scratch)
	for (int j = 0; j < NumberY; ++j)
	  for (int i = 0; i < NumberX; ++i)
	    Alloy[tmpIndice][j][i] = false;
      ++tmpIndice;
      this->AbovePotentialValue[k] = tmp1;
      this->AboveSize[k] = cellSizeZ;

    }

  for (int i = 0; i < NumberY; ++i)
    {
      delete[] base[i];
      delete[] interface[i];
    }
  delete[] base;
  delete[] interface;

  double percent =  inNCounter / dotCounter;
  this->Concentration = percent;
  if (fileName != 0)
    {
      ofstream File;
      File.open(fileName, ios::out | ios::app);
      File << "Pyramid quantum dot parameters\n";
      File << "Dimension (cells): X = " << this->NumberX << ", Y = " << this->NumberY << ", Z = " << this->NumberZ << '\n';
      File << "Number of cells under and above the dot: Under = " << this->Under << ", Above = " << this->Above << '\n';
      File << "Probability given: " << noInNProbability << ", " << withInNProbability << '\n';
      File << "Real proportion: " << percent << '\n';
      File << "Base and top radius: Base = " << this->BaseRadius << ", Top = " << this->TopRadius << '\n';
      File << "Wetting layer thickness: " << this->WettingWidth << '\n';
      File << "Electric fields under the dot, in the wetting layers, in and above the dot: " << downField << ", " << wettingField << ", " << dotField << ", " << upField << '\n';
      File << "Band offset: " << offset << '\n';
      File << "Growth direction lattice constant: " << cellSizeZ << endl;
      File.close();
    }

  cout << "Real proportion: " << percent << endl;  
}

// shift the potential with a given quantity
//
// delta = shift value

void HardBoxPyramidQuantumDotThreeDConstantCellPotential::ShiftPotential(double delta)
{
  for (int k = 0; k < this->Under; ++k)
    this->UnderPotentialValue[k] += delta;
  for (int k = 0; k < (this->NumberZ - this->Under - this->Above); ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	this->PotentialValue[k][j][i] += delta;
  for (int k = 0; k < this->Above; ++k)
    this->AbovePotentialValue[k] += delta;
}

// determine if there is any In in the first neigbor region (6 possibilities)
//
// m, n, h = three coordinations of the considered cell
// return = true if there is any InN neighborhood, false otherwise

bool HardBoxPyramidQuantumDotThreeDConstantCellPotential::Neighborhood(int m, int n, int h)
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

bool HardBoxPyramidQuantumDotThreeDConstantCellPotential::InTheDot(int x, int y, int z)
{
  //position of base center
  int Cx = this->NumberX/2;
  int Cy = this->NumberY/2;

  if ((z >= this->Under) && (z < (this->Under + this->WettingWidth)))
    return true; // in the wetting layer

  if (z >= (this->Under + this->WettingWidth) && (z < (this->NumberZ - this->Above)))
    {
      /* the code in the if segment below is used for the old potential, before February 2004 */
      /*
      if (z == this->Under + this->WettingWidth)
	{
	  double v1 = fabs((double)(x - Cx));
	  double v2 = fabs((double)(y - Cy));
	  double tin1 = v1 - this->BaseRadius + v2 / sqrt(3.0);
	  double tin2 = v2 - this->BaseRadius * sqrt(3.0) / 2.0;
	  // in the dot
	  if ((tin1 < 0) && (tin2 < 0))
	    return true;
	  else 
	    return false;
	}
      */

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

void HardBoxPyramidQuantumDotThreeDConstantCellPotential::SaveDiagram(char* fileName)
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

void HardBoxPyramidQuantumDotThreeDConstantCellPotential::LoadDiagram(char* fileName)
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

// save the potential in a file as a reduced form
//
// fileName = name of the file to stock the potential

void HardBoxPyramidQuantumDotThreeDConstantCellPotential::SavePotentialWithConstantField(char* fileName)
{
  ofstream file(fileName);
  if (!file.is_open())
    cout << "Error when open the file: " << fileName << " to write diagram." << endl;
  file << this->Under << " ";
  for (int k = 0; k < this->Under; ++k)
    file << this->UnderPotentialValue[k] << " " << this->UnderSize[k] << " ";
  file << '\n';
  file << this->Above << " ";
  for (int k = 0; k < this->Above; ++k)
    file << this->AbovePotentialValue[k] << " " << this->AboveSize[k] << " ";
  file << '\n';
  for (int k = 0; k < (this->NumberZ - this->Above - this->Under); ++k)
    {
      for (int j = 0; j < this->NumberY; ++j)
	{
	  for (int i = 0; i < this->NumberX; ++i)
	    file << this->PotentialValue[k][j][i] << '\t';
	  file << '\n';
	}
      file << '\n';
    }
  file.close();
}

// load the potential from a file as a reduced form
//
// fileName = name of the potential file

void HardBoxPyramidQuantumDotThreeDConstantCellPotential::LoadPotentialWithConstantField(char* fileName)
{
  ifstream file(fileName);
  if (! file.is_open())
    {
      cout << "Error when open the potential file. Exit now" << endl;
      exit(1);
    }
  file >> this->Under;
  for (int k = 0; k < this->Under; ++k)
    file >> this->UnderPotentialValue[k] >> this->UnderSize[k];
  file >> this->Above;
  for (int k = 0; k < this->Above; ++k)
    file >> this->AbovePotentialValue[k] >> this->AboveSize[k];

  for (int k = 0; k < (this->NumberZ - this->Above - this->Under); ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	file >> this->PotentialValue[k][j][i];
  file.close();
}

// save the whole diagram presentation in a bitmap file
//
// fileName = name of the file to stock the diagram presentation

void HardBoxPyramidQuantumDotThreeDConstantCellPotential::SaveBmpPicture(char* fileName)
{
}
