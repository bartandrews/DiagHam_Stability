////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//      class of potential in three directions with constant cylinders        //
//                                                                            //
//                        last modification : 04/22/2004                      //
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


#include "Tools/Potential/ThreeDConstantCylinderPotential.h"
#include "Tools/Potential/OneDConstantCellPotential.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

// destructor
//

ThreeDConstantCylinderPotential::~ThreeDConstantCylinderPotential()
{
  delete[] this->PotentialValue;
}

// shift the potential with a given quantity
//
// delta = shift value

void ThreeDConstantCylinderPotential::ShiftPotential(double delta)
{
  for (int k = 0; k < this->NumberZ; ++k)
    this->PotentialValue[k] += delta;
}

// save the diagram of atoms in a file
//
// fileName = name of the file to stock the diagram

void ThreeDConstantCylinderPotential::SaveDiagram(char* fileName)
{
}

// load the diagram of atoms from a file
//
// fileName = name of the file in which the diagram is stocked

void ThreeDConstantCylinderPotential::LoadDiagram(char* fileName)
{
}

// save the potential in a file
//
// fileName = name of the file to stock the potential
// the file contains NumberZ blocks, each block has NumberY lines and NumberX columns

void ThreeDConstantCylinderPotential::SavePotential(char* fileName)
{
}

// load the potential from a file
//
// fileName = name of the file in which the potential is stocked

void ThreeDConstantCylinderPotential::LoadPotential(char* fileName)
{
}

// get a 1-D potential by averaging through a gaussian weight
//
// sigma = the variance of the gaussian function
// m = value of angular momentum
// return = pointer to the one dimension potential

OneDConstantCellPotential* ThreeDConstantCylinderPotential::GaussianReductionOneDimension(double sigma, int m)
{
  double* cellWidth = new double [this->NumberZ];
  double* potentialValue = new double [this->NumberZ];
  if (m < 0)
    m = -m;
  for (int k = 0; k < this->NumberZ; ++k)
    {
      double tmp = this->CylinderRadius[k] * this->CylinderRadius[k] / (sigma * sigma);
      double tmp1 = 1.0, tmp2 = 1.0;
      for (int i = 0; i < m; ++i)
	{
	  tmp2 *= tmp;
	  tmp1 += tmp2;
	}
      cellWidth[k] = this->CylinderHeight[k];      
      if (this->CylinderRadius[k] >= 0.0)
	potentialValue[k] = this->PotentialValue[k] * (1.0 - exp(-tmp) * tmp1);
      else
	potentialValue[k] = this->PotentialValue[k];
      //cout << potentialValue[k] << endl;
    }
  OneDConstantCellPotential* potential = new OneDConstantCellPotential (this->NumberZ, cellWidth, potentialValue);
  return potential;
}

// get a 1-D potential by averaging through a weight of C(r^2 - a^2) exp(-r^2/(2s^2))
//
// lambda = the variance of the gaussian function
// sigma = the variance of the gaussian function for the 1S state
// m = value of angular momentum
// return = pointer to the one dimension potential

OneDConstantCellPotential* ThreeDConstantCylinderPotential::SecondDegreeReductionOneDimension(double lambda, double sigma, int m)
{
  double* cellWidth = new double [this->NumberZ];
  double* potentialValue = new double [this->NumberZ];
  if (m < 0)
    m = -m;

  // only the case m = 0
  for (int k = 0; k < this->NumberZ; ++k)
    {
      if (this->CylinderRadius[k] >= 0.0)
	{
	  double tmp = this->CylinderRadius[k] * this->CylinderRadius[k] / (lambda * lambda);
	  
	  double l = lambda * lambda; double s = sigma * sigma;
	  double a = 2 * l * s / (l + s);
	  
	  double numerator = l * l * (tmp * tmp + 2 * tmp + 2) - 2 * l * a * (tmp + 1) + a * a;
	  double denominator = 2 * l * l - 2 * l * a + a * a;

	  potentialValue[k] = this->PotentialValue[k] * (1 - exp(-tmp) * numerator / denominator);
	}
      else
	potentialValue[k] = this->PotentialValue[k];
      cellWidth[k] = this->CylinderHeight[k];      

      //cout << potentialValue[k] << endl;
    }  

  OneDConstantCellPotential* potential = new OneDConstantCellPotential (this->NumberZ, cellWidth, potentialValue);
  return potential;
}
