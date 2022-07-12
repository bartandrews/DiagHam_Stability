////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//    class of dot potential in three directions with constant cylinders      //
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


#include "Tools/Potential/QuantumDotThreeDConstantCylinderPotential.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;


// constructor from data
//
// belowHeight = height of cylinder below the wetting layer
// wettingWidth = width of the wetting layer
// nbrCylinderDot = number of cylinders in the dot
// dotHeight = height of the dot
// baseRadius = radius of the dot base
// topRadius = radius of the dot top
// aboveHeight = height of cylinder above the dot
// underBarrier = height of the layer serving as well barrier
// superCylinderRadius = radius of the supercylinder 

QuantumDotThreeDConstantCylinderPotential::QuantumDotThreeDConstantCylinderPotential( double belowHeight, double wettingWidth, int nbrCylinderDot, double dotHeight, double baseRadius, double topRadius, double aboveHeight, double underBarrier, double superCylinderRadius)
{
  this->SuperCylinderRadius = superCylinderRadius;
  this->NbrCylinderDot = nbrCylinderDot;
  this->NumberZ = this->NbrCylinderDot + 4;
  this->CylinderHeight = new double[this->NumberZ];
  this->CylinderRadius = new double[this->NumberZ];
  this->PotentialValue =  new double[this->NumberZ];

  this->UnderBarrier = underBarrier;
  this->CylinderHeight[0] = this->UnderBarrier;
  this->CylinderRadius[0] = -1.0;

  this->BelowHeight = belowHeight;
  this->CylinderHeight[1] = this->BelowHeight;
  this->CylinderRadius[1] = -1.0;

  this->WettingWidth = wettingWidth;
  this->CylinderHeight[2] = this->WettingWidth;
  this->CylinderRadius[2] = -1.0;

  this->DotHeight = dotHeight;
  this->BaseRadius = baseRadius;
  this->TopRadius = topRadius;
  for (int k = 3; k < (3 + this->NbrCylinderDot); ++k)
    {
      this->CylinderHeight[k] = this->DotHeight / double(this->NbrCylinderDot);
      this->CylinderRadius[k] = this->BaseRadius - double(k - 3) * (this->BaseRadius - this->TopRadius) / double(this->NbrCylinderDot - 1);
    }
  
  this->AboveHeight = aboveHeight;
  this->CylinderHeight[this->NumberZ - 1] = this->AboveHeight;
  this->CylinderRadius[this->NumberZ - 1] = -1.0;

}

// destructor
//

QuantumDotThreeDConstantCylinderPotential::~QuantumDotThreeDConstantCylinderPotential()
{
}

// construct the potential profile from potential values
//
// dotPotential = dot potential, in comparision with the bulk
// wellPotential = barrier potential, in comparison with the bulk

void QuantumDotThreeDConstantCylinderPotential::ConstructPotential(double dotPotential, double wellPotential)
{
  this->PotentialValue[0] = wellPotential;
  this->PotentialValue[1] = 0.0;
  this->PotentialValue[2] = dotPotential;
  for (int k = 0; k < this->NbrCylinderDot; ++k)
    this->PotentialValue[k + 3] = dotPotential;
  this->PotentialValue[this->NumberZ - 1] = 0.0;
}

// add two barriers of the same potential just below and above the dot + WL (for Vienna experimenters)
//
// belowBarrier = width of barrier layer just below the WL (in Angstrom unit)
// aboveBarrier = width of barrier layer just above  the dot (in Angstrom unit)
// potential = potential to add in these barriers (in eV unit)  

void QuantumDotThreeDConstantCylinderPotential::AddBarrierPotential (double belowBarrier, double aboveBarrier, double potential)
{
  if ((belowBarrier > this->UnderBarrier) || (aboveBarrier > this->AboveHeight))
    {
      cout << "The barrier potential inserted is too thick. Try another thinner one" << endl;
      return;
    }

  delete[] this->CylinderHeight; delete[] this->CylinderRadius; 
  double* tmpPotential = this->PotentialValue;

  this->NumberZ = this->NbrCylinderDot + 6;
  this->CylinderHeight = new double[this->NumberZ];
  this->CylinderRadius = new double[this->NumberZ];
  this->PotentialValue =  new double[this->NumberZ];
  
  // geometry
  this->CylinderHeight[0] = this->UnderBarrier;
  this->CylinderRadius[0] = -1.0;

  this->CylinderHeight[1] = this->BelowHeight - belowBarrier;
  this->CylinderRadius[1] = -1.0;

  this->CylinderHeight[2] = belowBarrier;
  this->CylinderRadius[2] = -1.0;

  this->CylinderHeight[3] = this->WettingWidth;
  this->CylinderRadius[3] = -1.0;

  for (int k = 4; k < (4 + this->NbrCylinderDot); ++k)
    {
      this->CylinderHeight[k] = this->DotHeight / double(this->NbrCylinderDot);
      this->CylinderRadius[k] = this->BaseRadius - double(k - 4) * (this->BaseRadius - this->TopRadius) / double(this->NbrCylinderDot - 1);
    }

  this->CylinderHeight[this->NumberZ - 2] = aboveBarrier;
  this->CylinderRadius[this->NumberZ - 2] = -1.0;

  this->CylinderHeight[this->NumberZ - 1] = this->AboveHeight - aboveBarrier;
  this->CylinderRadius[this->NumberZ - 1] = -1.0;
  
  // potential
  this->PotentialValue[0] = tmpPotential[0];
  this->PotentialValue[1] = tmpPotential[1];
  this->PotentialValue[2] = this->PotentialValue[1] + potential;
  this->PotentialValue[3] = tmpPotential[2];
  for (int k = 0; k < this->NbrCylinderDot; ++k)
    this->PotentialValue[k + 4] = tmpPotential[k + 3];
  this->PotentialValue[this->NumberZ - 2] = tmpPotential[1] + potential;
  this->PotentialValue[this->NumberZ - 1] = this->PotentialValue[1];

  delete[] tmpPotential;

}

// save the whole diagram presentation in a bitmap file
//
// fileName = name of the file to stock the diagram presentation

void QuantumDotThreeDConstantCylinderPotential::SaveBmpPicture(char* fileName)
{
}

