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


#ifndef THREEDCONSTANTCYLINDERPOTENTIAL_H
#define THREEDCONSTANTCYLINDERPOTENTIAL_H


#include "config.h"
#include "Tools/Potential/AbstractPotential.h"

class OneDConstantCellPotential;

class ThreeDConstantCylinderPotential : public AbstractPotential
{

 protected:

  // radius of the supercylinder
  double SuperCylinderRadius;

  // number of cylinders in Z direction
  int NumberZ;

  // height of each cylinder
  double* CylinderHeight;

  // radius of each cylinder, equal to -1.0 if infinite
  double* CylinderRadius;

  // value of constant potential in each cylinder
  double* PotentialValue;

 public:

  // destructor
  //
  virtual ~ThreeDConstantCylinderPotential();

  // get the super cylinder radius
  //
  double GetSuperCylinderRadius();

  // get the number of cylinders in Y direction
  //
  int GetNbrCylinderZ();

  // shift the potential with a given quantity
  //
  // delta = shift value
  virtual void ShiftPotential(double delta);

  // save the diagram of atoms in a file
  //
  // fileName = name of the file to stock the diagram
  virtual void SaveDiagram(char* fileName);

  // load the diagram of atoms from a file
  //
  // fileName = name of the file in which the diagram is stocked
  virtual void LoadDiagram(char* fileName);

  // save the potential in a file
  //
  // fileName = name of the file to stock the potential
  virtual void SavePotential(char* fileName);

  // load the potential from a file
  //
  // fileName = name of the file in which the potential is stocked
  virtual void LoadPotential(char* fileName);

  // assign the potential a value at a given position 
  //
  // k = z coordinate of the considered cylinder 
  // value = value of potential
  virtual void SetPotential(int k, double value);

  // get the potential at a given position
  //
  // k = z coordinate of the considered cylinder
  // return value = the potential in the cylinder
  virtual double GetPotential(int k);

  // get the cylinder height of a given cylinder
  //
  //  k = z coordinate of the considered cylinder
  // return value = the height of the cylinder
  virtual double GetHeight(int k);

  // get the cylinder radius of a given cylinder
  //
  //  k = z coordinate of the considered cylinder
  // return value = the radius of the cylinder
  virtual double GetRadius(int k);

  // get a 1-D potential by averaging through a gaussian weight
  //
  // sigma = the variance of the gaussian function
  // m = value of angular momentum
  // return = pointer to the one dimension potential
  OneDConstantCellPotential* GaussianReductionOneDimension (double sigma, int m = 0);

  // get a 1-D potential by averaging through a weight of C(r^2 - a^2) exp(-r^2/(2s^2))
  //
  // lambda = the variance of the gaussian function
  // sigma = the variance of the gaussian function for the 1S state
  // m = value of angular momentum
  // return = pointer to the one dimension potential  
  OneDConstantCellPotential* SecondDegreeReductionOneDimension (double lambda, double sigma, int m);

};

// get the super cylinder radius
//
  
inline double ThreeDConstantCylinderPotential::GetSuperCylinderRadius()
{
  return this->SuperCylinderRadius;
}

// get the number of cylinders in Y direction
//

inline int ThreeDConstantCylinderPotential::GetNbrCylinderZ()
{
  return this->NumberZ;
}

// assign the potential a value at a given position 
//
// k = z coordinate of the considered cylinder 
// value = value of potential

inline void ThreeDConstantCylinderPotential::SetPotential(int k, double value)
{
  this->PotentialValue[k] = value;
}

// get the potential at a given position
//
// k = z coordinate of the considered cylinder
// return value = the potential in the cylinder

inline double ThreeDConstantCylinderPotential::GetPotential(int k)
{
  return this->PotentialValue[k]; 
}

// get the cylinder height of a given cylinder
//
//  k = z coordinate of the considered cylinder
// return value = the height of the cylinder

inline double ThreeDConstantCylinderPotential::GetHeight(int k)
{
  return this->CylinderHeight[k];
}

// get the cylinder radius of a given cylinder
//
//  k = z coordinate of the considered cylinder
// return value = the radius of the cylinder

inline double ThreeDConstantCylinderPotential::GetRadius(int k)
{
  return this->CylinderRadius[k]; 
}

#endif
