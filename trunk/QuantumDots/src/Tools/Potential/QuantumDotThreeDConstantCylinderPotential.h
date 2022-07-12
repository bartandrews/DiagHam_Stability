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


#ifndef QUANTUMDOTTHREEDCONSTANTCYLINDERPOTENTIAL_H
#define QUANTUMDOTTHREEDCONSTANTCYLINDERPOTENTIAL_H


#include "config.h"
#include "Tools/Potential/ThreeDConstantCylinderPotential.h"


class QuantumDotThreeDConstantCylinderPotential : public ThreeDConstantCylinderPotential
{

 protected:

  // height of well barrier below the dot
  double UnderBarrier;

  // height of cylinder below the wetting layer
  double BelowHeight;

  // width of the wetting layer
  double WettingWidth;

  // number of cylinders in the dot
  int NbrCylinderDot;

  // height of the dot
  double DotHeight;

  // radius of the dot base
  double BaseRadius;

  // radius of the dot top
  double TopRadius;

  // height of cylinder above the dot
  double AboveHeight;

 public:

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
  QuantumDotThreeDConstantCylinderPotential(double belowHeight, double wettingWidth, int nbrCylinderDot, double dotHeight, double baseRadius, double topRadius, double aboveHeight, double underBarrier = 0.0, double superCylinderRadius = -1.0);
  
  // destructor
  //
  ~QuantumDotThreeDConstantCylinderPotential();

  // construct the potential profile from potential values
  //
  // dotPotential = dot potential, in comparision with the bulk
  // wellPotential = barrier potential, in comparison with the bulk
  void ConstructPotential(double dotPotential, double wellPotential = 0.0);

  // add two barriers of the same potential just below and above the dot + WL (for Vienna experimenters)
  //
  // belowBarrier = width of barrier layer just below the WL (in Angstrom unit)
  // aboveBarrier = width of barrier layer just above  the dot (in Angstrom unit)
  // potential = potential to add in these barriers (in eV unit)  
  void AddBarrierPotential (double belowBarrier, double aboveBarrier, double potential);  

  // save the whole diagram presentation in a bitmap file
  //
  // fileName = name of the file to stock the diagram presentation
  virtual void SaveBmpPicture(char* fileName);  
};


#endif
