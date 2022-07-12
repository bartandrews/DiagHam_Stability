////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//       class for periodic average spectra with XY reflexion symmetry        //
//                                                                            //
//                        last modification : 04/04/2004                      //
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



#ifndef HARDBOXSPECTRA_H
#define HARDBOXSPECTRA_H

#include "config.h"

#include "HilbertSpace/ThreeDOneParticle.h"


class HardBoxSpectra
{
 protected:
  // wave function basis dimension in the x direction
  int NbrStateX;
  
  // wave function basis dimension in the y direction
  int NbrStateY;
  
  // wave function basis dimension in the z direction
  int NbrStateZ;

  double*** Coefficients;

 public:

  // constructor from a Hilbert space and a file
  //
  // space = Hilbert space describing the particle
  // fileName = name of the state file
  HardBoxSpectra(ThreeDOneParticle* space, char* fileName);

  // get the overlap of derived functions
  //
  // space = Hilbert space describing the other particle
  // fileName = the file to stock the other function
  // sizeX, sizeY, sizeZ = size of sample in X, Y and Z directions
  // overlapX, overlapY = reference to the return values
  void GetDerivedOverlap (ThreeDOneParticle* space, char* fileName, double sizeX, double sizeY, double sizeZ, double &Overlap, double &OverlapX, double &OverlapY);

};

#endif
