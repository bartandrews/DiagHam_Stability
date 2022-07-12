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



#ifndef XYREFLEXIONSYMMETRICPERIODICSPECTRA_H
#define XYREFLEXIONSYMMETRICPERIODICSPECTRA_H

#include "config.h"

#include "HilbertSpace/PeriodicXYReflexionZPeriodicThreeDOneParticle.h"


class XYReflexionSymmetricPeriodicSpectra
{
 protected:
  // wave function basis dimension in the x direction
  int NbrStateX;
  int LowerImpulsionX;
  
  // wave function basis dimension in the y direction
  int NbrStateY;
  int LowerImpulsionY;
  
  // wave function basis dimension in the z direction
  int NbrStateZ;
  int LowerImpulsionZ;

  double*** RealCoefficients;
  double*** ImaginaryCoefficients;

 public:

  // constructor from a Hilbert space and a file
  //
  // space = Hilbert space describing the particle
  // fileName = name of the state file
  XYReflexionSymmetricPeriodicSpectra(PeriodicXYReflexionZPeriodicThreeDOneParticle* space, char* fileName);

  // get the wave function value of a state at a given point
  //
  // x, y, z : the position of the point
  // SizeX, SizeY, SizeZ : the 3D-sizes of the sample
  // Real, Imaginary : references to the real and imaginary components of the wave function
  void WaveFunctionValue(double x, double SizeX, double y, double SizeY, double z, double SizeZ, double& Real, double& Imaginary);

  // get the integrated probability density over z direction
  //
  // x, y = coordinates of the point
  // sizeX, sizeY = dimensions of the system in x and y directions
  // return = value of the probability density  
  double PlanarProbabilityDensity(double x, double sizeX, double y, double sizeY);

  // get the value of impulsion operators with another wave function <this|p|another>
  //
  // space = Hilbert space describing the other particle
  // fileName = the file to stock the other function
  // sizeX, sizeY, sizeZ = size of sample in X, Y and Z directions
  // impulsionX, impulsionY, impulsionZ = reference to the return values
  void GetImpulsion(PeriodicXYReflexionZPeriodicThreeDOneParticle* space, char* fileName, double sizeX, double sizeY, double sizeZ, double &realImpulsionX, double &imaginaryImpulsionX, double &realImpulsionY, double &imaginaryImpulsionY, double &realImpulsionZ, double &imaginaryImpulsionZ);

  // get the value of impulsion  with a continuum state, described by Z plane waves
  //
  // impulsionX, impulsionY = value of the planer impulsion
  // nbrStateZ, lowerImpulsionZ = Hilbert space characteristics 
  // fileName =  the file to stock the other function
  // sizeZ = size of the sample in Z direction
  // impulsionZ = reference to the return values  
  void GetImpulsionWithContinuum (int impulsionX, int impulsionY, int nbrStateZ, int lowerImpulsionZ, char* fileName, double sizeZ, double &realImpulsionZ, double &imaginaryImpulsionZ);

  // get the overlap of derived functions
  //
  // space = Hilbert space describing the other particle
  // fileName = the file to stock the other function
  // sizeX, sizeY, sizeZ = size of sample in X, Y and Z directions
  // overlapX, overlapY = reference to the return values
  void GetDerivedOverlap (PeriodicXYReflexionZPeriodicThreeDOneParticle* space, char* fileName, double sizeX, double sizeY, double sizeZ, double &realOverlap, double &imaginaryOverlap, double &realOverlapX, double &imaginaryOverlapX, double &realOverlapY, double &imaginaryOverlapY);

  // get the probability of finding the particle in a cube
  //
  // minX, maxX, minY, maxY, minZ, maxZ = bounds of the cube in unit of proportion compared to the whole length
  // return = value of the probability  
  double GetCubeProbability (double minX, double maxX, double minY, double maxY, double minZ, double maxZ);

};

#endif
