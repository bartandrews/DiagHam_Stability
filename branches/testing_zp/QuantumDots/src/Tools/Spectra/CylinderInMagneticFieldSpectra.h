////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//          class for periodic average spectra with Landau states             //
//                                                                            //
//                        last modification : 05/04/2004                      //
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


#ifndef CYLINDERINMAGNETICFIELDSPECTRA_H
#define CYLINDERINMAGNETICFIELDSPECTRA_H


#include "config.h"

#include "HilbertSpace/PlanarRotationSymmetryZPeriodicOneParticle.h"

class CylinderInMagneticFieldSpectra
{
 protected:

  // quantum number of kinetic momentum in Z direction
  int NumberM;
  
  // wave function basis dimension in the plane
  int NbrStateR;
  
  // wave function basis dimension in the z direction
  int NbrStateZ;
  int LowerImpulsionZ;
  
  // magnetic field in Z direction
  double Bz;

  double** RealCoefficients;
  double** ImaginaryCoefficients;

 public:

  // constructor from a Hilbert space and a file
  //
  // space = Hilbert space describing the particle
  // fileName = name of the state file
  // bz = magnetic field in Z direction
  CylinderInMagneticFieldSpectra(PlanarRotationSymmetryZPeriodicOneParticle* space, char* fileName, double bz);

  // get the wave function value of a state at a given point
  //
  // x, y, z : the position of the point
  // Real, Imaginary : references to the real and imaginary components of the wave function
  void WaveFunctionValue(double x, double y, double z, double SizeZ, double& Real, double& Imaginary);

  // get the probability density in z direction (i.e. to sum the probability in the plane)
  //
  // z = z position
  // sizeZ = sample size in Z direction
  // return = probability density in Z direction at the given point
  double ZProbabilityDensity(double z, double sizeZ);

  // get the value of impulsion operators with another wavefunction <this|p|another>
  //
  // space = Hilbert space describing the other particle
  // fileName = the file to stock the other function
  // sizeZ = size of sample in Z direction
  // impulsionX, impulsionY, impulsionZ = reference to the return values
  void GetImpulsion(PlanarRotationSymmetryZPeriodicOneParticle* space, char* fileName, double sizeZ, double &realImpulsionX, double &imaginaryImpulsionX, double &realImpulsionY, double &imaginaryImpulsionY, double &realImpulsionZ, double &imaginaryImpulsionZ);

  // get mean value of position operator
  //
  // space = Hilbert space describing the other particle
  // fileName = the file to stock the other function
  // sizeZ = size of sample in Z direction
  // positionX, positionY, positionZ = reference to the return values
  void GetMeanPosition(PlanarRotationSymmetryZPeriodicOneParticle* space, char* fileName, double sizeZ, double &realPositionX, double &imaginaryPositionX, double &realPositionY, double &imaginaryPositionY, double &realPositionZ, double &imaginaryPositionZ);

  // get the value of <phi|r²|phi>
  //
  // return = the value of <phi|r²|phi> in Angstrom unit
  double GetSquaredRadius ();
  
};

#endif
