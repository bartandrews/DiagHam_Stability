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


#include "Tools/Spectra/CylinderInMagneticFieldSpectra.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;

#define LENGTH_FACTOR 256.8

// constructor from a Hilbert space and a file
//
// space = Hilbert space describing the particle
// fileName = name of the state file
// bz = magnetic field in Z direction
CylinderInMagneticFieldSpectra::CylinderInMagneticFieldSpectra(PlanarRotationSymmetryZPeriodicOneParticle* space, char* fileName, double bz)
{
  this->NumberM = space->GetLz();
  this->NbrStateR = space->GetNbrStateR();
  this->NbrStateZ = space->GetNbrStateZ();
  this->LowerImpulsionZ = space->GetLowerImpulsionZ();
  this->Bz = bz;

  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Error in open the file: " << fileName << endl;
      exit(0);
    }

  this->RealCoefficients = new double* [this->NbrStateR];
  this->ImaginaryCoefficients = new double* [this->NbrStateR];
  for (int i = 0; i < this->NbrStateR; ++i)
    {
      this->RealCoefficients[i] = new double [this->NbrStateZ];
      this->ImaginaryCoefficients[i] = new double [this->NbrStateZ];         
      for (int k = 0; k < this->NbrStateZ; ++k)		
	File >> this->RealCoefficients[i][k] >> this->ImaginaryCoefficients[i][k];	    	
    } 
  File.close();
}

// get the wave function value of a state at a given point
//
// x, y, z : the position of the point
// SizeZ : the 3D-sizes of the sample
// Real, Imaginary : references to the real and imaginary components of the wave function

void CylinderInMagneticFieldSpectra::WaveFunctionValue(double x, double y, double z, double SizeZ, double& Real, double& Imaginary)
{
}

// get the probability density in z direction (i.e. to sum the probability in the plane)
//
// z = z position
// sizeZ = sample size in Z direction
// return = probability density in Z direction at the given point

double CylinderInMagneticFieldSpectra::ZProbabilityDensity(double z, double sizeZ)
{
  int LengthZ = (this->NbrStateZ - 1) * 2 + 1;
  double* real = new double [LengthZ];
  double* imaginary = new double [LengthZ];
  int OriginZ = this->NbrStateZ - 1;
  double tmp = 2.0 * M_PI * z / sizeZ;
  for (int delta = 0; delta < LengthZ; ++delta)
    {
      real[delta] = cos (tmp * double (delta - OriginZ));
      imaginary[delta] = sin (tmp * double (delta - OriginZ));
    }
  double Tmp = 0.0; int IndexZ = 0;
  double* TmpRealCoefficient; double* TmpImaginaryCoefficient; 
  double TmpRe1 = 0.0, TmpIm1 = 0.0;
  for (int n = 0; n < this->NbrStateR; ++n)
    {
      TmpRealCoefficient = this->RealCoefficients[n]; TmpImaginaryCoefficient = this->ImaginaryCoefficients[n];
      for (int p1 = 0; p1 < this->NbrStateZ; ++p1)
	{
	  IndexZ = -p1 + OriginZ;
	  
	  for (int p2 = 0; p2 < this->NbrStateZ; ++p2) 
	    {
	      TmpRe1 = TmpRealCoefficient[p1] * TmpRealCoefficient[p2] + TmpImaginaryCoefficient[p1] * TmpImaginaryCoefficient[p2];
	      TmpIm1 = TmpRealCoefficient[p1] * TmpImaginaryCoefficient[p2] - TmpImaginaryCoefficient[p1] * TmpRealCoefficient[p2]; 
	      Tmp += (TmpRe1 * real[IndexZ] - TmpIm1 * imaginary[IndexZ]);
	      ++IndexZ;
	    }
	}
    }
  Tmp /= sizeZ; 
  return Tmp;
}

// get the value of impulsion operators with another wavefunction <this|p|another>
//
// space = Hilbert space describing the other particle
// fileName = the file to stock the other function
// sizeZ = size of sample in Z direction
// impulsionX, impulsionY, impulsionZ = reference to the return values

void CylinderInMagneticFieldSpectra::GetImpulsion(PlanarRotationSymmetryZPeriodicOneParticle* space, char* fileName, double sizeZ, double &realImpulsionX, double &imaginaryImpulsionX, double &realImpulsionY, double &imaginaryImpulsionY, double &realImpulsionZ, double &imaginaryImpulsionZ)
{
  double OrbitRadius = LENGTH_FACTOR / sqrt(this->Bz);
  int numberM = space->GetLz();
  int nbrStateR = space->GetNbrStateR();
  int nbrStateZ = space->GetNbrStateZ();
  int lowerImpulsionZ = space->GetLowerImpulsionZ();

  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Error in open the file: " << fileName << endl;
      exit(0);
    }
  double** realCoefficients = new double* [nbrStateR];
  double** imaginaryCoefficients = new double* [nbrStateR];
  for (int i = 0; i < nbrStateR; ++i)
    {
      realCoefficients[i] = new double [nbrStateZ];
      imaginaryCoefficients[i] = new double [nbrStateZ];         
      for (int k = 0; k < nbrStateZ; ++k)		
	File >> realCoefficients[i][k] >> imaginaryCoefficients[i][k]; 	
    }
  File.close();  

  int MaxLowerZ = this->LowerImpulsionZ;
  if (this->LowerImpulsionZ < lowerImpulsionZ)
    MaxLowerZ = lowerImpulsionZ;

  int MinUpperR = this->NbrStateR, MinUpperZ = this->LowerImpulsionZ + this->NbrStateZ;
  if (MinUpperR > (nbrStateR))
    MinUpperR = nbrStateR;
  if (MinUpperZ > (lowerImpulsionZ + nbrStateZ))
    MinUpperZ = lowerImpulsionZ + nbrStateZ;

  realImpulsionX = 0.0;
  realImpulsionY = 0.0;
  realImpulsionZ = 0.0;
  imaginaryImpulsionX = 0.0;
  imaginaryImpulsionY = 0.0;
  imaginaryImpulsionZ = 0.0;
  double TmpRe = 0.0, TmpIm = 0.0;
  double* Re1; double* Im1; double* Re2; double* Im2; double* Re3; double* Im3;
  double TmpRe1 = 0.0, TmpIm1 = 0.0;
  // x and y directions
  if (this->NumberM == numberM)
    {
      realImpulsionX = 0.0;
      imaginaryImpulsionX = 0.0;
      realImpulsionY = 0.0;
      imaginaryImpulsionY = 0.0;
    }
  else
    {
      TmpRe = 0.0; TmpIm = 0.0;
      for (int n = 0; n < (MinUpperR - 1); ++n)
	{	      
	  Re1 = this->RealCoefficients[n];
	  Im1 = this->ImaginaryCoefficients[n];
	  Re2 = realCoefficients[n];
	  Im2 = imaginaryCoefficients[n];
	  Re3 = this->RealCoefficients[n + 1];
	  Im3 = this->ImaginaryCoefficients[n + 1];
	  TmpRe1 = 0.0; TmpIm1 = 0.0;
	  for (int p = MaxLowerZ; p < MinUpperZ; ++p)
	    {
	      TmpRe1 += (Re1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]);
	      TmpIm1 += (-Re1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]);
	      TmpRe1 += (Re3[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im3[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]);
	      TmpIm1 += (-Re3[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] + Im3[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]);	      
	    }
	  TmpRe += (TmpRe1 * sqrt(double(n + 1)));
	  TmpIm += (TmpIm1 * sqrt(double(n + 1)));
	}
      Re1 = this->RealCoefficients[MinUpperR - 1];
      Im1 = this->ImaginaryCoefficients[MinUpperR - 1];
      Re2 = realCoefficients[MinUpperR - 1];
      Im2 = imaginaryCoefficients[MinUpperR - 1];
      TmpRe1 = 0.0; TmpIm1 = 0.0;
      for (int p = MaxLowerZ; p < MinUpperZ; ++p)
	{
	  TmpRe1 += (Re1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]);
	  TmpIm1 += (-Re1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]);
		    
	}
      TmpRe += (TmpRe1 * sqrt(double(MinUpperR)));
      TmpIm += (TmpIm1 * sqrt(double(MinUpperR)));
      TmpRe /= (sqrt(8.0) * OrbitRadius); TmpIm /= (sqrt(8.0) * OrbitRadius); 
      realImpulsionX = -TmpIm; realImpulsionY = -TmpIm;
      imaginaryImpulsionX = TmpRe; imaginaryImpulsionY = TmpRe; 
    }
  // z direction
  if (this->NumberM != numberM)
    {
      realImpulsionZ = 0.0;
      imaginaryImpulsionZ = 0.0;
    }
  else
    {          
      for (int n = 0; n < MinUpperR; ++n)
	{	      
	  Re1 = this->RealCoefficients[n];
	  Im1 = this->ImaginaryCoefficients[n];
	  Re2 = realCoefficients[n];
	  Im2 = imaginaryCoefficients[n];
	  for (int p = MaxLowerZ; p < MinUpperZ; ++p)
	    {
	      realImpulsionZ += ((Re1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]) * double(p));
	      imaginaryImpulsionZ += ((Re1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] - Im1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]) * double(p));
	    }
	}
      realImpulsionZ = realImpulsionZ * 2.0 * M_PI / sizeZ; imaginaryImpulsionZ = imaginaryImpulsionZ * 2.0 * M_PI / sizeZ;	    
    }  
}

// get mean value of position operator
//
// space = Hilbert space describing the other particle
// fileName = the file to stock the other function
// sizeZ = size of sample in Z direction
// positionX, positionY, positionZ = reference to the return values

void CylinderInMagneticFieldSpectra::GetMeanPosition(PlanarRotationSymmetryZPeriodicOneParticle* space, char* fileName, double sizeZ, double &realPositionX, double &imaginaryPositionX, double &realPositionY, double &imaginaryPositionY, double &realPositionZ, double &imaginaryPositionZ)
{
  double OrbitRadius = LENGTH_FACTOR / sqrt(this->Bz);
  int numberM = space->GetLz();
  int nbrStateR = space->GetNbrStateR();
  int nbrStateZ = space->GetNbrStateZ();
  int lowerImpulsionZ = space->GetLowerImpulsionZ();

  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Error in open the file: " << fileName << endl;
      exit(0);
    }
  double** realCoefficients = new double* [nbrStateR];
  double** imaginaryCoefficients = new double* [nbrStateR];
  for (int i = 0; i < nbrStateR; ++i)
    {
      realCoefficients[i] = new double [nbrStateZ];
      imaginaryCoefficients[i] = new double [nbrStateZ];         
      for (int k = 0; k < nbrStateZ; ++k)		
	File >> realCoefficients[i][k] >> imaginaryCoefficients[i][k]; 	
    }
  File.close();  

  int MaxLowerZ = this->LowerImpulsionZ;
  if (this->LowerImpulsionZ < lowerImpulsionZ)
    MaxLowerZ = lowerImpulsionZ;

  int MinUpperR = this->NbrStateR, MinUpperZ = this->LowerImpulsionZ + this->NbrStateZ;
  if (MinUpperR > (nbrStateR))
    MinUpperR = nbrStateR;
  if (MinUpperZ > (lowerImpulsionZ + nbrStateZ))
    MinUpperZ = lowerImpulsionZ + nbrStateZ;

  realPositionX = 0.0;
  realPositionY = 0.0;
  realPositionZ = 0.0;
  imaginaryPositionX = 0.0;
  imaginaryPositionY = 0.0;
  imaginaryPositionZ = 0.0;
  double TmpRe = 0.0, TmpIm = 0.0;
  double* Re1; double* Im1; double* Re2; double* Im2; double* Re3; double* Im3;
  double TmpRe1 = 0.0, TmpIm1 = 0.0;

   // x and y directions
  if (this->NumberM == numberM)
    {
      realPositionX = 0.0;
      imaginaryPositionX = 0.0;
      realPositionY = 0.0;
      imaginaryPositionY = 0.0;
    } 
  else
    {
      TmpRe = 0.0; TmpIm = 0.0;
      for (int n = 0; n < (MinUpperR - 1); ++n)
	{	      
	  Re1 = this->RealCoefficients[n];
	  Im1 = this->ImaginaryCoefficients[n];
	  Re2 = realCoefficients[n];
	  Im2 = imaginaryCoefficients[n];
	  Re3 = this->RealCoefficients[n + 1];
	  Im3 = this->ImaginaryCoefficients[n + 1];
	  TmpRe1 = 0.0; TmpIm1 = 0.0;
	  for (int p = MaxLowerZ; p < MinUpperZ; ++p)
	    {
	      TmpRe1 += (Re1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]);
	      TmpIm1 += (-Re1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]);
	      TmpRe1 -= (Re3[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im3[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]);
	      TmpIm1 -= (-Re3[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] + Im3[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]);	      
	    }
	  TmpRe += (TmpRe1 * sqrt(double(n + 1)));
	  TmpIm += (TmpIm1 * sqrt(double(n + 1)));
	}
      Re1 = this->RealCoefficients[MinUpperR - 1];
      Im1 = this->ImaginaryCoefficients[MinUpperR - 1];
      Re2 = realCoefficients[MinUpperR - 1];
      Im2 = imaginaryCoefficients[MinUpperR - 1];
      TmpRe1 = 0.0; TmpIm1 = 0.0;
      for (int p = MaxLowerZ; p < MinUpperZ; ++p)
	{
	  TmpRe1 += (Re1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]);
	  TmpIm1 += (-Re1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]);
		    
	}
      TmpRe += (TmpRe1 * sqrt(double(MinUpperR)));
      TmpIm += (TmpIm1 * sqrt(double(MinUpperR)));
      TmpRe *= (OrbitRadius / sqrt(2.0)); TmpIm *= (OrbitRadius / sqrt(2.0)); 
      realPositionX = -TmpIm; realPositionY = -TmpIm;
      imaginaryPositionX = TmpRe; imaginaryPositionY = TmpRe;      
    }
}

// get the value of <phi|r²|phi>
//
// return = the value of <phi|r²|phi> in Angstrom unit

double CylinderInMagneticFieldSpectra::GetSquaredRadius ()
{
  double OrbitRadius = LENGTH_FACTOR / sqrt(this->Bz);
  if (this->NumberM == 0)
    {
      double Res = 0.0;
      for (int p = 0; p < this->NbrStateZ; ++p)
	{
	  Res += (this->RealCoefficients[0][p] * this->RealCoefficients[0][p] + this->ImaginaryCoefficients[0][p] * this->ImaginaryCoefficients[0][p]);
	  Res -= (this->RealCoefficients[0][p] * this->RealCoefficients[1][p] + this->ImaginaryCoefficients[0][p] * this->ImaginaryCoefficients[1][p]);
	}
      for (int n = 1; n < (this->NbrStateR - 1); ++n)
	{
	  for (int p = 0; p < this->NbrStateZ; ++p)
	    { 
	      Res += (2 * n + 1) * (this->RealCoefficients[n][p] * this->RealCoefficients[n][p] + this->ImaginaryCoefficients[n][p] * this->ImaginaryCoefficients[n][p]);
	      Res -= n * (this->RealCoefficients[n][p] * this->RealCoefficients[n - 1][p] + this->ImaginaryCoefficients[n][p] * this->ImaginaryCoefficients[n - 1][p]);
	      Res -= (n + 1) * (this->RealCoefficients[n][p] * this->RealCoefficients[n + 1][p] + this->ImaginaryCoefficients[n][p] * this->ImaginaryCoefficients[n + 1][p]);    
	    }
	}
      int N = this->NbrStateR - 1;
      for (int p = 0; p < this->NbrStateZ; ++p)
	{
	  Res += (2 * N + 1) * (this->RealCoefficients[N][p] * this->RealCoefficients[N][p] + this->ImaginaryCoefficients[N][p] * this->ImaginaryCoefficients[N][p]);
	  Res -= N * (this->RealCoefficients[N][p] * this->RealCoefficients[N - 1][p] + this->ImaginaryCoefficients[N][p] * this->ImaginaryCoefficients[N - 1][p]);	  
	} 
      return Res * 2.0 * OrbitRadius * OrbitRadius;
    }

  return 0;
}
