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


#include "Tools/Spectra/XYReflexionSymmetricPeriodicSpectra.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;


// constructor from a Hilbert space and a file
//
// space = Hilbert space describing the particle
// fileName = name of the state file

XYReflexionSymmetricPeriodicSpectra::XYReflexionSymmetricPeriodicSpectra(PeriodicXYReflexionZPeriodicThreeDOneParticle* space, char* fileName)
{
  this->NbrStateX = space->GetNbrStateX();
  this->NbrStateY = space->GetNbrStateY();
  this->NbrStateZ = space->GetNbrStateZ();
  this->LowerImpulsionX = space->GetLowerImpulsionX();
  this->LowerImpulsionY = space->GetLowerImpulsionY();
  this->LowerImpulsionZ = space->GetLowerImpulsionZ();

  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Error in open the file: " << fileName << endl;
      exit(0);
    }

  this->RealCoefficients = new double** [this->NbrStateX];
  this->ImaginaryCoefficients = new double** [this->NbrStateX];
  for (int i = 0; i < this->NbrStateX; ++i)
    {
      this->RealCoefficients[i] = new double* [this->NbrStateY];
      this->ImaginaryCoefficients[i] = new double* [this->NbrStateY];   
      for (int j = 0; j < this->NbrStateY; ++j)
	{
	  this->RealCoefficients[i][j] = new double [this->NbrStateZ];
	  this->ImaginaryCoefficients[i][j] = new double [this->NbrStateZ];
	  for (int k = 0; k < this->NbrStateZ; ++k)	    
	    File >> this->RealCoefficients[i][j][k] >> this->ImaginaryCoefficients[i][j][k];	    
	}
    } 
  File.close();
}

// get the wave function value of a state at a given point
//
// x, y, z : the position of the point
// SizeX, SizeY, SizeZ : the 3D-sizes of the sample
// Real, Imaginary : references to the real and imaginary components of the wave function

void XYReflexionSymmetricPeriodicSpectra::WaveFunctionValue(double x, double SizeX, double y, double SizeY, double z, double SizeZ, double& Real, double& Imaginary)
{
  // pairX pairY
  if ((this->LowerImpulsionX == 0) && (this->LowerImpulsionY == 0))
    {
      double TmpX = 0.0; double TmpY = 0.0; double TmpZ = 0.0;
      double TmpRe = 0.0; double TmpIm = 0.0;
      double* TmpRealCoefficients; double* TmpImaginaryCoefficients;
      double TmpRe2 = 0.0; double TmpIm2 = 0.0;
      double FactorX = 2.0 * M_PI * x / SizeX; double FactorY = 2.0 * M_PI * y / SizeY; double FactorZ = 2.0 * M_PI * z / SizeZ;
      double TmpReZ, TmpImZ, TmpReYZ, TmpImYZ;
      // m = 0
      TmpX = sqrt(0.5);
      TmpReYZ = 0.0; TmpImYZ = 0.0;
      // m = 0; n = 0
      TmpY = sqrt(0.5);
      TmpRealCoefficients = this->RealCoefficients[0][0];
      TmpImaginaryCoefficients = this->ImaginaryCoefficients[0][0];
      TmpReZ = 0.0; TmpImZ = 0.0;
      for (int p = 0; p < this->NbrStateZ; ++p)
	{
	  TmpZ = double(p + this->LowerImpulsionZ) * FactorZ;
	  TmpRe2 = cos(TmpZ); TmpIm2 = sin(TmpZ);
	  TmpReZ += (TmpRealCoefficients[p] * TmpRe2 - TmpImaginaryCoefficients[p] * TmpIm2);
	  TmpImZ += (TmpRealCoefficients[p] * TmpIm2 + TmpImaginaryCoefficients[p] * TmpRe2);
	}
      TmpReYZ += TmpY * TmpReZ; 
      TmpImYZ += TmpY * TmpImZ;
      TmpRe += TmpX * TmpReYZ;
      TmpIm += TmpX * TmpImYZ;

      // m = 0; n != 0   
      TmpReYZ = 0.0; TmpImYZ = 0.0;   
      for (int n = 1; n < this->NbrStateY; ++n)
	{
	  TmpY = cos(double(n) * FactorY);
	  TmpRealCoefficients = this->RealCoefficients[0][n];
	  TmpImaginaryCoefficients = this->ImaginaryCoefficients[0][n];
	  TmpReZ = 0.0; TmpImZ = 0.0;
	  for (int p = 0; p < this->NbrStateZ; ++p)
	    {
	      TmpZ = double(p + this->LowerImpulsionZ) * FactorZ;
	      TmpRe2 = cos(TmpZ); TmpIm2 = sin(TmpZ);
	      TmpReZ += (TmpRealCoefficients[p] * TmpRe2 - TmpImaginaryCoefficients[p] * TmpIm2);
	      TmpImZ += (TmpRealCoefficients[p] * TmpIm2 + TmpImaginaryCoefficients[p] * TmpRe2);
	    }
	  TmpReYZ += TmpY * TmpReZ; 
	  TmpImYZ += TmpY * TmpImZ;
	}
      TmpRe += TmpX * TmpReYZ;
      TmpIm += TmpX * TmpImYZ;

      // m != 0
      for (int m = 1; m < this->NbrStateX; ++m)
	{
	  TmpX = cos(double(m) * FactorX);
	  TmpReYZ = 0.0; TmpImYZ = 0.0;
	  // m != 0; n = 0
	  TmpY = sqrt(0.5);
	  TmpRealCoefficients = this->RealCoefficients[m][0];
	  TmpImaginaryCoefficients = this->ImaginaryCoefficients[m][0];
	  TmpReZ = 0.0; TmpImZ = 0.0;
	  for (int p = 0; p < this->NbrStateZ; ++p)
	    {
	      TmpZ = double(p + this->LowerImpulsionZ) * FactorZ;
	      TmpRe2 = cos(TmpZ); TmpIm2 = sin(TmpZ);
	      TmpReZ += (TmpRealCoefficients[p] * TmpRe2 - TmpImaginaryCoefficients[p] * TmpIm2);
	      TmpImZ += (TmpRealCoefficients[p] * TmpIm2 + TmpImaginaryCoefficients[p] * TmpRe2);
	    }
	  TmpReYZ += TmpY * TmpReZ; 
	  TmpImYZ += TmpY * TmpImZ;
	  TmpRe += TmpX * TmpReYZ;
	  TmpIm += TmpX * TmpImYZ;

	  // m != 0; n!= 0
	  TmpReYZ = 0.0; TmpImYZ = 0.0;
	  for (int n = 1; n < this->NbrStateY; ++n)
	    {
	      TmpY = cos(double(n) * FactorY);
	      TmpRealCoefficients = this->RealCoefficients[m][n];
	      TmpImaginaryCoefficients = this->ImaginaryCoefficients[m][n];
	      TmpReZ = 0.0; TmpImZ = 0.0;
	      for (int p = 0; p < this->NbrStateZ; ++p)
		{
		  TmpZ = double(p + this->LowerImpulsionZ) * FactorZ;
		  TmpRe2 = cos(TmpZ); TmpIm2 = sin(TmpZ);
		  TmpReZ += (TmpRealCoefficients[p] * TmpRe2 - TmpImaginaryCoefficients[p] * TmpIm2);
		  TmpImZ += (TmpRealCoefficients[p] * TmpIm2 + TmpImaginaryCoefficients[p] * TmpRe2);
		}
	      TmpReYZ += TmpY * TmpReZ; 
	      TmpImYZ += TmpY * TmpImZ;
	    }
	  TmpRe += TmpX * TmpReYZ;
	  TmpIm += TmpX * TmpImYZ;
	}
      Real = TmpRe * sqrt(4.0 / (SizeX * SizeY * SizeZ));
      Imaginary = TmpIm * sqrt(4.0 / (SizeX * SizeY * SizeZ));   
    }  

  // pairX impairY
  if ((this->LowerImpulsionX == 0) && (this->LowerImpulsionY == 1))
    {
      double TmpX = 0.0; double TmpY = 0.0; double TmpZ = 0.0;
      double TmpRe = 0.0; double TmpIm = 0.0;
      double* TmpRealCoefficients; double* TmpImaginaryCoefficients;
      double TmpRe2 = 0.0; double TmpIm2 = 0.0;
      double FactorX = 2.0 * M_PI * x / SizeX; double FactorY = 2.0 * M_PI * y / SizeY; double FactorZ = 2.0 * M_PI * z / SizeZ;
      double TmpReZ, TmpImZ, TmpReYZ, TmpImYZ;
      // m = 0
      TmpX = sqrt(0.5); 
      TmpReYZ = 0.0; TmpImYZ = 0.0;   
      for (int n = 0; n < this->NbrStateY; ++n)
	{
	  TmpY = sin(double(n + 1) * FactorY);
	  TmpRealCoefficients = this->RealCoefficients[0][n];
	  TmpImaginaryCoefficients = this->ImaginaryCoefficients[0][n];
	  TmpReZ = 0.0; TmpImZ = 0.0;
	  for (int p = 0; p < this->NbrStateZ; ++p)
	    {
	      TmpZ = double(p + this->LowerImpulsionZ) * FactorZ;
	      TmpRe2 = cos(TmpZ); TmpIm2 = sin(TmpZ);
	      TmpReZ += (TmpRealCoefficients[p] * TmpRe2 - TmpImaginaryCoefficients[p] * TmpIm2);
	      TmpImZ += (TmpRealCoefficients[p] * TmpIm2 + TmpImaginaryCoefficients[p] * TmpRe2);
	    }
	  TmpReYZ += TmpY * TmpReZ; 
	  TmpImYZ += TmpY * TmpImZ;
	}
      TmpRe += TmpX * TmpReYZ;
      TmpIm += TmpX * TmpImYZ;

      // m != 0
      for (int m = 1; m < this->NbrStateX; ++m)
	{
	  TmpX = cos(double(m) * FactorX);
	  TmpReYZ = 0.0; TmpImYZ = 0.0;
	  for (int n = 0; n < this->NbrStateY; ++n)
	    {
	      TmpY = sin(double(n + 1) * FactorY);
	      TmpRealCoefficients = this->RealCoefficients[m][n];
	      TmpImaginaryCoefficients = this->ImaginaryCoefficients[m][n];
	      TmpReZ = 0.0; TmpImZ = 0.0;
	      for (int p = 0; p < this->NbrStateZ; ++p)
		{
		  TmpZ = double(p + this->LowerImpulsionZ) * FactorZ;
		  TmpRe2 = cos(TmpZ); TmpIm2 = sin(TmpZ);
		  TmpReZ += (TmpRealCoefficients[p] * TmpRe2 - TmpImaginaryCoefficients[p] * TmpIm2);
		  TmpImZ += (TmpRealCoefficients[p] * TmpIm2 + TmpImaginaryCoefficients[p] * TmpRe2);
		}
	      TmpReYZ += TmpY * TmpReZ; 
	      TmpImYZ += TmpY * TmpImZ;
	    }
	  TmpRe += TmpX * TmpReYZ;
	  TmpIm += TmpX * TmpImYZ;
	}
      Real = TmpRe * sqrt(4.0 / (SizeX * SizeY * SizeZ));
      Imaginary = TmpIm * sqrt(4.0 / (SizeX * SizeY * SizeZ));   
    }  

  // impairX pairY
  if ((this->LowerImpulsionX == 1) && (this->LowerImpulsionY == 0))
    {
      double TmpX = 0.0; double TmpY = 0.0; double TmpZ = 0.0;
      double TmpRe = 0.0; double TmpIm = 0.0;
      double* TmpRealCoefficients; double* TmpImaginaryCoefficients;
      double TmpRe2 = 0.0; double TmpIm2 = 0.0;
      double FactorX = 2.0 * M_PI * x / SizeX; double FactorY = 2.0 * M_PI * y / SizeY; double FactorZ = 2.0 * M_PI * z / SizeZ;
      double TmpReZ, TmpImZ, TmpReYZ, TmpImYZ;

      for (int m = 0; m < this->NbrStateX; ++m)
	{
	  TmpX = sin(double(m + 1) * FactorX);
	  TmpReYZ = 0.0; TmpImYZ = 0.0;
	  // n = 0
	  TmpY = sqrt(0.5);
	  TmpRealCoefficients = this->RealCoefficients[m][0];
	  TmpImaginaryCoefficients = this->ImaginaryCoefficients[m][0];
	  TmpReZ = 0.0; TmpImZ = 0.0;
	  for (int p = 0; p < this->NbrStateZ; ++p)
	    {
	      TmpZ = double(p + this->LowerImpulsionZ) * FactorZ;
	      TmpRe2 = cos(TmpZ); TmpIm2 = sin(TmpZ);
	      TmpReZ += (TmpRealCoefficients[p] * TmpRe2 - TmpImaginaryCoefficients[p] * TmpIm2);
	      TmpImZ += (TmpRealCoefficients[p] * TmpIm2 + TmpImaginaryCoefficients[p] * TmpRe2);
	    }
	  TmpReYZ += TmpY * TmpReZ; 
	  TmpImYZ += TmpY * TmpImZ;
	  TmpRe += TmpX * TmpReYZ;
	  TmpIm += TmpX * TmpImYZ;

	  // n!= 0
	  TmpReYZ = 0.0; TmpImYZ = 0.0;
	  for (int n = 1; n < this->NbrStateY; ++n)
	    {
	      TmpY = cos(double(n) * FactorY);
	      TmpRealCoefficients = this->RealCoefficients[m][n];
	      TmpImaginaryCoefficients = this->ImaginaryCoefficients[m][n];
	      TmpReZ = 0.0; TmpImZ = 0.0;
	      for (int p = 0; p < this->NbrStateZ; ++p)
		{
		  TmpZ = double(p + this->LowerImpulsionZ) * FactorZ;
		  TmpRe2 = cos(TmpZ); TmpIm2 = sin(TmpZ);
		  TmpReZ += (TmpRealCoefficients[p] * TmpRe2 - TmpImaginaryCoefficients[p] * TmpIm2);
		  TmpImZ += (TmpRealCoefficients[p] * TmpIm2 + TmpImaginaryCoefficients[p] * TmpRe2);
		}
	      TmpReYZ += TmpY * TmpReZ; 
	      TmpImYZ += TmpY * TmpImZ;
	    }
	  TmpRe += TmpX * TmpReYZ;
	  TmpIm += TmpX * TmpImYZ;
	}
      Real = TmpRe * sqrt(4.0 / (SizeX * SizeY * SizeZ));
      Imaginary = TmpIm * sqrt(4.0 / (SizeX * SizeY * SizeZ));   
    } 

  // impairX impairY
  if ((this->LowerImpulsionX == 1) && (this->LowerImpulsionY == 1))
    {
      double TmpX = 0.0; double TmpY = 0.0; double TmpZ = 0.0;
      double TmpRe = 0.0; double TmpIm = 0.0;
      double* TmpRealCoefficients; double* TmpImaginaryCoefficients;
      double TmpRe2 = 0.0; double TmpIm2 = 0.0;
      double FactorX = 2.0 * M_PI * x / SizeX; double FactorY = 2.0 * M_PI * y / SizeY; double FactorZ = 2.0 * M_PI * z / SizeZ;
      double TmpReZ, TmpImZ, TmpReYZ, TmpImYZ;

      for (int m = 0; m < this->NbrStateX; ++m)
	{
	  TmpX = sin(double(m + 1) * FactorX);
	  TmpReYZ = 0.0; TmpImYZ = 0.0;
	  for (int n = 0; n < this->NbrStateY; ++n)
	    {
	      TmpY = sin(double(n + 1) * FactorY);
	      TmpRealCoefficients = this->RealCoefficients[m][n];
	      TmpImaginaryCoefficients = this->ImaginaryCoefficients[m][n];
	      TmpReZ = 0.0; TmpImZ = 0.0;
	      for (int p = 0; p < this->NbrStateZ; ++p)
		{
		  TmpZ = double(p + this->LowerImpulsionZ) * FactorZ;
		  TmpRe2 = cos(TmpZ); TmpIm2 = sin(TmpZ);
		  TmpReZ += (TmpRealCoefficients[p] * TmpRe2 - TmpImaginaryCoefficients[p] * TmpIm2);
		  TmpImZ += (TmpRealCoefficients[p] * TmpIm2 + TmpImaginaryCoefficients[p] * TmpRe2);
		}
	      TmpReYZ += TmpY * TmpReZ; 
	      TmpImYZ += TmpY * TmpImZ;
	    }
	  TmpRe += TmpX * TmpReYZ;
	  TmpIm += TmpX * TmpImYZ;
	}
      Real = TmpRe * sqrt(4.0 / (SizeX * SizeY * SizeZ));
      Imaginary = TmpIm * sqrt(4.0 / (SizeX * SizeY * SizeZ));   
    } 
  else
    return;
}

// get the integrated probability density over z direction
//
// x, y = coordinates of the point
// sizeX, sizeY = dimensions of the system in x and y directions
// return = value of the probability density

double XYReflexionSymmetricPeriodicSpectra::PlanarProbabilityDensity(double x, double sizeX, double y, double sizeY)
{
  double* FunctionX = new double [this->NbrStateX];
  double* FunctionY = new double [this->NbrStateY];
  double FactorX = 2.0 * M_PI * x / sizeX; double FactorY = 2.0 * M_PI * y / sizeY; 
  for (int i = 0; i < this->NbrStateX; ++i)
    if (this->LowerImpulsionX == 0)
      FunctionX[i] = cos ((i + this->LowerImpulsionX) * FactorX);
    else
      FunctionX[i] = sin ((i + this->LowerImpulsionX) * FactorX);
  for (int j = 0; j < this->NbrStateY; ++j)
    if (this->LowerImpulsionY == 0)
      FunctionY[j] = cos ((j + this->LowerImpulsionY) * FactorY);
    else
      FunctionY[j] = sin ((j + this->LowerImpulsionY) * FactorY);

  double* TmpReal1; double* TmpImaginary1;  double* TmpReal2; double* TmpImaginary2;
  double TmpRe = 0.0, TmpIm = 0.0; double TmpRe1 = 0.0, TmpIm1 = 0.0; double TmpRe2 = 0.0, TmpIm2 = 0.0;
  for (int m1 = 0; m1 < this->NbrStateX; ++m1)
    for (int m2 = 0; m2 < this->NbrStateX; ++m2)
      {
	TmpRe1 = 0.0; TmpIm1 = 0.0;
	for (int n1 = 0; n1 < this->NbrStateY; ++n1)
	  {
	    TmpReal1 = this->RealCoefficients[m1][n1];
	    TmpImaginary1 = this->ImaginaryCoefficients[m1][n1];
	    for (int n2 = 0; n2 < this->NbrStateY; ++n2)
	      {
		TmpReal2 = this->RealCoefficients[m2][n2];
		TmpImaginary2 = this->ImaginaryCoefficients[m2][n2];
		TmpRe2 = 0.0; TmpIm2 = 0.0;
		for (int p = 0; p < this->NbrStateZ; ++p)
		  {
		    TmpRe2 += (TmpReal1[p] * TmpReal2[p] + TmpImaginary1[p] * TmpImaginary2[p]);
		    TmpIm2 += (TmpReal1[p] * TmpImaginary2[p] - TmpImaginary1[p] * TmpReal2[p]);
		  }
		TmpRe1 += (TmpRe2 * FunctionY[n1] * FunctionY[n2]);
		TmpRe1 += (TmpIm2 * FunctionY[n1] * FunctionY[n2]);
	      }
	  }
	TmpRe += (TmpRe1 * FunctionX[m1] * FunctionX[m2]);
	TmpIm += (TmpIm1 * FunctionX[m1] * FunctionX[m2]);
      }
  return (TmpRe * TmpRe + TmpIm * TmpIm) / sqrt (sizeX * sizeY) ;
}

// get the value of impulsion operators with another wavefunction <this|p|another>
//
// space = Hilbert space describing the other particle
// fileName = the file to stock the other function
// sizeX, sizeY, sizeZ = size of sample in X, Y and Z directions
// impulsionX, impulsionY, impulsionZ = reference to the return values

void XYReflexionSymmetricPeriodicSpectra::GetImpulsion(PeriodicXYReflexionZPeriodicThreeDOneParticle* space, char* fileName, double sizeX, double sizeY, double sizeZ, double &realImpulsionX, double &imaginaryImpulsionX, double &realImpulsionY, double &imaginaryImpulsionY, double &realImpulsionZ, double &imaginaryImpulsionZ)
{
  int nbrStateX = space->GetNbrStateX();
  int nbrStateY = space->GetNbrStateY();
  int nbrStateZ = space->GetNbrStateZ();
  int lowerImpulsionX = space->GetLowerImpulsionX();
  int lowerImpulsionY = space->GetLowerImpulsionY();
  int lowerImpulsionZ = space->GetLowerImpulsionZ();

  int MaxLowerX = this->LowerImpulsionX, MaxLowerY = this->LowerImpulsionY, MaxLowerZ = this->LowerImpulsionZ;
  if (this->LowerImpulsionX < lowerImpulsionX)
    MaxLowerX = lowerImpulsionX;
  if (this->LowerImpulsionY < lowerImpulsionY)
    MaxLowerY = lowerImpulsionY;
  if (this->LowerImpulsionZ < lowerImpulsionZ)
    MaxLowerZ = lowerImpulsionZ;

  int MinUpperX = this->LowerImpulsionX + this->NbrStateX, MinUpperY = this->LowerImpulsionY + this->NbrStateY, MinUpperZ = this->LowerImpulsionZ + this->NbrStateZ;
  if (MinUpperX > (lowerImpulsionX + nbrStateX))
    MinUpperX = lowerImpulsionX + nbrStateX;
  if (MinUpperY > (lowerImpulsionY + nbrStateY))
    MinUpperY = lowerImpulsionY + nbrStateY;
  if (MinUpperZ > (lowerImpulsionZ + nbrStateZ))
    MinUpperZ = lowerImpulsionZ + nbrStateZ;

  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Error in open the file: " << fileName << endl;
      exit(0);
    }
  double*** realCoefficients = new double** [nbrStateX];
  double*** imaginaryCoefficients = new double** [nbrStateX];
  for (int i = 0; i < nbrStateX; ++i)
    {
      realCoefficients[i] = new double* [nbrStateY];
      imaginaryCoefficients[i] = new double* [nbrStateY];   
      for (int j = 0; j < nbrStateY; ++j)
	{
	  realCoefficients[i][j] = new double [nbrStateZ];
	  imaginaryCoefficients[i][j] = new double [nbrStateZ];
	  for (int k = 0; k < nbrStateZ; ++k)	      
	    File >> realCoefficients[i][j][k] >> imaginaryCoefficients[i][j][k];	    
	}
    }
  File.close();  


  realImpulsionX = 0.0;
  realImpulsionY = 0.0;
  realImpulsionZ = 0.0;
  imaginaryImpulsionX = 0.0;
  imaginaryImpulsionY = 0.0;
  imaginaryImpulsionZ = 0.0;
  double TmpRe = 0.0, TmpIm = 0.0;
  double* Re1; double* Im1; double* Re2; double* Im2;
  // x direction
  if ((this->LowerImpulsionY != lowerImpulsionY) || (this->LowerImpulsionX == lowerImpulsionX))
    {
      realImpulsionX = 0.0;
      imaginaryImpulsionX = 0.0;
    }
  else
    {
      for (int m = MaxLowerX; m < MinUpperX; ++m)
	{
	  TmpRe = 0.0; TmpIm = 0.0;
	  for (int n = MaxLowerY; n < MinUpperY; ++n)
	    {
	      Re1 = this->RealCoefficients[m - this->LowerImpulsionX][n - this->LowerImpulsionY];
	      Im1 = this->ImaginaryCoefficients[m - this->LowerImpulsionX][n - this->LowerImpulsionY];
	      Re2 = realCoefficients[m - lowerImpulsionX][n - lowerImpulsionY];
	      Im2 = imaginaryCoefficients[m - lowerImpulsionX][n - lowerImpulsionY];
	      for (int p = MaxLowerZ; p < MinUpperZ; ++p)
		{
		  TmpRe += (Re1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]);
		  TmpIm += (Re1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] - Im1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]);
		}
	    }
	  realImpulsionX += ((double) m * TmpRe); imaginaryImpulsionX += ((double) m * TmpIm);
	}
      realImpulsionX = realImpulsionX * 2.0 * M_PI / sizeX; imaginaryImpulsionX = imaginaryImpulsionX * 2.0 * M_PI / sizeX;
    }

  // y direction
  if ((this->LowerImpulsionX != lowerImpulsionX) || (this->LowerImpulsionY == lowerImpulsionY))
    {
      realImpulsionY = 0.0;
      imaginaryImpulsionY = 0.0;
    }
  else
    {
      for (int m = MaxLowerX; m < MinUpperX; ++m)
	{
	  for (int n = MaxLowerY; n < MinUpperY; ++n)
	    {
	      TmpRe = 0.0; TmpIm = 0.0;
	      Re1 = this->RealCoefficients[m - this->LowerImpulsionX][n - this->LowerImpulsionY];
	      Im1 = this->ImaginaryCoefficients[m - this->LowerImpulsionX][n - this->LowerImpulsionY];
	      Re2 = realCoefficients[m - lowerImpulsionX][n - lowerImpulsionY];
	      Im2 = imaginaryCoefficients[m - lowerImpulsionX][n - lowerImpulsionY];
	      for (int p = MaxLowerZ; p < MinUpperZ; ++p)
		{
		  TmpRe += (Re1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]);
		  TmpIm += (Re1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] - Im1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]);
		}
	      realImpulsionY += (double(n) * TmpRe); imaginaryImpulsionY += (double(n) * TmpIm);
	    }
	}
      realImpulsionY = realImpulsionY * 2.0 * M_PI / sizeY; imaginaryImpulsionY = imaginaryImpulsionY * 2.0 * M_PI / sizeY;
    }

  // z direction
  if ((this->LowerImpulsionX != lowerImpulsionX) || (this->LowerImpulsionY != lowerImpulsionY))
    {
      realImpulsionZ = 0.0;
      imaginaryImpulsionZ = 0.0;
    }
  else
    {    
      for (int m = MaxLowerX; m < MinUpperX; ++m)
	{
	  for (int n = MaxLowerY; n < MinUpperY; ++n)
	    {	      
	      Re1 = this->RealCoefficients[m - this->LowerImpulsionX][n - this->LowerImpulsionY];
	      Im1 = this->ImaginaryCoefficients[m - this->LowerImpulsionX][n - this->LowerImpulsionY];
	      Re2 = realCoefficients[m - lowerImpulsionX][n - lowerImpulsionY];
	      Im2 = imaginaryCoefficients[m - lowerImpulsionX][n - lowerImpulsionY];
	      for (int p = MaxLowerZ; p < MinUpperZ; ++p)
		{
		  realImpulsionZ += ((Re1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]) * double(p));
		  imaginaryImpulsionZ += ((Re1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] - Im1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]) * double(p));
		}
	    }
	}
      realImpulsionZ = realImpulsionZ * 2.0 * M_PI / sizeZ; imaginaryImpulsionZ = imaginaryImpulsionZ * 2.0 * M_PI / sizeZ;	    
    }  
}

// get the value of impulsion  with a continuum state, described by Z plane waves
//
// impulsionX, impulsionY = value of the planer impulsion
// nbrStateZ, lowerImpulsionZ = Hilbert space characteristics 
// fileName =  the file to stock the other function
// sizeZ = size of the sample in Z direction
// impulsionZ = reference to the return values

void XYReflexionSymmetricPeriodicSpectra::GetImpulsionWithContinuum(int impulsionX, int impulsionY, int nbrStateZ, int lowerImpulsionZ, char* fileName, double sizeZ, double &realImpulsionZ, double &imaginaryImpulsionZ)
{
  int MaxLowerZ = this->LowerImpulsionZ;
  if (this->LowerImpulsionZ < lowerImpulsionZ)
    MaxLowerZ = lowerImpulsionZ;

  int MinUpperZ = this->LowerImpulsionZ + this->NbrStateZ;
  if (MinUpperZ > (lowerImpulsionZ + nbrStateZ))
    MinUpperZ = lowerImpulsionZ + nbrStateZ;

  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Error in open the file: " << fileName << endl;
      exit(0);
    }
  double* realCoefficients = new double [nbrStateZ];
  double* imaginaryCoefficients = new double [nbrStateZ];
  for (int k = 0; k < nbrStateZ; ++k)	      
    File >> realCoefficients[k] >> imaginaryCoefficients[k];	    
  File.close();   

  int deltaX = impulsionX - this->LowerImpulsionX; int deltaY = impulsionY - this->LowerImpulsionY; 
  if ((deltaX < 0) || (deltaX > this->NbrStateX) || (deltaY < 0) || (deltaY >= this->NbrStateY))
    {
      realImpulsionZ = 0.0; imaginaryImpulsionZ = 0.0;
      return;
    }
  double* Re1 = this->RealCoefficients[deltaX][deltaY];
  double* Im1 = this->ImaginaryCoefficients[deltaX][deltaY];
  realImpulsionZ = 0.0; imaginaryImpulsionZ = 0.0;
  for (int p = MaxLowerZ; p < MinUpperZ; ++p)
    { 
      realImpulsionZ += ((Re1[p - this->LowerImpulsionZ] * realCoefficients[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * imaginaryCoefficients[p - lowerImpulsionZ]) * double(p));
      imaginaryImpulsionZ += ((Re1[p - this->LowerImpulsionZ] * imaginaryCoefficients[p - lowerImpulsionZ] - Im1[p - this->LowerImpulsionZ] * realCoefficients[p - lowerImpulsionZ]) * double(p));
    }
  realImpulsionZ = realImpulsionZ * 2.0 * M_PI / sizeZ; imaginaryImpulsionZ = imaginaryImpulsionZ * 2.0 * M_PI / sizeZ;	   
}

// get the overlap of derived functions
//
// space = Hilbert space describing the other particle
// fileName = the file to stock the other function
// sizeX, sizeY, sizeZ = size of sample in X, Y and Z directions
// overlap, overlapX, overlapY = reference to the return values

void XYReflexionSymmetricPeriodicSpectra::GetDerivedOverlap (PeriodicXYReflexionZPeriodicThreeDOneParticle* space, char* fileName, double sizeX, double sizeY, double sizeZ, double &realOverlap, double &imaginaryOverlap, double &realOverlapX, double &imaginaryOverlapX, double &realOverlapY, double &imaginaryOverlapY)
{
  int nbrStateX = space->GetNbrStateX();
  int nbrStateY = space->GetNbrStateY();
  int nbrStateZ = space->GetNbrStateZ();
  int lowerImpulsionX = space->GetLowerImpulsionX();
  int lowerImpulsionY = space->GetLowerImpulsionY();
  int lowerImpulsionZ = space->GetLowerImpulsionZ();

  int MaxLowerX = this->LowerImpulsionX, MaxLowerY = this->LowerImpulsionY, MaxLowerZ = this->LowerImpulsionZ;
  if (this->LowerImpulsionX < lowerImpulsionX)
    MaxLowerX = lowerImpulsionX;
  if (this->LowerImpulsionY < lowerImpulsionY)
    MaxLowerY = lowerImpulsionY;
  if (this->LowerImpulsionZ < lowerImpulsionZ)
    MaxLowerZ = lowerImpulsionZ;

  int MinUpperX = this->LowerImpulsionX + this->NbrStateX, MinUpperY = this->LowerImpulsionY + this->NbrStateY, MinUpperZ = this->LowerImpulsionZ + this->NbrStateZ;
  if (MinUpperX > (lowerImpulsionX + nbrStateX))
    MinUpperX = lowerImpulsionX + nbrStateX;
  if (MinUpperY > (lowerImpulsionY + nbrStateY))
    MinUpperY = lowerImpulsionY + nbrStateY;
  if (MinUpperZ > (lowerImpulsionZ + nbrStateZ))
    MinUpperZ = lowerImpulsionZ + nbrStateZ;

  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Error in open the file: " << fileName << endl;
      exit(0);
    }
  double*** realCoefficients = new double** [nbrStateX];
  double*** imaginaryCoefficients = new double** [nbrStateX];
  for (int i = 0; i < nbrStateX; ++i)
    {
      realCoefficients[i] = new double* [nbrStateY];
      imaginaryCoefficients[i] = new double* [nbrStateY];   
      for (int j = 0; j < nbrStateY; ++j)
	{
	  realCoefficients[i][j] = new double [nbrStateZ];
	  imaginaryCoefficients[i][j] = new double [nbrStateZ];
	  for (int k = 0; k < nbrStateZ; ++k)	      
	    File >> realCoefficients[i][j][k] >> imaginaryCoefficients[i][j][k];	    
	}
    }
  File.close();  

  realOverlap = 0.0; imaginaryOverlap = 0.0; realOverlapX = 0.0; imaginaryOverlapX = 0.0; realOverlapY = 0.0; imaginaryOverlapY = 0.0; 

  for (int m = 0; m < this->NbrStateX; ++m)
    for (int n = 0; n < this->NbrStateY; ++n)
      for (int p = 0; p < this->NbrStateZ; ++p)
	{
	  realOverlap += (this->RealCoefficients[m][n][p] * realCoefficients[m][n][p] + this->ImaginaryCoefficients[m][n][p] * imaginaryCoefficients[m][n][p]);
	  imaginaryOverlap += (this->RealCoefficients[m][n][p] * imaginaryCoefficients[m][n][p] - this->ImaginaryCoefficients[m][n][p] * realCoefficients[m][n][p]);
	}

  for (int m = 0; m < this->NbrStateX; ++m)
    for (int n = 0; n < this->NbrStateY; ++n)
      for (int p = 0; p < this->NbrStateZ; ++p)
	{
	  realOverlapX += (this->RealCoefficients[m][n][p] * realCoefficients[m][n][p] + this->ImaginaryCoefficients[m][n][p] * imaginaryCoefficients[m][n][p]) * double (m + this->LowerImpulsionX) * double (m + this->LowerImpulsionX);
	  imaginaryOverlapX += (this->RealCoefficients[m][n][p] * imaginaryCoefficients[m][n][p] - this->ImaginaryCoefficients[m][n][p] * realCoefficients[m][n][p]) * double (m + this->LowerImpulsionX) * double (m + this->LowerImpulsionX);
	} 

  for (int m = 0; m < this->NbrStateX; ++m)
    for (int n = 0; n < this->NbrStateY; ++n)
      for (int p = 0; p < this->NbrStateZ; ++p)
	{
	  realOverlapY += (this->RealCoefficients[m][n][p] * realCoefficients[m][n][p] + this->ImaginaryCoefficients[m][n][p] * imaginaryCoefficients[m][n][p]) * double (n + this->LowerImpulsionY) * double (n + this->LowerImpulsionY);
	  imaginaryOverlapY += (this->RealCoefficients[m][n][p] * imaginaryCoefficients[m][n][p] - this->ImaginaryCoefficients[m][n][p] * realCoefficients[m][n][p]) * double (n + this->LowerImpulsionY) * double (n + this->LowerImpulsionY);
	} 
  realOverlapX *= (4.0 * M_PI * M_PI / (sizeX * sizeX)); imaginaryOverlapX *= (4.0 * M_PI * M_PI / (sizeX * sizeX)); 
  realOverlapY *= (4.0 * M_PI * M_PI / (sizeY * sizeY)); imaginaryOverlapY *= (4.0 * M_PI * M_PI / (sizeY * sizeY)); 
}

// get the probability of finding the particle in a cube
//
// minX, maxX, minY, maxY, minZ, maxZ = bounds of the cube in unit of proportion compared to the whole length
// return = value of the probability  

double XYReflexionSymmetricPeriodicSpectra::GetCubeProbability (double minX, double maxX, double minY, double maxY, double minZ, double maxZ)
{
  
}
