////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                         class for periodic  spectra                        //
//                                                                            //
//                        last modification : 12/20/2003                      //
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


#include "config.h"
#include "Tools/Spectra/PeriodicSpectra.h"
#include "Tools/Potential/TetrapodThreeDConstantCellPotential.h"


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
// space: Hilbert space describing the particle
// fileName: name of the state file

PeriodicSpectra::PeriodicSpectra(PeriodicThreeDOneParticle* space, char* fileName)
{
  this->NbrStateX = space->GetNbrStateX();
  this->NbrStateY = space->GetNbrStateY();
  this->NbrStateZ = space->GetNbrStateZ();
  this->LowerImpulsionX = space->GetLowerImpulsionX();
  this->LowerImpulsionY = space->GetLowerImpulsionY();
  this->LowerImpulsionZ = space->GetLowerImpulsionZ();
  
  int LengthX = (this->NbrStateX - 1) * 2 + 1; int LengthY = (this->NbrStateY - 1) * 2 + 1; int LengthZ = (this->NbrStateZ - 1) * 2 + 1;
  int OriginX = this->NbrStateX - 1; int OriginY = this->NbrStateY - 1; int OriginZ = this->NbrStateZ - 1;
  this->RealWaveFunctionOverlapX = new double [LengthX];
  this->ImaginaryWaveFunctionOverlapX = new double [LengthX];
  this->RealSquareOverlapX = new double [LengthX];
  this->ImaginarySquareOverlapX = new double [LengthX];
  double Tmp = 0.0;
  for (int i = 0; i < LengthX; ++i)
    {
      if (i != OriginX)
	{
	  Tmp = 1.0 / (2.0 * M_PI * double(i - OriginX));
	  this->RealWaveFunctionOverlapX[i] = 0.0;
	  this->ImaginaryWaveFunctionOverlapX[i] = -Tmp;
	  this->RealSquareOverlapX[i] = Tmp * Tmp * 2.0;
	  this->ImaginarySquareOverlapX[i] = -Tmp;
	}
      else
	{
	  this->RealWaveFunctionOverlapX[i] = 0.5;
	  this->ImaginaryWaveFunctionOverlapX[i] = 0.0;
	  this->RealSquareOverlapX[i] = 1.0 / 3.0;
	  this->ImaginarySquareOverlapX[i] = 0.0;
	}
    }  
  this->RealWaveFunctionOverlapY = new double [LengthY];
  this->ImaginaryWaveFunctionOverlapY = new double [LengthY];
  this->RealSquareOverlapY = new double [LengthY];
  this->ImaginarySquareOverlapY = new double [LengthY];
  for (int i = 0; i < LengthY; ++i)
    {
      if (i != OriginY)
	{
	  Tmp =  1.0 / (2.0 * M_PI * double(i - OriginY));
	  this->RealWaveFunctionOverlapY[i] = 0.0;
	  this->ImaginaryWaveFunctionOverlapY[i] = -Tmp;
	  this->RealSquareOverlapY[i] = Tmp * Tmp * 2.0;
	  this->ImaginarySquareOverlapY[i] = -Tmp;
	}
      else
	{
	  this->RealWaveFunctionOverlapY[i] = 0.5;
	  this->ImaginaryWaveFunctionOverlapY[i] = 0.0;
	  this->RealSquareOverlapY[i] = 1.0 / 3.0;
	  this->ImaginarySquareOverlapY[i] = 0.0;
	}
    } 
  this->RealWaveFunctionOverlapZ = new double [LengthZ];
  this->ImaginaryWaveFunctionOverlapZ = new double [LengthZ];
  this->RealSquareOverlapZ = new double [LengthZ];
  this->ImaginarySquareOverlapZ = new double [LengthZ];
  for (int i = 0; i < LengthZ; ++i)
    {
      if (i != OriginZ)
	{
	  Tmp = 1.0 / (2.0 * M_PI * double(i - OriginZ));
	  this->RealWaveFunctionOverlapZ[i] = 0.0;
	  this->ImaginaryWaveFunctionOverlapZ[i] = -Tmp;
	  this->RealSquareOverlapZ[i] = Tmp * Tmp * 2.0;
	  this->ImaginarySquareOverlapZ[i] = -Tmp;
	}
      else
	{
	  this->RealWaveFunctionOverlapZ[i] = 0.5;
	  this->ImaginaryWaveFunctionOverlapZ[i] = 0.0;
 	  this->RealSquareOverlapZ[i] = 1.0 / 3.0;
	  this->ImaginarySquareOverlapZ[i] = 0.0;
	}
    } 

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

// get mean value in X direction
//
// return: position in 1.0 scale

double PeriodicSpectra::GetMeanValueX(double& squareX)
{
  int OriginX = this->NbrStateX - 1;
  double TmpRe = 0.0, TmpIm = 0.0;
  double Tmp = 0.0;
  squareX = 0.0;
  for (int m1 = 0; m1 < this->NbrStateX; ++m1)
    for (int m2 = 0; m2 < this->NbrStateX; ++m2)
      {
	TmpRe = 0.0; TmpIm = 0.0;
	for (int n = 0; n < this->NbrStateY; ++n)
	  for (int p = 0; p < this->NbrStateZ; ++p)
	    {
	      TmpRe += (RealCoefficients[m1][n][p] * RealCoefficients[m2][n][p] + ImaginaryCoefficients[m1][n][p] * ImaginaryCoefficients[m2][n][p]);
	      TmpIm += (RealCoefficients[m1][n][p] * ImaginaryCoefficients[m2][n][p] - ImaginaryCoefficients[m1][n][p] * RealCoefficients[m2][n][p]);
	    }
	Tmp += (TmpRe * this->RealWaveFunctionOverlapX[-m1 + m2 + OriginX] - TmpIm * this->ImaginaryWaveFunctionOverlapX[-m1 + m2 + OriginX]);   
  	squareX += (TmpRe * this->RealSquareOverlapX[m1 - m2 + OriginX] - TmpIm * this->ImaginarySquareOverlapX[-m1 + m2 + OriginX]);   
      }
  return Tmp;
}

// get mean value in Y direction
//
// return: position in 1.0 scale

double PeriodicSpectra::GetMeanValueY(double& squareY)
{
  int OriginY = this->NbrStateY - 1;
  double TmpRe = 0.0, TmpIm = 0.0;
  double Tmp = 0.0;
  squareY = 0.0;
  for (int n1 = 0; n1 < this->NbrStateY; ++n1)
    for (int n2 = 0; n2 < this->NbrStateY; ++n2)
      {
	TmpRe = 0.0; TmpIm = 0.0;
	for (int m = 0; m < this->NbrStateX; ++m)
	  for (int p = 0; p < this->NbrStateZ; ++p)
	    {
	      TmpRe += (RealCoefficients[m][n1][p] * RealCoefficients[m][n2][p] + ImaginaryCoefficients[m][n1][p] * ImaginaryCoefficients[m][n2][p]);
	      TmpIm += (RealCoefficients[m][n1][p] * ImaginaryCoefficients[m][n2][p] - ImaginaryCoefficients[m][n1][p] * RealCoefficients[m][n2][p]);
	    }
	Tmp += (TmpRe * this->RealWaveFunctionOverlapY[-n1 + n2 + OriginY] - TmpIm * this->ImaginaryWaveFunctionOverlapY[-n1 + n2 + OriginY]); 
 	squareY += (TmpRe * this->RealSquareOverlapY[-n1 + n2 + OriginY] - TmpIm * this->ImaginarySquareOverlapY[-n1 + n2 + OriginY]);        
      }
  return Tmp;
}

// get mean value in Z direction
//
// return: position in 1.0 scale

double PeriodicSpectra::GetMeanValueZ(double& squareZ)
{
  int OriginZ = this->NbrStateZ - 1;
  double TmpRe = 0.0, TmpIm = 0.0;
  double Tmp = 0.0;
  squareZ = 0.0;
  for (int p1 = 0; p1 < this->NbrStateZ; ++p1)
    for (int p2 = 0; p2 < this->NbrStateZ; ++p2)
      {
	TmpRe = 0.0; TmpIm = 0.0;
	for (int m = 0; m < this->NbrStateX; ++m)
	  for (int n = 0; n < this->NbrStateY; ++n)
	    {
	      TmpRe += (RealCoefficients[m][n][p1] * RealCoefficients[m][n][p2] + ImaginaryCoefficients[m][n][p1] * ImaginaryCoefficients[m][n][p2]);
	      TmpIm += (RealCoefficients[m][n][p1] * ImaginaryCoefficients[m][n][p2] - ImaginaryCoefficients[m][n][p1] * RealCoefficients[m][n][p2]);
	    }
	Tmp += (TmpRe * this->RealWaveFunctionOverlapZ[-p1 + p2 + OriginZ] - TmpIm * this->ImaginaryWaveFunctionOverlapZ[-p1 + p2 + OriginZ]);       
	squareZ += (TmpRe * this->RealSquareOverlapZ[-p1 + p2 + OriginZ] - TmpIm * this->ImaginarySquareOverlapZ[-p1 + p2 + OriginZ]);  
      }
  return Tmp;
}

// get the wave function value of a state at a given point
//
// x, y, z : the position of the point
// SizeX, SizeY, SizeZ : the 3D-sizes of the sample
// Real, Imaginary : references to the real and imaginary components of the wave function

void PeriodicSpectra::WaveFunctionValue(double x, double SizeX, double y, double SizeY, double z, double SizeZ, double& Real, double& Imaginary)
{
  double TmpX = 0.0; double TmpY = 0.0; double TmpZ = 0.0;
  double TmpRe = 0.0; double TmpIm = 0.0;
  double* TmpRealCoefficients; double* TmpImaginaryCoefficients;
  double TmpRe2 = 0.0; double TmpIm2 = 0.0;
  for (int m = 0; m < this->NbrStateX; ++m)
    {
      TmpX = double(m + this->LowerImpulsionX) * x / SizeX;
      for (int n = 0; n < this->NbrStateY; ++n)
	{
	  TmpY = double(n + this->LowerImpulsionY) * y /SizeY + TmpX;
	  TmpRealCoefficients = this->RealCoefficients[m][n];
	  TmpImaginaryCoefficients = this->ImaginaryCoefficients[m][n];
	  for (int p = 0; p < this->NbrStateZ; ++p)
	    {
	      TmpZ = 2 * M_PI * (TmpY + double(p + this->LowerImpulsionZ) * z / SizeZ);
	      TmpRe2 = cos(TmpZ); TmpIm2 = sin(TmpZ);
	      TmpRe += (TmpRealCoefficients[p] * TmpRe2 - TmpImaginaryCoefficients[p] * TmpIm2);
	      TmpIm += (TmpRealCoefficients[p] * TmpIm2 + TmpImaginaryCoefficients[p] * TmpRe2);
	    }
	}
    }
  Real = TmpRe / sqrt (SizeX * SizeY * SizeZ); Imaginary = TmpIm / sqrt (SizeX * SizeY * SizeZ);
}

// get the value of impulsion operators with another wavefunction <this|p|another>
//
// fileName = the file to stock the other function
// sizeX, sizeY, sizeZ = size of sample in X, Y and Z directions
// impulsionX, impulsionY, impulsionZ = reference to the return values

void PeriodicSpectra::GetImpulsion(char* fileName, double sizeX, double sizeY, double sizeZ, double &realImpulsionX, double &imaginaryImpulsionX, double &realImpulsionY, double &imaginaryImpulsionY, double &realImpulsionZ, double &imaginaryImpulsionZ)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Error in open the file: " << fileName << endl;
      exit(0);
    }
  realImpulsionX = 0.0;
  realImpulsionY = 0.0;
  realImpulsionZ = 0.0;
  imaginaryImpulsionX = 0.0;
  imaginaryImpulsionY = 0.0;
  imaginaryImpulsionZ = 0.0;

  double tmpRe = 0.0, tmpIm = 0.0;
  double tmpIm1 = 0.0, tmpRe1 = 0.0;
  double tmpIm2 = 0.0, tmpRe2 = 0.0;
  double tmpIm3 = 0.0, tmpRe3 = 0.0;
  double Im = 0.0, Re = 0.0;

  for (int m = 0; m < this->NbrStateX; ++m)
    {
      tmpRe1 = 0.0; tmpIm1 = 0.0;
      for (int n = 0; n < this->NbrStateY; ++n)
	{
	  tmpIm2 = 0.0; tmpRe2 = 0.0;
	  for (int p = 0; p < this->NbrStateZ; ++p)
	    {
	      File >> tmpRe >> tmpIm;
	      tmpRe3 = this->RealCoefficients[m][n][p] * tmpRe + this->ImaginaryCoefficients[m][n][p] * tmpIm;
	      tmpIm3 = this->RealCoefficients[m][n][p] * tmpIm - this->ImaginaryCoefficients[m][n][p] * tmpRe;
	      realImpulsionZ += (((double) p) * tmpRe3);
	      imaginaryImpulsionZ += (((double) p) * tmpIm3);
	      tmpRe2 += tmpRe3;
	      tmpIm2 += tmpIm3;
	    }
	  tmpRe1 += tmpRe2; tmpIm1 += tmpIm2;
	  realImpulsionY += double (n) * tmpRe2;
	  imaginaryImpulsionY += double (n) * tmpIm2;
	}
      Re += tmpRe1; Im += tmpIm1;
      realImpulsionX += double (m) * tmpRe1;
      imaginaryImpulsionX += double (m) * tmpIm1;
    }
  realImpulsionX = (realImpulsionX + double(this->LowerImpulsionX) * Re) * 2.0 * M_PI / sizeX;
  realImpulsionY = (realImpulsionY + double(this->LowerImpulsionY) * Re) * 2.0 * M_PI / sizeY;
  realImpulsionZ = (realImpulsionZ + double(this->LowerImpulsionZ) * Re) * 2.0 * M_PI / sizeZ;
  imaginaryImpulsionX = (imaginaryImpulsionX + double(this->LowerImpulsionX) * Im) * 2.0 * M_PI / sizeX;
  imaginaryImpulsionY = (imaginaryImpulsionY + double(this->LowerImpulsionY) * Im) * 2.0 * M_PI / sizeY;
  imaginaryImpulsionZ = (imaginaryImpulsionZ + double(this->LowerImpulsionZ) * Im) * 2.0 * M_PI / sizeZ;

  File.close();
}

// get the overlap of derived functions
//
// space = Hilbert space describing the other particle
// fileName = the file to stock the other function
// sizeX, sizeY, sizeZ = size of sample in X, Y and Z directions
// overlap, overlapX, overlapY = reference to the return values

void PeriodicSpectra::GetDerivedOverlap (PeriodicThreeDOneParticle* space, char* fileName, double sizeX, double sizeY, double sizeZ, double &realOverlap, double &imaginaryOverlap, double &realOverlapX, double &imaginaryOverlapX, double &realOverlapY, double &imaginaryOverlapY)
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

double PeriodicSpectra::GetCubeProbability (double minX, double maxX, double minY, double maxY, double minZ, double maxZ)
{
  Complex** overlapX = new Complex* [this->NbrStateX];
  Complex** overlapY = new Complex* [this->NbrStateY];
  Complex** overlapZ = new Complex* [this->NbrStateZ];

  for (int m1 = 0; m1 < this->NbrStateX; ++m1)    
    {
      overlapX [m1] = new Complex [this->NbrStateX];
      for (int m2 = 0; m2 < this->NbrStateX; ++m2)
	this->EvaluateOneDOverlap(m1, m2, minX, maxX, overlapX[m1][m2]);    
    }

  for (int n1 = 0; n1 < this->NbrStateY; ++n1)    
    {
      overlapY [n1] = new Complex [this->NbrStateY];
      for (int n2 = 0; n2 < this->NbrStateY; ++n2)
	this->EvaluateOneDOverlap(n1, n2, minY, maxY, overlapY[n1][n2]);    
    }

  for (int p1 = 0; p1 < this->NbrStateZ; ++p1)    
    {
      overlapZ [p1] = new Complex [this->NbrStateZ];
      for (int p2 = 0; p2 < this->NbrStateZ; ++p2)
	this->EvaluateOneDOverlap(p1, p2, minZ, maxZ, overlapZ[p1][p2]);    
    }

  Complex tmp1 = 0.0, tmp2 = 0.0; double value = 0.0;
  double* tmpRe1; double* tmpIm1; double* tmpRe2; double* tmpIm2; 
  Complex coef1 = 0.0, coef2 = 0.0; 
  for (int m1 = 0; m1 < this->NbrStateX; ++m1)    
    for (int m2 = 0; m2 <  this->NbrStateX; ++m2)
      {
	tmp1 = 0.0;
	for (int n1 = 0; n1 < this->NbrStateY; ++n1)    	    
	  {
	    tmpRe1 = this->RealCoefficients [m1][n1];
	    tmpIm1 = this->ImaginaryCoefficients [m1][n1];
	    for (int n2 = 0; n2 <  this->NbrStateY; ++n2)
	      {
		tmp2 = 0.0; 		
		tmpRe2 = this->RealCoefficients [m2][n2];
		tmpIm2 = this->ImaginaryCoefficients [m2][n2];
		for (int p1 = 0; p1 < this->NbrStateZ; ++p1)
		  {
		    coef1.Re = tmpRe1 [p1]; coef1.Im = tmpIm1 [p1];		 
		    for (int p2 = 0; p2 <  this->NbrStateZ; ++p2)
		      {
			coef2.Re = tmpRe2 [p2]; coef2.Im = tmpIm2 [p2];
			tmp2 += (coef1 * coef2 * overlapZ [p1][p2]);
		      }
		  }
		tmp1 += (tmp2 * overlapY [n1][n2]);	       
	      }
	  }
	value += (Real (tmp1) * Real (overlapX [m1][m2]) - Imag (tmp1) * Imag (overlapX [m1][m2]));
      }
  delete[] overlapX; delete[] overlapY; delete[] overlapZ; 

  return value;
}

// get the probability of the particle in different parts of the tetrapod 
//
// potential = pointer to the tetrapod potential
// SphereProbability = reference to the probability in the sphere
// ArmProbability = reference to the probabiliti in the four arms

void PeriodicSpectra::GetTetrapodProbability (TetrapodThreeDConstantCellPotential* potential, double& SphereProbability, double& ArmProbability)
{
  int NbrCellX = potential->GetNbrCellX ();
  int NbrCellY = potential->GetNbrCellY ();
  int NbrCellZ = potential->GetNbrCellZ ();

  double deltaX = 1.0 / ((double) NbrCellX);
  double deltaY = 1.0 / ((double) NbrCellY);
  double deltaZ = 1.0 / ((double) NbrCellZ);
  
  double minX, maxX, minY, maxY, minZ, maxZ;
  double sphere = 0.0; double arm = 0.0;
  for (int k = 0; k < NbrCellZ; ++k)
    {
      minZ = ((double) k) * deltaZ; maxZ = minZ + deltaZ;
      for (int j = 0; j < NbrCellY; ++j)
	{
	  minY = ((double) j) * deltaY; maxY = minY + deltaY;
	  for (int i = 0; i < NbrCellX; ++i)
	    {
	      minX = ((double) i) * deltaX; maxX = minX + deltaX;
	      if (potential->InTheDot (i, j, k))		
		sphere += this->GetCubeProbability (minX, maxX, minY, maxY, minZ, maxZ);
	      else
		if (potential->InTheArms (i, j, k))
		  arm += this->GetCubeProbability (minX, maxX, minY, maxY, minZ, maxZ);
	    }
	}
    }
  
  SphereProbability = sphere;; ArmProbability = arm;
}

// get the overlap value of one d plane wave functions
//
// m1, m2 = indices of the one d function
// min, max = bound of the integral, in unit of proportion compared to the whole length
// real, imaginary = real and imaginary values

void PeriodicSpectra::EvaluateOneDOverlap (int m1, int m2, double min, double max, Complex& c)
{
  if (m1 == m2)
    {
      c.Re = max - min;
      c.Im = 0.0;
    }
  else
    {
      double delta = 2.0 * M_PI * ((double) -m1 + m2);  
      c.Re = cos (delta * max) - cos (delta * min);
      c.Im = -sin (delta * max) + sin (delta * min);
      
      c.Re /= delta; c.Im /= delta; 
    }
  return;
}
