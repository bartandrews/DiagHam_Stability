
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//        class for periodic average spectra with Fourier-Bessel basis        //
//                                                                            //
//                        last modification : 07/05/2004                      //
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


#include "Tools/Spectra/CylinderQuantumDotSpectra.h"
#include "MathTools/BesselJZeros.h"
#include "Tools/Potential/ThreeDConstantCylinderPotential.h"
#include "Tools/Potential/QuantumDotThreeDConstantCylinderPotential.h"

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
// bz = magnetic field in Z direction
CylinderQuantumDotSpectra::CylinderQuantumDotSpectra(PlanarRotationSymmetryZPeriodicOneParticle* space, char* fileName, double bz)
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

// get the value of impulsion operators with another wavefunction <this|p|another>
//
// space = Hilbert space describing the other particle
// fileName = the file to stock the other function
// sizeR = size of the super-cylinder in plane
// sizeZ = size of sample in Z direction
// impulsionX, impulsionY, impulsionZ = reference to the return values

void CylinderQuantumDotSpectra::GetImpulsion(PlanarRotationSymmetryZPeriodicOneParticle* space, char* fileName, double sizeR, double sizeZ, double &realImpulsionX, double &imaginaryImpulsionX, double &realImpulsionY, double &imaginaryImpulsionY, double &realImpulsionZ, double &imaginaryImpulsionZ)
{
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
  double* Re1; double* Im1; double* Re2; double* Im2;
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
      for (int n1 = 0; n1 < MinUpperR; ++n1)
	{
	  Re1 = this->RealCoefficients[n1];
	  Im1 = this->ImaginaryCoefficients[n1];
	  for (int n2 = 0; n2 < MinUpperR; ++n2)
	    {
	      Re2 = realCoefficients[n2];
	      Im2 = imaginaryCoefficients[n2];
	      TmpRe1 = 0.0; TmpIm1 = 0.0;
	      for (int p = MaxLowerZ; p < MinUpperZ; ++p) 
		{
		  TmpRe1 += (Re1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]);
		  TmpIm1 += (-Re1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]);
		}
	      double fraction = BesselJZeros[0][n1] / BesselJZeros[1][n2];
	      TmpRe += (TmpRe1 / (fraction - 1.0 / fraction));
	      TmpIm += (TmpIm1 / (fraction - 1.0 / fraction));
	    }
	}
      realImpulsionX = TmpRe / sizeR; realImpulsionY = TmpRe / sizeR; 
      imaginaryImpulsionX = TmpIm / sizeR; imaginaryImpulsionY = TmpIm / sizeR; 
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

// get the probability integrated in the dot to find the particle
//
// potential = pointer to a 3D potential with constant value in a cylinder

double CylinderQuantumDotSpectra::GetDotProbability(QuantumDotThreeDConstantCylinderPotential* potential)
{
  int addition;
  if ((this->NumberM == 1) || (this->NumberM == -1))
    addition = 1;
  else
    if (this->NumberM == 0)
      addition = 0;
    else
      cout << "Attention! This quantum number of kinetic momentum has not been taken into account: " << this->NumberM << endl;

  double** RealWaveFunctionOverlapZ; double** ImaginaryWaveFunctionOverlapZ;
  if (!this->EvaluatePlaneWaveFunctionOverlap(potential, this->NbrStateZ, RealWaveFunctionOverlapZ, ImaginaryWaveFunctionOverlapZ))
    cout << "Error in evaluation of function overlap in Z direction. Stop!" << endl;
  
  double WettingRadius = potential->GetRadius(3) * 3; // wetting radius
  int nbrCylinder = potential->GetNbrCylinderZ();
  double RSize = potential->GetSuperCylinderRadius();
  double radius;
  double* Zeros = new double [this->NbrStateR];
  double* Normalization = new double [this->NbrStateR];
  double** Fraction = new double* [this->NbrStateR];
  double** BesselOne = new double* [this->NbrStateR];
  double** BesselTwo = new double* [this->NbrStateR];
  for (int n = 0; n < this->NbrStateR; ++n)
    {
      Zeros[n] = BesselJZeros[addition][n];
      if (addition == 0)
	Normalization[n] = 1.0 / j1(Zeros[n]);
      if (addition == 1)
	Normalization[n] = 1.0 / j0(Zeros[n]);

      Fraction[n] = new double [nbrCylinder];
      BesselOne[n] = new double [nbrCylinder];
      BesselTwo[n] = new double [nbrCylinder];
      for (int k = 0; k < nbrCylinder; ++k)
	{
	  if (k != 2)
	    radius = potential->GetRadius(k);
	  else 
	    radius = WettingRadius; // for the wetting layer!!!
	  if (radius > 0.0)
	    {
	      Fraction[n][k] = radius * Zeros[n] / RSize;
	      if (addition == 0)
		{
		  BesselOne[n][k] = j0(Fraction[n][k]);
		  BesselTwo[n][k] = -j1(Fraction[n][k]);
		}
	      if (addition == 1)
		{
		  BesselOne[n][k] = j1(Fraction[n][k]);
		  BesselTwo[n][k] = j0(Fraction[n][k]);
		}		
	    }
	}
    }
  double*** RadialOverlap = new double ** [this->NbrStateR];
  for (int n1 = 0; n1 < this->NbrStateR; ++n1)
    {
      RadialOverlap[n1] = new double* [this->NbrStateR];
      for (int n2 = 0; n2 < this->NbrStateR; ++n2)
	{ 
	  RadialOverlap[n1][n2] = new double [nbrCylinder];
	  for (int k = 0; k < nbrCylinder; ++k)
	    {
	      if (k != 2)
		radius = potential->GetRadius(k);
	      else 
		radius = WettingRadius; // for the wetting layer!!!
	      if (radius > 0.0)
		{
		  if (n1 != n2)		    
		    RadialOverlap[n1][n2][k] = 2.0 * (Fraction[n1][k] * BesselOne[n2][k] * BesselTwo[n1][k] - Fraction[n2][k] * BesselOne[n1][k] * BesselTwo[n2][k]) * (Normalization[n1] * Normalization[n2] / (Zeros[n2] * Zeros[n2] - Zeros[n1] * Zeros[n1]));
		  else
		    RadialOverlap[n1][n2][k] = ((radius * radius) / (RSize * RSize)) * (BesselTwo[n1][k] * BesselTwo[n1][k] - (2.0 * double(addition) * BesselOne[n1][k] * BesselTwo[n1][k] / Fraction[n1][k]) + BesselOne[n1][k] * BesselOne[n1][k]) * Normalization[n1] * Normalization[n1];
		}	      
	      else
		RadialOverlap[n1][n2][k] = 0;
	    }	  	  
	}
    }
  
  int IndexZ; int LimitZ = 0;
  int LengthZ = (this->NbrStateZ - 1) * 2 + 1; int OriginZ = this->NbrStateZ - 1;
  double* tmpRadial; double* tmpRealVertical; double* tmpImaginaryVertical;
  double TmpRe = 0.0; double TmpIm = 0.0; double TmpRe1 = 0.0; double TmpIm1 = 0.0; double TmpRe2 = 0.0; double TmpIm2 = 0.0;

  for (int n1 = 0; n1 < this->NbrStateR; ++n1)    
    for (int n2 = 0; n2 < this->NbrStateR; ++n2)
      { 
	tmpRadial =  RadialOverlap[n1][n2];
	for (int p1 = 0; p1 < this->NbrStateZ; ++p1)
	  {	   
	    IndexZ = -p1 + OriginZ;
	    LimitZ = LengthZ - p1;
	    for (int p2 = 0; IndexZ < LimitZ; ++IndexZ, ++p2)
	      {
		tmpRealVertical = RealWaveFunctionOverlapZ [IndexZ];
		tmpImaginaryVertical = ImaginaryWaveFunctionOverlapZ [IndexZ];
		TmpRe1 = 0.0; TmpIm1 = 0.0;
		for (int k = 0; k < nbrCylinder; ++k)		  
		  {  
		    TmpRe1 += (tmpRadial[k] * tmpRealVertical[k]);
		    TmpIm1 += (tmpRadial[k] * tmpImaginaryVertical[k]);
		  }
		TmpRe2 = (this->RealCoefficients[n1][p1] * this->RealCoefficients[n2][p2] + this->ImaginaryCoefficients[n1][p1] * this->ImaginaryCoefficients[n2][p2]);
		TmpIm2 = (this->RealCoefficients[n1][p1] * this->ImaginaryCoefficients[n2][p2] - this->ImaginaryCoefficients[n1][p1] * this->RealCoefficients[n2][p2]);
		TmpRe += (TmpRe2 * TmpRe1 - TmpIm2 * TmpIm1);
		//TmpIm += (TmpRe2 * TmpIm1 + TmpIm2 * TmpRe1);
	      }
	  }
      }    
  return TmpRe;
}

// evaluate the plane wave function overlap
//
// potential = pointer to the potential
// nbrState = number of states chosen for this direction
// realArray = 2D array containing the real elements of the overlap
// imaginaryArray = 2D array containing the imaginary elements of the overlap

bool CylinderQuantumDotSpectra::EvaluatePlaneWaveFunctionOverlap(QuantumDotThreeDConstantCylinderPotential* &potential, int nbrState, double** &realArray, double** &imaginaryArray)
{
  int nbrCylinder = potential->GetNbrCylinderZ();
  double* ZPosition = new double [nbrCylinder + 1];
  ZPosition[0] = 0.0;
  for (int k = 0; k < nbrCylinder; ++k)    
    ZPosition[k + 1] = ZPosition[k] + potential->GetHeight(k);      
  double ZSize = ZPosition[nbrCylinder];
      
  realArray = new double* [nbrState * 2 - 1];
  imaginaryArray = new double* [nbrState * 2 - 1];   
  int Origin = nbrState - 1;
  double Diff = 0.0, Tmp = 0.0;
  for (int p = 0; p < (nbrState * 2 - 1); ++p)
    {      
      realArray[p] = new double [nbrCylinder];
      imaginaryArray[p] = new double [nbrCylinder];
      if (p == Origin)	
	for (int k = 0; k < nbrCylinder; ++k)
	  {
	    realArray[Origin][k] = (ZPosition[k + 1] - ZPosition[k]) / ZSize;
	    imaginaryArray[Origin][k] = 0.0;     
	  }
      else
	{
	  Diff = 2.0 * M_PI * double(p - Origin);
	  Tmp = Diff / ZSize;
	  Diff = 1.0 / Diff;
	  for (int k = 0; k < nbrCylinder; ++k)
	    {
	      realArray[p][k] = Diff * (sin(Tmp * ZPosition[k + 1]) - sin(Tmp * ZPosition[k]));
	      imaginaryArray[p][k] = Diff * (cos(Tmp * ZPosition[k]) - cos(Tmp * ZPosition[k + 1]));     
	    }
	}
    }
  return true;
}

