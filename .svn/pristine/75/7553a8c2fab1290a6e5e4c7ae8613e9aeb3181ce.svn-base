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


#include "Tools/Spectra/HardBoxSpectra.h"

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

HardBoxSpectra::HardBoxSpectra(ThreeDOneParticle* space, char* fileName)
{
  this->NbrStateX = space->GetNbrStateX();
  this->NbrStateY = space->GetNbrStateY();
  this->NbrStateZ = space->GetNbrStateZ();

  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Error in open the file: " << fileName << endl;
      exit(0);
    }

  this->Coefficients = new double** [this->NbrStateX];
  for (int i = 0; i < this->NbrStateX; ++i)
    {
      this->Coefficients[i] = new double* [this->NbrStateY];
      for (int j = 0; j < this->NbrStateY; ++j)
	{
	  this->Coefficients[i][j] = new double [this->NbrStateZ];
	  for (int k = 0; k < this->NbrStateZ; ++k)	    
	    File >> this->Coefficients[i][j][k];	    
	}
    } 
  File.close();
}

void HardBoxSpectra::GetDerivedOverlap (ThreeDOneParticle* space, char* fileName, double sizeX, double sizeY, double sizeZ, double &Overlap, double &OverlapX, double &OverlapY)
{
  int nbrStateX = space->GetNbrStateX();
  int nbrStateY = space->GetNbrStateY();
  int nbrStateZ = space->GetNbrStateZ();

  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Error in open the file: " << fileName << endl;
      exit(0);
    }
  double*** Coefficients = new double** [nbrStateX];
  for (int i = 0; i < nbrStateX; ++i)
    {
      Coefficients[i] = new double* [nbrStateY];
      for (int j = 0; j < nbrStateY; ++j)
	{
	  Coefficients[i][j] = new double [nbrStateZ];
	  for (int k = 0; k < nbrStateZ; ++k)	      
	    File >> Coefficients[i][j][k];	    
	}
    }
  File.close();  

  Overlap = 0.0; OverlapX = 0.0; OverlapY = 0.0;

  for (int m = 0; m < this->NbrStateX; ++m)
    for (int n = 0; n < this->NbrStateY; ++n)
      for (int p = 0; p < this->NbrStateZ; ++p)	
	Overlap += (this->Coefficients[m][n][p] * Coefficients[m][n][p]);	

  for (int m = 0; m < this->NbrStateX; ++m)
    for (int n = 0; n < this->NbrStateY; ++n)
      for (int p = 0; p < this->NbrStateZ; ++p)	
	OverlapX += (this->Coefficients[m][n][p] * Coefficients[m][n][p]) * double (m + 1) * double (m + 1);	

  for (int m = 0; m < this->NbrStateX; ++m)
    for (int n = 0; n < this->NbrStateY; ++n)
      for (int p = 0; p < this->NbrStateZ; ++p)	
	OverlapY += (this->Coefficients[m][n][p] * Coefficients[m][n][p]) * double (n + 1) * double (n + 1);

  OverlapX *= (M_PI * M_PI / (sizeX * sizeX)); 
  OverlapY *= (M_PI * M_PI / (sizeY * sizeY)); 
}
