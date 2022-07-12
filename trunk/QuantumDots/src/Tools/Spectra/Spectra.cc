////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                            base class for spectra                          //
//                                                                            //
//                        last modification : 09/15/2003                      //
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


#include "Tools/Spectra/Spectra.h"
#include "Vector/RealVector.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using std::ofstream;
using std::ios;
using std::cout;
using std::endl;


// default constructor
//

Spectra::Spectra()
{
  AxeX = 0;
  AxeY = 0;
  PointNumber = 0;
}

// constructor
//
// N = number of points

Spectra::Spectra(int N)
{
  AxeX = new RealVector(N);
  AxeY = new RealVector(N);
  PointNumber = N;
}

// constructor from a set of files with 2 column data. Each peak is assimilated to a Lorentzian function with FWHM = Gamma
//
//
// FileNumber = number of files, Files = name of files
// LineNumber = number of lines, Gamma = FWHM
// Xmax, Xmin, dX = the two bounds and the step

Spectra::Spectra(int FileNumber, char** Files, int * LineNumber, double Gamma, double Xmin, double Xmax, double dX)
{
  int number = 0;
  int N = (Xmax - Xmin) / dX;
  double * Energy = new double [N];
  double * Absorption = new double [N]; double tmp1; double tmp2 = 0.0; double g = Gamma * Gamma / 4;

  for (int i = 0; i < N; ++i)
    {
      Energy[i] = Xmin + dX * i;
      Absorption[i] = 0.0;
    }

  for (int i = 0; i < FileNumber; ++i)
    {
      int n = LineNumber[i];
      number += n;
      double** tmp;
      tmp = new double* [2]; tmp[0] = new double [n]; tmp[1] = new double [n];
      ifstream file;
      file.open(Files[i],ios::out);
      if (!file.is_open())
        {
	  cout << "Error in open the file: " << Files[i] << "Exit now" << endl;
	  exit(0);
	}
      for (int j = 0; j < n; ++j)
	file >> tmp[0][j] >> tmp[1][j];
      file.close();

      for (int j = 0; j < N; ++j)
	{
	  tmp1 = Energy[j]; tmp2 = 0.0;
	  for (int k = 0; k < n; ++k)
	    tmp2 += tmp[1][k] / ((tmp1 - tmp[0][k]) * (tmp1 - tmp[0][k]) + g);
	  Absorption[j] += tmp2;
	}
      delete[] tmp;
      tmp = 0;
    }
  double tmp3 = Gamma/(2 * M_PI * number);
  for (int i = 0; i < N; ++i)
    Absorption[i] *= tmp3;

  this->AxeX = new RealVector(Energy, N);
  this->AxeY = new RealVector(Absorption, N);
  this->PointNumber = N;
}


// virtual destructor
//

Spectra::~Spectra()
{
  if ((AxeX != NULL) && (AxeY != NULL))
    {
      delete AxeX;
      delete AxeY;
    }
}

// virtual method to write the spectrum in a file in ASCII mode
//
// fileName = name of the file where the spectrum will be stored
// return = true if no error occurs

bool Spectra::WriteSpectra(char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  File.precision(14);
  for (int i = 0; i < this->PointNumber; ++i)
    {
      File << (*AxeX)[i] << '\t';
      File << (*AxeY)[i] << '\n';
    }
  File.close();
  return true;
}

// Add two string and an integer to a string
//
// s1, s2 = strings of component, n = integer
// s = destination string

void AddString(char* s, char* s1, int n, char* s2)
{
  strcpy(s, s1);
  char  temp[80];
  sprintf (temp, "%d", n);
  strcat (s, temp);
  strcat(s, s2);
}

// read spectrum raw data from a file
// 
// filename = name of  the file that conatins the spectrum (with optional relative/absolute path)
// energies = array where energy values will be stored
// nbrValues = number of energy values to retrieve from the file
// return value = true if no error occured

bool Spectra::ReadSpectrum(char* filename, double* energies, int nbrValues)
{
  ifstream File;
  File.open(filename, ios::out);
  if (!File.is_open())
    {
      cout << "error while opening file : " << filename << endl;
      return false;
    }
  for (int j = 0; j < nbrValues; ++j)
    {
      if (File.tellg() < 0)
	{
	  cout << filename <<  " has to few eigenvalues" << endl;
	  File.close();
	  return false;
	}
      else
	File >> energies[j];
    }
  File.close();
  return true;  
}
