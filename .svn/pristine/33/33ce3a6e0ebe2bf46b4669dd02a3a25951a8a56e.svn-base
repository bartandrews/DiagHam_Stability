////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                          class for average spectra                         //
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


#include "Tools/Spectra/AverageSpectra.h"
#include "Vector/RealVector.h"

#include <fstream>
#include <math.h>

using std::ifstream;
using std::ofstream;
using std::ios;
using std::cout;
using std::endl;


// default constructor
//

AverageSpectra::AverageSpectra()
{
  AxeX = 0;
  AxeY = 0;
  MeanXX = 0;
  MeanY = 0;
  MeanYY = 0;
  MeanZ = 0;
  MeanZZ = 0;
  PointNumber = 0;
}

// constructor from two files: energy and state
// 
// StateFile = file containing all states
// EnergyFile = file containing all energy
// Number = number of states
// M & N = number of cells in x and y directions
// a & b = lattice constants in x and y directions

AverageSpectra::AverageSpectra(char* StateFile, char* EnergyFile ,int Number, int M, int N, double a, double b)
{
  double Lx = M * a;
  double Ly = N * b;
  PointNumber = Number;
  double** state = new double* [N];

  ifstream Energy(EnergyFile); ifstream State(StateFile);

  AxeX = new RealVector(Number);
  AxeY = new RealVector(Number);
  MeanXX = new RealVector(Number);
  MeanY = new RealVector(Number);
  MeanYY = new RealVector(Number);

  double tmp1 = 0.0; double tmp2= 0.0; double tmp = 0.0;

  for (int i = 0; i < Number; ++i)
    {      
      Energy >> (*AxeX)[i];
      for (int n = 0; n < N; ++n)
	{
	  state[n] = new double [M];
	  for (int m = 0; m < M; ++m)
	    State >> state[n][m];
	}
      tmp1 = 0.0; tmp2 = 0.0;
      for (int m1 = 0; m1 < M; ++m1)
	{	  
	  for (int m2 = 0; m2 < m1; ++m2)
	    {
	      tmp = 0.0;
	      for (int n = 0; n < N; ++n)
		tmp += state[n][m1] * state[n][m2];
	      tmp1 += tmp * this->moyenX(m1 + 1, m2 + 1, Lx);
	      tmp2 += tmp * this->moyenXX(m1 + 1, m2 + 1, Lx);
	    }
	  tmp =0.0;
	  for (int n = 0; n < N; ++n)
	    tmp += state[n][m1] * state[n][m1];
	  tmp1 += (tmp * this->moyenX(m1 + 1, m1 + 1, Lx) / 2);
	  tmp2 += (tmp * this->moyenXX(m1 + 1, m1 + 1, Lx) / 2);	  
	}
      (*AxeY)[i] = tmp1 * 2;
      (*MeanXX)[i] = tmp2 * 2;
      tmp1 = 0.0; tmp2 = 0.0;
      for (int n1 = 0; n1 < N; ++n1)
	{	  
	  for (int n2 = 0; n2 < n1; ++n2)
	    {
	      tmp = 0.0;
	      for (int m = 0; m < M; ++m)
		tmp += state[n1][m] * state[n2][m];
	      tmp1 += tmp * this->moyenX(n1 + 1, n2 + 1, Ly);
	      tmp2 += tmp * this->moyenXX(n1 + 1, n2 + 1, Ly);
	    }
	  tmp =0.0;
	  for (int m = 0; m < M; ++m)
	    tmp += state[n1][m] * state[n1][m];
	  tmp1 += (tmp * this->moyenX(n1 + 1, n1 + 1, Ly) / 2);
	  tmp2 += (tmp * this->moyenXX(n1 + 1, n1 + 1, Ly) / 2);
	}
      (*MeanY)[i] = tmp1 * 2;
      (*MeanYY)[i] = tmp2 * 2;
    }
  Energy.close(); State.close();
}

// constructor from an energy file and a set of state files
// 
// StateFile = a set containing all states
// EnergyFile = file containing all energy
// Number: = number of states
// M, N and H = number of cells in x, y and z directions
// a, b and c = lattice constants in x, y and z directions

AverageSpectra::AverageSpectra(char** StateFile, char* EnergyFile ,int Number, int M, int N, int H, double a, double b, double c)
{
  double Lx = M * a;
  double Ly = N * b;
  double Lz = H * c;
  PointNumber = Number;
  AxeX = new RealVector(Number);
  AxeY = new RealVector(Number);
  MeanXX = new RealVector(Number);
  MeanY = new RealVector(Number);
  MeanYY = new RealVector(Number);
  MeanZ = new RealVector(Number);
  MeanZZ = new RealVector(Number);
  double tmp1 = 0.0; double tmp2 = 0.0; double tmp = 0.0;
  double*** state = new double** [M];
  ifstream Energy(EnergyFile);
  for (int i = 0; i < Number; ++i)
    {        
      Energy >> (*AxeX)[i]; ifstream State(StateFile[i]); 
      for (int m = 0; m < M; ++m)
	{
	  state[m] = new double* [N];
	  for (int n = 0; n < N; ++n)
	    {
	      state[m][n] = new double [H];
	      for (int h = 0; h < H; ++h)
		State >> state[m][n][h];
	    }
	}      
      tmp1 = 0.0; tmp2 = 0.0;
      for (int m1 = 0; m1 < M; ++m1)
	{	  
	  for (int m2 = 0; m2 < m1; ++m2)
	    {
	      tmp = 0.0;
	      for (int n = 0; n < N; ++n)
		for (int h = 0; h < H; ++h)
		  tmp += state[m1][n][h] * state[m2][n][h];
	      tmp1 += tmp * this->moyenX(m1 + 1, m2 + 1, Lx);
	      tmp2 += tmp * this->moyenXX(m1 + 1, m2 + 1, Lx);
	    }
	  tmp =0.0;
	  for (int n = 0; n < N; ++n)
	    for (int h = 0; h < H; ++h)
	      tmp += state[m1][n][h] * state[m1][n][h];
	  tmp1 += (tmp * this->moyenX(m1 + 1, m1 + 1, Lx) / 2);
	  tmp2 += (tmp * this->moyenXX(m1 + 1, m1 + 1, Lx) / 2);	  
	}     
      (*AxeY)[i] = tmp1 * 2;
      (*MeanXX)[i] = tmp2 * 2;
      tmp1 = 0.0; tmp2 = 0.0;
      for (int n1 = 0; n1 < N; ++n1)
	{	  
	  for (int n2 = 0; n2 < n1; ++n2)
	    {
	      tmp = 0.0;
	      for (int m = 0; m < M; ++m)
		for (int h = 0; h < H; ++h)
		  tmp += state[m][n1][h] * state[m][n2][h];
	      tmp1 += tmp * this->moyenX(n1 + 1, n2 + 1, Ly);
	      tmp2 += tmp * this->moyenXX(n1 + 1, n2 + 1, Ly);
	    }
	  tmp =0.0;
	  for (int m = 0; m < M; ++m)
	    for (int h = 0; h < H; ++h)
	      tmp += state[m][n1][h] * state[m][n1][h];
	  tmp1 += (tmp * this->moyenX(n1 + 1, n1 + 1, Ly) / 2);
	  tmp2 += (tmp * this->moyenXX(n1 + 1, n1 + 1, Ly) / 2);
	}      
      (*MeanY)[i] = tmp1 * 2;
      (*MeanYY)[i] = tmp2 * 2;
      tmp1 = 0.0; tmp2 = 0.0;
      for (int h1 = 0; h1 < H; ++h1)
	{	  
	  for (int h2 = 0; h2 < h1; ++h2)
	    {
	      tmp = 0.0;
	      for (int m = 0; m < M; ++m)
		for (int n = 0; n < N; ++n)
		  tmp += state[m][n][h1] * state[m][n][h2];
	      tmp1 += tmp * this->moyenX(h1 + 1, h2 + 1, Lz);
	      tmp2 += tmp * this->moyenXX(h1 + 1, h2 + 1, Lz);
	    }
	  tmp =0.0;
	  for (int m = 0; m < M; ++m)
	    for (int n = 0; n < N; ++n)
	      tmp += state[m][n][h1] * state[m][n][h1];
	  tmp1 += (tmp * this->moyenX(h1 + 1, h1 + 1, Lz) / 2);
	  tmp2 += (tmp * this->moyenXX(h1 + 1, h1 + 1, Lz) / 2);
	}
      (*MeanZ)[i] = tmp1 * 2;
      (*MeanZZ)[i] = tmp2 * 2; State.close();      
    }
  Energy.close();
}

// write data in a file with 5 columns: energy, <x>, <x²>, <y>, <y²>
//
// fileName = name of the file where the spectrum will be stored
// return = true if no error occurs

bool AverageSpectra::WriteXY(char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out | ios::app);
  File.precision(14);
  for (int i = 0; i < this->PointNumber; ++i)
    {
      File << (*AxeX)[i] << '\t';
      File << (*AxeY)[i] << '\t';
      File << (*MeanXX)[i] << '\t';
      File << (*MeanY)[i] << '\t';   
      File << (*MeanYY)[i] << '\n';
    }
  File.close();
  return true;
}


// write data in a file with 7 columns: energy, <x>, <x²>, <y>, <y²>, <z> and <z²>
//
// fileName = name of the file where the spectrum will be stored
// return = true if no error occurs

bool AverageSpectra::WriteXYZ(char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out | ios::app);
  File.precision(14);
  for (int i = 0; i < this->PointNumber; ++i)
    {
      File << (*AxeX)[i] << '\t';
      File << (*AxeY)[i] << '\t';
      File << (*MeanXX)[i] << '\t';
      File << (*MeanY)[i] << '\t';   
      File << (*MeanYY)[i] << '\t';
      File << (*MeanZ)[i] << '\t';   
      File << (*MeanZZ)[i] << '\n';
    }
  File.close();
  return true;
}

// write data in a file with 2 columns: <x> and var(x)
//
// fileName = name of the file to store data
// return = true if no error occurs

bool AverageSpectra::WriteVarMean(char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out | ios::app);
  File.precision(14);
  for (int i = 0; i < this->PointNumber; ++i)
    {
      File << (*AxeY)[i] << '\t';
      File << sqrt((*MeanXX)[i] - (((*AxeY)[i]) * ((*AxeY)[i]))) << '\n';
    }
  File.close();
  return true;
}

