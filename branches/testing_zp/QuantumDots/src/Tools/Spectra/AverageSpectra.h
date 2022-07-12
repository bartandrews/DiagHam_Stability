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


#ifndef AVERAGESPECTRA_H
#define AVERAGESPECTRA_H

#include "config.h"

#include "Tools/Spectra/Spectra.h"

#include <math.h>

class RealVector;

class AverageSpectra : public Spectra
{
 protected:

  // average and square average values in X, Y and Z directions
  RealVector* MeanXX;
  RealVector* MeanY;
  RealVector* MeanYY;
  RealVector* MeanZ;
  RealVector* MeanZZ;

 public:

  // default constructor
  //
  AverageSpectra();
  
  // constructor from two files: energy and state
  // 
  // StateFile = file containing all states
  // EnergyFile = file containing all energy
  // Number = number of states
  // M & N = number of cells in x and y directions
  // a & b = lattice constants in x and y directions
  AverageSpectra(char* StateFile, char* EnergyFile ,int Number, int M, int N, double a, double b);

  // constructor from an energy file and a set of state files
  // 
  // StateFile = a set containing all states
  // EnergyFile = file containing all energy
  // Number: = number of states
  // M, N and H = number of cells in x, y and z directions
  // a, b and c = lattice constants in x, y and z directions
  AverageSpectra(char** StateFile, char* EnergyFile ,int Number, int M, int N, int H, double a, double b, double c);

  // write data in a file with 5 columns: energy, <x>, <x²>, <y>, <y²>
  //
  // fileName = name of the file where the spectrum will be stored
  // return = true if no error occurs
  bool WriteXY(char* fileName);

  // write data in a file with 7 columns: energy, <x>, <x²>, <y>, <y²>, <z> and <z²>
  //
  // fileName = name of the file where the spectrum will be stored
  // return = true if no error occurs
  bool WriteXYZ(char* fileName);

  // write data in a file with 2 columns: <x> and var(x)
  //
  // fileName = name of the file to store data
  // return = true if no error occurs
  bool WriteVarMean(char* FileName);

  // calculate the term: <m|x|n> where m, n are sinus functions
  //
  // m, n = indices of sinus functions
  // L = length of the sample
  double moyenX(int m, int n, double L);

  // calculate the term: <m|x²|n> where m, n are sinus functions
  //
  // m, n = indices of sinus functions
  // L = length of the sample
  double moyenXX(int m, int n, double L);

};

// calculate the term: <m|x|n> where m, n are sinus functions
//
// m, n = indices of sinus functions
// L = length of the sample

inline double AverageSpectra::moyenX(int m, int n, double L)
{
  double tmp = M_PI*(m*m - n*n);
  if (((m + n)%2) != 0)
    return -(8*m)*(n*L)/(tmp * tmp);
  else
    { 
      if (m == n) 
	return L/2;
      else 
	return 0;
    }
}

// calculate the term: <m|x²|n> where m, n are sinus functions
//
// m, n = indices of sinus functions
// L = length of the sample

inline double AverageSpectra::moyenXX(int m, int n, double L)
{
  double tmp = M_PI*(m*m - n*n);
  if (m != n)
    {
      if (((m + n)%2) == 0)
        return 8*m*n*L*L/(tmp * tmp);
      else
        return -8*m*n*L*L/(tmp * tmp);
    }
  else
    return L*L*(2*M_PI*M_PI*m*m-3)/(6*M_PI*M_PI*m*m);
}

#endif
