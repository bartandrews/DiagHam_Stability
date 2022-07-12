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


#ifndef SPECTRA_H
#define SPECTRA_H

#include "config.h"

#include <iostream>

class RealVector;

class Spectra
{
 protected:

  // the x-axe, can be time of energy or anything else
  RealVector* AxeX;
  
  // the y-axe
  RealVector* AxeY;

  // number of discrete points
  int PointNumber;

 public:

  // default constructor
  //
  Spectra();

  // constructor
  //
  // N = number of points
  Spectra(int N);
 
  // constructor from a set of files with 2 column data. Each peak is assimilated to a Lorentzian function with FWHM = Gamma
  //
  //
  // FileNumber = number of files, Files = name of files
  // LineNumber = number of lines, Gamma = FWHM
  // Xmax, Xmin, dX = the two bounds and the step
  Spectra(int FileNumber, char** Files, int * StateNumber, double Gamma, double Emin, double Emax, double dE);

  // virtual destructor
  //
  virtual ~Spectra();

  // virtual method to write the spectrum in a file in ASCII mode
  //
  // fileName = name of the file where the spectrum will be stored
  // return = true if no error occurs
  virtual bool WriteSpectra(char * fileName);

 protected:

  // read spectrum raw data from a file
  // 
  // filename = name of  the file that conatins the spectrum (with optional relative/absolute path)
  // energies = array where energy values will be stored
  // nbrValues = number of energy values to retrieve from the file
  // return value = true if no error occured
  bool ReadSpectrum(char* filename, double* energies, int nbrValues);

};

// Add two string and an integer to a string
//
// s1, s2 = strings of component, n = integer
// s = destination string
void AddString(char* s, char* s1, int n, char* s2);

#endif
