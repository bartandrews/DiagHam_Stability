////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                            base class for overlap                          //
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


#ifndef OVERLAPSPECTRA_H
#define OVERLAPSPECTRA_H

#include "config.h"
#include "Tools/Spectra/Spectra.h"

#include <iostream>

class OverlapSpectra : public Spectra
{
 public:

  // default constructor
  // 
  OverlapSpectra();

  // constructor from a BINARY data file which contains the dimension at each new line
  //
  // ElectronStateFile, ElectronEnergyFile, ElectronNumber = state file, energy file and number of states for electrons
  // HoleStateFile, HoleEnergyFile, HoleNumber = state file, energy file and number of states for holes 
  OverlapSpectra(char* ElectronStateFile, char* ElectronEnergyFile ,int ElectronNumber, char* HoleStateFile, char* HoleEnergyFile, int HoleNumber);

  // constructor from ASCII data file
  //
  // ElectronStateFile, ElectronEnergyFile, ElectronNumber = state file, energy file and number of states for electrons
  // HoleStateFile, HoleEnergyFile, HoleNumber = state file, energy file and number of states for holes   
  // NumberState = number of states
  OverlapSpectra(char* ElectronStateFile, char* ElectronEnergyFile ,int ElectronNumber, char* HoleStateFile, char* HoleEnergyFile, int HoleNumber, int NumberState);

  // method to write in ASCII mode with 2 columns: recombination energy and square overlap
  //
  // fileName =  name of the file where the spectrum will be stored
  // return = true if no error occurs
  bool WriteSquareOverlap(char* fileName);

  // evaluate the overlap between two functions in periodic basis
  //
  // ElectronStateFile & HoleStateFile: state files
  // Dimension: dimension of the Hilbert space
  // Real: reference to the real part of overlap
  // Imaginary: reference to the imaginary part of overlap
  void GetOverlap(char* ElectronStateFile, char* HoleStateFile, int Dimension, double& Real, double& Imaginary);
};


#endif
