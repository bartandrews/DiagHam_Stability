////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                  class for spectra of time resolved experiment             //
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


#ifndef TIMERESOLVEDPLSPECTRA_H
#define TIMERESOLVEDPLSPECTRA_H

#include "config.h"

#include "Tools/Spectra/Spectra.h"

class TimeResolvedPLSpectra : public Spectra
{
 protected:

  // the time where the signal decrease 10 times
  double T10;

 public:

  // constructor from a set of squared overlap files
  //
  // FileNumber = number of squared overlap files
  // FileName = name of square overlap files
  // ElectronState, HoleState = number of states in each file
  // Emin, Emax = range of energy observed
  // T, NT = time and number of time intervals
  TimeResolvedPLSpectra(int FileNumber, char** FileName, int* ElectronState, int* HoleState, double Emin, double Emax, double T, int NT);


};

#endif
