////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class for calculation of Pseudopotential Coefficients  //
//                                                                            //
//                        last modification : 19/11/2007                      //
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


#ifndef ABSTRACTZDENSITYPROFILE_H
#define ABSTRACTZDENSITYPROFILE_H


#include "config.h"



class AbstractZDensityProfile {
 private:
  

 public:
  
  enum ZDensityProfileTypes
    {
      TabulatedProfile = 0x00000,
      InfiniteWell = 0x00001,
      FangHoward = 0x00002,
      InfiniteWellExc = 0x00003,
    };
  
  // virtual destructor
  virtual ~AbstractZDensityProfile();

  // get minimum and maximum value of density profile where the probability density is larger than precision
  // min = minimum value of z offset
  // max = maximum value of z offset
  // precision = requested precision
  virtual void GetSupport(double &min, double &max, double precision=1e-10) = 0;

  // evaluate the density for a given offset
  // z = offset of distribution
  virtual double GetValue(double z) = 0;

  // get type of the density profile
  virtual int GetType() = 0;

  // a static class function to return an actual DensityProfile object of some type
  static AbstractZDensityProfile* CreateZDensityProfile (char *type, double width);

  // a static class function to return the name of a type of DensityProfile
  // type = a string carrying either the number of the density profile, or a filename (if tabulated)
  static char *DensityProfileName(char *type);
  


};


#endif
