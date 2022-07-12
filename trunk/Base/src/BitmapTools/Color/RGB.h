////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of RGB Encoding Color                          //
//                                                                            //
//                        last modification : 07/06/2000                      //
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


#ifndef RGB_H
#define RGB_H


#define EPSILONCOLOR 0.002


class HLS;
class PicRGB;

class RGB
{

public:
  
  // intensity in standard RGB
  double RedIntensity;
  double GreenIntensity;
  double BlueIntensity;

  //constructors
  RGB ();
  RGB (const double& Int);
  RGB (const double& RedInt, const double& GreenInt, const double& BlueInt);
  RGB (const RGB& Col);
  RGB (const PicRGB& Col);
  RGB (HLS& Col);
  
  // assignement
  RGB& operator = (const RGB& Col);
  RGB& operator = (const PicRGB& Col);
  
  //destructor
  ~RGB () {};
  
  // conversions
  RGB& operator = (HLS& Col);
  
  // sum of two RGB colors
  RGB& operator += (const RGB& Col2);
  friend RGB operator + (const RGB& Col1, const RGB& Col2);
  
  // difference of two RGB colors
  RGB& operator -= (const RGB& Col2);
  friend RGB operator - (const RGB& Col1, const RGB& Col2);
  
  // product of two RGB colors
  RGB& operator *= (const RGB& Col2);
  friend RGB operator * (const RGB& Col1, const RGB& Col2);
  
  // product of a RGB color with a real
  RGB& operator *= (double& x);
  friend RGB operator * (const RGB& Col1, double x);
  friend RGB operator * (double x, const RGB& Col1);
  
  // division of a RGB color by a real
  RGB& operator /= (double& x);
  friend RGB operator / (const RGB& Col1, double& x);
  
  // truncate all wrong components of a RGB color
  friend RGB RegularizeRGBColor (const RGB& Col);

  // compare two colors versus their intensity
  bool operator < (const RGB& Col);
  bool operator > (const RGB& Col);
  
  bool operator <= (const RGB& Col);
  bool operator >= (const RGB& Col);

};

#endif
#include "config.h"
