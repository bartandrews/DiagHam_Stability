////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of HLS Encoding Color                          //
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


#ifndef HLS_H
#define HLS_H


#include "config.h"
#include "BitmapTools/Color/RGB.h"


class HLS
{

public:
  
  // intensity in standard HLS
  double Hue;
  double Lightness;
  double Saturation;
  
  //constructors
  HLS ();
  HLS (const double& Hue, const double& Lightness, const double& Saturation);
  HLS (const HLS& Col);
  HLS (const RGB& Col);
  
  // assignement
  HLS& operator = (const HLS& Col);
  
  //destructor
  ~HLS () {};
  
  // conversions
  HLS& operator = (const RGB& Col);
  
  // product of a HLS color with a real
  HLS& operator *= (double& x);
  friend HLS operator * (HLS& Col1, double& x);
  friend HLS operator * (double& x, HLS& Col1);
  
};

#endif
