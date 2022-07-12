////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of RGB Encoding Color for Picture Recording               //
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


#ifndef PICRGB_H
#define PICRGB_H


#include "config.h"
#include <fstream>
#include "Color/RGB.h"


using std::ofstream;
using std::ifstream;


class PicBGR;


class PicRGB
{

public:
  
  // intensity in standard RGB
  unsigned char Red;
  unsigned char Green;
  unsigned char Blue;
  
  //constructors
  PicRGB ();
  PicRGB (const unsigned char& RedInt, const unsigned char& GreenInt, const unsigned char& BlueInt);
  PicRGB (const PicRGB& Col);
  PicRGB (const PicBGR& Col);
  PicRGB (const RGB& Col);
  
  // assignement
  PicRGB& operator = (const PicRGB& Col);
  PicRGB& operator = (const PicBGR& Col);
  PicRGB& operator = (const RGB& Col);
  
  //destructor
  ~PicRGB () {};
  
  //compare two colors
  friend bool operator == (PicRGB& Col1,  PicRGB& Col2);
  friend bool operator != (PicRGB& Col1,  PicRGB& Col2);
  
  // storage in a file
  friend ofstream& operator << (ofstream& FileOut, PicRGB& Col);
  friend ifstream& operator >> (ifstream& FileOut, PicRGB& Col);
  
};

#endif
