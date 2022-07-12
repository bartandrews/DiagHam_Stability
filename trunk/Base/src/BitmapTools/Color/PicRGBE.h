////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of RGBE Encoding Color for Picture Recording              //
//                                                                            //
//                        last modification : 02/06/2002                      //
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


#ifndef PICRGBE_H
#define PICRGBE_H


#include "config.h"
#include <fstream>
#include "BitmapTools/Color/RGB.h"


using std::ofstream;
using std::ifstream;


class PicRGB;


class PicRGBE
{

public:
  
  // red intensity using 8 bits coding
  unsigned char Red;
  // green intensity using 8 bits coding
  unsigned char Green;
  // blue intensity using 8 bits coding
  unsigned char Blue;
  // global exponent factor for the three color components 
  unsigned char Exponent;
  
  // default constructor
  //
  PicRGBE ();

  // basic constructor
  //
  // red = red component 
  // green = green component 
  // blue = blue component 
  // exponent = global exponent for all color components
  PicRGBE (const unsigned char& red, const unsigned char& green, const unsigned char& blue, 
	   const unsigned char& exponent);

  // copy constructor
  //
  // color = color to copy
  PicRGBE (const PicRGBE& color);

  
  // copy constructor from a PicRGB (assuming no global exponent factor)
  //
  // color = color to copy
  PicRGBE (const PicRGB& color);

  // copy constructor from a RGB value
  //
  // color = color to copy  
  PicRGBE (const RGB& color);
  
  // destructor
  //
  ~PicRGBE ();
  
  // assignement 
  // 
  // color = color to copy
  // return value = reference on current PicRGBE
  PicRGBE& operator = (const PicRGBE& color);

  // assignement 
  // 
  // color = color to copy
  // return value = reference on current PicRGBE
  PicRGBE& operator = (const PicRGB& color);

  // assignement 
  // 
  // color = color to copy
  // return value = reference on current PicRGBE
  PicRGBE& operator = (const RGB& color);
  
  // get RGB value associated to the current PicRGBE value
  // 
  // return value = corresponding RGB value
  RGB GetRGB ();
  
  // get RGB value associated to the current PicRGBE value
  // 
  // color = reference on RGB value where result has to be stored
  // return value = reference on RGB value
  RGB& GetRGB (RGB& color);
    
  // compare if two colors are identical
  //
  // color1 = reference on color to the left hand side of the operator
  // color1 = reference on color to the right hand side of the operator
  // return value = true if colors are identical
  friend bool operator == (PicRGBE& color1,  PicRGBE& color2);

  // compare if two colors are different
  //
  // color1 = reference on color to the left hand side of the operator
  // color1 = reference on color to the right hand side of the operator
  // return value = true if colors are different
  friend bool operator != (PicRGBE& color1,  PicRGBE& color2);
  
  // store PicRGBE value in a file
  //
  // file = reference on output file stream
  // color = reference on color to store
  // return value = reference on output file stream
  friend ofstream& operator << (ofstream& file, PicRGBE& color);

  // read PicRGBE from a file a file
  //
  // file = reference on input file stream
  // color = reference on color where read value has to be stored
  // return value = reference on intput file stream
  friend ifstream& operator >> (ifstream& file, PicRGBE& color);
  
};

#endif
