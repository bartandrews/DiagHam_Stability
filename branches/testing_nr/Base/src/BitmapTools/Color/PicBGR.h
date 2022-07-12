////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of BGR Encoding Color for Picture Recording             //
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


#ifndef PICBGR_H
#define PICBGR_H

#include "config.h"
#include <fstream>
#include "BitmapTools/Color/RGB.h"
#include "BitmapTools/Color/PicRGB.h"


using std::ofstream;
using std::ifstream;


class PicBGR
{

public:

  // intensity in standard BGR
  unsigned char Red;
  unsigned char Green;
  unsigned char Blue;

  //constructors
  PicBGR ();
  PicBGR (const unsigned char& RedInt, const unsigned char& GreenInt, const unsigned char& BlueInt);
  PicBGR (const PicBGR& Col);
  PicBGR (const PicRGB& Col);
  PicBGR (const RGB& Col);
  
  // assignement
  PicBGR& operator = (const PicBGR& Col);
  PicBGR& operator = (const PicRGB& Col);
  PicBGR& operator = (const RGB& Col);

  //destructor
  ~PicBGR () {};
  
  //compare two colors
  friend bool operator == (PicBGR& Col1,  PicBGR& Col2);
  friend bool operator != (PicBGR& Col1,  PicBGR& Col2);
  
  // storage in a file
  friend ofstream& operator << (ofstream& FileOut, PicBGR& Col);
  friend ifstream& operator >> (ifstream& FileOut, PicBGR& Col);

};

#endif
