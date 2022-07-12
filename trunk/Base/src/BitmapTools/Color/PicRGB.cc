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


#include "BitmapTools/Color/PicRGB.h"
#include "BitmapTools/Color/PicBGR.h"


//constructors

PicRGB::PicRGB ()
{
  this->Red = 255;
  this->Green = 255;
  this->Blue = 255;
}

PicRGB::PicRGB (const unsigned char& RedInt, const unsigned char& GreenInt, const unsigned char& BlueInt)
{
  this->Red = RedInt;
  this->Green = GreenInt;
  this->Blue = BlueInt;
}

PicRGB::PicRGB (const PicRGB& Col)
{
  this->Red = Col.Red;
  this->Green = Col.Green;
  this->Blue = Col.Blue;
}

PicRGB::PicRGB (const PicBGR& Col)
{
  this->Red = Col.Red;
  this->Green = Col.Green;
  this->Blue = Col.Blue;
}

PicRGB::PicRGB (const RGB& Col)
{
  if (Col.RedIntensity >= 1)
    this->Red = (unsigned char) 255;
  else
    this->Red = (unsigned char) (255 * Col.RedIntensity);
  if (Col.GreenIntensity >= 1)
    this->Green = (unsigned char) 255;
  else
    this->Green = (unsigned char) (255 * Col.GreenIntensity);
  if (Col.BlueIntensity >= 1)
    this->Blue = (unsigned char) 255;
  else
    this->Blue = (unsigned char) (255 * Col.BlueIntensity);
}

// assignement

PicRGB& PicRGB::operator = (const PicRGB& Col)
{
  this->Red = Col.Red;
  this->Green = Col.Green;
  this->Blue = Col.Blue;
  return *this;
}

PicRGB& PicRGB::operator = (const PicBGR& Col)
{
  this->Red = Col.Red;
  this->Green = Col.Green;
  this->Blue = Col.Blue;
  return *this;
}

PicRGB& PicRGB::operator = (const RGB& Col)
{
  if (Col.RedIntensity >= 1)
    this->Red = (unsigned char) 255;
  else
    this->Red = (unsigned char) (255 * Col.RedIntensity);
  if (Col.GreenIntensity >= 1)
    this->Green = (unsigned char) 255;
  else
    this->Green = (unsigned char) (255 * Col.GreenIntensity);
  if (Col.BlueIntensity >= 1)
    this->Blue = (unsigned char) 255;
  else
    this->Blue = (unsigned char) (255 * Col.BlueIntensity);
  return *this;
}

//compare two colors

bool operator == (PicRGB& Col1,  PicRGB& Col2)
{
  return ((Col1.Red == Col2.Red) && (Col1.Green == Col2.Green) && (Col1.Blue == Col2.Blue));
}

bool operator != (PicRGB& Col1,  PicRGB& Col2)
{
  return ((Col1.Red != Col2.Red) || (Col1.Green != Col2.Green) || (Col1.Blue != Col2.Blue));
}

// storage in a file

ofstream& operator << (ofstream& FileOut, PicRGB& Col)
{
  FileOut.put((char&) Col.Red);
  FileOut.put((char&) Col.Green);
  FileOut.put((char&) Col.Blue);
  return FileOut;
}

ifstream& operator >> (ifstream& FileOut, PicRGB& Col)
{
  FileOut.get((char&) Col.Red);
  FileOut.get((char&) Col.Green);
  FileOut.get((char&) Col.Blue);
  return FileOut;
}

