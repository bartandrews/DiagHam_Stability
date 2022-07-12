////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of BGR Encoding Color for Picture Recording               //
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


#include "BitmapTools/Color/PicBGR.h"


//constructors

PicBGR::PicBGR ()
{
  this->Red = 255;
  this->Green = 255;
  this->Blue = 255;
}

PicBGR::PicBGR (const unsigned char& RedInt, const unsigned char& GreenInt, const unsigned char& BlueInt)
{
  this->Red = RedInt;
  this->Green = GreenInt;
  this->Blue = BlueInt;
}

PicBGR::PicBGR (const PicBGR& Col)
{
  this->Red = Col.Red;
  this->Green = Col.Green;
  this->Blue = Col.Blue;
}

PicBGR::PicBGR (const PicRGB& Col)
{
  this->Red = Col.Red;
  this->Green = Col.Green;
  this->Blue = Col.Blue;
}

PicBGR::PicBGR (const RGB& Col)
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

PicBGR& PicBGR::operator = (const PicBGR& Col)
{
  this->Red = Col.Red;
  this->Green = Col.Green;
  this->Blue = Col.Blue;
  return *this;
}

PicBGR& PicBGR::operator = (const PicRGB& Col)
{
  this->Red = Col.Red;
  this->Green = Col.Green;
  this->Blue = Col.Blue;
  return *this;
}

PicBGR& PicBGR::operator = (const RGB& Col)
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

bool operator == (PicBGR& Col1,  PicBGR& Col2)
{
  return ((Col1.Red == Col2.Red) && (Col1.Green == Col2.Green) && (Col1.Blue == Col2.Blue));
}

bool operator != (PicBGR& Col1,  PicBGR& Col2)
{
  return ((Col1.Red != Col2.Red) || (Col1.Green != Col2.Green) || (Col1.Blue != Col2.Blue));
}

// storage in a file

ofstream& operator << (ofstream& FileOut, PicBGR& Col)
{
  FileOut.put((char&) Col.Blue);
  FileOut.put((char&) Col.Green);
  FileOut.put((char&) Col.Red);
  return FileOut;
}

ifstream& operator >> (ifstream& FileOut, PicBGR& Col)
{
  FileOut.get((char&) Col.Blue);
  FileOut.get((char&) Col.Green);
  FileOut.get((char&) Col.Red);
  return FileOut;
}

