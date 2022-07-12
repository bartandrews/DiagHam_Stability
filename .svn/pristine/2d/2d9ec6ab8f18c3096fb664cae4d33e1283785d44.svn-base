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


#include "Color/PicRGBE.h"
#include "Color/PicRGB.h"

#include <math.h>


// default constructor
//

PicRGBE::PicRGBE ()
{
  this->Red = 0;
  this->Green = 0;
  this->Blue = 0;
  this->Exponent = 0;
}

// basic constructor
//
// red = red component 
// green = green component 
// blue = blue component 
// exponent = global exponent for all color components

PicRGBE::PicRGBE (const unsigned char& red, const unsigned char& green, const unsigned char& blue, 
		  const unsigned char& exponent)
{
  this->Red = red;
  this->Green = green;
  this->Blue = blue;
  this->Exponent = exponent;
}

// copy constructor
//
// color = color to copy

PicRGBE::PicRGBE (const PicRGBE& color)
{
  this->Red = color.Red;
  this->Green = color.Green;
  this->Blue = color.Blue;
  this->Exponent = color.Exponent;
}

  
// copy constructor from a PicRGB (assuming no global exponent factor)
//
// color = color to copy

PicRGBE::PicRGBE (const PicRGB& color)
{
  this->Red = color.Red;
  this->Green = color.Green;
  this->Blue = color.Blue;
  this->Exponent = 128;
}

// copy constructor from a RGB value
//
// color = color to copy  

PicRGBE::PicRGBE (const RGB& color)
{
  double Max = color.RedIntensity;
  if (Max < color.GreenIntensity)
    Max = color.GreenIntensity;
  if (Max < color.BlueIntensity)
    Max = color.BlueIntensity;
  int TmpExponent;
  Max = frexp (Max, &TmpExponent) * 255.9999999999 / Max;
  this->Red = (unsigned char) (Max * color.RedIntensity);
  this->Green = (unsigned char) (Max * color.GreenIntensity);
  this->Blue = (unsigned char) (Max * color.BlueIntensity);
  this->Exponent = (unsigned char)  (128 + TmpExponent);
}

// destructor
//

PicRGBE::~PicRGBE ()
{
}

// assignement 
// 
// color = color to copy
// return value = reference on current PicRGBE

PicRGBE& PicRGBE::operator = (const PicRGBE& color)
{
  this->Red = color.Red;
  this->Green = color.Green;
  this->Blue = color.Blue;
  this->Exponent = color.Exponent;
  return *this;
}

// assignement 
// 
// color = color to copy
// return value = reference on current PicRGBE

PicRGBE& PicRGBE::operator = (const PicRGB& color)
{
  this->Red = color.Red;
  this->Green = color.Green;
  this->Blue = color.Blue;
  this->Exponent = 128;
  return *this;
}

// assignement 
// 
// color = color to copy
// return value = reference on current PicRGBE

PicRGBE& PicRGBE::operator = (const RGB& color)
{
  double Max = color.RedIntensity;
  if (Max < color.GreenIntensity)
    Max = color.GreenIntensity;
  if (Max < color.BlueIntensity)
    Max = color.BlueIntensity;
  int TmpExponent;
  Max = frexp (Max, &TmpExponent) * 255.9999999999 / Max;
  this->Red = (unsigned char) (Max * color.RedIntensity);
  this->Green = (unsigned char) (Max * color.GreenIntensity);
  this->Blue = (unsigned char) (Max * color.BlueIntensity);
  this->Exponent = (unsigned char)  (128 + TmpExponent);
  return *this;
}

// get RGB value associated to the current PicRGBE value
// 
// return value = corresponding RGB value

RGB PicRGBE::GetRGB ()
{
  if (this->Exponent == 0)
    {
      RGB TmpColor (0.0, 0.0, 0.0);
      return TmpColor;
    }
  double TmpExponent = ldexp(1.0, ((int) this->Exponent)-(136));
  RGB TmpColor ((((double) this->Red) + 0.5) * TmpExponent, (((double) this->Green) + 0.5) * TmpExponent, 
		(((double) this->Blue) + 0.5) * TmpExponent);
  return TmpColor;
}
  
// get RGB value associated to the current PicRGBE value
// 
// color = reference on RGB value where result has to be stored
// return value = reference on RGB value

RGB& PicRGBE::GetRGB (RGB& color)
{
  if (this->Exponent == 0)
    {
      color.RedIntensity = 0.0;
      color.GreenIntensity = 0.0;
      color.BlueIntensity = 0.0;
      return color;
    }
  double TmpExponent = ldexp(1.0, ((int) this->Exponent)-(136));
  color.RedIntensity = (((double) this->Red) + 0.5) * TmpExponent;
  color.GreenIntensity = (((double) this->Green) + 0.5) * TmpExponent;
  color.BlueIntensity = (((double) this->Blue) + 0.5) * TmpExponent;
  return color;
}
  
// compare if two colors are identical
//
// color1 = reference on color to the left hand side of the operator
// color1 = reference on color to the right hand side of the operator
// return value = true if colors are identical

bool operator == (PicRGBE& color1,  PicRGBE& color2)
{
  return ((color1.Red == color2.Red) && (color1.Green == color2.Green) && (color1.Blue == color2.Blue) && 
	  (color1.Exponent == color2.Exponent));
}

// compare if two colors are different
//
// color1 = reference on color to the left hand side of the operator
// color1 = reference on color to the right hand side of the operator
// return value = true if colors are different

bool operator != (PicRGBE& color1,  PicRGBE& color2)
{
  return ((color1.Red != color2.Red) || (color1.Green != color2.Green) || (color1.Blue != color2.Blue) || 
	  (color1.Exponent != color2.Exponent));
}

// store PicRGBE value in a file
//
// file = reference on output file stream
// color = reference on color to store
// return value = reference on output file stream

ofstream& operator << (ofstream& file, PicRGBE& color)
{
  file.write((char*) &(color.Red), sizeof(unsigned char));
  file.write((char*) &(color.Green), sizeof(unsigned char));
  file.write((char*) &(color.Blue), sizeof(unsigned char));
  file.write((char*) &(color.Exponent), sizeof(unsigned char));
  return file;
}

// read PicRGBE from a file a file
//
// file = reference on input file stream
// color = reference on color where read value has to be stored
// return value = reference on intput file stream

ifstream& operator >> (ifstream& file, PicRGBE& color)
{
  file.read((char*) &(color.Red), sizeof(unsigned char));
  file.read((char*) &(color.Green), sizeof(unsigned char));
  file.read((char*) &(color.Blue), sizeof(unsigned char));
  file.read((char*) &(color.Exponent), sizeof(unsigned char));
  return file;
}

