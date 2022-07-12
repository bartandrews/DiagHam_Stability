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


#include "BitmapTools/Color/HLS.h"


//constructors

HLS::HLS ()
{
  this->Hue = 0;
  this->Lightness = 1.0;
  this->Saturation = 0;
}

HLS::HLS (const double& Hue, const double& Lightness, const double& Saturation)
{
  this->Hue = Hue;
  this->Lightness = Lightness;
  this->Saturation = Saturation;
}

HLS::HLS (const HLS& Col)
{
  this->Hue = Col.Hue;
  this->Lightness = Col.Lightness;
  this->Saturation = Col.Saturation;
}

HLS::HLS (const RGB& Col)
{
  double Max = Col.RedIntensity;
  double Min = Col.RedIntensity;
  if (Col.GreenIntensity > Max)
    Max = Col.GreenIntensity;
  if (Col.BlueIntensity > Max)
    Max = Col.BlueIntensity;
  if (Col.GreenIntensity < Min)
    Min = Col.GreenIntensity;
  if (Col.BlueIntensity < Min)
    Min = Col.BlueIntensity;
  if (Max > 1)
    Max = 1;
  if (Min < 0)
    Min = 0;
  this->Lightness = (Max + Min)/2;
  if ((Max >= Min - EPSILONCOLOR) && (Max <= Min + EPSILONCOLOR))
    {
      this->Saturation = 0;
      this->Hue = 0;
    }
  else
    {
      if (this->Lightness <= 0.5)
	this->Saturation = (Max - Min)/(Max + Min);
      else
	this->Saturation = (Max - Min)/(2- (Max + Min));
      if (Col.RedIntensity >= Max)
      	this->Hue = (Col.GreenIntensity - Col.BlueIntensity) / (Max - Min);
      if (Col.GreenIntensity >= Max)
      	this->Hue = 2 + (Col.BlueIntensity - Col.RedIntensity) / (Max - Min);
      if (Col.BlueIntensity >= Max)
      	this->Hue = 4 + (Col.RedIntensity - Col.GreenIntensity) / (Max - Min);
      this->Hue *= 60;
      if (this->Hue < 0)
      	this->Hue += 360;

    }
}

// assignement

HLS& HLS::operator = (const HLS& Col)
{
  this->Hue = Col.Hue;
  this->Lightness = Col.Lightness;
  this->Saturation = Col.Saturation;
  return *this;
}

// conversions

HLS& HLS::operator = (const RGB& Col)
{
  double Max = Col.RedIntensity;
  double Min = Col.RedIntensity;
  if (Col.GreenIntensity > Max)
    Max = Col.GreenIntensity;
  if (Col.BlueIntensity > Max)
    Max = Col.BlueIntensity;
  if (Col.GreenIntensity < Min)
    Min = Col.GreenIntensity;
  if (Col.BlueIntensity < Min)
    Min = Col.BlueIntensity;
  if (Max > 1)
    Max = 1;
  if (Min < 0)
    Min = 0;
  this->Lightness = (Max + Min)/2;
  if ((Max >= Min - EPSILONCOLOR) && (Max <= Min + EPSILONCOLOR))
    {
      this->Saturation = 0;
      this->Hue = 0;
    }
  else
    {
      if (this->Lightness <= 0.5)
	this->Saturation = (Max - Min)/(Max + Min);
      else
	this->Saturation = (Max - Min)/(2- Max + Min);
      if (Col.RedIntensity >= Max)
      	this->Hue = (Col.GreenIntensity - Col.BlueIntensity) / (Max - Min);
      if (Col.GreenIntensity >= Max)
      	this->Hue = 2 + (Col.BlueIntensity - Col.RedIntensity) / (Max - Min);
      if (Col.BlueIntensity >= Max)
      	this->Hue = 4 + (Col.RedIntensity - Col.GreenIntensity) / (Max - Min);
      this->Hue *= 60;
      if (this->Hue < 0)
      	this->Hue += 360;      
    }
  return *this;
}

// product of a HLS color with a real

HLS& HLS::operator *= (double& x)
{
  double Max;
  if (this->Lightness <= 0.5)
    Max = this->Lightness * (1 + this->Saturation);
  else
    Max = this->Lightness + this->Saturation * (1 - this->Lightness);
  double Min = 2 * this->Lightness - Max;
  this->Lightness *= x;
  if (this->Lightness <= 0.5)
    this->Saturation = (Max - Min)/(Max + Min);
  else
    this->Saturation = (Max - Min)/(2- (Max + Min));
  return *this;
}

HLS operator * (HLS& Col1, double& x)
{
  double Max;
  if (Col1.Lightness <= 0.5)
    Max = Col1.Lightness * (1 + Col1.Saturation);
  else
    Max = Col1.Lightness + Col1.Saturation  * (1 - Col1.Lightness);
  double Min = 2 * Col1.Lightness - Max;
  return HLS(Col1.Hue, Col1.Lightness * x, Col1.Saturation);
}

HLS operator * (double& x, HLS& Col1)
{
  double Max;
  if (Col1.Lightness <= 0.5)
    Max = Col1.Lightness * (1 + Col1.Saturation);
  else
    Max = Col1.Lightness + Col1.Saturation  * (1 - Col1.Lightness);
  double Min = 2 * Col1.Lightness - Max;
  return HLS(Col1.Hue, Col1.Lightness * x, Col1.Saturation);
}
