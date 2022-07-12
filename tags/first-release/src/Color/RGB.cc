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


#include "Color/RGB.h"
#include "Color/HLS.h"
#include "Color/PicRGB.h"


//constructors

RGB::RGB ()
{
  this->RedIntensity = 1;
  this->GreenIntensity = 1;
  this->BlueIntensity = 1;
}

RGB::RGB (const double& Int)
{
  this->RedIntensity = Int;
  this->GreenIntensity = Int;
  this->BlueIntensity = Int;
}

RGB::RGB (const double& RedInt, const double& GreenInt, const double& BlueInt)
{
  this->RedIntensity = RedInt;
  this->GreenIntensity = GreenInt;
  this->BlueIntensity = BlueInt;
}

RGB::RGB (const RGB& Col)
{
  this->RedIntensity = Col.RedIntensity;
  this->GreenIntensity = Col.GreenIntensity;
  this->BlueIntensity = Col.BlueIntensity;
}

RGB::RGB (const PicRGB& Col)
{
  this->RedIntensity = ((double)Col.Red) / 255.0;
  this->GreenIntensity = ((double)Col.Green) / 255.0;
  this->BlueIntensity = ((double)Col.Blue) / 255.0;
}

RGB::RGB (HLS& Col)
{
  if (Col.Lightness > 1)
    Col.Lightness = 1;
  if (Col.Lightness < 0)
    Col.Lightness = 0;
  if (Col.Saturation > 1)
    Col.Saturation = 1;
  if (Col.Saturation < 0)
    Col.Saturation = 0;
  double Max;
  if (Col.Lightness <= 0.5)
    Max = Col.Lightness * (1 + Col.Saturation);
  else
    Max = Col.Lightness + Col.Saturation * (1 - Col.Lightness);
  double Min = 2 * Col.Lightness - Max;
  if (Col.Saturation == 0)
    {
      this->RedIntensity = Col.Lightness;
      this->GreenIntensity = Col.Lightness;
      this->BlueIntensity = Col.Lightness;
    }
  else
    {
      if (Col.Hue < 0)
	Col.Hue += 360;
      if (Col.Hue > 360)
	Col.Hue -= 360;
      if (Col.Hue < 60)
	{
	  this->RedIntensity = Max;
	  this->GreenIntensity = Min + (Max - Min) * Col.Hue / 60;
	  this->BlueIntensity = Min;
	}
      else
	if (Col.Hue < 120)
	  {
	    this->RedIntensity = Min + (Max - Min) * (120 - Col.Hue) / 60;
	    this->GreenIntensity = Max;
	    this->BlueIntensity = Min;
	  }
	else
	  if (Col.Hue < 180)
	    {
	      this->RedIntensity = Min;
	      this->GreenIntensity = Max;
	      this->BlueIntensity =  Min + (Max - Min) * (Col.Hue - 120) / 60;
	    }
	  else
	    if (Col.Hue < 240)
	      {
		this->RedIntensity = Min;
		this->GreenIntensity =  Min + (Max - Min) * (240 - Col.Hue) / 60;
		this->BlueIntensity = Max;
	      }
	    else
	      if (Col.Hue < 300)
		{
		  this->RedIntensity =  Min + (Max - Min) * (Col.Hue - 240) / 60;
		  this->GreenIntensity = Min;
		  this->BlueIntensity = Max;
		}
	      else
		{
		  this->RedIntensity = Max;
		  this->GreenIntensity = Min;
		  this->BlueIntensity =  Min + (Max - Min) * (360 - Col.Hue) / 60;
		}
    }
}

// conversions

RGB& RGB::operator = (HLS& Col)
{
  if (Col.Lightness > 1)
    Col.Lightness = 1;
  if (Col.Lightness < 0)
    Col.Lightness = 0;
  if (Col.Saturation > 1)
    Col.Saturation = 1;
  if (Col.Saturation < 0)
    Col.Saturation = 0;
  double Max;
  if (Col.Lightness <= 0.5)
    Max = Col.Lightness * (1 - Col.Saturation);
  else
    Max = Col.Lightness + Col.Saturation / (1 + Col.Saturation);
  double Min = 2 * Col.Lightness - Max;
  if (Col.Saturation == 0)
    {
      this->RedIntensity = Col.Lightness;
      this->GreenIntensity = Col.Lightness;
      this->BlueIntensity = Col.Lightness;
    }
  else
    {
      double Value;
      if (Col.Hue < 0)
	Col.Hue += 360;
      if (Col.Hue > 360)
	Col.Hue -= 360;
      if (Col.Hue < 60)
      	Value = Max + (Min - Max) * Col.Hue / 60;
      else
      	if (Col.Hue < 180)
	  Value = Min;
	else
	  if (Col.Hue < 240)
	    Value = Max + (Min - Max) * (240 - Col.Hue) / 60;
	  else
	    Value = Max;
      this->RedIntensity = Value + 2;
      this->GreenIntensity = Value;
      this->BlueIntensity = Value - 2;
    }
  return *this;
}

// assignement

RGB& RGB::operator = (const RGB& Col)
{
  this->RedIntensity = Col.RedIntensity;
  this->GreenIntensity = Col.GreenIntensity;
  this->BlueIntensity = Col.BlueIntensity;
  return *this;
}

RGB& RGB::operator = (const PicRGB& Col)
{
  this->RedIntensity = ((double)Col.Red) / 255;
  this->GreenIntensity = ((double)Col.Green) / 255;
  this->BlueIntensity = ((double)Col.Blue) / 255;
  return *this;
}

// sum of two colors

RGB& RGB::operator += (const RGB& Col2)
{
  this->RedIntensity += Col2.RedIntensity;
  this->GreenIntensity += Col2.GreenIntensity;
  this->BlueIntensity += Col2.BlueIntensity;
  return *this;
}

RGB operator + (const RGB& Col1, const RGB& Col2)
{
  return RGB(Col1.RedIntensity + Col2.RedIntensity, Col1.GreenIntensity + Col2.GreenIntensity, Col1.BlueIntensity + Col2.BlueIntensity);
}

// difference of two RGB colors

RGB& RGB::operator -= (const RGB& Col2)
{
  this->RedIntensity -= Col2.RedIntensity;
  this->GreenIntensity -= Col2.GreenIntensity;
  this->BlueIntensity -= Col2.BlueIntensity;
  return *this;
}

RGB operator - (const RGB& Col1, const RGB& Col2)
{
  return RGB(Col1.RedIntensity - Col2.RedIntensity, Col1.GreenIntensity - Col2.GreenIntensity, Col1.BlueIntensity - Col2.BlueIntensity);
}

// product of two colors

RGB& RGB::operator *= (const RGB& Col2)
{
  this->RedIntensity *= Col2.RedIntensity;
  this->GreenIntensity *= Col2.GreenIntensity;
  this->BlueIntensity *= Col2.BlueIntensity;
  return *this;
}

RGB operator * (const RGB& Col1, const RGB& Col2)
{
  return RGB(Col1.RedIntensity * Col2.RedIntensity, Col1.GreenIntensity * Col2.GreenIntensity, Col1.BlueIntensity * Col2.BlueIntensity);
}

// product of a RGB color with a real

RGB& RGB::operator *= (double& x)
{
  this->RedIntensity *= x;
  this->GreenIntensity *= x;
  this->BlueIntensity *= x;
  return *this;
}

RGB operator * (const RGB& Col1, double x)
{
  return RGB(Col1.RedIntensity * x, Col1.GreenIntensity * x, Col1.BlueIntensity * x);
}

RGB operator * (double x, const RGB& Col1)
{
  return RGB(Col1.RedIntensity * x, Col1.GreenIntensity * x, Col1.BlueIntensity * x);
}

// division of a RGB color by a real

RGB& RGB::operator /= (double& x)
{
  this->RedIntensity /= x;
  this->GreenIntensity /= x;
  this->BlueIntensity /= x;
  return *this;
}

RGB operator / (const RGB& Col1, double& x)
{
  return RGB(Col1.RedIntensity / x, Col1.GreenIntensity / x, Col1.BlueIntensity / x);
}

// truncate all wrong components of a RGB color
RGB RegularizeRGBColor (const RGB& Col)
{
  RGB TmpRGB = Col;
  if (TmpRGB.RedIntensity > 1)
    TmpRGB.RedIntensity = 1;
  if (TmpRGB.RedIntensity < 0)
    TmpRGB.RedIntensity = 0;
  if (TmpRGB.GreenIntensity > 1)
    TmpRGB.GreenIntensity = 1;
  if (TmpRGB.GreenIntensity < 0)
    TmpRGB.GreenIntensity = 0;
  if (TmpRGB.BlueIntensity > 1)
    TmpRGB.BlueIntensity = 1;
  if (TmpRGB.BlueIntensity < 0)
    TmpRGB.BlueIntensity = 0;
  return TmpRGB;
}
// compare two colors versus their intensity

bool RGB::operator < (const RGB& Col)
{
  double min1 = this->RedIntensity;
  double max1 = this->RedIntensity;
  if (min1 > this->GreenIntensity)
    min1 = this->GreenIntensity;
  if (min1 > this->BlueIntensity)
    min1 = this->BlueIntensity;
  if (max1 < this->GreenIntensity)
    max1 = this->GreenIntensity;
  if (max1 < this->BlueIntensity)
    max1 = this->BlueIntensity;
  double min2 = Col.RedIntensity;
  double max2 = Col.RedIntensity;
  if (min2 > Col.GreenIntensity)
    min2 = Col.GreenIntensity;
  if (min2 > Col.BlueIntensity)
    min2 = Col.BlueIntensity;
  if (max2 < Col.GreenIntensity)
    max2 = Col.GreenIntensity;
  if (max2 < Col.BlueIntensity)
    max2 = Col.BlueIntensity;
  if ((min1 + max1) < (min2 + max2))
    return true;
  else
    return false;
}

bool RGB::operator > (const RGB& Col)
{
  double min1 = this->RedIntensity;
  double max1 = this->RedIntensity;
  if (min1 > this->GreenIntensity)
    min1 = this->GreenIntensity;
  if (min1 > this->BlueIntensity)
    min1 = this->BlueIntensity;
  if (max1 < this->GreenIntensity)
    max1 = this->GreenIntensity;
  if (max1 < this->BlueIntensity)
    max1 = this->BlueIntensity;
  double min2 = Col.RedIntensity;
  double max2 = Col.RedIntensity;
  if (min2 > Col.GreenIntensity)
    min2 = Col.GreenIntensity;
  if (min2 > Col.BlueIntensity)
    min2 = Col.BlueIntensity;
  if (max2 < Col.GreenIntensity)
    max2 = Col.GreenIntensity;
  if (max2 < Col.BlueIntensity)
    max2 = Col.BlueIntensity;
  if ((min1 + max1) > (min2 + max2))
    return true;
  else
    return false;
}

bool RGB::operator <= (const RGB& Col)
{
  double min1 = this->RedIntensity;
  double max1 = this->RedIntensity;
  if (min1 > this->GreenIntensity)
    min1 = this->GreenIntensity;
  if (min1 > this->BlueIntensity)
    min1 = this->BlueIntensity;
  if (max1 < this->GreenIntensity)
    max1 = this->GreenIntensity;
  if (max1 < this->BlueIntensity)
    max1 = this->BlueIntensity;
  double min2 = Col.RedIntensity;
  double max2 = Col.RedIntensity;
  if (min2 > Col.GreenIntensity)
    min2 = Col.GreenIntensity;
  if (min2 > Col.BlueIntensity)
    min2 = Col.BlueIntensity;
  if (max2 < Col.GreenIntensity)
    max2 = Col.GreenIntensity;
  if (max2 < Col.BlueIntensity)
    max2 = Col.BlueIntensity;
  if ((min1 + max1) <= (min2 + max2))
    return true;
  else
    return false;
}

bool RGB::operator >= (const RGB& Col)
{
  double min1 = this->RedIntensity;
  double max1 = this->RedIntensity;
  if (min1 > this->GreenIntensity)
    min1 = this->GreenIntensity;
  if (min1 > this->BlueIntensity)
    min1 = this->BlueIntensity;
  if (max1 < this->GreenIntensity)
    max1 = this->GreenIntensity;
  if (max1 < this->BlueIntensity)
    max1 = this->BlueIntensity;
  double min2 = Col.RedIntensity;
  double max2 = Col.RedIntensity;
  if (min2 > Col.GreenIntensity)
    min2 = Col.GreenIntensity;
  if (min2 > Col.BlueIntensity)
    min2 = Col.BlueIntensity;
  if (max2 < Col.GreenIntensity)
    max2 = Col.GreenIntensity;
  if (max2 < Col.BlueIntensity)
    max2 = Col.BlueIntensity;
  if ((min1 + max1) >= (min2 + max2))
    return true;
  else
    return false;
}

