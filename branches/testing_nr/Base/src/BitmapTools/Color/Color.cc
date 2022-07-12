////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                              class of Color                                //
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


#include "config.h"
#include "BitmapTools/Color/Color.h"
#include <math.h>

#ifdef __OPENGL__
#include <GL/gl.h>
#endif


//constructors

Color::Color ()
{
  this->RGBCol = RGB();
  this->AlphaChannel = 1.0l;
}

Color::Color (const double& Int)
{
  this->RGBCol = RGB(Int);
  this->AlphaChannel = 1.0l;
}

Color::Color (const double& RedInt, const double& GreenInt, const double& BlueInt)
{
  this->RGBCol = RGB(RedInt, GreenInt, BlueInt);
  this->AlphaChannel = 1.0l;
}

Color::Color (const RGB& Col)
{
  this->RGBCol = Col;
  this->AlphaChannel = 1.0l;
}


Color::Color (const Color& Col)
{
  this->RGBCol = Col.RGBCol;
  this->AlphaChannel = Col.AlphaChannel;
}

/*Color::Color (const PicRGB& Col)
{
  this->RGBCol = (RGB) Col;
}*/

Color::Color (const double& RedInt, const double& GreenInt, const double& BlueInt, double Alpha)
{
  this->RGBCol = RGB(RedInt, GreenInt, BlueInt);
  this->AlphaChannel = Alpha;
}

Color::Color (const RGB& Col, double Alpha)
{
  this->RGBCol = Col;
  this->AlphaChannel = Alpha;
}
 
/*Color::Color (const PicRGB& Col, double Alpha)
{
  this->RGBCol = (RGB) Col;
  this->AlphaChannel = Alpha;
}*/

// assignement

Color& Color::operator = (const Color& Col)
{
  this->RGBCol = Col.RGBCol;
  this->AlphaChannel = Col.AlphaChannel;
  return *this;
}

Color& Color::operator = (const RGB& Col)
{
  this->RGBCol = Col;
  this->AlphaChannel = 1.0l;
  return *this;
}

Color& Color::operator = (const PicRGB& Col)
{
  this->RGBCol = (RGB) Col;
  this->AlphaChannel = 1.0l;
  return *this;
}

// Get Color Component

RGB Color::GetRGBComponent()
{
  return this->RGBCol;
}

double Color::GetAlphaComponent()
{
  return this->AlphaChannel;
}

// Set Color Component

void Color::SetRGBComponent(const RGB& Col)
{
  this->RGBCol = Col;
}

void Color::SetAlphaComponent(double Alpha)
{
  this->AlphaChannel = Alpha;
}

// sum of two colors

Color& Color::operator += (const Color& Col2)
{
  this->AlphaChannel += Col2.AlphaChannel;
  this->RGBCol = this->RGBCol + Col2.RGBCol;
  return *this;
}

Color operator + (const Color& Col1, const Color& Col2)
{
  double TmpAlpha = Col1.AlphaChannel + Col2.AlphaChannel;
  return Color(Col1.RGBCol + Col2.RGBCol, TmpAlpha);
}

// difference of two colors

Color& Color::operator -= (const Color& Col2)
{
  this->AlphaChannel -= Col2.AlphaChannel;
  this->RGBCol = this->RGBCol - Col2.RGBCol;
  return *this;
}

Color operator - (const Color& Col1, const Color& Col2)
{
  double TmpAlpha = Col1.AlphaChannel - Col2.AlphaChannel;
  return Color(Col1.RGBCol - Col2.RGBCol, TmpAlpha);
}

// product of two colors

Color& Color::operator *= (const Color& Col2)
{
  this->RGBCol *= Col2.RGBCol;
  this->AlphaChannel *= Col2.AlphaChannel;  
  return *this;
}

Color operator * (const Color& Col1, const Color& Col2)
{
  return Color(Col1.RGBCol * Col2.RGBCol, Col1.AlphaChannel * Col2.AlphaChannel);
}

// product of a color with a real

Color& Color::operator *= (double x)
{
  this->RGBCol *= x;
  this->AlphaChannel *= x;
  return *this;
}

Color operator * (double x, const Color& Col1)
{
  return Color(x * Col1.RGBCol, x * Col1.AlphaChannel);
}

Color operator * (const Color& Col1, double x)
{
  return Color(x * Col1.RGBCol, x * Col1.AlphaChannel);
}

// division of a color by a real

Color& Color::operator /= (double x)
{
  this->RGBCol /= x;
  this->AlphaChannel /= x;
  return *this;
}

Color operator / (const Color& Col1, double x)
{
  return Color(Col1.RGBCol / x, Col1.AlphaChannel / x);
}


// truncate all wrong components of a color

Color RegularizeColor (const Color& Col)
{
  return Color (RegularizeRGBColor(Col.RGBCol), Col.AlphaChannel);
}

// compare two colors versus their intensity

bool Color::operator < (const Color& Col)
{
  return (this->RGBCol < Col.RGBCol);
}

bool Color::operator > (const Color& Col)
{
  return (this->RGBCol > Col.RGBCol);
}

bool Color::operator <= (const Color& Col)
{
  return (this->RGBCol <= Col.RGBCol);
}

bool Color::operator >= (const Color& Col)
{
  return (this->RGBCol >= Col.RGBCol);
}

// return color norm

double Color::norm ()
{
  return sqrt(0.3333333333 * (this->RGBCol.RedIntensity * this->RGBCol.RedIntensity
			      + this->RGBCol.GreenIntensity * this->RGBCol.GreenIntensity
			      + this->RGBCol.BlueIntensity * this->RGBCol.BlueIntensity));
}

// OpenGL methods

#ifdef __OPENGL__  

// Push Color on current OpenGL stack
//  

void Color::GLPushColor ()
{
  glColor4d(this->RGBCol.RedIntensity, this->RGBCol.GreenIntensity, this->RGBCol.BlueIntensity, this->AlphaChannel);
}

#endif

