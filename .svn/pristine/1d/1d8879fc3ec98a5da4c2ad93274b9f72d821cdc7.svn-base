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


#ifndef COLOR_H
#define COLOR_H


#include "config.h"
#include "BitmapTools/Color/RGB.h"
#include "BitmapTools/Color/PicRGB.h"

class Color
{

public:

  // intensity of each color component
  RGB RGBCol;
  double AlphaChannel;
  
  //constructors
  Color ();
  Color (const double& Int);
  Color (const double& RedInt, const double& GreenInt, const double& BlueInt);
  Color (const RGB& Col);
//  Color (const PicRGB& Col);
  Color (const Color& Col);
  Color (const double& RedInt, const double& GreenInt, const double& BlueInt, double Alpha);
  Color (const RGB& Col, double Alpha);
//  Color (const PicRGB& Col, double Alpha);
  
  // assignement
  Color& operator = (const Color& Col);
  Color& operator = (const RGB& Col);
  Color& operator = (const PicRGB& Col);
  Color& operator = (double Alpha); 

  //destructor
  ~Color () {};
  
  // Get Color Component
  RGB GetRGBComponent();
  double GetAlphaComponent();

  // Set Color Component
  void SetRGBComponent(const RGB& Col);
  void SetAlphaComponent(double Alpha);

  // sum of two colors
  Color& operator += (const Color& Col2);
  friend Color operator + (const Color& Col1, const Color& Col2);
  
  // difference of two colors
  Color& operator -= (const Color& Col2);
  friend Color operator - (const Color& Col1, const Color& Col2);
  
  // product of two colors
  Color& operator *= (const Color& Col2);
  friend Color operator * (const Color& Col1, const Color& Col2);
  
  // product of a color with a real
  Color& operator *= (double x);
  friend Color operator * (double x, const Color& Col1);
  friend Color operator * (const Color& Col1, double x);
  
  // division of a color by a real
  Color& operator /= (double x);
  friend Color operator / (const Color& Col1, double x);
  
  // truncate all wrong components of a color
  friend Color RegularizeColor (const Color& Col);

  // compare two colors versus their intensity
  bool operator < (const Color& Col);
  bool operator > (const Color& Col);
  
  bool operator <= (const Color& Col);
  bool operator >= (const Color& Col);
    
  // return color norm
  double norm ();

  // OpenGL methods

#ifdef __OPENGL__  

  // Push Color on current OpenGL stack
  //  
  void GLPushColor ();

#endif

};

#endif

