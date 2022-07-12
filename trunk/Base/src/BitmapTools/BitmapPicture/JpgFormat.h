////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of Bitmap Picture Using Jpg Format                //
//                     (interface to Thomas G. Lane's jpeglib)                //
//                                                                            //
//                        last modification : 09/04/2001                      //
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
  

#ifndef JPGFORMAT_H
#define JPGFORMAT_H


#include "config.h"
#include <fstream>
#include "BitmapTools/BitmapPicture/AbstractBitmapPicture.h"

#ifdef __OPENGL__
#include <GL/gl.h>
#endif


using std::ofstream;
using std::ifstream;


class JpgFormat : public AbstractBitmapPicture
{


protected:

  int Height;
  int Width;

  int Quality;
  
  PicRGB** Picture;
  
  unsigned char** Picture256Colors;
  PicRGB* ColorMap;
  int NbrColor;

#ifdef __OPENGL__

  unsigned char* GLTexture;
  GLuint GLTextureID;

  int GLWidth; 
  int GLHeight;
  
#endif
  
public:

  // default constructor
  //
  JpgFormat ();

  // constructor for an empty black picture
  //
  // width = picture width
  // height = picture height
  JpgFormat (int width, int height);

  // constructor for an empty picture using user background color 
  //
  // width = picture width
  // height = picture height
  // color = background color
  JpgFormat (int width, int height, PicRGB& color);
  
  // destructor
  //
  ~JpgFormat ();
  
  // Return picture height
  // 
  // return value = picture height
  int PictureHeight() {return this->Height;};

  // Return picture width
  // 
  // return value = picture width
  int PictureWidth() {return this->Width;};

  // Set Pixel Value
  //
  // x = pixel coordinate along x axis
  // y = pixel coordinate along y axis
  // color = pixel color
  void SetPixel(int x, int y, PicRGB& color);
  
  // Get Pixel Value
  //
  // x = pixel coordinate along x axis
  // y = pixel coordinate along y axis
  // reurn value = pixel color
  PicRGB GetPixel(int x, int y);
  
  // Save Picture
  // 
  // fileName = picture filename
  // return value = true if operation is successfull
  bool SavePicture (char* fileName);

  // Load Picture
  // 
  // fileName = picture filename
  // return value = true if operation is successfull
  bool LoadPicture (char* fileName);

  // set picture quality for destructive compression algorithm
  //
  // quality = picture quality (range from 0 to 100)
  void SetQuality (int quality);
 
  // clone a bitmap picture type with a new size
  // 
  // width = new width (0 if old width has to be kept)
  // height = new height (0 if old height has to be kept)
  // return value = pointer to resulting bitmap picture 
  AbstractBitmapPicture* Clone (int width = 0, int height = 0);

#ifdef __OPENGL__

  // Load bitmap picture as an open gl texture
  //
  // alpha = alpha factor to use if no alpha channel is defined in bitmap file (if greater than 1.0, use RGB texture only)
  void GLLoadBitmapAsTexture (double alpha = 2.0);

  // use bitmap picture as current open gl texture
  //
  void GLUseBitmapAsTexture ();

#endif

};

#endif
