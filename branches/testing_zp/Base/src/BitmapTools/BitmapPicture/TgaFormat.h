////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of Bitmap Picture Using Targa Format 24 and 32 bits         //
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


#ifndef TGAFORMAT_H
#define TGAFORMAT_H


#include "config.h"
#include "BitmapTools/BitmapPicture/AbstractBitmapPicture.h"

#ifdef __OPENGL__
#include <GL/gl.h>
#endif


class TgaFormat : public AbstractBitmapPicture
{

protected:

  int Height;
  int Width;
  
  PicRGB** Picture;
  
  unsigned char** Picture256Colors;
  PicRGB* ColorMap;
  int NbrColor;
  
  char* ImageID;
  int ImageIDLenght;
  
#ifdef __OPENGL__

  unsigned char* GLTexture;
  GLuint GLTextureID;

  int GLWidth; 
  int GLHeight;
  
#endif
  
public:

  enum TgaMode
  {
    UnPackedColorMapped = 0x01,
    UnPacked = 0x02,
    PackedColorMapped = 0x09,
    Packed = 0x0a,
    BottomLeft = 0x00,
    BottomRight = 0x10,
    TopLeft = 0x20,
    TopRight = 0x30
  };

  // default constructor
  //
  TgaFormat ();

  // constructor for an empty black picture
  //
  // width = picture width
  // height = picture height
  TgaFormat (int L, int H);

  // constructor for an empty picture using user background color 
  //
  // width = picture width
  // height = picture height
  // color = background color
  TgaFormat (int L, int H, PicRGB& Col);

  // destructor
  //
  ~TgaFormat ();
  
  // Return size of the picture
  // 
  // return value = picture height
  int PictureHeight() {return Height;};

  // Return picture width
  // 
  // return value = picture width
  int PictureWidth() {return Width;};

  // Set Pixel Value
  //
  // x = pixel coordinate along x axis
  // y = pixel coordinate along y axis
  // color = pixel color
  void SetPixel(int x, int y, PicRGB& Col);
  
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
  bool SavePicture (char* FileName);

  // Save Picture
  // 
  // fileName = picture filename
  // mode = TGA mode to use
  // return value = true if operation is successfull
  bool SavePicture (char* FileName, int mode);

  // Load Picture
  // 
  // fileName = picture filename
  // return value = true if operation is successfull
  bool LoadPicture (char* FileName);
  
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
