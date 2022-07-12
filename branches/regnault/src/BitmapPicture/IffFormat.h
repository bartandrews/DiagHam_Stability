////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of Bitmap Picture using iff file format              //
//                                                                            //
//                        last modification : 08/09/2002                      //
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
  

#ifndef IFFFORMAT_H
#define IFFFORMAT_H


#include "config.h"
#include <fstream>
#include "BitmapPicture/AbstractBitmapPicture.h"

#ifdef __OPENGL__
#include <GL/gl.h>
#endif


class IffFormat : public AbstractBitmapPicture
{


protected:

  // picture height
  int Height;
  // picture width
  int Width;
  
  // number of bitplanes
  int NbrBitplanes;
  // masking mode
  int MaskingMode;
  // compression mode;
  int CompressionMode;

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
  IffFormat ();

  // constructor for an empty black picture
  //
  // width = picture width
  // height = picture height
  IffFormat (int L, int H);

  // constructor for an empty picture using user background color 
  //
  // width = picture width
  // height = picture height
  // color = background color
  IffFormat (int L, int H, PicRGB& Col);
  
  //destructor
  //
  ~IffFormat ();
  
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

 private:

  // read tag
  //
  // file = reference on current input file stream
  // return value = tag
  unsigned long ReadTag(ifstream& file);

  // proceed header chunk
  //
  // file = reference on current input file stream
  // chunkSize = chunk size
  void ProceedHeader(ifstream& file, unsigned long chunkSize);

  // proceed color map chunk
  //
  // file = reference on current input file stream
  // chunkSize = chunk size
  void ProceedColorMap(ifstream& file, unsigned long chunkSize);

  // proceed body chunk
  //
  // file = reference on current input file stream
  // chunkSize = chunk size
  void ProceedBody(ifstream& file, unsigned long chunkSize);

  // unpack a scanline 
  //
  // file = reference on current input file stream
  // scanline = pointer to array where unpacked scanline has to be stored
  // scanlineSize = number of byte per scanline
  void UnpackScanline (ifstream& file, unsigned char* scanline, int scanlineSize);

};

#endif
