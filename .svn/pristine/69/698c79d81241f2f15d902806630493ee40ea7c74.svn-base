////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of Bitmap Picture Using Hdr Format (Radiance)            //
//                                                                            //
//                        last modification : 27/05/2002                      //
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
  

#ifndef HDRFORMAT_H
#define HDRFORMAT_H


#include "config.h"
#include <fstream>
#include "BitmapTools/BitmapPicture/AbstractBitmapPicture.h"
#include "BitmapTools/Color/RGB.h"
#include "BitmapTools/Color/PicRGBE.h"
#include "GeneralTools/List.h"

#ifdef __OPENGL__
#include <GL/gl.h>
#endif


using std::ofstream;
using std::ifstream;


class HdrFormat : public AbstractBitmapPicture
{


protected:

  int Height;
  int Width;

  // horizontal orientation 
  bool LeftToRight;
  // vertical orientation 
  bool TopToBottom;
  
  // temporary scanline
  PicRGBE* Scanline;

  RGB** Picture;

  // list of all extracted commands
  List<char*> Commands;
  // list of the values associated to each command
  List<char*> CommandArgument;

#ifdef __OPENGL__

  unsigned char* GLTexture;
  GLuint GLTextureID;

  int GLWidth; 
  int GLHeight;
  
#endif
  
public:

  // default constructor
  //
  HdrFormat ();

  // constructor for an empty black picture
  //
  // width = picture width
  // height = picture height
  HdrFormat (int L, int H);

  // constructor for an empty picture using user background color 
  //
  // width = picture width
  // height = picture height
  // color = background color
  HdrFormat (int L, int H, RGB& Col);
  
  //destructor
  //
  ~HdrFormat ();
  
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
  
   // Set Pixel Value
  //
  // x = pixel coordinate along x axis
  // y = pixel coordinate along y axis
  // color = pixel color
  void SetPixel(int x, int y, Color& color);
  
 // Get Pixel Value
  //
  // x = pixel coordinate along x axis
  // y = pixel coordinate along y axis
  // reurn value = pixel color
  PicRGB GetPixel(int x, int y);
  
  // Get Pixel Value
  //
  // x = pixel coordinate along x axis
  // y = pixel coordinate along y axis
  // color = reference on color that will contain pixel color
  void GetPixel(int x, int y, Color& color);

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

  // read command in a hdr file 
  //
  // file = reference on input file
  // fileSize = input file size (minus number of already read bytes)
  // string = temporary string that has to be used for storing command
  // return value = true if it was the last command
  bool ReadCommand (ifstream& file, int fileSize, char* string);

  // read a scanline in a hdr file 
  //
  // file = reference on input file
  // pos = poistion where to store the scanline
  // return value = true if the scanline has beeen read
  bool ReadScanline (ifstream& file, int pos);

  // read a scanline in a hdr file using old scanline format 
  //
  // file = reference on input file
  // pos = poistion where to store the scanline
  // return value = true if the scanline has beeen read
  bool ReadOldScanline (ifstream& file, int pos);

  // write a scanline in a hdr file 
  //
  // file = reference on input file
  // return value = true if the scanline has been written
  bool WriteScanline (ofstream& file);

};

#endif
