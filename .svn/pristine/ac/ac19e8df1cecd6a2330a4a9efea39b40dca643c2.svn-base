////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of Abstract Bitmap Picture                      //
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


#ifndef _ABSTRACTBITMAPPICTURE_H
#define _ABSTRACTBITMAPPICTURE_H


#include "config.h"
#include "Color/PicRGB.h"
#include "Color/Color.h"


class AbstractBitmapPicture;


// read a picture from a file and return corresponding AbstractBitmapPicture object
//
// fileName = file name 
// return value = pointer on corresponding AbstractBitmapPicture
AbstractBitmapPicture* ReadBitmapPicture (char* fileName);

// save a picture in a given format
//
// fileName = file name 
// inputPicture = picture to save
// quality = picture quality for output format (if not in the range 0 to 100, take input picture quality)
// return value = true if picture have been correctly saved
bool SaveBitmapPicture (char* fileName, AbstractBitmapPicture* inputPicture, int quality = -1);

// create a picture in a given format
//
// fileName = file name with an extension corresponding to the chosen file format
// width = picture width
// height = picture height
// return value = pointer to the created bitmap picture
AbstractBitmapPicture* CreateBitmapPicture (char* fileName, int width, int height);


class AbstractBitmapPicture
{

 public:

  // Return picture height
  //
  virtual int PictureHeight() = 0;

  // Return picture width
  //
  virtual int PictureWidth() = 0;
  
  // Set Pixel Value
  //
  // x = pixel coordinate along x axis
  // y = pixel coordinate along y axis
  // color = pixel color
  virtual void SetPixel(int x, int y, PicRGB& Col) = 0;
  
  // Set Pixel Value
  //
  // x = pixel coordinate along x axis
  // y = pixel coordinate along y axis
  // color = pixel color
  virtual void SetPixel(int x, int y, Color& color);
  
  // Get Pixel Value in 24 bits precision
  //
  // x = pixel coordinate along x axis
  // y = pixel coordinate along y axis
  // return value = pixel color
  virtual PicRGB GetPixel(int x, int y) = 0;

  // Get Pixel Value
  //
  // x = pixel coordinate along x axis
  // y = pixel coordinate along y axis
  // color = reference on color that will contain pixel color
  virtual void GetPixel(int x, int y, Color& color);

  // Save Picture
  //
  // fileName = picture filename
  // return value = true if operation is successfull
  virtual bool SavePicture (char* FileName) = 0;
  
  // Load Picture
  //
  // fileName = picture filename
  // return value = true if operation is successfull
  virtual bool LoadPicture (char* fileName) = 0;
  
  // virtual destructor
  //
  virtual ~AbstractBitmapPicture ();

  // resample a picture and store result in a new one
  //
  // width = new width
  // height = new height (0 if ratio has to be kept)
  // return value = pointer to resulting bitmap picture 
  virtual AbstractBitmapPicture* Resample (int width, int height = 0);

  // rotate a picture and store result in a new one (counter clockwise orientation)
  //
  // angle = rotation angle in degree
  // return value = pointer on the new resulting picture
  virtual AbstractBitmapPicture* Rotate (double angle);

  // set picture quality for destructive compression algorithm
  //
  // quality = picture quality (range from 0 to 100)
  virtual void SetQuality (int quality) {};
 
  // get picture quality for destructive compression algorithm
  //
  // return value = picture quality (range from 0 to 100)
  virtual int GetQuality ();
 
  // clone a bitmap picture type with a new size
  // 
  // width = new width (0 if old width has to be kept)
  // height = new height (0 if old height has to be kept)
  // return value = pointer to resulting bitmap picture 
  virtual AbstractBitmapPicture* Clone (int width = 0, int height = 0) = 0;

#ifdef __OPENGL__

  // Load bitmap picture as an open gl texture
  //
  // alpha = alpha factor to use if no alpha channel is defined in bitmap file (if greater than 1.0, use RGB texture only)
  virtual void GLLoadBitmapAsTexture (double alpha = 2.0) = 0;

  // use bitmap picture as current open gl texture
  //
  virtual void GLUseBitmapAsTexture () = 0;

#endif

};

#endif
