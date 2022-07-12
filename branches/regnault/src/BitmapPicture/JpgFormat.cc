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

 
#include "BitmapPicture/JpgFormat.h"
#include <stdio.h>
#include <setjmp.h>
#ifdef __JPEG__
extern "C"
{
#include <jpeglib.h>
#include <jerror.h>
}
#endif
#include <stdlib.h>


// default constructor
//

JpgFormat::JpgFormat () 
{
  this->Quality = 100;
  this->Height = 0;
  this->Width = 0;
  this->Picture = 0;
  this->Picture256Colors = 0;
  this->ColorMap = 0;
  this->NbrColor = 0;
#ifdef __OPENGL__
  this->GLTextureID = 0;
#endif
}

// constructor for an empty black picture
//
// width = picture width
// height = picture height

JpgFormat::JpgFormat (int width, int height) 
{
  this->Quality = 100;
  this->Height = height;
  this->Width = width;
  this->Picture = new PicRGB* [this->Height];
  for (int i = 0; i < this->Height; i++)
    {
      this->Picture[i] = new PicRGB [this->Width];
    }
  this->ColorMap = 0;
  this->Picture256Colors = 0;
  this->NbrColor = 0;
#ifdef __OPENGL__
  this->GLTextureID = 0;
#endif
}

// constructor for an empty picture using user background color 
//
// width = picture width
// height = picture height
// color = background color

JpgFormat::JpgFormat (int width, int height, PicRGB& color) 
{
  this->Quality = 100;
  this->Height = height;
  this->Width = width;
  this->Picture = new PicRGB* [this->Height];
  for (int i = 0; i < this->Height; i++)
    {
     this->Picture[i] = new PicRGB [this->Width];
      for (int j = 0; j < this->Width; j++)
      this->Picture[i][j] = color;
    }
  this->Picture256Colors = 0;
  this->ColorMap = 0;
  this->NbrColor = 0;
#ifdef __OPENGL__
  this->GLTextureID = 0;
#endif
}
  
// destructor
//

JpgFormat::~JpgFormat () 
{
  if (this->Picture != 0)
    {
      for (int i = 0; i < this->Height; i++)
	delete[] this->Picture[i];
      delete[] this->Picture;
    }
  if (this->ColorMap != 0)
    {
      delete[] this->ColorMap;
    }
  if (this->Picture256Colors != 0)
    {
      for (int i = 0; i < this->Height; i++)
	delete[] this->Picture256Colors[i];
      delete[] this->Picture256Colors;
    }  
#ifdef __OPENGL__
  if (this->GLTexture != 0)
    delete[] this->GLTexture;
#endif
}
  
// Set Pixel Value
//
// x = pixel coordinate along x axis
// y = pixel coordinate along y axis
// color = pixel color

void JpgFormat::SetPixel(int x, int y, PicRGB& color) 
{
  if (this->Picture != 0)
    {
      this->Picture[y][x] = color;
    }
  if (this->Picture256Colors != 0)
    {
      int Dif = 765;
      int Index = 0;
      int CurrentDif;
      for (int i = 0; i < this->NbrColor; i++)
	{
	  CurrentDif = abs((int)this->ColorMap[i].Red - (int)color.Red) +
	    abs((int)this->ColorMap[i].Green - (int)color.Green) +
	    abs((int)this->ColorMap[i].Blue - (int)color.Blue);
	  if (CurrentDif < Dif)
	    {
	      Index = i;
	      Dif = CurrentDif;
	    }
	}
      this->Picture256Colors [y][x] = (unsigned char) Index;
    }
  return;
}
  
// Get Pixel Value
//
// x = pixel coordinate along x axis
// y = pixel coordinate along y axis
// reurn value = pixel color

PicRGB JpgFormat::GetPixel(int x, int y) 
{
  if (this->Picture != 0)
    {
      return this->Picture[y][x];
    }
  if (this->Picture256Colors != 0)
    {
      return this->ColorMap[this->Picture256Colors[y][x]];
    }
  return PicRGB();
}
  
// set picture quality for destructive compression algorithm
//
// quality = picture quality (range from 0 to 100)

void JpgFormat::SetQuality (int quality)
{
  this->Quality = quality;
  if (this->Quality > 100)
    this->Quality = 100;
  else
    if (this->Quality < 0)
      this->Quality = 0;
}
 
// Save Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool JpgFormat::SavePicture (char* fileName) 
{
#ifdef __JPEG__
  struct jpeg_compress_struct JpegDef;
  struct jpeg_error_mgr JpegError;
  FILE* File; 
  if ((File= fopen (fileName, "wb")) == 0)
    return false;
  JpegDef.err = jpeg_std_error(&JpegError);
  jpeg_create_compress(&JpegDef);                                               
  JpegDef.image_height = this->Height;
  JpegDef.image_width = this->Width;
  if (Picture256Colors != 0)
    {
      JpegDef.input_components = 1;
      JpegDef.in_color_space = JCS_GRAYSCALE;
      jpeg_set_defaults(&JpegDef);
      jpeg_set_quality(&JpegDef, this->Quality, TRUE);
      jpeg_stdio_dest(&JpegDef, File);
      jpeg_start_compress(&JpegDef, TRUE);
      while (JpegDef.next_scanline < JpegDef.image_height)
	{
	  jpeg_write_scanlines(&JpegDef, &(this->Picture256Colors[JpegDef.next_scanline]), 1);
	}      
    }
  else
    {
      JpegDef.input_components = 3;
      JpegDef.in_color_space = JCS_RGB;
      jpeg_set_defaults(&JpegDef);
      jpeg_set_quality(&JpegDef, this->Quality, TRUE);
      jpeg_stdio_dest(&JpegDef, File);
      jpeg_start_compress(&JpegDef, TRUE);
      unsigned char* TmpBuffer = new unsigned char [3 * this->Width];
      while (JpegDef.next_scanline < JpegDef.image_height)
	{
	  for (int i = 0; i < this->Width; i++)
	    {
	      TmpBuffer[3 * i] = this->Picture[JpegDef.next_scanline][i].Red;
	      TmpBuffer[3 * i + 1] = this->Picture[JpegDef.next_scanline][i].Green;
	      TmpBuffer[3 * i + 2] = this->Picture[JpegDef.next_scanline][i].Blue; 
	    }
	  jpeg_write_scanlines(&JpegDef, &TmpBuffer, 1);
	}      
      delete[] TmpBuffer;
    }
  jpeg_finish_compress(&JpegDef);
  jpeg_destroy_compress(&JpegDef);
  fclose (File);
  return true;
#else
  return false;
#endif
}

// Load Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool JpgFormat::LoadPicture (char* fileName) 
{
#ifdef __JPEG__
  struct jpeg_decompress_struct JpegDef;
  struct jpeg_error_mgr JpegError;
  FILE* File; 
  if ((File= fopen (fileName, "rb")) == 0)
    return false;
  JpegDef.err = jpeg_std_error(&JpegError);
  jpeg_create_decompress(&JpegDef);                                               
  jpeg_stdio_src(&JpegDef, File);
  jpeg_read_header(&JpegDef, TRUE);
  jpeg_start_decompress(&JpegDef);
  this->Height = JpegDef.output_height;
  this->Width = JpegDef.output_width;

  if (JpegDef.output_components == 1)
    {
      this->Picture = 0;
      this->Picture256Colors = new unsigned char* [this->Height];
      this->NbrColor = JpegDef.actual_number_of_colors;
      this->NbrColor = 256;
      this->ColorMap = new PicRGB [this->NbrColor];
      unsigned char col = 0; 
      for (int i = 0; i < this->NbrColor; i++)
	{
	  this->ColorMap[i] = PicRGB(col, col, col);
	  col++;
	}
      while (JpegDef.output_scanline < JpegDef.output_height)
	{
	  this->Picture256Colors[JpegDef.output_scanline] = new unsigned char [this->Width];
	  jpeg_read_scanlines(&JpegDef, &(this->Picture256Colors[JpegDef.output_scanline]), 1);
	}
    }
  else
    {
      this->Picture256Colors = 0;
      this->NbrColor = 0;
      this->ColorMap = 0;
      this->Picture = new PicRGB* [this->Height];
      unsigned char* TmpBuffer = new unsigned char [3 * this->Width];
      while (JpegDef.output_scanline < JpegDef.output_height)
	{
	  this->Picture[JpegDef.output_scanline] = new PicRGB [this->Width];
	  jpeg_read_scanlines(&JpegDef, &TmpBuffer, 1);
	  for (int i = 0; i < this->Width; i++)
	    this->Picture[JpegDef.output_scanline - 1][i] = PicRGB (TmpBuffer[3 * i], TmpBuffer[3 * i + 1], 
								TmpBuffer[3 * i + 2]);
	}      
      delete[] TmpBuffer;
    }
  jpeg_finish_decompress(&JpegDef);
  jpeg_destroy_decompress(&JpegDef);
  fclose (File);
  return true;
#else
  return false;
#endif
}

// clone a bitmap picture type with a new size
// 
// width = new width (0 if old width has to be kept)
// height = new height (0 if old height has to be kept)
// return value = pointer to resulting bitmap picture 

AbstractBitmapPicture* JpgFormat::Clone (int width, int height)
{
  if (width == 0)
    width = this->Width;
  if (height == 0)
    height = this->Height;
  return new JpgFormat (width, height);
}

#ifdef __OPENGL__


// Load bitmap picture as an open gl texture
//
// alpha = alpha factor to use if no alpha channel is defined in bitmap file (if greater than 1.0, use RGB texture only)

void JpgFormat::GLLoadBitmapAsTexture (double alpha)
{
  if (this->GLTextureID != 0)
    return;

  glPixelStorei (GL_UNPACK_ALIGNMENT, 1);
  glGenTextures(1, &(this->GLTextureID));
  glBindTexture (GL_TEXTURE_2D, this->GLTextureID);

  this->GLWidth = 1;
  while (this->GLWidth < this->Width)
    this->GLWidth <<=  1;
  this->GLHeight = 1;
  while (this->GLHeight < this->Height)
    this->GLHeight <<= 1;
  if ((alpha > 1.0) || (alpha < 0.0))
    {
      this->GLTexture = new unsigned char [this->GLWidth * this->GLHeight * 3];
      int Index = 0;
      if (this->Picture256Colors == 0)
	{
	  int j = 0;
	  for (; j < this->Height; j++)
	    {
	      int i = 0;
	      for (; i < this->Width; i++)
		{
		  this->GLTexture[Index++] = this->Picture[j][i].Red;
		  this->GLTexture[Index++] = this->Picture[j][i].Green;
		  this->GLTexture[Index++] = this->Picture[j][i].Blue;
		}
	      for (; i < this->GLWidth; i++)
		{
		  this->GLTexture[Index++] = 0;
		  this->GLTexture[Index++] = 0;
		  this->GLTexture[Index++] = 0;
		}
	    }
	  for (; j < this->GLHeight; j++)
	    for (int i = 0; i < this->GLWidth; i++)
	      {
		this->GLTexture[Index++] = 0;
		this->GLTexture[Index++] = 0;
		this->GLTexture[Index++] = 0;
	      }
	}
      else
	{
	  int j = 0;
	  for (; j < this->Height; j++)
	    {
	      int i = 0;
	      for (; i < this->Width; i++)
		{
		  this->GLTexture[Index++] = this->ColorMap[this->Picture256Colors[j][i]].Red;
		  this->GLTexture[Index++] = this->ColorMap[this->Picture256Colors[j][i]].Green;
		  this->GLTexture[Index++] = this->ColorMap[this->Picture256Colors[j][i]].Blue;
		}
	      for (; i < this->GLWidth; i++)
		{
		  this->GLTexture[Index++] = 0;
		  this->GLTexture[Index++] = 0;
		  this->GLTexture[Index++] = 0;
		}
	    }
	  for (; j < this->GLHeight; j++)
	    for (int i = 0; i < this->GLWidth; i++)
	      {
		this->GLTexture[Index++] = 0;
		this->GLTexture[Index++] = 0;
		this->GLTexture[Index++] = 0;
	      }
	}      
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, this->GLWidth, this->GLHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, this->GLTexture);
    }
  else
    {
      this->GLTexture = new unsigned char [this->GLWidth * this->GLHeight * 4];
      int Index = 0;
      unsigned char AlphaFactor = (unsigned char) (alpha * 255);
      if (this->Picture256Colors == 0)
	{
	  int j = 0;
	  for (; j < this->Height; j++)
	    {
	      int i = 0;
	      for (; i < this->Width; i++)
		{
		  this->GLTexture[Index++] = this->Picture[j][i].Red;
		  this->GLTexture[Index++] = this->Picture[j][i].Green;
		  this->GLTexture[Index++] = this->Picture[j][i].Blue;
		  this->GLTexture[Index++] = AlphaFactor;
		}
	      for (; i < this->GLWidth; i++)
		{
		  this->GLTexture[Index++] = 0;
		  this->GLTexture[Index++] = 0;
		  this->GLTexture[Index++] = 0;
		  this->GLTexture[Index++] = 0;
		}
	    }
	  for (; j < this->GLHeight; j++)
	    for (int i = 0; i < this->GLWidth; i++)
	      {
		this->GLTexture[Index++] = 0;
		this->GLTexture[Index++] = 0;
		this->GLTexture[Index++] = 0;
		this->GLTexture[Index++] = 0;
	      }
	}
      else
	{
	  int j = 0;
	  for (; j < this->Height; j++)
	    {
	      int i = 0;
	      for (; i < this->Width; i++)
		{
		  this->GLTexture[Index++] = this->ColorMap[this->Picture256Colors[j][i]].Red;
		  this->GLTexture[Index++] = this->ColorMap[this->Picture256Colors[j][i]].Green;
		  this->GLTexture[Index++] = this->ColorMap[this->Picture256Colors[j][i]].Blue;
		  this->GLTexture[Index++] = AlphaFactor;
		}
	      for (; i < this->GLWidth; i++)
		{
		  this->GLTexture[Index++] = 0;
		  this->GLTexture[Index++] = 0;
		  this->GLTexture[Index++] = 0;
		  this->GLTexture[Index++] = 0;
		}
	    }
	  for (; j < this->GLHeight; j++)
	    for (int i = 0; i < this->GLWidth; i++)
	      {
		this->GLTexture[Index++] = 0;
		this->GLTexture[Index++] = 0;
		this->GLTexture[Index++] = 0;
		this->GLTexture[Index++] = 0;
	      }
	}  
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, this->GLWidth, this->GLHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, this->GLTexture);
    }
}

// use bitmap picture as current open gl texture
//

void JpgFormat::GLUseBitmapAsTexture ()
{
  glBindTexture (GL_TEXTURE_2D, this->GLTextureID);
  glScaled(((double) this->Width) / ((double) this->GLWidth), ((double) this->Height) / ((double) this->GLHeight), 0.0);
}

#endif
