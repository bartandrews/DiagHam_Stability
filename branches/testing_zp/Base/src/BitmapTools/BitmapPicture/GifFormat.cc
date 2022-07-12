////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of Bitmap Picture Using Gif Format               //
//                     (interface to Toshio Kuratomi libungif)                //
//                                                                            //
//                        last modification : 28/02/2002                      //
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

 
#include "BitmapTools/BitmapPicture/GifFormat.h"
#include <stdio.h>
#ifdef __GIF__
extern "C"
{
#include <gif_lib.h>
}
#endif
#include <stdlib.h>


// default constructor
//

GifFormat::GifFormat () 
{
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

GifFormat::GifFormat (int width, int height) 
{
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

GifFormat::GifFormat (int width, int height, PicRGB& color) 
{
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

GifFormat::~GifFormat () 
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

void GifFormat::SetPixel(int x, int y, PicRGB& color) 
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

PicRGB GifFormat::GetPixel(int x, int y) 
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

void GifFormat::SetQuality (int quality)
{
}
 
// Save Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool GifFormat::SavePicture (char* fileName) 
{
  return false;
}

// Load Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool GifFormat::LoadPicture (char* fileName) 
{
#ifdef __GIF__
  GifFileType* File = DGifOpenFileName(fileName);
  if (File == 0) 
    return false;
  if (DGifSlurp(File) != GIF_OK)
    {
      DGifCloseFile(File);
      return false;
    }
  this->Width = File->SavedImages[0].ImageDesc.Width;
  this->Height = File->SavedImages[0].ImageDesc.Height;
  ColorMapObject* ColorMap = File->SavedImages[0].ImageDesc.ColorMap;
  if (ColorMap == 0)
    ColorMap = File->SColorMap;
  if (ColorMap->BitsPerPixel <= 8)
    {
      int ColorMapSize = ColorMap->ColorCount;	
      this->ColorMap = new PicRGB [ColorMapSize];
      int Pos = 0;
      for (int i = 0; i < ColorMapSize; i ++)
	{
	  this->ColorMap[i] = PicRGB (ColorMap->Colors[Pos].Red,
				      ColorMap->Colors[Pos].Green, 
				      ColorMap->Colors[Pos++].Blue);
	}
      this->Picture256Colors = new unsigned char* [2 * this->Height];
      Pos = 0;
      if (ColorMap->BitsPerPixel <= 8)
	{
	  unsigned char Mask = 0xff >> (8 - ColorMap->BitsPerPixel);
	  if (File->SavedImages[0].ImageDesc.Interlace == 0)
	    {
	      for (int j = 0; j < this->Height; j++)
		{
		  this->Picture256Colors[j] = new unsigned char [this->Width];
		  for (int i = 0; i < this->Width; i++)
		    {
		      this->Picture256Colors[j][i] = File->SavedImages[0].RasterBits[Pos++] & Mask;
		    }
		}      
	    }
	  else
	    {
	      int Pos2;
	      int TmpHeight = this->Height >> 3;
	      if ((TmpHeight << 3) == this->Height)
		TmpHeight--;
	      for (int j = 0; j <= TmpHeight; j++)
		{
		  Pos2 = (j << 3);
		  this->Picture256Colors[Pos2] = new unsigned char [this->Width];
		  for (int i = 0; i < this->Width; i++)
		    {
		      this->Picture256Colors[Pos2][i] = File->SavedImages[0].RasterBits[Pos++] & Mask;
		    }
		}
	      TmpHeight = this->Height >> 3;
	      if (((TmpHeight << 3) + 4) == this->Height)
		TmpHeight--;
	      for (int j = 0; j <= TmpHeight; j++)
		{
		  Pos2 = (j << 3) + 4;
		  this->Picture256Colors[Pos2] = new unsigned char [this->Width];
		  for (int i = 0; i < this->Width; i++)
		    {
		      this->Picture256Colors[Pos2][i] = File->SavedImages[0].RasterBits[Pos++] & Mask;
		    }
		}
	      TmpHeight = this->Height >> 2;
	      if (((TmpHeight << 2) + 2) == this->Height)
		TmpHeight--;
	      for (int j = 0; j <= TmpHeight; j++)
		{
		  Pos2 = (j << 2) + 2;
		  this->Picture256Colors[Pos2] = new unsigned char [this->Width];
		  for (int i = 0; i < this->Width; i++)
		    {
		      this->Picture256Colors[Pos2][i] = File->SavedImages[0].RasterBits[Pos++] & Mask;
		    }
		}
	      TmpHeight = this->Height >> 1;
	      if (((TmpHeight << 1) + 1) == this->Height)
		TmpHeight--;
	      for (int j = 0; j <= TmpHeight; j++)
		{
		  Pos2 = (j << 1) + 1;
		  this->Picture256Colors[Pos2] = new unsigned char [this->Width];
		  for (int i = 0; i < this->Width; i++)
		    {
		      this->Picture256Colors[Pos2][i] = File->SavedImages[0].RasterBits[Pos++] & Mask;
		    }
		}
	    }
	}
    }
  else
    if (File->SavedImages[0].ImageDesc.ColorMap->BitsPerPixel == 24)
      {
      }

  DGifCloseFile(File);
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

AbstractBitmapPicture* GifFormat::Clone (int width, int height)
{
  if (width == 0)
    width = this->Width;
  if (height == 0)
    height = this->Height;
  return new GifFormat (width, height);
}

#ifdef __OPENGL__


// Load bitmap picture as an open gl texture
//
// alpha = alpha factor to use if no alpha channel is defined in bitmap file (if greater than 1.0, use RGB texture only)

void GifFormat::GLLoadBitmapAsTexture (double alpha)
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
		  if (this->Picture256Colors[j][i] != 0)
		    this->GLTexture[Index++] = AlphaFactor;
		  else
		    this->GLTexture[Index++] = 0;
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

void GifFormat::GLUseBitmapAsTexture ()
{
  glBindTexture (GL_TEXTURE_2D, this->GLTextureID);
  glScaled(((double) this->Width) / ((double) this->GLWidth), ((double) this->Height) / ((double) this->GLHeight), 0.0);
}

#endif
