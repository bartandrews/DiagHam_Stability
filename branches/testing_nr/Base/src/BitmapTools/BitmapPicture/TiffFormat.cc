////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of Bitmap Picture Using Tiff Format               //
//                      (interface to Sam Leffler-SGI libtiff)                //
//                                                                            //
//                        last modification : 25/08/2001                      //
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

 
#include "BitmapTools/BitmapPicture/TiffFormat.h"
#include <stdio.h>
#ifdef __TIFF__
extern "C"
{
#include <tiffio.h>
}
#endif
#include <stdlib.h>


// default constructor
//

TiffFormat::TiffFormat () 
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

TiffFormat::TiffFormat (int width, int height) 
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

TiffFormat::TiffFormat (int width, int height, PicRGB& color) 
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

TiffFormat::~TiffFormat () 
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

void TiffFormat::SetPixel(int x, int y, PicRGB& color) 
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

PicRGB TiffFormat::GetPixel(int x, int y) 
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

void TiffFormat::SetQuality (int quality)
{
}
 
// Save Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool TiffFormat::SavePicture (char* fileName) 
{
  return false;
}

// Load Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool TiffFormat::LoadPicture (char* fileName) 
{
#ifdef __TIFF__
  TIFF* TiffHandle = TIFFOpen (fileName, "r");
  if (TiffHandle == 0)
    return false;
  uint32 TmpWidth;
  uint32 TmpHeight;
  TIFFGetField(TiffHandle, TIFFTAG_IMAGELENGTH, &TmpHeight);
  TIFFGetField(TiffHandle, TIFFTAG_IMAGEWIDTH, &TmpWidth);
  this->Width = TmpWidth;
  this->Height = TmpHeight;
  unsigned char* TmpBuffer = (unsigned char*) _TIFFmalloc(TIFFScanlineSize(TiffHandle)); 
  short PictureType = 0;
  TIFFGetField(TiffHandle, TIFFTAG_PHOTOMETRIC, &PictureType);
  switch (PictureType)
    {
    case 2:
      {
	int Pos = 0;
	this->Picture = new PicRGB* [this->Height];
	for (int j = 0; j < this->Height; j++)
	  {
	    Pos = 0;
	    TIFFReadScanline(TiffHandle, TmpBuffer, j, 2);
	    this->Picture[j] = new PicRGB [this->Width];
	    for (int i = 0; i < this->Width; i++)
	      {
		this->Picture[j][i] = PicRGB (TmpBuffer[Pos], TmpBuffer[Pos + 1], TmpBuffer[Pos + 2]);
		Pos += 3;
	      }
	  }      
      }
      break;
    case 3:
      {
	short BitsPerSample = 0;
	TIFFGetField(TiffHandle, TIFFTAG_BITSPERSAMPLE, &BitsPerSample);
	int ColorMapSize = 1 << BitsPerSample;	
	this->ColorMap = new PicRGB [ColorMapSize];
	int16* TmpRedColorMap;
	int16* TmpGreenColorMap;
	int16* TmpBlueColorMap ;
	TIFFGetField(TiffHandle, TIFFTAG_COLORMAP, &TmpRedColorMap, &TmpGreenColorMap, &TmpBlueColorMap);	
	int Pos = 0;
	for (int i = 0; i < ColorMapSize; i ++)
	  {
	    this->ColorMap[i] = PicRGB ((unsigned char) ((TmpRedColorMap[Pos] >> 8) & 0xff),
					(unsigned char) ((TmpGreenColorMap[Pos] >> 8) & 0xff),
					(unsigned char) ((TmpBlueColorMap[Pos] >> 8) & 0xff));
	    Pos++;
	  }
	if (BitsPerSample == 4)
	  {
	    this->Picture256Colors = new unsigned char* [this->Height];
	    for (int j = 0; j < this->Height; j++)
	      {
		Pos = 0;
		TIFFReadScanline(TiffHandle, TmpBuffer, j, 2);
		this->Picture256Colors[j] = new unsigned char [this->Width];
		for (int i = 0; i < this->Width; i++)
		  {
		    this->Picture256Colors[j][i] = (unsigned char) (TmpBuffer[Pos] & 0x7);
		    this->Picture256Colors[j][i++] = (unsigned char) ((TmpBuffer[Pos] >> 4) & 0x7);
		  }
	      }      
	  }
	else
	  {
	    this->Picture256Colors = new unsigned char* [this->Height];
	    for (int j = 0; j < this->Height; j++)
	      {
		Pos = 0;
		TIFFReadScanline(TiffHandle, TmpBuffer, j, 2);
		this->Picture256Colors[j] = new unsigned char [this->Width];
		for (int i = 0; i < this->Width; i++)
		  {
		    this->Picture256Colors[j][i] = TmpBuffer[Pos++];
		  }
	      }      
	  }
      }
      break;
    case 1:
      {
	short BitsPerSample = 0;
	TIFFGetField(TiffHandle, TIFFTAG_BITSPERSAMPLE, &BitsPerSample);
	int ColorMapSize = 1 << BitsPerSample;
	unsigned char ColorInc = 128 >> (BitsPerSample - 1);
	this->ColorMap = new PicRGB [ColorMapSize];
	unsigned short* TmpColorMap = (unsigned short*) _TIFFmalloc(6 * ColorMapSize);
	unsigned short Col = 0;
	TIFFGetField(TiffHandle, TIFFTAG_COLORMAP, &TmpColorMap);	
	int Pos = 0;
	for (int i = 0; i < ColorMapSize; i ++)
	  {
	    this->ColorMap[i] = PicRGB (Col, Col, Col);
	    Col += ColorInc;
	  }
	_TIFFfree(TmpColorMap);
	if (BitsPerSample == 4)
	  {
	    this->Picture256Colors = new unsigned char* [this->Height];
	    for (int j = 0; j < this->Height; j++)
	      {
		Pos = 0;
		TIFFReadScanline(TiffHandle, TmpBuffer, j, 2);
		this->Picture256Colors[j] = new unsigned char [this->Width];
		for (int i = 0; i < this->Width; i++)
		  {
		    this->Picture256Colors[j][i] = (unsigned char) (TmpBuffer[Pos] & 0x7);
		    this->Picture256Colors[j][i++] = (unsigned char) ((TmpBuffer[Pos] >> 4) & 0x7);
		  }
	      }      
	  }
	else
	  {
	    this->Picture256Colors = new unsigned char* [this->Height];
	    for (int j = 0; j < this->Height; j++)
	      {
		Pos = 0;
		TIFFReadScanline(TiffHandle, TmpBuffer, j, 2);
		this->Picture256Colors[j] = new unsigned char [this->Width];
		for (int i = 0; i < this->Width; i++)
		  {
		    this->Picture256Colors[j][i] = TmpBuffer[Pos++];
		  }
	      }      
	  }
      }
      break;
    }
  TIFFClose(TiffHandle);
  _TIFFfree(TmpBuffer);
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

AbstractBitmapPicture* TiffFormat::Clone (int width, int height)
{
  if (width == 0)
    width = this->Width;
  if (height == 0)
    height = this->Height;
  return new TiffFormat (width, height);
}

#ifdef __OPENGL__


// Load bitmap picture as an open gl texture
//
// alpha = alpha factor to use if no alpha channel is defined in bitmap file (if greater than 1.0, use RGB texture only)

void TiffFormat::GLLoadBitmapAsTexture (double alpha)
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

void TiffFormat::GLUseBitmapAsTexture ()
{
  glBindTexture (GL_TEXTURE_2D, this->GLTextureID);
  glScaled(((double) this->Width) / ((double) this->GLWidth), ((double) this->Height) / ((double) this->GLHeight), 0.0);
}

#endif
