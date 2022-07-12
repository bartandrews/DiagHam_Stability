////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of Bitmap Picture Using PCX Format 24 bits              //
//                                                                            //
//                        last modification : 31/08/2001                      //
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


#include "BitmapPicture/PcxFormat.h"
#include "Color/PicBGR.h"
#include <stdlib.h>


using std::ios;


// default constructor
//

PcxFormat::PcxFormat ()
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

PcxFormat::PcxFormat (int L, int H)
{
  this->Height = H;
  this->Width = L;
  this->Picture256Colors = new unsigned char* [this->Height];
  for (int i = 0; i < Width; i++)
    this->Picture256Colors[i] = new unsigned char [this->Width];
  this->ColorMap = new PicRGB[256];
  this->NbrColor = 256;
#ifdef __OPENGL__
  this->GLTextureID = 0;
#endif
}

// constructor for an empty picture using user background color 
//
// width = picture width
// height = picture height
// color = background color

PcxFormat::PcxFormat (int L, int H, PicRGB& Col)
{
  this->Height = H;
  this->Width = L;
  this->Picture = new PicRGB* [this->Width];
  for (int i = 0; i < this->Width; i++)
    {
     this->Picture[i] = new PicRGB [this->Height];
      for (int j = 0; j < this->Height; j++)
      this->Picture[i][j] = Col;
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

PcxFormat::~PcxFormat ()
{
  if (this->Picture256Colors != 0)
    {
      for (int i = 0; i < this->Height; i++)
	delete[] this->Picture256Colors[i];
      delete[] this->Picture256Colors;
      delete[] this->ColorMap;
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

void PcxFormat::SetPixel(int x, int y, PicRGB& Col)
{
  int Dif = 765;
  int Index = 0;
  int CurrentDif;
  for (int i = 0; i < this->NbrColor; i++)
    {
      CurrentDif = abs((int)this->ColorMap[i].Red - (int)Col.Red) +
	abs((int)this->ColorMap[i].Green - (int)Col.Green) +
	abs((int)this->ColorMap[i].Blue - (int)Col.Blue);
      if (CurrentDif < Dif)
	{
	  Index = i;
	  Dif = CurrentDif;
	}
    }
  this->Picture256Colors [y][x] = (unsigned char) Index;
}

// Get Pixel Value
//
// x = pixel coordinate along x axis
// y = pixel coordinate along y axis
// reurn value = pixel color

PicRGB	PcxFormat::GetPixel(int x, int y)
{
  if ((x < 0) || (x >= this->Width) || (y < 0) || (y >= this->Height))
    return PicRGB(0, 0, 0);
  return this->ColorMap[this->Picture256Colors [y][x]];
}

// Save Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool PcxFormat::SavePicture (char* FileName)
{
  return true;
}

// Load Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool PcxFormat::LoadPicture (char* FileName)
{
  ifstream File;
  unsigned char Tmp;
  File.open(FileName, ios::binary | ios::in);
  File.read ((char*) &Tmp, sizeof(unsigned char));
  unsigned char Version;
  File.read ((char*) &Version, sizeof(unsigned char));
  unsigned char Encoding;
  File.read ((char*) &Encoding, sizeof(unsigned char));
  unsigned char BitPerPixel;
  File.read ((char*) &BitPerPixel, sizeof(unsigned char));
  short XMin;
  short XMax;
  short YMin;
  short YMax;
  File.read ((char*) (&XMin), sizeof(short));
  File.read ((char*) (&YMin), sizeof(short));
  File.read ((char*) (&XMax), sizeof(short));
  File.read ((char*) (&YMax), sizeof(short));
  this->Width = XMax - XMin + 1;
  this->Height = YMax - YMin + 1;
  File.seekg(4, ios::cur);
  File.seekg(65, ios::beg);
  unsigned char NbrBitplane;
  File.read ((char*) &NbrBitplane, sizeof(unsigned char));
  unsigned short ScanlineSize;
  File.read ((char*) (&ScanlineSize), sizeof(unsigned short));
  unsigned short PaletteInfo;
  File.read ((char*) (&PaletteInfo), sizeof(unsigned short));
  this->NbrColor = 1 << BitPerPixel;
  if (BitPerPixel == 8)
    {
      File.seekg(128, ios::beg);
      this->Picture256Colors = new unsigned char* [this->Height];
      int Pos;
      unsigned char Times;  
      for (int j = 0; j < this->Height; j++)
	{
	  Pos = 0;
	  this->Picture256Colors[j] = new unsigned char [this->Width];
	  while (Pos < this->Width)
	    {
	      File.read ((char*) &Times, sizeof(unsigned char));
	      if ((Times >> 6) == 0x03)
		{
		  Times &= 0x3f;
		  File.read ((char*) &Tmp, sizeof(unsigned char));
		  while ((Times-- > 0) && (Pos < this->Width))
		    this->Picture256Colors[j][Pos++] = Tmp;
		}
	      else
		{
		  this->Picture256Colors[j][Pos++] = Times;
		}
	    }
	}
    }
  File.seekg(-769, ios::end);
  File.read ((char*) &Tmp, sizeof(unsigned char));
  this->ColorMap = new PicRGB [256];
  for (int i = 0; i < 256; i++)
    {
      File.read ((char*) &Tmp, sizeof(unsigned char));
      this->ColorMap[i].Red = Tmp;
      File.read ((char*) &Tmp, sizeof(unsigned char));
      this->ColorMap[i].Green = Tmp;
      File.read ((char*) &Tmp, sizeof(unsigned char));
      this->ColorMap[i].Blue = Tmp;
    }
  File.close();
  return true;
}

// clone a bitmap picture type with a new size
// 
// width = new width (0 if old width has to be kept)
// height = new height (0 if old height has to be kept)
// return value = pointer to resulting bitmap picture 

AbstractBitmapPicture* PcxFormat::Clone (int width, int height)
{
  if (width == 0)
    width = this->Width;
  if (height == 0)
    height = this->Height;
  return new PcxFormat (width, height);
}

#ifdef __OPENGL__

// Load bitmap picture as an open gl texture
//
// alpha = alpha factor to use if no alpha channel is defined in bitmap file (if greater than 1.0, use RGB texture only)

void PcxFormat::GLLoadBitmapAsTexture (double alpha)
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
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, this->GLWidth, this->GLHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, this->GLTexture);
    }
  else
    {
      this->GLTexture = new unsigned char [this->GLWidth * this->GLHeight * 4];
      int Index = 0;
      unsigned char AlphaFactor = (unsigned char) (alpha * 255);
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
}

// use bitmap picture as current open gl texture
//

void PcxFormat::GLUseBitmapAsTexture ()
{
  glBindTexture (GL_TEXTURE_2D, this->GLTextureID);
  glScaled(((double) this->Width) / ((double) this->GLWidth), ((double) this->Height) / ((double) this->GLHeight), 0.0);
}

#endif

