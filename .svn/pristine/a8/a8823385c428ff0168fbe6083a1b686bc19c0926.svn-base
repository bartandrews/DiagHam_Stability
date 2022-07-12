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


#include "BitmapPicture/TgaFormat.h"
#include "Color/PicBGR.h"
#include <fstream>
#include <stdlib.h>


using std::ios;


// default constructor
// 

TgaFormat::TgaFormat ()
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
TgaFormat::TgaFormat (int L, int H)
{
  this->Height = H;
  this->Width = L;
  this->Picture = new PicRGB* [this->Width];
  for (int i = 0; i < Width; i++)
    this->Picture[i] = new PicRGB [this->Height];
  this->Picture256Colors = 0;
  this->ColorMap = 0;
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

TgaFormat::TgaFormat (int L, int H, PicRGB& Col)
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

TgaFormat::~TgaFormat ()
{
  if (this->Picture != 0)
    {
      for (int i = 0; i < this->Width; i++)
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

void TgaFormat::SetPixel(int x, int y, PicRGB& Col)
{
  if (this->Picture != 0)
    this->Picture [x][y] = Col;
  else
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
      this->Picture256Colors [x][y] = (unsigned char) Index;
   }
}

// Get Pixel Value
//
// x = pixel coordinate along x axis
// y = pixel coordinate along y axis
// reurn value = pixel color

PicRGB TgaFormat::GetPixel(int x, int y)
{
  if (this->Picture != 0)
    return this->Picture [x][y];
  else
    return this->ColorMap[this->Picture256Colors [x][y]];
}

// Save Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool TgaFormat::SavePicture (char* FileName)
{
  return this->SavePicture (FileName, 0x02);
}

// Save Picture
// 
// fileName = picture filename
// mode = TGA mode to use
// return value = true if operation is successfull

bool TgaFormat::SavePicture (char* FileName, int mode)
{
  if (this->Picture == 0)
    return false;
  ofstream File;
  File.open(FileName, ios::binary | ios::out);
  unsigned char IDLenght = 0;
  unsigned char ColorMapType = 0;
  unsigned char ImageType = (unsigned char) (mode & 0x0f);
  File.put(IDLenght);
  File.put(ColorMapType);
  File.put(ImageType);
  unsigned short FirstEntryIndex = 0;
  unsigned short ColorMapLength = 0;
  unsigned char ColorMapEntrySize = 0;
  File.write ((char*) (&FirstEntryIndex), sizeof(unsigned short));
  File.write ((char*) (&ColorMapLength), sizeof(unsigned short));
  File.put(ColorMapEntrySize);
  unsigned short XOriginImage = 0;
  unsigned short YOriginImage = 0;
  unsigned short Width = (unsigned short) this->Width;
  unsigned short Height = (unsigned short) this->Height;
  unsigned char PixelDepth = 24;
  unsigned char ImageDescriptor = 0;
  File.write ((char*) (&XOriginImage), sizeof(unsigned short));
  File.write ((char*) (&YOriginImage), sizeof(unsigned short));
  File.write ((char*) (&Width), sizeof(unsigned short));
  File.write ((char*) (&Height), sizeof(unsigned short));
  File.put(PixelDepth);
  File.put(ImageDescriptor);
  unsigned char Order = (unsigned char) ((mode & 0x0f) / 0x10);
  unsigned char AlphaChannel = (unsigned char) (ImageDescriptor & 31);

  if (ImageType == TgaFormat::UnPacked)
    {
      PicBGR TmpBGR;
      for (int j = this->Height -1; j >= 0; j--)
	{
	  for (int i = 0; i < this->Width; i++)
	    {
	      TmpBGR = this->Picture[i][j];
	      File << TmpBGR;
	    }
	}
    }
  else
    {
      PicBGR TmpBGR;
      for (int j = this->Height -1; j >= 0; j--)
	{
	  int i = 0;
	  while (i < this->Width)
	    {
	      int k = i;
	      unsigned char RepetitionCountField = 0;
	      while ((k < (this->Width-1)) && (this->Picture[k][j] != this->Picture[k+1][j]) && (RepetitionCountField != 127))
		{
		  k++;
		   RepetitionCountField = (unsigned char) (RepetitionCountField + 1);
		}
	       if (k == (this->Width - 1))
		 {
		   File.put(RepetitionCountField);
		   TmpBGR = this->Picture[i++][j]; 
		   File << TmpBGR;
		   while (RepetitionCountField != 0)
		     {
		       TmpBGR = this->Picture[i++][j]; 
		       File << TmpBGR;
	               RepetitionCountField = (unsigned char) (RepetitionCountField - 1);
		     }
		 }
	       if ((RepetitionCountField == 127) && (k < (this->Width-1)) && (this->Picture[k][j] != this->Picture[k+1][j]))
		 {
		   File.put(RepetitionCountField);
		   while (RepetitionCountField != 0)
		     {
		       TmpBGR = this->Picture[i++][j]; 
		       File << TmpBGR;
	               RepetitionCountField = (unsigned char) (RepetitionCountField - 1);
		     }
		 }
	       if ((k < (this->Width-1)) && (this->Picture[k][j] == this->Picture[k+1][j]))
		 {
		   if (RepetitionCountField != 0)
		     {
		       RepetitionCountField = (unsigned char) (RepetitionCountField - 1);
		       File.put(RepetitionCountField);
		       RepetitionCountField = (unsigned char) (RepetitionCountField + 1);
		       while (RepetitionCountField != 0)
			 {
			   TmpBGR = this->Picture[i++][j]; 
			   File << TmpBGR;
			   RepetitionCountField = (unsigned char) (RepetitionCountField - 1);
		   }
		     }
		   RepetitionCountField = 0;
		   while ((k < (this->Width-1)) && (this->Picture[k][j] == this->Picture[k+1][j]) && (RepetitionCountField != 127))
		     {
		       k++;
	               RepetitionCountField = (unsigned char) (RepetitionCountField + 1);
		     }
		   if (k == (this->Width - 1))
		     {
		       RepetitionCountField = (unsigned char) (RepetitionCountField + 128);
		       File.put(RepetitionCountField);
		       TmpBGR = this->Picture[i][j]; 
		       File << TmpBGR;
		       i += RepetitionCountField - 127;
		     }
		   if ((RepetitionCountField == 127) && (k < (this->Width-1)) && (this->Picture[k][j] == this->Picture[k+1][j]))
		     {
		       RepetitionCountField = (unsigned char) (RepetitionCountField + 128);
		       File.put(RepetitionCountField);
		       TmpBGR = this->Picture[i][j]; 
		       File << TmpBGR;
		       i += RepetitionCountField - 127;
		     }
		   if ((k < (this->Width-1)) && (this->Picture[k][j] != this->Picture[k+1][j]))
		     {
		       RepetitionCountField = (unsigned char) (RepetitionCountField + 128);
		       File.put(RepetitionCountField);
		       TmpBGR = this->Picture[i][j]; 
		       File << TmpBGR;
		       i += RepetitionCountField - 127;
		     }
		 }
	    }
	}
    }
  File.close();
  return true;
}

// Load Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool TgaFormat::LoadPicture (char* FileName)
{
  ifstream File;
  File.open(FileName, ios::binary | ios::in);
  unsigned char IDLenght;
  unsigned char ColorMapType;
  unsigned char ImageType;
  File.get((char&) IDLenght);
  File.get((char&) ColorMapType);
  File.get((char&) ImageType);
  if ((ImageType != TgaFormat::UnPacked) && (ImageType != TgaFormat::Packed)
      && (ImageType != TgaFormat::UnPackedColorMapped) && (ImageType != TgaFormat::PackedColorMapped))
    {
      File.close();
      return false;
    }
  unsigned short FirstEntryIndex;
  unsigned short ColorMapLength;
  unsigned char ColorMapEntrySize;
  File.read ((char*) (&FirstEntryIndex), sizeof(unsigned short));
  File.read ((char*) (&ColorMapLength), sizeof(unsigned short));
  File.get((char&) ColorMapEntrySize);
  this->NbrColor = (int) ColorMapLength - (int) FirstEntryIndex;
  if ((ColorMapType == 1) && (ColorMapEntrySize != 24))
    {
      File.close();
      return false;
    }
  unsigned short XOriginImage;
  unsigned short YOriginImage;
  unsigned short Width;
  unsigned short Height;
  unsigned char PixelDepth;
  unsigned char ImageDescriptor;
  File.read ((char*) (&XOriginImage), sizeof(unsigned short));
  File.read ((char*) (&YOriginImage), sizeof(unsigned short));
  File.read ((char*) (&Width), sizeof(unsigned short));
  File.read ((char*) (&Height), sizeof(unsigned short));
  File.get((char&) PixelDepth);
  if (((PixelDepth != 24) && (ColorMapType ==0)) || ((PixelDepth != 8) && (ColorMapType ==1)))
    {
      File.close();
      return false;
    }
  long TmpWidth = Width;
  long TmpHeight = Height;
  File.get((char&) ImageDescriptor);
  unsigned char Order = (unsigned char) ((ImageDescriptor & 48) / 16);
  unsigned char AlphaChannel = (unsigned char) (ImageDescriptor & 31);
  
  this->ImageIDLenght = IDLenght;
  if (IDLenght != 0)
    {
      this->ImageID = new char [IDLenght];
      File.read (this->ImageID, this->ImageIDLenght);
    }
  else
    this->ImageID = 0;
  if (ColorMapType == 1)
    {
      if (this->Picture != 0)
	{
	  for (int i = 0; i < this->Width; i++)
	    delete[] this->Picture[i];
	  delete[] this->Picture;
	  this->Picture = 0;
	}
      if (this->ColorMap != 0)
	{
	  delete[] this->ColorMap;
	  this->ColorMap = 0;
	}
      else
      	this->ColorMap = new PicRGB [this->NbrColor];
      if (this->Picture256Colors == 0)
	{
	  this->Height = TmpHeight;
	  this->Width = TmpWidth;
	  this->Picture256Colors = new unsigned char* [this->Width];
	  for (int i = 0; i < this->Width; i++)
	    this->Picture256Colors[i] = new unsigned char [this->Height];
	}
      else
	if ((this->Height != TmpHeight) || (this->Width != TmpWidth))
	  {
	    for (int i = 0; i < this->Width; i++)
	      delete[] this->Picture[i];
	      delete[] this->Picture256Colors;
	      this->Height = TmpHeight;
	      this->Width = TmpWidth;
	      this->Picture256Colors = new unsigned char* [this->Width];
	      for (int i = 0; i < this->Width; i++)
		this->Picture256Colors[i] = new unsigned char [this->Height];
	  }
      PicBGR Col;
      File.seekg((int) FirstEntryIndex, ios::cur);
      for (int i = 0; i < this->NbrColor; i++)
	{
	  File >> Col;
	  this->ColorMap[i] = Col;
	}
      if (ImageType == TgaFormat::UnPackedColorMapped)
	{
	  for (int j = this->Height -1; j >= 0; j--)
	    for (int i = 0; i < this->Width; i++)
	      File.get((char&) Picture256Colors[i][j]);
	}
      else
	{
	  unsigned char RepetitionCountField;
	  unsigned char ColIndex;
	  int NbrPixel;
	  for (int j = this->Height -1; j >= 0; j--)
	    {
	      int i = 0;
	      while (i < this->Width)
	      	{
		  File.get((char&) RepetitionCountField);
		  NbrPixel = (RepetitionCountField & 127) + 1;
		  if (RepetitionCountField & 128)
	            {
		      File.get((char&) ColIndex);
		      for (int k = 0; ((k < NbrPixel) && (i < this->Width)); k++)
			this->Picture256Colors[i++][j] = ColIndex;
	            }
		  else
	            {
		      for (int k = 0; ((k < NbrPixel) && (i < this->Width)); k++)
			{
			  File.get((char&) ColIndex);
			  this->Picture256Colors[i++][j] = ColIndex;
			}
	            }
	      	}
	    }
	}
    }
  else
    {
      this->NbrColor = 0;
      if (this->Picture256Colors != 0)
	{
	  for (int i = 0; i < this->Width; i++)
	    delete[] this->Picture256Colors[i];
	  delete[] this->Picture256Colors;
	  this->Picture256Colors = 0;
	}
      if (this->ColorMap != 0)
	{
	  delete[] this->ColorMap;
	  this->ColorMap = 0;
	}
      if (this->Picture == 0)
	{
	  this->Height = TmpHeight;
	  this->Width = TmpWidth;
	  this->Picture = new PicRGB* [this->Width];
	  for (int i = 0; i < this->Width; i++)
	    this->Picture[i] = new PicRGB [this->Height];
	}
      else
	if ((this->Height != TmpHeight) || (this->Width != TmpWidth))
	  {
	    for (int i = 0; i < this->Width; i++)
	      delete[] this->Picture[i];
	    delete[] this->Picture;
	    this->Height = TmpHeight;
	    this->Width = TmpWidth;
	    this->Picture = new PicRGB* [this->Width];
	    for (int i = 0; i < this->Width; i++)
	      this->Picture[i] = new PicRGB [this->Height];
	   }
      PicBGR Col;
      if (ImageType == TgaFormat::UnPacked)
	{
	  for (int j = this->Height -1; j >= 0; j--)
	    {
	      for (int i = 0; i < this->Width; i++)
	      	{
		  File >> Col;
		  this->Picture[i][j] = Col;
	      	}
	    }
	}
      else
	{
	  unsigned char RepetitionCountField;
	  int NbrPixel;
	  for (int j = this->Height -1; j >= 0; j--)
	    {
	      int i = 0;
	      while (i < this->Width)
	      	{
		  File.get((char&) RepetitionCountField);
		  NbrPixel = (RepetitionCountField & 127) + 1;
		  if (RepetitionCountField & 128)
	            {
		      File >> Col;
		      for (int k = 0; ((k < NbrPixel) && (i < this->Width)); k++)
			this->Picture[i++][j] = Col;
	            }
		  else
	            {
		      for (int k = 0; ((k < NbrPixel) && (i < this->Width)); k++)
			{
			  File >> Col;
			  this->Picture[i++][j] = Col;
			}
	            }
	      	}
	    }
	}
    }
  File.close();
  return true;
}

// clone a bitmap picture type with a new size
// 
// width = new width (0 if old width has to be kept)
// height = new height (0 if old height has to be kept)
// return value = pointer to resulting bitmap picture 

AbstractBitmapPicture* TgaFormat::Clone (int width, int height)
{
  if (width == 0)
    width = this->Width;
  if (height == 0)
    height = this->Height;
  return new TgaFormat (width, height);
}

#ifdef __OPENGL__

// Load bitmap picture as an open gl texture
//
// alpha = alpha factor to use if no alpha channel is defined in bitmap file (if greater than 1.0, use RGB texture only)
  
void TgaFormat::GLLoadBitmapAsTexture (double alpha)
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

// use bitmap picture as current open gl texture
//

void TgaFormat::GLUseBitmapAsTexture ()
{
  glBindTexture (GL_TEXTURE_2D, this->GLTextureID);
  glScaled(((double) this->Width) / ((double) this->GLWidth), ((double) this->Height) / ((double) this->GLHeight), 0.0);
}

#endif
