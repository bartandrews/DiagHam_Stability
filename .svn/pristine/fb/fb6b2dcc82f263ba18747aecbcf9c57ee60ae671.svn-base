////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of Bitmap Picture Using Bmp Format 24 bits              //
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


#include "BitmapTools/BitmapPicture/BmpFormat.h"
#include "BitmapTools/Color/PicBGR.h"
#include <stdlib.h>


using std::ios;


// default constructor
//

BmpFormat::BmpFormat ()
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

BmpFormat::BmpFormat (int L, int H)
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

BmpFormat::BmpFormat (int L, int H, PicRGB& Col)
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

BmpFormat::~BmpFormat ()
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

void BmpFormat::SetPixel(int x, int y, PicRGB& Col)
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

PicRGB	BmpFormat::GetPixel(int x, int y)
{
  if ((x < 0) || (x >= this->Width) || (y < 0) || (y >= this->Height))
    return PicRGB(0, 0, 0);
  if (this->Picture != 0)
    return this->Picture [x][y];
  else
    return this->ColorMap[this->Picture256Colors [x][y]];
}

// Save Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool BmpFormat::SavePicture (char* FileName)
{
  if (this->Picture == 0)
    return false;
  ofstream File;
  File.open(FileName, ios::binary | ios::out);
  char c;
  c = 'B';
  File.write ((char*) (&c), sizeof(char));
  c = 'M';
  File.write ((char*) (&c), sizeof(char));
  int mod = 4 - (3 * this->Width) % 4;
  if (mod == 4)
    mod = 0;
  unsigned long PicSize = (3 * this->Width + mod) * this->Height;
  unsigned long FileSize = 54 + PicSize;
  unsigned long Offset = 54;
  unsigned long Size = 40;
  unsigned long Zero = 0;
  unsigned char ZeroC = 0;
  File.write ((char*) (&FileSize), sizeof(unsigned long));
  File.write ((char*) (&Zero), sizeof(unsigned long));
  File.write ((char*) (&Offset), sizeof(unsigned long));
  File.write ((char*) (&Size), sizeof(unsigned long));
  File.write ((char*) (&this->Width), sizeof(long));
  File.write ((char*) (&this->Height), sizeof(long));
  unsigned short Plane = 1;
  File.write ((char*) (&Plane), sizeof(unsigned short));
  Plane = 24;
  File.write ((char*) (&Plane), sizeof(unsigned short));
  File.write ((char*) (&Zero), sizeof(unsigned long));
  File.write ((char*) (&PicSize), sizeof(unsigned long));
  File.write ((char*) (&Zero), sizeof(unsigned long));
  File.write ((char*) (&Zero), sizeof(unsigned long));
  File.write ((char*) (&Zero), sizeof(unsigned long));
  File.write ((char*) (&Zero), sizeof(unsigned long));
  for (int j = this->Height - 1; j >= 0; j--)
    {
      PicBGR TmpBGR;
      for (int i = 0; i < this->Width; i++)
	{
	  TmpBGR = this->Picture[i][j];
	  File << TmpBGR;
	}
      int i = 0;
      while (i++ < mod)
	File.write ((char*) (&ZeroC), sizeof(unsigned char));
    }
  File.close();
  return true;
}

// Load Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool BmpFormat::LoadPicture (char* FileName)
{
  ifstream File;
  char c;
  File.open(FileName, ios::binary | ios::in);
  File >> c;
  if (c != 'B')
    {
      File.close();
      return false;
    }
  File >> c;
  if (c != 'M')
    {
      File.close();
      return false;
    }
  unsigned long FileSize;
  unsigned long Offset;
  File.read ((char*) (&FileSize), sizeof(unsigned long));
  File.seekg(4, ios::cur);
  File.read ((char*) (&Offset), sizeof(unsigned long));
  
  unsigned long Size, Compression;
  long TmpWidth, TmpHeight;
  unsigned short Plane, BitCount;
  File.read ((char*) (&Size), sizeof(unsigned long));
  File.read ((char*) (&TmpWidth), sizeof(long));
  File.read ((char*) (&TmpHeight), sizeof(long));
  File.read ((char*) (&Plane), sizeof(unsigned short));
  if (Plane != 1)
    {
      File.close();
      return false;
    }
  File.read ((char*) (&BitCount), sizeof(unsigned short));
  if ((BitCount != 24) && (BitCount != 8))
    {
       File.close();
       return false;
    }
  File.read ((char*) (&Compression), sizeof(unsigned long));
  if ((BitCount == 24) && (Compression != 0))
    {
      File.close();
      return false;
    }
  if (BitCount == 24)
    {
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
      File.seekg(Offset, ios::beg);
      PicBGR Col;
      int mod = 4 - (3 * this->Width) % 4;
      if (mod == 4)
	 mod = 0;
      for (int j = this->Height -1; j >= 0; j--)
	{
	   for (int i = 0; i < this->Width; i++)
	     {
	       File >> Col;
	       this->Picture[i][j] = Col;
	     }
	   File.seekg(mod, ios::cur);
	   }
    }
   else
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
	   this->ColorMap = new PicRGB [256];
	 }
       else
	 this->ColorMap = new PicRGB [256];
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
       File.seekg(Size + 14, ios::beg);
       for (int i = 0; i < 256; i++)
	 {
	   File >> Col;
	   File.seekg(1, ios::cur);
	   this->ColorMap[i] = Col;
	 }
       if (Compression == 0)
	 {
	   File.seekg(Offset, ios::beg);
	   int mod = 4 - (this->Width) % 4;
	   if (mod == 4)
	     mod = 0;
	   for (int j = this->Height -1; j >= 0; j--)
	     {
	       for (int i = 0; i < this->Width; i++)
		 File.get((char&) Picture256Colors[i][j]);
	       File.seekg(mod, ios::cur);
	     }
	 }
       else
	 {
	   File.seekg(Offset, ios::beg);
	   unsigned char RepetitionCountField;
	   unsigned char ColIndex;
	   bool EndBitmap = false;
	   int j = this->Height -1;
	   while (EndBitmap == false)
	     {
	       bool EndLine = false;
	       int i = 0;
	       while (EndLine == false)
		  {
		    File.get((char&) RepetitionCountField);
		    if (RepetitionCountField != 0)
		      {
			File.get((char&) ColIndex);
	            	for (int k = 0; k < (int) RepetitionCountField; k++)
			  this->Picture256Colors[i++][j] = ColIndex;
		      }
		    else
		      {
			File.get((char&) RepetitionCountField);
			switch ((int) RepetitionCountField)
			  {
	                  case 0:
			    {
			      EndLine = true;
			      j--;
			      break;
			    }
			  case 1:
			    {
			      EndLine = true;
			      EndBitmap = true;
			      break;
			    }
			  case 2:
			    {
			      File.get((char&) RepetitionCountField);
			      i += (int) RepetitionCountField;
			      File.get((char&) RepetitionCountField);
			      j -= (int) RepetitionCountField;
			      break;
			    }
			  default:
			    {
			      for (int k = 0; k < (int) RepetitionCountField; k++)
				{
				  File.get((char&) ColIndex);
				  this->Picture256Colors[i++][j] = ColIndex;
				}
			      if ((((int) RepetitionCountField) % 2) != 0)
				File.seekg(1, ios::cur);
			    }
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

AbstractBitmapPicture* BmpFormat::Clone (int width, int height)
{
  if (width == 0)
    width = this->Width;
  if (height == 0)
    height = this->Height;
  return new BmpFormat (width, height);
}

#ifdef __OPENGL__

// Load bitmap picture as an open gl texture
//
// alpha = alpha factor to use if no alpha channel is defined in bitmap file (if greater than 1.0, use RGB texture only)

void BmpFormat::GLLoadBitmapAsTexture (double alpha)
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
		  this->GLTexture[Index++] = this->Picture[i][j].Red;
		  this->GLTexture[Index++] = this->Picture[i][j].Green;
		  this->GLTexture[Index++] = this->Picture[i][j].Blue;
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
	  for (; j < this->Height; ++j)
	    {
	      int i = 0;
	      for (; i < this->Width; i++)
		{
		  this->GLTexture[Index++] = this->Picture[i][j].Red;
		  this->GLTexture[Index++] = this->Picture[i][j].Green;
		  this->GLTexture[Index++] = this->Picture[i][j].Blue;
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

void BmpFormat::GLUseBitmapAsTexture ()
{
  glBindTexture (GL_TEXTURE_2D, this->GLTextureID);
  glScaled(((double) this->Width) / ((double) this->GLWidth), ((double) this->Height) / ((double) this->GLHeight), 0.0);
}

#endif

