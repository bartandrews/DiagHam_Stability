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


#include "BitmapTools/BitmapPicture/IffFormat.h"
#include "BitmapTools/BitmapPicture/IffTags.h"
#include "BitmapTools/Color/PicBGR.h"
#include <stdlib.h>
#include <iostream>


using std::ios;
using std::cout;
using std::endl;
using std::dec;
using std::hex;


// function to convert Big Endian encoded variable into Little Endian encoded variable if __LITTLEENDIAN__ is define
//
// File = reference on input file stream
// Var = reference on variable where result has to be stored

template<class ClassName>
void ReadBigEndian (ifstream& File, ClassName& Var)
{
  File.read ((char*) &Var, sizeof(ClassName));
#ifdef __LITTLEENDIAN__
  ClassName  TmpVar = Var;
  unsigned char* TmpBin1 = (unsigned char*) &Var;
  unsigned char* TmpBin2 = (unsigned char*) &TmpVar;
  int max = sizeof(ClassName) >> 1;
  for (int i = 0; i < max; i++)
    {
      TmpBin1[i] = TmpBin2[sizeof(ClassName) - i -1];
      TmpBin1[sizeof(ClassName) - i -1] = TmpBin2[i];
    }
#endif
};

// default constructor
//

IffFormat::IffFormat ()
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

IffFormat::IffFormat (int L, int H)
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

IffFormat::IffFormat (int L, int H, PicRGB& Col)
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

IffFormat::~IffFormat ()
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

void IffFormat::SetPixel(int x, int y, PicRGB& Col)
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

PicRGB	IffFormat::GetPixel(int x, int y)
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

bool IffFormat::SavePicture (char* FileName)
{
  return false;
}

// Load Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool IffFormat::LoadPicture (char* FileName)
{
  ifstream File;
  File.open(FileName, ios::binary | ios::in);

  //check file format
  unsigned long Tag = this->ReadTag(File);
  if (Tag != FORM_TAG)
    {
      File.close();
      return false;
    }
  long Size;
  ReadBigEndian (File, Size);
  Tag = this->ReadTag(File);
  if (Tag != ILBM_TAG)
    {
      File.close();
      return false;  
    }
  File.seekg(0, ios::end);
  long Size2 = File.tellg();
  Size2 -= 8;
  if (Size != Size2)
    {
      cout << "warning wrong size given by FORM chunk" << endl;
    }
  Size2 += 8;
  File.seekg(12, ios::beg);  
  Tag = this->ReadTag(File);
  ReadBigEndian (File, Size);
  if (Tag != BMHD_TAG)
    {
      File.close();
      return false;  
    }

  this->ProceedHeader(File, Size);
  this->ColorMap = 0;
  while (File.tellg() < Size2)
    {
      Tag = this->ReadTag(File);
      ReadBigEndian(File, Size);
      switch (Tag)
	{
	case CMAP_TAG:
	  {
	    this->ProceedColorMap(File, Size);
	  }
	  break;
	case BODY_TAG:
	  {
	    this->ProceedBody(File, Size);
	    File.close();
	    return true;
	  }
	  break;
	default:
	  {
	    cout << "unkonwn chunck " << hex << Tag <<  dec << endl;
	    File.seekg(Size, ios::cur);
	  }
	  break;
      }
    }
  File.close();
  return true;
}

// unpack a scanline 
//
// file = reference on current input file stream
// scanline = pointer to array where unpacked scanline has to be stored
// scanlineSize = number of byte per scanline

void IffFormat::UnpackScanline (ifstream& file, unsigned char* scanline, int scanlineSize)
{
  int Count = 0;
  char TmpByte;
  unsigned char TmpByte2;
  int Count2;
  while (Count < scanlineSize)
    {
      file.read(&TmpByte, sizeof(char));
      if (TmpByte != -128)
	{
	  if (TmpByte < 0)
	    {	      
	      file.read((char*) &TmpByte2, sizeof(unsigned char));
	      Count2 = 1 -  (int) TmpByte;
	      for (;(Count2 > 0) && (Count < scanlineSize); --Count2)
		{
		  scanline[Count] = TmpByte2;
		  ++Count;
		}
	    }
	  else
	    {
	      for (; (TmpByte >= 0) && (Count < scanlineSize); --TmpByte)
		{
		  file.read((char*) &TmpByte2, sizeof(unsigned char));
		  scanline[Count] = TmpByte2;
		  ++Count;
		}
	    }
	}
    }
}

// clone a bitmap picture type with a new size
// 
// width = new width (0 if old width has to be kept)
// height = new height (0 if old height has to be kept)
// return value = pointer to resulting bitmap picture 

AbstractBitmapPicture* IffFormat::Clone (int width, int height)
{
  if (width == 0)
    width = this->Width;
  if (height == 0)
    height = this->Height;
  return new IffFormat (width, height);
}

// read tag
//
// file = reference on current input file stream
// return value = tag

unsigned long IffFormat::ReadTag(ifstream& file) 
{
  unsigned long Var;
  ReadBigEndian(file, Var);
  return Var;
}
  
// proceed header chunk
//
// file = reference on current input file stream
// chunkSize = chunk size

void IffFormat::ProceedHeader(ifstream& file, unsigned long chunkSize)
{
  unsigned short TmpSize;
  short TmpPos;
  ReadBigEndian(file, TmpSize);
  this->Width = TmpSize;
  ReadBigEndian(file, TmpSize);
  this->Height = TmpSize;

  //unused pixel position
  ReadBigEndian(file, TmpPos);
  ReadBigEndian(file, TmpPos);

  unsigned char TmpNbrBitplanes;
  file.read ((char*) &TmpNbrBitplanes, sizeof(unsigned char));
  this->NbrBitplanes = TmpNbrBitplanes;
  file.read ((char*) &TmpNbrBitplanes, sizeof(unsigned char));
  this->MaskingMode = TmpNbrBitplanes;
  file.read ((char*) &TmpNbrBitplanes, sizeof(unsigned char));
  this->CompressionMode = TmpNbrBitplanes;

  // padding byte
  file.read ((char*) &TmpNbrBitplanes, sizeof(unsigned char));
  // transparent color (unused)
  ReadBigEndian(file, TmpSize);
  //pixel aspect ratio (unused)
  file.read ((char*) &TmpNbrBitplanes, sizeof(unsigned char));
  file.read ((char*) &TmpNbrBitplanes, sizeof(unsigned char));
  // page width and height (unused)
  ReadBigEndian(file, TmpSize);
  ReadBigEndian(file, TmpSize);
}

// proceed color map chunk
//
// file = reference on current input file stream
// chunkSize = chunk size

void IffFormat::ProceedColorMap(ifstream& file, unsigned long chunkSize)
{
  this->NbrColor = chunkSize / 3;
  this->ColorMap = new PicRGB [NbrColor];
  unsigned char TmpRed;
  unsigned char TmpGreen;
  unsigned char TmpBlue;
  for (int i = 0; i < NbrColor; ++i)
    {
      file.read((char*) &TmpRed, sizeof(unsigned char));
      this->ColorMap[i].Red = TmpRed;
      file.read((char*) &TmpGreen, sizeof(unsigned char));
      this->ColorMap[i].Green = TmpGreen;
      file.read((char*) &TmpBlue, sizeof(unsigned char));
      this->ColorMap[i].Blue = TmpBlue;
    }
  // padding byte
  if (NbrColor & 1)
    file.read((char*) &TmpRed, sizeof(unsigned char));    
}

// proceed body chunk
//
// file = reference on current input file stream
// chunkSize = chunk size

void IffFormat::ProceedBody(ifstream& file, unsigned long chunkSize)
{
  int TmpWidth = this->Width >> 3;
  int TmpCompWidth = this->Width & 0x7;
  if (this->NbrBitplanes <= 8)
    {
      // indexed picture
      // check if it is a gray scale picture
      if (this->ColorMap == 0)
	{
	  this->NbrColor = 1 << this->NbrBitplanes;
	  this->ColorMap = new PicRGB [NbrColor];
	  unsigned char TmpRed;
	  for (int i = 0; i < NbrColor; ++i)
	    {
	      TmpRed = i << (8 - this->NbrBitplanes);
	      this->ColorMap[i].Red = TmpRed;
	      this->ColorMap[i].Green = TmpRed;
	      this->ColorMap[i].Blue = TmpRed;
	    }
	}
      this->Picture256Colors = new unsigned char* [this->Width];
      for (int i = 0; i < this->Width; ++i)
	{
	  this->Picture256Colors[i] = new unsigned char [this->Height];
	  for (int j = 0; j < this->Height; ++j)
	    this->Picture256Colors[i][j] = 0;
	}

      unsigned char TmpByte;
     
      if (this->CompressionMode == 0)
	{
	  // uncompressed indexed picture
	  for (int j = 0; j < this->Height; ++j)
	    {
	      for (int k = 0; k < this->NbrBitplanes; ++k)
		{
		  int TmpWidth2 = 0;
		  for (int i = 0; i < TmpWidth; ++i)
		    {	    
		      file.read((char*) &TmpByte, sizeof(unsigned char));
		      this->Picture256Colors[TmpWidth2++][j] |= ((TmpByte & 0x80) >> 7) << k;
		      this->Picture256Colors[TmpWidth2++][j] |= ((TmpByte & 0x40) >> 6) << k;
		      this->Picture256Colors[TmpWidth2++][j] |= ((TmpByte & 0x20) >> 5) << k;
		      this->Picture256Colors[TmpWidth2++][j] |= ((TmpByte & 0x10) >> 4) << k;
		      this->Picture256Colors[TmpWidth2++][j] |= ((TmpByte & 0x08) >> 3) << k;
		      this->Picture256Colors[TmpWidth2++][j] |= ((TmpByte & 0x04) >> 2) << k;
		      this->Picture256Colors[TmpWidth2++][j] |= ((TmpByte & 0x02) >> 1) << k;
		      this->Picture256Colors[TmpWidth2++][j] |= (TmpByte & 0x01) << k;
		    }
		  if (TmpCompWidth != 0)
		    {
		      file.read((char*) &TmpByte, sizeof(unsigned char));
		      for (int l = 7; l > (7 - TmpCompWidth); --l)
			this->Picture256Colors[TmpWidth2++][j] |= ((TmpByte & (0x01 << l)) >> l) << k;
		    }
		if ((TmpWidth & 0x1) == 0)
		  file.read((char*) &TmpByte, sizeof(unsigned char));		  
		}
	      if (MaskingMode == 1)
		{
		  for (int i = 0; i < TmpWidth; ++i)
		    {	    
		      file.read((char*) &TmpByte, sizeof(unsigned char));
		    }
		  if (TmpCompWidth != 0)
		    file.read((char*) &TmpByte, sizeof(unsigned char));
		}
	    }
	}
      else
	// compressed indexed picture
	{
	  int TmpScanlineWidth = ((this->Width + 15) >> 3) & 0xfffffffe; 
	  unsigned char* TmpScanline = new unsigned char [TmpScanlineWidth * (2 + this->NbrBitplanes)];
	  unsigned char* TmpScanline2 = TmpScanline;
	  for (int j = 0; j < this->Height; ++j)
	    {
	      TmpScanline = TmpScanline2;
	      if (MaskingMode == 1)
		this->UnpackScanline(file, TmpScanline, (1 + this->NbrBitplanes) * TmpScanlineWidth);
	      else
		this->UnpackScanline(file, TmpScanline, this->NbrBitplanes * TmpScanlineWidth);	      
	      for (int k = 0; k < this->NbrBitplanes; ++k)
		{
		  int TmpWidth2 = 0;
		  int i = 0;
		  for (; i < TmpWidth; i++)
		    {
		      this->Picture256Colors[TmpWidth2++][j] |= ((TmpScanline[i] & 0x80) >> 7) << k;
		      this->Picture256Colors[TmpWidth2++][j] |= ((TmpScanline[i] & 0x40) >> 6) << k;
		      this->Picture256Colors[TmpWidth2++][j] |= ((TmpScanline[i] & 0x20) >> 5) << k;
		      this->Picture256Colors[TmpWidth2++][j] |= ((TmpScanline[i] & 0x10) >> 4) << k;
		      this->Picture256Colors[TmpWidth2++][j] |= ((TmpScanline[i] & 0x08) >> 3) << k;
		      this->Picture256Colors[TmpWidth2++][j] |= ((TmpScanline[i] & 0x04) >> 2) << k;
		      this->Picture256Colors[TmpWidth2++][j] |= ((TmpScanline[i] & 0x02) >> 1) << k;
		      this->Picture256Colors[TmpWidth2++][j] |= (TmpScanline[i] & 0x01) << k;
		    }
		  if (TmpCompWidth != 0)
		    {
		      for (int l = 7; l > (7 - TmpCompWidth); --l)
			this->Picture256Colors[TmpWidth2++][j] |= ((TmpScanline[i] & (0x01 << l)) >> l) << k;
		    }
		  TmpScanline += TmpScanlineWidth;
		}
	    }
	  delete[] TmpScanline2;
	}
      if (chunkSize & 0x1)
	file.read((char*) &TmpByte, sizeof(unsigned char));	
      return;
    }
  if (this->NbrBitplanes == 24)
    {
      this->Picture = new PicRGB* [this->Width];
      for (int i = 0; i < this->Width; i++)
	{
	  this->Picture[i] = new PicRGB [this->Height];
	  for (int j = 0; j < this->Height; ++j)
	    {
	      this->Picture[i][j].Red = 0;
	      this->Picture[i][j].Green = 0;
	      this->Picture[i][j].Blue = 0;
	    }
	}
      if (this->CompressionMode == 0)
	{
	  unsigned char TmpByte;
	  for (int j = 0; j < this->Height; ++j)
	    {
	      for (int k = 0; k < 8; ++k)
		{
		  int TmpWidth2 = 0;
		  for (int i = 0; i < TmpWidth; ++i)
		    {
		      file.read((char*) &TmpByte, sizeof(unsigned char));
		      this->Picture[TmpWidth2++][j].Red |= ((TmpByte & 0x80) >> 7) << k;
		      this->Picture[TmpWidth2++][j].Red |= ((TmpByte & 0x40) >> 6) << k;
		      this->Picture[TmpWidth2++][j].Red |= ((TmpByte & 0x20) >> 5) << k;
		      this->Picture[TmpWidth2++][j].Red |= ((TmpByte & 0x10) >> 4) << k;
		      this->Picture[TmpWidth2++][j].Red |= ((TmpByte & 0x08) >> 3) << k;
		      this->Picture[TmpWidth2++][j].Red |= ((TmpByte & 0x04) >> 2) << k;
		      this->Picture[TmpWidth2++][j].Red |= ((TmpByte & 0x02) >> 1) << k;
		      this->Picture[TmpWidth2++][j].Red |= (TmpByte & 0x01) << k;
		    }
		  if (TmpCompWidth != 0)
		    {
		      file.read((char*) &TmpByte, sizeof(unsigned char));
		      for (int l = 7; l > (7 - TmpCompWidth); --l)
			this->Picture[TmpWidth2++][j].Red |= ((TmpByte & (0x01 << l)) >> l) << k;
		    }
		}
	      for (int k = 0; k < 8; ++k)
		{
		  int TmpWidth2 = 0;
		  for (int i = 0; i < TmpWidth; ++i)
		    {
		      file.read((char*) &TmpByte, sizeof(unsigned char));
		      this->Picture[TmpWidth2++][j].Green |= ((TmpByte & 0x80) >> 7) << k;
		      this->Picture[TmpWidth2++][j].Green |= ((TmpByte & 0x40) >> 6) << k;
		      this->Picture[TmpWidth2++][j].Green |= ((TmpByte & 0x20) >> 5) << k;
		      this->Picture[TmpWidth2++][j].Green |= ((TmpByte & 0x10) >> 4) << k;
		      this->Picture[TmpWidth2++][j].Green |= ((TmpByte & 0x08) >> 3) << k;
		      this->Picture[TmpWidth2++][j].Green |= ((TmpByte & 0x04) >> 2) << k;
		      this->Picture[TmpWidth2++][j].Green |= ((TmpByte & 0x02) >> 1) << k;
		      this->Picture[TmpWidth2++][j].Green |= (TmpByte & 0x01) << k;
		    }
		  if (TmpCompWidth != 0)
		    {
		      file.read((char*) &TmpByte, sizeof(unsigned char));
		      for (int l = 7; l > (7 - TmpCompWidth); --l)
			this->Picture[TmpWidth2++][j].Green |= ((TmpByte & (0x01 << l)) >> l) << k;
		    }
		}
	      for (int k = 0; k < 8; ++k)
		{
		  int TmpWidth2 = 0;
		  for (int i = 0; i < TmpWidth; ++i)
		    {
		      file.read((char*) &TmpByte, sizeof(unsigned char));
		      this->Picture[TmpWidth2++][j].Blue |= ((TmpByte & 0x80) >> 7) << k;
		      this->Picture[TmpWidth2++][j].Blue |= ((TmpByte & 0x40) >> 6) << k;
		      this->Picture[TmpWidth2++][j].Blue |= ((TmpByte & 0x20) >> 5) << k;
		      this->Picture[TmpWidth2++][j].Blue |= ((TmpByte & 0x10) >> 4) << k;
		      this->Picture[TmpWidth2++][j].Blue |= ((TmpByte & 0x08) >> 3) << k;
		      this->Picture[TmpWidth2++][j].Blue |= ((TmpByte & 0x04) >> 2) << k;
		      this->Picture[TmpWidth2++][j].Blue |= ((TmpByte & 0x02) >> 1) << k;
		      this->Picture[TmpWidth2++][j].Blue |= (TmpByte & 0x01) << k;
		    }
		  if (TmpCompWidth != 0)
		    {
		      file.read((char*) &TmpByte, sizeof(unsigned char));
		      for (int l = 7; l > (7 - TmpCompWidth); --l)
			this->Picture[TmpWidth2++][j].Blue |= ((TmpByte & (0x01 << l)) >> l) << k;
		    }
		}
	      if (MaskingMode == 1)
		{
		  for (int i = 0; i < TmpWidth; ++i)
		    {	    
		      file.read((char*) &TmpByte, sizeof(unsigned char));
		    }
		  if (TmpCompWidth != 0)
		    file.read((char*) &TmpByte, sizeof(unsigned char));
		}
	    }
	}
      else
	{
	  int TmpScanlineWidth = ((this->Width + 15) >> 3) & 0xfffffffe;
	  unsigned char* TmpScanline = new unsigned char [TmpScanlineWidth * 27];
	  unsigned char* TmpScanline2 = TmpScanline;
	  for (int j = 0; j < this->Height; ++j)
	    {
	      TmpScanline = TmpScanline2;
	      if (MaskingMode == 1)
		this->UnpackScanline(file, TmpScanline, 25 * TmpScanlineWidth);
	      else
		this->UnpackScanline(file, TmpScanline, 24 * TmpScanlineWidth);	      
	      for (int k = 0; k < 8; ++k)
		{
		  int TmpWidth2 = 0;
		  int i = 0;
		  for (; i < TmpWidth; i++)
		    {
		      this->Picture[TmpWidth2++][j].Red |= ((TmpScanline[i] & 0x80) >> 7) << k;
		      this->Picture[TmpWidth2++][j].Red |= ((TmpScanline[i] & 0x40) >> 6) << k;
		      this->Picture[TmpWidth2++][j].Red |= ((TmpScanline[i] & 0x20) >> 5) << k;
		      this->Picture[TmpWidth2++][j].Red |= ((TmpScanline[i] & 0x10) >> 4) << k;
		      this->Picture[TmpWidth2++][j].Red |= ((TmpScanline[i] & 0x08) >> 3) << k;
		      this->Picture[TmpWidth2++][j].Red |= ((TmpScanline[i] & 0x04) >> 2) << k;
		      this->Picture[TmpWidth2++][j].Red |= ((TmpScanline[i] & 0x02) >> 1) << k;
		      this->Picture[TmpWidth2++][j].Red |= (TmpScanline[i] & 0x01) << k;
		    }
		  if (TmpCompWidth != 0)
		    {
		      for (int l = 7; l > (7 - TmpCompWidth); --l)
			this->Picture[TmpWidth2++][j].Red |= ((TmpScanline[i] & (0x01 << l)) >> l) << k;
		    }
		  TmpScanline += TmpScanlineWidth;
		}
	      for (int k = 0; k < 8; ++k)
		{
		  int TmpWidth2 = 0;
		  int i = 0;;
		  for (; i < TmpWidth; i++)
		    {
		      this->Picture[TmpWidth2++][j].Green |= ((TmpScanline[i] & 0x80) >> 7) << k;
		      this->Picture[TmpWidth2++][j].Green |= ((TmpScanline[i] & 0x40) >> 6) << k;
		      this->Picture[TmpWidth2++][j].Green |= ((TmpScanline[i] & 0x20) >> 5) << k;
		      this->Picture[TmpWidth2++][j].Green |= ((TmpScanline[i] & 0x10) >> 4) << k;
		      this->Picture[TmpWidth2++][j].Green |= ((TmpScanline[i] & 0x08) >> 3) << k;
		      this->Picture[TmpWidth2++][j].Green |= ((TmpScanline[i] & 0x04) >> 2) << k;
		      this->Picture[TmpWidth2++][j].Green |= ((TmpScanline[i] & 0x02) >> 1) << k;
		      this->Picture[TmpWidth2++][j].Green |= (TmpScanline[i] & 0x01) << k;
		    }
		  if (TmpCompWidth != 0)
		    {
		      for (int l = 7; l > (7 - TmpCompWidth); --l)
			this->Picture[TmpWidth2++][j].Green |= ((TmpScanline[i] & (0x01 << l)) >> l) << k;
		    }
		  TmpScanline += TmpScanlineWidth;
		}
	      for (int k = 0; k < 8; ++k)
		{
		  int TmpWidth2 = 0;
		  int i = 0;
		  for (; i < TmpWidth; i++)
		    {
		      this->Picture[TmpWidth2++][j].Blue |= ((TmpScanline[i] & 0x80) >> 7) << k;
		      this->Picture[TmpWidth2++][j].Blue |= ((TmpScanline[i] & 0x40) >> 6) << k;
		      this->Picture[TmpWidth2++][j].Blue |= ((TmpScanline[i] & 0x20) >> 5) << k;
		      this->Picture[TmpWidth2++][j].Blue |= ((TmpScanline[i] & 0x10) >> 4) << k;
		      this->Picture[TmpWidth2++][j].Blue |= ((TmpScanline[i] & 0x08) >> 3) << k;
		      this->Picture[TmpWidth2++][j].Blue |= ((TmpScanline[i] & 0x04) >> 2) << k;
		      this->Picture[TmpWidth2++][j].Blue |= ((TmpScanline[i] & 0x02) >> 1) << k;
		      this->Picture[TmpWidth2++][j].Blue |= (TmpScanline[i] & 0x01) << k;
		    }
		  if (TmpCompWidth != 0)
		    {
		      for (int l = 7; l > (7 - TmpCompWidth); --l)
			this->Picture[TmpWidth2++][j].Blue |= ((TmpScanline[i] & (0x01 << l)) >> l) << k;
		    }
		  TmpScanline += TmpScanlineWidth;
		}
	    }
	  delete[] TmpScanline2;
	}
      if (chunkSize & 0x1)
	{
	  unsigned char TmpByte;
	  file.read((char*) &TmpByte, sizeof(unsigned char));	
	}
      return;
    }
}

#ifdef __OPENGL__

// Load bitmap picture as an open gl texture
//
// alpha = alpha factor to use if no alpha channel is defined in bitmap file (if greater than 1.0, use RGB texture only)

void IffFormat::GLLoadBitmapAsTexture (double alpha)
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

void IffFormat::GLUseBitmapAsTexture ()
{
  glBindTexture (GL_TEXTURE_2D, this->GLTextureID);
  glScaled(((double) this->Width) / ((double) this->GLWidth), ((double) this->Height) / ((double) this->GLHeight), 0.0);
}

#endif

