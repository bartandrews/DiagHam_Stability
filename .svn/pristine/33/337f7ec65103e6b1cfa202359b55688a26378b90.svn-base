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


#include "BitmapTools/BitmapPicture/HdrFormat.h"
#include "BitmapTools/Color/PicBGR.h"
#include "GeneralTools/ListIterator.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>


using std::ios;
using std::cout;
using std::endl;


// default constructor
//

HdrFormat::HdrFormat ()
{
  this->Height = 0;
  this->Width = 0;
  this->Picture = 0;
#ifdef __OPENGL__
  this->GLTextureID = 0;
#endif
}

// constructor for an empty black picture
//
// width = picture width
// height = picture height

HdrFormat::HdrFormat (int L, int H)
{
  this->Height = H;
  this->Width = L;
  this->Picture = new RGB* [this->Height];
  for (int i = 0; i < this->Height; i++)
    this->Picture[i] = new RGB [this->Width];
#ifdef __OPENGL__
  this->GLTextureID = 0;
#endif
}

// constructor for an empty picture using user background color 
//
// width = picture width
// height = picture height
// color = background color

HdrFormat::HdrFormat (int L, int H, RGB& Col)
{
  this->Height = H;
  this->Width = L;
  this->Picture = new RGB* [this->Height];
  for (int i = 0; i < this->Height; i++)
    {
      this->Picture[i] = new RGB [this->Width];
      for (int j = 0; j < this->Width; j++)
      this->Picture[i][j] = Col;
    }
#ifdef __OPENGL__
  this->GLTextureID = 0;
#endif
}

// destructor
//

HdrFormat::~HdrFormat ()
{
  if (this->Picture != 0)
    {
      for (int i = 0; i < this->Height; i++)
	delete[] this->Picture[i];
    }  
  char** TmpString;
  ListIterator<char*> CommandIterator(this->Commands);
  while ((TmpString = CommandIterator()))
    delete[] *TmpString;
  CommandIterator.DefineList(this->CommandArgument);
  while ((TmpString = CommandIterator()))
    delete[] *TmpString;

#ifdef __OPENGL__
  if (this->GLTexture != 0)
    delete[] this->GLTexture;
#endif
}

// Set Pixel Value
//
// x = pixel coordinate along x axis
// y = pixel coordinate along y axis

void HdrFormat::SetPixel(int x, int y, PicRGB& Col)
{
  this->Picture [y][x] = Col;
}

// Set Pixel Value
//
// x = pixel coordinate along x axis
// y = pixel coordinate along y axis
// color = pixel color

void HdrFormat::SetPixel(int x, int y, Color& color)
{
  this->Picture [y][x] = color.GetRGBComponent();
}
  
// Get Pixel Value
//
// x = pixel coordinate along x axis
// y = pixel coordinate along y axis
// reurn value = pixel color

PicRGB	HdrFormat::GetPixel(int x, int y)
{
  if ((x < 0) || (x >= this->Width) || (y < 0) || (y >= this->Height))
    return PicRGB(0, 0, 0);
  return this->Picture [y][x];
}

// Get Pixel Value
//
// x = pixel coordinate along x axis
// y = pixel coordinate along y axis
// color = reference on color that will contain pixel color

void HdrFormat::GetPixel(int x, int y, Color& color)
{
  if ((x < 0) || (x >= this->Width) || (y < 0) || (y >= this->Height))
    return;
  Color TmpColor (this->Picture [y][x]);
  color = TmpColor;
}


// Save Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool HdrFormat::SavePicture (char* FileName)
{
  ofstream File;
  File.open(FileName, ios::binary | ios::out);
  File.write("#?RADIANCE\n# Made with DarkRay\nFORMAT=32-bit_rle_rgbe\nEXPOSURE=1.0000000000000\n\n", 80);
  char* TmpString = new char [24];
  sprintf (TmpString, "-Y %d +X %d\n", this->Height, this->Width);
  File.write(TmpString, strlen(TmpString) * sizeof(char));
  delete[] TmpString;

  this->Scanline = new PicRGBE [this->Width];
  RGB* TmpScanline;
  int j = 0;
  for (; j < this->Height; ++j)//1; ++j)//
    {
      TmpScanline = this->Picture[j];
      for (int i = 0; i < this->Width; ++i)
	this->Scanline[i] = TmpScanline[i];
      if ((this->WriteScanline(File)) == false)
	j = this->Height + 1;
    }
  delete[] this->Scanline;
  File.close();
  if (j > this->Height)
    return false;
  return true;
}

// write a scanline in a hdr file 
//
// file = reference on input file
// return value = true if the scanline has been written

bool HdrFormat::WriteScanline (ofstream& file)
{
  unsigned char TmpChar = 2;
  file.write((char*) &TmpChar, sizeof(unsigned char));
  file.write((char*) &TmpChar, sizeof(unsigned char));
  TmpChar = (unsigned char) ((this->Width >> 8) & 0x7f);
  file.write((char*) &TmpChar, sizeof(unsigned char));
  TmpChar = (unsigned char) (this->Width & 0xff);
  file.write((char*) &TmpChar, sizeof(unsigned char)); 
  int Pos = 0;
  int Pos2 = 0;
  unsigned char Length = 0;
  int TotalLength = 0;
  int ReducedWidth = this->Width - 1;
  while (Pos2 < this->Width)
    {
      Length = 1;
      while ((Length < 128) && (Pos2 < ReducedWidth) && (this->Scanline[Pos2].Red != this->Scanline[Pos2 + 1].Red))
	{
	  ++Length;
	  ++Pos2;
	}
      if (Length == 1)
	{
	  while ((Length < 127) && (Pos2 < ReducedWidth) && (this->Scanline[Pos2].Red == this->Scanline[Pos2 + 1].Red))
	    {
	      ++Length;
	      ++Pos2;
	    }
	  ++Pos2;
	  TotalLength += Length;
	  Length |= 0x80;
	  file.write((char*) &Length, sizeof(unsigned char));
	  file.write((char*) &(this->Scanline[Pos].Red), sizeof(unsigned char));
	}
      else
	{
	  if ((Length == 128) || (Pos2 == ReducedWidth))
	    {
	      TotalLength += Length;
	      file.write((char*) &Length, sizeof(unsigned char));
	      for (int j = Pos; j <= Pos2; ++j)
		file.write((char*) &(this->Scanline[j].Red), sizeof(unsigned char));
	      ++Pos2;
	    }
	  else
	    {
	      --Length;
	      TotalLength += Length;
	      file.write((char*) &Length, sizeof(unsigned char));
	      for (int j = Pos; j < Pos2; ++j)
		file.write((char*) &(this->Scanline[j].Red), sizeof(unsigned char));
	    }
	}
      Pos = Pos2;
    }
  Pos = 0;
  Pos2 = 0;
  while (Pos2 < this->Width)
    {
      Length = 1;
      while ((Length < 128) && (Pos2 < ReducedWidth) && (this->Scanline[Pos2].Green != this->Scanline[Pos2 + 1].Green))
	{
	  ++Length;
	  ++Pos2;
	}
      if (Length == 1)
	{
	  while ((Length < 127) && (Pos2 < ReducedWidth) && (this->Scanline[Pos2].Green == this->Scanline[Pos2 + 1].Green))
	    {
	      ++Length;
	      ++Pos2;
	    }
	  ++Pos2;
	  TotalLength += Length;
	  Length |= 0x80;
	  file.write((char*) &Length, sizeof(unsigned char));
	  file.write((char*) &(this->Scanline[Pos].Green), sizeof(unsigned char));
	}
      else
	{
	  if ((Length == 128) || (Pos2 == ReducedWidth))
	    {
	      TotalLength += Length;
	      file.write((char*) &Length, sizeof(unsigned char));
	      for (int j = Pos; j <= Pos2; ++j)
		file.write((char*) &(this->Scanline[j].Green), sizeof(unsigned char));
	      ++Pos2;
	    }
	  else
	    {
	      --Length;
	      TotalLength += Length;
	      file.write((char*) &Length, sizeof(unsigned char));
	      for (int j = Pos; j < Pos2; ++j)
		file.write((char*) &(this->Scanline[j].Green), sizeof(unsigned char));
	    }
	}
      Pos = Pos2;
    }
  Pos = 0;
  Pos2 = 0;
  while (Pos2 < this->Width)
    {
      Length = 1;
      while ((Length < 128) && (Pos2 < ReducedWidth) && (this->Scanline[Pos2].Blue != this->Scanline[Pos2 + 1].Blue))
	{
	  ++Length;
	  ++Pos2;
	}
      if (Length == 1)
	{
	  while ((Length < 127) && (Pos2 < ReducedWidth) && (this->Scanline[Pos2].Blue == this->Scanline[Pos2 + 1].Blue))
	    {
	      ++Length;
	      ++Pos2;
	    }
	  ++Pos2;
	  TotalLength += Length;
	  Length |= 0x80;
	  file.write((char*) &Length, sizeof(unsigned char));
	  file.write((char*) &(this->Scanline[Pos].Blue), sizeof(unsigned char));
	}
      else
	{
	  if ((Length == 128) || (Pos2 == ReducedWidth))
	    {
	      TotalLength += Length;
	      file.write((char*) &Length, sizeof(unsigned char));
	      for (int j = Pos; j <= Pos2; ++j)
		file.write((char*) &(this->Scanline[j].Blue), sizeof(unsigned char));
	      ++Pos2;
	    }
	  else
	    {
	      --Length;
	      TotalLength += Length;
	      file.write((char*) &Length, sizeof(unsigned char));
	      for (int j = Pos; j < Pos2; ++j)
		file.write((char*) &(this->Scanline[j].Blue), sizeof(unsigned char));
	    }
	}
      Pos = Pos2;
    }
  Pos = 0;
  Pos2 = 0;
  while (Pos2 < this->Width)
    {
      Length = 1;
      while ((Length < 128) && (Pos2 < ReducedWidth) && (this->Scanline[Pos2].Exponent != this->Scanline[Pos2 + 1].Exponent))
	{
	  ++Length;
	  ++Pos2;
	}
      if (Length == 1)
	{
	  while ((Length < 127) && (Pos2 < ReducedWidth) && (this->Scanline[Pos2].Exponent == this->Scanline[Pos2 + 1].Exponent))
	    {
	      ++Length;
	      ++Pos2;
	    }
	  ++Pos2;
	  TotalLength += Length;
	  Length |= 0x80;
	  file.write((char*) &Length, sizeof(unsigned char));
	  file.write((char*) &(this->Scanline[Pos].Exponent), sizeof(unsigned char));
	}
      else
	{
	  if ((Length == 128) || (Pos2 == ReducedWidth))
	    {
	      TotalLength += Length;
	      file.write((char*) &Length, sizeof(unsigned char));
	      for (int j = Pos; j <= Pos2; ++j)
		file.write((char*) &(this->Scanline[j].Exponent), sizeof(unsigned char));
	      ++Pos2;
	    }
	  else
	    {
	      --Length;
	      TotalLength += Length;
	      file.write((char*) &Length, sizeof(unsigned char));
	      for (int j = Pos; j < Pos2; ++j)
		file.write((char*) &(this->Scanline[j].Exponent), sizeof(unsigned char));
	    }
	}
      Pos = Pos2;
    }
  return true;
}

// Load Picture
// 
// fileName = picture filename
// return value = true if operation is successfull

bool HdrFormat::LoadPicture (char* FileName)
{
  ifstream File;
  int FileSize = 0;
  File.open(FileName, ios::binary | ios::in);
  File.seekg(0, ios::end);
  FileSize = File.tellg();
  File.seekg(0, ios::beg);
  if (FileSize < 10)
    {
      File.close();
      return false;
    }  
  char* TestString = new char [16];
  File.read (TestString, sizeof(char) * 11);
  if (strncmp(TestString, "#?RADIANCE\n", 11) != 0)
    {
      delete[] TestString;
      return false;
    }
  delete[] TestString;
  char* TmpString = new char [FileSize - 11];
  while (this->ReadCommand(File, FileSize - 11, TmpString) == false);
  delete[] TmpString;
  this->Scanline = new PicRGBE [2 * this->Width];
  this->Picture = new RGB* [this->Height];
  int j;
  if (this->TopToBottom == true)
    for (j = 0; j < this->Height; ++j)
      {
	this->Picture[j] = new RGB [this->Width];
	if (this->ReadScanline(File, j) == false)
	  {
	    cout << "error while reading " << FileName << endl;
	    for (int k = j + 1; k < this->Height; ++k)
	      this->Picture[k] = new RGB [this->Width];
	    j = this->Height + 1;
	  }
      }
  else
    for (j = this->Height - 1; j >= 0; --j)
      {
	this->Picture[j] = new RGB [this->Width];
	if (this->ReadScanline(File, j) == false)
	  {
	    cout << "error while reading " << FileName << endl;
	    for (int k = j - 1; k >= 0; --k)
	      this->Picture[k] = new RGB [this->Width];	      
	    j = -1;
	  }
      }
  File.close();
  delete[] this->Scanline;
  if ((j < 0) || (j > this->Height))
    return true;
  return true;
}

// read command in a hdr file 
//
// file = reference on input file
// fileSize = input file size (minus number of already read bytes)
// string = temporary string that has to be used for storing command
// return value = true if it was the last command

bool HdrFormat::ReadCommand (ifstream& file, int fileSize, char* string)
{
  int TmpNbrReadChar = 0;
  if (fileSize == 0)
    return true;
  file.read (string, sizeof(char));
  while ((TmpNbrReadChar <= fileSize) && (string[TmpNbrReadChar] != '\n'))
    {
      ++TmpNbrReadChar;
      file.read (string + TmpNbrReadChar, sizeof(char));
    }
  string[TmpNbrReadChar + 1] = '\0';
  if ((string[0] == '#') || (string[0] == '\n'))
    return false;
  if (((string[0] == '+') || (string[0] == '-')) && ((string[1] == 'X') || (string[1] == 'Y')))
    {
      string[TmpNbrReadChar + 1] = '\0';
      int SignPosition = 2;
      while ((SignPosition < TmpNbrReadChar) && ((string[SignPosition] != '+') && (string[SignPosition] != '-')))
	++SignPosition;
      int Pos = 2;
      while ((Pos < TmpNbrReadChar) && (string[Pos] == ' '))
	++Pos;
      if (string[1] == 'X')
	{
	  this->Width = atoi(string + Pos);
	  if (string[0] == '+')
	    this->LeftToRight = true;
	  else
	    this->LeftToRight = false;
	}
      else
	{
	  this->Height = atoi(string + Pos);
	  if (string[0] == '+')
	    this->TopToBottom = false;
	  else
	    this->TopToBottom = true;
	}
      Pos = SignPosition + 2;
      while ((Pos < TmpNbrReadChar) && (string[Pos] == ' '))
	++Pos;
      if (string[SignPosition + 1] == 'X')
	{
	  this->Width = atoi(string + Pos);
	  if (string[SignPosition] == '+')
	    this->LeftToRight = true;
	  else
	    this->LeftToRight = false;
	}
      else
	{
	  this->Height = atoi(string + Pos);
	  if (string[SignPosition] == '+')
	    this->TopToBottom = false;
	  else
	    this->TopToBottom = true;
	}
      return true;
    }
  int EqualPosition = 0;
  while ((EqualPosition < TmpNbrReadChar) && (string[EqualPosition] != '=')) 
    ++EqualPosition;
  char* TmpString = new char [EqualPosition + 1];
  strncpy (TmpString, string, EqualPosition);
  TmpString[EqualPosition] = '\0';
  this->Commands += TmpString;
  TmpString = new char [TmpNbrReadChar - EqualPosition - 1];
  strncpy (TmpString, string + EqualPosition + 1, TmpNbrReadChar - EqualPosition - 2);
  TmpString[TmpNbrReadChar - EqualPosition - 2] = '\0';
  this->CommandArgument += TmpString;
  return false;
}

// read a scanline in a hdr file 
//
// file = reference on input file
// pos = poistion where to store the scanline
// return value = true if the scanline has beeen read

bool HdrFormat::ReadScanline (ifstream& file, int pos)
{
  unsigned char TmpChar;
  file.read((char*) &TmpChar, sizeof(unsigned char));
  if (TmpChar != 2)
    {
      file.seekg(-1, ios::cur);
      return this->ReadOldScanline(file, pos);
    }
  file.read((char*) &TmpChar, sizeof(unsigned char));
  if (TmpChar != 2)
    {
      file.seekg(-2, ios::cur);
      return this->ReadOldScanline(file, pos);      
    }
  file.read((char*) &TmpChar, sizeof(unsigned char));
  if ((TmpChar & 0x80) != 0)
    {
      file.seekg(-3, ios::cur);
      return this->ReadOldScanline(file, pos);      
    }    
  int TmpWidth = ((unsigned int) TmpChar) << 8;
  file.read((char*) &TmpChar, sizeof(unsigned char));
  TmpWidth |= (unsigned int) TmpChar;
  if (this->Width != TmpWidth)
    return false;

  unsigned char Count;
  for (int i = 0; i < TmpWidth;)
    {
      file.read((char*) &Count, sizeof(unsigned char));
      if (Count > 128)
	{
	  Count &= 127;
	  file.read((char*) &TmpChar, sizeof(unsigned char));	  
	  while (Count--)
	    {
	      this->Scanline[i].Red = TmpChar;
	      ++i;
	    }
	}
      else
	{
	  while (Count--)
	    {
	      file.read((char*) &TmpChar, sizeof(unsigned char));	  
	      this->Scanline[i].Red = TmpChar;
	      ++i;
	    }
	}
    }
  for (int i = 0; i < TmpWidth;)
    {
      file.read((char*) &Count, sizeof(unsigned char));
      if (Count > 128)
	{
	  Count &= 127;
	  file.read((char*) &TmpChar, sizeof(unsigned char));	  
	  while (Count--)
	    {
	      this->Scanline[i].Green = TmpChar;
	      ++i;
	    }
	}
      else
	{
	  while (Count--)
	    {
	      file.read((char*) &TmpChar, sizeof(unsigned char));	  
	      this->Scanline[i].Green = TmpChar;
	      ++i;
	    }
	}
    }
  for (int i = 0; i < TmpWidth;)
    {
      file.read((char*) &Count, sizeof(unsigned char));
      if (Count > 128)
	{
	  Count &= 127;
	  file.read((char*) &TmpChar, sizeof(unsigned char));	  
	  while (Count--)
	    {
	      this->Scanline[i].Blue = TmpChar;
	      ++i;
	    }
	}
      else
	{
	  while (Count--)
	    {
	      file.read((char*) &TmpChar, sizeof(unsigned char));	  
	      this->Scanline[i].Blue = TmpChar;
	      ++i;
	    }
	}
    }
  for (int i = 0; i < TmpWidth;)
    {
      file.read((char*) &Count, sizeof(unsigned char));
      if (Count > 128)
	{
	  Count &= 127;
	  file.read((char*) &TmpChar, sizeof(unsigned char));	  
	  while (Count--)
	    {
	      this->Scanline[i].Exponent = TmpChar;
	      ++i;
	    }
	}
      else
	{
	  while (Count--)
	    {
	      file.read((char*) &TmpChar, sizeof(unsigned char));	  
	      this->Scanline[i].Exponent = TmpChar;
	      ++i;
	    }
	}
    }
  for (int i = 0; i < TmpWidth; ++i)
    this->Scanline[i].GetRGB(this->Picture[pos][i]);
  return true;
}

// read a scanline in a hdr file using old scanline format 
//
// file = reference on input file
// pos = poistion where to store the scanline
// return value = true if the scanline has beeen read

bool HdrFormat::ReadOldScanline (ifstream& file, int pos)
{
  file.read((char*) &(this->Scanline[0].Red), sizeof(unsigned char));
  file.read((char*) &(this->Scanline[0].Green), sizeof(unsigned char));
  file.read((char*) &(this->Scanline[0].Blue), sizeof(unsigned char));
  file.read((char*) &(this->Scanline[0].Exponent), sizeof(unsigned char));
  int Shift = 0;
  int i;
  for (i = 1; i < this->Width;)
    {
      file.read((char*) &(this->Scanline[i].Red), sizeof(unsigned char));
      file.read((char*) &(this->Scanline[i].Green), sizeof(unsigned char));
      file.read((char*) &(this->Scanline[i].Blue), sizeof(unsigned char));
      file.read((char*) &(this->Scanline[i].Exponent), sizeof(unsigned char));
      if ((this->Scanline[i].Red == 1) && (this->Scanline[i].Green == 1) && (this->Scanline[i].Blue == 1))
	{
	  for (int j = (((unsigned int) this->Scanline[i].Exponent) << Shift); j > 0; --j)
	    {
	      this->Scanline[i] = this->Scanline[i - 1];
	      ++i;
	    }
	  Shift += 8;
	}
      else
	{	  
	  ++i;
	  Shift = 0;
	}
    }
  for (i = 0; i < this->Width; ++i)
    this->Scanline[i].GetRGB(this->Picture[pos][i]);
  return true;
}

// clone a bitmap picture type with a new size
// 
// width = new width (0 if old width has to be kept)
// height = new height (0 if old height has to be kept)
// return value = pointer to resulting bitmap picture 

AbstractBitmapPicture* HdrFormat::Clone (int width, int height)
{
  if (width == 0)
    width = this->Width;
  if (height == 0)
    height = this->Height;
  return new HdrFormat (width, height);
}

#ifdef __OPENGL__

// Load bitmap picture as an open gl texture
//
// alpha = alpha factor to use if no alpha channel is defined in bitmap file (if greater than 1.0, use RGB texture only)

void HdrFormat::GLLoadBitmapAsTexture (double alpha)
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

void HdrFormat::GLUseBitmapAsTexture ()
{
  glBindTexture (GL_TEXTURE_2D, this->GLTextureID);
  glScaled(((double) this->Width) / ((double) this->GLWidth), ((double) this->Height) / ((double) this->GLHeight), 0.0);
}

#endif

