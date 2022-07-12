////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.06                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
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


#include "BitmapTools/BitmapPicture/AbstractBitmapPicture.h"
#include "BitmapTools/BitmapPicture/JpgFormat.h"
#include "BitmapTools/BitmapPicture/TgaFormat.h"
#include "BitmapTools/BitmapPicture/BmpFormat.h"
#include "BitmapTools/BitmapPicture/PcxFormat.h"
#include "BitmapTools/BitmapPicture/TiffFormat.h"
#include "BitmapTools/BitmapPicture/GifFormat.h"
#include "BitmapTools/BitmapPicture/HdrFormat.h"
#include "BitmapTools/BitmapPicture/IffFormat.h"
#include "BitmapTools/Color/Color.h"
#include "BitmapTools/Color/RGB.h"

#include <string.h>
#include <math.h>

// virtual destructor
//

AbstractBitmapPicture::~AbstractBitmapPicture() {};

// Set Pixel Value
//
// x = pixel coordinate along x axis
// y = pixel coordinate along y axis
// color = pixel color

void AbstractBitmapPicture::SetPixel(int x, int y, Color& color)
{
  PicRGB TmpRGB2 = color.GetRGBComponent();
  this->SetPixel(x, y, TmpRGB2);
}
  
// Get Pixel Value
//
// x = pixel coordinate along x axis
// y = pixel coordinate along y axis
// color = reference on color that will contain pixel color

void AbstractBitmapPicture::GetPixel(int x, int y, Color& color)
{
  color = this->GetPixel(x, y);
}

// read a picture from a file and return corresponding AbstractBitmapPicture object
//
// fileName = file name 
// return value = pointer on corresponding AbstractBitmapPicture

AbstractBitmapPicture* ReadBitmapPicture (char* fileName)
{
  AbstractBitmapPicture* Picture = 0;
  char* Extension = strrchr (fileName, (int) '.');
  if (strcasecmp (Extension, ".tga") == 0)
    Picture = (AbstractBitmapPicture*) new TgaFormat();
  else
    if (strcasecmp (Extension, ".bmp") == 0)
      Picture = (AbstractBitmapPicture*) new BmpFormat();
    else
      if (strcasecmp (Extension, ".jpg") == 0)
	Picture = (AbstractBitmapPicture*) new JpgFormat();
      else
	if (strcasecmp (Extension, ".tif") == 0)
	  Picture = (AbstractBitmapPicture*) new TiffFormat();
	else
	  if (strcasecmp (Extension, ".pcx") == 0)
	    Picture = (AbstractBitmapPicture*) new PcxFormat();
	  else
	    if (strcasecmp (Extension, ".gif") == 0)
	      Picture = (AbstractBitmapPicture*) new GifFormat();
	    else
	      if (strcasecmp (Extension, ".gif") == 0)
		Picture = (AbstractBitmapPicture*) new GifFormat();
	      else
		if (strcasecmp (Extension, ".hdr") == 0)
		  Picture = (AbstractBitmapPicture*) new HdrFormat();
		else
		  if (strcasecmp (Extension, ".iff") == 0)
		    Picture = (AbstractBitmapPicture*) new IffFormat();
  if ((Picture != 0) && (Picture->LoadPicture (fileName) == false))
    {
      delete Picture;
      Picture = 0;
    }
  return Picture;
}

// save a picture in a given format
//
// fileName = file name 
// inputPicture = picture to save
// quality = picture quality for output format (if not in the range 0 to 100, take input picture quality)
// return value = true if picture have been correctly saved

bool SaveBitmapPicture (char* fileName, AbstractBitmapPicture* inputPicture, int quality)
{
  AbstractBitmapPicture* Picture = 0;
  char* Extension = strrchr (fileName, (int) '.');
  if (strcasecmp (Extension, ".tga") == 0)
    Picture = (AbstractBitmapPicture*) new TgaFormat(inputPicture->PictureWidth(), inputPicture->PictureHeight());
  else
    if (strcasecmp (Extension, ".bmp") == 0)
      Picture = (AbstractBitmapPicture*) new BmpFormat(inputPicture->PictureWidth(), inputPicture->PictureHeight());
    else
      if (strcasecmp (Extension, ".jpg") == 0)
	Picture = (AbstractBitmapPicture*) new JpgFormat(inputPicture->PictureWidth(), inputPicture->PictureHeight());
      else
	if (strcasecmp (Extension, ".tif") == 0)
	  Picture = (AbstractBitmapPicture*) new TiffFormat(inputPicture->PictureWidth(), inputPicture->PictureHeight());
	else
	  if (strcasecmp (Extension, ".pcx") == 0)
	    Picture = (AbstractBitmapPicture*) new PcxFormat(inputPicture->PictureWidth(), inputPicture->PictureHeight());
	  else
	    if (strcasecmp (Extension, ".gif") == 0)
	      Picture = (AbstractBitmapPicture*) new GifFormat(inputPicture->PictureWidth(), inputPicture->PictureHeight());
	    else
	      if (strcasecmp (Extension, ".hdr") == 0)
		Picture = (AbstractBitmapPicture*) new HdrFormat(inputPicture->PictureWidth(), inputPicture->PictureHeight());
	      else
		if (strcasecmp (Extension, ".iff") == 0)
		  Picture = (AbstractBitmapPicture*) new IffFormat(inputPicture->PictureWidth(), inputPicture->PictureHeight());
  bool Flag = false;
  if (Picture != 0)
    {  
      if ((quality < 0) || (quality > 100))
	Picture->SetQuality(inputPicture->GetQuality());
      else
	Picture->SetQuality(quality);
      PicRGB TmpPicRGB;
      for (int i = 0; i < inputPicture->PictureWidth(); i++)
	for (int j = 0; j < inputPicture->PictureHeight(); j++)
	  {
	    TmpPicRGB = inputPicture->GetPixel(i, j);
	    Picture->SetPixel(i, j, TmpPicRGB);
	  }
      Flag = Picture->SavePicture(fileName);
      delete Picture;
    }
  return Flag;
}

// create a picture in a given format
//
// fileName = file name with an extension corresponding to the chosen file format
// width = picture width
// height = picture height
// return value = pointer to the created bitmap picture

AbstractBitmapPicture* CreateBitmapPicture (char* fileName, int width, int height)
{
  AbstractBitmapPicture* Picture = 0;
  char* Extension = strrchr (fileName, (int) '.');
  if ((strcmp (Extension, ".tga") == 0) || (strcmp (Extension, ".TGA") == 0))
    Picture = (AbstractBitmapPicture*) new TgaFormat(width, height);
  else
    if ((strcmp (Extension, ".bmp") == 0) || (strcmp (Extension, ".BMP") == 0))
      Picture = (AbstractBitmapPicture*) new BmpFormat(width, height);
    else
      if ((strcmp (Extension, ".jpg") == 0) || (strcmp (Extension, ".JPG") == 0))
	Picture = (AbstractBitmapPicture*) new JpgFormat(width, height);
      else
	if ((strcmp (Extension, ".tif") == 0) || (strcmp (Extension, ".TIF") == 0))
	  Picture = (AbstractBitmapPicture*) new TiffFormat(width, height);
	else
	  if ((strcmp (Extension, ".pcx") == 0) || (strcmp (Extension, ".PCX") == 0))
	    Picture = (AbstractBitmapPicture*) new PcxFormat(width, height);
	  else
	    if ((strcmp (Extension, ".gif") == 0) || (strcmp (Extension, ".GIF") == 0))
	      Picture = (AbstractBitmapPicture*) new GifFormat(width, height);
	    else
	      if ((strcmp (Extension, ".hdr") == 0) || (strcmp (Extension, ".HDR") == 0))
		Picture = (AbstractBitmapPicture*) new HdrFormat(width, height);
  return Picture;
}

// get picture quality for destructive compression algorithm
//
// return value = picture quality (range from 0 to 100)

int AbstractBitmapPicture::GetQuality ()
{
  return 100;
}
 
// resample a picture and store result in a new one
//
// width = new width
// height = new height (0 if ratio has to be kept)

AbstractBitmapPicture* AbstractBitmapPicture::Resample (int width, int height)
{
  if (height == 0)
    height = (width * this->PictureHeight()) /this->PictureWidth();
  AbstractBitmapPicture* TmpPicture = this->Clone(width, height);
  if ((width <= this->PictureWidth()) && (height <= this->PictureHeight()))
    {
      double XInc = (double) this->PictureWidth() / (double) width;
      double YInc = (double) this->PictureHeight() / (double) height;
      double YWeight1 = 1.0;
      double YWeight2 = 1.0;
      double XWeight1 = 1.0;
      double XWeight2 = 1.0;
      int XStart = 0;
      int YStart = 0;
      int XEnd = 0;
      int YEnd = 0;
      PicRGB TmpRGB;
      for (int i = 0; i < height; i++)
	{
	  YWeight1 = YInc * (double) i;
	  YWeight2 = YWeight1 + YInc;
	  YStart = (int) YWeight1;
	  YWeight1 = 1.0 - (YWeight1 - (double) YStart);
	  YEnd = (int) YWeight2;
	  YWeight2 = (YWeight2 - (double) YEnd);
	  if (YWeight2 < 1e-12)
	    YWeight2 = 0;
	  for (int j = 0; j < width; j++)
	    {
	      XWeight1 = XInc * (double) j;
	      XWeight2 = XWeight1 + XInc;
	      XStart = (int) XWeight1;
	      XWeight1 = 1.0 - (XWeight1 - (double) XStart);
	      XEnd = (int) XWeight2;
	      XWeight2 = (XWeight2 - (double) XEnd);
	      Color TmpCol = ((RGB) this->GetPixel(XStart, YStart)) * XWeight1;
	      for (int k = XStart + 1; k < XEnd; k++)
		TmpCol += ((RGB) this->GetPixel(k, YStart));
	      if (XWeight2 > 0)
		TmpCol += ((RGB) this->GetPixel(XEnd, YStart)) * XWeight2;
	      TmpCol = TmpCol.GetRGBComponent() * YWeight1;
	      for (int k = YStart + 1; k < YEnd; k++)
		{
		  TmpCol += ((RGB) this->GetPixel(XStart, k)) * XWeight1;
		  for (int l = XStart + 1; l < XEnd; l++)
		    TmpCol += ((RGB) this->GetPixel(l, k));
		  if (XWeight2 > 0)
		    TmpCol += ((RGB) this->GetPixel(XEnd, k)) * XWeight2;
		}
	      if (YWeight2 > 0)
		{
		  Color TmpCol2 = ((RGB) this->GetPixel(XStart, YEnd)) * XWeight1;
		  for (int k = XStart + 1; k < XEnd; k++)
		    TmpCol2 += ((RGB) this->GetPixel(k, YEnd));
		  if (XWeight2 > 0)
		    TmpCol2 += ((RGB) this->GetPixel(XEnd, YEnd)) * XWeight2;
		  TmpCol2 = TmpCol2.GetRGBComponent() * YWeight2;
		  TmpCol += TmpCol2;
		}
	      TmpCol /= (XInc * YInc);
	      TmpRGB = TmpCol.GetRGBComponent();
	      TmpPicture->SetPixel(j, i, TmpRGB);
	    }
	}
    }
  return TmpPicture;
}

// rotate a picture and store result in a new one (counter clockwise orientation)
//
// angle = rotation angle in degree
// return value = pointer on the new resulting picture

AbstractBitmapPicture* AbstractBitmapPicture::Rotate (double angle)
{
  angle = fmod (angle, 360.0);
  if (angle == 0)
    return this;
  if (angle == 90.0)
    {
      PicRGB TmpRGB; 
      AbstractBitmapPicture* TmpPicture = this->Clone(this->PictureHeight(), this->PictureWidth());
      for (int i = 0; i < this->PictureWidth(); i++)
	for (int j = 0; j < this->PictureHeight(); j++)
	  {
	    TmpRGB = this->GetPixel(i, j);
	    TmpPicture->SetPixel(j, this->PictureWidth() - 1 - i, TmpRGB);
	  }
      return TmpPicture;
    }
  if (angle == 270.0)
    {
      PicRGB TmpRGB; 
      AbstractBitmapPicture* TmpPicture = this->Clone(this->PictureHeight(), this->PictureWidth());
      for (int i = 0; i < this->PictureWidth(); i++)
	for (int j = 0; j < this->PictureHeight(); j++)
	  {
	    TmpRGB = this->GetPixel(i, j);
	    TmpPicture->SetPixel(this->PictureHeight() - 1 - j, i,  TmpRGB);
	  }
      return TmpPicture;
    }
  return 0;
}
