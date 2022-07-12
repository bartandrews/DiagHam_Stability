////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of generic 3D tight binding model                 //
//                                                                            //
//                        last modification : 03/10/2012                      //
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


#include "config.h"
#include "Tools/FTITightBinding/Generic3DTightBindingModel.h"
#include "GeneralTools/Endian.h"

#include <fstream>


using std::cout;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::endl;


// default constructor
//
// fileName = name of the binary file that contains the band structure information

Generic3DTightBindingModel::Generic3DTightBindingModel(char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return;
    }
  File.seekg (0l, ios::end);
  unsigned long FileSize = File.tellg ();
  File.close();

  File.open(fileName, ios::binary | ios::in);
  ReadLittleEndian(File, this->NbrBands);
  ReadLittleEndian(File, this->NbrStatePerBand);
  this->Inversion.Resize(this->NbrBands, this->NbrBands);
  for (int i = 0; i < this->NbrBands; ++i)
      for (int j = 0; j < this->NbrBands; ++j)
      {
          Complex TmpInversion = 0.0;
          ReadLittleEndian(File, TmpInversion);
          this->Inversion[i][j] = TmpInversion;
      }
  int HeaderSize = -1;
  ReadLittleEndian(File, HeaderSize);
  int CorrectDimension = 3;
  int CorrectHeaderSize = (((this->NbrBands + 2) * CorrectDimension) * sizeof(double)) + ((CorrectDimension + 1) * sizeof(int));
  if (HeaderSize >= CorrectHeaderSize)
    {
      int TmpDimension = -1;
      ReadLittleEndian(File, TmpDimension);
      HeaderSize -= sizeof(int);
      if (TmpDimension != CorrectDimension)
	{
	  cout << "error : " << fileName << " does not contain a valid 3D band structure" << endl;
	  this->NbrBands = 0;
	  this->NbrStatePerBand = 0;
	  File.close();
	  return;
	}
      ReadLittleEndian(File, this->NbrSiteX);
      ReadLittleEndian(File, this->KxFactor);
      ReadLittleEndian(File, this->GammaX);	 
      this->EmbeddingX.Resize(this->NbrBands);
      for (int i = 0; i < this->NbrBands; ++i)
      {
          double Tmp = 0.0;
          ReadLittleEndian(File, Tmp);
          this->EmbeddingX[i] = Tmp;
      }
      ReadLittleEndian(File, this->NbrSiteY);
      ReadLittleEndian(File, this->KyFactor);
      ReadLittleEndian(File, this->GammaY);	  
      this->EmbeddingY.Resize(this->NbrBands);
      for (int i = 0; i < this->NbrBands; ++i)
      {
          double Tmp = 0.0;
          ReadLittleEndian(File, Tmp);
          this->EmbeddingY[i] = Tmp;
      }
      ReadLittleEndian(File, this->NbrSiteZ);
      ReadLittleEndian(File, this->KzFactor);
      ReadLittleEndian(File, this->GammaZ);	  
      this->EmbeddingZ.Resize(this->NbrBands);
      for (int i = 0; i < this->NbrBands; ++i)
      {
          double Tmp = 0.0;
          ReadLittleEndian(File, Tmp);
          this->EmbeddingZ[i] = Tmp;
      }
      HeaderSize -= (CorrectHeaderSize - sizeof(int));
      File.seekg (HeaderSize, ios::cur);
    }
  else
    {
      cout << "error : " << fileName << " does not contain a valid 3D band structure" << endl;
      this->NbrBands = 0;
      this->NbrStatePerBand = 0;
      File.close();
      return;
    }
  this->EnergyBandStructure = new double*[this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
      for (int j = 0; j < this->NbrStatePerBand; ++j)
	{
	  ReadLittleEndian(File, this->EnergyBandStructure[i][j]);
	}
    }
  if (FileSize == ((sizeof(double) * this->NbrStatePerBand * this->NbrBands) + sizeof(long) + sizeof(int) + sizeof(int) + (this->NbrBands * this->NbrBands * sizeof(Complex)) + HeaderSize))
    {
      this->OneBodyBasis = 0;
    }
  else
    {
      this->OneBodyBasis = new ComplexMatrix [this->NbrStatePerBand];
      for (int j = 0; j < this->NbrStatePerBand; ++j)	
	{
	  this->OneBodyBasis[j].ReadMatrix(File);
	}     
    }
  File.close();
}

// destructor
//

Generic3DTightBindingModel::~Generic3DTightBindingModel()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute
void Generic3DTightBindingModel::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
}
