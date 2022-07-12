////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of generic 2D tight binding model                 //
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
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
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

Generic2DTightBindingModel::Generic2DTightBindingModel(char* fileName)
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

  // data written before header5
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

  // processing generic read instructions
  int HeaderSize = this->ReadHeader(File);
  if (HeaderSize<0)
    {
      cout << "error : " << fileName << " does not contain a valid 2D band structure" << endl;
      this->NbrBands = 0;
      this->NbrStatePerBand = 0;
      File.close();
      return;
    }
  this->ReadEigensystem(File, HeaderSize, FileSize);
  File.close();

  this->Nx1 = this->NbrSiteX;
  this->Ny1 = 0;
  this->Nx2 = 0;
  this->Ny2 = this->NbrSiteY;  
  this->ProjectedMomenta = new double* [this->NbrStatePerBand];
  for (int i = 0; i < this->NbrStatePerBand; ++i)
    this->ProjectedMomenta[i] = new double [2];
  this->ComputeAllProjectedMomenta();
}

// destructor
//

Generic2DTightBindingModel::~Generic2DTightBindingModel()
{
}


// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute
void Generic2DTightBindingModel::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
}
