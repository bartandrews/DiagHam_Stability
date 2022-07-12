////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of lz-sorted spectrum for QHE on sphere              //
//                                                                            //
//                        last modification : 18/04/2005                      //
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
#include "Tools/FQHESpectrum/QHEOnSphereLzSortedSpectrum.h"

#include <string.h>
#include <fstream>
#include <math.h>


using std::ifstream;
using std::ios;
using std::endl;
using std::cout;


// constructor
// 
// nbrParticles = number of particles
// lzMax = maximum Lz value reached by a particle
// fermionicFlag = fermionic statistics flag (true if the system if fermionic)
// fileName = name of the file that contains the spectrum datas
// error = relative error that has to be used to test if two energies are degenerated 

QHEOnSphereLzSortedSpectrum::QHEOnSphereLzSortedSpectrum (int nbrParticles, int lzMax, bool fermionicFlag, char* fileName, double error)
{
  this->NbrParticles = nbrParticles;
  this->LzMax = lzMax;
  this->FermionicFlag = fermionicFlag;
  this->Error = error;
  if (this->FermionicFlag == true)
    this->MaxTotalLz = (this->NbrParticles * (this->LzMax - this->NbrParticles + 1)) >> 1;	 
  else
    this->MaxTotalLz = (this->NbrParticles * this->LzMax) >> 1;
  this->NbrEnergies = new int [this->MaxTotalLz + 1];
  this->NbrDistinctEnergies = new int [this->MaxTotalLz + 1];
  this->Spectrum = new double* [this->MaxTotalLz + 1];
  this->Degeneracy = new int* [this->MaxTotalLz + 1];
  this->ConvertionTable = new int* [this->MaxTotalLz + 1];
  for (int i = 0; i <= this->MaxTotalLz; ++i)
    {
      this->NbrEnergies[i] = 0;
      this->NbrDistinctEnergies[i] = 0;
    }
  if (this->ParseSpectrumFile(fileName) == false)
    {
      this->NbrParticles = 0;
      delete[] this->Spectrum;
      delete[] this->Degeneracy;
      delete[] this->ConvertionTable;
      delete[] this->NbrEnergies;	  
      delete[] this->NbrDistinctEnergies;
    }
}

// constructor from a file, retrieving other informations from its name
// 
// fileName = name of the file that contains the spectrum datas
// error = relative error that has to be used to test if two energies are degenerated 

QHEOnSphereLzSortedSpectrum::QHEOnSphereLzSortedSpectrum (char* fileName, double error)
{
  if (this->RetrieveInformationFromName(fileName, this->NbrParticles, this->LzMax, this->FermionicFlag) == true)
    {
      this->Error = error;
      if (this->FermionicFlag == true)
	this->MaxTotalLz = (this->NbrParticles * (this->LzMax - this->NbrParticles + 1)) >> 1;	 
      else
	this->MaxTotalLz = (this->NbrParticles * this->LzMax) >> 1;
      this->NbrEnergies = new int [this->MaxTotalLz + 1];
      this->NbrDistinctEnergies = new int [this->MaxTotalLz + 1];
      this->Spectrum = new double* [this->MaxTotalLz + 1];
      this->Degeneracy = new int* [this->MaxTotalLz + 1];
      this->ConvertionTable = new int* [this->MaxTotalLz + 1];
      for (int i = 0; i <= this->MaxTotalLz; ++i)
	{
	  this->NbrEnergies[i] = 0;
	  this->NbrDistinctEnergies[i] = 0;
	}
      if (this->ParseSpectrumFile(fileName) == false)
	{
	  this->NbrParticles = 0;
	  delete[] this->Spectrum;
	  delete[] this->Degeneracy;
	  delete[] this->ConvertionTable;
	  delete[] this->NbrEnergies;	
	  delete[] this->NbrDistinctEnergies;
	}
    }
  else
    {
      this->NbrParticles = 0;
    }
}
  
// destructor
//

QHEOnSphereLzSortedSpectrum::~QHEOnSphereLzSortedSpectrum ()
{
  if (this->NbrParticles != 0)
    {
      for (int i = 0; i <= this->MaxTotalLz; ++i)
	{
	  if (this->NbrEnergies[i] > 0)
	    {
	      delete[] this->Spectrum[i];
	      delete[] this->Degeneracy[i];
	      delete[] this->ConvertionTable[i];
	    }
	}
      delete[] this->Spectrum;
      delete[] this->Degeneracy;
      delete[] this->ConvertionTable;
      delete[] this->NbrEnergies;
      delete[] this->NbrDistinctEnergies;
    }
}

// test if read spectrum is valid
//
// retur value = true if the spectrum has been read and is valid

bool QHEOnSphereLzSortedSpectrum::IsSpectrumValid()
{
  if (this->NbrParticles != 0)
    return true;
  else
    return false;
}

// parse spectrum content from a file 
//
// fileName = name of the file that contains spectrum datas
// return value = true if no error occurs

bool QHEOnSphereLzSortedSpectrum::ParseSpectrumFile(char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      return false;
    } 

  File.seekg(0, ios::end);
  long FileSize = File.tellg();
  File.seekg(0, ios::beg);
  int TmpLzValue;
  double Dummy;
  while ((File.tellg() < FileSize) && (File.tellg() >= 0))
    {
      File >> TmpLzValue >> Dummy;
      if ((TmpLzValue <= this->MaxTotalLz) && (TmpLzValue >= 0))
	this->NbrEnergies[TmpLzValue]++;
      TmpLzValue = -1;
    } 
  File.close();

  for (int i = 0; i <= this->MaxTotalLz; ++i)
    {
      if (this->NbrEnergies[i] > 0)
	{
	  this->Spectrum[i] = new double [this->NbrEnergies[i]];
	  this->ConvertionTable[i] = new int [this->NbrEnergies[i]];
	  this->Degeneracy[i] = new int [this->NbrEnergies[i]];
	  this->NbrEnergies[i] = 0;
	}
    }
  
  ifstream File2;
  File2.open(fileName, ios::binary | ios::in);
  while ((File2.tellg() < FileSize) && (File2.tellg() >= 0))
    {
      File2 >> TmpLzValue >> Dummy;
      if ((TmpLzValue <= this->MaxTotalLz) && (TmpLzValue >= 0))
	{
	  if (this->NbrEnergies[TmpLzValue] == 0)
	    {
	      this->Spectrum[TmpLzValue][0] = Dummy;
	      this->Degeneracy[TmpLzValue][0] = 1;
	      this->ConvertionTable[TmpLzValue][0] = 0;
	      this->NbrEnergies[TmpLzValue]++;
	    }
	  else
	    {
	      double Diff = fabs(this->Spectrum[TmpLzValue][this->NbrDistinctEnergies[TmpLzValue]] - Dummy);
	      if ((Diff < (this->Error * fabs(Dummy))) || (Diff < this->Error))
		{
		  this->Degeneracy[TmpLzValue][this->NbrDistinctEnergies[TmpLzValue]]++;
		  this->ConvertionTable[TmpLzValue][this->NbrEnergies[TmpLzValue]] = this->NbrDistinctEnergies[TmpLzValue];
		}
	      else
		{
		  this->NbrDistinctEnergies[TmpLzValue]++;
		  this->Spectrum[TmpLzValue][this->NbrDistinctEnergies[TmpLzValue]] = Dummy;
		  this->Degeneracy[TmpLzValue][this->NbrDistinctEnergies[TmpLzValue]] = 1;
		  this->ConvertionTable[TmpLzValue][this->NbrEnergies[TmpLzValue]] = this->NbrDistinctEnergies[TmpLzValue];
		}
	      this->NbrEnergies[TmpLzValue]++;
	    }
	}
      TmpLzValue = -1;
    } 
 
  File2.close();

  for (int i = 0; i <= this->MaxTotalLz; ++i)
    if (this->NbrEnergies[i] > 0)
      this->NbrDistinctEnergies[i]++;

  return true;
}

// get system information from a formatted spectrum file name 
//
// fileName = name of the file that contains spectrum datas (can include partial of full path)
// nbrParticles = reference on the variable where the number of particles has to be stored
// lzMax = reference on the variable where the maximum Lz value reached by a particle has to be stored
// fermionicFlag = reference on the variable where the fermionic statistics flag (true if the system if fermionic) has to be stored
// return value = true if no error occurs (aka all values have been successfully retrieved)

bool QHEOnSphereLzSortedSpectrum::RetrieveInformationFromName(char* fileName, int& nbrParticles, int& lzMax, bool& fermionicFlag)
{
  int End = strlen(fileName);
  int Start = End;
  while ((Start > 0) && (fileName[Start] != '/'))
    --Start;
  fileName += Start;
  char* TmpString = strstr (fileName, "fermions");
  if (TmpString == 0)
    fermionicFlag = false;
  else
    fermionicFlag = true;
  TmpString = strstr (fileName, "_n_");
  if (TmpString == 0)
    return false;
  TmpString += 3;
  char* TmpError;
  nbrParticles = strtol(TmpString, &TmpError, 0);
  if (TmpError == TmpString)
    {
      return false;
    }
  TmpString = strstr (fileName, "_2s_");
  if (TmpString == 0)
    return false;
  TmpString += 4;
  lzMax = strtol(TmpString, &TmpError, 0);
  if (TmpError == TmpString)
    {
      return false;
    }
  return true;
}

// print spectrum
//
// str = reference on the output stream
// showDegeneracy = true if degeneracy has to be written
// return value = reference on the output stream

ostream& QHEOnSphereLzSortedSpectrum::PrintSpectrum (ostream& str, bool showDegeneracy)
{
  if (this->NbrParticles != 0)
    {
      if (showDegeneracy == false)
	{
	  for (int i = 0; i <= this->MaxTotalLz; ++i)
	    if (this->NbrEnergies[i] > 0) 
	      for (int j = 0; j < this->NbrEnergies[i]; ++j)
		str << i << " " << this->Spectrum[i][this->ConvertionTable[i][j]] << endl;
	}
      else
	{
	  for (int i = 0; i <= this->MaxTotalLz; ++i)
	    if (this->NbrEnergies[i] > 0) 
	      for (int j = 0; j < this->NbrDistinctEnergies[i]; ++j)
		str << i << " " << this->Spectrum[i][j] << " "  << this->Degeneracy[i][j] << endl;
	}
    }
  return str;
}
