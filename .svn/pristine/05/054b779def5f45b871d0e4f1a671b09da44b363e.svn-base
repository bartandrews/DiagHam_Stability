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
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>


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

FQHEOnSphereLzSortedSpectrum::FQHEOnSphereLzSortedSpectrum (int nbrParticles, int lzMax, bool fermionicFlag, char* fileName, double error)
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

FQHEOnSphereLzSortedSpectrum::FQHEOnSphereLzSortedSpectrum (char* fileName, double error)
{
  this->NbrParticles = 0;
  this->LzMax = 0;
  this->FermionicFlag = false;
  if (FQHEOnSphereFindSystemInfoFromFileName(fileName, this->NbrParticles, this->LzMax, this->FermionicFlag) == true)
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

FQHEOnSphereLzSortedSpectrum::~FQHEOnSphereLzSortedSpectrum ()
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

bool FQHEOnSphereLzSortedSpectrum::IsSpectrumValid()
{
  if (this->NbrParticles != 0)
    return true;
  else
    return false;
}

// get the highest Lz value avalailable within the spectrum
// 
// return value = twice the highest Lz value

int FQHEOnSphereLzSortedSpectrum::GetMaxLzValue ()
{
  int TmpLz = this->MaxTotalLz;
  while ((TmpLz > 0) && (this->NbrEnergies[TmpLz] == 0))
    {
      --TmpLz;  
    }
  TmpLz <<= 1;
  if (((this->NbrParticles * this->LzMax) & 1) != 0)
    ++TmpLz;
  return TmpLz;
}

// parse spectrum content from a file 
//
// fileName = name of the file that contains spectrum datas
// return value = true if no error occurs

bool FQHEOnSphereLzSortedSpectrum::ParseSpectrumFile(char* fileName)
{
  MultiColumnASCIIFile Spectrum;
  if (Spectrum.Parse(fileName) == false)
    {
      Spectrum.DumpErrors(cout);
      return false;
    }

  long TotalNbrEnergies = Spectrum.GetNbrLines();
  int* TmpLzValues = Spectrum.GetAsIntegerArray(0);
  double* TmpEnergies = Spectrum.GetAsDoubleArray(1);

  int TmpLzValue = 0;
  double TmpEnergy = 0.0;

  for (long i = 0l; i < TotalNbrEnergies; ++i)
    {
      TmpLzValue = TmpLzValues[i];
      if ((TmpLzValue <= this->MaxTotalLz) && (TmpLzValue >= 0))
	this->NbrEnergies[TmpLzValue]++;
    } 

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
  
  for (long i = 0l; i < TotalNbrEnergies; ++i)
    {
      TmpLzValue = TmpLzValues[i];
      TmpEnergy = TmpEnergies[i];
      if ((TmpLzValue <= this->MaxTotalLz) && (TmpLzValue >= 0))
	{
	  if (this->NbrEnergies[TmpLzValue] == 0)
	    {
	      this->Spectrum[TmpLzValue][0] = TmpEnergy;
	      this->Degeneracy[TmpLzValue][0] = 1;
	      this->ConvertionTable[TmpLzValue][0] = 0;
	      this->NbrEnergies[TmpLzValue]++;
	    }
	  else
	    {
	      double Diff = fabs(this->Spectrum[TmpLzValue][this->NbrDistinctEnergies[TmpLzValue]] - TmpEnergy);
	      if ((Diff < (this->Error * fabs(TmpEnergy))) || (Diff < this->Error))
		{
		  this->Degeneracy[TmpLzValue][this->NbrDistinctEnergies[TmpLzValue]]++;
		  this->ConvertionTable[TmpLzValue][this->NbrEnergies[TmpLzValue]] = this->NbrDistinctEnergies[TmpLzValue];
		}
	      else
		{
		  this->NbrDistinctEnergies[TmpLzValue]++;
		  this->Spectrum[TmpLzValue][this->NbrDistinctEnergies[TmpLzValue]] = TmpEnergy;
		  this->Degeneracy[TmpLzValue][this->NbrDistinctEnergies[TmpLzValue]] = 1;
		  this->ConvertionTable[TmpLzValue][this->NbrEnergies[TmpLzValue]] = this->NbrDistinctEnergies[TmpLzValue];
		}
	      this->NbrEnergies[TmpLzValue]++;
	    }
	}
      TmpLzValue = -1;
    } 
 
  for (int i = 0; i <= this->MaxTotalLz; ++i)
    if (this->NbrEnergies[i] > 0)
      this->NbrDistinctEnergies[i]++;

  delete[] TmpLzValues;
  delete[] TmpEnergies;

  return true;
}

// print spectrum
//
// str = reference on the output stream
// showDegeneracy = true if degeneracy has to be written
// return value = reference on the output stream

ostream& FQHEOnSphereLzSortedSpectrum::PrintSpectrum (ostream& str, bool showDegeneracy)
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
