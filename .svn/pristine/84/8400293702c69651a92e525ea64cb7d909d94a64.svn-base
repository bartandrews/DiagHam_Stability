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


#ifndef FQHEONSPHERELZSORTEDSPECTRUM_H
#define FQHEONSPHERELZSORTEDSPECTRUM_H


#include "config.h"

#include <iostream>


using std::ostream;


class FQHEOnSphereLzSortedSpectrum
{

 protected:

  // number of particles
  int NbrParticles;
  // twice the maximum Lz value reached by a particle
  int LzMax;
  // maximum total Lz value that can be reached by the system
  int MaxTotalLz;
  // fermionic statistics flag 
  bool FermionicFlag;

  // relative error that has to be used to test if two energies are degenerated 
  double Error;

  // number of energies per total Lz value
  int* NbrEnergies;
  // number of distinct energies per total Lz value
  int* NbrDistinctEnergies;
  // array containing the spectrum (without duplicating degenerate values)
  double** Spectrum;
  // degeneracy of each energy value
  int** Degeneracy;
  // convertion table that gives index in the Spectrum table from the absolute index (i.e. the one that doesn't take degeneracy into account)
  int** ConvertionTable;
  
 public:

  // constructor
  // 
  // nbrParticles = number of particles
  // lzMax = twice the maximum Lz value reached by a particle
  // fermionicFlag = fermionic statistics flag (true if the system if fermionic)
  // fileName = name of the file that contains the spectrum datas
  // error = relative error that has to be used to test if two energies are degenerated 
  FQHEOnSphereLzSortedSpectrum (int nbrParticles, int lzMax, bool fermionicFlag, char* fileName, double error = 1e-13);

  // constructor from a file, retrieving other informations from its name
  // 
  // fileName = name of the file that contains the spectrum datas
  // error = relative error that has to be used to test if two energies are degenerated 
  FQHEOnSphereLzSortedSpectrum (char* fileName, double error = 1e-13);
  
  // destructor
  //
  ~FQHEOnSphereLzSortedSpectrum ();

  // test if read spectrum is valid
  //
  // retur value = true if the spectrum has been read and is valid
  bool IsSpectrumValid();

  // get degeneracy of a given energy
  //
  // lz = twice the value of the total Lz
  // index = index of the state corresponding energy  (absolute index i.e. the one that doesn't take degeneracy into account)
  // return value = degeneracy (-1 if an error occured)
  int GetDegeneracy (int lz, int index);

  // get energy value
  //
  // lz = twice the value of the total Lz
  // index = index of the state corresponding energy  (absolute index i.e. the one that doesn't take degeneracy into account)
  // return value = energy value
  double GetEnergy (int lz, int index);

  // get the highest Lz value avalailable within the spectrum
  // 
  // return value = twice the highest Lz value
  int GetMaxLzValue ();

  // print spectrum
  //
  // str = reference on the output stream
  // showDegeneracy = true if degeneracy has to be written
  // return value = reference on the output stream  
  ostream& PrintSpectrum (ostream& str, bool showDegeneracy = false);

 protected:

  // parse spectrum content from a file 
  //
  // fileName = name of the file that contains spectrum datas
  // return value = true if no error occurs
  bool ParseSpectrumFile(char* fileName);

};

// get degeneracy of a given energy
//
// lz = twice the value of the total Lz
// index = index of the state corresponding energy  (absolute index i.e. the one that doesn't take degeneracy into account)
// return value = degeneracy (-1 if an error occured)

inline int FQHEOnSphereLzSortedSpectrum::GetDegeneracy (int lz, int index)
{
  if (((lz >> 1) <= MaxTotalLz) && (index < this->NbrEnergies[lz >> 1]))
    return this->Degeneracy[lz >> 1][this->ConvertionTable[lz >> 1][index]];
  else
    return -1;
}

// get energy value
//
// lz = twice the value of the total Lz
// index = index of the state corresponding energy  (absolute index i.e. the one that doesn't take degeneracy into account)
// return value = energy value

inline double FQHEOnSphereLzSortedSpectrum::GetEnergy (int lz, int index)
{
  return this->Spectrum[lz >> 1][this->ConvertionTable[lz >> 1][index]];
}

#endif
