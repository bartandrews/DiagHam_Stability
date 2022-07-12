////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of Landau level spectrum on sphere                 //
//                                                                            //
//                        last modification : 14/01/2005                      //
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
#include "Tools/FQHEWaveFunction/LandauSpectrumOnSphere.h"

#include <string.h>
#include <iostream>


using std::cout;
using std::endl;


// default constructor
//

LandauSpectrumOnSphere::LandauSpectrumOnSphere()
{
  this->NbrLandauLevels = 0;
  this->TwiceS = 0;
  this->LevelOccupation = 0;
  this->NbrParticles = 0;
}

// constructor from datas
//
// nbrLandauLevels = number of Landau levels
// momentum = twice the value of the momentum in the lowest Landau level
// description = pointer to the string containing the description

LandauSpectrumOnSphere::LandauSpectrumOnSphere(int nbrLandauLevels, int momentum, char* description)
{
  this->NbrLandauLevels = nbrLandauLevels;
  this->TwiceS = momentum;
  this->NbrParticles = this->ParseOccupationDescription(description, this->LevelOccupation);
  if (this->NbrParticles != 0)
    {
      this->Flag.Initialize();
    }
  else
    {
      this->NbrLandauLevels = 0;
    }
}

// copy constructor
//
// spectrum = reference on the spectrum to copy

LandauSpectrumOnSphere::LandauSpectrumOnSphere(const LandauSpectrumOnSphere& spectrum)
{
  this->NbrLandauLevels = spectrum.NbrLandauLevels;
  this->TwiceS = spectrum.TwiceS;
  this->LevelOccupation = spectrum.LevelOccupation;
  this->NbrParticles = spectrum.NbrParticles;
  this->Flag = spectrum.Flag;
}

// destructor
//

LandauSpectrumOnSphere::~LandauSpectrumOnSphere()
{
  if ((this->NbrLandauLevels != 0) && (this->Flag.Shared() == false))
    {
      for (int i = 0; i < this->NbrLandauLevels; ++i)
	delete[] this->LevelOccupation[i];
      delete[] this->LevelOccupation;
    }
}

// assignement
//
// spectrum = reference on the spectrum to assign
// return value = reference on current spectrum

LandauSpectrumOnSphere& LandauSpectrumOnSphere::operator = (const LandauSpectrumOnSphere& spectrum)
{
  if ((this->NbrLandauLevels != 0) && (this->Flag.Shared() == false))
    {
      for (int i = 0; i < this->NbrLandauLevels; ++i)
	delete[] this->LevelOccupation[i];
      delete[] this->LevelOccupation;
    }
  this->NbrLandauLevels = spectrum.NbrLandauLevels;
  this->TwiceS = spectrum.TwiceS;
  this->LevelOccupation = spectrum.LevelOccupation;
  this->NbrParticles = spectrum.NbrParticles;
  this->Flag = spectrum.Flag;
  return *this;
}

// get occupation information from a formatted string
//
// descriptionString = pointer to the string containing the description
// descriptionArray = reference on the array where description has to be stored
// return value = number of particles (0 if an error occured)

int LandauSpectrumOnSphere::ParseOccupationDescription (char* descriptionString, int**& descriptionArray)
{
  int TmpNbrParticles = 0;
  descriptionArray = new int* [this->NbrLandauLevels];
  int Pos = 0;
  for (int i = 0; i < this->NbrLandauLevels; ++i)
    {
      int Lim = this->TwiceS + 1 + (2 * i);
      descriptionArray[i] = new int [Lim];
      int Pos2 = 0;
      while ((Pos2 < Lim) && (descriptionString[Pos] != '\0')) 
	{
	  if ((descriptionString[Pos] == ' ') || (descriptionString[Pos] == '\t') || (descriptionString[Pos] == '\n'))
	    {
	      ++Pos;
	    }
	  else
	    {
	      if ((descriptionString[Pos] != '0') &&  (descriptionString[Pos] != '1'))
		{
		  for (int j = 0; j <= i; ++j)
		    delete[] descriptionArray[j];
		  delete[] descriptionArray;
		  return 0;
		}
	      else
		{
		  descriptionArray[i][Pos2] = (int) (descriptionString[Pos] - '0');
		  TmpNbrParticles += descriptionArray[i][Pos2];
		  ++Pos2;
		  ++Pos;
		}
	    }
	}
      while ((descriptionString[Pos] != '\0') && ((descriptionString[Pos] == ' ') || (descriptionString[Pos] == '\t') || (descriptionString[Pos] == '\n')))
	{
	  ++Pos;
	}
      if (((i + 1) != this->NbrLandauLevels) && ((Pos2 != Lim) || (descriptionString[Pos] != '|')))
	{
	  for (int j = 0; j <= i; ++j)
	    delete[] descriptionArray[j];
	  delete[] descriptionArray;
	  return 0;
	}
      else
	{
	  ++Pos;
	}
    }
  return TmpNbrParticles;
}

// print Landau level spectrum
//
// str = reference on the output stream
// return value = reference on the output stream

ostream& LandauSpectrumOnSphere::PrintSpectrum(ostream& str)
{
  int MaxMomentum = this->TwiceS + 2 * (this->NbrLandauLevels - 1);
  for (int j = this->NbrLandauLevels - 1; j >= 0; --j)
    {
      for (int i = 0; i < (this->NbrLandauLevels - 1 - j); ++i)
	str << "  ";
      for (int i = 0; i <= MaxMomentum; ++i)
	str << this->LevelOccupation[j][i] << " ";
      str << endl;
      MaxMomentum -= 2;
    }
  return str;
}

// evaluate number of changes needed to go from one spectrum to another
//
// spectrum = reference on the spectrum to compare with
// return value = number of changes (-1 if an error occurs)

int LandauSpectrumOnSphere::EvaluateDistance(LandauSpectrumOnSphere& spectrum)
{
  if (((this->NbrLandauLevels != spectrum.NbrLandauLevels) || (this->TwiceS != spectrum.TwiceS) ||
       (this->NbrParticles != spectrum.NbrParticles)) && (this->NbrLandauLevels != 0))
    {
      return -1;
    }
  int Distance = 0;
  int MaxMomentum = this->TwiceS;
  for (int j = 0; j < this->NbrLandauLevels; ++j)
    {
      for (int i = 0; i <= MaxMomentum; ++i)
	{
	  if (this->LevelOccupation[j][i] > spectrum.LevelOccupation[j][i])
	    Distance += (this->LevelOccupation[j][i] - spectrum.LevelOccupation[j][i]);
	  else
	    if (this->LevelOccupation[j][i] < spectrum.LevelOccupation[j][i])
	      Distance += (spectrum.LevelOccupation[j][i] - this->LevelOccupation[j][i]);
	}
       MaxMomentum += 2;
   }
  return (Distance >> 1);
}
