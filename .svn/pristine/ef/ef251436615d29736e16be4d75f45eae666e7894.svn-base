////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     set of functions used to managed files related to QHE on torus         //
//                                                                            //
//                        last modification : 12/04/2010                      //
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


#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"
#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cstring>
#include <cstdlib>


using std::cout;
using std::endl;


// try to guess system information from file name
//
// filename = file name
// nbrParticles = reference to the number of particles
// kyMax = reference to the maximum momentum for a single particle
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured

bool FQHEOnTorusFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& kyMax, bool& statistics)
{
  char* StrNbrParticles;

  StrNbrParticles = strstr(filename, "_n_");
  if (StrNbrParticles == 0)
    StrNbrParticles = strstr(filename, "_p_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  nbrParticles = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess number of particles from file name " << filename << endl;
      return false;            
    }

  StrNbrParticles = strstr(filename, "_2s_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  kyMax = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
	  else
	    StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess maximum momentum from file name " << filename << endl;
      return false;            
    }

  if (strstr(filename, "fermion") == 0)
    {
      if (strstr(filename, "boson") == 0)
	{
	  cout << "can't guess particle statistics from file name " << filename << endl;
	  return false;	  
	}
      else
	{
	  statistics = false;
	}
    }
  else
    {
      statistics = true;
    }
  return true;
}


// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// kyMax = reference to the maximum momentum for a single particle
// ky = reference to the y projection of the angular momentum
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured

bool FQHEOnTorusFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& ky, bool& statistics)
{
  if (FQHEOnTorusFindSystemInfoFromFileName(filename, nbrParticles, kyMax, statistics) == false)
    return false;
  char* StrNbrParticles;

  StrNbrParticles = strstr(filename, "_ky_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      if (StrNbrParticles[SizeString] == '-')
	++SizeString;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '.') && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  ky = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '.';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess ky momentum from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// kyMax = reference to the momentum for a single particle
// kx = reference to the x momentum
// ky = reference to the y momentum
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured

bool FQHEOnTorusFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& kx, int& ky, bool& statistics)
{
  if (FQHEOnTorusFindSystemInfoFromVectorFileName(filename, nbrParticles, kyMax, ky, statistics) == false)
    return false;
  char* StrNbrParticles;

  StrNbrParticles = strstr(filename, "_kx_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      if (StrNbrParticles[SizeString] == '-')
	++SizeString;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '.') || (StrNbrParticles[SizeString] == '_')) && (SizeString != 0))
	{
	  char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  kx = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess kx momentum from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// kyMax = reference to the momentum for a single particle
// kx = reference to the x momentum
// ky = reference to the y momentum
// ratio = reference on the aspect ratio
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured

bool FQHEOnTorusFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& kx, int& ky, double& ratio, bool& statistics)
{
  if (FQHEOnTorusFindSystemInfoFromVectorFileName(filename, nbrParticles, kyMax, kx, ky, statistics) == false)
    {
      return false;
    }
  char* StrRatio = strstr(filename, "_ratio_");
  if (StrRatio != 0)
    {
      StrRatio += 7;
      ratio = atof(StrRatio);
    }
  else
    return false;
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// kyMax = reference to the momentum for a single particle
// ky = reference to the y momentum
// ratio = reference on the aspect ratio
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured

bool FQHEOnTorusFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& ky, double& ratio, bool& statistics)
{
  if (FQHEOnTorusFindSystemInfoFromVectorFileName(filename, nbrParticles, kyMax, ky, statistics) == false)
    {
      return false;
    }
  char* StrRatio = strstr(filename, "_ratio_");
  if (StrRatio != 0)
    {
      StrRatio += 7;
      ratio = atof(StrRatio);
    }
  else
    return false;
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// NbrParticles = reference to the number of particles 
// NbrFluxQuanta (kyMax) = reference to the momentum for a single particle
// Momentum (ky) = reference to the y momentum
// Ratio = reference on the aspect ratio
// Statistics = reference to flag for fermionic statistics (if set to false on input, do not parse from filename)
// return value = true if no error occured

bool FQHEOnTorusFindSystemInfoFromVectorFileName_SpectralResponse(char* filename, int& NbrParticles, int& NbrFluxQuanta, int& Momentum, double& Ratio, bool& Statistics)
{
  FilenameStatisticsCheck(Statistics, filename);
  FilenameIntegerSearch(NbrParticles, filename, "_n_");
  FilenameIntegerSearch(NbrFluxQuanta, filename, "_2s_");
  FilenameIntegerSearch(Momentum, filename, "_ky_");
  FilenameDoubleSearch(Ratio, filename, "_ratio_");
  return true;
}

// try to guess system information from file name for system suth an SU(2) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles
// kyMax = reference to the maximum momentum for a single particle
// ky = reference to the y projection of the angular momentum
// sz = reference to twice the z projection of the total spin
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured

bool FQHEOnTorusWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& ky, int& sz, bool& statistics)
{
  if (FQHEOnTorusFindSystemInfoFromVectorFileName(filename, nbrParticles, kyMax, ky, statistics) == false)
    return false;
  char* StrNbrParticles;

  StrNbrParticles = strstr(filename, "_sz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      if (StrNbrParticles[SizeString] == '-')
	++SizeString;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  sz = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess sz from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name for system suth an SU(2) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles
// kyMax = reference to the maximum momentum for a single particle
// kx = reference to the x momentum
// ky = reference to the y projection of the angular momentum
// sz = reference to twice the z projection of the total spin
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnTorusWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& kx, int& ky, int& sz, bool& statistics)
{
  if (FQHEOnTorusFindSystemInfoFromVectorFileName(filename, nbrParticles, kyMax, kx, ky, statistics) == false)
    return false;
  char* StrNbrParticles;

  StrNbrParticles = strstr(filename, "_sz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      if (StrNbrParticles[SizeString] == '-')
	++SizeString;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  sz = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess sz from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name for system suth an SU(3) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles
// kyMax = reference to the maximum momentum for a single particle
// ky = reference to the y projection of the angular momentum
// totalTz = reference to twice the total Tz value
// totalY = reference to three time the total Y value
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured

bool FQHEOnTorusWithSU3SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& ky, int& totalTz, int& totalY, bool& statistics)
{
  if (FQHEOnTorusFindSystemInfoFromVectorFileName(filename, nbrParticles, kyMax, ky, statistics) == false)
    return false;
  char* StrNbrParticles;

  StrNbrParticles = strstr(filename, "_tz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      if (StrNbrParticles[SizeString] == '-')
	++SizeString;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  totalTz = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess sz from file name " << filename << endl;
      return false;            
    }
  StrNbrParticles = strstr(filename, "_y_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      if (StrNbrParticles[SizeString] == '-')
	++SizeString;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  totalY = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess y from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name for system suth an SU(7) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles
// kyMax = reference to the maximum momentum for a single particle
// ky = reference to the y projection of the angular momentum
// totalSz = reference to twice the total spin value
// totalIz = reference to the total isospin value
// totalPz = reference to the total entanglement value
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured

bool FQHEOnTorusWithSU4SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& ky, int& totalSz, int& totalIz, int& totalPz, bool& statistics)
{
  if (FQHEOnTorusWithSpinFindSystemInfoFromVectorFileName(filename, nbrParticles, kyMax, ky, totalSz, statistics) == false)
    return false;
  char* StrNbrParticles;

  StrNbrParticles = strstr(filename, "_iz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      if (StrNbrParticles[SizeString] == '-')
	++SizeString;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  totalIz = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess iz from file name " << filename << endl;
      return false;            
    }
  StrNbrParticles = strstr(filename, "_pz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      if (StrNbrParticles[SizeString] == '-')
	++SizeString;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  totalPz = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess pz from file name " << filename << endl;
      return false;            
    }
  return true;
}

