////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  set of functions used to manage file names                //
//                   related to the FQHE on cylinder geometry                 //
//                                                                            //
//                        last modification : 08/06/2016                      //
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


#include "Tools/FQHEFiles/FQHEOnCylinderFileTools.h"

#include <iostream>
#include <cstring>
#include <cstdlib>


using std::cout;
using std::endl;


// try to guess system information from file name
//
// filename = file name
// nbrParticles = reference to the number of particles 
// lzMax = reference to twice the maximum momentum for a single particle
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons)
// ratio = reference to the cylinder aspect ratio
// perimeter = reference to the cylinder perimeter
// return value = true if no error occured

bool FQHEOnCylinderFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& lzMax, bool& statistics,
					      double& ratio, double& perimeter)
{
  char* StrNbrParticles;
  StrNbrParticles = strstr(filename, "_n_");
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
      cout << "can't guess the number of particles from file name " << filename << endl;
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
	  lzMax = atoi(StrNbrParticles);
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
	  cout << "can't guess the particle statistics from file name " << filename << endl;
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
  ratio = 0.0;
  perimeter = 0.0;
  StrNbrParticles = strstr(filename, "_perimeter_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 11;
      perimeter = atof(StrNbrParticles);
    }
  StrNbrParticles = strstr(filename, "_ratio_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 7;
      ratio = atof(StrNbrParticles);
    }
  if ((ratio == 0.0) && (perimeter == 0.0))
    {
      cout << "can't guess the cylinder perimeter of aspect ratio from file name " << filename << endl;    
      return false;
    }
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles
// lzMax = reference to twice the maximum momentum for a single particle
// ky = reference to the momentum along the cylinder perimeter
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// ratio = reference to the cylinder aspect ratio
// perimeter = reference to the cylinder perimeter
// return value = true if no error occured

bool FQHEOnCylinderFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, bool& statistics,
						    double& ratio, double& perimeter)
{
  if (FQHEOnCylinderFindSystemInfoFromFileName(filename, nbrParticles, lzMax, statistics, ratio, perimeter) == false)
    return false;
  char* StrNbrParticles;
  StrNbrParticles = strstr(filename, "_lz_");
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
	  lz = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '.';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess the momentum from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name for system with an SU(2) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles
// lzMax = reference to twice the maximum momentum for a single particle
// ky = reference to the momentum along the cylinder perimeter
// sz = reference to twice the z projection of the total spin
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// ratio = reference to the cylinder aspect ratio
// perimeter = reference to the cylinder perimeter
// return value = true if no error occured

bool FQHEOnCylinderWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& ky, int& sz, bool& statistics,
							    double& ratio, double& perimeter)
{
  if (FQHEOnCylinderFindSystemInfoFromVectorFileName(filename, nbrParticles, lzMax, ky, statistics, ratio, perimeter) == false)
    {
      return false;
    }
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
      cout << "can't guess the z projection of the total spin from file name " << filename << endl;
      return false;            
    }
  return true;
}

