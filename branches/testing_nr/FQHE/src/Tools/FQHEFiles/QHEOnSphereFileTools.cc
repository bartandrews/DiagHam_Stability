////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    set of functions used to managed files related to QHE on sphere         //
//                                                                            //
//                        last modification : 06/06/2006                      //
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


#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>


using std::cout;
using std::endl;


// try to guess system information from file name
//
// filename = file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSphereFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& lzMax, bool& statistics)
{
  char* StrNbrParticles;
  if (nbrParticles == 0)
    {
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
	  cout << "can't guess number of particles from file name " << filename << endl;
	  return false;            
	}
    }
  if (lzMax == 0)
    {
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
    }
  if (statistics == true)
    {
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
    }
  return true;
}


// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSphereFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, bool& statistics)
{
  if (FQHEOnSphereFindSystemInfoFromFileName(filename, nbrParticles, lzMax, statistics) == false)
    return false;
  char* StrNbrParticles;
  if (lz == 0)
    {
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
	  cout << "can't guess z projection of the angular momentum from file name " << filename << endl;
	  return false;            
	}
    }
  return true;
}

// try to guess system information from file name for system suth an SU(2) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// sz = reference to twice the z projection of the total spin (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, int& sz, bool& statistics)
{
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(filename, nbrParticles, lzMax, lz, statistics) == false)
    return false;
  char* StrNbrParticles;
  if (sz == 0)
    {
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
	  cout << "can't guess z projection of the total spin from file name " << filename << endl;
	  return false;            
	}
    }
  return true;
}


// try to guess system information from file name for system with an SU(2) degree of freedom and discrete symmetries
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// sz = reference to twice the z projection of the total spin (grab it only if initial value is 0)
// szSymmetry = reference on the flag for the Sz<->-Sz symmetry
// szSymmetryMinusParity = reference on the flag for the minus parity sector of the Sz<->-Sz symmetry
// lzSymmetry = reference on the flag for the Lz<->-Lz symmetry
// lzSymmetryMinusParity = reference on the flag for the minus parity sector of the Lz<->-Lz symmetry
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, int& sz, 
							  bool& szSymmetry, bool& szSymmetryMinusParity, bool& lzSymmetry, bool& lzSymmetryMinusParity, bool& statistics)
{
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(filename, nbrParticles, lzMax, lz, sz, statistics) == false)
    return false;
  szSymmetry = false;
  szSymmetryMinusParity = false;
  lzSymmetry = false;
  lzSymmetryMinusParity = false;
  if (strstr(filename, "_szpsym_") != 0)
    {
      szSymmetry = true;
      return true;
    }
  if (strstr(filename, "_szmsym_") != 0)
    {
      szSymmetry = true;
      szSymmetryMinusParity = true;
      return true;
    }
  if (strstr(filename, "_lzpsym_") != 0)
    {
      lzSymmetry = true;
      return true;
    }
  if (strstr(filename, "_lzmsym_") != 0)
    {
      lzSymmetry = true;
      lzSymmetryMinusParity = true;
      return true;
    }
  if (strstr(filename, "_lzpszpsym_") != 0)
    {
      szSymmetry = true;
      lzSymmetry = true;
      return true;
    }
  if (strstr(filename, "_lzmszpsym_") != 0)
    {
      szSymmetry = true;
      lzSymmetry = true;
      lzSymmetryMinusParity = true;
      return true;
    }
  if (strstr(filename, "_lzpszmsym_") != 0)
    {
      szSymmetry = true;
      szSymmetryMinusParity = true;
      lzSymmetry = true;
      return true;
    }
  if (strstr(filename, "_lzmszmsym_") != 0)
    {
      szSymmetry = true;
      szSymmetryMinusParity = true;
      lzSymmetry = true;
      lzSymmetryMinusParity = true;
      return true;
    }
  return true;
  
}


// try to guess system information from file name for system with an SU(3) degree of freedom and discrete symmetries
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// tz = reference to twice the Tz quantum number (grab it only if initial value is 0)
// y = reference to three time the Y quantum number (grab it only if initial value is 0)
// tzSymmetry = reference on the flag for the Tz<->-Tz symmetry
// tzSymmetryMinusParity = reference on the flag for the minus parity sector of the Y<->-Y symmetry
// ySymmetry = reference on the flag for the Z3 symmetry
// lzSymmetry = reference on the flag for the Lz<->-Lz symmetry
// lzSymmetryMinusParity = reference on the flag for the minus parity sector of the Lz<->-Lz symmetry
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FQHEOnSphereWithSU3SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, int& tz, int& y,
							     bool& tzSymmetry, bool& tzSymmetryMinusParity,  bool& ySymmetry, 
							     bool& lzSymmetry, 
							     bool& lzSymmetryMinusParity, bool& statistics)
{
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(filename, nbrParticles, lzMax, lz, statistics) == false)
    return false;
  char* StrNbrParticles;
  if (tz == 0)
    {
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
	      tz = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
    }
  if (y == 0)
    {
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
	      y = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess z projection of the total spin from file name " << filename << endl;
	  return false;            
	}
    }
  tzSymmetry = false;
  tzSymmetryMinusParity = false;
  ySymmetry = false;
  lzSymmetry = false;
  lzSymmetryMinusParity = false;
  if (strstr(filename, "tzp") != 0)
    {
      tzSymmetry = true;
    }
  if (strstr(filename, "tzm") != 0)
    {
      tzSymmetry = true;
      tzSymmetryMinusParity = true;
    }
  if (strstr(filename, "z3sym") != 0)
    {
      ySymmetry = true;
    }
  if (strstr(filename, "lzp") != 0)
    {
      lzSymmetry = true;
    }
  if (strstr(filename, "lzm") != 0)
    {
      lzSymmetry = true;
      lzSymmetryMinusParity = true;
    }
  return true;
}

// try to guess system information from file name for system suth an SU(4) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// sz = reference to twice the z projection of the total spin (grab it only if initial value is 0)
// iz = reference to twice the z projection of the total isospin (grab it only if initial value is 0)
// pz = reference to twice the entanglement (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSphereWithSU4SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, int& sz, int& iz, int& pz, bool& statistics)
{
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(filename, nbrParticles, lzMax, lz, sz, statistics) == false)
    return false;
  char* StrNbrParticles;
  if (iz == 0)
    {
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
	      sz = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess z projection of the total isospin from file name " << filename << endl;
	  return false;            
	}
    }
  if (pz == 0)
    {
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
	      sz = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess the entropy from file name " << filename << endl;
	  return false;            
	}
    }
  return true;
}

