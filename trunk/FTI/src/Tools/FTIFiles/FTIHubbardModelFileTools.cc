////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     set of functions used to managed files related to Hubbard models       //
//                                                                            //
//                        last modification : 19/06/2014                      //
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


#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include <iostream>
#include <cstring>
#include <cstdlib>


using std::cout;
using std::endl;


// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference to the number sites
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// return value = true if no error occured

bool FTIHubbardModelFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& nbrSites, bool& statistics)
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
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  nbrParticles = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
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
  StrNbrParticles = strstr(filename, "_ns_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  nbrSites = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      nbrSites = 0;
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
// nbrSites = reference to the number sites
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference to flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured

bool FTIHubbardModelFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, bool& statistics, bool& gutzwiller)
{
  if (FTIHubbardModelFindSystemInfoFromFileName(filename, nbrParticles, nbrSites, statistics) == false)
    {
      cout << "basic" << endl;
      return false;
    }
  char* GutzwillerFlag = strstr(filename, "_gutzwiller_");
  if (GutzwillerFlag != 0)
    gutzwiller = true;
  else
    gutzwiller = false;
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference to the number sites
// szValue = reference on the value of the total spin
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference to flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured

bool FTIHubbardModelFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, bool& statistics, bool& gutzwiller)
{
  if (FTIHubbardModelFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSites, statistics, gutzwiller) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  int SizeString;
  StrNbrParticles = strstr(filename, "_sz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (((StrNbrParticles[SizeString] >= '0') 
																	&& (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] != '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  szValue = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
    
  if (StrNbrParticles == 0)
    {
      return false;            
    }
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrSites = reference to the number sites
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// parity = reference to the particle number parity
// return value = true if no error occured

bool FTIHubbardModelFindSystemInfoFromVectorFileName(char* filename, int& nbrSites, bool& statistics, int& parity)
{
  char* StrNbrParticles;
  StrNbrParticles = strstr(filename, "_par_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 5;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.'))&& (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  parity = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess particle number parity from file name " << filename << endl;
      return false;            
    }
  StrNbrParticles = strstr(filename, "_ns_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  nbrSites = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      nbrSites = 0;
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
// nbrSites = reference on the number sites
// xMomentum = reference on the momentum sector in the x direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWith1DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& xMomentum, int& xPeriodicity, bool& statistics, bool& gutzwiller)
{
  if (FTIHubbardModelFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSites, statistics, gutzwiller) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  StrNbrParticles = strstr(filename, "_kx_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  xMomentum = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess momentum sector from file name " << filename << endl;
      return false;            
    }
  StrNbrParticles = strstr(filename, "_x_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  xPeriodicity = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
//       cout << "can't guess periodicity from file name " << filename << endl;
      return false;            
    }
  return true;
}



// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// yPeriodicity = reference on the periodicity in the y direction with respect to site numbering 
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& xMomentum, int& yMomentum, int& xPeriodicity, int& yPeriodicity, bool& statistics, bool& gutzwiller)
{
  if (FTIHubbardModelWith1DTranslationFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSites, xMomentum, xPeriodicity, statistics, gutzwiller) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  int SizeString;
    
  StrNbrParticles = strstr(filename, "_ky_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  yMomentum = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess y momentum sector from file name " << filename << endl;
      return false;            
    }
  StrNbrParticles = strstr(filename, "_y_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  yPeriodicity = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
//       cout << "can't guess periodicity from file name " << filename << endl;
      return false;            
    }
  return true;
}



// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// szValue = reference on the value of the total spin
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, bool& statistics, bool& gutzwiller)
{    
  if (FTIHubbardModelFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSites, statistics, gutzwiller) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  int SizeString;
  StrNbrParticles = strstr(filename, "_sz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (((StrNbrParticles[SizeString] >= '0') 
																	&& (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] != '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  szValue = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
    
  if (StrNbrParticles == 0)
    {
//       cout << "can't guess sz value sector from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// szValue = reference on the value of the total spin
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, int& szSymmetry, bool& statistics, bool& gutzwiller)
{
  if (FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSites, szValue, szSymmetry, statistics, gutzwiller) == false)
    {
      return false;
    } 
  char* StrNbrParticles;
  int SizeString;
  StrNbrParticles = strstr(filename, "_sz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (((StrNbrParticles[SizeString] >= '0') 
																	&& (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] != '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  szValue = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
    
  if (StrNbrParticles == 0)
    {
      return false;            
    }
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// szValue = reference on the value of the total spin
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// yPeriodicity = reference on the periodicity in the y direction with respect to site numbering
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured

bool FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, 
								      int& xMomentum, int& yMomentum, int& xPeriodicity, int& yPeriodicity, 
								      bool& statistics, bool& gutzwiller)
{
  if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSites, xMomentum, yMomentum, 
								       xPeriodicity, yPeriodicity, statistics, gutzwiller) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  int SizeString;
  StrNbrParticles = strstr(filename, "_sz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (((StrNbrParticles[SizeString] >= '0') 
																	&& (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] != '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  szValue = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
    
  if (StrNbrParticles == 0)
    {
      return false;            
    }
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// szValue = reference on the value of the total spin
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// yPeriodicity = reference on the periodicity in the y direction with respect to site numbering
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, int& szSymmetry, int& xMomentum, int& yMomentum, int& xPeriodicity, int& yPeriodicity, bool& statistics, bool& gutzwiller)
{
  if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSites, szValue, xMomentum, yMomentum, 
								       xPeriodicity, yPeriodicity, statistics, gutzwiller) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  szSymmetry = 0;
  StrNbrParticles = strstr(filename, "_szsym_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 7;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && 
	     (((StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] == '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  szSymmetry = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	{
	  cout << "error while retrieving the Sz<->-Sz parity" << endl;
	  return false;
	}
    }
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// szValue = reference on the value of the total spin
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// minNbrSinglets = minimum number of on-site singlets
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// yPeriodicity = reference on the periodicity in the y direction with respect to site numbering
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, int& szSymmetry, int& minNbrSinglets, 
								      int& xMomentum, int& yMomentum, int& xPeriodicity, int& yPeriodicity, bool& statistics, bool& gutzwiller)
{
  if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSites, szValue, szSymmetry, xMomentum, yMomentum, 
								       xPeriodicity, yPeriodicity, statistics, gutzwiller) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  minNbrSinglets = 0;
  StrNbrParticles = strstr(filename, "_minnbrsinglet_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 15;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && 
	     (((StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] == '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  minNbrSinglets = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	{
	  cout << "error while retrieving the minimum number of on-site singlet" << endl;
	  return false;
	}
    }
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrUnitCellX = reference on the number of unit cell  in the x direction
// nbrUnitCellY = reference on the number of unit cell  in the y direction
// nbrSiteInUnitCellX  = reference on the number of site in each unit cell  in the x direction
// nbrSiteInUnitCellY  = reference on the number of site in each unit cell  in the y direction
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured

bool FTIHofstadterdModelFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrUnitCellX, int& nbrUnitCellY, int& nbrSiteInUnitCellX,int& nbrSiteInUnitCellY, bool& statistics, bool& gutzwiller)
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
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  nbrParticles = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
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



  StrNbrParticles = strstr(filename, "_x_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  nbrUnitCellX = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess the number of unit cell in x direction " << filename << endl;
      return false;            
    }

  StrNbrParticles = strstr(filename, "_y_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  nbrUnitCellY = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess the number of unit cell in y direction " << filename << endl;
      return false;            
    }



  StrNbrParticles = strstr(filename, "_X_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  nbrSiteInUnitCellX = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess number of site in the magnetic cell in x direction  from file name " << filename << endl;
      return false;            
    }

  StrNbrParticles = strstr(filename, "_Y_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  nbrSiteInUnitCellY = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess number of site in the magnetic cell in y direction  from file name " << filename << endl;
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
  
  char* GutzwillerFlag = strstr(filename, "_gutzwiller_");
  if (GutzwillerFlag != 0)
    gutzwiller = true;
  else
    gutzwiller = false;
  return true;
  
}




// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// nbrUnitCellX = reference on the number of unit cell  in the x direction
// nbrUnitCellY = reference on the number of unit cell  in the y direction
// nbrSiteInUnitCellX  = reference on the number of site in each unit cell  in the x direction
// nbrSiteInUnitCellY  = reference on the number of site in each unit cell  in the y direction
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured

bool FTIHofstadterdModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& xMomentum, int& yMomentum,  int& nbrUnitCellX, int& nbrUnitCellY, int& nbrSiteInUnitCellX,int& nbrSiteInUnitCellY, bool& statistics, bool& gutzwiller)
{
  if (FTIHofstadterdModelFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrUnitCellX, nbrUnitCellY, nbrSiteInUnitCellX, nbrSiteInUnitCellY, statistics, gutzwiller) == false) 
    {
      return false;
    }

  char* StrNbrParticles;
  int SizeString;
    
  StrNbrParticles = strstr(filename, "_kx_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  xMomentum = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess x momentum sector from file name " << filename << endl;
      return false;            
    }

  StrNbrParticles = strstr(filename, "_ky_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  yMomentum = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess y momentum sector from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// fluxInUnitCell, = reference on the number of flux in the unit cell
// nbrUnitCellX = reference on the number of unit cell  in the x direction
// nbrUnitCellY = reference on the number of unit cell  in the y direction
// nbrSiteInUnitCellX  = reference on the number of site in each unit cell  in the x direction
// nbrSiteInUnitCellY  = reference on the number of site in each unit cell  in the y direction
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// translationFlag = bool indicating if translations are used
// return value = true if no error occured

bool FTIHofstadterdModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& xMomentum, int& yMomentum, int & fluxInUnitCell, int& nbrUnitCellX, int& nbrUnitCellY, int& nbrSiteInUnitCellX,int& nbrSiteInUnitCellY, bool& statistics, bool& gutzwiller,  bool& translationFlag)
{
  if (FTIHofstadterdModelFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrUnitCellX, nbrUnitCellY, nbrSiteInUnitCellX, nbrSiteInUnitCellY, statistics, gutzwiller) == false) 
    {
      return false;
    }

  char* StrNbrParticles;
  int SizeString;

  StrNbrParticles = strstr(filename, "_q_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  fluxInUnitCell = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "Impossible to obtain the number of flux per unit cell" << filename << endl;
      return false;            
    }

  StrNbrParticles = strstr(filename, "_kx_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  xMomentum = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	  translationFlag = true ;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      translationFlag = false;
      cout << "No translations from file name " << filename << endl;
      return true;            
    }

  StrNbrParticles = strstr(filename, "_ky_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  yMomentum = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess y momentum sector from file name " << filename << endl;
      return false;            
    }
  return true;
}


// try to guess system information from file name
//
// filename = vector file name
// szValue = reference on the value of the total spin
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// minNbrSinglets = minimum number of on-site singlets
// return value = true if no error occured

bool FTIHofstadterModelWithSzFindSystemInfoFromVectorFileName(char* filename,int & szValue, int& szSymmetry, int& minNbrSinglets)
{    
  char* StrNbrParticles;
  int SizeString;
  StrNbrParticles = strstr(filename, "_sz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (((StrNbrParticles[SizeString] >= '0') 
																	&& (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] != '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  szValue = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
    
  if (StrNbrParticles == 0)
    {
      cout << "can't guess sz value sector from file name " << filename << endl;
      return false;            
    }
  szSymmetry = 0;
  StrNbrParticles = strstr(filename, "_szsym_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 7;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && 
	     (((StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] == '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  szSymmetry = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	{
	  cout << "error while retrieving the Sz<->-Sz parity" << endl;
	  return false;
	}
    }
  minNbrSinglets = 0;
  StrNbrParticles = strstr(filename, "_minnbrsinglet_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 15;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && 
	     (((StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] == '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  minNbrSinglets = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	{
	  cout << "error while retrieving the minimum number of on-site singlet" << endl;
	  return false;
	}
    }
  return true;
}



// try to guess system information from file name
//
// filename = vector file name
// szValue = reference on the value of the total spin
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// minNbrSinglets = minimum number of on-site singlets
// usingSzSymmetryFlag = bool indicating if the SzSymmetry is used
// usingNbrSingletConstraintFlag = bool indicating if the constraint on the number of singlets is used
// return value = true if no error occured

bool FTIHofstadterModelWithSzFindSystemInfoFromVectorFileName(char* filename,int & szValue, int& szSymmetry, int& minNbrSinglets, bool & usingSzSymmetryFlag, bool & usingNbrSingletConstraintFlag)
{    
  char* StrNbrParticles;
  int SizeString;
  StrNbrParticles = strstr(filename, "_sz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (((StrNbrParticles[SizeString] >= '0') 
																	&& (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] != '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  szValue = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  
  if (StrNbrParticles == 0)
    {
      cout << "can't guess sz value sector from file name " << filename << endl;
      return false;            
    }
  szSymmetry = 0;
  StrNbrParticles = strstr(filename, "_szsym_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 7;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && 
	     (((StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] == '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  szSymmetry = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	  usingSzSymmetryFlag = true;
	}
      else
	{
	  cout << "No Sz<->-Sz parity" << endl;
	  usingSzSymmetryFlag = false;
	}
    }
  minNbrSinglets = 0;
  StrNbrParticles = strstr(filename, "_minnbrsinglet_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 15;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && 
	     (((StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] == '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  minNbrSinglets = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	  usingNbrSingletConstraintFlag = true;
	}
      else
	{
	  cout << "No minimum number of on-site singlet" << endl;
	  usingNbrSingletConstraintFlag = false;
	}
    }
  return true;
}


// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// szValue = reference on the value of the total spin
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// yPeriodicity = reference on the periodicity in the y direction with respect to site numbering
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, int& szSymmetry, int& xMomentum, int& yMomentum, int& xPeriodicity, int& yPeriodicity, bool& statistics, bool& gutzwiller, bool& clusterExclusion)
{
  if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSites, szValue, xMomentum, yMomentum, 
								       xPeriodicity, yPeriodicity, statistics, gutzwiller, clusterExclusion) == false)
    {
      return false;
    }
    
  char* StrNbrParticles;
  int SizeString;
  StrNbrParticles = strstr(filename, "_szp_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 5;
      SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (( (StrNbrParticles[SizeString] <= '9'))))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  szSymmetry = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
      {
	cout << "error while retrieving the Sz<->-Sz parity" << endl;
	return false;
      }
    }
  else
  {
    cout << "Not using Sz symmetry" << endl;
    return false;
  }
    
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// szValue = reference on the value of the total spin
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// xMomentum = reference on the momentum sector in the x direction
// yMomentum = reference on the momentum sector in the y direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// yPeriodicity = reference on the periodicity in the y direction with respect to site numbering
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& szValue, int& xMomentum, int& yMomentum, int& xPeriodicity, int& yPeriodicity, bool& statistics, bool& gutzwiller, bool& clusterExclusion)
{
  if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSites, szValue, xMomentum, yMomentum, 
								       xPeriodicity, yPeriodicity, statistics, gutzwiller) == false)
    {
      return false;
    }
  char* ClusterExclusionFlag = strstr(filename, "_clustercharging_exclusion_");
  if (ClusterExclusionFlag != 0)
    clusterExclusion = true;
  else
    clusterExclusion = false;
  return true;
}