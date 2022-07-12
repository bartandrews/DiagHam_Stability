////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    set of functions used to managed files related to QHE on lattice         //
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


#include "Tools/FQHEFiles/QHEOnLatticeFileTools.h"

#include <iostream>


using std::cout;
using std::endl;



// try to guess system information from file name
//
// filename = file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lx, ly = indication of lattice geometry
// interaction = strength of interaction parameter
// flux = number of flux quanta
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// hardcore = returns true if hardcore bosons encountered
// return value = true if no error occured
bool FQHEOnLatticeFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& lx, int& ly, double &interaction, int &flux, bool& statistics, bool& hardcore)
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
  if (lx == 0)
    {
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
	      StrNbrParticles[SizeString] = '\0';
	      lx = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess length lx from file name " << filename << endl;
	  return false;            
	}
    }
  if (ly == 0)
    {
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
	      StrNbrParticles[SizeString] = '\0';
	      ly = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess length ly from file name " << filename << endl;
	  return false;            
	}
    }
  StrNbrParticles = strstr(filename, "hardcore");
  if (StrNbrParticles != 0)
    {
      interaction=0.0;
      hardcore=true;
    }
  if ((interaction == 0.0)&&(hardcore=false))
    {
      StrNbrParticles = strstr(filename, "_u_");
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
	      interaction = strtod(StrNbrParticles,NULL);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess interaction u from file name " << filename << endl;
	  return false;            
	}
    }
  
  if (flux == 0)
    {
      StrNbrParticles = strstr(filename, "_q_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 3;
	  int SizeString = 0;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if (SizeString != 0)
	    {
	      char tmpC = StrNbrParticles[SizeString];
	      StrNbrParticles[SizeString] = '\0';
	      flux = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = tmpC;
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess flux from file name " << filename << endl;
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
// lx, ly = indication of lattice geometry
// interaction = strength of interaction parameter
// flux = number of flux quanta
// nbrState = number of eigenstate
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// hardcore = returns true if hardcore bosons encountered
// return value = true if no error occured
bool FQHEOnLatticeFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lx, int& ly, double &interaction, int &flux, int &nbrState, bool& statistics, bool &hardcore)
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
  if (lx == 0)
    {
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
	      StrNbrParticles[SizeString] = '\0';
	      lx = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess length lx from file name " << filename << endl;
	  return false;            
	}
    }
  if (ly == 0)
    {
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
	      StrNbrParticles[SizeString] = '\0';
	      ly = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess length ly from file name " << filename << endl;
	  return false;            
	}
    }
  StrNbrParticles = strstr(filename, "hardcore");
  if (StrNbrParticles != 0)
    {
      interaction=0.0;
      hardcore=true;
    }
  if ((interaction == 0)&& (hardcore==false))
    {
      StrNbrParticles = strstr(filename, "_u_");
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
	      interaction = strtod(StrNbrParticles,NULL);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess interaction from file name " << filename << endl;
	  return false;            
	}
    }
  // definitely search for _q_ even if q!=0
    {
      StrNbrParticles = strstr(filename, "_q_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 3;
	  int SizeString = 0;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9') && (StrNbrParticles[SizeString] >= '0'))
	    ++SizeString;
	  if (SizeString != 0)
	    {
	      char tmpC=StrNbrParticles[SizeString];
	      StrNbrParticles[SizeString] = '\0';
	      if (flux == 0)
		{
		  flux = atoi(StrNbrParticles);
		}
	      StrNbrParticles[SizeString] = tmpC;
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess flux from file name " << filename << endl;
	  return false;
	}
      if (nbrState==0)
	{
	  // parse file number
	  StrNbrParticles++;
	  
	  int SizeString = 0;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.')
		 && (StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '.') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      nbrState = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '.';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	  
	  if (StrNbrParticles == 0)
	    {
	      cout << "can't guess state number from file name " << filename << endl;
	      return false;
	    }
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
// lx, ly = indication of lattice geometry
// interaction = strength of interaction parameter
// kyMomentum = momentum in y-direction when using symmetries
// flux = number of flux quanta
// nbrState = number of eigenstate
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// hardcore = returns true if hardcore bosons encountered
// return value = true if no error occured
bool FQHEOnLatticeFindSystemInfoWithKyFromVectorFileName(char* filename, int& nbrParticles, int& lx, int& ly, int &kyMomentum, double &interaction, int &flux, int &nbrState, bool& statistics, bool &hardcore)
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
  if (lx == 0)
    {
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
	      StrNbrParticles[SizeString] = '\0';
	      lx = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess length lx from file name " << filename << endl;
	  return false;            
	}
    }
  if (ly == 0)
    {
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
	      StrNbrParticles[SizeString] = '\0';
	      ly = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess length ly from file name " << filename << endl;
	  return false;            
	}
    }
  StrNbrParticles = strstr(filename, "hardcore");
  if (StrNbrParticles != 0)
    {
      interaction=0.0;
      hardcore=true;
    }
  if ((interaction == 0) && (hardcore==false))
    {
      StrNbrParticles = strstr(filename, "_u_");
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
	      interaction = strtod(StrNbrParticles,NULL);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess interaction from file name " << filename << endl;
	  return false;            
	}
    }
  if (kyMomentum == 0)
    {
      StrNbrParticles = strstr(filename, "_k_");
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
	      kyMomentum = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "could not guess momentum ky from file name " << filename << endl;
	  return false;
	}
    }
  // definitely search for _q_ even if q!=0
    {
      StrNbrParticles = strstr(filename, "_q_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 3;
	  int SizeString = 0;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if (SizeString != 0)
	    {
	      char tmpC=StrNbrParticles[SizeString];
	      StrNbrParticles[SizeString] = '\0';
	      if (flux == 0)
		{
		  flux = atoi(StrNbrParticles);
		}
	      StrNbrParticles[SizeString] = tmpC;
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess flux from file name " << filename << endl;
	  return false;
	}
      if (nbrState==0)
	{
	  // parse file number
	  StrNbrParticles++;
	  
	  int SizeString = 0;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.')
		 && (StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '.') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      nbrState = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '.';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	  
	  if (StrNbrParticles == 0)
	    {
	      cout << "can't guess state number from file name " << filename << endl;
	      return false;
	    }
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
