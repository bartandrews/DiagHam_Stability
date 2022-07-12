////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             set of functions used to managed files related to FQHE         //
//                          on square lattice                                 //
//                                                                            //
//                        last modification : 14/10/2016                      //
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


#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"
#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cstring>
#include <cstdlib>

using std::cout;
using std::endl;


// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSquareLatticeFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, bool& statistics)
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
	  nbrSiteX = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      nbrSiteX = 0;
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
	  nbrSiteY = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      nbrSiteY = 0;
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
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumX = reference to the momentum along the y direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& momentumX, int& momentumY, bool& statistics)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromFileName(filename, nbrParticles, nbrSiteX, nbrSiteY, statistics) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  StrNbrParticles = strstr(filename, "_kx_");
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
	  momentumX = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      momentumX = 0;
    }
  StrNbrParticles = strstr(filename, "_ky_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      char TmpChar = StrNbrParticles[SizeString];
      if (((StrNbrParticles[SizeString] == '.') || (StrNbrParticles[SizeString] == '_')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  momentumY = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      momentumY = 0;
    }
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction (=-1 if undefined because there are off block diagonal elements. If we are in the block diagonal case of Wannier, then momentumX gives the good sector of the vector)
// momentumY = reference to the momentum along the y direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSquareLatticeWannierFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& momentumX, int& momentumY, bool& statistics)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromFileName(filename, nbrParticles, nbrSiteX, nbrSiteY, statistics) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  StrNbrParticles = strstr(filename, "_kx_");
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
	  momentumX = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      momentumX = -1;
    }
  StrNbrParticles = strstr(filename, "_ky_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      char TmpChar = StrNbrParticles[SizeString];
      if (((StrNbrParticles[SizeString] == '.') || (StrNbrParticles[SizeString] == '_')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  momentumY = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      momentumY = 0;
    }
  return true;
}


// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumX = reference to the momentum along the y direction
// totalSz = reference to the Sz value
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& momentumX, int& momentumY, int& totalSz, bool& statistics)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSiteX, nbrSiteY, momentumX, momentumY, statistics) == false)
    {
      return false;
    }
  char* StrNbrParticles = strstr(filename, "_sz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      char TmpChar = StrNbrParticles[SizeString];
      if (((StrNbrParticles[SizeString] == '.') || (StrNbrParticles[SizeString] == '_')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  totalSz = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      totalSz  = 0;
    }
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumX = reference to the momentum along the y direction
// totalSz = reference to the Sz value
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& totalSz, bool& statistics)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromFileName(filename, nbrParticles, nbrSiteX, nbrSiteY, statistics) == false)
    {
      return false;
    }
  char* StrNbrParticles = strstr(filename, "_sz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      char TmpChar = StrNbrParticles[SizeString];
      if (((StrNbrParticles[SizeString] == '.') || (StrNbrParticles[SizeString] == '_')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  totalSz = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      totalSz  = 0;
    }
  return true;
}

// try to guess system information from file name 
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumX = reference to the momentum along the y direction
// mass = reference to the mass term
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& momentumX, int& momentumY, double& mass, bool& statistics)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSiteX, nbrSiteY, momentumX, momentumY, statistics) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  StrNbrParticles = strstr(filename, "_m_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (SizeString != 0)
	{
	  char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  mass = atof(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      mass = 0.0;
    }
  return true;
}

// try to guess system information from file name for a cubic lattice
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// nbrSiteZ = reference to the number sites along the z direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnCubicLatticeFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& nbrSiteZ, bool& statistics)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromFileName(filename, nbrParticles, nbrSiteX, nbrSiteY, statistics) == false)
    {
      return false;
    }
  char* StrNbrParticles = strstr(filename, "_z_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '.') || (StrNbrParticles[SizeString] == '_')) && (SizeString != 0))
	{
	  char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  nbrSiteZ = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      nbrSiteZ = 0;
    }
  return true;
}

// try to guess system information from file name for a cubic lattice
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// nbrSiteZ = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumY = reference to the momentum along the y direction
// momentumZ = reference to the momentum along the z direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnCubicLatticeFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& nbrSiteZ, int& momentumX, int& momentumY, int& momentumZ, bool& statistics)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSiteX, nbrSiteY, momentumX, momentumY, statistics) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  StrNbrParticles = strstr(filename, "_ky_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '.') || (StrNbrParticles[SizeString] == '_')) && (SizeString != 0))
	{
	  char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  momentumY = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      momentumY = 0;
    }
  StrNbrParticles = strstr(filename, "_kz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '.') || (StrNbrParticles[SizeString] == '_')) && (SizeString != 0))
	{
	  char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  momentumZ = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      momentumZ = 0;
    }
  StrNbrParticles = strstr(filename, "_z_");
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
	  nbrSiteZ = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      nbrSiteZ = 0;
    }
  return true;
}

// try to guess system information from file name for a Hofstadter lattice model
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// NbrCellX = number of magnetic unit cells along X
// NbrCellY = number of magnetic unit cells along Y
// Interaction = onsite interaction
// FluxPerCell = number of flux quanta per unit cell
// NbrState = number of the eigenstate
// Statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// Hardcore = flag indicating hard-core bosons
// EmbeddingFlag = flag indicating whether embedding is used
// Axis = character indicating axis of Landau gauge
// GammaX = periodic boundary conditions along X
// GammaY = periodic boundary conditions along Y
// MomentumX = momentum along x-direction
// MomentumY = momentum along y-direction
// UnitCellX = size of magnetic unit cell along x
// UnitCellY = size of magnetic unit cell along y
// 
bool FQHEOnSquareLatticeFindSystemInfoFromVectorFileName_Hofstadter(char* filename, int& NbrParticles, int& NbrCellX, int& NbrCellY, double& Interaction, int& FluxPerCell, int& NbrState, bool& Statistics, bool& Hardcore, bool& EmbeddingFlag, char& Axis, double& GammaX, double& GammaY, int& MomentumX, int& MomentumY, int& UnitCellX, int& UnitCellY)
{
  FilenameStatisticsCheck(Statistics, filename);
  FilenameIntegerSearch(UnitCellX, filename, "_X_");
  FilenameIntegerSearch(UnitCellY, filename, "_Y_");
  FilenameIntegerSearch(FluxPerCell, filename, "_q_");
  FilenameCharacterSearch(Axis, filename, "landau-");
  if (Axis != 'x' && Axis != 'y')
    {
      cout << "Did not find a valid Landau-gauge axis in file name - defaulting to landau-y" << endl;
      Axis = 'y';
    }
  FilenameIntegerSearch(NbrParticles, filename, "_n_");
  FilenameIntegerSearch(NbrCellX, filename, "_x_");
  FilenameIntegerSearch(NbrCellY, filename, "_y_");
  FilenameBooleanSearch(Hardcore, filename, "hardcore");
  if (Hardcore==false)
    FilenameDoubleSearch(Interaction, filename, "_u_");
  else
    Interaction=0;
  FilenameDoubleSearch(GammaX, filename, "_gx_");
  FilenameDoubleSearch(GammaY, filename, "_gy_");
  FilenameBooleanSearch(EmbeddingFlag, filename, "emb");
  FilenameIntegerSearch(MomentumX, filename, "_kx_");
  FilenameIntegerSearch(MomentumY, filename, "_ky_");
  FilenamePenultimateDotIntegerSearch(NbrState, filename);
  return true;
}


// try to guess system information from file name for a Hofstadter lattice model
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// NbrCellX = number of magnetic unit cells along X
// NbrCellY = number of magnetic unit cells along Y
// Interaction = onsite interaction
// FluxPerCell = number of flux quanta per unit cell
// NbrState = number of the eigenstate
// Statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// Hardcore = flag indicating hard-core bosons
// EmbeddingFlag = flag indicating whether embedding is used
// Axis = character indicating axis of Landau gauge
// GammaX = periodic boundary conditions along X
// GammaY = periodic boundary conditions along Y
// MomentumX = momentum along x-direction
// MomentumY = momentum along y-direction
// UnitCellX = size of magnetic unit cell along x
// UnitCellY = size of magnetic unit cell along y
// enlargeCell = true if unit cell contains 2 * FluxPerCell magnetic fluxes
// muS = amplitude of the symmetry-breaking on-site potential
// nbrBands = number of bands that have to be considered
// return value = true if no error occured
// 
bool FQHEOnSquareLatticeFindSystemInfoFromVectorFileName_Hofstadter(char* filename, int& NbrParticles, int& NbrCellX, int& NbrCellY, double& Interaction, int& FluxPerCell, int& NbrState, bool& Statistics, bool& Hardcore, bool& EmbeddingFlag, char& Axis, double& GammaX, double& GammaY, int& MomentumX, int& MomentumY, int& UnitCellX, int& UnitCellY, bool& enlargeCell, double& muS, int& nbrBands)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName_Hofstadter(filename, NbrParticles, NbrCellX, NbrCellY, Interaction, FluxPerCell, NbrState, Statistics, Hardcore, EmbeddingFlag, Axis, GammaX, GammaY, MomentumX, MomentumY, UnitCellX, UnitCellY) == false)
    return false;
    
  FilenameBooleanSearch (enlargeCell, filename, "_Xeff_");
  if (enlargeCell == false)
    return false;
  FilenameDoubleSearch (muS, filename, "_mus_");
  bool multipleBands = false;
  FilenameBooleanSearch (multipleBands, filename, "_b_0-");
  if (multipleBands)
    FilenameIntegerSearch (nbrBands, filename, "_b_0-");
  else
    nbrBands = 0;
  ++nbrBands;
  return true;  
}