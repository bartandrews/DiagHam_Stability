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


#include "Tools/SUNSpinFileTools.h"

#include <iostream>
#include <cstring>
#include <cstdlib>


using std::cout;
using std::endl;


// try to guess system information from file name
//
// filename = file name
// levelN = level of SU(N) symmetry
// descriptor = description of lattice
// nbrSpins = total number of spins
// cartanQuantumNumbers = quantum numbers of subspace, returned only if present in filename
// return value = true if no error occured
bool SUNSpinFindSystemInfoFromFileName(char* filename, int& levelN, char* &descriptor, int &nbrSpins, int* &cartanQuantumNumbers)
{
  char* StrNbrParticles;
  char* DescStart = NULL;
  char* DescEnd = NULL;
  if (levelN == 0)
    {
      StrNbrParticles = strstr(filename, "_SU");
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
	      levelN = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	      DescStart = StrNbrParticles+1;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess level of SU(N) from file name " << filename << endl;
	  return false;            
	}
    }
  if (nbrSpins == 0)
    {
      StrNbrParticles = strstr(filename, "_n_");
      if (StrNbrParticles != 0)
	{
	  DescEnd = StrNbrParticles-1;
	  StrNbrParticles += 3;
	  int SizeString = 0;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      nbrSpins = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess number of spins from file name " << filename << endl;
	  return false;            
	}
    }

  if (descriptor == NULL)
    {
      if (DescStart == NULL)
	DescStart = filename;
      if (DescEnd == NULL)
	DescEnd = filename+strlen(filename)-1;
      descriptor = new char[DescEnd-DescStart+2];
      strncpy(descriptor,DescStart,DescEnd-DescStart+1);
    }

  
  return true;
}


