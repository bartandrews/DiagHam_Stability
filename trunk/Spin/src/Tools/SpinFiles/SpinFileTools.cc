////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     set of functions used to managed files related to spin system          //
//                                                                            //
//                        last modification : 22/06/2010                      //
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


#include "Tools/SpinFiles/SpinFileTools.h"

#include <iostream>
#include <cstring>
#include <cstdlib>


using std::cout;
using std::endl;


// try to guess system information from file name
//
// filename = file name
// nbrSpins = reference to the number of spins
// return value = true if no error occured

bool SpinFindSystemInfoFromFileName(char* filename, int& nbrSpins)
{
  char* StrNbrSpins;
  StrNbrSpins = strstr(filename, "_n_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 3;
      int SizeString = 0;
      while ((StrNbrSpins[SizeString] != '\0') && (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if (SizeString != 0)
	{
	  char Tmp = StrNbrSpins[SizeString];
	  StrNbrSpins[SizeString] = '\0';
	  nbrSpins = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = Tmp;
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess number of spins from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name
//
// filename = file name
// nbrSpins = reference to the number of spins
// spin = reference to twice the spin value per site
// return value = true if no error occured

bool SpinFindSystemInfoFromFileName(char* filename, int& nbrSpins, int& spin)
{
  char* StrNbrSpins;
  StrNbrSpins = strstr(filename, "_n_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 3;
      int SizeString = 0;
      while ((StrNbrSpins[SizeString] != '\0') && (StrNbrSpins[SizeString] != '_') && (StrNbrSpins[SizeString] >= '0') 
	     && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrSpins[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrSpins[SizeString] = '\0';
	  nbrSpins = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = '_';
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess number of spins from file name " << filename << endl;
      return false;            
    }

  StrNbrSpins = strstr(filename, "spin_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 5;
      int SizeString = 0;
      while ((StrNbrSpins[SizeString] != '\0') && (StrNbrSpins[SizeString] != '_') && (StrNbrSpins[SizeString] >= '0') 
	     && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrSpins[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrSpins[SizeString] = '\0';
	  spin = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = '_';
	  StrNbrSpins += SizeString;
	  if (strstr(StrNbrSpins, "_2_") != StrNbrSpins)
	    {
	      spin *= 2;
	    }
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess spin value from file name " << filename << endl;
      return false;            
    }

  return true;
}

// try to guess system information from file name
//
// filename = file name
// nbrSpins = reference to the number of spins
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// return value = true if no error occured

bool SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int& spin)
{
  if (SpinFindSystemInfoFromFileName(filename, nbrSpins, spin) == false)
    return false;
  char* StrNbrSpins;

  StrNbrSpins = strstr(filename, "_sz_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 4;
      int SizeString = 0;
      if (StrNbrSpins[SizeString] == '-')
	++SizeString;
      while ((StrNbrSpins[SizeString] != '\0') && (StrNbrSpins[SizeString] != '.') && (StrNbrSpins[SizeString] != '_') && 
	     (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrSpins[SizeString] == '.') && (SizeString != 0))
	{
	  StrNbrSpins[SizeString] = '\0';
	  sz = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = '.';
	  StrNbrSpins += SizeString;
	}
      else
	{
	  if ((StrNbrSpins[SizeString] == '_') && (SizeString != 0))
	    {
	      StrNbrSpins[SizeString] = '\0';
	      sz = atoi(StrNbrSpins);
	      StrNbrSpins[SizeString] = '_';
	  StrNbrSpins += SizeString;
	    }
	  else
	    {
	      StrNbrSpins = 0;
	    }
	}
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess total Sz projection from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name
//
// filename = file name
// nbrSpins = reference to the number of spins
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// momentum = reference on the momentum
// return value = true if no error occured

bool SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int& spin, int& momentum)
{
  if (SpinFindSystemInfoFromVectorFileName(filename, nbrSpins, sz, spin) == false)
    return false;
  char* StrMomentum;
 
  StrMomentum = strstr(filename, "_k_");
  if (StrMomentum != 0)
    {
      StrMomentum += 3;
      int SizeString = 0;
      if (StrMomentum[SizeString] == '-')
	++SizeString;
      while ((StrMomentum[SizeString] != '\0') && (StrMomentum[SizeString] != '_') && (StrMomentum[SizeString] != '.') && (StrMomentum[SizeString] >= '0') 
	     && (StrMomentum[SizeString] <= '9'))
	++SizeString;
      if ((StrMomentum[SizeString] == '_') && (SizeString != 0))
	{
	  StrMomentum[SizeString] = '\0';
	  momentum = atoi(StrMomentum);
	  StrMomentum[SizeString] = '_';
	  StrMomentum += SizeString;
	}
      else
	{
	  if ((StrMomentum[SizeString] == '.') && (SizeString != 0))
	    {
	      StrMomentum[SizeString] = '\0';
	      momentum = atoi(StrMomentum);
	      StrMomentum[SizeString] = '.';
	      StrMomentum += SizeString;
	    }
	  else
	    {
	      StrMomentum = 0;
	    }
	}
    }
  if (StrMomentum == 0)
    {
      cout << "can't guess momentum from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name without a fixed total Sz value
//
// filename = file name
// nbrSpins = reference to the number of spins
// spin = reference to twice the spin value per site
// momentum = reference on the momentum
// return value = true if no error occured

bool SpinAllSzFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& spin, int& momentum)
{
  if (SpinFindSystemInfoFromFileName(filename, nbrSpins, spin) == false)
    return false;
  char* StrMomentum;
 
  StrMomentum = strstr(filename, "_k_");
  if (StrMomentum != 0)
    {
      StrMomentum += 3;
      int SizeString = 0;
      if (StrMomentum[SizeString] == '-')
	++SizeString;
      while ((StrMomentum[SizeString] != '\0') && (StrMomentum[SizeString] != '_') && (StrMomentum[SizeString] != '.') && (StrMomentum[SizeString] >= '0') 
	     && (StrMomentum[SizeString] <= '9'))
	++SizeString;
      if ((StrMomentum[SizeString] == '_') && (SizeString != 0))
	{
	  StrMomentum[SizeString] = '\0';
	  momentum = atoi(StrMomentum);
	  StrMomentum[SizeString] = '_';
	  StrMomentum += SizeString;
	}
      else
	{
	  if ((StrMomentum[SizeString] == '.') && (SizeString != 0))
	    {
	      StrMomentum[SizeString] = '\0';
	      momentum = atoi(StrMomentum);
	      StrMomentum[SizeString] = '.';
	      StrMomentum += SizeString;
	    }
	  else
	    {
	      StrMomentum = 0;
	    }
	}
    }
  if (StrMomentum == 0)
    {
      cout << "can't guess momentum from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name
//
// filename = file name
// nbrSpins = reference to the number of spins 
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// inversion =  reference on the inversion parity, will be non-zero only if the vector is encoded with the inversion symmetry
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// return value = true if no error occured

bool SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int& spin, int& inversion, int& szSymmetry)
{
  if (SpinFindSystemInfoFromVectorFileName(filename, nbrSpins, sz, spin) == false)
    return false;
  inversion = 0;
  szSymmetry = 0;
  char* StrNbrParticles = strstr(filename, "_invsym_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 8;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && 
	     (((StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] == '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  inversion = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	{
	  cout << "error while retrieving the inversion parity" << endl;
	  return false;
	}
    }
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

// try to guess system information from file name without a fixed total Sz value
//
// filename = file name
// nbrSpins = reference to the number of spins 
// spin = reference to twice the spin value per site
// momentum = reference on the momentum
// inversion =  reference on the inversion parity
// szSymmetry =  reference on the Sz<->-Sz parity
// return value = true if no error occured

bool SpinAllSzFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& spin, int& inversion, int& szSymmetry)
{
  if (SpinFindSystemInfoFromFileName(filename, nbrSpins, spin) == false)
    return false;
  inversion = 0;
  szSymmetry = 0;
  char* StrNbrParticles = strstr(filename, "_invsym_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 8;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && 
	     (((StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] == '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  inversion = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	{
	  cout << "error while retrieving the inversion parity" << endl;
	  return false;
	}
    }
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
// filename = file name
// nbrSpins = reference to the number of spins 
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// momentum = reference on the momentum
// inversion =  reference on the inversion parity, will be non-zero only if the vector is encoded with the inversion symmetry
// szSymmetry =  reference on the Sz<->-Sz parity, will be non-zero only if the vector is encoded with the Sz<->-Sz symmetry
// return value = true if no error occured

bool SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int& spin, int& momentum, int& inversion, int& szSymmetry)
{
  if (SpinFindSystemInfoFromVectorFileName(filename, nbrSpins, sz, spin, momentum) == false)
    return false;
  inversion = 0;
  szSymmetry = 0;
  char* StrNbrParticles = strstr(filename, "_invsym_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 8;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && 
	     (((StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] == '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  inversion = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	{
	  cout << "error while retrieving the inversion parity" << endl;
	  return false;
	}
    }
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

// try to guess system information from file name without a fixed total Sz value
//
// filename = file name
// nbrSpins = reference to the number of spins 
// spin = reference to twice the spin value per site
// momentum = reference on the momentum
// inversion =  reference on the inversion parity
// szSymmetry =  reference on the Sz<->-Sz parity
// return value = true if no error occured

bool SpinAllSzFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& spin, int& momentum, int& inversion, int& szSymmetry)
{
  if (SpinAllSzFindSystemInfoFromVectorFileName(filename, nbrSpins, spin, momentum) == false)
    return false;
  inversion = 0;
  szSymmetry = 0;
  char* StrNbrParticles = strstr(filename, "_invsym_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 8;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && 
	     (((StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] == '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  inversion = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	{
	  cout << "error while retrieving the inversion parity" << endl;
	  return false;
	}
    }
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

// try to guess system information from file name for a 2d spin system with translations
//
// filename = file name
// nbrSpins = reference to the number of spins 
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// xMomentum = reference on the momentum along the x direction
// xPeriodicity = reference on the number of sites along the x direction
// yMomentum = reference on the momentum along the y direction
// yPeriodicity = reference on the number of sites along the y direction
// return value = true if no error occured

bool SpinWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int& spin, int& xMomentum, int& xPeriodicity,
							   int& yMomentum, int& yPeriodicity)
{
  if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(filename, nbrSpins, spin, xMomentum, xPeriodicity, yMomentum, yPeriodicity) == false)
    return false;
  char* StrNbrSpins;

  StrNbrSpins = strstr(filename, "_sz_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 4;
      int SizeString = 0;
      if (StrNbrSpins[SizeString] == '-')
	++SizeString;
      while ((StrNbrSpins[SizeString] != '\0') && (StrNbrSpins[SizeString] != '.') && 
	     (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9') && (StrNbrSpins[SizeString] != '_'))
	++SizeString;
      if ((StrNbrSpins[SizeString] == '.') && (SizeString != 0))
	{
	  StrNbrSpins[SizeString] = '\0';
	  sz = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = '.';
	  StrNbrSpins += SizeString;
	}
      else
	{
	  if ((StrNbrSpins[SizeString] == '_') && (SizeString != 0))
	    {
	      StrNbrSpins[SizeString] = '\0';
	      sz = atoi(StrNbrSpins);
	      StrNbrSpins[SizeString] = '_';
	  StrNbrSpins += SizeString;
	    }
	  else
	    {
	      StrNbrSpins = 0;
	    }
	}
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess total Sz projection from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name for a 2d spin system with translations
//
// filename = file name
// nbrSpins = reference to the number of spins 
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// xMomentum = reference on the momentum along the x direction
// xPeriodicity = reference on the number of sites along the x direction
// yMomentum = reference on the momentum along the y direction
// yPeriodicity = reference on the number of sites along the y direction
// inversion =  reference on the inversion parity
// return value = true if no error occured

bool SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int& spin, int& xMomentum, int& xPeriodicity,
							   int& yMomentum, int& yPeriodicity, int& inversion)
{
  if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(filename, nbrSpins, sz, spin, xMomentum, xPeriodicity, yMomentum, yPeriodicity) == false)
    return false;
  char* StrNbrParticles = strstr(filename, "_i_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && 
	     (((StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] == '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  inversion = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess inversion sector from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name for a 2d spin system with translations
//
// filename = file name
// nbrSpins = reference to the number of spins 
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// xMomentum = reference on the momentum along the x direction
// xPeriodicity = reference on the number of sites along the x direction
// yMomentum = reference on the momentum along the y direction
// yPeriodicity = reference on the number of sites along the y direction
// inversion =  reference on the inversion parity
// szSymmetry =  reference on the Sz<->-Sz parity
// return value = true if no error occured

bool SpinWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int& spin, int& xMomentum, int& xPeriodicity,
							   int& yMomentum, int& yPeriodicity, int& inversion, int& szSymmetry)
{
  if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(filename, nbrSpins, sz, spin, xMomentum, xPeriodicity,
							    yMomentum, yPeriodicity) == false)
    {
      return false;
    }
  inversion = 0;
  szSymmetry = 0;
  char* StrNbrParticles = strstr(filename, "_invsym_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 8;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && 
	     (((StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] == '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  inversion = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	{
	  cout << "error while retrieving the inversion parity" << endl;
	  return false;
	}
    }
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

// try to guess system information from file name for a 2d spin system with translations
//
// filename = file name
// nbrSpins = reference to the number of spins 
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// xMomentum = reference on the momentum along the x direction
// xPeriodicity = reference on the number of sites along the x direction
// yMomentum = reference on the momentum along the y direction
// yPeriodicity = reference on the number of sites along the y direction
// return value = true if no error occured

bool SpinWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& spin, int& xMomentum, int& xPeriodicity,
							   int& yMomentum, int& yPeriodicity)
{
  if (SpinFindSystemInfoFromFileName(filename, nbrSpins, spin) == false)
    return false;
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
      cout << "can't guess x-momentum sector from file name " << filename << endl;
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
      cout << "can't guess x-periodicity from file name " << filename << endl;
      return false;            
    }
  StrNbrParticles = strstr(filename, "_ky_");
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
	  yMomentum = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess y-momentum sector from file name " << filename << endl;
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
      cout << "can't guess y-periodicity from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name for a 2d spin system with translations
//
// filename = file name
// nbrSpins = reference to the number of spins 
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// xMomentum = reference on the momentum along the x direction
// xPeriodicity = reference on the number of sites along the x direction
// yMomentum = reference on the momentum along the y direction
// yPeriodicity = reference on the number of sites along the y direction
// return value = true if no error occured

bool SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int& spin, int& momentum, int& inversion, int& szSymmetry, int& offset)
{
  SpinFindSystemInfoFromVectorFileName (filename, nbrSpins, sz, spin, momentum, inversion, szSymmetry);
  char* StrNbrParticles;
  StrNbrParticles = strstr(filename, "_off_");
  if (StrNbrParticles != 0)
    {
      cout << StrNbrParticles << endl;
      StrNbrParticles += 5;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  offset = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  else
    offset = 0;
  return true;
}

// try to guess system information from file name for a 2d spin system with translations
//
// filename = file name
// nbrSpins = reference to the number of spins 
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// xMomentum = reference on the momentum along the x direction
// xPeriodicity = reference on the number of sites along the x direction
// yMomentum = reference on the momentum along the y direction
// yPeriodicity = reference on the number of sites along the y direction
// inversion =  reference on the inversion parity
// return value = true if no error occured
bool SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& spin, int& xMomentum, int& xPeriodicity,
							   int& yMomentum, int& yPeriodicity, int& inversion)
{
  if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(filename, nbrSpins, spin, xMomentum, xPeriodicity, yMomentum, yPeriodicity) == false)
    return false;
  char* StrNbrParticles = strstr(filename, "_i_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && 
	     (((StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] == '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  inversion = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess inversion sector from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name
//
// filename = file name
// nbrSites = reference to the number of sites
// qValue = reference to Zn charge
// return value = true if no error occured

bool PottsFindSystemInfoFromVectorFileName(char* filename, int& nbrSites, int& qValue)
{
  if (SpinFindSystemInfoFromFileName(filename, nbrSites) == false)
    return false;
  char* StrNbrSpins;
  StrNbrSpins = strstr(filename, "_q_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 3;
      int SizeString = 0;
      if (StrNbrSpins[SizeString] == '-')
	++SizeString;
      while ((StrNbrSpins[SizeString] != '\0') && 
	     (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((SizeString > 1) || ((SizeString == 1) && (StrNbrSpins[0] != '-')))
	{
	  char TmpChar = StrNbrSpins[SizeString];
	  StrNbrSpins[SizeString] = '\0';
	  qValue = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = TmpChar;
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess the total Zn charge from file name " << filename << endl;
      return false;            
    }
  return true; 
}

// try to guess system information from file name without a fixed total Q value
//
// filename = file name
// nbrSpins = reference to the number of spins 
// nValue = reference to the n of the Zn
// qValue = reference to Zn charge
// return value = true if no error occured

bool PottsFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& nValue, int& qValue)
{
  if (PottsFindSystemInfoFromVectorFileName(filename, nbrSpins, qValue) == false)
    return false;
  char* StrNbrSpins;
  StrNbrSpins = strstr(filename, "potts");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 5;
      if ((StrNbrSpins[0] >= '0') && (StrNbrSpins[0] <= '9'))
	{
	  char TmpChar = StrNbrSpins[1];
	  StrNbrSpins[1] = '\0';
	  nValue = atoi(StrNbrSpins);
	  StrNbrSpins[1] = TmpChar;
	}
      else
	{
	  StrNbrSpins = 0;
	}
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess Zn from file name " << filename << endl;
      return false;
    }
  return true;  
}

// try to guess system information from file name without a fixed total Q value
//
// filename = file name
// nbrSpins = reference to the number of spins 
// nValue = reference to the n of the Zn
// qValue = reference to Zn charge
// momentum = reference on the momentum
// return value = true if no error occured

bool PottsFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& nValue, int& qValue, int& momentum)
{
  if (PottsFindSystemInfoFromVectorFileName(filename, nbrSpins, nValue, qValue) == false)
    return false;
  char* StrMomentum = strstr(filename, "_k_");
  if (StrMomentum != 0)
    {
      StrMomentum += 3;
      int SizeString = 0;
      if (StrMomentum[SizeString] == '-')
	++SizeString;
      while ((StrMomentum[SizeString] != '\0') && (StrMomentum[SizeString] != '_') && (StrMomentum[SizeString] != '.') && (StrMomentum[SizeString] >= '0') 
	     && (StrMomentum[SizeString] <= '9'))
	++SizeString;
      if ((StrMomentum[SizeString] == '_') && (SizeString != 0))
	{
	  StrMomentum[SizeString] = '\0';
	  momentum = atoi(StrMomentum);
	  StrMomentum[SizeString] = '_';
	  StrMomentum += SizeString;
	}
      else
	{
	  if ((StrMomentum[SizeString] == '.') && (SizeString != 0))
	    {
	      StrMomentum[SizeString] = '\0';
	      momentum = atoi(StrMomentum);
	      StrMomentum[SizeString] = '.';
	      StrMomentum += SizeString;
	    }
	  else
	    {
	      StrMomentum = 0;
	    }
	}
    }
  if (StrMomentum == 0)
    {
      cout << "can't guess momentum from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name without a fixed total Q value
//
// filename = file name
// nbrSpins = reference to the number of spins 
// nValue = reference to the n of the Zn
// qValue = reference to Zn charge
// momentum = reference on the momentum
// inversion =  reference on the inversion parity
// return value = true if no error occured

bool PottsFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& nValue, int& qValue, int& momentum, int& inversion)
{
  if (PottsFindSystemInfoFromVectorFileName(filename, nbrSpins, nValue, qValue, momentum) == false)
    return false;
  inversion = 0;
  char* StrNbrParticles = strstr(filename, "_invsym_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 8;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && 
	     (((StrNbrParticles[SizeString] >= '0') && (StrNbrParticles[SizeString] <= '9')) || (StrNbrParticles[SizeString] == '-')))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  inversion = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	{
	  cout << "error while retrieving the inversion parity" << endl;
	  return false;
	}
    }
  return true;
}

bool PEPSFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz)
{
  char* StrNbrSpins;
  StrNbrSpins = strstr(filename, "_l_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 3;
      int SizeString = 0;
      while ((StrNbrSpins[SizeString] != '\0') && (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if (SizeString != 0)
	{
	  char Tmp = StrNbrSpins[SizeString];
	  StrNbrSpins[SizeString] = '\0';
	  nbrSpins = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = Tmp;
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess number of spins from file name " << filename << endl;
      return false;            
    }
  
  
  StrNbrSpins = strstr(filename, "_sz_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 4;
      int SizeString = 0;
      if (StrNbrSpins[SizeString] == '-')
	++SizeString;
      while ((StrNbrSpins[SizeString] != '\0') && (StrNbrSpins[SizeString] != '.') && 
	     (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((SizeString != 0))
	{
	  char Tmp = StrNbrSpins[SizeString];
	  StrNbrSpins[SizeString] = '\0';
	  sz = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = Tmp;
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess total Sz projection from file name " << filename << endl;
      return false;            
    }
  return true;
}


bool PEPSFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int & momentum)
{
  if (PEPSFindSystemInfoFromVectorFileName(filename, nbrSpins, sz) == false)
    return false;

  char* StrMomentum;
 
  StrMomentum = strstr(filename, "_k_");
  if (StrMomentum != 0)
    {
      StrMomentum += 3;
      int SizeString = 0;
      if (StrMomentum[SizeString] == '-')
	++SizeString;
      while ((StrMomentum[SizeString] != '\0') && (StrMomentum[SizeString] != '_') && (StrMomentum[SizeString] != '.') && (StrMomentum[SizeString] >= '0') 
	     && (StrMomentum[SizeString] <= '9'))
	++SizeString;
      if ((StrMomentum[SizeString] == '_') && (SizeString != 0))
	{
	  StrMomentum[SizeString] = '\0';
	  momentum = atoi(StrMomentum);
	  StrMomentum[SizeString] = '_';
	  StrMomentum += SizeString;
	}
      else
	{
	  if ((StrMomentum[SizeString] == '.') && (SizeString != 0))
	    {
	      StrMomentum[SizeString] = '\0';
	      momentum = atoi(StrMomentum);
	      StrMomentum[SizeString] = '.';
	      StrMomentum += SizeString;
	    }
	  else
	    {
	      StrMomentum = 0;
	    }
	}
    }
  if (StrMomentum == 0)
    {
      cout << "can't guess momentum from file name " << filename << endl;
      return false;            
    }
  return true;
}



bool PEPSFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int & valueOfZBra, int & valueOfZKet)
{
  if (PEPSFindSystemInfoFromVectorFileName(filename, nbrSpins, sz) == false)
    return false;
  
  char* StrNbrSpins;
  StrNbrSpins = strstr(filename, "_zbra_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 6;
      int SizeString = 0;
      if (StrNbrSpins[SizeString] == '-')
	++SizeString;
      while ((StrNbrSpins[SizeString] != '\0') && 
	     (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((SizeString > 1) || ((SizeString == 1) && (StrNbrSpins[0] != '-')))
	{
	  StrNbrSpins[SizeString] = '\0';
	  valueOfZBra = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = '_';
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess the value of Z symmetry bra  from file name " << filename << endl;
      return false;            
    }

  StrNbrSpins = strstr(filename, "_zket_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 6;
      int SizeString = 0;
      if (StrNbrSpins[SizeString] == '-')
	++SizeString;
      while ((StrNbrSpins[SizeString] != '\0') && 
	     (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((SizeString > 1) || ((SizeString == 1) && (StrNbrSpins[0] != '-')))
	{
	  char tmp = StrNbrSpins[SizeString];
	  StrNbrSpins[SizeString] = '\0';
	  valueOfZKet = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = tmp;
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess the value of Z symmetry ket  from file name " << filename << endl;
      return false;            
    }
  
  return true;
}

bool PEPSFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int & momentum, int & valueOfZBra, int & valueOfZKet)
{
  if (PEPSFindSystemInfoFromVectorFileName(filename, nbrSpins, sz, momentum) == false)
    return false;

  char* StrNbrSpins;
  StrNbrSpins = strstr(filename, "_zbra_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 6;
      int SizeString = 0;
      if (StrNbrSpins[SizeString] == '-')
	++SizeString;
      while ((StrNbrSpins[SizeString] != '\0') && 
	     (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((SizeString > 1) || ((SizeString == 1) && (StrNbrSpins[0] != '-')))
	{
	  StrNbrSpins[SizeString] = '\0';
	  valueOfZBra = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = '_';
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess the value of Z symmetry bra  from file name " << filename << endl;
      return false;            
    }

  StrNbrSpins = strstr(filename, "_zket_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 6;
      int SizeString = 0;
      if (StrNbrSpins[SizeString] == '-')
	++SizeString;
      while ((StrNbrSpins[SizeString] != '\0') && 
	     (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((SizeString > 1) || ((SizeString == 1) && (StrNbrSpins[0] != '-')))
	{
	  char tmp = StrNbrSpins[SizeString];
	  StrNbrSpins[SizeString] = '\0';
	  valueOfZKet = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = tmp;
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess the value of Z symmetry ket  from file name " << filename << endl;
      return false;            
    }
  
  return true;
}



bool PEPSFindSubLatticeNumbersFromVectorFileName(char* filename, int& nbrSpins, int& sz, int & momentum, int & valueSubLatticeZeroBra, int & valueSubLatticeZeroKet )
{
  if (PEPSFindSystemInfoFromVectorFileName(filename, nbrSpins, sz, momentum) == false)
   return false;
   
  char* StrNbrSpins;
  StrNbrSpins = strstr(filename, "_sbra_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 6;
      int SizeString = 0;
      if (StrNbrSpins[SizeString] == '-')
	++SizeString;
      while ((StrNbrSpins[SizeString] != '\0') && 
	     (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((SizeString > 1) || ((SizeString == 1) && (StrNbrSpins[0] != '-')))
	{
	  StrNbrSpins[SizeString] = '\0';
	  valueSubLatticeZeroBra = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = '_';
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess the value of S symmetry bra  from file name " << filename << endl;
      return false;            
    }

  StrNbrSpins = strstr(filename, "_sket_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 6;
      int SizeString = 0;
      if (StrNbrSpins[SizeString] == '-')
	++SizeString;
      while ((StrNbrSpins[SizeString] != '\0') && 
	     (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((SizeString > 1) || ((SizeString == 1) && (StrNbrSpins[0] != '-')))
	{
	  char tmp = StrNbrSpins[SizeString];
	  StrNbrSpins[SizeString] = '\0';
	  valueSubLatticeZeroKet = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = tmp;
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess the value of S symmetry ket  from file name " << filename << endl;
      return false;            
    }
  return true;
}

bool PEPSFindSubLatticeNumbersFromVectorFileName(char* filename, int& nbrSpins, int& sz, int & momentum, int & valueSubLatticeZeroBra, int & valueSubLatticeZeroKet, int & valueSubLatticeZeroProduct) 
{
  if (PEPSFindSubLatticeNumbersFromVectorFileName(filename, nbrSpins, sz, momentum, valueSubLatticeZeroBra, valueSubLatticeZeroKet) == false)
    return false;
  
  char* StrNbrSpins;  
  StrNbrSpins = strstr(filename, "_sprod_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 7;
      int SizeString = 0;
      if (StrNbrSpins[SizeString] == '-')
	++SizeString;
      while ((StrNbrSpins[SizeString] != '\0') && 
	     (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((SizeString > 1) || ((SizeString == 1) && (StrNbrSpins[0] != '-')))
	{
	  char tmp = StrNbrSpins[SizeString];
	  StrNbrSpins[SizeString] = '\0';
	  valueSubLatticeZeroProduct = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = tmp;
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess the value of the product of S symmetry from file name " << filename << endl;
      return false;            
    }
 
  return true;
}


bool PEPSFindSystemInfoFromVectorFileNameUndescribedSpace(char* filename, int& nbrSpins, int& bondDimension)
{
  char* StrNbrSpins;
  StrNbrSpins = strstr(filename, "_l_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 3;
      int SizeString = 0;
      while ((StrNbrSpins[SizeString] != '\0') && (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if (SizeString != 0)
	{
	  char Tmp = StrNbrSpins[SizeString];
	  StrNbrSpins[SizeString] = '\0';
	  nbrSpins = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = Tmp;
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess number of spins from file name " << filename << endl;
      return false;            
    }
  
  
  StrNbrSpins = strstr(filename, "_d_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 3;
      int SizeString = 0;
      if (StrNbrSpins[SizeString] == '-')
	++SizeString;
      while ((StrNbrSpins[SizeString] != '\0') && (StrNbrSpins[SizeString] != '.') && 
	     (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((SizeString != 0))
	{
	  char Tmp = StrNbrSpins[SizeString];
	  StrNbrSpins[SizeString] = '\0';
	  bondDimension = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = Tmp;
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess the bond dimension  file name " << filename << endl;
      return false;            
    }
  return true;
}


bool PEPSFindSystemInfoFromVectorFileNameUndescribedSpace(char* filename, int& nbrSpins, int& bondDimension, int & momentum)
{
  if (PEPSFindSystemInfoFromVectorFileNameUndescribedSpace(filename, nbrSpins, bondDimension) == false)
    return false;

  char* StrMomentum;
 
  StrMomentum = strstr(filename, "_k_");
  if (StrMomentum != 0)
    {
      StrMomentum += 3;
      int SizeString = 0;
      if (StrMomentum[SizeString] == '-')
	++SizeString;
      while ((StrMomentum[SizeString] != '\0') && (StrMomentum[SizeString] != '_') && (StrMomentum[SizeString] != '.') && (StrMomentum[SizeString] >= '0') 
	     && (StrMomentum[SizeString] <= '9'))
	++SizeString;
      if ((StrMomentum[SizeString] == '_') && (SizeString != 0))
	{
	  StrMomentum[SizeString] = '\0';
	  momentum = atoi(StrMomentum);
	  StrMomentum[SizeString] = '_';
	  StrMomentum += SizeString;
	}
      else
	{
	  if ((StrMomentum[SizeString] == '.') && (SizeString != 0))
	    {
	      StrMomentum[SizeString] = '\0';
	      momentum = atoi(StrMomentum);
	      StrMomentum[SizeString] = '.';
	      StrMomentum += SizeString;
	    }
	  else
	    {
	      StrMomentum = 0;
	    }
	}
    }
  if (StrMomentum == 0)
    {
      cout << "can't guess momentum from file name " << filename << endl;
      return false;            
    }
  return true;
}

