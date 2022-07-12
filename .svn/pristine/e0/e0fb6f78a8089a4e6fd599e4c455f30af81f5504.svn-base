////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          set of functions used to handle torus pseudopotential files       //
//                                                                            //
//                        last modification : 17/04/2012                      //
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
#include "Tools/FQHEFiles/FQHETorusPseudopotentialTools.h"
#include "GeneralTools/ConfigurationParser.h"


#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cstdio>

using std::cout;
using std::endl;
using std::string;


// get pseudopototentials for particles on torus from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrPseudoPotentials = reference on the number of pseudopotentials
// pseudoPotentials = reference on the array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
// return value = true if no error occured

bool FQHETorusGetPseudopotentials (char* fileName, int& nbrPseudoPotentials, double*& pseudoPotentials)
{
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }
  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', pseudoPotentials, nbrPseudoPotentials) == false)
    {
      cout << "Pseudopotentials has a wrong value in " << fileName << endl;
      return false;
    }
  return true;
}


// get pseudopototentials for particles on torus from file with an interaction name
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrPseudoPotentials = reference on the number of pseudopotentials
// pseudoPotentials = reference on the array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
// interactionName = naming convention read from definition, or generated from LL index
// return value = true if no error occured

bool FQHETorusGetPseudopotentials (char* fileName, int& nbrPseudoPotentials, double*& pseudoPotentials, char*& interactionName)
{
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }
  if (InteractionDefinition["Name"] == NULL)
    {
      cout << "Attention, using unnamed interaction! Please include a line 'Name = ...'" << endl;
      interactionName = new char[10];
      sprintf(interactionName,"unnamed");
    }
  else
    {
      interactionName = new char[strlen(InteractionDefinition["Name"])+1];
      strcpy(interactionName, InteractionDefinition["Name"]);
    }
  InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', pseudoPotentials, nbrPseudoPotentials);
  
  return true;
}


// get pseudopototentials for particles on torus from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrFluxQuanta = number of flux quanta
// nbrPseudoPotentials = reference on the number of pseudopotentials
// pseudoPotentials = reference on the array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
// oneBodyPseudoPotentials  = array with the one-body pseudo-potentials
// return value = true if no error occured

bool FQHETorusGetPseudopotentials (char* fileName, int nbrFluxQuanta, int& nbrPseudoPotentials, double*& pseudoPotentials, double*& oneBodyPotentials)
{
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }
  if (InteractionDefinition["Pseudopotentials"] != NULL)
    {
      if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', pseudoPotentials, nbrPseudoPotentials) == false)
	{
	  cout << "Pseudopotentials has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  else
    {
      pseudoPotentials = 0;
      nbrPseudoPotentials = 0;
    }
  int TmpNbrPseudoPotentials = 0;
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentials", ' ', oneBodyPotentials, TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != nbrFluxQuanta)
	{
	  cout << "OneBodyPotentials has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  return true;
}

// get pseudopototentials for particles on torus from file
// 
// fileName = name of the file that contains the pseudopotantial description
// landauLevel = index of Coulomb Landau-level
// nbrPseudoPotentials = reference on the number of pseudopotentials
// pseudoPotentials = reference on the array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
// interactionName = naming convention read from definition, or generated from LL index
// return value = true if no error occured

bool FQHETorusGetPseudopotentials (char* fileName, bool haveCoulomb, int &landauLevel, int& nbrPseudoPotentials, double*& pseudoPotentials, char*& interactionName)
{
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }

  if (InteractionDefinition["CoulombLandauLevel"] != NULL)
    {
      landauLevel = atoi(InteractionDefinition["CoulombLandauLevel"]);
      haveCoulomb=true;
    }
  else
    {
      haveCoulomb=false;
    }
  if (InteractionDefinition["Name"] == NULL)
    {
      if ((InteractionDefinition["CoulombLandauLevel"] != NULL) && (InteractionDefinition["Pseudopotentials"] == NULL))
	{
	  interactionName = new char[18];
	  if (landauLevel>=0)
	    sprintf(interactionName,"coulomb_l_%d",landauLevel);
	  else
	    sprintf(interactionName,"graphene_l_%d",-landauLevel);
	}
      else
	{
	  cout << "Attention, using unnamed interaction! Please include a line 'Name = ...'" << endl;
	  interactionName = new char[10];
	  sprintf(interactionName,"unnamed");
	}
    }
  else
    {
      interactionName = new char[strlen(InteractionDefinition["Name"])+1];
      strcpy(interactionName, InteractionDefinition["Name"]);
    }
  InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', pseudoPotentials, nbrPseudoPotentials);
  
  return true;
}


// get pseudopototentials for particles on torus with SU(2) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// return value = true if no error occured

bool FQHETorusSU2GetPseudopotentials (char* fileName, int* nbrPseudoPotentials, double** pseudoPotentials)
{
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }
  int TmpNbrPseudoPotentials;
  double* TmpPseudoPotentials;
  bool Flag = false;
  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      for (int i = 0; i < 3; ++i)
	{
	  nbrPseudoPotentials[i] = TmpNbrPseudoPotentials;
	  pseudoPotentials[i] = new double[nbrPseudoPotentials[i]];
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    pseudoPotentials[i][j] = TmpPseudoPotentials[j];
	}
    }
  else
    if (InteractionDefinition["Pseudopotentials"] != 0)
      {
	cout << "Pseudopotentials has a wrong value in " << fileName << endl;
	return false;
      }
  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpUp", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[0] = TmpNbrPseudoPotentials;
      pseudoPotentials[0] = new double[nbrPseudoPotentials[0]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[0][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsUpUp"] != 0)
      {
	cout << "PseudopotentialsUpUp has a wrong value in " << fileName << endl;
	return false;
      }
  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[1] = TmpNbrPseudoPotentials;
      pseudoPotentials[1] = new double[nbrPseudoPotentials[1]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[1][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsDownDown"] != 0)
      {
	cout << "PseudopotentialsDownDown has a wrong value in " << fileName << endl;
	return false;
      }
  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[2] = TmpNbrPseudoPotentials;
      pseudoPotentials[2] = new double[nbrPseudoPotentials[2]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[2][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsUpDown"] != 0)
      {
	cout << "PseudopotentialsUpDown has a wrong value in " << fileName << endl;
	return false;
      }
  return true;
}

// get pseudopototentials for particles on torus with SU(2) spin from file, including the one-body pseudopotentials if any
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrFluxQuanta = number of flux quanta
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// oneBobyPseudoPotentials  = array with the one body pseudo-potentials (sorted from the ky=0 to the ky=nphi-1)
//                            first index refered to the spin sector (sorted as up-up, down-down, up-down)
// return value = true if no error occured

bool FQHETorusSU2GetPseudopotentials (char* fileName, int nbrFluxQuanta, int* nbrPseudoPotentials, double** pseudoPotentials, double** oneBobyPseudoPotentials)
{
  if (FQHETorusSU2GetPseudopotentials (fileName, nbrPseudoPotentials, pseudoPotentials) == false)
    {
      return false;
    }
  return FQHETorusSU2GetOneBodyPseudopotentials (fileName, nbrFluxQuanta, oneBobyPseudoPotentials[0], oneBobyPseudoPotentials[1], oneBobyPseudoPotentials[2]);
}

// get pseudopototentials for particles on torus with SU(2) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrFluxQuanta = number of flux quanta
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// return value = true if no error occured

bool FQHETorusSU2GetOneBodyPseudopotentials (char* fileName, int nbrFluxQuanta, double*& oneBodyPotentialUpUp, double*& oneBodyPotentialDownDown, double*& oneBodyPotentialUpDown)
{
  int TmpNbrPseudoPotentials;
  ConfigurationParser InteractionDefinition;
  oneBodyPotentialUpUp = 0;
  oneBodyPotentialDownDown = 0;
  oneBodyPotentialUpDown = 0;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialUpUp", ' ', oneBodyPotentialUpUp, TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != nbrFluxQuanta)
	{
	  cout << "OneBodyPotentialUpUp has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialDownDown", ' ', oneBodyPotentialDownDown, TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != nbrFluxQuanta)
	{
	  cout << "OneBodyPotentialDownDown has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialUpDown", ' ', oneBodyPotentialUpDown, TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != nbrFluxQuanta)
	{
	  cout << "OneBodyPotentialUpDown has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  return true;
}

// get pseudopototentials for particles on torus with SU(3) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as 11, 12, 13, 22, 23, 33)
// return value = true if no error occured

bool FQHETorusSU3GetPseudopotentials (char* fileName, int* nbrPseudoPotentials, double** pseudoPotentials)
{
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }
  int TmpNbrPseudoPotentials;
  double* TmpPseudoPotentials;
  bool Flag = false;
  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      for (int i = 0; i < 6; ++i)
	{
	  nbrPseudoPotentials[i] = TmpNbrPseudoPotentials;
	  pseudoPotentials[i] = new double[nbrPseudoPotentials[i]];
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    pseudoPotentials[i][j] = TmpPseudoPotentials[j];
	}
    }
  else
    if (InteractionDefinition["Pseudopotentials"] != 0)
      {
	cout << "Pseudopotentials has a wrong value in " << fileName << endl;
	return false;
      }
  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials11", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[0] = TmpNbrPseudoPotentials;
      pseudoPotentials[0] = new double[nbrPseudoPotentials[0]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[0][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["Pseudopotentials11"] != 0)
      {
	cout << "Pseudopotentials11 has a wrong value in " << fileName << endl;
	return false;
      }
  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials12", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[1] = TmpNbrPseudoPotentials;
      pseudoPotentials[1] = new double[nbrPseudoPotentials[1]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[1][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["Pseudopotentials12"] != 0)
      {
	cout << "Pseudopotentials12 has a wrong value in " << fileName << endl;
	return false;
      }
  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials13", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[2] = TmpNbrPseudoPotentials;
      pseudoPotentials[2] = new double[nbrPseudoPotentials[2]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[2][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["Pseudopotentials13"] != 0)
      {
	cout << "Pseudopotentials13 has a wrong value in " << fileName << endl;
	return false;
      }
  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials22", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[3] = TmpNbrPseudoPotentials;
      pseudoPotentials[3] = new double[nbrPseudoPotentials[3]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[3][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["Pseudopotentials22"] != 0)
      {
	cout << "Pseudopotentials22 has a wrong value in " << fileName << endl;
	return false;
      }
  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials23", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[4] = TmpNbrPseudoPotentials;
      pseudoPotentials[4] = new double[nbrPseudoPotentials[4]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[4][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["Pseudopotentials23"] != 0)
      {
	cout << "Pseudopotentials23 has a wrong value in " << fileName << endl;
	return false;
      }
  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials33", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[5] = TmpNbrPseudoPotentials;
      pseudoPotentials[5] = new double[nbrPseudoPotentials[5]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[5][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["Pseudopotentials33"] != 0)
      {
	cout << "Pseudopotentials33 has a wrong value in " << fileName << endl;
	return false;
      }
  return true;
}

// get pseudopototentials for particles on torus with SU(3) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrFluxQuanta = number of flux quanta
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as 11, 12, 13, 22, 23, 33)
// oneBodyPseudoPotentials  = array with the one body pseudo-potentials (sorted from the ky=0 to the ky=nphi-1)
//                            first index refered to the spin sector (sorted as 11, 12, 13, 22, 23, 33)
// return value = true if no error occured

bool FQHETorusSU3GetPseudopotentials (char* fileName, int nbrFluxQuanta, int* nbrPseudoPotentials, double** pseudoPotentials, double** oneBodyPseudoPotentials)
{
  if (FQHETorusSU3GetPseudopotentials(fileName, nbrPseudoPotentials, pseudoPotentials) == false)
    return  false;
  return FQHETorusSU3GetOneBodyPseudopotentials (fileName, nbrFluxQuanta, oneBodyPseudoPotentials);
}
// get pseudopototentials for particles on torus with SU(3) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrFluxQuanta = number of flux quanta
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// oneBodyPseudoPotentials  = array with the one body pseudo-potentials (sorted from the ky=0 to the ky=nphi-1)
//                            first index refered to the spin sector (sorted as 11, 12, 13, 22, 23, 33)
// return value = true if no error occured

bool FQHETorusSU3GetOneBodyPseudopotentials (char* fileName, int nbrFluxQuanta, double** oneBodyPseudoPotentials)
{
  int TmpNbrPseudoPotentials;
  ConfigurationParser InteractionDefinition;
  for (int i = 0; i < 6; ++i)
    oneBodyPseudoPotentials[i] = 0;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotential11", ' ', oneBodyPseudoPotentials[0], TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != nbrFluxQuanta)
	{
	  cout << "OneBodyPotential11 has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotential12", ' ', oneBodyPseudoPotentials[1], TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != nbrFluxQuanta)
	{
	  cout << "OneBodyPotential12 has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotential13", ' ', oneBodyPseudoPotentials[2], TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != nbrFluxQuanta)
	{
	  cout << "OneBodyPotential13 has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotential22", ' ', oneBodyPseudoPotentials[3], TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != nbrFluxQuanta)
	{
	  cout << "OneBodyPotential22 has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotential23", ' ', oneBodyPseudoPotentials[4], TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != nbrFluxQuanta)
	{
	  cout << "OneBodyPotential23 has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotential33", ' ', oneBodyPseudoPotentials[5], TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != nbrFluxQuanta)
	{
	  cout << "OneBodyPotential33 has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  return true;
}

// get pseudopototentials for particles on torus with SU(4) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as pplus-upplus, upplus-upminus, upplus-downplus, upplus-downminus, 
//                                                                     upminus-upminus, upminus-downplus, upminus-downminus, downplus-downplus, downplus-downminus, downminus-downminus)
// return value = true if no error occured

bool FQHETorusSU4GetPseudopotentials (char* fileName, int* nbrPseudoPotentials, double** pseudoPotentials)
{
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }
  int TmpNbrPseudoPotentials;
  double* TmpPseudoPotentials;
  bool Flag = false;
  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      for (int i = 0; i < 10; ++i)
	{
	  nbrPseudoPotentials[i] = TmpNbrPseudoPotentials;
	  pseudoPotentials[i] = new double[nbrPseudoPotentials[i]];
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    pseudoPotentials[i][j] = TmpPseudoPotentials[j];
	}
    }
  else
    if (InteractionDefinition["Pseudopotentials"] != 0)
      {
	cout << "Pseudopotentials has a wrong value in " << fileName << endl;
	return false;
      }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpPlusUpPlus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[0] = TmpNbrPseudoPotentials;
      pseudoPotentials[0] = new double[nbrPseudoPotentials[0]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[0][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsUpPlusUpPlus"] != 0)
      {
	cout << "PseudopotentialsUpPlusUpPlus has a wrong value in " << fileName << endl;
	return false;
      }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpPlusUpMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[1] = TmpNbrPseudoPotentials;
      pseudoPotentials[1] = new double[nbrPseudoPotentials[1]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[1][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsUpPlusUpMinus"] != 0)
      {
	cout << "PseudopotentialsUpPlusUpMinus has a wrong value in " << fileName << endl;
	return false;
      }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpPlusDownPlus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[2] = TmpNbrPseudoPotentials;
      pseudoPotentials[2] = new double[nbrPseudoPotentials[2]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[2][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsUpPlusDownPlus"] != 0)
      {
	cout << "PseudopotentialsUpPlusDownPlus has a wrong value in " << fileName << endl;
	return false;
      }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpPlusDownMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[3] = TmpNbrPseudoPotentials;
      pseudoPotentials[3] = new double[nbrPseudoPotentials[3]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[3][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsUpPlusDownMinus"] != 0)
      {
	cout << "PseudopotentialsUpPlusDownMinus has a wrong value in " << fileName << endl;
	return false;
      }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpMinusUpMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[4] = TmpNbrPseudoPotentials;
      pseudoPotentials[4] = new double[nbrPseudoPotentials[4]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[4][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsUpMinusUpMinus"] != 0)
      {
	cout << "PseudopotentialsUpMinusUpMinus has a wrong value in " << fileName << endl;
	return false;
      }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpMinusDownPlus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[5] = TmpNbrPseudoPotentials;
      pseudoPotentials[5] = new double[nbrPseudoPotentials[5]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[5][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsUpMinusDownPlus"] != 0)
      {
	cout << "PseudopotentialsUpMinusDownPlus has a wrong value in " << fileName << endl;
	return false;
      }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpMinusDownMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[6] = TmpNbrPseudoPotentials;
      pseudoPotentials[6] = new double[nbrPseudoPotentials[6]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[6][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsUpMinusDownMinus"] != 0)
      {
	cout << "PseudopotentialsUpMinusDownMinus has a wrong value in " << fileName << endl;
	return false;
      }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownPlusDownPlus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[7] = TmpNbrPseudoPotentials;
      pseudoPotentials[7] = new double[nbrPseudoPotentials[7]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[7][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsDownPlusDownPlus"] != 0)
      {
	cout << "PseudopotentialsDownPlusDownPlus has a wrong value in " << fileName << endl;
	return false;
      }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownPlusDownMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[8] = TmpNbrPseudoPotentials;
      pseudoPotentials[8] = new double[nbrPseudoPotentials[8]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[8][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsDownPlusDownMinus"] != 0)
      {
	cout << "PseudopotentialsDownPlusDownMinus has a wrong value in " << fileName << endl;
	return false;
      }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownMinusDownMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[9] = TmpNbrPseudoPotentials;
      pseudoPotentials[9] = new double[nbrPseudoPotentials[9]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[9][j] = TmpPseudoPotentials[j];
    }
  else
    {
      if (InteractionDefinition["PseudopotentialsDownMinusDownMinus"] != 0)
	{
	  cout << "PseudopotentialsDownMinusDownMinus has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  return true;
}

