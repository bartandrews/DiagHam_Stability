////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         set of functions used to handle sphere pseudopotential files       //
//                                                                            //
//                        last modification : 26/02/2009                      //
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
#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"
#include "GeneralTools/ConfigurationParser.h"


#include <iostream>


using std::cout;
using std::endl;


// get pseudopototentials for particles on sphere with SU(2) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMax = reference on twice the maximum Lz value
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// return value = true if no error occured

bool FQHESphereSU2GetPseudopotentials (char* fileName, int lzMax, double** pseudoPotentials,
				       double* oneBodyPotentialUpUp, double* oneBodyPotentialDownDown)
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
      if (TmpNbrPseudoPotentials > (lzMax +1))
	{
	  cout << "Invalid number of pseudo-potentials in Pseudopotentials" << endl;
	  return false;	  
	}
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	  pseudoPotentials[i][j] = TmpPseudoPotentials[j];
      if (TmpNbrPseudoPotentials <= lzMax)
	{
	  cout << "warning : number of pseudo-potentials in Pseudopotentials is lower than expected, padding with zeros" << endl;
	  for (int i = 0; i < 3; ++i)
	    for (int j = TmpNbrPseudoPotentials; j <= lzMax; ++j)
	      pseudoPotentials[i][j] = 0.0;	  
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
      if (TmpNbrPseudoPotentials != (lzMax +1))
	{
	  cout << "Invalid number of pseudo-potentials in PseudopotentialsUpUp" << endl;
	  return false;	  
	}
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[0][j] = TmpPseudoPotentials[j];
      if (TmpNbrPseudoPotentials <= lzMax)
	{
	  cout << "warning : number of pseudo-potentials in PseudopotentialsUpUp is lower than expected, padding with zeros" << endl;
	  for (int j = TmpNbrPseudoPotentials; j <= lzMax; ++j)
	    pseudoPotentials[0][j] = 0.0;	  
	}
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
      if (TmpNbrPseudoPotentials > (lzMax +1))
	{
	  cout << "Invalid number of pseudo-potentials in PseudopotentialsDownDown" << endl;
	  return false;	  
	}
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[1][j] = TmpPseudoPotentials[j];
      if (TmpNbrPseudoPotentials <= lzMax)
	{
	  cout << "warning : number of pseudo-potentials in PseudopotentialsDownDown is lower than expected, padding with zeros" << endl;
	  for (int j = TmpNbrPseudoPotentials; j <= lzMax; ++j)
	    pseudoPotentials[1][j] = 0.0;
	}
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
      if (TmpNbrPseudoPotentials > (lzMax +1))
	{
	  cout << "Invalid number of pseudo-potentials in PseudopotentialsUpDown" << endl;
	  return false;	  
	}
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[2][j] = TmpPseudoPotentials[j];
      if (TmpNbrPseudoPotentials <= lzMax)
	{
	  cout << "warning : number of pseudo-potentials in PseudopotentialsUpDown is lower than expected, padding with zeros" << endl;
	  for (int j = TmpNbrPseudoPotentials; j <= lzMax; ++j)
	    pseudoPotentials[2][j] = 0.0;
	}
    }
  else
    if (InteractionDefinition["PseudopotentialsUpDown"] != 0)
      {
	cout << "PseudopotentialsUpDown has a wrong value in " << fileName << endl;
	return false;
      }
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialUpUp", ' ', oneBodyPotentialUpUp, TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != (lzMax + 1))
	{
	  cout << "OneBodyPotentialUpUp has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialDownDown", ' ', oneBodyPotentialDownDown, TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != (lzMax + 1))
	{
	  cout << "OneBodyPotentialDownDown has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  return true;
}

