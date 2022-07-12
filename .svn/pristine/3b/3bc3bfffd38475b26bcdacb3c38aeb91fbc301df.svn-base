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
#include <string>


using std::cout;
using std::endl;
using std::string;


// get pseudopototentials for spinless particles on sphere from file
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMax = reference on twice the maximum Lz value
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
// onebodyPotential =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state), null pointer if none
// return value = true if no error occured

bool FQHESphereGetPseudopotentials (char* fileName, int lzMax, double* pseudoPotentials, double*& oneBodyPotential)
{
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }
  if (InteractionDefinition["Pseudopotentials"] == 0)
    {
      cout << "Pseudopotentials are not defined in " << fileName << endl;
      return false;
    }  
  int TmpNbrPseudoPotentials;
  double* TmpPseudoPotentials;
  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax < (lzMax + 1))
	{
	  cout << "warning, Pseudopotentials has less entries than the number of orbitals, padding with zeroes" << endl;
	  for (int j = 0; j < TmpMax; ++j)
	    pseudoPotentials[j] = TmpPseudoPotentials[j];
 	  for (int j = TmpMax; j <= lzMax; ++j)
	    pseudoPotentials[j] = 0.0; 
	}      
      else
	{
	  if (TmpMax > (lzMax + 1))
	    {
	      TmpMax = lzMax + 1;
	      cout << "warning, Pseudopotentials has more entries than the number of orbitals and will be truncated" << endl;
	    }
	  for (int j = 0; j <= lzMax; ++j)
	    pseudoPotentials[j] = TmpPseudoPotentials[j];
	}
    }
  else
    {
      cout << "Pseudopotentials has a wrong value in " << fileName << endl;
      return false;
    }
  oneBodyPotential = 0;
  if (InteractionDefinition.GetAsDoubleArray("Onebodypotentials", ' ', oneBodyPotential, TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != (lzMax + 1))
	{
	  cout << "OneBodyPotential has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  else
    {
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentials", ' ', oneBodyPotential, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (lzMax + 1))
	    {
	      cout << "OneBodyPotential has a wrong number of components or has a wrong value in " << fileName << endl;
	      return false;
	    }
	}
    }
  return true;
}

// get the one-body potential for spinless particles on sphere from file
// 
// fileName = name of the file that contains the one-body potential description
// lzMax = reference on twice the maximum Lz value
// onebodyPotential =  reference on the one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state)
// return value = true if no error occured

bool FQHESphereGetOneBodyPotentials (char* fileName, int& lzMax, double*& oneBodyPotential)
{
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }
  oneBodyPotential = 0;
  if (InteractionDefinition.GetAsDoubleArray("Onebodypotentials", ' ', oneBodyPotential, lzMax) == false)
    {
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentials", ' ', oneBodyPotential, lzMax) == false)
	{
	  cout << "OneBodyPotentials Onebodypotentials or  is not defined in " << fileName << endl;
	  return false;
	}
    }
  --lzMax;
  return true;
}

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
				       double*& oneBodyPotentialUpUp, double*& oneBodyPotentialDownDown)
{
  double *TmpOneBodyPseudopotentialUpDown=NULL;

  bool Rst = FQHESphereSU2GetPseudopotentials (fileName, lzMax, pseudoPotentials,
					       oneBodyPotentialUpUp, oneBodyPotentialDownDown, TmpOneBodyPseudopotentialUpDown);
  if (TmpOneBodyPseudopotentialUpDown!=NULL)
    delete[] TmpOneBodyPseudopotentialUpDown;
  return Rst;
}


// get pseudopototentials for particles on sphere with SU(2) spin from file, including a tunneling term
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMax = reference on twice the maximum Lz value
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyPotentialUpDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// return value = true if no error occured

bool FQHESphereSU2GetPseudopotentials (char* fileName, int lzMax, double** pseudoPotentials,
				       double*& oneBodyPseudopotentialUpUp, double*& oneBodyPseudopotentialDownDown, double*& oneBodyPseudopotentialUpDown)
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
      if (TmpNbrPseudoPotentials > (lzMax +1))
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
  // section needed only in all-sz modes
  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsPairTunneling", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      if (TmpNbrPseudoPotentials > (lzMax +1))
	{
	  cout << "Invalid number of pseudo-potentials in PseudopotentialsPairTunneling" << endl;
	  return false;
	}
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[3][j] = TmpPseudoPotentials[j];
      if (TmpNbrPseudoPotentials <= lzMax)
	{
	  cout << "warning : number of pseudo-potentials in PseudopotentialsPairTunneling is lower than expected, padding with zeros" << endl;
	  for (int j = TmpNbrPseudoPotentials; j <= lzMax; ++j)
	    pseudoPotentials[3][j] = 0.0;
	}
    }
  else
    if (InteractionDefinition["PseudopotentialsPairTunneling"] != 0)
      {
	cout << "PseudopotentialsPairTunneling has a wrong value in " << fileName << endl;
	return false;
      }
  // end all-sz insertion
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialUpUp", ' ', oneBodyPseudopotentialUpUp, TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != (lzMax + 1))
	{
	  cout << "OneBodyPotentialUpUp has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialDownDown", ' ', oneBodyPseudopotentialDownDown, TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != (lzMax + 1))
	{
	  cout << "OneBodyPotentialDownDown has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialUpDown", ' ', oneBodyPseudopotentialUpDown, TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != (lzMax + 1))
	{
	  cout << "OneBodyPotentialUpDown has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  return true;
}

// get pseudopototentials for particles on sphere with SU(2) spin from file, including a tunneling term and pairing
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMax = reference on twice the maximum Lz value
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyPotentialUpDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyPotentialUpDown = reference on the one body potential for the pairing term (sorted from component on the lowest Lz state to component on the highest Lz state), null pointer if none
// return value = true if no error occured

bool FQHESphereSU2GetPseudopotentialsWithPairing (char* fileName, int lzMax, double** pseudoPotentials,
						  double*& oneBodyPseudopotentialUpUp, double*& oneBodyPseudopotentialDownDown, double*& oneBodyPseudopotentialUpDown,
						  double*& oneBodyPseudopotentialPairing)
{
  if (FQHESphereSU2GetPseudopotentials(fileName, lzMax, pseudoPotentials,
				       oneBodyPseudopotentialUpUp, oneBodyPseudopotentialDownDown, oneBodyPseudopotentialUpDown) == false)
    {
      return false;
    }
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }
  int TmpNbrPseudoPotentials = 0;
  if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialPairing", ' ', oneBodyPseudopotentialPairing, TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != (lzMax + 1))
	{
	  cout << "OneBodyPotentialPairing has a wrong number of components or has a wrong value in " << fileName << endl;
	  return false;
	}
    }
  return true;
}

// get pseudopototentials for particles on sphere with two landau levels
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMax = reference on twice the maximum Lz value of the LLL
// pseudoPotentials = array or arrays of pseudo-potentials. 9 in total which go in ascending order of index p = 0-8 where label p = l*3 + r where l and r label the ll indices on the left and right of the 
//                   interaction respectively and take values 0:up-up, 1: down-down, 2: up-down. Assumed that space is already allocated.
// return value = true if no error occured

bool FQHESphereTwoLandauLevelGetPseudopotentials (char* fileName, int lzMax, double** pseudoPotentials)
{
  // these are the labels of the arrays as they will be in the file.
  string PseudoLabels[10] = {"PseudopotentialsUpUpUpUp","PseudopotentialsUpUpDownDown","PseudopotentialsUpUpUpDown",
			    "PseudopotentialsDownDownUpUp","PseudopotentialsDownDownDownDown","PseudopotentialsDownDownUpDown",
			    "PseudopotentialsUpDownUpUp","PseudopotentialsUpDownDownDown","PseudopotentialsUpDownUpDown","PseudopotentialsUpDownDownUp"};
  // these are the lenghts of the arrays corresponding to the labels above. 			    
  int PseudoLengths[10] = { lzMax+3, lzMax+1, lzMax+1, lzMax+1, lzMax+1, lzMax, lzMax+1, lzMax, lzMax+1, lzMax+1}; 
  
  
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }
  int TmpNbrPseudoPotentials;
  double* TmpPseudoPotentials;
  for ( int Idx = 0 ; Idx < 10; Idx++ ) 
    {
      if (InteractionDefinition.GetAsDoubleArray(PseudoLabels[Idx].c_str(), ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials > PseudoLengths[Idx])
	    {
	      cout << "Invalid number of pseudo-potentials in " << PseudoLabels[Idx]  << endl;
	      return false;	  
	    }
	  if (TmpNbrPseudoPotentials < PseudoLengths[Idx])
	    {
	      cout << "warning : number of pseudo-potentials in " << PseudoLabels[Idx] << " is lower than expected, padding with zeros" << endl;	      
	      for (int j = TmpNbrPseudoPotentials; j < PseudoLengths[Idx]; ++j)
		  pseudoPotentials[Idx][j] = 0.0;	  
	    }
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    pseudoPotentials[Idx][j] = TmpPseudoPotentials[j];
	  delete [] TmpPseudoPotentials;
  
	}
  else if (InteractionDefinition[PseudoLabels[Idx].c_str()] != 0)
      {
	cout << PseudoLabels[Idx] << " has a wrong value in " << fileName << endl;
	return false;
      }
    }
  
  return true;
}

// get pseudopototentials for particles on sphere with SU(4) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMax = twice the maximum Lz value of the LLL
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                    first index refered to the spin sector (sorted as 1-1, 1-2, 1-3, 2-2, 2-3, 3-3)
// return value = true if no error occured

bool FQHESphereSU3GetPseudopotentials (char* fileName, int lzMax, double** pseudoPotentials)
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
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, Pseudopotentials has more entries than the number of orbitals and will be truncated" << endl;
	}
      for (int i = 0; i < 6; ++i)
	{
	  pseudoPotentials[i] = new double[lzMax + 1];
	  for (int j = 0; j < TmpMax; ++j)
	    pseudoPotentials[i][j] = TmpPseudoPotentials[j];
	  for (int j = TmpMax; j <= lzMax; ++j)
	    pseudoPotentials[i][j] = 0.0;	  
	}
    }
  else
    {
      if (InteractionDefinition["Pseudopotentials"] != 0)
	{
	  cout << "Pseudopotentials has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials11", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, Pseudopotentials11 has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[0] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[0][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[0][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["Pseudopotentials11"] != 0)
	{
	  cout << "Pseudopotentials11 has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials12", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, Pseudopotentials12 has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[1] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[1][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[1][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["Pseudopotentials12"] != 0)
	{
	  cout << "Pseudopotentials12 has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials13", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, Pseudopotentials13 has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[2] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[2][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[2][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["Pseudopotentials13"] != 0)
	{
	  cout << "Pseudopotentials13 has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials22", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, Pseudopotentials22 has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[3] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[3][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[3][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["Pseudopotentials22"] != 0)
	{
	  cout << "Pseudopotentials22 has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials23", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, Pseudopotentials23 has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[4] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[4][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[4][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["Pseudopotentials23"] != 0)
	{
	  cout << "Pseudopotentials23 has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials33", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, Pseudopotentials33 has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[5] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[5][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[5][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["Pseudopotentials33"] != 0)
	{
	  cout << "Pseudopotentials33 has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  return true;
}

// get pseudopototentials for particles on sphere with SU(4) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// lzMax = twice the maximum Lz value of the LLL
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as pplus-upplus, upplus-upminus, upplus-downplus, upplus-downminus, 
//                                                                     upminus-upminus, upminus-downplus, upminus-downminus, downplus-downplus, downplus-downminus, downminus-downminus)
// return value = true if no error occured

bool FQHESphereSU4GetPseudopotentials (char* fileName, int lzMax, double** pseudoPotentials)
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
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, Pseudopotentials has more entries than the number of orbitals and will be truncated" << endl;
	}
      for (int i = 0; i < 10; ++i)
	{
	  pseudoPotentials[i] = new double[lzMax + 1];
	  for (int j = 0; j < TmpMax; ++j)
	    pseudoPotentials[i][j] = TmpPseudoPotentials[j];
	  for (int j = TmpMax; j <= lzMax; ++j)
	    pseudoPotentials[i][j] = 0.0;	  
	}
    }
  else
    {
      if (InteractionDefinition["Pseudopotentials"] != 0)
	{
	  cout << "Pseudopotentials has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpPlusUpPlus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, PseudopotentialsUpPlusUpPlus has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[0] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[0][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[0][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["PseudopotentialsUpPlusUpPlus"] != 0)
	{
	  cout << "PseudopotentialsUpPlusUpPlus has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpPlusUpMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, PseudopotentialsUpPlusUpMinus has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[1] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[1][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[1][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["PseudopotentialsUpPlusUpMinus"] != 0)
	{
	  cout << "PseudopotentialsUpPlusUpMinus has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpPlusDownPlus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, PseudopotentialsUpPlusDownPlus has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[2] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[2][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[2][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["PseudopotentialsUpPlusDownPlus"] != 0)
	{
	  cout << "PseudopotentialsUpPlusDownPlus has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpPlusDownMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, PseudopotentialsUpPlusDownMinus has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[3] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[3][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[3][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["PseudopotentialsUpPlusDownMinus"] != 0)
	{
	  cout << "PseudopotentialsUpPlusDownMinus has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpMinusUpMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, PseudopotentialsUpMinusUpMinus has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[4] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[4][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[4][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["PseudopotentialsUpMinusUpMinus"] != 0)
	{
	  cout << "PseudopotentialsUpMinusUpMinus has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpMinusDownPlus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, PseudopotentialsUpMinusDownPlus has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[5] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[5][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[5][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["PseudopotentialsUpMinusDownPlus"] != 0)
	{
	  cout << "PseudopotentialsUpMinusDownPlus has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpMinusDownMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, PseudopotentialsUpMinusDownMinus has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[6] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[6][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[6][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["PseudopotentialsUpMinusDownMinus"] != 0)
	{
	  cout << "PseudopotentialsUpMinusDownMinus has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownPlusDownPlus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, PseudopotentialsDownPlusDownPlus has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[7] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[7][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[7][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["PseudopotentialsDownPlusDownPlus"] != 0)
	{
	  cout << "PseudopotentialsDownPlusDownPlus has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownPlusDownMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, PseudopotentialsDownPlusDownMinus has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[8] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[8][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[8][j] = 0.0;	  
    }
  else
    {
      if (InteractionDefinition["PseudopotentialsDownPlusDownMinus"] != 0)
	{
	  cout << "PseudopotentialsDownPlusDownMinus has a wrong value in " << fileName << endl;
	  return false;
	}
    }

  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownMinusDownMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      int TmpMax = TmpNbrPseudoPotentials;
      if (TmpMax > (lzMax + 1))
	{
	  TmpMax = lzMax + 1;
	  cout << "warning, PseudopotentialsDownMinusDownMinus has more entries than the number of orbitals and will be truncated" << endl;
	}
      pseudoPotentials[9] = new double[lzMax + 1];
      for (int j = 0; j < TmpMax; ++j)
	pseudoPotentials[9][j] = TmpPseudoPotentials[j];
      for (int j = TmpMax; j <= lzMax; ++j)
	pseudoPotentials[9][j] = 0.0;	  
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
