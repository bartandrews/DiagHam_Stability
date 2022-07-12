////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             set of functions used to squeezed (aka Haldane) basis          //
//                                                                            //
//                        last modification : 24/02/2009                      //
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
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>
#include <cstdlib>


using std::cout;
using std::endl;


// get the root parition from a file
// 
// rootFileName = name of the file that contains the root description
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the maximum Lz value
// referenceState = array where the root partition description will be stored
// return value = true if no error occured

bool FQHEGetRootPartition (char* rootFileName, int& nbrParticles, int& lzMax, int*& referenceState)
{
  ConfigurationParser ReferenceStateDefinition;
  if (ReferenceStateDefinition.Parse(rootFileName) == false)
    {
      ReferenceStateDefinition.DumpErrors(cout) << endl;
      return false;
    }
  if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", nbrParticles) == false) || (nbrParticles <= 0))
    {
      cout << "NbrParticles is not defined or as a wrong value" << endl;
      return false;
    }
  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", lzMax) == false) || (lzMax < 0))
    {
      cout << "LzMax is not defined or as a wrong value" << endl;
      return false;
    }
  int MaxNbrLz;
  if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', referenceState, MaxNbrLz) == false)
    {
      cout << "error while parsing ReferenceState in " << rootFileName << endl;
      return false;     
    }
  if (MaxNbrLz != (lzMax + 1))
    {
      cout << "wrong LzMax value in ReferenceState" << endl;
      return false;     
    }
  return true;
}

// auxiliary routine used in FQHEGetRootPartitionSU2
// lzMax = maxLz
// referenceState = string valued reference state with extended syntax
int CountExtendedRootStates(int lzMax, char **referenceState)
{
  int NbrUp=0;
  int NbrDown=0;
  int NbrSinglet=0;
  for (int j = 0; j <= lzMax; ++j)
    {
      if ((referenceState[j][0] == 'x') || (referenceState[j][0] == 'X'))
	{
	  ++NbrUp;
	  ++NbrDown;
	}
	else
	  if ((referenceState[j][0] == 'u') || (referenceState[j][0] == 'U'))
	    ++NbrUp;
	  else
	    if ((referenceState[j][0] == 'd') || (referenceState[j][0] == 'D'))
	      ++NbrDown;
	    else
	      if ((referenceState[j][0] == 's') || (referenceState[j][0] == 'S'))
		++NbrSinglet;
    }
  if (NbrSinglet&1)
    {
      cout << "Number of electrons in a singlet state must be an even number for all roots";
      exit(-1);
    }
  FactorialCoefficient Bico(1);
  Bico.FactorialMultiply(NbrSinglet);
  Bico.FactorialDivide(NbrSinglet/2);
  Bico.FactorialDivide(NbrSinglet/2);
  return (int)Bico.GetIntegerValue();
}

// recursively write all root states for known positions of singlets
// nbrSinglets = numbers of singlets present
// singlets = positions of singlets
// lzMax = maximum Lz
// background = invariant part of state
// referenceStates = array to write reference states
// currentNbrSinglet = singlets remaining to be filled
// currentNbrUp = up-spins to be filled
// pos = number of first root state in array referenceStates to be written to
// return value = position from which new roots have to be stored
int ShiftedGenerateExtendedRootStates(int nbrSinglets, int* singlets, int lzMax, int *background, int **&referenceStates, int currentNbrSinglet, int currentNbrUp, int pos)
{
  if (currentNbrSinglet == 0) return pos;

  if (currentNbrSinglet == 1)
    {
      referenceStates[pos] = new int [lzMax + 1];
      for (int j=0; j<=lzMax; ++j)
	referenceStates[pos][j]=background[j];
      if (currentNbrUp>0)	
	referenceStates[pos][singlets[0]]=2;
      else
	referenceStates[pos][singlets[0]]=1;
      return pos+1;
    }

  if (currentNbrSinglet-currentNbrUp>0)
    {
      // place down spin:
      int TmpPos = ShiftedGenerateExtendedRootStates(nbrSinglets, singlets, lzMax, background, referenceStates, currentNbrSinglet-1, currentNbrUp, pos);
      for (int i = pos; i < TmpPos; i++)
	referenceStates[i][singlets[currentNbrSinglet-1]] = 1;
//       if (currentNbrSinglet==nbrSinglets)
// 	{
// 	  for (int i = pos; i < TmpPos; i++)
// 	    {
// 	      cout << "Root "<<i<<": ";
// 	      for (int j=0; j<=lzMax; ++j)
// 		cout <<referenceStates[i][j]<<" ";
// 	      cout << endl;
// 	    }
// 	}
      pos=TmpPos;
    }  
  if (currentNbrUp>0)
    {
      // place up spin:
      int TmpPos = ShiftedGenerateExtendedRootStates(nbrSinglets, singlets, lzMax, background, referenceStates, currentNbrSinglet-1, currentNbrUp-1, pos);
      for (int i = pos; i < TmpPos; i++)
	referenceStates[i][singlets[currentNbrSinglet-1]] = 2;
//       if (currentNbrSinglet==nbrSinglets)
// 	{
// 	  for (int i = pos; i < TmpPos; i++)
// 	    {
// 	      cout << "Root "<<i<<": ";
// 	      for (int j=0; j<=lzMax; ++j)
// 		cout <<referenceStates[i][j]<<" ";
// 	      cout << endl;
// 	    }
// 	}
      pos=TmpPos;
    }
  return pos;	  
}

// generate all root states from a generalized form
// lzMax = maximum Lz
// pos = number of first root state in array referenceStates
// referenceStates = array to write reference states
// extReferenceState = extended reference state
void GenerateExtendedRootStates(int lzMax, int& pos, int **&referenceStates, char **extReferenceState)
{
  int *Singlets = new int[lzMax+1];
  int NbrSinglets = 0;
  int *Background = new int[lzMax+1];
  for (int j = 0; j <= lzMax; ++j)
    {
      if (extReferenceState[j][0] == '0')
	Background[j] = 0;
      else
	if ((extReferenceState[j][0] == 'x') || (extReferenceState[j][0] == 'X'))
	  Background[j] = 3;
	else
	  if ((extReferenceState[j][0] == 'u') || (extReferenceState[j][0] == 'U'))
	    Background[j] = 2;
	  else
	    if ((extReferenceState[j][0] == 'd') || (extReferenceState[j][0] == 'D'))
	      Background[j] = 1;
	    else
	      if ((extReferenceState[j][0] == 's') || (extReferenceState[j][0] == 'S'))
		{
		  Singlets[NbrSinglets]=j;
		  ++NbrSinglets;
		}
    }
  pos+=ShiftedGenerateExtendedRootStates(NbrSinglets, Singlets, lzMax, Background, referenceStates, NbrSinglets, NbrSinglets/2, pos);
  delete [] Singlets;
  delete [] Background;
}

// get the root partition from a file in the SU2 case
// 
// rootFileName = name of the file that contains the root description
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the maximum Lz value
// referenceStates = array where the root partition descriptions will be stored
// nbrReferenceStates = number of root partitions that have been extracted
// return value = true if no error occured

bool FQHEGetRootPartitionSU2 (char* rootFileName, int& nbrParticles, int& lzMax, 
			      int**& referenceStates, int& nbrReferenceStates)
{
  ConfigurationParser ReferenceStateDefinition;
  if (ReferenceStateDefinition.Parse(rootFileName) == false)
    {
      ReferenceStateDefinition.DumpErrors(cout) << endl;
      return false;
    }
  if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", nbrParticles) == false) || (nbrParticles <= 0))
    {
      cout << "NbrParticles is not defined or as a wrong value" << endl;
      return false;
    }
  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", lzMax) == false) || (lzMax < 0))
    {
      cout << "LzMax is not defined or as a wrong value" << endl;
      return false;
    }
  char*** TmpReferenceStates;
  int* NbrLzMaxPerStates;
  if (ReferenceStateDefinition.GetAsStringMultipleArray("ReferenceStates", '|', ' ', TmpReferenceStates, nbrReferenceStates, NbrLzMaxPerStates) == false)
    {
      cout << "error while parsing ReferenceStates in " << rootFileName << endl;
      return false;     
    }
  bool HaveExtendedFormat=false;
  bool *StateHaveExtendedFormat=new bool[nbrReferenceStates];
  int TotalReferenceStates=0;
  for (int i = 0; i < nbrReferenceStates; ++i)
    {
      StateHaveExtendedFormat[i]=false;
      if (NbrLzMaxPerStates[i] != (lzMax + 1))
	{
	  cout << "wrong LzMax value in ReferenceState " << i << endl;
	  return false;     
	}
      for (int j = 0; j <= lzMax; ++j)
	{
	  if ((TmpReferenceStates[i][j][0] == 's') || (TmpReferenceStates[i][j][0] == 'S'))
	    {
	      StateHaveExtendedFormat[i]=true;
	      HaveExtendedFormat=true;
	    }
	}
      if (StateHaveExtendedFormat[i])
	TotalReferenceStates+=CountExtendedRootStates(lzMax, TmpReferenceStates[i]);
      else
	++TotalReferenceStates;
    }
  cout << "total: "<<TotalReferenceStates<<" root states"<<endl;
  referenceStates = new int*[TotalReferenceStates];
  int pos=0;
  for (int i = 0; i < nbrReferenceStates; ++i)
    {
      if (StateHaveExtendedFormat[i])
	{
	  GenerateExtendedRootStates(lzMax, pos, referenceStates, TmpReferenceStates[i]);
	}
      else
	{
	  referenceStates[pos] = new int [lzMax + 1];
	  for (int j = 0; j <= lzMax; ++j)
	    {
	      if (TmpReferenceStates[i][j][0] == '0')
		referenceStates[pos][j] = 0;
	      else
		if ((TmpReferenceStates[i][j][0] == 'x') || (TmpReferenceStates[i][j][0] == 'X'))
		  referenceStates[pos][j] = 3;
		else
		  if ((TmpReferenceStates[i][j][0] == 'u') || (TmpReferenceStates[i][j][0] == 'U'))
		    referenceStates[pos][j] = 2;
		  else
		    if ((TmpReferenceStates[i][j][0] == 'd') || (TmpReferenceStates[i][j][0] == 'D'))
		      referenceStates[pos][j] = 1;	  
	    }
	  ++pos;
	}
    }
  nbrReferenceStates = TotalReferenceStates;
  return true;
}


// get the root partition from a file in the SU2 case
// 
// rootFileName = name of the file that contains the root description
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the maximum Lz value
// referenceStates = array where the root partition descriptions will be stored
// nbrReferenceStates = number of root partitions that have been extracted
// texturelessFlag = flag to indicate whether or not to consider spin texture when performing squeezing
// return value = true if no error occured

bool FQHEGetRootPartitionSU2 (char* rootFileName, int& nbrParticles, int& lzMax, 
			      int**& referenceStates, int& nbrReferenceStates, bool &texturelessFlag)
{
  ConfigurationParser ReferenceStateDefinition;
  if (ReferenceStateDefinition.Parse(rootFileName) == false)
    {
      ReferenceStateDefinition.DumpErrors(cout) << endl;
      return false;
    }
  if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", nbrParticles) == false) || (nbrParticles <= 0))
    {
      cout << "NbrParticles is not defined or as a wrong value" << endl;
      return false;
    }
  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", lzMax) == false) || (lzMax < 0))
    {
      cout << "LzMax is not defined or as a wrong value" << endl;
      return false;
    }
  if (ReferenceStateDefinition.GetAsBoolean("Textureless", texturelessFlag) == false)
    {
      texturelessFlag = false;
    }  
  char*** TmpReferenceStates;
  int* NbrLzMaxPerStates;
  if (ReferenceStateDefinition.GetAsStringMultipleArray("ReferenceStates", '|', ' ', TmpReferenceStates, nbrReferenceStates, NbrLzMaxPerStates) == false)
    {
      cout << "error while parsing ReferenceStates in " << rootFileName << endl;
      return false;     
    }
  bool HaveExtendedFormat=false;
  bool *StateHaveExtendedFormat=new bool[nbrReferenceStates];
  int TotalReferenceStates=0;
  for (int i = 0; i < nbrReferenceStates; ++i)
    {
      StateHaveExtendedFormat[i]=false;
      if (NbrLzMaxPerStates[i] != (lzMax + 1))
	{
	  cout << "wrong LzMax value in ReferenceState " << i << endl;
	  return false;     
	}
      for (int j = 0; j <= lzMax; ++j)
	{
	  if ((TmpReferenceStates[i][j][0] == 's') || (TmpReferenceStates[i][j][0] == 'S'))
	    {
	      StateHaveExtendedFormat[i]=true;
	      HaveExtendedFormat=true;
	    }
	}
      if (StateHaveExtendedFormat[i])
	TotalReferenceStates+=CountExtendedRootStates(lzMax, TmpReferenceStates[i]);
      else
	++TotalReferenceStates;
    }
  cout << "total: "<<TotalReferenceStates<<" root states"<<endl;
  referenceStates = new int*[TotalReferenceStates];
  int pos=0;
  for (int i = 0; i < nbrReferenceStates; ++i)
    {
      if (StateHaveExtendedFormat[i])
	{
	  GenerateExtendedRootStates(lzMax, pos, referenceStates, TmpReferenceStates[i]);
	}
      else
	{
	  referenceStates[pos] = new int [lzMax + 1];
	  for (int j = 0; j <= lzMax; ++j)
	    {
	      if (TmpReferenceStates[i][j][0] == '0')
		referenceStates[pos][j] = 0;
	      else
		if ((TmpReferenceStates[i][j][0] == 'x') || (TmpReferenceStates[i][j][0] == 'X'))
		  referenceStates[pos][j] = 3;
		else
		  if ((TmpReferenceStates[i][j][0] == 'u') || (TmpReferenceStates[i][j][0] == 'U'))
		    referenceStates[pos][j] = 2;
		  else
		    if ((TmpReferenceStates[i][j][0] == 'd') || (TmpReferenceStates[i][j][0] == 'D'))
		      referenceStates[pos][j] = 1;	  
	    }
	  ++pos;
	}
    }
  nbrReferenceStates = TotalReferenceStates;
  return true;
}



