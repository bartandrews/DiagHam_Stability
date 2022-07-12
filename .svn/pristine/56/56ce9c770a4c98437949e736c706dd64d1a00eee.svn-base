#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/StringTools.h"

#include "Options/Options.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <cstring> 


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHETorusBosonsFilterOccupation" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  //  ArchitectureManager Architecture;

  Manager += SystemGroup;
  //  Architecture.AddOptionGroup(&Manager);
//  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
//  Manager += ToolsGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "vector file that corresponds to the input state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-occ", "maximum occupation per obital", 1);
  (*OutputGroup) += new BooleanOption  ('\n', "discard-vector", "do not store the filtered state");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusBosonsFilterOccupation -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles;
  int NbrFluxQuanta;
  int KxMomentum;
  int KyMomentum;
  double Ratio = 1.0;
  bool Statistics = true;
  bool MagneticTranslationFlag = true;
  if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						  NbrParticles, NbrFluxQuanta, KxMomentum, KyMomentum, Statistics) == false)
    {
      MagneticTranslationFlag = false;
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						  NbrParticles, NbrFluxQuanta, KxMomentum, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state") << endl;
	  return -1;
	}
    }

  cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << NbrFluxQuanta << " Ky=" << KyMomentum << " ";
  if (Statistics == true)
    {
      cout << "FQHETorusBosonsFilterOccupation is only relevant for bosonic states" << endl;
    }

  if (MagneticTranslationFlag == true)
    {
      BosonOnTorusWithMagneticTranslationsShort* SpaceWithTranslations = 0;
      SpaceWithTranslations = new BosonOnTorusWithMagneticTranslationsShort(NbrParticles, NbrFluxQuanta, KxMomentum, KyMomentum);
      ComplexVector InputStateWithTranslations;
      if (InputStateWithTranslations.ReadVector(Manager.GetString("input-state")) == false)
	{
	  cout << "error while reading " << Manager.GetString("input-state") << endl;
	  return -1;
	}
      if (InputStateWithTranslations.GetVectorDimension() != SpaceWithTranslations->GetHilbertSpaceDimension())
	{
	  cout << "error: vector and Hilbert-space have unequal dimensions " << InputStateWithTranslations.GetVectorDimension() 
	       << " " << SpaceWithTranslations->GetHilbertSpaceDimension() << endl;
	  return -1;
	}
      unsigned long* TmpState = new unsigned long [NbrFluxQuanta];
      ComplexVector OutputStateWithTranslations (SpaceWithTranslations->GetLargeHilbertSpaceDimension(), true);
      double TotalWeight = 0.0;
      unsigned long MaxOccupation = (unsigned long) Manager.GetInteger("max-occ");
      for (long i = 0l; i < SpaceWithTranslations->GetLargeHilbertSpaceDimension(); ++i)
	{
	  SpaceWithTranslations->GetOccupationNumber(i, TmpState);
	  int Pos = 0;
	  while ((Pos < NbrFluxQuanta) && (TmpState[Pos] <= MaxOccupation))
	    ++Pos;
	  if (Pos == NbrFluxQuanta)
	    {
	      OutputStateWithTranslations[i] = InputStateWithTranslations[i];
	      TotalWeight += SqrNorm(InputStateWithTranslations[i]);
	    }
	}
      cout << "total weight with occupation <= " << MaxOccupation << " : " <<  TotalWeight << endl;
      if (Manager.GetBoolean("discard-vector") == false)
	{
	  OutputStateWithTranslations /= OutputStateWithTranslations.Norm();
	  char* TmpString = new char [32];
	  sprintf (TmpString, "torus_maxocc_%ld_", Manager.GetInteger("max-occ"));  
	  char* OutputFileName = ReplaceString(Manager.GetString("input-state"), "torus_", TmpString);
	  if (OutputStateWithTranslations.WriteVector(OutputFileName) == false)
	    {
	      cout << "error while writing " << OutputFileName << endl;
	      return -1;
	    }  
	}
    }
  else
    {
      BosonOnTorusShort* Space = 0;
      Space = new BosonOnTorusShort(NbrParticles, NbrFluxQuanta, KxMomentum, KyMomentum);
      ComplexVector InputState;
      if (InputState.ReadVectorTest(Manager.GetString("input-state")) == true)
	{
	  if (InputState.ReadVector(Manager.GetString("input-state")) == false)
	    {
	      cout << "error while reading " << Manager.GetString("input-state") << endl;
	      return -1;
	    }
	  if (InputState.GetVectorDimension() != Space->GetHilbertSpaceDimension())
	    {
	      cout << "error: vector and Hilbert-space have unequal dimensions " << InputState.GetVectorDimension() 
		   << " " << Space->GetHilbertSpaceDimension() << endl;
	      return -1;
	    }
	  unsigned long* TmpState = new unsigned long [NbrFluxQuanta];
	  ComplexVector OutputState (Space->GetLargeHilbertSpaceDimension(), true);
	  double TotalWeight = 0.0;
	  unsigned long MaxOccupation = (unsigned long) Manager.GetInteger("max-occ");
	  for (long i = 0l; i < Space->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      Space->GetOccupationNumber(i, TmpState);
	      int Pos = 0;
	      while ((Pos < NbrFluxQuanta) && (TmpState[Pos] <= MaxOccupation))
		++Pos;
	      if (Pos == NbrFluxQuanta)
		{
		  OutputState[i] = InputState[i];
		  TotalWeight += SqrNorm(InputState[i]);
		}
	    }
	  cout << "total weight with occupation <= " << MaxOccupation << " : " <<  TotalWeight << endl;
	  if (Manager.GetBoolean("discard-vector") == false)
	    {
	      OutputState /= OutputState.Norm();
	      char* TmpString = new char [64];
	      sprintf (TmpString, "torus_kysym_maxocc_%ld_", Manager.GetInteger("max-occ"));  
	      char* OutputFileName = ReplaceString(Manager.GetString("input-state"), "torus_kysym_", TmpString);
	      if (OutputState.WriteVector(OutputFileName) == false)
		{
		  cout << "error while writing " << OutputFileName << endl;
		  return -1;
		}  
	    }
	}
      else
	{
	  RealVector RealInputState;
	  if (RealInputState.GetVectorDimension() != Space->GetHilbertSpaceDimension())
	    {
	      cout << "error: vector and Hilbert-space have unequal dimensions " << RealInputState.GetVectorDimension() 
		   << " " << Space->GetHilbertSpaceDimension() << endl;
	      return -1;
	    }
	  unsigned long* TmpState = new unsigned long [NbrFluxQuanta];
	  RealVector OutputState (Space->GetLargeHilbertSpaceDimension(), true);
	  double TotalWeight = 0.0;
	  unsigned long MaxOccupation = (unsigned long) Manager.GetInteger("max-occ");
	  for (long i = 0l; i < Space->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      Space->GetOccupationNumber(i, TmpState);
	      int Pos = 0;
	      while ((Pos < NbrFluxQuanta) && (TmpState[Pos] <= MaxOccupation))
		++Pos;
	      if (Pos == NbrFluxQuanta)
		{
		  OutputState[i] = RealInputState[i];
		  TotalWeight += RealInputState[i] * RealInputState[i];
		}
	    }
	  cout << "total weight with occupation <= " << MaxOccupation << " : " <<  TotalWeight << endl;
	  if (Manager.GetBoolean("discard-vector") == false)
	    {
	      OutputState /= OutputState.Norm();
	      char* TmpString = new char [64];
	      sprintf (TmpString, "torus_kysym_maxocc_%ld_", Manager.GetInteger("max-occ"));  
	      char* OutputFileName = ReplaceString(Manager.GetString("input-state"), "torus_kysym_", TmpString);
	      if (OutputState.WriteVector(OutputFileName) == false)
		{
		  cout << "error while writing " << OutputFileName << endl;
		  return -1;
		}  
	    }
	}
    }
  return 0;
}
