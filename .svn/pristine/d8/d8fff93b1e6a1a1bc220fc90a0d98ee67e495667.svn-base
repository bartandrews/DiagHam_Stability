#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/StringTools.h"

#include "Options/Options.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

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
  OptionManager Manager ("FQHESphereBosonsFilterOccupation" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  //  ArchitectureManager Architecture;

  Manager += SystemGroup;
  //  Architecture.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
//  Manager += ToolsGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "vector file that corresponds to the input state");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-occ", "maximum occupation per obital", 1);
  (*SystemGroup) += new BooleanOption ('\n', "maxocc-percomponent", "for multicomponent systems, fix the maximum occupation per obital and per component (the default behavior only fixes the maximum occupation per obital)");
  (*OutputGroup) += new BooleanOption  ('\n', "discard-vector", "do not store the filtered state");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonsFilterOccupation -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  

  int NbrParticles = 0;
  int LzMax = 0;
  int TotalLz = 0;
  bool Statistics = true;
  int TotalSz = 0;
  bool SU2Symmetry = true;
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
							   NbrParticles, LzMax, TotalLz, TotalSz, Statistics) == false)
    {
      SU2Symmetry = false;
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						       NbrParticles, LzMax, TotalLz, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state") << endl;
	  return -1;
	}
    }

  cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << LzMax << " Lz=" << TotalLz << endl;
  if (Statistics == true)
    {
      cout << "FQHESphereBosonsFilterOccupation is only relevant for bosonic states" << endl;
    }

  RealVector InputState;
  if (InputState.ReadVector(Manager.GetString("input-state")) == false)
    {
      cout << "error while reading " << Manager.GetString("input-state") << endl;
      return -1;
    }

  if (SU2Symmetry == false)
    {
      BosonOnSphereShort* Space = 0;
      if ((Manager.GetBoolean("haldane") == false) && (strstr(Manager.GetString("input-state"), "_haldane_") == 0))
	{
	  Space = new BosonOnSphereShort(NbrParticles, TotalLz, LzMax);
	}
      else
	{
	  int* ReferenceState = 0;
	  if (Manager.GetString("reference-file") == 0)
	    {
	      cout << "error, a reference file is needed" << endl;
	      return 0;
	    }
	  ConfigurationParser ReferenceStateDefinition;
	  if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
	    {
	      ReferenceStateDefinition.DumpErrors(cout) << endl;
	      return 0;
	    }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
	    {
	      cout << "NbrParticles is not defined or as a wrong value" << endl;
	      return 0;
	    }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax < 0))
	    {
	      cout << "LzMax is not defined or as a wrong value" << endl;
	      return 0;
	    }
	  int MaxNbrLz;
	  if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
	    {
	      cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
	      return 0;     
	    }
	  if (MaxNbrLz != (LzMax + 1))
	    {
	      cout << "wrong LzMax value in ReferenceState" << endl;
	      return 0;     
	    }
	  if (Manager.GetString("load-hilbert") != 0)
	    Space = new BosonOnSphereHaldaneBasisShort(Manager.GetString("load-hilbert"));
	  else
	    {
	      Space = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, LzMax, ReferenceState);	  
	    }
	}
      if (InputState.GetVectorDimension() != Space->GetHilbertSpaceDimension())
	{
	  cout << "error: vector and Hilbert-space have unequal dimensions " << InputState.GetVectorDimension() 
	       << " " << Space->GetHilbertSpaceDimension() << endl;
	  return -1;
	}
      unsigned long* TmpState = new unsigned long [LzMax + 1];
      RealVector OutputState (Space->GetLargeHilbertSpaceDimension(), true);
      double TotalWeight = 0.0;
      unsigned long MaxOccupation = (unsigned long) Manager.GetInteger("max-occ");
      for (long i = 0l; i < Space->GetLargeHilbertSpaceDimension(); ++i)
	{
	  Space->GetOccupationNumber(i, TmpState);
	  int Pos = 0;
	  while ((Pos <= LzMax) && (TmpState[Pos] <= MaxOccupation))
	    ++Pos;
	  if (Pos > LzMax)
	    {
	      OutputState[i] = InputState[i];
	      TotalWeight += SqrNorm(InputState[i]);
	    }
	}
      cout << "total weight with occupation <= " << MaxOccupation << " : " <<  TotalWeight << endl;
      if (Manager.GetBoolean("discard-vector") == false)
	{
	  OutputState /= OutputState.Norm();
	  char* TmpString = new char [32];
	  sprintf (TmpString, "bosons_maxocc_%ld_", Manager.GetInteger("max-occ"));  
	  char* OutputFileName = ReplaceString(Manager.GetString("input-state"), "bosons_", TmpString);
	  if (OutputState.WriteVector(OutputFileName) == false)
	    {
	      cout << "error while writing " << OutputFileName << endl;
	      return -1;
	    }  
	}
    }
  else
    {
      BosonOnSphereWithSpin* Space = new BosonOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz);
      if (InputState.GetVectorDimension() != Space->GetHilbertSpaceDimension())
	{
	  cout << "error: vector and Hilbert-space have unequal dimensions " << InputState.GetVectorDimension() 
	       << " " << Space->GetHilbertSpaceDimension() << endl;
	  return -1;
	}
      unsigned long* TmpStateUp = new unsigned long [LzMax + 1];
      unsigned long* TmpStateDown = new unsigned long [LzMax + 1];
      RealVector OutputState (Space->GetLargeHilbertSpaceDimension(), true);
      double TotalWeight = 0.0;
      unsigned long MaxOccupation = (unsigned long) Manager.GetInteger("max-occ");

      if (Manager.GetBoolean("maxocc-percomponent") == true)
	{
	  for (long i = 0l; i < Space->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      Space->GetBosonicDescription(i, TmpStateUp, TmpStateDown);
	      int Pos = 0;
	      while ((Pos <= LzMax) && (TmpStateUp[Pos] <= MaxOccupation) && (TmpStateDown[Pos] <= MaxOccupation))
		++Pos;
	      if (Pos > LzMax)
		{
		  OutputState[i] = InputState[i];
		  TotalWeight += SqrNorm(InputState[i]);
		}
	    }
	}
      else
	{
	  for (long i = 0l; i < Space->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      Space->GetBosonicDescription(i, TmpStateUp, TmpStateDown);
	      int Pos = 0;
	      while ((Pos <= LzMax) && ((TmpStateUp[Pos] + TmpStateDown[Pos]) <= MaxOccupation))
		++Pos;
	      if (Pos > LzMax)
		{
		  OutputState[i] = InputState[i];
		  TotalWeight += SqrNorm(InputState[i]);
		}
	    }
	}
      cout << "total weight with occupation <= " << MaxOccupation << " : " <<  TotalWeight << endl;
      if (Manager.GetBoolean("discard-vector") == false)
	{
	  OutputState /= OutputState.Norm();
	  char* TmpString = new char [64];
	  if (Manager.GetBoolean("maxocc-percomponent") == true)
	    {
	      sprintf (TmpString, "bosons_su2_maxoccup_%ld_maxoccdown_%ld_", Manager.GetInteger("max-occ"), Manager.GetInteger("max-occ"));  
	    }
	  else
	    {
	      sprintf (TmpString, "bosons_su2_maxocc_%ld_", Manager.GetInteger("max-occ"));  
	    }
	  char* OutputFileName = ReplaceString(Manager.GetString("input-state"), "bosons_", TmpString);
	  if (OutputState.WriteVector(OutputFileName) == false)
	    {
	      cout << "error while writing " << OutputFileName << endl;
	      return -1;
	    }  
	}
    }
  return 0;
}
