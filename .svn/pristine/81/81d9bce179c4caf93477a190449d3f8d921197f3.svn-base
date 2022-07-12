#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairing.h"
#include "HilbertSpace/QuasiholeOnSphereWithSpin.h"

#include "Architecture/ArchitectureManager.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHEOnCylinderFileTools.h"

#include "Vector/RealVector.h"


#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates" , "0.01");
  OptionGroup* SystemGroup  = new OptionGroup("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('\n', "input-state", "name of the file containing the state to convert");
  (*SystemGroup) += new SingleStringOption('\n', "degenerate-states", "name of the file containing a list of states (override input-state)");
  (*SystemGroup) += new SingleStringOption ('\n', "directory", "use a specific directory for the input data instead of the current one");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int LzMax = 0;
  int TotalSz = 0;
  int TotalLz = 0;
  bool Statistics = true;  
  int KValue = 1;
  int RValue = 2;
  bool UseCylinderFlag = false;
  double Ratio = 0.0;
  double Perimeter = 0.0;
  

  if ((Manager.GetString("input-state") == 0) && (Manager.GetString("degenerate-states") == 0))
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity -h" << endl;
      return -1;
    }

  int NbrInputStates = 0;

  if (Manager.GetString("input-state") != 0)
    {
      if (FQHEOnCylinderWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, LzMax, TotalLz, TotalSz, Statistics, Ratio, Perimeter) == false)
	{
	  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, LzMax, TotalLz, TotalSz, Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state")  << endl;
	      return -1;
	    }
	}
      else
	{
	  UseCylinderFlag = true;
	}
       NbrInputStates = 1;
   }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      NbrInputStates = DegenerateFile.GetNbrLines();
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  if (FQHEOnCylinderWithSpinFindSystemInfoFromVectorFileName(DegenerateFile(0, i), NbrParticles, LzMax, TotalLz, TotalSz, Statistics, Ratio, Perimeter) == false)
	    {
	      if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(DegenerateFile(0, i), NbrParticles, LzMax, TotalLz, TotalSz, Statistics) == false)
		{
		  cout << "error while retrieving system parameters from file name " << DegenerateFile(0, i)  << endl;
		  return -1;
		}
	    }
	  else
	    {
	      UseCylinderFlag = true;
	    }
 	}
   }


  char* FilePrefix = new char[512];

  if (UseCylinderFlag == true)
    {
      if (Perimeter > 0.0)	
	{
	  if (Statistics == true)
	    {
	      sprintf (FilePrefix, "fermions_cylinder_perimeter_%.6f", Perimeter);
	    }
	  else
	    {
	      sprintf (FilePrefix, "bosons_cylinder_perimeter_%.6f", Perimeter);
	    }
	}
      else
	{
	  if (Statistics == true)
	    {
	      sprintf (FilePrefix, "fermions_cylinder_ratio_%.6f", Ratio);
	    }
	  else
	    {
	      sprintf (FilePrefix, "bosons_cylinder_ratio_%.6f", Ratio);
	    }      
	}
    }
  else
    {
      if (Statistics == true)
	{
	  sprintf (FilePrefix, "fermions");
	}
      else
	{
	  sprintf (FilePrefix, "bosons");
	}
    }

  RealVector* InputStates = new RealVector[NbrInputStates];
  char** InputStateNames = new char*[NbrInputStates];
  if (Manager.GetString("input-state") != 0)
    {
      InputStateNames[0] = new char [strlen(Manager.GetString("input-state")) + 1];
      strcpy (InputStateNames[0], Manager.GetString("input-state"));
      if (InputStates[0].ReadVector (Manager.GetString("input-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("input-state") << endl;
	  return -1;      
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      if (InputStates[0].ReadVector (DegenerateFile(0, 0)) == false)
	{
	  cout << "can't open vector file " << DegenerateFile(0, 0) << endl;
	  return -1;      
	}	  
      InputStateNames[0] = new char [strlen(DegenerateFile(0, 0)) + 1];
      strcpy (InputStateNames[0], DegenerateFile(0, 0));
      for (int i = 1; i < NbrInputStates; ++i)
	{
	  if (InputStates[i].ReadVector (DegenerateFile(0, i)) == false)
	    {
	      cout << "can't open vector file " << DegenerateFile(0, i) << endl;
	      return -1;      
	    }	  
	  if (InputStates[0].GetVectorDimension() != InputStates[i].GetVectorDimension())
	    {
	      cout << "error, " << DegenerateFile(0, 0) << " and " <<  DegenerateFile(0, i) << " don't have the same  dimension (" 
		   << InputStates[0].GetVectorDimension() << " and " << InputStates[i].GetVectorDimension()<< ")" << endl;
	      return -1;
	    }
	  InputStateNames[i] = new char [strlen(DegenerateFile(0, i)) + 1];
	  strcpy (InputStateNames[i], DegenerateFile(0, i));
	}
    }




  QuasiholeOnSphereWithSpinAndPairing* InputSpace;
  InputSpace = new QuasiholeOnSphereWithSpin (KValue, RValue, TotalLz, LzMax, NbrParticles, TotalSz, Manager.GetString("directory"), FilePrefix);
  QuasiholeOnSphereWithSpinAndPairing* OutputSpace;
  OutputSpace = new QuasiholeOnSphereWithSpinAndPairing (KValue, RValue, TotalLz, LzMax, TotalSz, Manager.GetString("directory"), FilePrefix);
  char* OriginalNbrParticleString = new char[32];
  char* FixedNbrParticleString = new char[32];
  sprintf (OriginalNbrParticleString, "n_%d", NbrParticles);
  sprintf (FixedNbrParticleString, "fixedn_%d_n_0", NbrParticles);
  for (int i = 0; i < NbrInputStates; ++i)
    {
      char* VectorOutputName = ReplaceString(InputStateNames[i], OriginalNbrParticleString, FixedNbrParticleString);
      cout << "input file : " << InputStateNames[i] << endl;
      RealVector TmpVector = OutputSpace->ConvertToNbodyBasis(InputStates[i], InputSpace);
      if (TmpVector.WriteVector(VectorOutputName) == false)
	{
	  cout << "error, can't write vector " << VectorOutputName << endl;
	}
      cout << "output file : " << VectorOutputName << endl;
      delete[] VectorOutputName;     
    }
  delete InputSpace;
  delete OutputSpace;

  return 0;
}
