#include "Vector/ComplexVector.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"


#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"



#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FTIProjectionIntheLowestBand" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");


  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the state to be projected");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name to store the projected state");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIProjectionIntheLowestBand -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int TotalKx = 0;
  int TotalKy = 0;
  int NbrParticles = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  bool Statistics = true;
  int TotalSpin = 0;

  
  if ( Manager.GetString("ground-file") == 0)
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type  FTIProjectionIntheLowestBand -h" << endl;
      return -1;
    }
  if ((Manager.GetString("ground-file") != 0) &&  (IsFile(Manager.GetString("ground-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("ground-file") << endl;
      return -1;
    }
  char* OutputFileName = 0;
  if (Manager.GetString("output-file") == 0)
    {
      OutputFileName = ReplaceString(Manager.GetString("ground-file"), "twoband", "singleband_projected");
      if (OutputFileName == 0)
	{
	  cout << "can't guess output name from " << Manager.GetString("ground-file") << "(should contain twoband)" << endl;
	  return 0;
	}
    }
  else
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }

 if (FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName( Manager.GetString("ground-file"), NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy, TotalSpin, Statistics) == false)
  {
     cout << "error while retrieving system parameters from file name " <<  Manager.GetString("ground-file") << endl;
     return -1;
 }

  if (Statistics == true ) 
  {
 	cout <<" Fermion are not yet supported"<<endl;
	return -1;
  }

  BosonOnSquareLatticeWithSU2SpinMomentumSpace * InitialSpace = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy); 
  BosonOnSquareLatticeMomentumSpace * FinalSpace = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy); 

  ComplexVector InitialState;
  if (InitialState.ReadVector(Manager.GetString("ground-file")) == false)
  {
      cout << "can't open vector file " <<  Manager.GetString("ground-file") << endl;
      return -1;      
  }

  ComplexVector FinalState (FinalSpace->GetHilbertSpaceDimension()); 
  InitialSpace->ProjectIntoTheLowestBand(&InitialState, FinalSpace, &FinalState);
  FinalState.WriteVector(OutputFileName);

  return 0;
}
