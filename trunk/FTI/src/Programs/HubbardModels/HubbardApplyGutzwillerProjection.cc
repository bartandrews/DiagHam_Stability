#include "Vector/ComplexVector.h"

#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "Options/Options.h"

#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

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
  OptionManager Manager ("HubbardApplyGutzwillerProjection" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption('i', "input-state", "name of the file containing the state whose Kx momentum has to be computed");
  (*SystemGroup) += new SingleStringOption('\n', "degenerate-states", "name of the file containing a list of states (override input-state)");
  (*SystemGroup) += new BooleanOption ('\n', "unnormalized", "store the Gutzwiller projected states without normalizing them");
  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardApplyGutzwillerProjection -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int NbrSites = 0;
  int* XMomentum = 0;
  int* YMomentum = 0;
  int* SzValue = 0;
  bool TotalSpinConservedFlag;
  int XPeriodicity = 0;
  int YPeriodicity = 0;
  bool Statistics = true;
  bool GutzwillerFlag = false;
  int NbrInputStates = 0;

  if ((Manager.GetString("input-state") == 0) && (Manager.GetString("degenerate-states") == 0))
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type HubbardApplyGutzwillerProjection -h" << endl;
      return -1;
    }

  if (Manager.GetString("input-state") != 0)
    {
      NbrInputStates = 1;
      XMomentum = new int[1];
      YMomentum = new int[1];
      SzValue = new int [1];
      if (IsFile(Manager.GetString("input-state")) == false)
	{
	  cout << "can't open file " << Manager.GetString("input-state") << endl;
	  return -1;
	}
      if ((FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites, XMomentum[0], YMomentum[0], XPeriodicity, YPeriodicity, Statistics, GutzwillerFlag) == false) && (FTIHubbardModelWith1DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites, XMomentum[0], XPeriodicity, Statistics, GutzwillerFlag) == false))
	{
	  cout << "error while retrieving system parameters from file name " <<Manager.GetString("input-state")  << endl;
	  return -1;
	}
      TotalSpinConservedFlag = FTIHubbardModelWithSzFindSystemInfoFromVectorFileName (Manager.GetString("input-state"), NbrParticles, NbrSites, SzValue[0], Statistics, GutzwillerFlag);
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
      XMomentum = new int [NbrInputStates];
      YMomentum = new int [NbrInputStates];
      SzValue = new int [NbrInputStates];
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  if ((FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0,i), NbrParticles, NbrSites, XMomentum[i], YMomentum[i],  XPeriodicity, YPeriodicity, Statistics, GutzwillerFlag) == false) && (FTIHubbardModelWith1DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0,i), NbrParticles, NbrSites, XMomentum[i], XPeriodicity, Statistics, GutzwillerFlag) == false))
	    {
	      cout << "error while retrieving system parameters from file name " << DegenerateFile(0, i) << endl;
	      return -1;
	    }
	  TotalSpinConservedFlag = FTIHubbardModelWithSzFindSystemInfoFromVectorFileName (DegenerateFile(0, i), NbrParticles, NbrSites, SzValue[i], Statistics, GutzwillerFlag);
	}
      
    }
  if (YPeriodicity == 0)
    {
      for (int i = 0; i < NbrInputStates; ++i)
	cout << "Nbr particles=" << NbrParticles << " Nbr sites=" << NbrSites << " Kx=" << XMomentum[i] << " Tx=" << XPeriodicity << endl;
    }
  else
    {
      for (int i = 0; i < NbrInputStates; ++i)
	cout << "Nbr particles=" << NbrParticles << " Nbr sites=" << NbrSites << " Kx=" << XMomentum[i] << " Tx=" << XPeriodicity << " Ky=" << YMomentum[i] << " Ty=" << YPeriodicity << endl;
    }
  ComplexVector* InputStates = new ComplexVector[NbrInputStates];
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


  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation** InputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* [NbrInputStates];
  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (Statistics == true)
	{
	  if (TotalSpinConservedFlag == false)
	    {
	      InputSpace[i] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum[i], XPeriodicity, YMomentum[i], YPeriodicity);
	    }
	  else
	    {
	      InputSpace[i] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, SzValue[i], NbrSites, XMomentum[i], XPeriodicity, YMomentum[i], YPeriodicity, 10000000);
	    }
	}
      else
	{
	  cout << "not available for bosons" << endl;
	  return -1;
	}
    }
  

  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation** OutputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* [NbrInputStates];
  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (Statistics == true)
	{
	  if (TotalSpinConservedFlag == false)
	    {
	      OutputSpace[i] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum[i], XPeriodicity, YMomentum[i], YPeriodicity);
	    }
	  else
	    {
	      OutputSpace[i] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, SzValue[i], NbrSites, XMomentum[i], XPeriodicity, YMomentum[i], YPeriodicity, 10000000);
	    }
	}
      else
	{
	  cout << "not available for bosons" << endl;
	  return -1;
	}
    }
  
  char* UnprojectedString = new char[256];
  char* GutzwillerProjectedString = new char[256];
  sprintf (UnprojectedString, "hubbard");
  if (Manager.GetBoolean("unnormalized") == false)
    {
      sprintf (GutzwillerProjectedString, "hubbard_gutzwiller_projected");
    }
  else
    {
      sprintf (GutzwillerProjectedString, "hubbard_unnormalized_gutzwiller_projected");
    }
 
  
  for (int i = 0; i < NbrInputStates; ++i)
    {
      char* VectorOutputName = ReplaceString(InputStateNames[i], UnprojectedString, GutzwillerProjectedString);
      ComplexVector TmpVector = OutputSpace[i]->ConvertToNbodyBasis(InputStates[i], InputSpace[i]);
      if (Manager.GetBoolean("unnormalized") == false)
	{
	  double TmpNorm = TmpVector.Norm();
	  cout << "norm of projected state " << i << " : " << TmpNorm << endl;      
	  TmpVector /= TmpNorm;
	  TmpVector *= Phase(-Arg(TmpVector[0]));
	}
      if (TmpVector.WriteVector(VectorOutputName) == false)
	{
	  cout << "error, can't write vector " << VectorOutputName << endl;
	}
      delete[] VectorOutputName;
    }
}
