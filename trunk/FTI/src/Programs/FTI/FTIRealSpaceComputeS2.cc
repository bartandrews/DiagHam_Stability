#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnHoneycombLatticeWithSpinSzSymmetryRealSpacePlaquetteExclusionAnd2DTranslation.h"

#include "Operator/ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator.h"


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
  cout.precision(14);
  OptionManager Manager ("FTIRealSpaceComputeS2" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "name of the file ");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerate-states", "single column file describing a set of degenerate states");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIRealSpaceComputeS2 -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int NbrSites = 0;
  bool Statistics = true;
  bool GutzwillerFlag = false;
  bool ClusterExclusionFlag = false;
  int NbrStates = 0;
  int MomentumFlag = false;
  int KxMomentum = 0;
  int XPeriodicity = 0;
  int KyMomentum = 0;
  int YPeriodicity = 0;
  bool FixedSzFlag = false;
  int TotalSz = 0;
  bool SzSymmetryFlag = false;
  int SzParity = 0;

  if ((Manager.GetString("input-state") == 0) && (Manager.GetString("degenerate-states") == 0))
    {
      cout << "error, a state file should be provided. See man page for option syntax or type FTIRealSpaceComputeS2 -h" << endl;
      return -1;
    }
  if (Manager.GetString("input-state") != 0)
    {
      NbrStates = 1;
      if (IsFile(Manager.GetString("input-state")) == false)
	{
	  cout << "can't open file " << Manager.GetString("input-state") << endl;
	  return -1;
	}
      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites, TotalSz, SzParity,
									   KxMomentum, KyMomentum, XPeriodicity, YPeriodicity, 
									   Statistics, GutzwillerFlag, ClusterExclusionFlag) == false)
	{
	  if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites, TotalSz,
									   KxMomentum, KyMomentum, XPeriodicity, YPeriodicity, 
									   Statistics, GutzwillerFlag, ClusterExclusionFlag) == false)
	    {
	      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites, TotalSz,
									   KxMomentum, KyMomentum, XPeriodicity, YPeriodicity, 
									   Statistics, GutzwillerFlag) == false)
	      {
		if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites,
									       KxMomentum, KyMomentum, XPeriodicity, YPeriodicity, 
									       Statistics, GutzwillerFlag) == false)
		{	  
		  if (FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites, TotalSz, Statistics, GutzwillerFlag) == false)
		  {
		    if (FTIHubbardModelFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites, Statistics, GutzwillerFlag) == 	false)
		      {
			cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state") << endl;
			return -1;
		      }
		  }
		    else
		    {
		      FixedSzFlag = true;
		    }
		}
	      else
		{
		  MomentumFlag = true;
		}
	    }
	  else
	    {
	      FixedSzFlag = true;
	      MomentumFlag = true;
	    }
	  }
	  else
	  {
	    MomentumFlag = true;
	    FixedSzFlag = true;
	    ClusterExclusionFlag = true;
	  }
	}
	else
	{
	  MomentumFlag = true;
	  FixedSzFlag = true;
	  ClusterExclusionFlag = true;
	  SzSymmetryFlag = true;
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
      NbrStates = DegenerateFile.GetNbrLines();
      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), NbrParticles, NbrSites, TotalSz,
									   KxMomentum, KyMomentum, XPeriodicity, YPeriodicity, 
									   Statistics, GutzwillerFlag) == false)
	{
	  if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), NbrParticles, NbrSites,
									       KxMomentum, KyMomentum, XPeriodicity, YPeriodicity, 
									       Statistics, GutzwillerFlag) == false)
	    {	  
	      if (FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), NbrParticles, NbrSites, TotalSz, Statistics, GutzwillerFlag) == false)
		{
		  if (FTIHubbardModelFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), NbrParticles, NbrSites, Statistics, GutzwillerFlag) == false)
		    {
		      cout << "error while retrieving system parameters from file name " << DegenerateFile(0, 0) << endl;
		      return -1;
		    }
		}
	      else
		{
		  FixedSzFlag = true;
		}
	    }
	  else
	    {
	      MomentumFlag = true;
	    }
	}
      else
	{
	  FixedSzFlag = true;
	  MomentumFlag = true;
	}
    }

      
  ComplexVector* InputStates = new ComplexVector[NbrStates];
  char** InputStateNames = new char*[NbrStates];
  if (Manager.GetString("input-state") != 0)
    {
      InputStateNames[0] = new char[strlen(Manager.GetString("input-state")) + 1];
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
      InputStateNames[0] = new char[strlen(DegenerateFile(0, 0)) + 1];
      strcpy (InputStateNames[0], DegenerateFile(0, 0));
      for (int i = 1; i < NbrStates; ++i)
	{
	  InputStateNames[i] = new char[strlen(DegenerateFile(0, i)) + 1];
	  strcpy (InputStateNames[i], DegenerateFile(0, i));
	  if (InputStates[i].ReadVector (DegenerateFile(0, i)) == false)
	    {
	      cout << "can't open vector file " << DegenerateFile(0, i) << endl;
	      return -1;      
	    }	  
	  if (InputStates[0].GetVectorDimension() != InputStates[i].GetVectorDimension())
	    {
	      cout << "error, " << DegenerateFile(0, 0) << " and " <<  DegenerateFile(0, i) << "don't have the same  dimension (" << InputStates[0].GetVectorDimension() << " and " << InputStates[i].GetVectorDimension()<< ")" << endl;
	      return -1;
	    }
	}
    }
  ParticleOnSphereWithSpin* InputSpace = 0;
  if (MomentumFlag == false)
    {
      if (Statistics == true)
	{
	  if (FixedSzFlag == true)
	    {
	      if (GutzwillerFlag == false)
		InputSpace = new FermionOnLatticeWithSpinRealSpace (NbrParticles, TotalSz, NbrSites);
	      else
		InputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, TotalSz, NbrSites);
	    }
	  else
	    {
	      if (GutzwillerFlag == false)
		InputSpace = new FermionOnLatticeWithSpinRealSpace (NbrParticles, NbrSites);
	      else
		InputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
	    }
	}
      else
	{
	  cout << "not available for bosons" << endl;
	  return -1;
	}
    }
  else
    {
      if (FixedSzFlag == true)
	{
	  if (Statistics == true)
	    {
	      if (GutzwillerFlag == false)
	      {
		if (ClusterExclusionFlag == false)
		  InputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, TotalSz, NbrSites, KxMomentum, XPeriodicity, KyMomentum, YPeriodicity);
		else
		{
		  bool MinusParitySector = true;
		  if (SzParity == 1)
		    MinusParitySector = false;
		  InputSpace = new FermionOnHoneycombLatticeWithSpinSzSymmetryRealSpacePlaquetteExclusionAnd2DTranslation (NbrParticles, TotalSz, NbrSites, MinusParitySector, KxMomentum, XPeriodicity, KyMomentum, YPeriodicity);
		}
		
	      }
	      else
		InputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, TotalSz, NbrSites, KxMomentum, XPeriodicity, KyMomentum, YPeriodicity);
	    }
	  else
	    {
	      cout << "not available for bosons" << endl;
	      return -1;
	    }
	}
      else
	{
	  if (Statistics == true)
	    {
	      if (GutzwillerFlag == false)
		InputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, NbrSites, KxMomentum, XPeriodicity, KyMomentum, YPeriodicity);
	      else
		InputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, KxMomentum, XPeriodicity, KyMomentum, YPeriodicity);
	    }
	  else
	    {
	      cout << "not available for bosons" << endl;
	      return -1;
	    }
	}
    }

  if (InputSpace->GetHilbertSpaceDimension() != InputStates[0].GetVectorDimension())
    {
      cout << "error, " << Manager.GetString("input-state")  << " has a wrong dimension (" <<InputStates[0].GetVectorDimension() << ", should be " << InputSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  

  if (MomentumFlag == false)
    {
    }
  else
    {
      ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator OperatorS2 ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) InputSpace, FixedSzFlag);
      for (int i = 0; i < NbrStates; ++i)
	{
	  OperatorMatrixElementOperation OperationS2 (&OperatorS2, InputStates[i], InputStates[i], InputStates[i].GetVectorDimension());
	  OperationS2.ApplyOperation(Architecture.GetArchitecture());
	  cout << InputStateNames[i] << " : " << " <S^2>=" << OperationS2.GetScalar() << " <S>=" << (sqrt(Norm(OperationS2.GetScalar()) + 0.25) - 0.5) << endl;
	}
    }
  return 0;
}

