#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"


#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <cstring>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereConvertSymmetrizedState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "input-file", "input state file name");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new BooleanOption  ('f', "fermion", "use fermionic statistic (override autodetection from input file name)");
  (*SystemGroup) += new BooleanOption  ('b', "boson", "use bosonic statistics (override autodetection from input file name)");
  (*SystemGroup) += new BooleanOption  ('r', "symmetrize", "symmetrize state (instead of unsymmetrizing it)");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (removing any occurence of lzsym_)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereConvertSymmetrizedState -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-file") == 0)
    {
      cout << "error, one input file should be provided. See man page for option syntax or type FQHESphereConvertSymmetrizedState -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("input-file")) == false)
    {
      cout << "can't open file " << Manager.GetString("input-file") << endl;
    }

  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      char* InputFileName = Manager.GetString("input-file"); 
      char* TagPosition = strcasestr(InputFileName, "lzsym_");
      if (TagPosition == 0)
	{
	  cout << "no default output name can be built from " << InputFileName << endl;
	  return -1;
	}
      OutputFileName = new char [strlen(InputFileName) - 5];
      strncpy (OutputFileName, InputFileName, TagPosition - InputFileName);
      strcpy (OutputFileName + (TagPosition - InputFileName), TagPosition + 6);
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrFluxQuanta = Manager.GetInteger("nbr-flux"); 
  bool SymmetrizeFlag = Manager.GetBoolean("symmetrize");
  bool Statistics = true;
  int TotalLz = 0;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("input-file"),
						   NbrParticles, NbrFluxQuanta, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-file") << endl;
      return -1;
    }
  if ((Manager.GetBoolean("boson") == true) || (Manager.GetBoolean("fermion") == true))
    {
      if (Manager.GetBoolean("boson") == true)
	Statistics = false;
      else
	Statistics = true;
    }
  if (((NbrParticles * NbrFluxQuanta) & 1) != (TotalLz != 0))
    {
      cout << "incompatible values for nbr-particles, nbr-flux and total-lz" << endl;
      return -1;
    }

  RealVector State;
  if (State.ReadVector (Manager.GetString("input-file")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-file") << endl;
      return -1;      
    }


  if (Statistics == true)
    {
      if (Manager.GetBoolean("haldane") == false)
	{
	  RealVector OutputState;
#ifdef __64_BITS__
	  if (NbrFluxQuanta <= 63)
#else
	    if (NbrFluxQuanta <= 31)
#endif
	      {
		FermionOnSphereSymmetricBasis InitialSpace(NbrParticles, NbrFluxQuanta);
                cout << "Initial dimension: "<<InitialSpace.GetHilbertSpaceDimension()<<endl;
		FermionOnSphere TargetSpace(NbrParticles, TotalLz, NbrFluxQuanta);
		cout << "Target dimension: "<<TargetSpace.GetHilbertSpaceDimension()<<endl;
		if (SymmetrizeFlag)
		  {
		    if (TargetSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
		      {
			cout << "dimension mismatch between Hilbert space and input state" << endl;
			return -1;
		      }
		    OutputState = InitialSpace.ConvertToSymmetricNbodyBasis(State, TargetSpace);
		  }
		else
		  {
		    if (InitialSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
		      {
			cout << "dimension mismatch between Hilbert space and input state" << endl;
			return -1;
		      }
		    OutputState = InitialSpace.ConvertToNbodyBasis(State, TargetSpace);
		  }
		if (OutputState.WriteVector(OutputFileName) == false)
		  {
		    cout << "error while writing output state " << OutputFileName << endl;
		    return -1;
		  }
	      }
	    else
#ifdef __128_BIT_LONGLONG__
	    if (NbrFluxQuanta <= 126)
#else
	      if (NbrFluxQuanta <= 62)
#endif
		{
		  FermionOnSphereSymmetricBasisLong InitialSpace(NbrParticles, NbrFluxQuanta);
		  FermionOnSphereLong TargetSpace(NbrParticles, TotalLz, NbrFluxQuanta);
		  if (SymmetrizeFlag)
		    {
		      if (TargetSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
			{
			  cout << "dimension mismatch between Hilbert space and input state" << endl;
			  return -1;
			}
		      OutputState = InitialSpace.ConvertToSymmetricNbodyBasis(State, TargetSpace);
		    }
		  else
		    {
		      if (InitialSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
			{
			  cout << "dimension mismatch between Hilbert space and input state" << endl;
			  return -1;
			}
		      OutputState = InitialSpace.ConvertToNbodyBasis(State, TargetSpace);
		    }
		  if (OutputState.WriteVector(OutputFileName) == false)
		    {
		      cout << "error while writing output state " << OutputFileName << endl;
		      return -1;
		    }
		}
	}
      else
	{
	  int* ReferenceState = 0;
	  if (Manager.GetString("reference-file") == 0)
	    {
	      ReferenceState = new int[NbrFluxQuanta + 1];
	      for (int i = 0; i <= NbrFluxQuanta; ++i)
		ReferenceState[i] = 0;
	      if (strcasecmp(Manager.GetString("reference-state"), "laughlin") == 0)
		for (int i = 0; i <= NbrFluxQuanta; i += 3)
		  ReferenceState[i] = 1;
	      else
		if (strcasecmp(Manager.GetString("reference-state"), "pfaffian") == 0)
		  for (int i = 0; i <= NbrFluxQuanta; i += 4)
		    {
		      ReferenceState[i] = 1;
		      ReferenceState[i + 1] = 1;
		    }
		else
		  if (strcasecmp(Manager.GetString("reference-state"), "readrezayi3") == 0)
		    for (int i = 0; i <= NbrFluxQuanta; i += 5)
		      {
			ReferenceState[i] = 1;
			ReferenceState[i + 1] = 1;
			ReferenceState[i + 2] = 1;
		      }
		  else
		    {
		      cout << "unknown reference state " << Manager.GetString("reference-state") << endl;
		      return -1;
		    }
	    }
	  else
	    {
	      ConfigurationParser ReferenceStateDefinition;
	      if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
		{
		  ReferenceStateDefinition.DumpErrors(cout) << endl;
		  return -1;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
		{
		  cout << "NbrParticles is not defined or as a wrong value" << endl;
		  return -1;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", NbrFluxQuanta) == false) || (NbrFluxQuanta <= 0))
		{
		  cout << "LzMax is not defined or as a wrong value" << endl;
		  return -1;
		}
	      int MaxNbrLz;
	      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
		{
		  cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
		  return -1;     
		}
	      if (MaxNbrLz != (NbrFluxQuanta + 1))
		{
		  cout << "wrong LzMax value in ReferenceState" << endl;
		  return -1;     
		}
	    }
	  RealVector OutputState;
#ifdef __64_BITS__
	  if (NbrFluxQuanta <= 63)
#else
	    if (NbrFluxQuanta <= 31)
#endif
	      {
		FermionOnSphereHaldaneSymmetricBasis InitialSpace(NbrParticles, NbrFluxQuanta, ReferenceState);
		cout << "Initial dimension: "<<InitialSpace.GetHilbertSpaceDimension()<<endl;
		FermionOnSphereHaldaneBasis TargetSpace(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);
		cout << "Target dimension: "<<TargetSpace.GetHilbertSpaceDimension()<<endl;
		if (SymmetrizeFlag)
		  {
		    if (TargetSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
		      {
			cout << "dimension mismatch between Hilbert space and input state" << endl;
			return -1;
		      }
		    OutputState = InitialSpace.ConvertToSymmetricHaldaneNbodyBasis(State, TargetSpace);
		  }
		else
		  {
		    if (InitialSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
		      {
			cout << "dimension mismatch between Hilbert space and input state" << endl;
			return -1;
		      }
		    OutputState = InitialSpace.ConvertToHaldaneNbodyBasis(State, TargetSpace);
		  }
		if (OutputState.WriteVector(OutputFileName) == false)
		  {
		    cout << "error while writing output state " << OutputFileName << endl;
		    return -1;
		  }
	      }
	    else
#ifdef __128_BIT_LONGLONG__
	    if (NbrFluxQuanta <= 126)
#else
	      if (NbrFluxQuanta <= 62)
#endif
		{
		  FermionOnSphereHaldaneSymmetricBasisLong InitialSpace(NbrParticles, NbrFluxQuanta, ReferenceState);
		  cout << "Initial dimension: "<<InitialSpace.GetHilbertSpaceDimension()<<endl;
		  FermionOnSphereHaldaneBasisLong TargetSpace(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);
		  cout << "Target dimension: "<<TargetSpace.GetHilbertSpaceDimension()<<endl;
		  if (SymmetrizeFlag)
		    {
		      if (TargetSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
			{
			  cout << "dimension mismatch between Hilbert space and input state" << endl;
			  return -1;
			}
		      OutputState = InitialSpace.ConvertToSymmetricHaldaneNbodyBasis(State, TargetSpace);
		    }
		  else
		    {
		      if (InitialSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
			{
			  cout << "dimension mismatch between Hilbert space and input state" << endl;
			  return -1;
			}
		      OutputState = InitialSpace.ConvertToHaldaneNbodyBasis(State, TargetSpace);
		    }
		  if (OutputState.WriteVector(OutputFileName) == false)
		    {
		      cout << "error while writing output state " << OutputFileName << endl;
		      return -1;
		    }
		}
	}
    }
  else
    {
      RealVector OutputState;
#ifdef  __64_BITS__
      if ((NbrFluxQuanta + NbrParticles - 1) < 63)
#else
        if ((NbrFluxQuanta + NbrParticles - 1) < 31)
#endif
          {
	    BosonOnSphereSymmetricBasisShort InitialSpace(NbrParticles, NbrFluxQuanta);
	    cout << "Initial dimension: "<<InitialSpace.GetHilbertSpaceDimension()<<endl;
	    BosonOnSphereShort TargetSpace(NbrParticles, TotalLz, NbrFluxQuanta);
	    cout << "Target dimension: "<<TargetSpace.GetHilbertSpaceDimension()<<endl;
	    if (SymmetrizeFlag)
	      {
		if (TargetSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
		  {
		    cout << "dimension mismatch between Hilbert space and input state" << endl;
		    return -1;
		  }
                cout << "ConvertToSymmetricNbodyBasis not implemented" << endl;
		//		OutputState = InitialSpace.ConvertToSymmetricNbodyBasis(State, TargetSpace);
	      }
	    else
	      {
		if (InitialSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
		  {
		    cout << "dimension mismatch between Hilbert space and input state" << endl;
		    return -1;
		  }
		OutputState = InitialSpace.ConvertToNbodyBasis(State, TargetSpace);
	      }
	    if (OutputState.WriteVector(OutputFileName) == false)
	      {
		cout << "error while writing output state " << OutputFileName << endl;
		return -1;
	      }
	  }
	else
	  {
            BosonOnSphereSymmetricBasis InitialSpace(NbrParticles, NbrFluxQuanta);
            cout << "Initial dimension: "<<InitialSpace.GetHilbertSpaceDimension()<<endl;
            BosonOnSphere TargetSpace(NbrParticles, TotalLz, NbrFluxQuanta);
	    cout << "Target dimension: "<<TargetSpace.GetHilbertSpaceDimension()<<endl;
	    if (SymmetrizeFlag)
	      {
		if (TargetSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
		  {
		    cout << "dimension mismatch between Hilbert space and input state" << endl;
		    return -1;
		  }
		cout << "ConvertToSymmetricNbodyBasis not implemented" << endl;
		//		OutputState = InitialSpace.ConvertToSymmetricNbodyBasis(State, TargetSpace);
	      }
	    else
	      {
		if (InitialSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
		  {
		    cout << "dimension mismatch between Hilbert space and input state" << endl;
		    return -1;
		  }
		OutputState = InitialSpace.ConvertToNbodyBasis(State, TargetSpace);
	      }
	    if (OutputState.WriteVector(OutputFileName) == false)
	      {
		cout << "error while writing output state " << OutputFileName << endl;
		return -1;
	      }
	  }      
    }
}

