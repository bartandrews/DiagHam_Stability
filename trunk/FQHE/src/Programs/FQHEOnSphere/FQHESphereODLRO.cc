#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string.h>


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereODLRO" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "vector file that describes the state that hasd to be truncated");
  (*SystemGroup) += new SingleStringOption  ('o', "output-state", "vector file that describes the smaller system");
  (*SystemGroup) += new BooleanOption  ('\n', "input-unnormalized", "indicates that the input state is written in the unnormalized basis");
  (*SystemGroup) += new BooleanOption  ('\n', "input-haldane", "use the squeezed basis or the input state");
  (*SystemGroup) += new BooleanOption  ('\n', "output-haldane", "use the squeezed basis or the output state");  
  (*SystemGroup) += new SingleStringOption  ('\n', "input-reference", "use a file as the definition of the reference state of the input state");
  (*SystemGroup) += new SingleStringOption  ('\n', "output-reference", "use a file as the definition of the reference state of the output state");
  (*SystemGroup) += new BooleanOption  ('\n', "inputhuge-basis", "use huge Hilbert space support for the input state");
  (*SystemGroup) += new BooleanOption  ('\n', "outputhuge-basis", "use huge Hilbert space support for the output state");
  (*SystemGroup) += new SingleStringOption  ('\n', "pattern", "pattern that has to be shared between the two n-body states");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "shift-pattern", "shift the pattern away from the pole from a given number of orbitals", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "no-unnormalized", "do not use the unnormalized basis as an intermediate step");
  (*SystemGroup) += new BooleanOption  ('\n', "north-south", " compute the ODLRO between north pole and south pole");

  (*OutputGroup) += new BooleanOption ('\n', "save-truncated", "save the truncated state");
  (*OutputGroup) += new SingleStringOption ('\n', "truncated-name", "output file name used to store the truncated state (default name uses input-state and add odlro_pattern to the interaction name)");

  (*PrecalculationGroup) += new SingleStringOption  ('\n', "inputload-hilbert", "load Hilbert space description from the input file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "outputload-hilbert", "load Hilbert space description from the output file",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereODLRO -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int OutputNbrParticles = 0;
  int OutputLzMax = 0;
  int OutputTotalLz = 0;
  
  int InputNbrParticles = 0;
  int InputLzMax = 0;
  int InputTotalLz = 0;
  bool Statistics = true;
  bool NoUnnormalization = Manager.GetBoolean("no-unnormalized");

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						   InputNbrParticles, InputLzMax, InputTotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from input state name " << Manager.GetString("input-state") << endl;
      return -1;
    }
  if ((Manager.GetBoolean("save-truncated") == false) && (Manager.GetBoolean("north-south") == false) &&
      (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("output-state"),
						    OutputNbrParticles, OutputLzMax, OutputTotalLz, Statistics) == false))
    {
      cout << "error while retrieving system parameters from input state name " << Manager.GetString("input-state") << endl;
      return -1;
    }
  
  int PatternNbrParticles = 0;
  int PatternLzMax = 0;
  int* Pattern = 0;
  if (FQHEGetRootPartition(Manager.GetString("pattern"), PatternNbrParticles, PatternLzMax, Pattern) == false)
    return -1;
  
  if ((Manager.GetBoolean("save-truncated") == true) || (Manager.GetBoolean("north-south") == true))
    {
      OutputNbrParticles = InputNbrParticles - PatternNbrParticles;
      OutputLzMax = InputLzMax - PatternLzMax - 1;
      OutputTotalLz = 0;
      for (int i = 0; i <= PatternLzMax; ++i)
	OutputTotalLz -= Pattern[i] * (2 * i);
      OutputTotalLz += InputTotalLz + (InputNbrParticles * InputLzMax);
      OutputTotalLz -= OutputNbrParticles * (2 * (PatternLzMax + 1));
      OutputTotalLz -= OutputNbrParticles * OutputLzMax;
    }

//   cout << OutputNbrParticles << " " << OutputLzMax << " " << OutputTotalLz << endl;

//   int MaxOutputTotalLz = (OutputLzMax - OutputNbrParticles + 1) * OutputNbrParticles;
//   cout << "MaxOutputTotalLz = " << MaxOutputTotalLz << endl;
//   if ((OutputTotalLz > MaxOutputTotalLz) || (OutputTotalLz < -MaxOutputTotalLz))
//     {
//       cout << "ODLRO=0" << endl;
//       return 0;
//     }

  RealVector InputState;
  if (InputState.ReadVector (Manager.GetString("input-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
      return -1;      
    }

  ParticleOnSphere* InputBasis = 0;
  if (Statistics == false)
    {
      if (Manager.GetBoolean("inputhuge-basis") == true)
	{
	  if (Manager.GetString("inputload-hilbert") == 0)
	    {
	      cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
	      return -1;
	    }
	  InputBasis = new BosonOnSphereHaldaneHugeBasisShort (Manager.GetString("inputload-hilbert"), Manager.GetInteger("memory"));
	}
      else
	{
	  if (Manager.GetBoolean("input-haldane") == true)
	    {
	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(Manager.GetString("input-reference"), InputNbrParticles, InputLzMax, ReferenceState) == false)
		return -1;
	      if (Manager.GetString("inputload-hilbert") != 0)
		InputBasis = new BosonOnSphereHaldaneBasisShort(Manager.GetString("inputload-hilbert"));	  
	      else
		InputBasis = new BosonOnSphereHaldaneBasisShort(InputNbrParticles, InputTotalLz, InputLzMax, ReferenceState);	  
	    }
	  else
	    {
	      InputBasis = new BosonOnSphereShort(InputNbrParticles, InputTotalLz, InputLzMax);
	    }
	}
    }
  else
    {
      if (Manager.GetBoolean("input-haldane") == true)
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("input-reference"), InputNbrParticles, InputLzMax, ReferenceState) == false)
	    return -1;
	  if (Manager.GetString("inputload-hilbert") != 0)
	    InputBasis = new FermionOnSphereHaldaneBasis(Manager.GetString("inputload-hilbert"));	  
	  else
	    InputBasis = new FermionOnSphereHaldaneBasis(InputNbrParticles, InputTotalLz, InputLzMax, ReferenceState);
	}
      else
	{
	  InputBasis = new FermionOnSphere(InputNbrParticles, InputTotalLz, InputLzMax);
	}
    }

  ParticleOnSphere* OutputBasis = 0;
  if (Statistics == false)
    {
      if (Manager.GetBoolean("outputhuge-basis") == true)
	{
	  if (Manager.GetString("outputload-hilbert") == 0)
	    {
	      cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
	      return -1;
	    }
	  OutputBasis = new BosonOnSphereHaldaneHugeBasisShort (Manager.GetString("outputload-hilbert"), Manager.GetInteger("memory"));
	}
      else
	{
	  if (Manager.GetBoolean("output-haldane") == true)
	    {
	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(Manager.GetString("output-reference"), OutputNbrParticles, OutputLzMax, ReferenceState) == false)
		return -1;
	      if (Manager.GetString("outputload-hilbert") != 0)
		OutputBasis = new BosonOnSphereHaldaneBasisShort(Manager.GetString("outputload-hilbert"));	  
	      else
		OutputBasis = new BosonOnSphereHaldaneBasisShort(OutputNbrParticles, OutputTotalLz, OutputLzMax, ReferenceState);	  
	    }
	  else
	    {
	      OutputBasis = new BosonOnSphereShort(OutputNbrParticles, OutputTotalLz, OutputLzMax);	  
	    }
	}
    }
  else
    {
      if (Manager.GetBoolean("output-haldane") == true)
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("output-reference"), OutputNbrParticles, OutputLzMax, ReferenceState) == false)
	    return -1;
	  if (Manager.GetString("outputload-hilbert") != 0)
	    OutputBasis = new FermionOnSphereHaldaneBasis(Manager.GetString("outputload-hilbert"));	  
	  else
	    OutputBasis = new FermionOnSphereHaldaneBasis(OutputNbrParticles, OutputTotalLz, OutputLzMax, ReferenceState);
	}
      else
	{
	  OutputBasis = new FermionOnSphere(OutputNbrParticles, OutputTotalLz, OutputLzMax);
	}
    }


  if ((Manager.GetBoolean("input-unnormalized") == false) && (NoUnnormalization == false))
    {
      InputBasis->ConvertToUnnormalizedMonomial(InputState, -1l);
    }
  RealVector TruncatedState = InputBasis->TruncateStateWithPatternConstraint(InputState, OutputBasis, Pattern, PatternLzMax + 1, Manager.GetInteger("shift-pattern"));
  
  double NorthNorm = TruncatedState.Norm();
  if (NorthNorm > 1e-10)
    {
      if (NoUnnormalization == false)
	OutputBasis->ConvertFromUnnormalizedMonomial(TruncatedState, -1l);
      NorthNorm = TruncatedState.Norm();
      cout << "NorthNorm = " << NorthNorm  << endl;
      if ((Manager.GetBoolean("save-truncated") == false) || (Manager.GetBoolean("north-south") == true))
	TruncatedState /= NorthNorm;
    }
  else
    {
      if (Manager.GetBoolean("north-south") == true)
	{
	  cout << "ODLRO=0" << endl;  
	  return 0;
	}
    }
  RealVector SouthPoleTruncatedState;
  if (Manager.GetBoolean("north-south") == true)
    {
      int* SouthPattern = new int [PatternLzMax + 1];
      for (int i = 0; i <= PatternLzMax; ++i)
	SouthPattern[i] = Pattern[PatternLzMax - i];
      int TmpShift = InputLzMax - PatternLzMax;
      SouthPoleTruncatedState = InputBasis->TruncateStateWithPatternConstraint(InputState, OutputBasis, SouthPattern, PatternLzMax + 1, TmpShift);      
      if (NoUnnormalization == false)
	{
	  double SouthNorm = SouthPoleTruncatedState.Norm();
	  if (SouthNorm > 1e-10)
	    {
	      if (NoUnnormalization == false)
		OutputBasis->ConvertFromUnnormalizedMonomial(SouthPoleTruncatedState, -1l);
	      SouthNorm = SouthPoleTruncatedState.Norm();
	      SouthPoleTruncatedState /= SouthNorm;
	    }
	  else
	    {
	      if (Manager.GetBoolean("north-south") == true)
		{
		  cout << "ODLRO=0" << endl;  
		  return 0;
		}
	    }
	}
    }
  
  if (Manager.GetBoolean("save-truncated") == true)
    {
      if (Manager.GetString("truncated-name") != 0)
	{
	  if (TruncatedState.WriteVector(Manager.GetString("truncated-name")) == false)
	    {
	      cout << "can't write vector " << Manager.GetString("truncated-name") << endl;
	      return -1;
	    }
	}
      else
	{
	  char* OutputName = new char [strlen(Manager.GetString("input-state")) + 8 + ((PatternLzMax + 1) * 2)];
	  long Size = strstr (Manager.GetString("input-state"), "_n_") - Manager.GetString("input-state");
	  strncpy (OutputName, Manager.GetString("input-state"), Size);
	  sprintf (OutputName + Size, "_odlro_");
	  Size += 7;
	  for (int i = 0; i <= PatternLzMax; ++i)
	    {
	      sprintf (OutputName + Size, "%d", Pattern[i]);
	      if (Pattern[i] > 9)
		Size += 2;
	      else
		++Size;
	    }
	  sprintf (OutputName + Size, "_n_%d_2s_%d_lz_%d.0.vec", OutputNbrParticles, OutputLzMax, OutputTotalLz);
	  if (TruncatedState.WriteVector(OutputName) == false)
	    {
	      cout << "can't write vector " << OutputName << endl;
	      return -1;
	    }
	  delete[] OutputName;
	}
      return 0;
    }
  if (Manager.GetBoolean("north-south") == true)
    {
      cout.precision(14); 
      cout << "ODLRO=" << (SouthPoleTruncatedState * TruncatedState) << endl;  
      return 0;
    }

  RealVector OutputState;
  if (OutputState.ReadVector(Manager.GetString("output-state")) == false)
    {
      cout << "can't open " << Manager.GetString("output-state") << endl;
      return -1;
    }

  cout.precision(14); 
  cout << "ODLRO=" << fabs(OutputState * TruncatedState) << " " << NorthNorm << endl;



  return 0;
}


