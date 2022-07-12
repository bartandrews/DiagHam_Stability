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
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereFuseStates" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('i', "input-states", "file that describes states to fuse");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "padding", "number of empty one-body states to insert between two fused Hilbert spaces", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state of the output state");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "name of the fused vector that will be generated");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the vector into a text file");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFuseStates -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int LzMax = 0;
  int TotalLz = 0;
  int Padding = Manager.GetInteger("padding");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  char* OutputTxtFileName = Manager.GetString("txt-output");
  bool SymmetrizedBasis = ((BooleanOption*) Manager["symmetrized-basis"])->GetBoolean();
  MultiColumnASCIIFile InputVectors;
  if (InputVectors.Parse(Manager.GetString("input-states")) == false)
    {
      InputVectors.DumpErrors(cout) << endl;
      return -1;
    }
  
  int LeftNbrParticles = 0;
  int LeftLzMax = 0;
  int LeftTotalLz = 0;
  bool Statistics = true;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(InputVectors(0, 0),
						   LeftNbrParticles, LeftLzMax, LeftTotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from left state name " << InputVectors(0, 0) << endl;
      return -1;
    }
  int RightNbrParticles = 0;
  int RightLzMax = 0;
  int RightTotalLz = 0;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(InputVectors(1, 0),
						   RightNbrParticles, RightLzMax, RightTotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from left state name " << InputVectors(1, 0) << endl;
      return -1;
    }

  NbrParticles = RightNbrParticles + LeftNbrParticles;
  LzMax = RightLzMax + LeftLzMax + Padding;
  TotalLz = 0;
  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = new char [256];
      if (Statistics == false)
	sprintf (OutputFileName, "bosons_fused_n_%d_2s_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalLz);
      else
	sprintf (OutputFileName, "fermions_fused_n_%d_2s_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalLz);
    }

  ParticleOnSphere* OutputBasis = 0;
  if (Statistics == false)
    {
      if (Manager.GetBoolean("huge-basis") == true)
	{
	  if (Manager.GetString("load-hilbert") == 0)
	    {
	      cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
	      return -1;
	    }
	  OutputBasis = new BosonOnSphereHaldaneHugeBasisShort (Manager.GetString("load-hilbert"), Manager.GetInteger("memory"));
	}
      else
	{
	  if (HaldaneBasisFlag == false)
	    OutputBasis = new BosonOnSphereShort(NbrParticles, TotalLz, LzMax);
	  else
	    {
	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceState) == false)
		return -1;
	      if (Manager.GetString("load-hilbert") != 0)
		OutputBasis = new BosonOnSphereHaldaneBasisShort(Manager.GetString("load-hilbert"));	  
	      else
		OutputBasis = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, LzMax, ReferenceState);	  
	    }
	}
    }
  else
    {
      if (HaldaneBasisFlag == false)
	OutputBasis = new FermionOnSphere(NbrParticles, TotalLz, LzMax);
      else
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceState) == false)
	    return -1;
	  if (Manager.GetString("load-hilbert") != 0)
	    OutputBasis = new FermionOnSphereHaldaneBasis(Manager.GetString("load-hilbert"));	  
	  else
	    OutputBasis = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, LzMax, ReferenceState);
	}
    }
  RealVector OutputState(OutputBasis->GetLargeHilbertSpaceDimension(), true);

  for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
    {
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(InputVectors(0, i),
						       LeftNbrParticles, LeftLzMax, LeftTotalLz, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from left state name " << InputVectors(0, i) << endl;
	  return -1;
	}
      ParticleOnSphere* LeftBasis = 0;
      if (Statistics == false)
	{
	  if ((InputVectors(2, i) == 0) || (strcmp("none", InputVectors(2, i)) == 0))
	    LeftBasis = new BosonOnSphereShort(LeftNbrParticles, LeftTotalLz, LeftLzMax);
	  else
	    {
	      int* LeftReferenceState = 0;
	      if (FQHEGetRootPartition(InputVectors(2, i), LeftNbrParticles, LeftLzMax, LeftReferenceState) == false)
		return -1;
	      if ((InputVectors(4, i) == 0) || (strcmp("none", InputVectors(4, i)) == 0))
		LeftBasis = new BosonOnSphereHaldaneBasisShort(LeftNbrParticles, LeftTotalLz, LeftLzMax, LeftReferenceState);	  
	      else
		LeftBasis = new BosonOnSphereHaldaneBasisShort(InputVectors(4, i));
	    }
	}
      else
	{
	  if ((InputVectors(2, i) == 0) || (strcmp("none", InputVectors(2, i)) == 0))
	    LeftBasis = new FermionOnSphere(LeftNbrParticles, LeftTotalLz, LeftLzMax);
	  else
	    {
	      int* LeftReferenceState = 0;
	      if (FQHEGetRootPartition(InputVectors(2, i), LeftNbrParticles, LeftLzMax, LeftReferenceState) == false)
		return -1;
	      if ((InputVectors(4, i) == 0) || (strcmp("none", InputVectors(4, i)) == 0))
		LeftBasis = new FermionOnSphereHaldaneBasis(LeftNbrParticles, LeftTotalLz, LeftLzMax, LeftReferenceState);	  
	      else
		LeftBasis = new FermionOnSphereHaldaneBasis(InputVectors(4, i));
	    }
	}
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(InputVectors(1, i),
						       RightNbrParticles, RightLzMax, RightTotalLz, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from left state name " << InputVectors(1, i) << endl;
	  return -1;
	}
      ParticleOnSphere* RightBasis = 0;
      if (Statistics == false)
	{
	  if ((InputVectors(3, i) == 0) || (strcmp("none", InputVectors(3, i)) == 0))
	    RightBasis = new BosonOnSphereShort(RightNbrParticles, RightTotalLz, RightLzMax);
	  else
	    {
	      int* RightReferenceState = 0;
	      if (FQHEGetRootPartition(InputVectors(3, i), RightNbrParticles, RightLzMax, RightReferenceState) == false)
		return -1;
	      if ((InputVectors(5, i) == 0) || (strcmp("none", InputVectors(5, i)) == 0))
		RightBasis = new BosonOnSphereHaldaneBasisShort(RightNbrParticles, RightTotalLz, RightLzMax, RightReferenceState);	  
	      else
		RightBasis = new BosonOnSphereHaldaneBasisShort(InputVectors(5, i));
	    }
	}
      else
	{
	  if ((InputVectors(3, i) == 0) || (strcmp("none", InputVectors(3, i)) == 0))
	    RightBasis = new FermionOnSphere(RightNbrParticles, RightTotalLz, RightLzMax);
	  else
	    {
	      int* RightReferenceState = 0;
	      if (FQHEGetRootPartition(InputVectors(3, i), RightNbrParticles, RightLzMax, RightReferenceState) == false)
		return -1;
	      if ((InputVectors(5, i) == 0) || (strcmp("none", InputVectors(5, i)) == 0))
		RightBasis = new FermionOnSphereHaldaneBasis(RightNbrParticles, RightTotalLz, RightLzMax, RightReferenceState);	  
	      else
		RightBasis = new FermionOnSphereHaldaneBasis(InputVectors(5, i));
	    }
	}
      RealVector LeftVector;
      if (LeftVector.ReadVector (InputVectors(0, i)) == false)
	{
	  cout << "can't open vector file " << InputVectors(0, i) << endl;
	  return -1;      
	}
      RealVector RightVector;
      if (RightVector.ReadVector (InputVectors(1, i)) == false)
	{
	  cout << "can't open vector file " << InputVectors(1, i) << endl;
	  return -1;      
	}


      OutputBasis->FuseStates(OutputState, LeftVector, RightVector, Padding, LeftBasis, RightBasis, SymmetrizedBasis);
      delete RightBasis;
      delete LeftBasis;      
    }

  if (OutputTxtFileName != 0)
    {
      ofstream File;
      File.open(OutputTxtFileName, ios::binary | ios::out);
      File.precision(14);
      for (long i = 0; i < OutputBasis->GetLargeHilbertSpaceDimension(); ++i)
	{
	  File << OutputState[i] << " ";
	  OutputBasis->PrintStateMonomial(File, i) << endl;
	}
      File.close();
    }
  if (OutputFileName != 0)
    {
      if (OutputState.WriteVector(OutputFileName) == false)
	{
	  cout << "error while writing output state " << OutputFileName << endl;
	  return -1;
	}	  
    }

  return 0;
}
