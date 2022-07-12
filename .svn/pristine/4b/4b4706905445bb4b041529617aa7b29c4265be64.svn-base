#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnDisk.h"
#include "HilbertSpace/FermionOnDiskHaldaneBasis.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/BosonOnDiskHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnDiskHaldaneHugeBasisShort.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/FQHEOnDiskFileTools.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHEDiskUnnormalizeState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "file that describes states to fuse");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state of the output state");
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "normalization", "indicates which component should be set to one", 0l);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "shift-orbitals", "shift the angular momentum of all orbitals", 0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "name of the unnormalized vector that will be generated");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the vector into a text file");
  (*OutputGroup) += new BooleanOption ('\n', "normalize", "normalize the state instead of unnormalizing");  
  (*OutputGroup) += new SingleDoubleOption  ('\n', "hide-component", "in the test output, hide state components whose absolute value is lower than a given error (0 if all components have to be shown", 0.0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEDiskUnnormalizeState -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int TotalLz = 0;
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  char* OutputTxtFileName = Manager.GetString("txt-output");
  int ForceMaxMomentum = 0;
  double Error = ((SingleDoubleOption*) Manager["hide-component"])->GetDouble();
	   
  bool Statistics = true;
  if (FQHEOnDiskFindSystemInfoFromFileName(Manager.GetString("input-state"), NbrParticles, ForceMaxMomentum, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from " << Manager.GetString("input-state") << endl;
      return -1;
    }

  int TmpMaxMomentum = 0;
  if (Statistics == false)
    {
      TmpMaxMomentum = TotalLz;
      if ((ForceMaxMomentum >= 0) && (ForceMaxMomentum < TmpMaxMomentum))
	TmpMaxMomentum = ForceMaxMomentum; 
    }
  else
    {
      TmpMaxMomentum = (TotalLz - (((NbrParticles - 1) * (NbrParticles - 2)) / 2));
      if ((ForceMaxMomentum >= 0) && (ForceMaxMomentum < TmpMaxMomentum))
	TmpMaxMomentum = ForceMaxMomentum;
    }     

  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = new char [512];
      if (Statistics == true)
	{
	  if (Manager.GetBoolean("normalize"))	
	    sprintf (OutputFileName, "fermions_disk_normalized_n_%d_lzmax_%d_lz_%d.0.vec", NbrParticles, TmpMaxMomentum, TotalLz);
	  else
	    sprintf (OutputFileName, "fermions_disk_unnormalized_n_%d_lzmax_%d_lz_%d.0.vec", NbrParticles, TmpMaxMomentum, TotalLz);
	}
      else
	{
	  if (Manager.GetBoolean("normalize"))	
	    sprintf (OutputFileName, "bosons_disk_normalized_n_%d_lzmax_%d_lz_%d.0.vec", NbrParticles, TmpMaxMomentum, TotalLz);
	  else
	    sprintf (OutputFileName, "bosons_disk_unnormalized_n_%d_lzmax_%d_lz_%d.0.vec", NbrParticles, TmpMaxMomentum, TotalLz);
	}
    }

  RealVector OutputState;
  if (OutputState.ReadVector (Manager.GetString("input-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
      return -1;      
    }
  
  ParticleOnSphere* OutputBasis = 0;
  if (Statistics == true)
    {
      if (HaldaneBasisFlag == false)
	OutputBasis = new FermionOnDisk(NbrParticles, TotalLz, TmpMaxMomentum);
      else
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, TmpMaxMomentum, ReferenceState) == false)
	    return -1;
	  if (Manager.GetString("load-hilbert") != 0)
	    OutputBasis = new FermionOnDiskHaldaneBasis(Manager.GetString("load-hilbert"));	  
	  else
	    OutputBasis = new FermionOnDiskHaldaneBasis(NbrParticles, TotalLz, TmpMaxMomentum, ReferenceState);
	}
    }
  else
    {
      if (HaldaneBasisFlag == false)
	OutputBasis = new BosonOnDiskShort(NbrParticles, TotalLz, TmpMaxMomentum);
      else
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, TmpMaxMomentum, ReferenceState) == false)
	    return -1;
	  if (Manager.GetBoolean("huge-basis") == true)
	    {
	      if (Manager.GetString("load-hilbert") == 0)
		{
		  cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
		  return -1;
		}
	      OutputBasis = new BosonOnDiskHaldaneHugeBasisShort (Manager.GetString("load-hilbert"), Manager.GetInteger("memory"));
	    }
	  else
	    {
	      if (Manager.GetString("load-hilbert") != 0)
		OutputBasis = new BosonOnDiskHaldaneBasisShort(Manager.GetString("load-hilbert"));	  
	      else
		OutputBasis = new BosonOnDiskHaldaneBasisShort(NbrParticles, TotalLz, TmpMaxMomentum, ReferenceState);	  
	    }
	}
    }

   if (Manager.GetBoolean("normalize"))
     OutputBasis->ShiftedConvertFromUnnormalizedMonomial(OutputState, Manager.GetInteger("shift-orbitals"), Manager.GetInteger("normalization"));
   else
     OutputBasis->ConvertToUnnormalizedMonomial(OutputState, Manager.GetInteger("normalization"));
  

  cout << OutputBasis->GetHilbertSpaceDimension() << " " << OutputState.GetVectorDimension() << endl;
  if (OutputTxtFileName != 0)
    {
      ofstream File;
      File.open(OutputTxtFileName, ios::binary | ios::out);
      File.precision(14);
      if (Error == 0.0)
	for (long i = 0; i < OutputBasis->GetLargeHilbertSpaceDimension(); ++i)
	  {
	    File << OutputState[i] << " ";
	    OutputBasis->PrintStateMonomial(File, i) << endl;
	  }
      else
	for (long i = 0; i < OutputBasis->GetLargeHilbertSpaceDimension(); ++i)
	  {
	    if (fabs(OutputState[i]) > Error)
	      {
		File << OutputState[i] << " ";
		OutputBasis->PrintStateMonomial(File, i) << endl;
	      }
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

