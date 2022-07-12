#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/BosonOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnSpherePTruncated.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnCP2.h"
#include "HilbertSpace/BosonOnCP2TzZ3Symmetry.h"
#include "HilbertSpace/FermionOnCP2.h"
#include "HilbertSpace/BosonOnS2xS2.h"
#include "HilbertSpace/FermionOnS2xS2.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

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
  OptionManager Manager ("FQHESphereUnnormalizeState" , "0.01");
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
  (*SystemGroup) += new BooleanOption  ('\n', "p-truncated", "use a p-truncated basis instead of the full n-body basis");
  (*SystemGroup) += new SingleIntegerOption ('\n', "p-truncation", "p-truncation for the p-truncated basis (if --p-truncated is used)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "boson-truncation", "maximum occupation for a given orbital for the p-truncated basis", 1);
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "tzZ3symmetrized-basis", "use fully symmetrized version of the basis on CP2 (only valid if total tz=0, y = 0, override detection from input vector name)");
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "normalization", "indicates which component should be set to one", 0l);
  (*SystemGroup) += new BooleanOption  ('\n', "symmetry-factor", "do not remove(add) the symmetry factor when (un)normalizing");
  (*SystemGroup) += new BooleanOption  ('\n', "conformal-limit", "indicate that the input state is in the conformal limit basis instead of the unnormalized basis");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "name of the unnormalized vector that will be generated");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the vector into a text file");
  (*OutputGroup) += new BooleanOption ('\n', "txt-separatespin", "for the text output, use the sign convention which separates spins");
  (*OutputGroup) += new BooleanOption ('\n', "normalize", "normalize the state instead of unnormalizing");  
  (*OutputGroup) += new BooleanOption ('\n', "noglobal-normalization", "when normalizing, just restoire the geometry/occupation factors without enforcing the full state to be normalized to one");  
  (*OutputGroup) += new SingleDoubleOption  ('\n', "hide-component", "in the text output, hide state components whose absolute value is lower than a given error (0 if all components have to be shown", 0.0);
  (*OutputGroup) += new SingleDoubleOption  ('\n', "hide-inputcomponent", "in the text output, hide state components whose absolute value is lower than a given error in the original normalized state (0 if all components have to be shown", 0.0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereUnnormalizeState -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int LzMax = 0;
  int LzMax2 = 0;
  int TotalLz = 0;
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  char* OutputTxtFileName = Manager.GetString("txt-output");
  bool SymmetrizedBasis = Manager.GetBoolean("symmetrized-basis");
  bool SU2Flag = false;
  bool S2xS2Flag = false;
  bool CP2Flag = false;
  int TotalTz = 0;
  int TotalY = 0;
  int TotalSz = 0;
  double Error = Manager.GetDouble("hide-component");
  double InputError = Manager.GetDouble("hide-inputcomponent");
  bool SymmetryFactor = !(Manager.GetBoolean("symmetry-factor"));
  
  bool SymFlagTz = false;
  bool SymFlagTzMinusParity = false;
  bool SymFlagTzZ3 = false;
	   
  bool Statistics = true;
  if ((strstr(Manager.GetString("input-state"), "su2") == 0) && (strstr(Manager.GetString("input-state"), "cp2") == 0) 
      && (strstr(Manager.GetString("input-state"), "s2xs2") == 0))
    {
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						       NbrParticles, LzMax, TotalLz, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from " << Manager.GetString("input-state") << endl;
	  return -1;
	}
    }
  else
    {
      if ((strstr(Manager.GetString("input-state"), "cp2") == 0) && (strstr(Manager.GetString("input-state"), "s2xs2") == 0))
	{
	  SU2Flag = true;
	  bool SymFlagLz = false;
	  bool SymFlagSz = false;
	  bool SymFlagLzParity = false;
	  bool SymFlagSzParity = false;
	  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
								   NbrParticles, LzMax, TotalLz, TotalSz, SymFlagLz, SymFlagSz, SymFlagLzParity, SymFlagSzParity, Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from " << Manager.GetString("input-state") << endl;
	      return -1;
	    }
	}
      else
	{
	  if (strstr(Manager.GetString("input-state"), "s2xs2") == 0)
	    {
	      CP2Flag = true;
	      
	      if (FQHEOnCP2FindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, LzMax, TotalTz, TotalY,SymFlagTz, SymFlagTzMinusParity, SymFlagTzZ3, Statistics) == false)
		{
		  cout << "error while retrieving system parameters from " << Manager.GetString("input-state") << endl;
		  return -1;
		}
	    }
	  else
	    {
	      S2xS2Flag = true;	      
	      if (FQHEOnS2xS2FindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, LzMax, LzMax2, TotalTz, TotalY, Statistics) == false)
		{
		  cout << "error while retrieving system parameters from " << Manager.GetString("input-state") << endl;
		  return -1;
		}
	    }
	}
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
	    {
	      sprintf (OutputFileName, "fermions_normalized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalLz);
	    }
	  else
	    {
	      if (CP2Flag == false)
		{
		  if (S2xS2Flag == false)
		    {
		      sprintf (OutputFileName, "fermions_unnormalized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalLz);
		    }
		  else
		    {
		      sprintf (OutputFileName, "fermions_unnormalized_s2xs2_n_%d_2s1_%d_2s2_%d_lz_%d_kz_%d.0.vec", NbrParticles, LzMax, LzMax2, TotalTz, TotalY);
		    }
		}
	      else
		{
		  sprintf(OutputFileName, "fermions_unnormalized_cp2_n_%d_2s_%d_tz_%d_y_%d.0.vec", NbrParticles, LzMax, TotalTz, TotalY);
		}
	    }
	}
      else
	{
	  if (Manager.GetBoolean("normalize"))
	    {	
	      sprintf (OutputFileName, "bosons_normalized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalLz);
	    }
	  else
	    {
	      if (CP2Flag == false)
		{
		  if (S2xS2Flag == false)
		    {
		      sprintf (OutputFileName, "bosons_unnormalized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalLz);
		    }
		  else
		    {
		      sprintf (OutputFileName, "bosons_unnormalized_s2xs2_n_%d_2s1_%d_2s2_%d_lz_%d_kz_%d.0.vec", NbrParticles, LzMax, LzMax2, TotalTz, TotalY);
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("tzZ3symmetrized-basis") == true || SymFlagTzZ3 == true)
		    {
		      sprintf(OutputFileName, "bosons_unnormalized_cp2_tzZ3sym_n_%d_2s_%d_tz_%d_y_%d.0.vec", NbrParticles, LzMax, TotalTz, TotalY);
		    }
		  else
		    {
		      sprintf(OutputFileName, "bosons_unnormalized_cp2_n_%d_2s_%d_tz_%d_y_%d.0.vec", NbrParticles, LzMax, TotalTz, TotalY);
		    }
		}
	    }	      	  
	}
    }

  RealVector OutputState;
  if (OutputState.ReadVector (Manager.GetString("input-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
      return -1;      
    }
  RealVector InputState;
  if (InputError != 0.0)
    {
      InputState.Copy(OutputState);
    }
  ParticleOnSphere* OutputBasis = 0;
  if (Statistics == true)
    {
      if (SU2Flag == false)
	{
	  if (Manager.GetBoolean("huge-basis") == true)
	    {
	      if (Manager.GetString("load-hilbert") == 0)
		{
		  cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
		  return -1;
		}
	      OutputBasis = new FermionOnSphereHaldaneHugeBasis (Manager.GetString("load-hilbert"), Manager.GetInteger("memory"));
	    }
	  else
	    {
	      if (HaldaneBasisFlag == false)
		{
		  if (Manager.GetBoolean("p-truncated") == false)
		    {
		      if (CP2Flag == false)
			{
			  if (S2xS2Flag == false)
			    {
			      OutputBasis = new FermionOnSphere(NbrParticles, TotalLz, LzMax);
			    }
			  else
			    {
			      OutputBasis = new FermionOnS2xS2(NbrParticles, LzMax, LzMax2, TotalTz, TotalY);
			    }
			}
		      else
			{
			  
			  OutputBasis = new FermionOnCP2(NbrParticles, LzMax, TotalTz, TotalY);
			}
		    }
		  else
		    {
		      int* ReferenceState = 0;
		      if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceState) == false)
			return -1;
		      OutputBasis = new FermionOnSpherePTruncated(NbrParticles, TotalLz, LzMax, Manager.GetInteger("p-truncation"), ReferenceState);
		    }
		}
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
	}
      else
	{
	  if (HaldaneBasisFlag == false)
	    OutputBasis = new FermionOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz);
	  else
	    {
	      int** ReferenceStates = 0;
	      int NbrReferenceStates;
	      if (Manager.GetString("reference-file") == 0)
		{
		  cout << "error, a reference file is needed for fermions in Haldane basis" << endl;
		  return 0;
		}
	      if (FQHEGetRootPartitionSU2(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceStates, NbrReferenceStates) == false)
		{
		  cout << "error while parsing " << Manager.GetString("reference-file") << endl;	      
		  return 0;
		}
#ifdef __64_BITS__
	      if (LzMax <= 31)
#else
		if (LzMax <= 15)
#endif
		  {
		    OutputBasis = new FermionOnSphereWithSpinHaldaneBasis(NbrParticles, TotalLz, LzMax, TotalSz, ReferenceStates, NbrReferenceStates);
		  }
		else
		  {
#ifdef __128_BIT_LONGLONG__
		    if (LzMax <= 63)
#else
		      if (LzMax <= 31)
#endif
			{
			  OutputBasis = new FermionOnSphereWithSpinHaldaneBasisLong (NbrParticles, TotalLz, LzMax, TotalSz, ReferenceStates, NbrReferenceStates);
			}
		      else
			{
			  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
			  return 0;
		    }	
		  }
	    }
	}
    }
  else
    { 
      if (SU2Flag == false)
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
		{
		  if (Manager.GetBoolean("p-truncated") == false)
		    {
		      if (CP2Flag == false)
			{
			  if (S2xS2Flag == false)
			    {
			      OutputBasis = new BosonOnSphereShort(NbrParticles, TotalLz, LzMax);
			    }
			  else
			    {
			      OutputBasis = new BosonOnS2xS2(NbrParticles, LzMax, LzMax2, TotalTz, TotalY);
			    }
			}
		      else
			{
			  if (Manager.GetBoolean("tzZ3symmetrized-basis") == true || SymFlagTzZ3 == true)
			    OutputBasis = new BosonOnCP2TzZ3Symmetry(NbrParticles, LzMax, TotalTz, TotalY, SymFlagTzMinusParity);
			  else
			    OutputBasis = new BosonOnCP2(NbrParticles, LzMax, TotalTz, TotalY);
			  
			}
		    }
		  else
		    {
		      int* ReferenceState = 0;
		      if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceState) == false)
			return -1;
		  OutputBasis = new BosonOnSpherePTruncated(NbrParticles, TotalLz, LzMax, Manager.GetInteger("p-truncation"), (int) Manager.GetInteger("boson-truncation"), ReferenceState);
		    }
		}
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
	    OutputBasis = new BosonOnSphereWithSU2Spin(NbrParticles, TotalLz, LzMax, TotalSz);
	}
    }

   if (Manager.GetBoolean("normalize"))
     {
       if (Manager.GetBoolean("conformal-limit") == false)
	 {
	   if (Manager.GetBoolean("noglobal-normalization") == false)
	     {
	       OutputBasis->ConvertFromUnnormalizedMonomial(OutputState, Manager.GetInteger("normalization"), SymmetryFactor);
	     }
	   else
	     {
	       OutputBasis->ConvertFromUnnormalizedMonomialNoGlobalNormalization(OutputState);
	     }
	 }
       else
	 {
	   OutputBasis->ConvertFromConformalLimit(OutputState, Manager.GetInteger("normalization"));
	 }
     }
   else
     {
       OutputBasis->ConvertToUnnormalizedMonomial(OutputState, Manager.GetInteger("normalization"), SymmetryFactor);
     }
  

  cout << OutputBasis->GetLargeHilbertSpaceDimension() << " " << OutputState.GetLargeVectorDimension() << endl;
  if (OutputTxtFileName != 0)
    {
      ofstream File;
      File.open(OutputTxtFileName, ios::binary | ios::out);
      File.precision(14);
      if ((SU2Flag == false) || (Manager.GetBoolean("txt-separatespin") == false))
	{
	  if (Error == 0.0)
	    {
	      if (InputError == 0.0)
		{
		  for (long i = 0; i < OutputBasis->GetLargeHilbertSpaceDimension(); ++i)
		    {
		      File << OutputState[i] << " ";
		      OutputBasis->PrintStateMonomial(File, i) << endl;
		    }
		}
	      else
		{
		  for (long i = 0; i < OutputBasis->GetLargeHilbertSpaceDimension(); ++i)
		    {
		      if (fabs(InputState[i]) > InputError)
			{
			  File << OutputState[i] << " ";
			  OutputBasis->PrintStateMonomial(File, i) << endl;
			}
		    }
		}
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
	}
      else
	{
	  FermionOnSphereWithSpin* TmpBasis = (FermionOnSphereWithSpin*)  OutputBasis;
	  double GlobalSign = TmpBasis->GetSpinSeparationSignFromIndex(0);
 	  if (Error == 0.0)
 	    for (long i = 0; i < OutputBasis->GetLargeHilbertSpaceDimension(); ++i)
 	      {
 		File << (OutputState[i] * GlobalSign * TmpBasis->GetSpinSeparationSignFromIndex(i)) << " ";
 		TmpBasis->PrintStateMonomialSeparatedSpin(File, i) << endl;
 	      }
 	  else
 	    for (long i = 0; i < OutputBasis->GetLargeHilbertSpaceDimension(); ++i)
 	      {
 		if (fabs(OutputState[i]) > Error)
 		  {
 		    File << (OutputState[i] * GlobalSign * TmpBasis->GetSpinSeparationSignFromIndex(i)) << " ";
 		    TmpBasis->PrintStateMonomialSeparatedSpin(File, i) << endl;
 		  }
 	      }
	}
      File.close();
    }
  else
    {
      if (OutputFileName != 0)
	{
	  if (OutputState.WriteVector(OutputFileName) == false)
	    {
	      cout << "error while writing output state " << OutputFileName << endl;
	      return -1;
	    }	  
	}
    }
  return 0;
}

