#include "config.h"

#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereFull.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnSpherePTruncated.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

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
  OptionManager Manager ("FQHESphereConvertHaldaneBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "input-file", "input state file name");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "p-truncated", "use a p-truncated basis instead of the full squeezed basis");
  (*SystemGroup) += new SingleIntegerOption ('\n', "p-truncation", "p-truncation for the p-truncated basis (if --p-truncated is used)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "boson-truncation", "maximum occupation for a given orbital for the p-truncated basis", 1);
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", -1);
  (*SystemGroup) += new BooleanOption  ('r', "reverse", "convert a state from the n-body basis to the Haldane basis");
  (*SystemGroup) += new BooleanOption  ('f', "fermion", "use fermionic statistic (override autodetection from input file name)");
  (*SystemGroup) += new BooleanOption  ('b', "boson", "use bosonic statistics (override autodetection from input file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "rational" , "use rational numbers instead of double precision floating point numbers");
  (*SystemGroup) += new BooleanOption  ('\n', "target-haldane", "target space is also a Haldane squeezed basis instead of the full n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "target-referencefile", "definition of the target reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "target-full", "target space is full Hilbert space (all Lz sectors)");
  (*SystemGroup) += new BooleanOption  ('\n', "target-ptruncated", "target space is also a p-truncated basis instead of the full n-body basis");

  (*SystemGroup) += new SingleIntegerOption ('\n', "target-ptruncation", "p-truncation for the target basis (if --target-ptruncated is used)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support (only available when both the source and target spaces are squeezed basis)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the non-symmetric Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "target-loadhilbert", "load the target Hilbert space description from the indicated file",0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (removing any occurence of haldane_)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereConvertHaldaneBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-file") == 0)
    {
      cout << "error, one input file should be provided. See man page for option syntax or type FQHESphereConvertHaldaneBasis -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("input-file")) == false)
    {
      cout << "can't open file " << Manager.GetString("input-file") << endl;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrFluxQuanta = Manager.GetInteger("nbr-flux"); 
  int TotalLz = 0;
  bool ReverseFlag = Manager.GetBoolean("reverse");
  bool Statistics = true;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("input-file"),
						   NbrParticles, NbrFluxQuanta, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-file") << endl;
      return -1;
    }
  if (Manager.GetInteger("total-lz") >= 0)
    TotalLz = Manager.GetInteger("total-lz"); 
  if ((Manager.GetBoolean("boson") == true) || (Manager.GetBoolean("fermion") == true))
    {
      if (Manager.GetBoolean("boson") == true)
	Statistics = false;
      else
	Statistics = true;
    }
  if (((NbrParticles * NbrFluxQuanta) & 1) != (TotalLz & 1))
    {
      cout << "incompatible values for nbr-particles, nbr-flux and total-lz" << endl;
      return -1;
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
      char* TagPosition = strcasestr(InputFileName, "haldane_");
      if (TagPosition == 0)
	{
	  cout << "no default output name can be built from " << InputFileName << endl;
	  return -1;
	}
      OutputFileName = new char [strlen(InputFileName) - 7];
      strncpy (OutputFileName, InputFileName, TagPosition - InputFileName);
      strcpy (OutputFileName + (TagPosition - InputFileName), TagPosition + 8);
    }

  if ((Manager.GetBoolean("huge-basis") == true) && (Manager.GetBoolean("target-haldane") == true))
    {
      if (Statistics == true)
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false)
	    return -1;
	  FermionOnSphereHaldaneHugeBasis* InitialSpace;
	  if (Manager.GetString("load-hilbert") != 0)
	    {
	      InitialSpace = new FermionOnSphereHaldaneHugeBasis(Manager.GetString("load-hilbert"), Manager.GetInteger("memory"));
	    }
	  else
	    {
	      cout << "initial space initialization requires to use the load-hilbert option in huge mode" << endl;
	      return -1;
	    }
	  FermionOnSphereHaldaneHugeBasis* TargetSpace = 0;
	  int* TargetReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("target-referencefile"), NbrParticles, NbrFluxQuanta, TargetReferenceState) == false)
	    return -1;
	  if (Manager.GetString("target-loadhilbert") != 0)
	    {
	      TargetSpace = new FermionOnSphereHaldaneHugeBasis(Manager.GetString("target-loadhilbert"), Manager.GetInteger("memory"));
	    }
	  else
	    {
	      cout <<"Calculating target space:"<<endl;
	      TargetSpace = new FermionOnSphereHaldaneHugeBasis (NbrParticles, TotalLz, NbrFluxQuanta, /*Manager.GetInteger("file-size")*/ 0, TargetReferenceState, /*((unsigned long) Manager.GetInteger("memory"))*/ 100ul << 20, false, /*Manager.GetInteger("huge-fulldim")*/ 0);
	      //cout << "target space initialization requires to use the target-loadhilbert option in huge mode" << endl;
	    }
	  if (Manager.GetBoolean("rational") == false)
	    {
	      RealVector State;
	      if (State.ReadVector (Manager.GetString("input-file")) == false)
		{
		  cout << "can't open vector file " << Manager.GetString("input-file") << endl;
		  return -1;      
		}
	      if ((ReverseFlag == false) && (InitialSpace->GetLargeHilbertSpaceDimension() != State.GetLargeVectorDimension()))
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;      
		}
	      if ((ReverseFlag == true) && (TargetSpace->GetLargeHilbertSpaceDimension() != State.GetLargeVectorDimension()))
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;      
		}	    
	      RealVector OutputState;
	      if (ReverseFlag == false)
		OutputState = InitialSpace->ConvertToNbodyBasis(State, *TargetSpace);
	      else
		OutputState = InitialSpace->ConvertFromNbodyBasis(State, *TargetSpace);
	      if (OutputState.WriteVector(OutputFileName) == false)
		{
		  cout << "error while writing output state " << OutputFileName << endl;
		  return -1;
		}
	    }
	  else
	    {
	      LongRationalVector State;
	      if (State.ReadVector (Manager.GetString("input-file")) == false)
		{
		  cout << "can't open vector file " << Manager.GetString("input-file") << endl;
		  return -1;      
		}
	      if ((ReverseFlag == false) && (InitialSpace->GetLargeHilbertSpaceDimension() != State.GetLargeVectorDimension()))
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;      
		}
	      if ((ReverseFlag == true) && (TargetSpace->GetLargeHilbertSpaceDimension() != State.GetLargeVectorDimension()))
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;      
		}	    
	      LongRationalVector OutputState;
	      if (ReverseFlag == false)
		OutputState = InitialSpace->ConvertToNbodyBasis(State, *TargetSpace);
	      else
		OutputState = InitialSpace->ConvertFromNbodyBasis(State, *TargetSpace);
	      if (OutputState.WriteVector(OutputFileName) == false)
		{
		  cout << "error while writing output state " << OutputFileName << endl;
		  return -1;
		}
	    }
	  delete InitialSpace;
	}
      else
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false)
	    return -1;
	  BosonOnSphereHaldaneHugeBasisShort* InitialSpace;
	  if (Manager.GetString("load-hilbert") != 0)
	    {
	      InitialSpace = new BosonOnSphereHaldaneHugeBasisShort(Manager.GetString("load-hilbert"), Manager.GetInteger("memory"));
	    }
	  else
	    {
	      cout << "initial space initialization requires to use the load-hilbert option in huge mode" << endl;
	      exit(1);
	    }
	  BosonOnSphereHaldaneHugeBasisShort* TargetSpace = 0;
	  int* TargetReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("target-referencefile"), NbrParticles, NbrFluxQuanta, TargetReferenceState) == false)
	    return -1;
	  if (Manager.GetString("target-loadhilbert") != 0)
	    {
	      TargetSpace = new BosonOnSphereHaldaneHugeBasisShort(Manager.GetString("target-loadhilbert"), Manager.GetInteger("memory"));
	    }
	  else
	    {
	      cout << "target space initialization requires to use the load-hilbert option in huge mode" << endl;
	    }
	  if (Manager.GetBoolean("rational") == false)
	    {
	      RealVector State;
	      if (State.ReadVector (Manager.GetString("input-file")) == false)
		{
		  cout << "can't open vector file " << Manager.GetString("input-file") << endl;
		  return -1;      
		}
	      if ((ReverseFlag == false) && (InitialSpace->GetLargeHilbertSpaceDimension() != State.GetLargeVectorDimension()))
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;      
		}
	      if ((ReverseFlag == true) && (TargetSpace->GetLargeHilbertSpaceDimension() != State.GetLargeVectorDimension()))
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;      
		}	    
	      RealVector OutputState;
	      if (ReverseFlag == false)
		OutputState = InitialSpace->ConvertToNbodyBasis(State, *TargetSpace);
	      else
		OutputState = InitialSpace->ConvertFromNbodyBasis(State, *TargetSpace);
	      if (OutputState.WriteVector(OutputFileName) == false)
		{
		  cout << "error while writing output state " << OutputFileName << endl;
		  return -1;
		}
	    }
	  else
	    {
	      LongRationalVector State;
	      if (State.ReadVector (Manager.GetString("input-file")) == false)
		{
		  cout << "can't open vector file " << Manager.GetString("input-file") << endl;
		  return -1;      
		}
	      if ((ReverseFlag == false) && (InitialSpace->GetHilbertSpaceDimension() != State.GetLargeVectorDimension()))
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;      
		}
	      if ((ReverseFlag == true) && (TargetSpace->GetHilbertSpaceDimension() != State.GetLargeVectorDimension()))
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;      
		}	    
	      LongRationalVector OutputState;
	      if (ReverseFlag == false)
		OutputState = InitialSpace->ConvertToNbodyBasis(State, *TargetSpace);
	      else
		OutputState = InitialSpace->ConvertFromNbodyBasis(State, *TargetSpace);
	      if (OutputState.WriteVector(OutputFileName) == false)
		{
		  cout << "error while writing output state " << OutputFileName << endl;
		  return -1;
		}
	    }
	  delete InitialSpace;
	}
    }
  else
    {
      if (Statistics == true)
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false)
	    return -1;
	  FermionOnSphere* InitialSpace;
	  if (Manager.GetString("load-hilbert") != 0)
	    InitialSpace = new FermionOnSphereHaldaneBasis(Manager.GetString("load-hilbert"));
	  else
            if (Manager.GetBoolean("p-truncated") == false)
	      {
	        InitialSpace = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);
	        if (Manager.GetString("save-hilbert") != 0)
		  {
		    InitialSpace->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		    return 0;
		  }
	      }
             else
              {
                InitialSpace = new FermionOnSpherePTruncated(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
              }
 
	  FermionOnSphere* TargetSpace = 0;
          if (Manager.GetBoolean("target-haldane") == false)
	    {
	      if (Manager.GetBoolean("target-ptruncated") == false)
               { 
                 if (Manager.GetBoolean("target-full") == false)
	           TargetSpace = new FermionOnSphere (NbrParticles, TotalLz, NbrFluxQuanta);
                 else
	           TargetSpace = new FermionOnSphereFull (NbrParticles, NbrFluxQuanta);
               }
	      else
		{
		  int* TargetReferenceState = 0;
		  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, TargetReferenceState) == false)
		    return -1;
		  TargetSpace = new FermionOnSpherePTruncated(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("target-ptruncation"), TargetReferenceState);
		}
	    }
	  else
	    {
	      int* TargetReferenceState = 0;
	      if (FQHEGetRootPartition(Manager.GetString("target-referencefile"), NbrParticles, NbrFluxQuanta, TargetReferenceState) == false)
		return -1;
	      if (Manager.GetString("load-hilbert") != 0)
		TargetSpace = new FermionOnSphereHaldaneBasis(Manager.GetString("target-loadhilbert"));
	      else
		{
		  TargetSpace = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, NbrFluxQuanta, TargetReferenceState);
		}
	    }
	  if (Manager.GetBoolean("rational") == false)
	    {
	      RealVector State;
	      if (State.ReadVector (Manager.GetString("input-file")) == false)
		{
		  cout << "can't open vector file " << Manager.GetString("input-file") << endl;
		  return -1;      
		}
	      if ((ReverseFlag == false) && (InitialSpace->GetHilbertSpaceDimension() != State.GetVectorDimension()))
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;      
		}
	      if ((ReverseFlag == true) && (TargetSpace->GetHilbertSpaceDimension() != State.GetVectorDimension()))
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;      
		}	    
	      RealVector OutputState;
              if (Manager.GetBoolean("p-truncated") == false)
                {  
	          if (ReverseFlag == false)
		    OutputState = ((FermionOnSphereHaldaneBasis*)InitialSpace)->ConvertToNbodyBasis(State, *TargetSpace);
	          else
		    OutputState = ((FermionOnSphereHaldaneBasis*)InitialSpace)->ConvertFromNbodyBasis(State, *TargetSpace);
	          if (OutputState.WriteVector(OutputFileName) == false)
		    {
		      cout << "error while writing output state " << OutputFileName << endl;
		      return -1;
		    }
                }
              else
                {
		    OutputState = ((FermionOnSpherePTruncated*)InitialSpace)->ConvertToHaldaneBasis(State, *((FermionOnSphereHaldaneBasis*) TargetSpace));

  	            if (OutputState.WriteVector(OutputFileName) == false)
		      {
		        cout << "error while writing output state " << OutputFileName << endl;
		        return -1;
		      }
                }
	    }
	  else
	    {
	      LongRationalVector State;
	      if (State.ReadVector (Manager.GetString("input-file")) == false)
		{
		  cout << "can't open vector file " << Manager.GetString("input-file") << endl;
		  return -1;      
		}
	      if ((ReverseFlag == false) && (InitialSpace->GetHilbertSpaceDimension() != State.GetVectorDimension()))
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;      
		}
	      if ((ReverseFlag == true) && (TargetSpace->GetHilbertSpaceDimension() != State.GetVectorDimension()))
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;      
		}	    
	      LongRationalVector OutputState;
	      if (ReverseFlag == false)
		OutputState = ((FermionOnSphereHaldaneBasis*)InitialSpace)->ConvertToNbodyBasis(State, *TargetSpace);
	      else
		OutputState =((FermionOnSphereHaldaneBasis*)InitialSpace)->ConvertFromNbodyBasis(State, *TargetSpace);
	      if (OutputState.WriteVector(OutputFileName) == false)
		{
		  cout << "error while writing output state " << OutputFileName << endl;
		  return -1;
		}
	    }
	  delete InitialSpace;
	}
      else
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false)
	    return -1;
#ifdef  __64_BITS__
	  if ((NbrFluxQuanta + NbrParticles - 1) < 63)
#else
	    if ((NbrFluxQuanta + NbrParticles - 1) < 31)	
#endif
	 {
	   BosonOnSphereHaldaneBasisShort *InitialSpace;
	   if (Manager.GetString("load-hilbert") != 0)
	     InitialSpace = new BosonOnSphereHaldaneBasisShort(Manager.GetString("load-hilbert"));
	   else
	     {
	       InitialSpace = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);	  
	       if (Manager.GetString("save-hilbert") != 0)
		 {
		   InitialSpace->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		   return 0;
		 }
	     }
	   BosonOnSphereShort* TargetSpace;
	   if (Manager.GetBoolean("target-haldane") == false)
	     {
	       if (Manager.GetBoolean("target-ptruncated") == false)
		 {
		   TargetSpace = new BosonOnSphereShort(NbrParticles, TotalLz, NbrFluxQuanta);
		 }
	      else
		{
		  int* TargetReferenceState = 0;
		  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, TargetReferenceState) == false)
		    return -1;
		  TargetSpace = new BosonOnSpherePTruncated(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("target-ptruncation"), (int) Manager.GetInteger("boson-truncation"), TargetReferenceState);
		}
	     }
	   else
	     {
	       int* TargetReferenceState = 0;
	       if (FQHEGetRootPartition(Manager.GetString("target-referencefile"), NbrParticles, NbrFluxQuanta, TargetReferenceState) == false)
		 return -1;
	       if (Manager.GetString("load-hilbert") != 0)
		 TargetSpace = new BosonOnSphereHaldaneBasisShort(Manager.GetString("target-loadhilbert"));
	       else
		 {
		   TargetSpace = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, NbrFluxQuanta, TargetReferenceState);
		 }
	     }
	   if (Manager.GetBoolean("rational") == false)
	     {
	       RealVector State;
	       if (State.ReadVector (Manager.GetString("input-file")) == false)
		 {
		   cout << "can't open vector file " << Manager.GetString("input-file") << endl;
		   return -1;      
		 }
	       if ((ReverseFlag == false) && (InitialSpace->GetHilbertSpaceDimension() != State.GetVectorDimension()))
		 {
		   cout << "dimension mismatch between Hilbert space and input state, have " << State.GetVectorDimension() << ", should be " << InitialSpace->GetHilbertSpaceDimension() << endl;
		   return -1;      
		 }
	       if ((ReverseFlag == true) && (TargetSpace->GetHilbertSpaceDimension() != State.GetVectorDimension()))
		 {
		   cout << "dimension mismatch between Hilbert space and input state, have " << State.GetVectorDimension() << ", should be " << TargetSpace->GetHilbertSpaceDimension() << endl;
		   return -1;      
		 }	    
	       RealVector OutputState;
	       if (ReverseFlag == false)
		 OutputState = InitialSpace->ConvertToNbodyBasis(State, *TargetSpace);
	       else
		 OutputState = InitialSpace->ConvertFromNbodyBasis(State, *TargetSpace);
	       if (OutputState.WriteVector(OutputFileName) == false)
		 {
		   cout << "error while writing output state " << OutputFileName << endl;
		   return -1;
		 }
	     }
	   else
	     {
	       LongRationalVector State;
	       if (State.ReadVector (Manager.GetString("input-file")) == false)
		 {
		   cout << "can't open vector file " << Manager.GetString("input-file") << endl;
		   return -1;      
		 }
	       if ((ReverseFlag == false) && (InitialSpace->GetHilbertSpaceDimension() != State.GetVectorDimension()))
		 {
		   cout << "dimension mismatch between Hilbert space and input state" << endl;
		   return -1;      
		 }
	       if ((ReverseFlag == true) && (TargetSpace->GetHilbertSpaceDimension() != State.GetVectorDimension()))
		 {
		   cout << "dimension mismatch between Hilbert space and input state" << endl;
		   return -1;      
		 }	    
	       LongRationalVector OutputState;
	       if (ReverseFlag == false)
		 OutputState = InitialSpace->ConvertToNbodyBasis(State, *TargetSpace);
	       else
		 OutputState = InitialSpace->ConvertFromNbodyBasis(State, *TargetSpace);
	       if (OutputState.WriteVector(OutputFileName) == false)
		 {
		   cout << "error while writing output state " << OutputFileName << endl;
		   return -1;
		 }
	     }
	   delete InitialSpace;
	 }
	}
}
}

