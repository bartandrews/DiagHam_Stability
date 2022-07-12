#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneLargeBasis.h"

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

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


// initialize the Hilbert space from the fuse data file
//
// inputVectors = reference on the fuse input data file
// rowIndex = index of the row to look into in the input data file
// fileNameIndex = column index of the file name in inputVectors
// rootIndex = column index of the root file description in inputVectors
// loadHilbertIndex  = column index of the hilbert space file  in inputVectors
// lzMax = total number of flux quanta
// return value = pointer to the Hilbert space (0 if an error occured)
ParticleOnSphereWithSpin* InitializeHilbertSpaceFromFuseData(MultiColumnASCIIFile& inputVectors, int rowIndex, int fileNameIndex, int rootIndex, int loadHilbertIndex, int lzMax);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereWithSpinFuseStates" , "0.01");
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
      cout << "see man page for option syntax or type FQHESphereWithSpinFuseStates -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int LzMax = 0;
  int TotalLz = 0;
  int TotalSz = 0;
  int Padding = Manager.GetInteger("padding");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  char* OutputTxtFileName = Manager.GetString("txt-output");
  bool SymmetrizedBasis = Manager.GetBoolean("symmetrized-basis");
  MultiColumnASCIIFile InputVectors;
  if (InputVectors.Parse(Manager.GetString("input-states")) == false)
    {
      InputVectors.DumpErrors(cout) << endl;
      return -1;
    }
  
  
  int* Paddings = 0;
  if (InputVectors(6, 0) != 0)
    {
      Paddings = InputVectors.GetAsIntegerArray(6);
    }
  else
    {
      Paddings = new int [InputVectors.GetNbrLines()];
      for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
	Paddings[i] = Padding;
    }
  int LeftNbrParticles = 0;
  int LeftLzMax = 0;
  int LeftTotalLz = 0;
  int LeftTotalSz = 0;
  bool Statistics = true;
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(InputVectors(0, 0),
							   LeftNbrParticles, LeftLzMax, LeftTotalLz, LeftTotalSz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from left state name " << InputVectors(0, 0) << endl;
      return -1;
    }
  int RightNbrParticles = 0;
  int RightLzMax = 0;
  int RightTotalLz = 0;
  int RightTotalSz = 0;
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(InputVectors(1, 0),
							   RightNbrParticles, RightLzMax, RightTotalLz, RightTotalSz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from left state name " << InputVectors(1, 0) << endl;
      return -1;
    }

  NbrParticles = RightNbrParticles + LeftNbrParticles;
  LzMax = RightLzMax + LeftLzMax + Paddings[0];
  TotalLz = 0;
  TotalSz = LeftTotalSz + RightTotalSz;
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
	sprintf (OutputFileName, "bosons_sphere_su2_fused_n_%d_2s_%d_sz_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalSz, TotalLz);
      else
	sprintf (OutputFileName, "fermions_sphere_su2_fused_n_%d_2s_%d_sz_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalSz, TotalLz);
    }

  ParticleOnSphereWithSpin* OutputBasis = 0;
  if (Statistics == false)
    {
      cout << "error : bosons are not supported" << endl;
      return 0;
    }
  else
    {
      if (HaldaneBasisFlag == false)
	{
#ifdef __64_BITS__
	  if (LzMax <= 31)
#else
	    if (LzMax <= 15)
#endif
	      {
		OutputBasis = new FermionOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz);
	      }
	    else
	      {
#ifdef __128_BIT_LONGLONG__
		if (LzMax <= 63)
#else
		  if (LzMax <= 31)
#endif
		    {
		      OutputBasis = new FermionOnSphereWithSpinLong(NbrParticles, TotalLz, LzMax, TotalSz);
		    }
		  else
		    {
		      cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		      return -1;
		    }	
	      }
	}
      else
	{
	  int** ReferenceStates = 0;
	  int NbrReferenceStates;
	  if (FQHEGetRootPartitionSU2(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceStates, NbrReferenceStates) == false)
	    {
	      cout << "error while parsing " << Manager.GetString("reference-file") << endl;	      
	      return -1;
	    }
	  if (Manager.GetString("load-hilbert") != 0)
	    {
#ifdef __64_BITS__
	      if (LzMax <= 31)
#else
		if (LzMax <= 15)
#endif
		  {
		    OutputBasis = 0;//new FermionOnSphereWithSpinHaldaneBasis(Manager.GetString("load-hilbert"));	  
		  }
		else
		  {
#ifdef __128_BIT_LONGLONG__
		    if (LzMax <= 63)
#else
		      if (LzMax <= 31)
#endif
			{
			  OutputBasis = 0;//new FermionOnSphereWithSpinHaldaneLargeBasis(Manager.GetString("load-hilbert"));	  
			}
		      else
			{
			  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
			  return -1;
			}	
		  }
	    }
	  else
	    {
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
			  return -1;
			}	
		  }
	    }
	}
    }

  RealVector OutputState(OutputBasis->GetLargeHilbertSpaceDimension(), true);

  for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
    {
      ParticleOnSphereWithSpin* LeftBasis = InitializeHilbertSpaceFromFuseData(InputVectors, i, 0, 2, 4, LzMax);
      if (LeftBasis == 0)
	return -1;
      ParticleOnSphere* RightBasis = InitializeHilbertSpaceFromFuseData(InputVectors, i, 1, 3, 5, LzMax);
      if (RightBasis == 0)
	return -1;
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

      OutputBasis->FuseStates(OutputState, LeftVector, RightVector, Paddings[i], LeftBasis, RightBasis, SymmetrizedBasis);
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

// initialize the Hilbert space from the fuse data file
//
// inputVectors = reference on the fuse input data file
// rowIndex = index of the row to look into in the input data file
// fileNameIndex = column index of the file name in inputVectors
// rootIndex = column index of the root file description in inputVectors
// loadHilbertIndex  = column index of the hilbert space file  in inputVectors
// lzMax = total number of flux quanta
// return value = pointer to the Hilbert space (0 if an error occured)

ParticleOnSphereWithSpin* InitializeHilbertSpaceFromFuseData(MultiColumnASCIIFile& inputVectors, int rowIndex, int fileNameIndex, int rootIndex, int loadHilbertIndex, int lzMax)
{
  int leftNbrParticles = 0;
  int leftLzMax = 0;
  int leftTotalLz = 0;
  int leftTotalSz = 0;
  bool Statistics = true;
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(inputVectors(fileNameIndex, rowIndex),
							   leftNbrParticles, leftLzMax, leftTotalLz, leftTotalSz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from left state name " << inputVectors(0, rowIndex) << endl;
      return 0;
    }
  ParticleOnSphereWithSpin* LeftBasis = 0;
  if (Statistics == false)
    {
      cout << "error : bosons are not supported (" << inputVectors(0, rowIndex) << ")" <<endl;
      return 0;
    }
  else
    {
      if ((inputVectors(rootIndex, rowIndex) == 0) || (strcmp("none", inputVectors(rootIndex, rowIndex)) == 0))
	{
#ifdef __64_BITS__
	  if (lzMax <= 31)
#else
	    if (lzMax <= 15)
#endif
	      {
		LeftBasis = new FermionOnSphereWithSpin(leftNbrParticles, leftTotalLz, leftLzMax, leftTotalSz);
	      }
	    else
	      {
#ifdef __128_BIT_LONGLONG__
		if (lzMax <= 63)
#else
		  if (lzMax <= 31)
#endif
		    {
		      LeftBasis = new FermionOnSphereWithSpinLong(leftNbrParticles, leftTotalLz, leftLzMax, leftTotalSz);
		    }
		  else
		    {
		      cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		      return 0;
		    }	
	      }
	}
      else
	{
	  int** LeftReferenceStates = 0;
	  int LeftNbrReferenceStates;
	  if (FQHEGetRootPartitionSU2(inputVectors(rootIndex, rowIndex), leftNbrParticles, leftLzMax, LeftReferenceStates, LeftNbrReferenceStates) == false)
	    return 0;
	  if ((inputVectors(loadHilbertIndex, rowIndex) == 0) || (strcmp("none", inputVectors(loadHilbertIndex, rowIndex)) == 0))
	    {
#ifdef __64_BITS__
	      if (lzMax <= 31)
#else
		if (lzMax <= 15)
#endif
		  {
		    LeftBasis = new FermionOnSphereWithSpinHaldaneBasis(leftNbrParticles, leftTotalLz, leftLzMax, leftTotalSz, LeftReferenceStates, LeftNbrReferenceStates);
		  }
		else
		  {
#ifdef __128_BIT_LONGLONG__
		    if (lzMax <= 63)
#else
		      if (lzMax <= 31)
#endif
			{
			  LeftBasis = new FermionOnSphereWithSpinHaldaneBasisLong (leftNbrParticles, leftTotalLz, leftLzMax, leftTotalSz, LeftReferenceStates, LeftNbrReferenceStates);
			}
		      else
			{
			  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
			  return 0;
			}	
		  }
	    }
	  else
	    {
#ifdef __64_BITS__
	      if (lzMax <= 31)
#else
		if (lzMax <= 15)
#endif
		  {
		    LeftBasis = 0;//new FermionOnSphereWithSpinHaldaneBasis(inputVectors(loadHilbertIndex, rowIndex));	  
		  }
		else
		  {
#ifdef __128_BIT_LONGLONG__
		    if (lzMax <= 63)
#else
		      if (lzMax <= 31)
#endif
			{
			  LeftBasis = 0;//new FermionOnSphereWithSpinHaldaneLargeBasis(inputVectors(loadHilbertIndex, rowIndex));	  
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
  return LeftBasis;
}
