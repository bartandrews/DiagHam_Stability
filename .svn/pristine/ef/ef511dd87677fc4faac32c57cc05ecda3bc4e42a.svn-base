#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"

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


// core part of the fuse state program
//
// inputFileName = name of the file that describes the states to fuse
// outputBasis = output basis
// outputState = output state (i.e. fused state)
// rationalOutputState = output state (i.e. fused state) in rational mode
// tmpOutputState = temporay state used in add mode
// rationalTmpOutputState = temporay state used in add and rational modes
// rationalFlag = use rational instead of double
// addFlag = use add mode instead of merge mode 
// symmetrizedBasisFlag = output basis is Lz symmetric
// defaultPadding = default padding to use
// return value = 0 if no error occured
int FQHESphereFuseStateCore(char* inputFileName, ParticleOnSphere* outputBasis, RealVector& outputState, LongRationalVector& rationalOutputState, RealVector& tmpOutputState, LongRationalVector& rationalTmpOutputState, bool rationalFlag, bool addFlag, bool symmetrizedBasisFlag, int defaultPadding);

// core part of the fuse multiple state program
//
// inputFileName = name of the file that describes the states to fuse
// outputBasis = output basis
// outputState = output state (i.e. fused state)
// rationalOutputState = output state (i.e. fused state) in rational mode
// tmpOutputState = temporay state used in add mode
// rationalTmpOutputState = temporay state used in add and rational modes
// rationalFlag = use rational instead of double
// addFlag = use add mode instead of merge mode 
// symmetrizedBasisFlag = output basis is Lz symmetric
// return value = 0 if no error occured
int FQHESphereFuseMultipleStateCore(char* inputFileName, ParticleOnSphere* outputBasis, RealVector& outputState, LongRationalVector& rationalOutputState, RealVector& tmpOutputState, LongRationalVector& rationalTmpOutputState, bool rationalFlag, bool addFlag, bool symmetrizedBasisFlag);

// get the target Hilbert space data for the fuse multiple state program
//
// inputFileName = name of the file that describes the states to fuse
// nbrParticles = reference on the number of particles
// totalLz = reference on twice the total Lz value
// statistics = reference on the particle statistics
// return value = 0 if no error occured
int FQHESphereFuseMultipleStateCoreGetTargetHilbertSpace(char* inputFileName, int& nbrParticles, int& lzMax, int& totalLz, bool& statistics);


// get the Hilbert space from file 
// 
// inputFileName = input vector
// referenceFile = name of the reference file for Haldane (squeezed) basis, can set to "none" or 0 if no Haldane basis is requested
// hilbertFile = name of the file where the Hilbert space is stored, can be set to "none" or 0 if it does not exist
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the maximum angular momentum for a single particle
// totalLz = reference on twice the total angular momentum
// return value = pointer to the Hilbert space
ParticleOnSphere* FQHESphereFuseStateGetHilbertSpace(char* inputFileName, char* referenceFile, char* hilbertFile, int& nbrParticles, int& lzMax, int& totalLz);


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
  (*SystemGroup) += new BooleanOption  ('\n', "rational" , "use rational numbers instead of double precision floating point numbers");
  (*SystemGroup) += new BooleanOption  ('\n', "add" , "add the different vector instead of merging them");
  (*SystemGroup) += new SingleStringOption  ('\n', "multiple-add", "use multiple description files, each of then in add mode, merging results of each description file");
  (*SystemGroup) += new SingleStringOption  ('\n', "multiple-fuse", "file describing multiple state fusion");

  (*OutputGroup) += new SingleStringOption ('o', "output-file", "name of the fused vector that will be generated");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the vector into a text file");
  
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFuseStates -h" << endl;
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
  bool Statistics = true;
  int Padding = Manager.GetInteger("padding");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  char* OutputTxtFileName = Manager.GetString("txt-output");
  bool SymmetrizedBasis = Manager.GetBoolean("symmetrized-basis");
  MultiColumnASCIIFile InputVectors;

  if (Manager.GetString("input-states") != 0) 
    {
      if (InputVectors.Parse(Manager.GetString("input-states")) == false)
	{
	  InputVectors.DumpErrors(cout) << endl;
	  return -1;
	}
    }
  else
    {
      if (Manager.GetString("multiple-add") != 0)
	{
	  MultiColumnASCIIFile MultipleAddFiles;
	  if (MultipleAddFiles.Parse(Manager.GetString("multiple-add")) == false)
	    {
	      MultipleAddFiles.DumpErrors(cout) << endl;
	      return -1;
	    }
	  if (InputVectors.Parse(MultipleAddFiles(0, 0)) == false)
	    {
	      InputVectors.DumpErrors(cout) << endl;
	      return -1;
	    }
	}
      else
	{
	  if (Manager.GetString("multiple-fuse") != 0)
	    {
	      MultiColumnASCIIFile MultipleFuseFiles;
	      if (MultipleFuseFiles.Parse(Manager.GetString("multiple-fuse")) == false)
		{
		  MultipleFuseFiles.DumpErrors(cout) << endl;
		  return -1;
		}
	      if (FQHESphereFuseMultipleStateCoreGetTargetHilbertSpace(MultipleFuseFiles(0, 0), NbrParticles, LzMax, TotalLz, Statistics) != 0)
		{
		  return -1;
		}
	      cout << "full system  : nbr particles=" << NbrParticles << "  lzmax=" << LzMax << "  totalLz = " << TotalLz << endl;
	    }
	  else
	    {
	      cout << "no input file, see man page for option syntax or type FQHESphereFuseStates -h" << endl;
	      return -1;
	    }
	}
    }
  
  
  int* Paddings = 0;
  double* Coefficients = 0;
  LongRational* RationalCoefficients = 0;
  if (Manager.GetString("multiple-fuse") == 0)
    {
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
      if (Manager.GetBoolean("rational") == false)
	{
	  if (InputVectors(7, 0) != 0)
	    {
	      Coefficients = InputVectors.GetAsDoubleArray(7);
	    }
	  else
	    {
	      Coefficients = new double [InputVectors.GetNbrLines()];
	      for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
		Coefficients[i] = 1.0;
	    }
	}
      else
	{
	  if (InputVectors(7, 0) != 0)
	    {
	      RationalCoefficients = InputVectors.GetAsLongRationalArray(7);
	    }
	  else
	    {
	      RationalCoefficients = new LongRational [InputVectors.GetNbrLines()];
	      for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
		RationalCoefficients[i] = 1l;
	    }
	}
    }

  int LeftNbrParticles = 0;
  int LeftLzMax = 0;
  int LeftTotalLz = 0;
  int RightNbrParticles = 0;
  int RightLzMax = 0;
  int RightTotalLz = 0;
  if (Manager.GetString("multiple-fuse") == 0)
    {
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(InputVectors(0, 0),LeftNbrParticles, LeftLzMax, LeftTotalLz, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from left state name " << InputVectors(0, 0) << endl;
	  return -1;
	}
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(InputVectors(1, 0), RightNbrParticles, RightLzMax, RightTotalLz, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from left state name " << InputVectors(1, 0) << endl;
	  return -1;
	}
      NbrParticles = RightNbrParticles + LeftNbrParticles;
      LzMax = RightLzMax + LeftLzMax  + Paddings[0] + 1;
      TotalLz = 0;
    }
  
  char* OutputFileName = 0;
  cout << "NbrParticles=" << NbrParticles << " LzMax=" << LzMax << " TotalLz=" << TotalLz << endl;
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
    }
  RealVector OutputState;
  LongRationalVector RationalOutputState;
	
  if (Manager.GetBoolean("rational") == false)
    OutputState = RealVector (OutputBasis->GetLargeHilbertSpaceDimension(), true);
  else
    RationalOutputState = LongRationalVector (OutputBasis->GetLargeHilbertSpaceDimension(),true);
  
  RealVector TmpOutputState;
  LongRationalVector RationalTmpOutputState;
  RealVector TmpOutputState2;
  LongRationalVector RationalTmpOutputState2;

  if ((Manager.GetBoolean("add") == true) || (Manager.GetString("multiple-add") != 0) || (Manager.GetString("multiple-fuse") != 0))
    {
      if (Manager.GetBoolean("rational") == false)
	{
	  TmpOutputState = RealVector (OutputBasis->GetLargeHilbertSpaceDimension(), true);
	}
      else
	{
	  RationalTmpOutputState = LongRationalVector (OutputBasis->GetLargeHilbertSpaceDimension(),true);
	}
    }
  if ((Manager.GetString("multiple-add") != 0) || (Manager.GetString("multiple-fuse") != 0))
    {
      if (Manager.GetBoolean("rational") == false)
	{
	  TmpOutputState2 = RealVector (OutputBasis->GetLargeHilbertSpaceDimension(), true);
	}
      else
	{
	  RationalTmpOutputState2 = LongRationalVector (OutputBasis->GetLargeHilbertSpaceDimension(),true);
	}
    }


  if (Manager.GetString("input-states") != 0)
    {
      if (FQHESphereFuseStateCore(Manager.GetString("input-states"), OutputBasis, OutputState, RationalOutputState, TmpOutputState, RationalTmpOutputState, Manager.GetBoolean("rational"), Manager.GetBoolean("add"), SymmetrizedBasis, Padding) != 0)
	return -1;  
    }
  else
    {
      if ((Manager.GetString("multiple-add") != 0) || (Manager.GetString("multiple-fuse") != 0))
	{
	  MultiColumnASCIIFile MultipleAddFiles;
	  if (Manager.GetString("multiple-add") != 0)
	    {
	      if (MultipleAddFiles.Parse(Manager.GetString("multiple-add")) == false)
		{
		  MultipleAddFiles.DumpErrors(cout) << endl;
		  return -1;
		}
	    }
	  else
	    {
	      if (MultipleAddFiles.Parse(Manager.GetString("multiple-fuse")) == false)
		{
		  MultipleAddFiles.DumpErrors(cout) << endl;
		  return -1;
		}
	    }
	  for (int i = 0; i < MultipleAddFiles.GetNbrLines(); ++i)
	    {
	      if (MultipleAddFiles.GetNbrColumns() == 1)
		{
		  if (Manager.GetString("multiple-add") != 0)
		    {
		      if (FQHESphereFuseStateCore(MultipleAddFiles(0, i), OutputBasis, TmpOutputState, RationalTmpOutputState, TmpOutputState2, RationalTmpOutputState2, Manager.GetBoolean("rational"), Manager.GetBoolean("add"), SymmetrizedBasis, Padding) != 0)
			return -1;  	  
		    }
		  else
		    {
		      if (FQHESphereFuseMultipleStateCore(MultipleAddFiles(0, i), OutputBasis, TmpOutputState, RationalTmpOutputState, TmpOutputState2, RationalTmpOutputState2, Manager.GetBoolean("rational"), Manager.GetBoolean("add"), SymmetrizedBasis) != 0)
			return -1;  	  
		    }
		}
	      else
		{
		  if (strcmp(MultipleAddFiles(1, i), "add") == 0)
		    {

		      if (Manager.GetString("multiple-add") != 0)
			{
 			  if (FQHESphereFuseStateCore(MultipleAddFiles(0, i), OutputBasis, TmpOutputState, RationalTmpOutputState, TmpOutputState2, RationalTmpOutputState2, Manager.GetBoolean("rational"), true, SymmetrizedBasis, Padding) != 0)
			    return -1;  	  
			}
		      else
			{
 			  if (FQHESphereFuseMultipleStateCore(MultipleAddFiles(0, i), OutputBasis, TmpOutputState, RationalTmpOutputState, TmpOutputState2, RationalTmpOutputState2, Manager.GetBoolean("rational"), true, SymmetrizedBasis) != 0)
			    return -1;  	  
			}
		    }
		  else
		    {
		      if (Manager.GetString("multiple-add") != 0)
			{
			  if (FQHESphereFuseStateCore(MultipleAddFiles(0, i), OutputBasis, TmpOutputState, RationalTmpOutputState, TmpOutputState2, RationalTmpOutputState2, Manager.GetBoolean("rational"), false, SymmetrizedBasis, Padding) != 0)
			    return -1;  	  
			}
		      else
			{
 			  if (FQHESphereFuseMultipleStateCore(MultipleAddFiles(0, i), OutputBasis, TmpOutputState, RationalTmpOutputState, TmpOutputState2, RationalTmpOutputState2, Manager.GetBoolean("rational"), false, SymmetrizedBasis) != 0)
			    return -1;  	  
			}
		    }	      		  
		}
	      if (Manager.GetBoolean("rational") == false)
		{
		  for (long j = 0; j < TmpOutputState.GetLargeVectorDimension(); ++j)
		    {
		      if (TmpOutputState[j] != 0.0)
			{
			  OutputState[j] = TmpOutputState[j];
			}
		    }
		  TmpOutputState.ClearVector();
		}
	      else
		{
		  for (long j = 0; j < RationalTmpOutputState.GetLargeVectorDimension(); ++j)
		    {
		      if (RationalTmpOutputState[j] != 0l)
			{
			  RationalOutputState[j] = RationalTmpOutputState[j];
			}
		    }
		  RationalTmpOutputState.ClearVector();
		}	      
	    }
	}
      else
	{
	  if (Manager.GetString("multiple-fuse") != 0)
	    {
	      MultiColumnASCIIFile MultipleFuseFiles;
	      if (MultipleFuseFiles.Parse(Manager.GetString("multiple-fuse")) == false)
		{
		  MultipleFuseFiles.DumpErrors(cout) << endl;
		  return -1;
		}
	    }
	}
    }

  if (Manager.GetBoolean("rational") == false)
    {
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
    }
  else
    {
      if (OutputTxtFileName != 0)
	{
	  ofstream File;
	  File.open(OutputTxtFileName, ios::binary | ios::out);
	  File.precision(14);
	  for (long i = 0; i < OutputBasis->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      File << RationalOutputState[i] << " ";
	      OutputBasis->PrintStateMonomial(File, i) << endl;
	    }
	  File.close();
	}
      if (OutputFileName != 0)
	{
	  if (RationalOutputState.WriteVector(OutputFileName) == false)
	    {
	      cout << "error while writing output state " << OutputFileName << endl;
	      return -1;
	    }	  
	}
    }

  return 0;
}


// core part of the fuse state program
//
// inputFileName = name of the file that describes the states to fuse
// outputBasis = output basis
// outputState = output state (i.e. fused state)
// rationalOutputState = output state (i.e. fused state) in rational mode
// tmpOutputState = temporay state used in add mode
// rationalTmpOutputState = temporay state used in add and rational modes
// rationalFlag = use rational instead of double
// addFlag = use add mode instead of merge mode 
// symmetrizedBasisFlag = output basis is Lz symmetric
// defaultPadding = default padding to use
// return value = 0 if no error occured

int FQHESphereFuseStateCore(char* inputFileName, ParticleOnSphere* outputBasis, RealVector& outputState, LongRationalVector& rationalOutputState, RealVector& tmpOutputState, LongRationalVector& rationalTmpOutputState, bool rationalFlag, bool addFlag, bool symmetrizedBasisFlag, int defaultPadding)
{
  MultiColumnASCIIFile InputVectors;
  if (InputVectors.Parse(inputFileName) == false)
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
	Paddings[i] = defaultPadding;
    }
  
  double* Coefficients = 0;
  LongRational* RationalCoefficients = 0;
  if (rationalFlag  == false)
    {
      if (InputVectors(7, 0) != 0)
	{
	  Coefficients = InputVectors.GetAsDoubleArray(7);
	}
      else
	{
	  Coefficients = new double [InputVectors.GetNbrLines()];
	  for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
	    Coefficients[i] = 1.0;
	}
    }
  else
    {
      if (InputVectors(7, 0) != 0)
	{
	  RationalCoefficients = InputVectors.GetAsLongRationalArray(7);
       }
      else
	{
	  RationalCoefficients = new LongRational [InputVectors.GetNbrLines()];
	  for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
	    RationalCoefficients[i] = 1l;
	}
    }

  for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
    {
      int LeftNbrParticles = 0;
      int LeftLzMax = 0;
      int LeftTotalLz = 0;
      int RightNbrParticles = 0;
      int RightLzMax = 0;
      int RightTotalLz = 0;
      ParticleOnSphere* LeftBasis = FQHESphereFuseStateGetHilbertSpace(InputVectors(0, i), InputVectors(2, i), InputVectors(4, i), LeftNbrParticles, LeftLzMax, LeftTotalLz);
      ParticleOnSphere* RightBasis = FQHESphereFuseStateGetHilbertSpace(InputVectors(1, i), InputVectors(3, i), InputVectors(5, i), RightNbrParticles, RightLzMax, RightTotalLz);
      if (rationalFlag == false)
	{
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
	  

	  cout << "local padding = " <<  Paddings[i] << "  local coeffcient = " << Coefficients[i] << endl;
	  if (addFlag == false)
	    outputBasis->FuseStates(outputState, LeftVector, RightVector, Paddings[i], LeftBasis, RightBasis, symmetrizedBasisFlag, Coefficients[i]);
	  else
	    {
	      outputBasis->FuseStates(tmpOutputState, LeftVector, RightVector, Paddings[i], LeftBasis, RightBasis, false, Coefficients[i]);
	      outputState += tmpOutputState;
	      tmpOutputState.ClearVector();
	    }
	}
      else
	{
	  LongRationalVector LeftVector;
	  if (LeftVector.ReadVector (InputVectors(0, i)) == false)
	    {
	      cout << "can't open vector file " << InputVectors(0, i) << endl;
	      return -1;      
	    }
	  LongRationalVector RightVector;
	  if (RightVector.ReadVector (InputVectors(1, i)) == false)
	    {
	      cout << "can't open vector file " << InputVectors(1, i) << endl;
	      return -1;      
	    }


	  cout << "local padding = " <<  Paddings[i] << "  local coeffcient = " << RationalCoefficients[i] << endl;
	  if (addFlag == false)
	    outputBasis->FuseStates(rationalOutputState, LeftVector, RightVector, Paddings[i], LeftBasis, RightBasis, symmetrizedBasisFlag, RationalCoefficients[i]);
	  else
	    {
	      outputBasis->FuseStates(rationalTmpOutputState, LeftVector, RightVector, Paddings[i], LeftBasis, RightBasis, false, RationalCoefficients[i]);	      
	      rationalOutputState += rationalTmpOutputState;
	      rationalTmpOutputState.ClearVector();
	    }
	}
      delete RightBasis;
      delete LeftBasis;      
    }
  if ((symmetrizedBasisFlag == true) && (addFlag == true))
    {
      if (rationalFlag == false)
	{
	  RealVector TmpState = outputBasis->GetLzSymmetricVector(outputBasis, outputState);
	  for (long j = 0; j < outputState.GetLargeVectorDimension(); ++j)
	    {
	      if (TmpState[j] != 0.0)
		{
		  outputState[j] = TmpState[j];
		}
	    } 	  
	}
      else
	{
	  LongRationalVector TmpState = outputBasis->GetLzSymmetricVector(outputBasis, rationalOutputState);
	  for (long j = 0; j < rationalOutputState.GetLargeVectorDimension(); ++j)
	    {
	      if (TmpState[j] != 0l)
		{
		  rationalOutputState[j] = TmpState[j];
		}
	    }
	}
    }
  return 0;
}

// core part of the fuse multiple state program
//
// inputFileName = name of the file that describes the states to fuse
// outputBasis = output basis
// outputState = output state (i.e. fused state)
// rationalOutputState = output state (i.e. fused state) in rational mode
// tmpOutputState = temporay state used in add mode
// rationalTmpOutputState = temporay state used in add and rational modes
// rationalFlag = use rational instead of double
// addFlag = use add mode instead of merge mode 
// symmetrizedBasisFlag = output basis is Lz symmetric
// defaultPadding = default padding to use
// return value = 0 if no error occured

int FQHESphereFuseMultipleStateCore(char* inputFileName, ParticleOnSphere* outputBasis, RealVector& outputState, LongRationalVector& rationalOutputState, RealVector& tmpOutputState, LongRationalVector& rationalTmpOutputState, bool rationalFlag, bool addFlag, bool symmetrizedBasisFlag)
{
  MultiColumnASCIIFile InputVectors;
  if (InputVectors.Parse(inputFileName) == false)
    {
      InputVectors.DumpErrors(cout) << endl;
      return -1;
    }
  int NbrFusedStates = InputVectors.GetNbrColumns();
  if (((NbrFusedStates % 4) != 0) || (NbrFusedStates < 8))
    {
      cout << "wrong number of columns in " << inputFileName << endl;
      return -1;
    }
  NbrFusedStates /= 4;
  int* Paddings = new int [NbrFusedStates];
  int* FusedNbrParticles = new int [NbrFusedStates];
  int* FusedLzMax = new int [NbrFusedStates];
  int* FusedTotalLz = new int [NbrFusedStates];
  ParticleOnSphere** FusedSpaces = new ParticleOnSphere* [NbrFusedStates];
  for (int j = 0; j < NbrFusedStates; ++j)
    {
      Paddings[j] = 0;
      FusedNbrParticles[j] = 0;
      FusedLzMax[j] = 0;
      FusedTotalLz[j] = 0;
    }
  double* Coefficients = 0;
  LongRational* RationalCoefficients = 0;
  if (rationalFlag  == false)
    {
      Coefficients = InputVectors.GetAsDoubleArray(0);
    }
  else
    {
      RationalCoefficients = InputVectors.GetAsLongRationalArray(0);
    }

  for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
    {
      int TmpNbrFusedStates = 0;
      RealVector* FusedVectors = new RealVector[NbrFusedStates];
      LongRationalVector* RationalFusedVectors = new LongRationalVector[NbrFusedStates];
      for (int j = 0; j < NbrFusedStates; ++j)
	{
	  cout << InputVectors(1 + (4 * j), i) << endl;
	  FusedSpaces[TmpNbrFusedStates] = FQHESphereFuseStateGetHilbertSpace(InputVectors(1 + (4 * j), i), InputVectors(2 + (4 * j), i), InputVectors(3 + (4 * j), i), FusedNbrParticles[TmpNbrFusedStates], FusedLzMax[TmpNbrFusedStates], FusedTotalLz[TmpNbrFusedStates]);
	  if (FusedSpaces[TmpNbrFusedStates] != 0)
	    {
	      if (j != (NbrFusedStates - 1))
		{
		  Paddings[TmpNbrFusedStates] = atoi(InputVectors(4 + (4 * j), i));
		}
	      if (rationalFlag == false)
		{
		  if (FusedVectors[TmpNbrFusedStates].ReadVector (InputVectors(1 + (4 * j), i)) == false)
		    {
		      cout << "can't open vector file " << InputVectors(1 + (4 * j), i) << endl;
		      return -1;      
		    }
		}
	      else
		{
		  if (RationalFusedVectors[TmpNbrFusedStates].ReadVector (InputVectors(1 + (4 * j), i)) == false)
		    {
		      cout << "can't open vector file " << InputVectors(1 + (4 * j), i) << endl;
		      return -1;      
		    }
		}
	      ++TmpNbrFusedStates;
	    }
	}
      if (rationalFlag == false)
	{
	  if (addFlag == false)
	    {
	      outputBasis->FuseMultipleStates(outputState, NbrFusedStates, FusedVectors, Paddings, FusedSpaces, symmetrizedBasisFlag, Coefficients[i]);
	    }
	  else
	    {
	      outputBasis->FuseMultipleStates(tmpOutputState, NbrFusedStates, FusedVectors, Paddings, FusedSpaces, false, Coefficients[i]);
	      outputState += tmpOutputState;
	      tmpOutputState.ClearVector();
	    }
	}
      else
	{
	  if (addFlag == false)
	    {
	      outputBasis->FuseMultipleStates(rationalOutputState, NbrFusedStates, RationalFusedVectors, Paddings, FusedSpaces, symmetrizedBasisFlag, RationalCoefficients[i]);
	    }
	  else
	    {
	      outputBasis->FuseMultipleStates(rationalTmpOutputState, NbrFusedStates, RationalFusedVectors, Paddings, FusedSpaces, false, RationalCoefficients[i]);	      
	      rationalOutputState += rationalTmpOutputState;
	      rationalTmpOutputState.ClearVector();
	    }
	}
      delete[] FusedVectors;
      delete[] RationalFusedVectors;      
      for (int j = 0; j < NbrFusedStates; ++j)
	{
	  delete FusedSpaces[j];
	}
    }
  if ((symmetrizedBasisFlag == true) && (addFlag == true))
    {
      if (rationalFlag == false)
	{
	  RealVector TmpState = outputBasis->GetLzSymmetricVector(outputBasis, outputState);
	  for (long j = 0; j < outputState.GetLargeVectorDimension(); ++j)
	    {
	      if (TmpState[j] != 0.0)
		{
		  outputState[j] = TmpState[j];
		}
	    } 	  
	}
      else
	{
	  LongRationalVector TmpState = outputBasis->GetLzSymmetricVector(outputBasis, rationalOutputState);
	  for (long j = 0; j < rationalOutputState.GetLargeVectorDimension(); ++j)
	    {
	      if (TmpState[j] != 0l)
		{
		  rationalOutputState[j] = TmpState[j];
		}
	    }
	}
    }
  delete[] Paddings;
  delete[] FusedNbrParticles;
  delete[] FusedLzMax;
  delete[] FusedTotalLz;
  delete[] FusedSpaces;
  if (rationalFlag  == false)
    {
      delete[] Coefficients;
    }
  else
    {
      delete[] RationalCoefficients;
    }
  return 0;
}

// get the target Hilbert space data for the fuse multiple state program
//
// inputFileName = name of the file that describes the states to fuse
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the angular momentum 
// totalLz = reference on twice the total Lz value
// statistics = reference on the particle statistics
// return value = 0 if no error occured

int FQHESphereFuseMultipleStateCoreGetTargetHilbertSpace(char* inputFileName, int& nbrParticles, int& lzMax, int& totalLz, bool& statistics)
{
  MultiColumnASCIIFile InputVectors;
  if (InputVectors.Parse(inputFileName) == false)
    {
      InputVectors.DumpErrors(cout) << endl;
      return -1;
    }
  int NbrFusedStates = InputVectors.GetNbrColumns();
  if (((NbrFusedStates % 4) != 0) || (NbrFusedStates < 8))
    {
      cout << "wrong number of columns in " << inputFileName << endl;
      return -1;
    }
  NbrFusedStates /= 4;
  nbrParticles = 0;
  lzMax = 0;
  totalLz = 0;
  for (int j = 0; j < NbrFusedStates; ++j)
    {
      int TmpNbrParticles = 0;
      int TmpLzMax = 0;
      int TmpTotalLz = 0;
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(InputVectors(1 + (4 * j), 0), TmpNbrParticles, TmpLzMax, TmpTotalLz, statistics) == false)
	{
	  cout << "error while retrieving system parameters from state name " << InputVectors(1 + (4 * j), 0) << endl;
	  return -1;
	}
      int LocalPadding = 0;
      if (j != (NbrFusedStates - 1))
	{
	  LocalPadding = atoi(InputVectors(4 + (4 * j), 0));
	}
      else
	{
	  LocalPadding = -1;
	}
      totalLz += ((TmpTotalLz +  (TmpLzMax * TmpNbrParticles)) >> 1) + (lzMax * TmpNbrParticles);
      nbrParticles += TmpNbrParticles;
      lzMax += TmpLzMax + LocalPadding + 1;
    }
  totalLz <<= 1;
  totalLz -= lzMax * nbrParticles;
  return 0;
}

// get the Hilbert space from file 
// 
// inputFileName = input vector
// referenceFile = name of the reference file for Haldane (squeezed) basis, can set to "none" or 0 if no Haldane basis is requested
// hilbertFile = name of the file where the Hilbert space is stored, can be set to "none" or 0 if it does not exist
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the maximum angular momentum for a single particle
// totalLz = reference on twice the total angular momentum
// return value = pointer to the Hilbert space

ParticleOnSphere* FQHESphereFuseStateGetHilbertSpace(char* inputFileName, char* referenceFile, char* hilbertFile, int& nbrParticles, int& lzMax, int& totalLz)
{
  if (strcmp("none" , inputFileName) == 0)
    return 0;
  bool Statistics = true;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(inputFileName, nbrParticles, lzMax, totalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from left state name " << inputFileName << endl;
      return 0;
    }
  ParticleOnSphere* TmpBasis = 0;
  if (Statistics == false)
    {
      if ((referenceFile == 0) || (strcmp("none", referenceFile) == 0))
	TmpBasis = new BosonOnSphereShort(nbrParticles, totalLz, lzMax);
      else
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(referenceFile, nbrParticles, lzMax, ReferenceState) == false)
	    return 0;
	  if ((hilbertFile == 0) || (strcmp("none", hilbertFile) == 0))
	    TmpBasis = new BosonOnSphereHaldaneBasisShort(nbrParticles, totalLz, lzMax, ReferenceState);	  
	  else
	    TmpBasis = new BosonOnSphereHaldaneBasisShort(hilbertFile);
	    }
    }
  else
    {
      if ((referenceFile == 0) || (strcmp("none", referenceFile) == 0))
	TmpBasis = new FermionOnSphere(nbrParticles, totalLz, lzMax);
      else
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(referenceFile, nbrParticles, lzMax, ReferenceState) == false)
	    return 0;
	  if ((hilbertFile == 0) || (strcmp("none", hilbertFile) == 0))
	    TmpBasis = new FermionOnSphereHaldaneBasis(nbrParticles, totalLz, lzMax, ReferenceState);	  
	  else
	    TmpBasis = new FermionOnSphereHaldaneBasis(hilbertFile);
	}
    }
  return TmpBasis;
}

