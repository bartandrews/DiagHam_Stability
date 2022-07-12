#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "Options/Options.h"
#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereMonomialsProductOperation.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"


#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <stdio.h>
#include <string.h>

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;
using std::ifstream;

int main(int argc, char** argv)
{
  cout.precision(14);
  
  OptionManager Manager ("FQHESphereLLLStatesProduct" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  
  ArchitectureManager Architecture;
  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
  Manager += OutputGroup;
	
  (*SystemGroup) += new MultipleStringOption ('\0', "states", "names of the vector files obtained using exact diagonalization");
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the Jack polynomial decomposition into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the Jack polynomial decomposition into a text file");
  // (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  // (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  // (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*OutputGroup) += new BooleanOption ('\n', "normalize", "the output vector will be normalize on the factory");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
	
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonsTimesFermions -h" << endl;
      return -1;
    }
	
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
	
  int NbrVectors;
  char** VectorFiles = Manager.GetStrings("states", NbrVectors);
  char* OutputFileName =  Manager.GetString("bin-output");
  char* OutputTxtFileName = ((SingleStringOption*) Manager["txt-output"])->GetString();
  bool FermionFlag = true;
  bool NormalizeFlag = Manager.GetBoolean("normalize");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");

  if (NbrVectors == 1)
    {
      int NbrParticles = 0;
      int LzMax = 0;
      int TotalLz = 0;

      if (NbrParticles == 0)
	if (FQHEOnSphereFindSystemInfoFromVectorFileName(VectorFiles[0], NbrParticles, LzMax, TotalLz, FermionFlag) == false)
	  {
	    return -1;
	  }
      
      int Parity = TotalLz & 1;
      if (Parity != ((NbrParticles * LzMax) & 1))
	{
	  cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
	  return -1;           
	}
      
      if (IsFile(VectorFiles[0]) == false)
	{
	  cout << "state " << VectorFiles[0] << " does not exist or can't be opened" << endl;
	  return -1;
	}
      
      RealVector GroundState;
      if (GroundState.ReadVector (VectorFiles[0]) == false)
	{
	  cout << "can't open vector file " << VectorFiles[0] << endl;
	  return -1;      
	}
      
      ParticleOnSphere * Space = 0;
      
      if (FermionFlag == false)
	{
	  if (HaldaneBasisFlag == false)
	    {
#ifdef  __64_BITS__
	      if ((LzMax + NbrParticles - 1) < 63)
#else
		if ((LzMax + NbrParticles - 1) < 31)	
#endif
		  {
		    Space = new BosonOnSphereShort (NbrParticles, TotalLz, LzMax);
		  }
		else
		  {
		    cout << "error, the Space needed required class BosonOnSphere" <<endl;
		  }
	    }
	  else
	    {
	      int* ReferenceState = 0;
	      if (((SingleStringOption*) Manager["reference-file"])->GetString() == 0)
		{
		  cout << "error, a reference file is needed for bosons in Haldane basis" << endl;
		  return -1;
		}
	      ConfigurationParser ReferenceStateDefinition;
	      if (ReferenceStateDefinition.Parse(((SingleStringOption*) Manager["reference-file"])->GetString()) == false)
		{
		  ReferenceStateDefinition.DumpErrors(cout) << endl;
		  return -1;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
		{
		  cout << "NbrParticles is not defined or as a wrong value" << endl;
		  return -1;
		}
	      int MaxNbrLz;
	      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
		{
		  cout << "error while parsing ReferenceState in " << ((SingleStringOption*) Manager["reference-file"])->GetString() << endl;
		  return -1;     
		}
	      if (MaxNbrLz != (LzMax + 1))
		{
		  cout << "wrong LzMax value in ReferenceState" << endl;
			cout <<MaxNbrLz <<" " <<LzMax<<endl;
		  return -1;
		}
#ifdef  __64_BITS__
	      if (LzMax  < 63)
#else
		if (LzMax  < 31)	
#endif
		  Space = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, LzMax, ReferenceState);
	    }
	}
      else
	{
	  if (HaldaneBasisFlag == false)
	    {
#ifdef __64_BITS__
	      if (LzMax <= 63)
		Space = new FermionOnSphere(NbrParticles, TotalLz, LzMax);
	      else
		{
		  //Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, LzMax);
		  cout <<"This type of Hilbert space is not supported yet"<<endl;
		  return -1;
		}
#else
	      if (LzMax <= 31)
		Space = new FermionOnSphere(NbrParticles, TotalLz, LzMax);
	      else
		{
		  //Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, LzMax);
		  cout <<"This type of Hilbert space is not supported yet"<<endl;
		  return -1;
		}
#endif
	    }
	  else
	    {
	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceState) == false)
		return -1;
	      Space = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, LzMax, ReferenceState);
	    }
	}
      
      ParticleOnSphere * FinalSpace = 0;
      
      FinalSpace = new BosonOnSphereShort(NbrParticles, 2*TotalLz, 2*LzMax);
      
      
      
      if (Space->GetHilbertSpaceDimension() != GroundState.GetVectorDimension())
	{
	  cout << "Number of rows of the vector is not equal to the Hilbert space dimension!";
	  return -1;
	}
      
      
      RealVector OutputVector(FinalSpace->GetHilbertSpaceDimension(),true);	
      MonomialsProductOperation Operation(Space,FinalSpace, &GroundState, &GroundState,&OutputVector,NormalizeFlag);
      Operation.ApplyOperation(Architecture.GetArchitecture());
      OutputVector.WriteVector(OutputFileName);
      
      if(OutputTxtFileName != 0)
	{
	  ofstream File;
	  File.open(OutputTxtFileName, ios::binary | ios::out);
	  File.precision(14);
	  for (long i = 0; i < FinalSpace->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      File << OutputVector[i] << " ";
	      FinalSpace->PrintStateMonomial(File, i) << endl;
	    }
	  File.close();
	}
      return 0;
    }
  
  if (NbrVectors != 2)
    {
      cout << "One or two vector files are required!"<<endl;
      return -1;
    }
  
  int NbrParticles[2] = {0,0};
  int TotalLz[2] = {0,0};
  int LzMax[2] = {0,0};
  
  bool FermionFlag1 = true;
  bool FermionFlag2 = true;
		
  ParticleOnSphere * OutputBasis1 = 0;
  ParticleOnSphere * OutputBasis2 = 0;
  ParticleOnSphere * FinalSpace = 0;
  
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(VectorFiles[0], NbrParticles[0], LzMax[0], TotalLz[0], FermionFlag1) == false)
    {
      return -1;
    }
  
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(VectorFiles[1], NbrParticles[1], LzMax[1], TotalLz[1], FermionFlag2) == false)
    {
      return -1;
    }
  
  if(NbrParticles[1] != NbrParticles[0])
    {
      cout << "NbrParticles in the two states must be the same"<<endl;
      return -1;
    }
  
  bool InverseFlag = false;
  
  
  if ((FermionFlag2 == false) && (FermionFlag1 == true))
    {
      OutputBasis1 = new BosonOnSphereShort(NbrParticles[1], TotalLz[1], LzMax[1]);
      OutputBasis2 = new FermionOnSphere (NbrParticles[0],TotalLz[0], LzMax[0]);
      FinalSpace = new FermionOnSphere(NbrParticles[0],TotalLz[0]+TotalLz[1],LzMax[0]+LzMax[1]);
      InverseFlag = true;
    }
  else
    {
      if(FermionFlag1 == true)
	{
	  OutputBasis1 = new FermionOnSphere (NbrParticles[0], TotalLz[0], LzMax[0]);
	  OutputBasis2 = new FermionOnSphere (NbrParticles[1], TotalLz[1], LzMax[1]);
	  FinalSpace = new BosonOnSphereShort(NbrParticles[0], TotalLz[0] + TotalLz[1], LzMax[0] + LzMax[1]);
	}
      else
	{
	  OutputBasis1 = new BosonOnSphereShort (NbrParticles[0],TotalLz[0], LzMax[0]);
	  
	  if(FermionFlag2 == false)
	    {
	      OutputBasis2 = new BosonOnSphereShort(NbrParticles[1], TotalLz[1], LzMax[1]);
	      FinalSpace = new BosonOnSphereShort(NbrParticles[0], TotalLz[0] + TotalLz[1], LzMax[0] + LzMax[1]);
	    }
	  else
	    {
	      OutputBasis2 = new FermionOnSphere(NbrParticles[1], TotalLz[1], LzMax[1]);
	      FinalSpace = new FermionOnSphere(NbrParticles[0], TotalLz[0] + TotalLz[1], LzMax[0] + LzMax[1]);
	    }
	}
    }
  
  RealVector GroundState1;
  RealVector GroundState2;
  
  if (InverseFlag == false)
    {	
      if (GroundState1.ReadVector (VectorFiles[0]) == false)
	{
	  cout << "can't open vector file " << VectorFiles[0] << endl;
	  return -1;      
	}
      
      if (GroundState2.ReadVector (VectorFiles[1]) == false)
	{
	  cout << "can't open vector file " << VectorFiles[1]<< endl;
	  return -1;
	}
    }
  else
    {
      if (GroundState1.ReadVector (VectorFiles[1]) == false)
	{
	  cout << "can't open vector file " << VectorFiles[1] << endl;
	  return -1;      
	}
      
      if (GroundState2.ReadVector (VectorFiles[0]) == false)
	{
	  cout << "can't open vector file " << VectorFiles[0]<< endl;
	  return -1;
	}
    }
  
  
  if (OutputBasis1->GetHilbertSpaceDimension() != GroundState1.GetVectorDimension())
    {
      cout << "Number of rows of the vector is not equal to the Hilbert space dimension!";
      return -1;
    }
  
  if (OutputBasis2->GetHilbertSpaceDimension() != GroundState2.GetVectorDimension())
    {
      cout << "Number of rows of the vector is not equal to the Hilbert space dimension!";
      return -1;
    }
  
  
  RealVector OutputVector(FinalSpace->GetHilbertSpaceDimension(),true);
  MonomialsProductOperation Operation(OutputBasis1,OutputBasis2,FinalSpace,&GroundState1,&GroundState2,&OutputVector,NormalizeFlag);
  Operation.ApplyOperation(Architecture.GetArchitecture());
  
  
  OutputVector.WriteVector(OutputFileName);	
  ofstream File;
  if(OutputTxtFileName != 0)
    {
      File.open(OutputTxtFileName, ios::binary | ios::out);
      File.precision(14);
      for (long i = 0; i < FinalSpace->GetLargeHilbertSpaceDimension(); ++i)
	{
	  File << OutputVector[i] << " ";
	  FinalSpace->PrintStateMonomial(File, i) << endl;
	}
    }
  File.close();
  return 0;
}
