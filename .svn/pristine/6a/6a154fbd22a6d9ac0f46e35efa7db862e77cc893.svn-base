#include "config.h"

#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"

#include "Operator/ParticleOnSphereWithSpinSMinusOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSpinAllSz.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereWithSpinSMinus" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleStringOption  ('i', "input-file", "name of the file corresponding to the input state");
  (*SystemGroup) += new  BooleanOption ('\n', "fully-unpolarized", "apply the S- minus as many times as need to reach the Sz=0 (or Sz=1/2) sector");
  (*SystemGroup) += new BooleanOption  ('\n', "all-sz", "use Hilbert-space comprising all sz sectors", 0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (changing the Sz value)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

  if (Manager.GetString("input-file")==NULL)
    {
      cout << "An input file is required!"<<endl;
      exit(-1);
    }

  int NbrParticles = 0;
  int LzMax = 0;
  int Lz = 0;
  int TotalSz = 0;
  bool Statistics = true;
  if (Manager.GetBoolean("all-sz"))
    TotalSz = -1;
  char* InputStateFile = new char [strlen(Manager.GetString("input-file")) + 1];
  strcpy(InputStateFile, Manager.GetString("input-file"));

  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(InputStateFile, NbrParticles, LzMax, Lz, 
							   TotalSz, Statistics) == false)
    {
      return -1;
    }

  int Parity = Lz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the same parity" << endl;
      return -1;           
    }
  if ((Manager.GetBoolean("all-sz")==false) && ((NbrParticles&1) != (TotalSz&1)))
    {
      cout << "Sz and NbrParticles must have the same parity" << endl;
      return -1;
    }

  RealVector InputState; 
  RealVector TargetVector; 
  if (InputState.ReadVector(Manager.GetString("input-file")) == false)
    {
      cout << "error while reading " << Manager.GetString("input-file") << endl;
      return -1;
    }
	
  unsigned long MemorySpace = 9ul << 20;
  ParticleOnSphereWithSpin* InputSpace;
  ParticleOnSphereWithSpin* OutputSpace;
  if (Statistics == true)
    {
      if (Manager.GetBoolean("all-sz")==false)
	{
#ifdef __64_BITS__
	  if (LzMax <= 31)
	    {
	      InputSpace = new FermionOnSphereWithSpin(NbrParticles, Lz, LzMax, TotalSz, MemorySpace);
	    }
	  else
	    {
	      // InputSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	      cout << "Fermions with Spin not defined yet for LzMax > 31"<<endl;
	      exit(-1);
	    }
#else
	  if (LzMax <= 15)
	    {
	      InputSpace = new FermionOnSphereWithSpin(NbrParticles, Lz, LzMax, TotalSz, MemorySpace);
	    }
	  else
	    {
	      // InputSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	      cout << "Fermions with Spin not defined yet for LzMax > 15, consider using a 64 bit machine!"<<endl;
	      exit(-1);
	    }
#endif
	}
      else
	{
	  InputSpace = new FermionOnSphereWithSpinAllSz (NbrParticles, Lz, LzMax, MemorySpace);
	}
    }
  else
    {
      if (Manager.GetBoolean("all-sz")==false)
	{
	  InputSpace = new BosonOnSphereWithSpin(NbrParticles, Lz, LzMax, TotalSz);
	}
      else
	{
	  InputSpace = new BosonOnSphereWithSpinAllSz(NbrParticles, Lz, LzMax, MemorySpace);
	}
    }

  int MinSz = TotalSz - 2;
  int CurrentSz = MinSz;
  if (Manager.GetBoolean("fully-unpolarized"))
    {
      if ((NbrParticles & 1) == 0)
	MinSz = 0;
      else
	MinSz = 1;
    }
  char* OutputFileName = 0;
  if (Manager.GetString("output-file") == 0)
    {
      char* TmpString = strstr(InputStateFile, "_sz_");
      if (TmpString == 0)
	{
	  cout << InputStateFile << " file name does not contain _sz_, can't create output file name" << endl;
	  return 0;
	}
      int SizeString = 4;
      while ((TmpString[SizeString] != '\0') && (TmpString[SizeString] != '_') && (TmpString[SizeString] != '.') &&
	     (((TmpString[SizeString] >= '0') && (TmpString[SizeString] <= '9')) || (TmpString[SizeString] != '-')))
	++SizeString;
       if (((TmpString[SizeString] == '_') || (TmpString[SizeString] == '.')) && (SizeString != 4))
	 {
	   char TmpChar = TmpString[4];
	   TmpString[4] = '\0';
	   OutputFileName = new char [strlen(InputStateFile) + 16 + strlen(TmpString + SizeString)];
	   sprintf (OutputFileName, "%s%d%s", InputStateFile, MinSz, TmpString + SizeString);
	   TmpString[4] = TmpChar;
	 }
       else
	 {
	   cout << "error while replacing sz value in " << InputStateFile << endl;
	   return 0;
	 }     
    }
  else
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }

  while (CurrentSz >= MinSz)
    {
      ParticleOnSphereWithSpin* OutputSpace;
      if (Statistics == true)
	{
	  if (Manager.GetBoolean("all-sz") == false)
	    {
#ifdef __64_BITS__
	      if (LzMax <= 31)
		{
		  OutputSpace = new FermionOnSphereWithSpin(NbrParticles, Lz, LzMax, CurrentSz, MemorySpace);
		  cout << NbrParticles << " " <<  Lz<< " " << LzMax<< " " << CurrentSz << endl;
		}
	      else
		{
		  // OutputSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax - (2 * i), TotalSz, MemorySpace);
		  cout << "Fermions with Spin not defined yet for LzMax > 31"<<endl;
		  exit(-1);
		}
#else
	      if (LzMax <= 15)
		{
		  OutputSpace = new FermionOnSphereWithSpin(NbrParticles, Lz, LzMax, CurrentSz, MemorySpace);
		}
	      else
		{
		  // OutputSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax - (2 * i), TotalSz, MemorySpace);
		  cout << "Fermions with Spin not defined yet for LzMax > 15, consider using a 64 bit machine!"<<endl;
		  exit(-1);
		}
#endif
	    }
	  else
	    {
	      OutputSpace = new FermionOnSphereWithSpinAllSz (NbrParticles, Lz, LzMax, MemorySpace);
	    }
	}
      else
	{
	  if (Manager.GetBoolean("all-sz")==false)
	    OutputSpace = new BosonOnSphereWithSpin(NbrParticles, Lz, LzMax, CurrentSz, MemorySpace);
	  else
	    {
	      OutputSpace = new BosonOnSphereWithSpinAllSz(NbrParticles, Lz, LzMax, MemorySpace);
	    }
	}
      InputSpace->SetTargetSpace(OutputSpace);
      TargetVector = RealVector(OutputSpace->GetHilbertSpaceDimension());
      if (OutputSpace->GetHilbertSpaceDimension() != InputSpace->GetTargetHilbertSpaceDimension())
	{
	  cout << "Problem with setting target space"<<endl;
	  exit(-1);
	}
      ParticleOnSphereWithSpinSMinusOperator SMinusOperator(InputSpace);
      RealVector TmpVector (OutputSpace->GetHilbertSpaceDimension());
      VectorOperatorMultiplyOperation Operation(&SMinusOperator, &InputState, &TmpVector);
      Operation.ApplyOperation(Architecture.GetArchitecture());
      InputState = TmpVector;
      InputState /= InputState.Norm();
      delete InputSpace;
      InputSpace = OutputSpace;
      CurrentSz -= 2;
    }
  if (InputState.WriteVector(OutputFileName) == false)
    {
      cout << "error while writing " <<  OutputFileName<< endl;
      return -1;
    }
  delete InputSpace;
  delete [] OutputFileName;
  return 0;
}

