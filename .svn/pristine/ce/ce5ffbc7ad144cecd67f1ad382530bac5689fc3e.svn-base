#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
//#include "HilbertSpace/BosonOnSphereWithSU2SpinAllSz.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"

#include "Operator/ParticleOnSphereWithSpinCreationOperator.h"
#include "Operator/ParticleOnSphereWithSpinAnnihilationOperator.h"

#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundState.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithEigenstates.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundStateFastDisk.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "GeneralTools/ListIterator.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <cstring> 


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereWithSpinApplyCreationOperatorToVector" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* DataGroup = new OptionGroup ("data options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  Manager += ToolsGroup;
  Manager += DataGroup;

 (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system for the initial state", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "spin projection of ground state", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "statistics", "particle statistics (boson or fermion, try to guess it from file name if not defined)");
  (*SystemGroup) += new BooleanOption  ('A', "all-sz", "assume a hilbert space including all sz values");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "orbital-number", "number of the orbital where a particle will be created (0, 1, 2, ..., lzmax)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "spin-number", "spin projection of the added/removed electron (0 for up, 1 for down)", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "input-reference", "use a haldane basis with the given reference file for the input file");
  (*SystemGroup) += new SingleStringOption  ('\n', "output-reference", "use a haldane basis with the given reference file for the output file");
  (*SystemGroup) += new BooleanOption  ('\n', "annihilate-particle", "annihilate particle instead of creating it");

  (*DataGroup) += new SingleStringOption  ('i', "input-file", "input vector file name");
  (*DataGroup) += new SingleStringOption  ('o', "output-file", "output vector file name");
  (*DataGroup) += new SingleStringOption  ('\n', "interaction-name", "interaction name for the output files when computing more than one state", "addremoveparticle");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereSingleModeApproximation -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz = Manager.GetInteger("total-lz");
  int TotalSz = Manager.GetInteger("total-sz");
  int OrbitalNumber = Manager.GetInteger("orbital-number");
  int SpinNumber = Manager.GetInteger("spin-number");
  bool FermionFlag = false;
  if (Manager.GetString("statistics") == 0)
    FermionFlag = true;
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-file"), NbrParticles, LzMax, TotalLz, TotalSz, 
                 FermionFlag) == false)
    {
      return -1;
    }
  if ((Manager.GetString("statistics")) != 0)
    {
      if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	{
	  FermionFlag = true;
	}
      else
	if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	  {
	    FermionFlag = false;
	  }
	else
	  {
	    cout << Manager.GetString("statistics") << " is an undefined statistics" << endl;
	  }  
    }

  int Parity = TotalLz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the same parity" << endl;
      return -1;           
    }
  if ((NbrParticles&1) != (TotalSz&1))
    {
      cout << "Sz and NbrParticles must have the same parity" << endl;
      return -1;
    }

  RealVector InitialVector; 
  RealVector TargetVector; 
  if (InitialVector.ReadVector(Manager.GetString("input-file")) == false)
    {
      cout << "error while reading " << Manager.GetString("input-file") << endl;
      return -1;
    }
	
  long MemorySpace = 9l << 20;
  char* OutputNamePrefix = new char [512];
  ParticleOnSphereWithSpin* InputSpace = 0;
  ParticleOnSphereWithSpin* OutputSpace = 0;

  int ResultingNbrParticles;
  if (Manager.GetBoolean("annihilate-particle") == false)
    ResultingNbrParticles = NbrParticles + 1; 
  else
    ResultingNbrParticles = NbrParticles - 1; 


  int ResultingTotalLz;
  if (Manager.GetBoolean("annihilate-particle") == false)
    ResultingTotalLz = TotalLz + (2 * OrbitalNumber - LzMax);
  else
    ResultingTotalLz = TotalLz - (2 * OrbitalNumber - LzMax);

  int ResultingTotalSz; 
  if ((SpinNumber == 0) && (Manager.GetBoolean("annihilate-particle") == false))
	ResultingTotalSz = TotalSz + 1; // c_{m,up}^+
  else if ((SpinNumber == 1) && (Manager.GetBoolean("annihilate-particle") == false))
	ResultingTotalSz = TotalSz -1;  // c_{m,down}^+
  else if ((SpinNumber == 0) && (Manager.GetBoolean("annihilate-particle") == true))
	ResultingTotalSz = TotalSz - 1; // c_{m,up}
  else 
	ResultingTotalSz = TotalSz + 1; // c_{m,down}

  int ResultingNbrUp = (ResultingNbrParticles + ResultingTotalSz)/2;
  int ResultingNbrDown = (ResultingNbrParticles - ResultingTotalSz)/2;
	
  if ( (ResultingNbrUp < 0) || (ResultingNbrUp > ResultingNbrParticles) || (ResultingNbrDown < 0) || (ResultingNbrDown > ResultingNbrParticles) )
    {
      cout << "Nbr up or down < 0 or > N" << endl; 
      return 0;
    }
	

  if (FermionFlag == true)
    {
      if (Manager.GetBoolean("all-sz")==false)
  {
#ifdef __64_BITS__
	      if (LzMax <= 31)
#else
		if (LzMax <= 15)
#endif
		  {
		        InputSpace = new FermionOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz, MemorySpace);
		        OutputSpace = new FermionOnSphereWithSpin(ResultingNbrParticles, ResultingTotalLz, LzMax, ResultingTotalSz, MemorySpace);

		  }
		else
		  {
#ifdef __128_BIT_LONGLONG__
		    if (LzMax <= 63)
#else
		      if (LzMax <= 31)
#endif
			{
			        InputSpace = new FermionOnSphereWithSpinLong(NbrParticles, TotalLz, LzMax, TotalSz, MemorySpace);
			        OutputSpace = new FermionOnSphereWithSpinLong(ResultingNbrParticles, ResultingTotalLz, LzMax, ResultingTotalSz, MemorySpace);
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
    InputSpace = new FermionOnSphereWithSpinAllSz (NbrParticles, TotalLz, LzMax, MemorySpace);
    OutputSpace = new FermionOnSphereWithSpinAllSz (ResultingNbrParticles, ResultingTotalLz, LzMax, MemorySpace);
  }
      sprintf (OutputNamePrefix, "fermions_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz_%d", Manager.GetString("interaction-name"), ResultingNbrParticles, LzMax, ResultingTotalSz, ResultingTotalLz);
    }
  else
    {
      if (Manager.GetBoolean("all-sz")==false)
  {
    InputSpace = new BosonOnSphereWithSU2Spin(NbrParticles, TotalLz, LzMax, TotalSz);
    OutputSpace = new BosonOnSphereWithSU2Spin(ResultingNbrParticles, ResultingTotalLz, LzMax, ResultingTotalSz);
  }
      else
  {
    InputSpace = new BosonOnSphereWithSU2Spin(NbrParticles, TotalLz, LzMax);
    OutputSpace = new BosonOnSphereWithSU2Spin(ResultingNbrParticles, ResultingTotalLz, LzMax);
  }
      sprintf (OutputNamePrefix, "bosons_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz_%d", Manager.GetString("interaction-name"), ResultingNbrParticles, LzMax, ResultingTotalSz, ResultingTotalLz);
    }

      cout << "2S= " << LzMax << endl;
      cout << "(N, Sz, Lz): " << NbrParticles << " " << TotalSz << " " << TotalLz << " ---> " <<  ResultingNbrParticles << " " << ResultingTotalSz << " " << ResultingTotalLz << endl;
      cout << "Initial dim= " << InputSpace->GetHilbertSpaceDimension() << " Target dim= " << OutputSpace->GetHilbertSpaceDimension() << endl;


      Architecture.GetArchitecture()->SetDimension(InputSpace->GetHilbertSpaceDimension());
      InputSpace->SetTargetSpace(OutputSpace);

      TargetVector = RealVector(OutputSpace->GetHilbertSpaceDimension(), true);
      if (OutputSpace->GetHilbertSpaceDimension()!=InputSpace->GetTargetHilbertSpaceDimension())
	{
	  cout << "Problem with setting target space"<<endl;
	  exit(-1);
	}


      if (Manager.GetBoolean("annihilate-particle") == false)
        {
          cout << "computing c^+_"<< OrbitalNumber << " sigma= " << SpinNumber << " |Psi>" << endl;
          ParticleOnSphereWithSpinCreationOperator TmpOperator(InputSpace, OrbitalNumber, SpinNumber);
          VectorOperatorMultiplyOperation Operation(&TmpOperator, &InitialVector, &TargetVector);
          Operation.ApplyOperation(Architecture.GetArchitecture());          
        } 
      else //annihilate
        {
          cout << "computing c_"<< OrbitalNumber << " sigma= " << SpinNumber << " |Psi>" << endl;
          ParticleOnSphereWithSpinAnnihilationOperator TmpOperator(InputSpace, OrbitalNumber, SpinNumber);
          VectorOperatorMultiplyOperation Operation(&TmpOperator, &InitialVector, &TargetVector);
          Operation.ApplyOperation(Architecture.GetArchitecture());          
        }  
      	

    char* OutputName = new char [strlen(OutputNamePrefix)+ 16];
    sprintf (OutputName, "%s.0.vec", OutputNamePrefix);

    TargetVector.WriteVector(OutputName); 

    cout << "Norm= " << TargetVector.Norm() << endl;  
    delete [] OutputName;
    delete InputSpace;
    return 0;
}
