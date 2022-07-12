#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"

#include "Operator/ParticleOnSphereDensityOperator.h"

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
  OptionManager Manager ("FQHESphereSingleModeApproximation" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
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

  (*SystemGroup) += new SingleStringOption('\n', "ground-state", "name of the file containing the ground state vector upon which rho_k acts");
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistics (default value))");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "lz-boost", "Lz momentum that has to be transfered via the SMA", 0);
  (*SystemGroup) += new BooleanOption ('\n', "compute-bilinears", "compute the action of all the bilinear operators on the ground state");
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "sma");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
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
  int LzBoost = Manager.GetInteger("lz-boost");
  bool Statistics = true;
  if (Manager.GetBoolean("boson") == true)
    {
      Statistics = false;
    }
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("ground-state"),
						   NbrParticles, LzMax, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("ground-state") << endl;
      return -1;
    }
  cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << LzMax << " Lz=" << TotalLz << " ";
  if (Statistics == false)
    {
      cout << "bosons";
    }
  else
    {
      cout << "fermions";
    }
  cout << endl;

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
 
  char* OutputNamePrefix = new char [512];
  ParticleOnSphere* TotalSpace = 0;
  ParticleOnSphere* TargetSpace = 0;
  int ResultingTotalLz = TotalLz + (2 * LzBoost);
  
  if (Statistics == false)
    {
#ifdef  __64_BITS__
      if ((LzMax + NbrParticles - 1) < 63)
#else
	if ((LzMax + NbrParticles - 1) < 31)	
#endif
	  {
	    TotalSpace = new BosonOnSphereShort(NbrParticles, TotalLz, LzMax);
	    TargetSpace = new BosonOnSphereShort(NbrParticles, ResultingTotalLz, LzMax);
            ((BosonOnSphereShort*)TotalSpace)->SetTargetSpace(TargetSpace);
	  }
	else
	  {
	    TotalSpace = new BosonOnSphere(NbrParticles, TotalLz, LzMax);
	    TargetSpace = new BosonOnSphereShort(NbrParticles, ResultingTotalLz, LzMax);
            ((BosonOnSphere*)TotalSpace)->SetTargetSpace(TargetSpace);
	  }
      sprintf (OutputNamePrefix, "bosons_%s_n_%d_2s_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, LzMax, ResultingTotalLz);
    }
  else
    {
      TotalSpace = new FermionOnSphere (NbrParticles, TotalLz, LzMax);
      TargetSpace = new FermionOnSphere (NbrParticles, ResultingTotalLz, LzMax);
      ((FermionOnSphere*)TotalSpace)->SetTargetSpace(TargetSpace);

      sprintf (OutputNamePrefix, "fermions_%s_n_%d_2s_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, LzMax, ResultingTotalLz);
    }


  Architecture.GetArchitecture()->SetDimension(TotalSpace->GetHilbertSpaceDimension());

  char* StateFileName = Manager.GetString("ground-state");
  if (IsFile(StateFileName) == false)
    {
      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }

  RealVector InputState;
  if (InputState.ReadVector(StateFileName) == false)
    {
      cout << "error while reading " << StateFileName << endl;
      return -1;
    }
  if (InputState.GetVectorDimension() != TotalSpace->GetHilbertSpaceDimension())
    {
      cout << "error: vector and Hilbert-space have unequal dimensions " << InputState.GetVectorDimension() << " "<< TotalSpace->GetHilbertSpaceDimension() << endl;
      return -1;
    }
 
  if (Manager.GetBoolean("compute-bilinears"))
    {
      RealVector TmpState(TargetSpace->GetHilbertSpaceDimension());
      int MinLzValue = 0;
      int MaxLzValue = LzMax;
      if (LzBoost >= 0)
	{
	  MaxLzValue = LzMax - LzBoost;	  
	}
      else
	{
	  MinLzValue = -LzBoost;	  
	}
      for (int m = MinLzValue; m <= MaxLzValue; ++m)
	{
	  cout << "computing c^+_"<< (m + LzBoost) << " c_" << m << " |Psi>" << endl;
	  ParticleOnSphereDensityOperator TmpOperator(TotalSpace, m + LzBoost, m);
	  VectorOperatorMultiplyOperation Operation(&TmpOperator, &InputState, &TmpState);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  char* OutputNameLz = new char [strlen(OutputNamePrefix)+ 16];
	  sprintf (OutputNameLz, "%s.%d.vec", OutputNamePrefix, (m - MinLzValue));
	  TmpState.WriteVector(OutputNameLz);
	}
      return 0;
    }
  else //proceed to explicitly compute the SMA wavefunction with given L = LzBoost
   {
      cout << "Computing SMA state "<<endl;
      cout << "Using the expression |psi_L> = sum_k (-1)^(S-k) <S,L+k; S,-k| L,L> c_{L+k}^+ c_{k} |psi_0> " << endl;

      ClebschGordanCoefficients Coefficients(LzMax, LzMax);
      
      int MinLzValue = 0;
      int MaxLzValue = LzMax;
      if (LzBoost >= 0)
        {
          MaxLzValue = LzMax - LzBoost;   
        }
      else
        {
          MinLzValue = -LzBoost;    
        }

      RealVector SMAState(TargetSpace->GetHilbertSpaceDimension(), true);
      RealVector TmpState(TargetSpace->GetHilbertSpaceDimension());
       for (int m = MinLzValue; m <= MaxLzValue; ++m)
        {
          cout << "computing c^+_"<< (m + LzBoost) << " c_" << m << " |Psi>" << endl;

          ParticleOnSphereDensityOperator* TmpOperator;
          TmpOperator = new ParticleOnSphereDensityOperator(TotalSpace, m + LzBoost, m);
          VectorOperatorMultiplyOperation Operation(TmpOperator, &InputState, &TmpState);
          Operation.ApplyOperation(Architecture.GetArchitecture());
    
          double TmpCoeff = Coefficients.GetCoefficient(2*m + 2*LzBoost - LzMax, LzMax-2*m, 2*LzBoost);
          
          if (m%2 == 1)
            TmpCoeff *= -1.0;

          cout << TmpCoeff << " " << TmpState.Norm() << endl; 
          SMAState.AddLinearCombination(TmpCoeff, TmpState);
        }
      cout << "Final state norm " << SMAState.Norm() << endl;
      if (SMAState.Norm() > 1e-10)
	      SMAState /= SMAState.Norm();  
      char* OutputNameLz = new char [strlen(OutputNamePrefix)+ 16];
      sprintf (OutputNameLz, "%s.0.vec", OutputNamePrefix);
      SMAState.WriteVector(OutputNameLz); 

   }


  cout << " Target Hilbert space dimension = " << TargetSpace->GetHilbertSpaceDimension() << endl;
  cout << " Groundstate Hilbert space dimension = " << TotalSpace->GetHilbertSpaceDimension() << endl;
 
  return 0;
}
