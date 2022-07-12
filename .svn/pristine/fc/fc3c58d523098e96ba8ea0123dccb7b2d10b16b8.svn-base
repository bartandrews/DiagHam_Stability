#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "HilbertSpace/FermionOnTorusWithSpin.h"
#include "HilbertSpace/BosonOnTorusWithSpin.h"

#include "Operator/ParticleOnSphereWithSpinDensityOperator.h"

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

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

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
  OptionManager Manager ("FQHETorusWithSpinFermionsSingleModeApproximation" , "0.01");
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

  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistics (default value))");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "max-momentum", "maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "y-momentum", "total momentum in the y direction of the ground state (negative if none or override autodetection from input file name if greater or equal to zero)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the total projection of spin", 0);
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "sma");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "kx", "momentum along x direction of the density operator (negative if none)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ky", "momentum along y direction of the density operator (negative if none)", -1);  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sma-spin", "0 to act with rho_down, 1 to evaluate rho_up", 0);
  (*SystemGroup) += new SingleDoubleOption   ('r', "ratio", 
					      "ratio between lengths along the x and y directions (-1 if has to be taken equal to nbr-particles/4)", -1);
  (*SystemGroup) += new BooleanOption ('\n', "compute-bilinears", "compute the action of all the bilinear operators on the ground state"); 
  (*MiscGroup) += new SingleStringOption('\n', "ground-state", "name of the file containing the ground state vector upon which rho_k acts");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusWithSU2SpinSMinus -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("ground-state") == 0)
    {
      cout << "error, an input state file should be provided. See man page for option syntax or type FQHETorusWithSU2SpinSMinus -h" << endl;
      return -1;
    }
  if ((Manager.GetString("ground-state") != 0) && 
      (IsFile(Manager.GetString("ground-state")) == false))
    {
      cout << "can't open file " << Manager.GetString("ground-state") << endl;
      return -1;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int MaxMomentum = Manager.GetInteger("max-momentum");
  int YMomentum = Manager.GetInteger("y-momentum");
  int TotalSz = Manager.GetInteger("total-sz");
  bool Statistics = true;
  if (Manager.GetBoolean("boson") == true)
    {
      Statistics = false;
    }
  if (FQHEOnTorusWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("ground-state"), NbrParticles, MaxMomentum, YMomentum, TotalSz, Statistics) == false)
   {
     cout << "error while retrieving system parameters from file name " << Manager.GetString("ground-state") << endl;
     return -1;
   }
  cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Ky=" << YMomentum << " ";
  if (Statistics == false)
    {
      cout << "bosons";
    }
  else
    {
      cout << "fermions";
    }
  cout << endl;

  int SMASpin = Manager.GetInteger("sma-spin");
  int Kx = Manager.GetInteger("kx");
  int MaxKx = Kx +1; 
  int Ky = Manager.GetInteger("ky");
  if (Manager.GetBoolean("compute-bilinears") == true)
    {
      if (Ky < 0)
	{
	  Ky = 0;
	}
    }
  else
    {
      cout << "Bilinear calculation only supported at the moment." << endl;
      //if (Kx < 0)
	//{
	//  Kx = 0;
	//  MaxKx = MaxMomentum;
	//}      
    }

  int ResultingYMomentum = (YMomentum + Ky) % MaxMomentum; 
  int MomentumModulo = FindGCD(NbrParticles, MaxMomentum);

  double XRatio = Manager.GetDouble("ratio");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
 
  char* OutputNamePrefix = new char [512];
  ParticleOnSphereWithSpin* TotalSpace = 0;
  ParticleOnSphereWithSpin* TargetSpace = 0;
  
  if (Statistics == false)
    {
#ifdef  __64_BITS__
      if ((MaxMomentum + NbrParticles - 1) < 63)
#else
	if ((MaxMomentum + NbrParticles - 1) < 31)	
#endif
	  {
	    TotalSpace = new BosonOnTorusWithSpin(NbrParticles, MaxMomentum, TotalSz, YMomentum);
	    TargetSpace = new BosonOnTorusWithSpin(NbrParticles, MaxMomentum, TotalSz, ResultingYMomentum);
            ((BosonOnTorusWithSpin*)TotalSpace)->SetTargetSpace(TargetSpace);
	  }
    sprintf (OutputNamePrefix, "bosons_torus_kysym_%s_n_%d_2s_%d_sz_%d_ky_%d", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum, TotalSz, ResultingYMomentum);
    }
  else
    {
      TotalSpace = new FermionOnTorusWithSpin (NbrParticles, MaxMomentum, TotalSz, YMomentum);
      TargetSpace = new FermionOnTorusWithSpin (NbrParticles, MaxMomentum, TotalSz, ResultingYMomentum);
      ((FermionOnTorusWithSpin*)TotalSpace)->SetTargetSpace(TargetSpace);

      sprintf (OutputNamePrefix, "fermions_torus_kysym_%s_n_%d_2s_%d_sz_%d_ky_%d", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum, TotalSz, ResultingYMomentum);
    }


  Architecture.GetArchitecture()->SetDimension(TotalSpace->GetHilbertSpaceDimension());

  char* StateFileName = Manager.GetString("ground-state");
  if (IsFile(StateFileName) == false)
    {
      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }

  RealVector InputState;
  ComplexVector ComplexInputState;
  bool ComplexVectorFlag = false;
  if (InputState.ReadVectorTest(StateFileName) == false)
    {
      ComplexVectorFlag = true;
      if (ComplexInputState.ReadVector(StateFileName) == false)
	{
	  cout << "error while reading " << StateFileName << endl;
	  return -1;
	}
      if (ComplexInputState.GetVectorDimension() != TotalSpace->GetHilbertSpaceDimension())
	{
	  cout << "error: vector and Hilbert-space have unequal dimensions " << ComplexInputState.GetVectorDimension() 
	       << " " << TotalSpace->GetHilbertSpaceDimension() << endl;
	  return -1;
	}
    }
  else
    {
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
    }

  if (Manager.GetBoolean("compute-bilinears"))
    {
      if (ComplexVectorFlag == false)
	{
	  RealVector TmpState(TargetSpace->GetHilbertSpaceDimension());
	  for (int m = 0; m < MaxMomentum; ++m)
	    {
	      cout << "computing c^+_"<< ((m + Ky) % MaxMomentum) << " c_" << m << " |Psi>" << endl;
	      ParticleOnSphereWithSpinDensityOperator TmpOperator(TotalSpace, (m + Ky) % MaxMomentum, SMASpin, m, SMASpin);
	      VectorOperatorMultiplyOperation Operation(&TmpOperator, &InputState, &TmpState);
	      Operation.ApplyOperation(Architecture.GetArchitecture());
	      char* OutputNameLz = new char [strlen(OutputNamePrefix)+ 16];
	      sprintf (OutputNameLz, "%s.%d.vec", OutputNamePrefix, m);
	      TmpState.WriteVector(OutputNameLz);
	    }
	}
      else
	{
	  ComplexVector TmpState(TargetSpace->GetHilbertSpaceDimension());
	  for (int m = 0; m < MaxMomentum; ++m)
	    {
	      cout << "computing c^+_"<< ((m + Ky) % MaxMomentum) << " c_" << m << " |Psi>" << endl;
	      ParticleOnSphereWithSpinDensityOperator TmpOperator(TotalSpace, (m + Ky) % MaxMomentum, SMASpin, m, SMASpin);
	      VectorOperatorMultiplyOperation Operation(&TmpOperator, &ComplexInputState, &TmpState);
	      Operation.ApplyOperation(Architecture.GetArchitecture());
	      char* OutputNameLz = new char [strlen(OutputNamePrefix)+ 16];
	      sprintf (OutputNameLz, "%s.%d.vec", OutputNamePrefix, m);
	      TmpState.WriteVector(OutputNameLz);
	    }
	}
      return 0;
    }

 return 0;
}
