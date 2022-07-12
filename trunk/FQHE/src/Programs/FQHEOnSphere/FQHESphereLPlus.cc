#include "config.h"

#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Operator/ParticleOnSphereLPlusOperator.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"

#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"

#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereLPlus" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* DataGroup = new OptionGroup ("data options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += DataGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz", "twice the total lz value of the system for the initial state", 0);
  (*SystemGroup) += new SingleStringOption  ('s', "statistics", "particle statistics (boson or fermion, try to guess it from file name if not defined)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lp", "number of time the L^+ operator has to be applied", 1);

  (*SystemGroup) += new SingleStringOption  ('\n', "input-reference", "use a haldane basis with the given reference file for the input file");
  (*SystemGroup) += new SingleStringOption  ('\n', "output-reference", "use a haldane basis with the given reference file for the output file");
  (*SystemGroup) += new BooleanOption  ('\n', "save-all", "save all states when applying L^+ multiple times");
  (*SystemGroup) += new BooleanOption  ('\n', "check-l2", "check the L^2 value after applying L^+");

  (*DataGroup) += new SingleStringOption  ('i', "input-file", "input vector file name");
  (*DataGroup) += new SingleStringOption  ('o', "output-file", "output vector file name");
  (*DataGroup) += new SingleStringOption  ('\n', "interaction-name", "interaction name for the output files when computing more than one state", "lplus");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereLPlus -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int Lz = Manager.GetInteger("lz");
  int NbrLPlus = Manager.GetInteger("nbr-lp");
  bool FermionFlag = false;
  if (Manager.GetString("statistics") == 0)
    FermionFlag = true;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("input-file"), NbrParticles, LzMax, Lz, FermionFlag) == false)
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
  int Parity = Lz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the same parity" << endl;
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
  ParticleOnSphere* IntialSpace;
  ParticleOnSphere* TargetSpace;
  cout << "Creating space for Lz="<<Lz<<"/2"<<endl;
  if (FermionFlag == true)
    {
#ifdef __64_BITS__
      if (LzMax <= 63)
	{
	  IntialSpace = new FermionOnSphere(NbrParticles, Lz, LzMax, MemorySpace);
	}
      else
	{
	  IntialSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	}
#else
      if (LzMax <= 31)
	{
	  IntialSpace = new FermionOnSphere(NbrParticles, Lz, LzMax, MemorySpace);
	}
      else
	{
	  IntialSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	}
#endif
    }
  else
    {
      IntialSpace = new BosonOnSphereShort(NbrParticles, Lz, LzMax);
    }
  for (int i = 1; i <= NbrLPlus; ++i)
    {
//       if (Lz + (2 * i) > LzMax)
// 	{
// 	  cout << "Cannot apply LPlus more than "<<i-1<<" times"<<endl;
// 	  exit(-1);
// 	}
      cout << "Creating space for Lz="<<Lz + (2 * i)<<"/2"<<endl;
      if (FermionFlag == true)
	{
#ifdef __64_BITS__
	  if (LzMax <= 63)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles, Lz + (2 * i), LzMax, MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles, Lz + (2 * i), LzMax, MemorySpace);
	    }
#else
	  if (LzMax <= 31)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles, Lz + (2 * i), LzMax, MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles, Lz + (2 * i), LzMax, MemorySpace);
	    }
#endif
	}
      else
	{
	  TargetSpace = new BosonOnSphereShort(NbrParticles, Lz + (2 * i), LzMax);
	}
      IntialSpace->SetTargetSpace(TargetSpace);
      TargetVector = RealVector(TargetSpace->GetHilbertSpaceDimension(), true);
      if (TargetSpace->GetHilbertSpaceDimension()!=IntialSpace->GetTargetHilbertSpaceDimension())
	{
	  cout << "Problem with setting target space"<<endl;
	  exit(-1);
	}
      ParticleOnSphereLPlusOperator LPlus(IntialSpace, Lz + (2 * i) - 2, LzMax);
      LPlus.Multiply(InitialVector, TargetVector);
      delete IntialSpace;
      IntialSpace = TargetSpace;
      InitialVector = TargetVector;
      double TmpNorm = InitialVector.Norm();
      cout << "Norm of the L^+|psi>= " << TmpNorm << endl;
      if (TmpNorm < 1e-12)
	 cout << "Very small norm, vector will not be normalized." << endl;
      else
        { 
	 InitialVector /= InitialVector.Norm();
         cout << "Normalizing the vector..." << endl;
        }  
      if (Manager.GetBoolean("check-l2") == true)
	{
	  ParticleOnSphereL2Hamiltonian Hamiltonian (IntialSpace, NbrParticles, LzMax, (Lz + (2 * i)), Architecture.GetArchitecture(), 1.0, 0);
	  RealVector TmpState(IntialSpace->GetHilbertSpaceDimension());
	  VectorHamiltonianMultiplyOperation Operation (&Hamiltonian, &InitialVector, &TmpState);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  double L2Value = TmpState * InitialVector;
	  double RawTmpAngularMomentum = 0.5 * (sqrt ((4.0 * L2Value) + 1.0) - 1.0);
	  cout << "checking L^2 value at Lz = " << (Lz + (2 * i)) << " : " << endl
	       << "<L^2> = " << L2Value << endl
	       << "<L> = " << RawTmpAngularMomentum << endl
	       << "-----------------------------------" << endl;
	}
      if (Manager.GetBoolean("save-all") == true)
	{
	  char* OutputName = new char [strlen(Manager.GetString("interaction-name")) + 256];
	  if (FermionFlag == true)
	    {
	      sprintf (OutputName, "fermions_%s_n_%d_2s_%d_lz_%d.0.vec", Manager.GetString("interaction-name"), NbrParticles, LzMax, (Lz + (2 * i)));
	    }
	  else
	    {
	      sprintf (OutputName, "bosons_%s_n_%d_2s_%d_lz_%d.0.vec", Manager.GetString("interaction-name"), NbrParticles, LzMax, (Lz + (2 * i)));
	    }
	  if (InitialVector.WriteVector(OutputName) == false)
	    {
	      cout << "error while writing " << OutputName << endl;
	      return -1;
	    }
	  delete[] OutputName;
	}
    }
  if (Manager.GetBoolean("save-all") == false)
    {
      char *OutputName;
      if (Manager.GetString("output-file")!=NULL)
	{
	  OutputName = new char[strlen(Manager.GetString("output-file"))+1];
	  strcpy(OutputName,Manager.GetString("output-file"));
	}
      else
	{
	  OutputName = new char[strlen(Manager.GetString("input-file"))+10];
	  sprintf(OutputName,"%s_L-",Manager.GetString("input-file"));
	}  
      if (InitialVector.WriteVector(OutputName) == false)
	{
	  cout << "error while writing " << OutputName << endl;
	  return -1;
	}
      delete [] OutputName;
    }
  delete IntialSpace;
  return 0;
}

