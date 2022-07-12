#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "Hamiltonian/AbstractQHEHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/InternalReorthogonalizedLanczosAlgorithm.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#ifdef __MPI__
#include <mpi.h>
#endif


using std::ios;
using std::cout;
using std::endl;

using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("QHEFermionsSphereWithSpin" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += ToolsGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "full-diag", 
						"maximum Hilbert space dimension for which full diagonalization is applied", 
						500, true, 100);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup)  += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "block-size", "size of the block used in the block Lanczos algorithm", 2);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "limit-time", "use limit in time instead of a number of lanczos iteration (0 if none, time in seconds)", 0);
  (*LanczosGroup)  += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*LanczosGroup) += new BooleanOption  ('\n', "force-reorthogonalize", 
					 "force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate", "evaluate eigenstates", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Lanczos iteration", false); 
  (*LanczosGroup) += new SingleStringOption  ('\n', "initial-vector", "use file as the initial vector for the Lanczos algorithm" , 0);
  (*LanczosGroup) += new  BooleanOption ('\n', "partial-lanczos", "only run a given number of Lanczos iterations" , false);
  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 15);
  (*SystemGroup) += new SingleIntegerOption  ('s', "SzTotal", "twice the z component of the total spin of the system", 0);
  (*SystemGroup) += new SingleDoubleOption  ('v', "V0-Interaction", "Interaction in s-wave channel", 0);
  (*SystemGroup) += new SingleDoubleOption  ('w', "V1-Interaction", "Interaction in p-wave channel", 1.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace","-g",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
  (*PrecalculationGroup) += new BooleanOption  ('\n', "allow-disk-storage", "expand memory for fast multiplication using disk storage",false);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsSphereWithSpin -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  bool GroundFlag = ((BooleanOption*) Manager["ground"])->GetBoolean();
  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int SzTotal = ((SingleIntegerOption*) Manager["SzTotal"])->GetInteger();
  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) Manager["fast-search"])->GetInteger()) << 20;
  int InitialLz = ((SingleIntegerOption*) Manager["initial-lz"])->GetInteger();
  int NbrLz = ((SingleIntegerOption*) Manager["nbr-lz"])->GetInteger();
  double V0 = ((SingleDoubleOption*) Manager["V0-Interaction"])->GetDouble();
  double V1 = ((SingleDoubleOption*) Manager["V1-Interaction"])->GetDouble();
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();
  char* SavePrecalculationFileName = ((SingleStringOption*) Manager["save-precalculation"])->GetString();
  bool onDiskCacheFlag = ((BooleanOption*) Manager["allow-disk-storage"])->GetBoolean();
  bool FirstRun = true;
  
#ifdef __MPI__
  char* OutputNameLz = Manager.GetFormattedString("fermions_sphere_spin_n_%nbr-particles%_2S_%lzmax%_Sz_%SzTotal%_V_%V0-Interaction%_W_%V1-Interaction%_lz%mpi%%ground%.dat");
#else 
  char* OutputNameLz = Manager.GetFormattedString("fermions_sphere_spin_n_%nbr-particles%_2S_%lzmax%_Sz_%SzTotal%_V_%V0-Interaction%_W_%V1-Interaction%_lz%ground%.dat");
#endif


  int NbrUp = (NbrFermions + SzTotal)/2;
  int NbrDown = (NbrFermions - SzTotal)/2;
  if ((NbrUp+NbrDown != NbrFermions) || ( NbrUp < 0 ) || (NbrDown < 0 ))
    {
      cout << "This value of Sz cannot be achieved with this particle number!" << endl;
      exit(5);
    }
  int Max = ((LzMax - NbrUp + 1) * NbrUp) + ((LzMax - NbrDown + 1) * NbrDown);

  int  L = 0;
  if ((abs(Max) & 1) != 0)
     L = 1;
  if (InitialLz >= 0)
    {
      L = InitialLz;
      if ((abs(Max) & 1) != 0)
	L |= 1;
      else
	L &= ~0x1;
    }
  if (GroundFlag == true)
      Max = L;
  else
    {
      if (NbrLz > 0)
	{
	  if (L + (2 * (NbrLz - 1)) < Max)
	    Max = L + (2 * (NbrLz - 1));
	}
    }
  for (; L <= Max; L += 2)
    {
      double Shift = -10.0;
      ParticleOnSphereWithSpin* Space;
#ifdef __64_BITS__
      if (LzMax <= 31)
        {
          Space = new FermionOnSphereWithSpin(NbrFermions, L, LzMax, SzTotal, MemorySpace);
        }
      else
	{
	  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
	  return -1;
	}	
#else
      if (LzMax <= 15)
        {
          Space = new FermionOnSphereWithSpin(NbrFermions, L, LzMax, SzTotal, MemorySpace);
	}
      else
	{
	  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
	  return -1;
	}	
#endif
      
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      /*
      // Some testing of the Hilbert Space...
      if (Space->GetHilbertSpaceDimension() < 100)
	for (int i=0; i<Space->GetHilbertSpaceDimension(); i++)
	  { Space->PrintState(cout, i); cout << endl; }
      else
	for (int i=0; i<100; i++)
	  { Space->PrintState(cout, i); cout << endl; }

      // test Four-Point Operators:
      ParticleOnSphere* Space2=0;
      if (NbrFermions==abs(SzTotal))
	{
	  Space2 = new FermionOnSphere(NbrFermions, L, LzMax, MemorySpace);
	  cout << "For comparison: space without spin:" <<endl;
	        if (Space->GetHilbertSpaceDimension() < 100)
		  {
		    if (Space->GetHilbertSpaceDimension() < 100)
		      for (int i=0; i<Space->GetHilbertSpaceDimension(); i++)
			{ Space2->PrintState(cout, i); cout << endl; }
		    else
		      for (int i=0; i<100; i++)
			{ Space2->PrintState(cout, i); cout << endl; }
		  }
	}
      for (int i=0; i<Space->GetHilbertSpaceDimension(); i++)
	{
	  cout << "Applying operators to: "; Space->PrintState(cout, i); cout << endl;
	  for (int m1=0; m1 <= LzMax; ++m1)
	    for (int m2=0; m2 <= m1; ++m2)
	      for (int m3=0; m3 <= LzMax; ++m3)
		{
		  int m4= m1+m2-m3, ret;
		  double coeff=0.0;
		  ret=Space->AddAddAdAd(i,m1,m2,m3,m4,coeff);
		  printf("Pair (m1=%d,d; m2=%d,d; m3=%d,d; m4=%d,d) : ",m1,m2,m3,m4);
		  if (ret < Space->GetHilbertSpaceDimension())
		    {
		      cout << coeff <<"*";
		      Space->PrintState(cout,ret);
		      cout << endl;
		    }
		  else cout << "void" << endl;
		}
	}
      */
      // introduce interaction and test Hamiltonian!
      AbstractQHEHamiltonian* Hamiltonian;

      Hamiltonian = new ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian(Space, NbrFermions, LzMax, V0, V1, Architecture.GetArchitecture(), Memory, onDiskCacheFlag, LoadPrecalculationFileName);
      Hamiltonian->ShiftHamiltonian(Shift);
      if (SavePrecalculationFileName != 0)
	{
	  Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
	}
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [120];
	  sprintf (EigenvectorName, "fermions_sphere_spin_n_%d_2S_%d_Sz_%d_lz_%d_V_%g_W_%g.ev",
		   NbrFermions, LzMax, SzTotal, L, V0, V1);
	}
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());

      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	  EigenvectorName = 0;
	}

      if ((Manager.GetString("initial-vector")!=0))
	{
	  InternalReorthogonalizedLanczosAlgorithm Lanczos(Architecture.GetArchitecture(),
							   /* nbrEigenvalue */ 1, 10, 1000);
	  Lanczos.SetHamiltonian(Hamiltonian);
	  RealVector InitialVector;
	  InitialVector.ReadVector(Manager.GetString("initial-vector"));	  
	  EigenvectorName = new char [120];
	  sprintf (EigenvectorName, "fermions_sphere_spin_n_%d_2S_%d_Sz_%d_lz_%d_V_%g_W_%g.ev-I",
		   NbrFermions, LzMax, SzTotal, L, V0, V1);
	  Lanczos.ProjectVector(InitialVector);
	  InitialVector.WriteVector(EigenvectorName);
	}

      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	  EigenvectorName = 0;
	}
      
      delete Hamiltonian;
      delete Space;
      
      if (FirstRun == true)
	FirstRun = false; 
    }
  delete [] OutputNameLz;
  return 0;
}
