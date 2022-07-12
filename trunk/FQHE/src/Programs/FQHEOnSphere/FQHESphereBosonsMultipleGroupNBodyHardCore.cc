#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "Hamiltonian/ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "GeneralTools/ConfigurationParser.h"

#include <iostream>
#include <cstring>
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

  // some running options and help
  OptionManager Manager ("FQHESphereBosonsMultipleGroupNBodyHardCore" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);
  ParticleOnSphereManager ParticleManager(false, true, 1);
  ParticleManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");

  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-nbody", "number of particle that can interact simultaneously through the n-body hard-core interaction", 2);
  (*SystemGroup) += new  SingleStringOption ('\n', "nbody-file", "file describing which n-body hard-core interactions have to be used");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files, default is nbody)");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");
  (*SystemGroup) += new BooleanOption  ('\n', "get-lvalue", "compute mean l value from <L^2> for each eigenvalue");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "full-diag", "maximum Hilbert space dimension for which full diagonalization is applied", 500, true, 100);

  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonsMultipleGroupNBodyHardCore -h" << endl;
      return -1;
    }
  if ( Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
    
    
  bool GroundFlag = ((BooleanOption*) Manager["ground"])->GetBoolean();
  int NbrBosons = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int NbrNBody = ((SingleIntegerOption*) Manager["nbr-nbody"])->GetInteger();
  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  int InitialLz = ((SingleIntegerOption*) Manager["initial-lz"])->GetInteger();
  int NbrLz = ((SingleIntegerOption*) Manager["nbr-lz"])->GetInteger();
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();  
  bool DiskCacheFlag = ((BooleanOption*) Manager["disk-cache"])->GetBoolean();
  bool FirstRun = true;
  double* NBodyWeightFactors = 0;
  
  int NbrGroups = 2;
  int * NbrBodys = new int [NbrGroups];
  NbrBodys[0] = 2;
  NbrBodys[1] = 2;
  
  /*if (((SingleStringOption*) Manager["nbody-file"])->GetString() != 0)
    {
      ConfigurationParser NBodyDefinition;
      if (NBodyDefinition.Parse(((SingleStringOption*) Manager["nbody-file"])->GetString()) == false)
	{
	  NBodyDefinition.DumpErrors(cout) << endl;
	  return -1;
	}
      if ((NBodyDefinition.GetAsSingleInteger("NbrNBody", NbrNBody) == false) || (NbrNBody < 2))
	{
	  cout << "NbrNBody is not defined or as a wrong value in " << ((SingleStringOption*) Manager["nbody-file"])->GetString() << endl;
	  return -1;
	}
      int TmpNbrNBody;
      if ((NBodyDefinition.GetAsDoubleArray("Weights", ' ', NBodyWeightFactors, TmpNbrNBody) == false) || ((TmpNbrNBody - NbrNBody) != 1))
	{
	  cout << "Weights is not defined or as a wrong value in " << ((SingleStringOption*) Manager["nbody-file"])->GetString() << endl;
	  return -1;
	}
    } */

  char* OutputNameLz;
  if (((SingleStringOption*) Manager["interaction-name"])->GetString() == 0)
    {
      OutputNameLz = Manager.GetFormattedString("bosons_hardcore_nbody_%nbr-nbody%_n_%nbr-particles%_2s_%lzmax%_lz.dat");
    }
  else
    {
      OutputNameLz = Manager.GetFormattedString("bosons_%interaction-name%_n_%nbr-particles%_2s_%lzmax%_lz.dat");
    }

  int Max = (LzMax * NbrBosons);
  int  L = 0;
  if ((abs(Max) & 1) != 0)
     L = 1;
  if (InitialLz >= 0)
    {
      L = InitialLz;
      if ((abs(Max) & 1) != (InitialLz & 1))
	L += 1;
    }
  if (GroundFlag == true)
      Max = L;
  else
    {
      if (NbrLz > 0)
	{
	  Max = L + (2 * (NbrLz - 1));
	}
    }
  for (; L <= Max; L += 2)
    {
      ParticleOnSphere* Space = (ParticleOnSphere*) ParticleManager.GetHilbertSpace(L);
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      AbstractQHEOnSphereHamiltonian* Hamiltonian = 0;
      if (NBodyWeightFactors == 0)
	{
	  Hamiltonian = new ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian(Space, NbrBosons, LzMax, NbrGroups, NbrBodys,
								     ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(), 
								     Architecture.GetArchitecture(), 
								     Memory, DiskCacheFlag,
								     LoadPrecalculationFileName);
	}
	
      double Shift = - 0.5 * ((double) (NbrBosons * NbrBosons)) / (0.5 * ((double) LzMax));
      cout <<"Shift = "<<Shift<<endl;
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  if (((SingleStringOption*) Manager["interaction-name"])->GetString() == 0)
	    {
	      EigenvectorName = new char [256];
	      sprintf (EigenvectorName, "bosons_hardcore_nbody_%d_n_%d_2s_%d_lz_%d", NbrNBody, NbrBosons, LzMax, L);
	    }
	  else
	    {
	      EigenvectorName = new char [256 + strlen(((SingleStringOption*) Manager["interaction-name"])->GetString())];
	      sprintf (EigenvectorName, "bosons_%s_n_%d_2s_%d_lz_%d", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrBosons, LzMax, L);
	    }
	}
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      delete Hamiltonian;
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      if (FirstRun == true)
	FirstRun = false;
      delete Space;
    }

  return 0;
}
