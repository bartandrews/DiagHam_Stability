#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"

#include "Hamiltonian/ParticleOnSphereTwoLandauLevelL2Hamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereTwoLandauLevelL2Diagonalize" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 5);
  (*SystemGroup) += new BooleanOption  ('\n', "restrict-polarization", "restrict number of particles in each Landau level (provided by nbrparticles-up and nbrparticles-down)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrparticles-up", "number of particles in N=1 LL", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrparticles-down", "number of particles in N=0 LL", 0);  
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 8);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "l2");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistic");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 1.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "energy-shift", "if non zero, override energy shift using the indicated value ", -10.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");

  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereTwoLandauLevelL2Diagonalize -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  bool PolarizationFlag = Manager.GetBoolean("restrict-polarization");  
  int NbrParticlesUp = Manager.GetInteger("nbrparticles-up");
  int NbrParticlesDown = Manager.GetInteger("nbrparticles-down");
  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz  = Manager.GetInteger("total-lz");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool FirstRun = true;
  char* OutputNameLz = new char [256 + strlen(Manager.GetString("interaction-name"))];
  if (Manager.GetBoolean("boson") == false)
    sprintf (OutputNameLz, "fermions_%s_n_%d_2s_%d_lz.dat", Manager.GetString("interaction-name"), NbrParticles, LzMax);
  else
    sprintf (OutputNameLz, "bosons_%s_n_%d_2s_%d_lz.dat", Manager.GetString("interaction-name"), NbrParticles, LzMax);
  int LzMaxUp = LzMax + 2;
  int LzMaxDown = LzMax;
  cout << "Total Lz = " << TotalLz << endl;
  cout << "Warning: this code will not run correctly if TotalLz is not defined!" << endl;

  ParticleOnSphereWithSpin* Space = 0;
  if (Manager.GetBoolean("boson") == false)
    {
      if (PolarizationFlag) 
          Space = new FermionOnSphereTwoLandauLevels(NbrParticlesUp, NbrParticlesDown, TotalLz, LzMaxUp, LzMaxDown);
      else
          Space = new FermionOnSphereTwoLandauLevels(NbrParticles, TotalLz, LzMaxUp, LzMaxDown);
    }
  else
    {
      cout << "bosons implementation not working" << endl;
      exit(2);
      Space = new BosonOnSphereTwoLandauLevels (NbrParticles, TotalLz, LzMaxUp, LzMaxDown);      
    }

  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    Memory = Architecture.GetArchitecture()->GetLocalMemory();
  AbstractQHEOnSphereHamiltonian* Hamiltonian = new ParticleOnSphereTwoLandauLevelL2Hamiltonian(Space, NbrParticles, LzMax+2, TotalLz,Manager.GetDouble("l2-factor"),
												Architecture.GetArchitecture(), 
												Memory, DiskCacheFlag,
												LoadPrecalculationFileName, false);

   if (Manager.GetString("energy-expectation") != 0 )
    {
      char* StateFileName = Manager.GetString("energy-expectation");
      if (IsFile(StateFileName) == false)
       {
         cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
         return -1;           
       }
      RealVector State;
      if (State.ReadVector(StateFileName) == false)
      {
        cout << "error while reading " << StateFileName << endl;
        return -1;
      }
     if (State.GetVectorDimension()!=Space->GetHilbertSpaceDimension())
      {
        cout << "error: vector and Hilbert-space have unequal dimensions"<<endl;
        return -1;
      }
     RealVector TmpState(Space->GetHilbertSpaceDimension());
     VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
     Operation.ApplyOperation(Architecture.GetArchitecture());
     double EnergyValue = State*TmpState;
     double Shift = 0.25 * TotalLz * TotalLz - 0.5 * TotalLz;
     cout << "< L^2 > = " << (EnergyValue + Shift) << endl;
     cout << "L = "<<0.5*(sqrt(4.0 * (EnergyValue + Shift) + 1.0) - 1.0)<<endl;
     return 0;
   }

  double Shift = -(0.25 * TotalLz * TotalLz - 0.5 * TotalLz); //Manager.GetDouble("energy-shift");
  //Hamiltonian->ShiftHamiltonian(Shift);
  char* EigenvectorName = 0;
  if (Manager.GetBoolean("eigenstate") == true)	
    {
      EigenvectorName = new char [64];
      if (Manager.GetBoolean("boson") == false)
	sprintf (EigenvectorName, "fermions_%s_n_%d_2s_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, LzMax, TotalLz);
      else
	sprintf (EigenvectorName, "bosons_%s_n_%d_2s_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, LzMax, TotalLz);
    }
  QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, TotalLz, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
  MainTaskOperation TaskOperation (&Task);
  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
  if (EigenvectorName != 0)
    delete[] EigenvectorName;
  delete Hamiltonian;
  if (FirstRun == true)
    FirstRun = false;
  delete Space;
  return 0;
}
