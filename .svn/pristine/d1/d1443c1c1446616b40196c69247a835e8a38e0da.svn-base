#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"

#include "Hamiltonian/ParticleOnSphereProjectorHamiltonian.h"
#include "Hamiltonian/LinearlySuperposedQHEOnSphereHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereGenericHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Options/Options.h"

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
  OptionManager Manager ("FQHESphereBosonsProjectorHamiltonian" , "0.01");
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

  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "projector-state", "file describing the projector state");
  (*SystemGroup) += new  SingleStringOption ('\n', "multiple-projectors", "one column formatted file giving the list of projector states");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing an optional 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");
  (*SystemGroup) += new BooleanOption  ('\n', "get-lvalue", "compute mean l value from <L^2> for each eigenvalue");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

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
      cout << "see man page for option syntax or type FQHESphereBosonsProjectorHamiltonian -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  bool GroundFlag = Manager.GetBoolean("ground");
  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  bool SymmetrizedBasis = Manager.GetBoolean("symmetrized-basis");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool FirstRun = true;
  double* PseudoPotentials = 0;
  // double* OneBodyPotentials = 0;
  double* FourBodyPotentials = 0;
  int TmpNbrFourBodyPseudoPotentials = 0;

  char* OutputNameLz = new char [256 + strlen(Manager.GetString("interaction-name"))];
  sprintf (OutputNameLz, "bosons_%s_n_%d_2s_%d_lz.dat", Manager.GetString("interaction-name"), NbrParticles, LzMax);

  if ((Manager.GetString("projector-state") == 0) && (Manager.GetString("multiple-projectors") == 0))
    {
      cout << "a projetor state has to be provided" << endl;
      return 0;
    }

  int NbrProjectors = 1;
  int ProjectorNbrParticles = 0;
  int ProjectorLzMax = 0;
  int ProjectorTotalLz = 0;
  BosonOnSphereShort** ProjectorSpaces;  
  RealVector* ProjectorStates;
  bool ProjectorFermionFlag = true;
  if (Manager.GetString("projector-state") != 0)
    {
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("projector-state"), ProjectorNbrParticles, ProjectorLzMax, ProjectorTotalLz, ProjectorFermionFlag) == false)
	{
	  cout << "can't parse information from " << Manager.GetString("projector-state") << endl;
	  return 0;
	}
      ProjectorSpaces = new BosonOnSphereShort*[NbrProjectors];
      ProjectorStates = new RealVector[NbrProjectors];
      ProjectorSpaces[0] = new BosonOnSphereShort(ProjectorNbrParticles, ProjectorTotalLz, ProjectorLzMax);
      if (ProjectorStates[0].ReadVector(Manager.GetString("projector-state")) == false)
	{
	  cout << "error while reading " << Manager.GetString("projector-state") << endl;
	  return -1;
	}
    }
  else
    {
      MultiColumnASCIIFile Description;
      if (Description.Parse(Manager.GetString("multiple-projectors")) == false)
	{
	  Description.DumpErrors(cout);
	  return -1;
	}
      NbrProjectors = Description.GetNbrLines();
      ProjectorSpaces = new BosonOnSphereShort*[NbrProjectors];
      ProjectorStates = new RealVector[NbrProjectors];
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(Description(0, 0), ProjectorNbrParticles, ProjectorLzMax, ProjectorTotalLz, ProjectorFermionFlag) == false)
	{
	  cout << "can't parse information from " << Description(0, 0) << endl;
	  return 0;
	}
      if (ProjectorStates[0].ReadVector(Description(0, 0)) == false)
	{
	  cout << "error while reading " << Description(0, 0) << endl;
	  return -1;
	}
      ProjectorSpaces[0] = new BosonOnSphereShort(ProjectorNbrParticles, ProjectorTotalLz, ProjectorLzMax);
      cout << ProjectorTotalLz << endl;
      for (int i = 1; i < NbrProjectors; ++i)
	{
	  int TmpProjectorNbrParticles = 0;
	  int TmpProjectorLzMax = 0;
	  int TmpProjectorTotalLz = 0;
	  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Description(0, i), TmpProjectorNbrParticles, TmpProjectorLzMax, TmpProjectorTotalLz, ProjectorFermionFlag) == false)
	    {
	      cout << "can't parse information from " << Description(0, i) << endl;
	      return 0;
	    }
	  cout << TmpProjectorTotalLz << endl;
	  if (TmpProjectorNbrParticles != ProjectorNbrParticles)
	    {
	      cout << "error number of particles in " << Description(0, i) << " does not match the one in " << Description(0, 0) << endl;	      
	    }
	  if (ProjectorLzMax != TmpProjectorLzMax)
	    {
	      cout << "error number of flux quanta in " << Description(0, i) << " does not match the one in " << Description(0, 0) << endl;	      
	    }
	  ProjectorSpaces[i] = new BosonOnSphereShort(ProjectorNbrParticles, TmpProjectorTotalLz, ProjectorLzMax);
 	  if (ProjectorStates[i].ReadVector(Description(0, i)) == false)
	    {
	      cout << "error while reading " << Description(0, i) << endl;
	      return -1;
	    }
	}
    }

  int Max = (LzMax * NbrParticles);
  int  L = InitialLz;

  if ((abs(Max) & 1) != (InitialLz & 1))
    L += 1;

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
      AbstractQHEOnSphereHamiltonian* Hamiltonian = 0;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      if (Manager.GetString("interaction-file") == 0)
	{
	  Hamiltonian = new ParticleOnSphereProjectorHamiltonian(Space, NbrParticles, LzMax, ProjectorStates, (ParticleOnSphere**) ProjectorSpaces, 
								 NbrProjectors, ProjectorNbrParticles,
								 ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(),
								 Architecture.GetArchitecture(), 
								 Memory, DiskCacheFlag,
								 LoadPrecalculationFileName);
	}
      else
	{
	  ConfigurationParser InteractionDefinition;
	  double* PseudoPotentials = 0;
	  if (InteractionDefinition.Parse(Manager.GetString("interaction-file")) == false)
	    {
	      InteractionDefinition.DumpErrors(cout) << endl;
	      return -1;
	    }
	  int TmpNbrPseudoPotentials;
	  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', PseudoPotentials, TmpNbrPseudoPotentials) == false)
	    {
	      cout << "Weights is not defined or as a wrong value in " << Manager.GetString("interaction-file") << endl;
	      return -1;
	    }
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials" << endl;
	      return -1;	  
	    }
	  AbstractQHEOnSphereHamiltonian** TmpHamiltonians = new AbstractQHEOnSphereHamiltonian*[2];
	  TmpHamiltonians[0] = new ParticleOnSphereProjectorHamiltonian(Space, NbrParticles, LzMax, ProjectorStates, 
									(ParticleOnSphere**) ProjectorSpaces, NbrProjectors, ProjectorNbrParticles,
									((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(),
									Architecture.GetArchitecture(), 
									Memory, DiskCacheFlag,
									LoadPrecalculationFileName);
	  TmpHamiltonians[1] = new ParticleOnSphereGenericHamiltonian(Space, NbrParticles, LzMax, PseudoPotentials,
								      ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(),
								      Architecture.GetArchitecture(), 
								      Memory);
	  double* TmpCoefficients = new double [2];
	  TmpCoefficients[0] = 1.0;
	  TmpCoefficients[1] = 1.0;	  
	  Hamiltonian = new LinearlySuperposedQHEOnSphereHamiltonian(2, TmpCoefficients, TmpHamiltonians);
	}

      double Shift = - 0.5 * ((double) (NbrParticles * NbrParticles)) / (0.5 * ((double) LzMax));
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [64];
	  sprintf (EigenvectorName, "bosons_%s_n_%d_2s_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, LzMax, L);
	}
      
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      delete Hamiltonian;
      if (FirstRun == true)
	FirstRun = false;
    }

  return 0;
}
