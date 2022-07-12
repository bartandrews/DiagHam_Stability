#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"

#include "Hamiltonian/ParticleOnSphereGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereGenericAnisotropicHamiltonian.h"
#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"

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
  OptionManager Manager ("FQHESphereBosonsTwoBodyGeneric" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("Tools options");
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
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new  SingleStringOption ('\n', "Vm2-file", "file describing the anisotropic pseudopotentials V_m^2, m = 1, 3, 5, etc.");
  (*SystemGroup) += new  SingleStringOption ('\n', "Vm4-file", "file describing the anisotropic pseudopotentials V_m^4, m = 1, 3, 5, etc.");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "chirality", "chirality of the spectral function (+/-1)", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "get-lvalue", "compute mean l value from <L^2> for each eigenvalue");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");

  (*PrecalculationGroup) += new BooleanOption ('\n', "no-hermitian", "do not use hermitian symmetry of the hamiltonian");
  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new SingleStringOption('\n', "conjugate-vector", "name of the file containing the state vector that will be used for inner product with H|psi_0>, where psi_0 is provided by --energy-expectation");
  (*MiscGroup) += new SingleStringOption('\n', "store-vector", "name of the file containing the vector obtained by acting with H on the ground state");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsNBodyHardCore -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  bool GroundFlag = Manager.GetBoolean("ground");
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  if (Manager.GetString("energy-expectation") != 0 ) Memory = 0x0l;
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool FirstRun = true;
  double* PseudoPotentials = new double[LzMax + 1];
  double* OneBodyPotentials = 0; 
  int Chirality = Manager.GetInteger("chirality");
  double* AnisotropicPseudoPotentialsAlpha2 = new double[LzMax + 1];
  double* OneBodyAnisotropicPotentialsAlpha2 = 0;
  double* AnisotropicPseudoPotentialsAlpha4 = new double[LzMax + 1];
  double* OneBodyAnisotropicPotentialsAlpha4 = 0;

  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "an isotropic interaction file not provided" << endl;
      return -1;
    }
  else
    {
      if (FQHESphereGetPseudopotentials(Manager.GetString("interaction-file"), LzMax, PseudoPotentials, OneBodyPotentials) == false)
	return -1;
    }
  char* InteractionName = 0;

  if (Manager.GetString("Vm2-file") == 0)
    {
      cout << "an anisotropic alpha=2 interaction file not provided" << endl;
      return -1;
    }
  else
    {
      if (FQHESphereGetPseudopotentials(Manager.GetString("Vm2-file"), LzMax, AnisotropicPseudoPotentialsAlpha2, OneBodyAnisotropicPotentialsAlpha2) == false)
        {
	  cout << "Problem reading in pseudopotentials." << endl;
	  return -1;
	}
    }

  cout << "Done reading in pseudopotentials alpha=2." << endl;

  if (Manager.GetString("Vm4-file") == 0)
    { 
      cout << "an anisotropic alpha=4 interaction not provided" << endl;
      return -1;
    }
  else
    {
      if (FQHESphereGetPseudopotentials(Manager.GetString("Vm4-file"), LzMax, AnisotropicPseudoPotentialsAlpha4, OneBodyAnisotropicPotentialsAlpha4) == false)
        {
	  cout << "Problem reading in pseudopotentials." << endl;
	  return -1;
	}
    }

  cout << "Done reading in pseudopotentials alpha=4." << endl;

  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(Manager.GetString("interaction-file")) == false)
    {
	  InteractionDefinition.DumpErrors(cout) << endl;
	  exit(-1);
    }
  if (InteractionDefinition["Name"] != 0)
    {
	  InteractionName = new char[strlen(InteractionDefinition["Name"]) + 1];
	  strcpy(InteractionName, InteractionDefinition["Name"]);
    }
  else
    {
	  InteractionName = new char[32 + strlen(Manager.GetString("interaction-name"))];
	  strcpy(InteractionName, Manager.GetString("interaction-name"));
	  InteractionName=new char[15 + strlen(Manager.GetString("interaction-name"))];
          strcpy (InteractionName, Manager.GetString("interaction-name")); 
    }



  char* OutputNameLz = new char [256 + strlen(InteractionName)];
  sprintf (OutputNameLz, "bosons_%s_n_%d_2s_%d_lz.dat", InteractionName, NbrBosons, LzMax);
  int Max = (LzMax * NbrBosons);
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
      ParticleOnSphere* Space = ParticleManager.GetHilbertSpace(L);
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      AbstractQHEOnSphereHamiltonian* Hamiltonian = 0;

      if (OneBodyPotentials == 0)
        Hamiltonian = new ParticleOnSphereGenericAnisotropicHamiltonian(Space, NbrBosons, LzMax, PseudoPotentials,
							     AnisotropicPseudoPotentialsAlpha2,
							     AnisotropicPseudoPotentialsAlpha4,
                                                             Chirality,
							   Manager.GetDouble("l2-factor"),
							   Architecture.GetArchitecture(), 
							   Memory, DiskCacheFlag,
							   LoadPrecalculationFileName, !Manager.GetBoolean("no-hermitian"));
      else
        Hamiltonian = new ParticleOnSphereGenericAnisotropicHamiltonian(Space, NbrBosons, LzMax, PseudoPotentials, OneBodyPotentials, 
							     AnisotropicPseudoPotentialsAlpha2,
							     AnisotropicPseudoPotentialsAlpha4,
                                                             Chirality,
							   Manager.GetDouble("l2-factor"),
							   Architecture.GetArchitecture(), 
							   Memory, DiskCacheFlag,
							   LoadPrecalculationFileName, !Manager.GetBoolean("no-hermitian"));
											 
      double Shift = - 0.5 * ((double) (NbrBosons * NbrBosons)) / (0.5 * ((double) LzMax));

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

	  //RealVector TmpState(Space->GetHilbertSpaceDimension());
	  //VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
	  //Operation.ApplyOperation(Architecture.GetArchitecture());
	  //double EnergyValue = State*TmpState;
	  //cout << "< Energy > = "<<EnergyValue<<endl;
	  //cout << "< shifted energy > = "<<EnergyValue + Shift<<endl;
	  //return 0;

      RealVector BraState;
      char* BraStateFileName = Manager.GetString("conjugate-vector");
      if (IsFile(BraStateFileName) == false)
          {
            cout << "state " << BraStateFileName << " does not exist or can't be opened" << endl;
            return -1;           
          }
      if (BraState.ReadVector(BraStateFileName) == false)
          {
            cout << "error while reading " << BraStateFileName << endl;
            return -1;
          }
      if (BraState.GetVectorDimension() != Space->GetHilbertSpaceDimension())
        {
        cout << "error: vector and Hilbert-space have unequal dimensions"<<endl;
        return -1;
        }

      double Overlap = BraState * State;
      cout << "< " << BraStateFileName << " | H | " << StateFileName << " > "<< " overlap= " << Overlap  << endl; 
      RealVector TmpState(Space->GetHilbertSpaceDimension());
      VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
      Operation.ApplyOperation(Architecture.GetArchitecture());
      double TmpNorm = TmpState * TmpState;

     char* StoreVectorFileName = Manager.GetString("store-vector");
      if (StoreVectorFileName != 0)
       {
         TmpState.WriteVector(StoreVectorFileName);
         cout << "Outputting H |GS>." << endl;           
       }


      double EnergyValue;
      if (BraStateFileName == StateFileName)  
        EnergyValue = State * TmpState;
          else
        EnergyValue = BraState * TmpState;        
      cout << "<Energy>= "<< (EnergyValue) << " Norm= " << TmpNorm << endl;
      return 0;
	 }


      
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [64];
	  sprintf (EigenvectorName, "bosons_%s_n_%d_2s_%d_lz_%d", InteractionName, NbrBosons, LzMax, L);
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
    }

  return 0;
}
