#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereTwoLandauLevels.h"

#include "Hamiltonian/ParticleOnSphereGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereTwoLandauLevelDeltaHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/QHEOnSphereMainTask.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

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
  OptionManager Manager ("FQHESphereFermionsTwoLandauLevelTwoBodyGeneric" , "0.01");
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

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new BooleanOption  ('\n', "restrict-polarization", "restrict number of particles in each Landau level (provided by nbrparticles-up and nbrparticles-down)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrparticles-up", "number of particles in N=1 LL", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrparticles-down", "number of particles in N=0 LL", 0);  
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "number of flux quanta", 6);
  (*SystemGroup) += new SingleIntegerOption  ('z', "initial-lz", "initial lz value", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "potential-lll", "one particle potential lll", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "potential-fll", "one particle potential sll", 0.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potentials");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");
  (*SystemGroup) += new BooleanOption  ('\n', "show-interaction-factors", "show the interaction factors being used");
  
  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
#endif
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");  

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsTwoLandauLevelTwoBodyGeneric -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
   {
      Manager.DisplayHelp (cout);
      return 0;
   }
  
  // Get options for the options manager.  
  bool GroundFlag = Manager.GetBoolean("ground");
  bool PolarizationFlag = Manager.GetBoolean("restrict-polarization");  
  int NbrFermions = Manager.GetInteger("nbr-particles");
  int NbrFermionsUp = Manager.GetInteger("nbrparticles-up");
  int NbrFermionsDown = Manager.GetInteger("nbrparticles-down");
  int NbrFluxQuanta = Manager.GetInteger("nbr-flux");
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  long Memory = Manager.GetInteger("memory") << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  double cyclotron_energies[2];
  cyclotron_energies[0] = (double)Manager.GetDouble("potential-lll");
  cyclotron_energies[1] = (double)Manager.GetDouble("potential-fll");
  bool ShowIntFactorsFlag = Manager.GetBoolean("show-interaction-factors");
  
  
  bool FirstRun = true;
  
  // Read in Pseudopotentials from the interaction file. 
  double** PseudoPotentials ;
  PseudoPotentials = new double*[10];
  for ( int i = 0 ; i < 10;  i++ ) {
      PseudoPotentials[i] = new double[NbrFluxQuanta+3]; //make enough space for each level and set all to zero
      for ( int j = 0 ; j < NbrFluxQuanta+3 ; j++ ) 
        {
	  PseudoPotentials[i][j] = 0.0;
        }
  }
  
  if ( Manager.GetString("interaction-file") == 0 ) 
    {
      cout << "No interaction file supplied. Will use delta interaction." << endl;
    }
  else
    {
      if ( Manager.GetString("interaction-file") != 0 )
        {
	  if ( !FQHESphereTwoLandauLevelGetPseudopotentials(Manager.GetString("interaction-file"), NbrFluxQuanta, PseudoPotentials) ) 
	    {
	      cout << "There were problems encountered when attempting to read the pseudopotential file: " << Manager.GetString("interaction-file") << endl;
	      return -1;
	    }	  
	}
    }	

  int Max;
  char* OutputNameLz = new char [256 + strlen(Manager.GetString("interaction-name"))];
  sprintf (OutputNameLz, "fermions_%s_n_%d_2s_%d_lz.dat", Manager.GetString("interaction-name"), NbrFermions, NbrFluxQuanta);
  if ( (NbrFermions % 2) == 1 )
    {
    	Max = NbrFluxQuanta+2;
	for ( int i = 0 ; i < (NbrFermions - 1) ; i+=2 )
	  Max += 2*(NbrFluxQuanta - i);
    } 
  else 
    {
	Max = NbrFluxQuanta+2;
	int i;
	for ( i = 2 ; i < (NbrFermions - 1) ; i+=2 )
	  Max += 2*(NbrFluxQuanta + 2 - i);    	
    	Max += (NbrFluxQuanta + 2 - NbrFermions);
    }
  
  int  L = InitialLz;

  if (L < -Max)
    L = -Max;
  else
    if (L > Max)
      L = Max;
  if ((abs(Max) & 1) != (abs(InitialLz) & 1))
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
      cout << "Lz="<<L<<endl;
      // create space.
      ParticleOnSphereWithSpin* Space = 0;
      if (PolarizationFlag) 
          Space = new FermionOnSphereTwoLandauLevels(NbrFermionsUp, NbrFermionsDown, L, NbrFluxQuanta+2, NbrFluxQuanta);
      else
          Space = new FermionOnSphereTwoLandauLevels(NbrFermions, L, NbrFluxQuanta+2, NbrFluxQuanta);
              
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
      
    //for ( int i = 0 ; i < Space->GetHilbertSpaceDimension() ; i++ ) 
	  // {
	  //   cout << i << ": " ;
	  //   Space->PrintState(cout, i);
	  //   cout << endl;
	  //  }
    // exit(2); 
      
      // Create data structure for the Hamiltonian.
      AbstractQHEOnSphereHamiltonian* Hamiltonian = 0;
      
      if ( Manager.GetString("interaction-file") == 0 ) 
	{
	  Hamiltonian = new ParticleOnSphereTwoLandauLevelDeltaHamiltonian(Space, NbrFermions, NbrFluxQuanta+2, NULL, cyclotron_energies,  
									   Architecture.GetArchitecture(), 
									   Memory, DiskCacheFlag,
									   LoadPrecalculationFileName, ShowIntFactorsFlag);
	}
      else 
	{
	  Hamiltonian = new ParticleOnSphereTwoLandauLevelDeltaHamiltonian(Space, NbrFermions, NbrFluxQuanta+2, PseudoPotentials, cyclotron_energies,  
									   Architecture.GetArchitecture(), 
									   Memory, DiskCacheFlag,
									   LoadPrecalculationFileName, ShowIntFactorsFlag);
	}

       double Shift = 0;

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
          cout << "< Energy > = "<<EnergyValue<<endl;
          cout << "< shifted energy > = "<<EnergyValue + Shift<<endl;
          return 0;
        }

      int LzMax = NbrFluxQuanta + 2;
      
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [128];
	  sprintf (EigenvectorName, "fermions_2ll_%s_n_%d_2s_%d_lz_%d", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrFermions, NbrFluxQuanta, L);
	}
      
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
      MainTaskOperation TaskOperation (&Task); TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      delete Hamiltonian;
      delete Space;
      
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      
      if (FirstRun == true)
	FirstRun = false;
    }
    
  for ( int i = 0 ; i < 10;  i++ ) 
    delete [] PseudoPotentials[i];
  delete []PseudoPotentials;
  delete []OutputNameLz;
  

  return 0;
}
