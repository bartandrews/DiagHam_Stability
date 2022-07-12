#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereLong.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOn4DSphere.h"
#include "HilbertSpace/BosonOn4DSphereLong.h"

#include "Hamiltonian/ParticleOn4DSphereDeltaHamiltonian.h"
#include "Hamiltonian/ParticleOn4DSphereThreeBodyDeltaHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

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
  OptionManager Manager ("FQHESphere4DBosonsDelta" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "number of flux quanta", 20);
  (*SystemGroup) += new BooleanOption ('\n', "three-body", "use three body delta interaction instead of two-body interaction");
  /*
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");*/
//   (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
//   (*SystemGroup) += new BooleanOption  ('\n', "get-lvalue", "compute mean l value from <L^2> for each eigenvalue");
//   (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-jz", " only evaluate one jz sector (negative if all sectors have to be computed) ", 1000);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kz", "only evaluate one kz sector (negative if all sectors have to be computed)", 1000);

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
      cout << "see man page for option syntax or type FQHESphere4DBosonsDelta -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int NbrFluxQuanta = Manager.GetInteger("nbr-flux");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool FirstRun = true;
  bool ThreeBodyFlag = Manager.GetBoolean("three-body");
  

  char* OutputName = new char [256];
  if (ThreeBodyFlag == false)
    sprintf (OutputName, "bosons_sphere4d_delta_n_%d_2s_%d.dat", NbrBosons, NbrFluxQuanta);
  else
    sprintf (OutputName, "bosons_sphere4d_threebody_delta_n_%d_2s_%d.dat", NbrBosons, NbrFluxQuanta);
  
  int MinJz = 0;
//   int MinJz = -NbrFluxQuanta*NbrBosons;
  int MaxJz = NbrFluxQuanta*NbrBosons;
  int NbrOrbitals = (NbrFluxQuanta + 1)*(NbrFluxQuanta + 2)*(NbrFluxQuanta + 3) / 6;
  
  if (Manager.GetInteger("only-jz") != 1000)
  {
   MinJz = Manager.GetInteger("only-jz");
   MaxJz = Manager.GetInteger("only-jz");
  }
    
  for (int jz = MinJz; jz <= MaxJz; ++jz)
    {
      int MinKz = 0;
//       int MinKz = -NbrFluxQuanta*NbrBosons;
      int MaxKz = NbrFluxQuanta*NbrBosons - jz;
//       int MaxKz = NbrFluxQuanta*NbrBosons;
      if (Manager.GetInteger("only-kz") != 1000)
	{
	  MinKz = Manager.GetInteger("only-kz");
	  MaxKz = Manager.GetInteger("only-kz");
	  if (((Manager.GetInteger("only-jz") + Manager.GetInteger("only-kz")) & 1 ) != ((NbrBosons * NbrFluxQuanta) & 1) or ((Manager.GetInteger("only-jz") + Manager.GetInteger("only-kz")) > NbrBosons * NbrFluxQuanta))
	    cout << "Incompatible values for jz, kz, nbr of particles and number of flux quanta" << endl;
	  if (Manager.GetInteger("only-jz") < Manager.GetInteger("only-kz"))
	    cout << "jz should be bigger than kz" << endl;
	}
      
      for (int kz = MinKz; kz <= MaxKz ; ++kz)
	{
	  if ((((jz + kz) & 1) == ((NbrBosons * NbrFluxQuanta) & 1)) && (abs(kz)<=abs(jz)) && (abs(kz)+abs(jz)<= (NbrBosons*NbrFluxQuanta)))
	    {
	      // 	   
	      cout << "(jz,kz) = (" << jz << "," << kz << ")" << endl; 
	      ParticleOnSphere* Space = 0;
	      if (NbrOrbitals + NbrBosons < 65)
		Space = new BosonOn4DSphere(NbrBosons, NbrFluxQuanta, jz, kz);
	      else
		Space = new BosonOn4DSphereLong(NbrBosons, NbrFluxQuanta, jz, kz);
	      
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      //       for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	      // 	Space->PrintState(cout, i);
	      
	      AbstractQHEHamiltonian* Hamiltonian = 0;
	      
	      if (ThreeBodyFlag == false)
		{
		  Hamiltonian = new ParticleOn4DSphereDeltaHamiltonian(Space, NbrBosons, NbrFluxQuanta,
								       Architecture.GetArchitecture(), 
								       Memory);
		}
	      else
		{
		  Hamiltonian = new ParticleOn4DSphereThreeBodyDeltaHamiltonian(Space, NbrBosons, NbrFluxQuanta, 0, Architecture.GetArchitecture(),  Memory, DiskCacheFlag, LoadPrecalculationFileName);
		}
	      //        double Shift = - 0.5 * ((double) (NbrBosons * NbrBosons)) / (0.5 * ((double) NbrFluxQuanta));
	      //        Hamiltonian->ShiftHamiltonian(Shift);
	      char* EigenvectorName = 0;
	      if (Manager.GetBoolean("eigenstate") == true)	
		{
		  EigenvectorName = new char [64];
		  if (ThreeBodyFlag == false)
		    sprintf (EigenvectorName, "bosons_sphere4d_delta_n_%d_2s_%d_jz_%d_kz_%d", NbrBosons, NbrFluxQuanta, jz, kz);
		  else 
		    sprintf (EigenvectorName, "bosons_sphere4d_threebody_delta_n_%d_2s_%d_jz_%d_kz_%d", NbrBosons, NbrFluxQuanta, jz, kz);
		}
	      
	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d", jz, kz);
	      
	      char* SubspaceLegend = new char[256];
	      sprintf (SubspaceLegend, "jz kz");
	      GenericRealMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, ContentPrefix, SubspaceLegend, 0, OutputName, FirstRun, EigenvectorName);
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
	  
	}
    }
  return 0;
}
