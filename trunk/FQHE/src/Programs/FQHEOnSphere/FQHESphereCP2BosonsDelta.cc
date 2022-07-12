#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereLong.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnCP2.h"
#include "HilbertSpace/BosonOnCP2TzSymmetry.h"
#include "HilbertSpace/BosonOnCP2TzZ3Symmetry.h"


#include "Hamiltonian/ParticleOnCP2DeltaHamiltonian.h"
#include "Hamiltonian/ParticleOnCP2ThreeBodyDeltaHamiltonian.h"

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
  OptionManager Manager ("FQHESphereCP2BosonsDelta" , "0.01");
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
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
//   (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
//   (*SystemGroup) += new BooleanOption  ('\n', "get-lvalue", "compute mean l value from <L^2> for each eigenvalue");
//   (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-tz", " only evaluate one jz sector (negative if all sectors have to be computed) ", 1000);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-y", "only evaluate one kz sector (negative if all sectors have to be computed)", 1000);
  (*SystemGroup) += new BooleanOption  ('\n', "tzsymmetrized-basis", "use Tz <-> -Tz symmetrized version of the basis (only valid if total-tz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "tzZ3symmetrized-basis", "use Tz <-> -Tz and Z3 permutations symmetrized version of the basis (only valid if total-tz=0 and total-y = 0)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-tzparity", "select the  Tz <-> -Tz symmetric sector with negative parity");
  (*SystemGroup) += new  SingleStringOption ('\n', "onebody-file", "file describing the one body potential in terms of pseudopotentials",0);
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
      cout << "see man page for option syntax or type FQHESphereCP2BosonsDelta -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int NbrFluxQuanta = Manager.GetInteger("nbr-flux");
  int NbrOrbitals = (NbrFluxQuanta + 1)*(NbrFluxQuanta + 2)/ 2;
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool FirstRun = true;
  bool ThreeBodyFlag = Manager.GetBoolean("three-body");
  bool TzSymmetrizedBasis = Manager.GetBoolean("tzsymmetrized-basis");
  bool TzZ3SymmetrizedBasis = Manager.GetBoolean("tzZ3symmetrized-basis");
  bool TzMinusParity = Manager.GetBoolean("minus-tzparity");
  double* OneBodyPotentials = 0;
  if (Manager.GetString("onebody-file") != 0)
  {
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(Manager.GetString("onebody-file")) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return -1;
     }
  int TmpNbrPseudoPotentials;
  if (InteractionDefinition.GetAsDoubleArray("Onebodypotentials", ' ', OneBodyPotentials, TmpNbrPseudoPotentials) == true)
    {
      if (TmpNbrPseudoPotentials != (NbrOrbitals))
	{
	   cout << "Invalid number of pseudo-potentials" << endl;
	   return -1;	  
	}
    }
  }
  char* OutputName = new char [256];
  if ((ThreeBodyFlag == false) && (TzSymmetrizedBasis == false) && (TzZ3SymmetrizedBasis == false))
    sprintf (OutputName, "bosons_cp2_delta_n_%d_2s_%d.dat", NbrBosons, NbrFluxQuanta);
  else
  {
    if (ThreeBodyFlag == false)
    {
      if (TzSymmetrizedBasis == false && TzZ3SymmetrizedBasis == false)
	sprintf (OutputName, "bosons_cp2_delta_n_%d_2s_%d.dat", NbrBosons, NbrFluxQuanta);
      if (TzSymmetrizedBasis == true && TzMinusParity == false)
	sprintf (OutputName, "bosons_cp2_tzpsym_delta_n_%d_2s_%d.dat", NbrBosons, NbrFluxQuanta);
      if (TzSymmetrizedBasis == true && TzMinusParity == true)
	sprintf (OutputName, "bosons_cp2_tzmsym_delta_n_%d_2s_%d.dat", NbrBosons, NbrFluxQuanta);
      if (TzZ3SymmetrizedBasis == true)
	sprintf (OutputName, "bosons_cp2_tzZ3sym_delta_n_%d_2s_%d.dat", NbrBosons, NbrFluxQuanta);
    }
    else
      if (TzSymmetrizedBasis == false && TzZ3SymmetrizedBasis == false)
	sprintf (OutputName, "bosons_cp2_threebody_delta_n_%d_2s_%d.dat", NbrBosons, NbrFluxQuanta);
      if (TzSymmetrizedBasis == true && TzMinusParity == false)
	sprintf (OutputName, "bosons_cp2_tzpsym_threebody_delta_n_%d_2s_%d.dat", NbrBosons, NbrFluxQuanta);
      if (TzSymmetrizedBasis == true && TzMinusParity == true)
	sprintf (OutputName, "bosons_cp2_tzmsym_threebody_delta_n_%d_2s_%d.dat", NbrBosons, NbrFluxQuanta);
      if (TzZ3SymmetrizedBasis == true)
	sprintf (OutputName, "bosons_cp2_tzZ3sym_threebody_delta_n_%d_2s_%d.dat", NbrBosons, NbrFluxQuanta);
  }
  
  if ((Manager.GetInteger("only-tz") == 1000) && (Manager.GetInteger("only-y") == 1000))
  {
    int MinR = 0;
    int MaxR = NbrFluxQuanta*NbrBosons;
  
    for (int r = MinR; r <= MaxR; ++r)
      {
	int MinS = 0;
	int MaxS = NbrFluxQuanta*NbrBosons - r;
	if (MaxS > r)
	  MaxS = r;
	for (int s = MinS; s <= MaxS ; ++s)
	{
	  int tz = r - s;
	  int y = 3*(r + s) - 2*NbrBosons*NbrFluxQuanta;
	  cout << "(tz,y) = (" << tz << "," << y << ")" << endl; 
	  if (TzSymmetrizedBasis == true)
	      {
		cout << "  Tz symmetrized basis ";
		if (TzMinusParity == true)
		  cout << "(minus parity) " << endl;
		else
		  cout << "(plus parity) " << endl;
	      }
	  ParticleOnSphere* Space = 0;
	  if (NbrOrbitals + NbrBosons < 65)
	  {
	    if (TzSymmetrizedBasis == false)
	      Space = new BosonOnCP2(NbrBosons, NbrFluxQuanta, tz, y);
	    else
	    {
	     if (tz != 0)
	       Space = new BosonOnCP2(NbrBosons, NbrFluxQuanta, tz, y);
	     else
	      Space = new BosonOnCP2TzSymmetry(NbrBosons, NbrFluxQuanta, tz, y, TzMinusParity);
	    }
	    
	  }
	  else
	    cout << " Warning : number of orbitals too big " << endl;
     
	Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	  Memory = Architecture.GetArchitecture()->GetLocalMemory();
//       for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
// 	Space->PrintState(cout, i);
       
	AbstractQHEHamiltonian* Hamiltonian = 0;
	if (ThreeBodyFlag == false)
	{
	  if (OneBodyPotentials == 0)
	    Hamiltonian = new ParticleOnCP2DeltaHamiltonian(Space, NbrBosons, NbrFluxQuanta, Architecture.GetArchitecture(), Memory);
	  else
	    Hamiltonian = new ParticleOnCP2DeltaHamiltonian(Space, NbrBosons, NbrFluxQuanta, OneBodyPotentials, Architecture.GetArchitecture(), Memory);
	}
	else
	{
	  Hamiltonian = new ParticleOnCP2ThreeBodyDeltaHamiltonian(Space, NbrBosons, NbrFluxQuanta, 0, Architecture.GetArchitecture(), Memory, DiskCacheFlag, LoadPrecalculationFileName);
	}
	
      
//        double Shift = - 0.5 * ((double) (NbrBosons * NbrBosons)) / (0.5 * ((double) NbrFluxQuanta));
//        Hamiltonian->ShiftHamiltonian(Shift);
	char* EigenvectorName = 0;
	if (Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [64];
	  if (ThreeBodyFlag == false)
	    if (TzSymmetrizedBasis == false && TzZ3SymmetrizedBasis == false)
	      sprintf (EigenvectorName, "bosons_cp2_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	    else
	    {
	      if(TzSymmetrizedBasis == true && TzMinusParity == false)
		sprintf (EigenvectorName, "bosons_cp2_tzpsym_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	      if(TzSymmetrizedBasis == true && TzMinusParity == true)
		sprintf (EigenvectorName, "bosons_cp2_tzmsym_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	      if (TzZ3SymmetrizedBasis == true)
		sprintf (EigenvectorName, "bosons_cp2_tzZ3sym_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	    }
	  else
	    if (TzSymmetrizedBasis == false && TzZ3SymmetrizedBasis == false)
	    sprintf (EigenvectorName, "bosons_cp2_threebody_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	  else
	    {
	      if(TzSymmetrizedBasis == true && TzMinusParity == false)
		sprintf (EigenvectorName, "bosons_cp2_tzpsym_threebody_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	      if(TzSymmetrizedBasis == true && TzMinusParity == true)
		sprintf (EigenvectorName, "bosons_cp2_tzmsym_threebody_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	      if (TzZ3SymmetrizedBasis == true)
		sprintf (EigenvectorName, "bosons_cp2_tzZ3sym_threebody_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	    }
	}
	
	char* ContentPrefix = new char[256];
	sprintf (ContentPrefix, "%d %d", tz, y);
      
	char* SubspaceLegend = new char[256];
	sprintf (SubspaceLegend, "tz y");
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
    else
    {
     int tz = Manager.GetInteger("only-tz");
     int y = Manager.GetInteger("only-y");
     if (((y + 3*tz + 2*NbrBosons*NbrFluxQuanta) % 6 != 0) || ((y - 3*tz + 2*NbrBosons*NbrFluxQuanta) % 6 != 0))
     {
       cout << "Y + 3Tz + 2N*Nphi and Y - 3Tz + 2N*Nphi should multiple of 6" << endl;
       return -1;
     }
     cout << "(tz,y) = (" << tz << "," << y << ")" << endl; 
     if (TzSymmetrizedBasis == true)
      {
	cout << "  Tz symmetrized basis ";
	if (TzMinusParity == true)
	  cout << "(minus parity) ";
	else
	  cout << "(plus parity) ";
      }
      if (TzZ3SymmetrizedBasis == true)
      {
	cout << "  Tz and Z3 symmetrized basis ";
	if (TzMinusParity == true)
	  cout << "(minus parity) ";
	else
	  cout << "(plus parity) ";
      }
     ParticleOnSphere* Space = 0;
     if (NbrOrbitals + NbrBosons < 65)
	{
	  if (TzSymmetrizedBasis == false && TzZ3SymmetrizedBasis == false)
	    Space = new BosonOnCP2(NbrBosons, NbrFluxQuanta, tz, y);
	  else
	  {
	    if (TzSymmetrizedBasis == true)
	    {
	      if (tz != 0)
	      {
		cout << "tzsymmetrized-basis mode only valid for tz = 0" << endl;
		return 0;
	      }
	      else
		Space = new BosonOnCP2TzSymmetry(NbrBosons, NbrFluxQuanta, tz, y, TzMinusParity);
	      }
	    if (TzZ3SymmetrizedBasis == true)
	    {
	      if ( tz == 0 && y ==0)
		Space = new BosonOnCP2TzZ3Symmetry(NbrBosons, NbrFluxQuanta, tz, y, TzMinusParity);
	      else
	      {
		cout << "tzZ3symmetrized-basis mode only valid for tz = 0" << endl;
		return 0;
	      }
	    }
	  }  
	 }
     else
	cout << " Warning : number of orbitals too big " << endl;
     Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
    if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
//       for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
// 	Space->PrintState(cout, i);
       
    AbstractQHEHamiltonian* Hamiltonian = 0;
    if (ThreeBodyFlag == false)
    {
      if (OneBodyPotentials == 0)
	Hamiltonian = new ParticleOnCP2DeltaHamiltonian(Space, NbrBosons, NbrFluxQuanta, Architecture.GetArchitecture(), Memory);
      else
	Hamiltonian = new ParticleOnCP2DeltaHamiltonian(Space, NbrBosons, NbrFluxQuanta, OneBodyPotentials, Architecture.GetArchitecture(), Memory);
    }
    else
	Hamiltonian = new ParticleOnCP2ThreeBodyDeltaHamiltonian(Space, NbrBosons, NbrFluxQuanta, 0, Architecture.GetArchitecture(),  Memory, DiskCacheFlag, LoadPrecalculationFileName);
      
//        double Shift = - 0.5 * ((double) (NbrBosons * NbrBosons)) / (0.5 * ((double) NbrFluxQuanta));
//        Hamiltonian->ShiftHamiltonian(Shift);
    char* EigenvectorName = 0;
    if (Manager.GetBoolean("eigenstate") == true)	
      {
	EigenvectorName = new char [64];
	  if (ThreeBodyFlag == false)
	    if (TzSymmetrizedBasis == false && TzZ3SymmetrizedBasis == false)
	      sprintf (EigenvectorName, "bosons_cp2_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	    else
	    {
	      if(TzSymmetrizedBasis == true && TzMinusParity == false)
		sprintf (EigenvectorName, "bosons_cp2_tzpsym_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	      if(TzSymmetrizedBasis == true && TzMinusParity == true)
		sprintf (EigenvectorName, "bosons_cp2_tzmsym_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	      if (TzZ3SymmetrizedBasis == true)
		sprintf (EigenvectorName, "bosons_cp2_tzZ3sym_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	    }
	else
	  if (TzSymmetrizedBasis == false && TzZ3SymmetrizedBasis == false)
	   sprintf (EigenvectorName, "bosons_cp2_threebody_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	  else
	    {
	      if(TzSymmetrizedBasis == true && TzMinusParity == false)
		sprintf (EigenvectorName, "bosons_cp2_tzpsym_threebody_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	      if(TzSymmetrizedBasis == true && TzMinusParity == true)
		sprintf (EigenvectorName, "bosons_cp2_tzmsym_threebody_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	      if (TzZ3SymmetrizedBasis == true)
		sprintf (EigenvectorName, "bosons_cp2_tzZ3sym_threebody_delta_n_%d_2s_%d_tz_%d_y_%d", NbrBosons, NbrFluxQuanta, tz, y);
	    }
      }
    char* ContentPrefix = new char[256];
    sprintf (ContentPrefix, "%d %d", tz, y);
      
    char* SubspaceLegend = new char[256];
    sprintf (SubspaceLegend, "tz y");
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

  return 0;
}
