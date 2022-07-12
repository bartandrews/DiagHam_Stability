#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereWithSU3Spin.h"

#include "Hamiltonian/ParticleOnSphereWithSU3SpinGenericHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "GeneralTools/ListIterator.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Options/Options.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);
    
  // some running options and help
  OptionManager Manager ("FQHESphereBosonsWithSU3Spin" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 15);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-tz", "twice the quantum number of the system associated to the Tz generator", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-y", "three time the quantum number of the system associated to the Y generator", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-n1", "number of type 1 particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-n2", "number of type 2 particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-n3", "number of type 3 particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new SingleDoubleOption ('\n', "spin1-flux", "inserted flux for particles with spin 1 (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "spin2-flux", "inserted flux for particles with spin 2 (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "spin3-flux", "inserted flux for particles with spin 3 (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "spin4-flux", "inserted flux for particles with spin 4 (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new  BooleanOption  ('\n', "redundant-kymomenta", "Calculate all subspaces up to Ky  = MaxMomentum-1", false);
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonsWithSU3Spin -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int TotalTz = Manager.GetInteger("total-tz");
  int TotalY = Manager.GetInteger("total-y");
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  if ((Manager.GetInteger("nbr-n1") + Manager.GetInteger("nbr-n2") + Manager.GetInteger("nbr-n3")) == NbrBosons)
    {
      TotalTz = (Manager.GetInteger("nbr-n1") - Manager.GetInteger("nbr-n2"));
      TotalY = (Manager.GetInteger("nbr-n1") + Manager.GetInteger("nbr-n2") - (2 * Manager.GetInteger("nbr-n3")));
    }
  else
    {
      int NbrN1 = (2 * NbrBosons) + TotalY + (3 * TotalTz);
      int NbrN2 = (2 * NbrBosons) + TotalY - (3 * TotalTz);
      int NbrN3 = NbrBosons - TotalY;
      if ((NbrN1 < 0 ) || (NbrN2 < 0 ) || (NbrN3 < 0) || ((NbrN1 % 6) != 0) || ((NbrN2 % 6) != 0) || ((NbrN3 % 3) != 0))
	{
	  cout << "These values of Tz and Y cannot be achieved with this particle number!" << endl;
	  return -1;
	}
      NbrN1 /= 6;
      NbrN2 /= 6;
      NbrN3 /= 3;
    }

  double** PseudoPotentials  = new double*[6];
  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      if (FQHESphereSU3GetPseudopotentials(Manager.GetString("interaction-file"), LzMax, PseudoPotentials) == false)
	return -1;
    }

  char* OutputFileName = new char [512];
  sprintf (OutputFileName, "bosons_sphere_su3_%s_n_%d_2s_%d_tz_%d_y_%d.dat", Manager.GetString("interaction-name"), NbrBosons, LzMax, TotalTz, TotalY);
  ofstream File;
  File.open(OutputFileName, ios::binary | ios::out);
  File.precision(14);

  int Max = LzMax * NbrBosons;
  cout << "maximum Lz value = " << Max << endl;
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
  if (NbrLz > 0)
    {
      if (L + (2 * (NbrLz - 1)) < Max)
	Max = L + (2 * (NbrLz - 1));
    }
  bool FirstRun = true;
  for (; L <= Max; L += 2)
    {
      cout << "----------------------------------------------------------------" << endl;
      BosonOnSphereWithSU3Spin Space (NbrBosons, L, LzMax,  TotalTz, TotalY);	

      Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());	
      AbstractQHEHamiltonian* Hamiltonian = new ParticleOnSphereWithSU3SpinGenericHamiltonian(&Space, NbrBosons, LzMax, PseudoPotentials, 0, 0, 0, Architecture.GetArchitecture(), Memory);
      double Shift = -10.0;
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if ( Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [100];
	  sprintf (EigenvectorName, "bosons_sphere_su3_%s_n_%d_2s_%d_tz_%d_y_%d_lz_%d", Manager.GetString("interaction-name"), NbrBosons, LzMax,
		   TotalTz, TotalY, L);
	}
      
      QHEOnSphereMainTask Task (&Manager, &Space, Hamiltonian, L, Shift, OutputFileName, FirstRun, EigenvectorName, LzMax);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      
      if (FirstRun == true)
	FirstRun = false;
      
      delete Hamiltonian;
    }
  File.close();
  delete[] OutputFileName;
  return 0;
}
