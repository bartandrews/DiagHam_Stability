#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereWithSU4Spin.h"

#include "Hamiltonian/ParticleOnSphereWithSU4SpinGenericHamiltonian.h"

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
  OptionManager Manager ("FQHESphereBosonsWithSU4Spin" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('i', "total-isosz", "twice the z component of the total isospin (i.e valley SU(2) degeneracy) of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('e', "total-entanglement", "twice the projection of the total spin-isopsin entanglement of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-n1", "number of up-plus particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-n2", "number of up-minus particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-n3", "number of down-plus particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-n4", "number of down-minus particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
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
      cout << "see man page for option syntax or type FQHESphereBosonsWithSU4Spin -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int TotalSz = Manager.GetInteger("total-sz");
  int TotalIz = Manager.GetInteger("total-isosz");
  int TotalPz = Manager.GetInteger("total-entanglement");
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  if ((Manager.GetInteger("nbr-n1") + Manager.GetInteger("nbr-n2") + Manager.GetInteger("nbr-n3") + Manager.GetInteger("nbr-n4")) == NbrBosons)
    {
      TotalSz = (Manager.GetInteger("nbr-n1") + Manager.GetInteger("nbr-n2")) - (Manager.GetInteger("nbr-n3") + Manager.GetInteger("nbr-n4"));
      TotalIz = (Manager.GetInteger("nbr-n1") + Manager.GetInteger("nbr-n3")) - (Manager.GetInteger("nbr-n2") + Manager.GetInteger("nbr-n4"));
      TotalPz = (Manager.GetInteger("nbr-n1") + Manager.GetInteger("nbr-n4")) - (Manager.GetInteger("nbr-n2") + Manager.GetInteger("nbr-n3"));
    }
  else
    {
      int NbrNUpPlus = (NbrBosons + TotalSz + TotalIz + TotalPz);
      int NbrNUpMinus = (NbrBosons + TotalSz - TotalIz - TotalPz);
      int NbrNDownPlus = (NbrBosons - TotalSz + TotalIz - TotalPz);
      int NbrNDownMinus = (NbrBosons - TotalSz - TotalIz + TotalPz);
      if ((NbrNUpPlus < 0 ) || (NbrNUpMinus < 0 ) || (NbrNDownPlus < 0) || (NbrNDownMinus < 0) || 
	  ((NbrNUpPlus & 3) != 0) || ((NbrNUpMinus & 3) != 0) || ((NbrNDownPlus & 3) != 0) || ((NbrNDownMinus & 3) != 0))
	{
	  cout << "These values of Sz, Iz and Pz cannot be achieved with this particle number!" << endl;
	  return -1;
	}
      NbrNUpPlus >>= 2;
      NbrNUpMinus >>= 2;
      NbrNDownPlus >>= 2;
      NbrNDownMinus >>= 2;
    }

  double** PseudoPotentials  = new double*[10];
  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      if (FQHESphereSU4GetPseudopotentials(Manager.GetString("interaction-file"), LzMax, PseudoPotentials) == false)
	return -1;
    }

  char* OutputFileName = new char [512];
  sprintf (OutputFileName, "bosons_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_pz_%d.dat", Manager.GetString("interaction-name"), NbrBosons, LzMax, TotalSz, TotalIz, TotalPz);
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
      BosonOnSphereWithSU4Spin Space (NbrBosons, L, LzMax,  TotalSz, TotalIz, TotalPz);	

      Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());	
      AbstractQHEHamiltonian* Hamiltonian = new ParticleOnSphereWithSU4SpinGenericHamiltonian(&Space, NbrBosons, LzMax, PseudoPotentials, 0.0, Architecture.GetArchitecture(), Memory);
      double Shift = -10.0;
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if ( Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [100];
	  sprintf (EigenvectorName, "bosons_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_pz_%d_lz_%d", Manager.GetString("interaction-name"), NbrBosons, LzMax,
		   TotalSz, TotalIz, TotalPz, L);
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
