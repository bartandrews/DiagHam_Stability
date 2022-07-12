#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnTorusWithSpin.h"
#include "HilbertSpace/BosonOnTorusShort.h"

#include "Hamiltonian/ParticleOnTorusWithSpinGenericHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Tools/FQHEFiles/FQHETorusPseudopotentialTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "GeneralTools/ListIterator.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "GeneralTools/FilenameTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Options/Options.h"

#include "MainTask/FQHEOnTorusMainTask.h"

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
  OptionManager Manager ("FQHETorusBosonsWithSpinTwoBodyGeneric" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('l', "max-momentum", "maximum momentum for a single particle", 18);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-spin", "total spin of the system", 0);
  (*SystemGroup) += new SingleIntegerOption ('y', "ky-momentum", "constraint on the total momentum modulo the maximum momentum (negative if none)", -1);
  (*SystemGroup) += new SingleDoubleOption ('r', "ratio", "ratio between the two torus lengths", 1.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "spinup-flux", "inserted flux for particles with spin up (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "spindown-flux", "inserted flux for particles with spin down (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");

  (*SystemGroup) += new BooleanOption  ('\n', "redundant-kymomenta", "Calculate all subspaces up to Ky  = MaxMomentum-1", false);
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
  (*ToolsGroup) += new BooleanOption  ('\n', "force-complex", "enforce calculation with complex numbers");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusBosonsWithSpinTwoBodyGeneric -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int TotalSpin = Manager.GetInteger("total-spin");
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int MaxMomentum = Manager.GetInteger("max-momentum");
  int YMomentum = Manager.GetInteger("ky-momentum");
  double XRatio = Manager.GetDouble("ratio");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  if ((TotalSpin & 1) != (NbrBosons & 1))
    {
      TotalSpin &= ~1;
      TotalSpin |= (NbrBosons & 1);
    }

  double** PseudoPotentials  = new double*[3];
  int* NbrPseudoPotentials  = new int[3];
  double** OneBodyPseudoPotentials  = new double*[3];
  double * OneBodyPotentialUpUp = 0;
  double * OneBodyPotentialDownDown = 0;
  double * OneBodyPotentialUpDown = 0;
  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      if (FQHETorusSU2GetPseudopotentials(Manager.GetString("interaction-file"), MaxMomentum, NbrPseudoPotentials, PseudoPotentials, OneBodyPseudoPotentials) == false)
	{
	  return -1;
	}
    }
    
  char* OutputName = new char [512];
  if (OneBodyPseudoPotentials[2] == 0)
    {
      if ((Manager.GetDouble("spinup-flux") == 0.0) && (Manager.GetDouble("spindown-flux") == 0.0))
	sprintf (OutputName, "bosons_torus_su2_kysym_%s_n_%d_2s_%d_sz_%d_ratio_%f.dat", Manager.GetString("interaction-name"), NbrBosons, MaxMomentum, TotalSpin, XRatio);
      else
	sprintf (OutputName, "bosons_torus_su2_kysym_%s_n_%d_2s_%d_sz_%d_ratio_%f_fluxup_%f_fluxdown_%f.dat", Manager.GetString("interaction-name"), NbrBosons, MaxMomentum, TotalSpin, XRatio, Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"));
    }
  else
    {
      if ((Manager.GetDouble("spinup-flux") == 0.0) && (Manager.GetDouble("spindown-flux") == 0.0))
	sprintf (OutputName, "bosons_torus_su2_kysym_%s_n_%d_2s_%d_ratio_%f.dat", Manager.GetString("interaction-name"), NbrBosons, MaxMomentum, XRatio);
      else
	sprintf (OutputName, "bosons_torus_su2_kysym_%s_n_%d_2s_%d_ratio_%f_fluxup_%f_fluxdown_%f.dat", Manager.GetString("interaction-name"), NbrBosons, MaxMomentum, XRatio, Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"));
    }
  ofstream File;
  File.open(OutputName, ios::binary | ios::out);
  File.precision(14);

  int MomentumModulo = FindGCD(NbrBosons, MaxMomentum);
  int YMaxMomentum;
  if (Manager.GetBoolean("redundant-kymomenta"))
    YMaxMomentum = (MaxMomentum - 1);
  else
    YMaxMomentum = (MomentumModulo - 1);
  if (YMomentum < 0)
    YMomentum = 0;
  else
    YMaxMomentum = YMomentum; 

  bool FirstRun = true;
  for (int YMomentum2 = YMomentum; YMomentum2 <= YMaxMomentum; ++YMomentum2)
    {
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;
      BosonOnTorusWithSpin * Space = 0;
      if (OneBodyPseudoPotentials[2] == 0)
	Space = new BosonOnTorusWithSpin (NbrBosons, MaxMomentum, TotalSpin, YMomentum2);
      else
	Space = new BosonOnTorusWithSpin (NbrBosons, MaxMomentum, YMomentum2);
	

      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      AbstractQHEHamiltonian* Hamiltonian = new ParticleOnTorusWithSpinGenericHamiltonian(Space, NbrBosons, MaxMomentum, XRatio,
											  NbrPseudoPotentials[0], PseudoPotentials[0],
											  NbrPseudoPotentials[1], PseudoPotentials[1],
											  NbrPseudoPotentials[2], PseudoPotentials[2],
											  Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"),
											  Architecture.GetArchitecture(), Memory, 0, OneBodyPseudoPotentials[0], 
											  OneBodyPseudoPotentials[1], OneBodyPseudoPotentials[2]);

      double Shift = -10.0;
      if (((BooleanOption*) Manager["force-complex"])->GetBoolean() == true)
        Shift = 0.0;

      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if ( Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [100];
	  char* TmpName = RemoveExtensionFromFileName(OutputName, ".dat");
	  sprintf (EigenvectorName, "%s_ky_%d", TmpName, YMomentum);
	}
      
      FQHEOnTorusMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, YMomentum2, Shift, OutputName, FirstRun, EigenvectorName);
      if (((BooleanOption*) Manager["force-complex"])->GetBoolean() == true)
	 Task.ForceComplex();
      	
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      
      if (FirstRun == true)
	FirstRun = false;
      
      delete Hamiltonian;
      delete Space;
    }
  File.close();
  delete[] OutputName;
  return 0;
}
