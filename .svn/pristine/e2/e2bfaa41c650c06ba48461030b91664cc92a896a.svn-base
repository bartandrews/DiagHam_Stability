#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinSzSymmetry.h"

#include "Hamiltonian/ParticleOnCylinderWithSpinGenericHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Tools/FQHEFiles/FQHETorusPseudopotentialTools.h"
#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

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


  OptionManager Manager ("FQHECylinderBosonsWithSpin" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);
  ParticleOnSphereManager ParticleManager(false, true, 2);
  ParticleManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleDoubleOption  ('r', "aspect-ratio", "aspect ratio of the cylinder", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "cylinder-perimeter", "if non zero, fix the cylinder perimeter (in magnetic length unit) instead of the aspect ratio", 0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "spinup-flux", "inserted flux for particles with spin up (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "spindown-flux", "inserted flux for particles with spin down (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "s2-factor", "multiplicative factor in front of an optional S^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new BooleanOption ('\n', "l2-s2-only", "compose Hamiltonian only of L2 and S2 terms");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "negative-lz", "calculate for negative, instead of positive Lz");
  (*SystemGroup) += new BooleanOption  ('\n', "all-lz", "calculate both negative and positive Lz");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "laughlin-exponent", "start the Haldane algorithm from Laughlin state with exponent m)", -1);
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleDoubleOption ('\n', "energy-shift", "apply a temporary energy shift used during the diagonalization", 0.0);

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "s2-memory", "amount of memory that can be allocated for fast multiplication of s2 term (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "l2-memory", "amount of memory that can be allocated for fast multiplication of l2 term (in Mbytes)", 500);
  (*PrecalculationGroup) += new BooleanOption  ('\n', "allow-disk-storage", "expand memory for fast multiplication using disk storage",false);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-space", "show detail of states in the Hilbert-space");
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the Hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECylinderBosonsWithSpin -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalSpin = Manager.GetInteger("total-sz");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  double Ratio = Manager.GetDouble("aspect-ratio");
  double Perimeter = Manager.GetDouble("cylinder-perimeter");
  if (Perimeter != 0.0)
    {
      Ratio = 2.0 * M_PI * (LzMax + 1) / (Perimeter * Perimeter);
    }
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;  
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
  bool onDiskCacheFlag = Manager.GetBoolean("allow-disk-storage");
  bool FirstRun = true;

  int NbrUp = (NbrBosons + TotalSpin) >> 1;
  int NbrDown = (NbrBosons - TotalSpin) >> 1;
  if ((NbrUp < 0) || (NbrDown < 0))
    {
      cout << "This value of the spin z projection cannot be achieved with this particle number!" << endl;
      return -1;
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
      if (FQHETorusSU2GetPseudopotentials(Manager.GetString("interaction-file"), LzMax + 1, NbrPseudoPotentials, PseudoPotentials, OneBodyPseudoPotentials) == false)
	{
	  return -1;
	}
    }

  char* DiscreteSymmetryName = new char[128];
  if (Manager.GetBoolean("szsymmetrized-basis") == false)
    {
      if (Manager.GetBoolean("lzsymmetrized-basis") == false)
	{
	  sprintf (DiscreteSymmetryName, "");
	}
      else
	{
	  if (Manager.GetBoolean("minus-lzparity") == false)
	    {
	      sprintf (DiscreteSymmetryName, "_lzsym_1");
	    }
	  else
	    {
	      sprintf (DiscreteSymmetryName, "_lzsym_-1");
	    }
	}
    }
  else
    {
      if (Manager.GetBoolean("lzsymmetrized-basis") == false)
	{
	  if (Manager.GetBoolean("minus-szparity") == false)
	    {
	      sprintf (DiscreteSymmetryName, "_szsym_1");
	    }
	  else
	    {
	      sprintf (DiscreteSymmetryName, "_szsym_-1");
	    }
	}
      else
	{
	  if (Manager.GetBoolean("minus-szparity") == false)
	    {
	      if (Manager.GetBoolean("minus-lzparity") == false)
		{
		  sprintf (DiscreteSymmetryName, "_lzsym_1_szsym_1");
		}
	      else
		{
		  sprintf (DiscreteSymmetryName, "_lzsym_-1_szsym_1");
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("minus-lzparity") == false)
		{
		  sprintf (DiscreteSymmetryName, "_lzsym_1_szsym_-1");
		}
	      else
		{
		  sprintf (DiscreteSymmetryName, "_lzsym_-1_szsym_-1");
		}
	    }
	}
    }
  char* GeometryName = new char[256];
  if (Manager.GetInteger("nbrspin-polarized") == 0)
    {
      if (Perimeter > 0.0)	
	{
	  sprintf (GeometryName, "cylinder_perimeter_%.6f_su2", Perimeter);
	}
      else
	{
	  sprintf (GeometryName, "cylinder_ratio_%.6f_su2", Ratio);
	}
    }
  else
    {
      if (Perimeter > 0.0)	
	{
	  sprintf (GeometryName, "cylinder_perimeter_%.6f_su2_polarized_%ld", Perimeter, Manager.GetInteger("nbrspin-polarized"));
	}
      else
	{
	  sprintf (GeometryName, "cylinder_ratio_%.6f_su2_polarized_%ld", Ratio, Manager.GetInteger("nbrspin-polarized"));
	}
    }
  char* OutputName = new char [512 + strlen(DiscreteSymmetryName) + strlen(GeometryName )+ strlen(Manager.GetString("interaction-name"))];
  if (OneBodyPseudoPotentials[2] == 0)
    {
      if ((Manager.GetDouble("spinup-flux") == 0.0) && (Manager.GetDouble("spindown-flux") == 0.0))
	{
	  sprintf (OutputName, "bosons_%s%s_%s_n_%d_2s_%d_sz_%d.dat", GeometryName, DiscreteSymmetryName, Manager.GetString("interaction-name"), 
		   NbrBosons, LzMax, TotalSpin);
	}
      else
	{
	  sprintf (OutputName, "bosons_%s%s_%s_n_%d_2s_%d_sz_%d_fluxup_%.6f_fluxdown_%.6f.dat", GeometryName, DiscreteSymmetryName, 
		   Manager.GetString("interaction-name"), 
		   NbrBosons, LzMax, TotalSpin, Ratio, Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"));
	}
    }
  else
    {
      if ((Manager.GetDouble("spinup-flux") == 0.0) && (Manager.GetDouble("spindown-flux") == 0.0))
	{
	  sprintf (OutputName, "bosons_%s%s_%s_n_%d_2s_%d.dat", GeometryName, DiscreteSymmetryName, Manager.GetString("interaction-name"), NbrBosons, LzMax);
	}
      else
	{
	  sprintf (OutputName, "bosons_%s%s_%s_n_%d_2s_%d_fluxup_%.6f_fluxdown_%.6f.dat", GeometryName, DiscreteSymmetryName, 
		   Manager.GetString("interaction-name"), NbrBosons, LzMax,
		   Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"));
	}
    }

  int Max = (LzMax * (NbrUp+NbrDown));
  cout << "maximum Ky value = " << Max << endl;
  int LSign = 1;
  if (Manager.GetBoolean("negative-lz"))
    LSign = -1;
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
  if (Manager.GetBoolean("all-lz"))
    {
      L = -Max;
    }
  for (; L <= Max; L += 2)
    {
      double Shift = Manager.GetDouble("energy-shift");
      ParticleOnSphereWithSpin* Space = 0;
      Space = (ParticleOnSphereWithSpin*) ParticleManager.GetHilbertSpace(L * LSign);
      if (Space->GetHilbertSpaceDimension() > 0)
	{
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  if (Manager.GetBoolean("show-space"))
	    for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	      Space->PrintState(cout, i) << endl;
	  
	  AbstractQHEOnSphereWithSpinHamiltonian* Hamiltonian;
	  Hamiltonian = new ParticleOnCylinderWithSpinGenericHamiltonian(Space, NbrBosons, LzMax, Ratio,
									 NbrPseudoPotentials[0], PseudoPotentials[0],
									 NbrPseudoPotentials[1], PseudoPotentials[1],
									 NbrPseudoPotentials[2], PseudoPotentials[2],
									 Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"),
									 Architecture.GetArchitecture(), Memory, 0, OneBodyPseudoPotentials[0], 
									 OneBodyPseudoPotentials[1], OneBodyPseudoPotentials[2]);
	  
	  
	  Hamiltonian->ShiftHamiltonian(Shift);
	  if (SavePrecalculationFileName != 0)
	    {
	      Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
	    }
	  char* EigenvectorName = 0;
	  if (Manager.GetBoolean("eigenstate") == true)	
	    {
	      char* TmpName = RemoveExtensionFromFileName(OutputName, ".dat");
	      EigenvectorName = new char [32 + strlen(TmpName)];
	      sprintf (EigenvectorName, "%s_lz_%d", TmpName, (L * LSign));
	      delete[] TmpName;
	    }
	  QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L*LSign, Shift, OutputName, FirstRun, EigenvectorName, LzMax);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  delete Hamiltonian;
	  if (EigenvectorName != 0)
	    {
	      delete[] EigenvectorName;
	      EigenvectorName = 0;
	    }
	  if (FirstRun == true)
	    FirstRun = false;
	}
      delete Space;      
    }
  delete[] OutputName;
  return 0;
}


