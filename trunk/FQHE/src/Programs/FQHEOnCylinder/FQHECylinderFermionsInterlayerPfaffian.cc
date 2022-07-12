#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetry.h"

#include "Hamiltonian/ParticleOnCylinderInterlayerPfaffian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/FQHEOnTorusMainTask.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "QuantumNumber/AbstractQuantumNumber.h"

#include "GeneralTools/ListIterator.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
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



  OptionManager Manager ("FQHECylinderFermionsInterlayerPfaffian" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption ('l', "max-momentum", "maximum momentum for a single particle (total l+1 orbitals, from -l/2 to l/2)", 15);
  (*SystemGroup) += new SingleIntegerOption ('s', "total-sz", "projection of spin", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use Sz <-> -Sz symmetrized version of the basis (only valid if total-sz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-szparity", "select the  Sz <-> -Sz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*SystemGroup) += new SingleIntegerOption ('y', "ky-momentum", "constraint on the total momentum along y-axis (negative if none)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-ky", "number of Ky values to evaluate", -1);
  (*SystemGroup) += new SingleDoubleOption ('r', "ratio", "ratio between the height and length of the cylinder (LH=2pi r N_{orb})", 1.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "3b");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");

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
      cout << "see man page for option syntax or type FQHETorusFermionsLaplacianDelta -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int SzTotal =((SingleIntegerOption*) Manager["total-sz"])->GetInteger();
  bool LzSymmetrizedBasis = Manager.GetBoolean("lzsymmetrized-basis");
  bool SzSymmetrizedBasis = Manager.GetBoolean("szsymmetrized-basis");
  int MaxMomentum = ((SingleIntegerOption*) Manager["max-momentum"])->GetInteger();
  int Momentum = ((SingleIntegerOption*) Manager["ky-momentum"])->GetInteger();
  int NbrKy = Manager.GetInteger("nbr-ky");
  double XRatio = ((SingleDoubleOption*) Manager["ratio"])->GetDouble();

  unsigned long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  bool FirstRun = true;
  
  char* OutputNameLz = new char [256];
  sprintf (OutputNameLz, "fermions_cylinder_ky_%s_n_%d_2s_%d_ratio_%f.dat", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum, XRatio);
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  File.precision(14);

  int NbrUp = (NbrParticles + SzTotal) >> 1;
  int NbrDown = (NbrParticles - SzTotal) >> 1;
  if ((NbrUp < 0 ) || (NbrDown < 0 ))
    {
      cout << "This value of the spin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }

  int Max = ((MaxMomentum - NbrUp + 1) * NbrUp) + ((MaxMomentum - NbrDown + 1) * NbrDown);

  int  Ky = Momentum;
  if (Ky < -Max)
    Ky = -Max;
  else
    if (Ky > Max)
      Ky = Max;
  if ((abs(Max) & 1) != (abs(Momentum) & 1))
    Ky += 1;
  if (NbrKy > 0)
   {
      if (Ky + (2 * (NbrKy - 1)) < Max)
        Max = Ky + (2 * (NbrKy - 1));
   }

  for (; Ky <= Max; Ky+=2 )
    {     
      cout<<" Ky = " << Ky << endl;
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;

     ParticleOnSphereWithSpin* Space = 0;
#ifdef __64_BITS__
	      if (MaxMomentum <= 31)
#else
		if (MaxMomentum <= 15)
#endif
		  {
                    if ((SzSymmetrizedBasis == false) && (LzSymmetrizedBasis == false))
	  	       Space = new FermionOnSphereWithSpin(NbrParticles, Momentum, MaxMomentum, SzTotal);
                    else //either Lz or Sz symmetrized basis
                      {
	                if ((SzSymmetrizedBasis == true)  && (SzTotal == 0) && (LzSymmetrizedBasis == true) && (Momentum == 0))
		           {
		               Space = new FermionOnSphereWithSpinLzSzSymmetry(NbrParticles, MaxMomentum, Manager.GetBoolean("minus-szparity"),
								    Manager.GetBoolean("minus-lzparity"));
		           }
		        else 
                           if ((SzSymmetrizedBasis == true)  && (SzTotal == 0))
                              {
		                 Space = new FermionOnSphereWithSpinSzSymmetry(NbrParticles, Momentum, MaxMomentum, Manager.GetBoolean("minus-szparity"));
                              }
		           else
			         Space = new FermionOnSphereWithSpinLzSymmetry(NbrParticles, MaxMomentum, SzTotal, Manager.GetBoolean("minus-lzparity"));
                       }	
     	           }

      cout << " Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;

      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();

      AbstractQHEHamiltonian* Hamiltonian = new ParticleOnCylinderInterlayerPfaffian (Space, NbrParticles, MaxMomentum, XRatio, Architecture.GetArchitecture(), Memory);

      double Shift = -10.0;
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [256];
	  sprintf (EigenvectorName, "fermions_cylinder_%s_n_%d_2s_%d_ratio_%f_ky_%d", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum, XRatio, Ky);
	}
      FQHEOnTorusMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, Ky, Shift, OutputNameLz, FirstRun, EigenvectorName);
      Task.SetKxValue(-1);
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

  return 0;
}
