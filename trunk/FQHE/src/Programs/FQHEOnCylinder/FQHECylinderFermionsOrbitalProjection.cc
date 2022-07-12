#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereFull.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereLong.h"

#include "Hamiltonian/ParticleOnCylinderOrbitalProjection.h"

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



  OptionManager Manager ("FQHECylinderFermionsOrbitalProjection" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption ('y', "ky-momentum", "constraint on the total momentum along y-axis (negative if none)", -1);
  (*SystemGroup) += new SingleDoubleOption ('r', "ratio", "ratio between the height and length of the cylinder (LH=2pi r N_{orb})", 1.0);
  (*SystemGroup) += new SingleIntegerOption ('\n', "orbital-index", "index of the orbital to be projected out", 0);
  (*SystemGroup) += new SingleDoubleOption ('a', "anisotropy", "shape (anisotropy) parameter", 1.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "x0", "x0 coordinate for the center of the orbital", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "y0", "y0 coordinate for the center of the orbital", 0.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "orbitalprojector");
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

  cout<<" ----------------------------------"<<endl;
  cout<<" Projector on mth coherent state centered at (x0,y0) "<<endl;
  cout<<" via formula n_m(x0,y0) = int_d^2q Lm(q^2/2) exp(-q^2/4) exp(-iqx.x0-iqy.y0) exp(iq.R)"<<endl;
  cout<<" ----------------------------------"<<endl;
  

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int MaxMomentum = ((SingleIntegerOption*) Manager["max-momentum"])->GetInteger();
  int Momentum = ((SingleIntegerOption*) Manager["ky-momentum"])->GetInteger();
  int NbrKy = 1;

  double XRatio = ((SingleDoubleOption*) Manager["ratio"])->GetDouble();

  double X0 = ((SingleDoubleOption*) Manager["x0"])->GetDouble();
  double Y0 = ((SingleDoubleOption*) Manager["y0"])->GetDouble();
  int OrbitalIndex = Manager.GetInteger("orbital-index");
  double Anisotropy = ((SingleDoubleOption*) Manager["anisotropy"])->GetDouble();

  cout<<"Projector onto orbital "<<OrbitalIndex<<" centered at x0= "<<X0<<" y0= "<<Y0<<endl;
  cout<<"Orbital shape: anisotropy = "<<Anisotropy<<endl;

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  bool FirstRun = true;
  
  char* OutputNameLz = new char [256];
  sprintf (OutputNameLz, "fermions_cylinder_ky_%s_n_%d_2s_%d_ratio_%f.dat", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum, XRatio);
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  File.precision(14);

  int Max = ((MaxMomentum - NbrParticles + 1) * NbrParticles);

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
     Max = Ky + (2 * (NbrKy - 1));
   }

  for (; Ky <= Max; Ky+=2 )
    {     
      cout<<" Ky = " << Ky << endl;
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;

      ParticleOnSphere* Space = 0;

#ifdef __64_BITS__
	  if (MaxMomentum <= 62)
#else
	  if (MaxMomentum <= 30)
#endif
  	    Space = new FermionOnSphereFull(NbrParticles, MaxMomentum);

	  else
#ifdef __128_BIT_LONGLONG__
	    if (MaxMomentum <= 126)
#else
	      if (MaxMomentum <= 62)
#endif
 	        Space = new FermionOnSphereLong(NbrParticles, Ky, MaxMomentum);
	      else
		Space = new FermionOnSphereUnlimited(NbrParticles, Ky, MaxMomentum);


      cout << " Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;

      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();

      AbstractQHEHamiltonian* Hamiltonian = new ParticleOnCylinderOrbitalProjection (Space, NbrParticles, MaxMomentum, XRatio,  OrbitalIndex, Anisotropy, X0, Y0, Architecture.GetArchitecture(), Memory);

      double Shift = 0.0;
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
