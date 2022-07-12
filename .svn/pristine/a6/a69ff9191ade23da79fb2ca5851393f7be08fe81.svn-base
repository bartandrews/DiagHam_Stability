#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereLong.h"

#include "Hamiltonian/ParticleOnCylinderGaffnianHamiltonian.h"
#include "Hamiltonian/ParticleOnCylinderHaffnianHamiltonian.h"


#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/FQHEOnTorusMainTask.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "QuantumNumber/AbstractQuantumNumber.h"

#include "GeneralTools/ListIterator.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

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



  OptionManager Manager ("FQHECylinderFermionsGaffnian" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-ky", "number of Ky values to evaluate", -1);
  (*SystemGroup) += new SingleDoubleOption ('r', "ratio", "ratio between the height and length of the cylinder (LH=2pi r N_{orb})", 1.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "3b-coupling", "amplitude of the 3-body interaction", 1.0);
  (*SystemGroup) += new BooleanOption  ('\n', "haffnian", "diagonalize Haffnian Hamiltonian instead of Gaffnian");
  (*SystemGroup) += new SingleDoubleOption ('\n', "confinement-potential", "amplitude of the quadratic confinement potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "electric-field", "parameter for the value of the electric field applied along the cylinder (a=eEl_B^2/hbar omega_c", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "b-field", "parameter for the value of the magnetic field [in T] when also the electric field is present (needed to set the scale for the kinetic term)", 0.0);
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
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
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
  int MaxMomentum = ((SingleIntegerOption*) Manager["max-momentum"])->GetInteger();
  int Momentum = ((SingleIntegerOption*) Manager["ky-momentum"])->GetInteger();
  int NbrKy = Manager.GetInteger("nbr-ky");
  double XRatio = ((SingleDoubleOption*) Manager["ratio"])->GetDouble();
  double Confinement = ((SingleDoubleOption*) Manager["confinement-potential"])->GetDouble();
  if (Confinement != 0.0)
    {
      cout << "Assuming quadratic confining potential sum_m (a*X_m^2) c_m^+ c_m " << endl;
      cout << "where X_m=2pi m/L and a = " << Confinement << endl;
    }
  double ThreeBodyCoupling = ((SingleDoubleOption*) Manager["3b-coupling"])->GetDouble();

  double ElectricFieldParameter = ((SingleDoubleOption*) Manager["electric-field"])->GetDouble();
  double BFieldParameter = ((SingleDoubleOption*) Manager["b-field"])->GetDouble();
  if (ElectricFieldParameter != 0)
   {
      cout << "Electric field is applied along the cylinder with magnitude a=eEl_B^2/hbar omega_c = " << ElectricFieldParameter << endl;
      cout << "The energies are given in units e^2/epsilon l_B, but there is a kinetic term " << endl;
      cout << "K= (hbar omega_c/2/e^2/epsilon l_B) [sqrt(1+a) Ne + sum_m a/(1+a) (2pi m/L)^2 l_B^2 c_m^+ c_m] " << endl;
      cout << "The code sets (hbar omega_c/2/e^2/epsilon l_B) = 0.194 sqrt{B[T]} " << endl;  
      cout << "and neglects the term sqrt(1+a)N_e (overall constant)." << endl;
   }

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  bool FirstRun = true;
  
  char* OutputNameLz = new char [256];
  sprintf (OutputNameLz, "fermions_cylinder_ky_%s_n_%d_2s_%d_ratio_%f.dat", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum, XRatio);
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  File.precision(14);

  int Max = (MaxMomentum * NbrParticles);
  int  Ky = Momentum;

  if ((abs(Max) & 1) != (Ky & 1))
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
  	    Space = new FermionOnSphere(NbrParticles, Ky, MaxMomentum);

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

      AbstractQHEHamiltonian* Hamiltonian;
      if (((BooleanOption*) Manager["haffnian"])->GetBoolean() == true)
          Hamiltonian = new ParticleOnCylinderHaffnianHamiltonian (Space, NbrParticles, MaxMomentum, XRatio, ThreeBodyCoupling, Confinement, ElectricFieldParameter, BFieldParameter, Architecture.GetArchitecture(), Memory);
      else
          Hamiltonian = new ParticleOnCylinderGaffnianHamiltonian (Space, NbrParticles, MaxMomentum, XRatio, ThreeBodyCoupling, Confinement, ElectricFieldParameter, BFieldParameter, Architecture.GetArchitecture(), Memory);


      double Shift = 0.0;
      Hamiltonian->ShiftHamiltonian(Shift);

     if (Manager.GetString("energy-expectation") != 0 )
	{
	  char* StateFileName = Manager.GetString("energy-expectation");
	  if (IsFile(StateFileName) == false)
	    {
	      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
	      return -1;           
	    }
	  ComplexVector State;
	  if (State.ReadVector(StateFileName) == false)
	    {
	      cout << "error while reading " << StateFileName << endl;
	      return -1;
	    }
	  if (State.GetVectorDimension() != Space->GetHilbertSpaceDimension())
	    {
	      cout << "error: vector and Hilbert-space have unequal dimensions"<<endl;
	      return -1;
	    }
	  ComplexVector TmpState(Space->GetHilbertSpaceDimension());
	  VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  Complex EnergyValue = State * TmpState;
	  cout << "< Energy > = " << (EnergyValue.Re - Shift) << " " << EnergyValue.Im << endl;
	  return 0;
	}

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
