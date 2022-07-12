#include "HilbertSpace/BosonOnLatticeGeneric.h"
#include "HilbertSpace/HardCoreBosonOnLatticeGeneric.h"
#include "HilbertSpace/SingleParticleOnLatticeGeneric.h"
#include "Hamiltonian/ParticleOnLatticeGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeExternalHamiltonian.h"
#include "Operator/ParticleOnLatticeVacancyOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "MainTask/QHEOnLatticeMainTask.h"
#include "MainTask/GenericRealMainTask.h"

#include "Tools/FQHESpectrum/LatticePhases.h"
#include "Tools/FQHEWaveFunction/GutzwillerOnLatticeWaveFunction.h"

#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


#include "Matrix/ComplexMatrix.h"




int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHELatticeBosonsImpurity" , "0.01");  
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  LatticePhases::AddOptionGroup(&Manager);
  Architecture.AddOptionGroup(&Manager);
  QHEOnLatticeMainTask::AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 8);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice (-1=all)", -1);
  (*SystemGroup) += new SingleDoubleOption  ('Q', "cont-flux", "multiples of flux quanta piercing the lattice", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('u', "contactU", "prefactor U of the contact interaction (kinetic term ~ 1)", 1.0);
  (*SystemGroup) += new MultipleDoubleOption  ('s', "solenoid-flux", "twist in periodic boundary conditions phi_x[,phi_y])",',');
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons");
  (*SystemGroup) += new SingleStringOption  ('e', "external-two-body", "use definition of two-body interactions from a file");
  (*SystemGroup) += new SingleStringOption  ('E', "external-name", "descriptor of external interaction (if in use)","ext");
  (*SystemGroup) += new SingleStringOption  ('\n', "one-body", "use definition of one-body potential/hoppings from a file (use -E flag to indicate presence in filename)");
  (*SystemGroup) += new BooleanOption  ('\n', "hopping-only", "evaluate only energy of hopping terms, excluding local potentials", false);
  (*SystemGroup) += new BooleanOption  ('\n', "all-flux", "calculate all values of the flux to test symmetry under n_phi->1-n_phi", false);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hamiltonian-shift", "additional shift of Hamiltonian",0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "shift-dice", "shift to zero of half-filled dice-lattice");
  (*SystemGroup) += new SingleStringOption('\n', "reference-state", "reference state vector, in which to create impurities");
  (*SystemGroup) += new BooleanOption  ('\n', "particle", "create particle rather than hole");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "impurity-sublattice", "sublattice on which impurity can be created", 0);
    
  (*PrecalculationGroup) += new BooleanOption ('\n', "no-hermitian", "do not use hermitian symmetry of the hamiltonian");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrBosons = Manager.GetInteger("nbr-particles");

  int TargetNbrBosons = NbrBosons-1;
  bool ParticleImpurity = Manager.GetBoolean("particle");
  if (ParticleImpurity) TargetNbrBosons+=2;
  if (TargetNbrBosons<=0)
    {
      cout << "Particle size too small" << endl;
      exit(1);
    }
  int NbrFluxQuanta = Manager.GetInteger("flux");  
  double SolenoidX=0.0, SolenoidY=0.0;
  {
    int tmpI;
    double *Fluxes=Manager.GetDoubles("solenoid-flux", tmpI);
    if (tmpI>0) SolenoidX=Fluxes[0];
    if (tmpI>1) SolenoidY=Fluxes[1];
    if (tmpI>0) delete [] Fluxes;
  }
  bool ReverseHopping = false;
  bool HardCore = Manager.GetBoolean("hard-core");
  double ContactU = Manager.GetDouble("contactU");
  double Shift = Manager.GetDouble("hamiltonian-shift");
  int ImpSublattice = Manager.GetInteger("impurity-sublattice");
  if (Manager.GetBoolean("shift-dice"))
    Shift = -NbrBosons*sqrt(6.0);
  if (HardCore) ContactU=0.0;
  if (ULONG_MAX>>20 < (unsigned long)Manager.GetInteger("memory"))
    cout << "Warning: integer overflow in memory request - you might want to use 64 bit code."<<endl;
  unsigned long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");

  if (Manager.GetString("reference-state") == 0 )
    {
      cout<<"Reference state is needed. Please use option --reference-state"<<endl;
    }

  // get the lattice geometry
  LatticePhases *Lattice = new LatticePhases();

  int NbrSites = Lattice->GetNbrSites();
  
  int NbrFluxValues = 1;
  if (Lattice->HavePredefinedFlux())
    {
      NbrFluxQuanta = Lattice->GetPredefinedFlux();
    }
  else
    {
      if (Manager.GetDouble("cont-flux")!=0.0)
	{
	  NbrFluxValues = 1;
	  NbrFluxQuanta = 0;
	}
      else
	{
	  if (NbrFluxQuanta == -1)
	    {
	      NbrFluxQuanta = 0;
	      if (Manager.GetBoolean("all-flux"))
		NbrFluxValues = NbrSites+1;
	      else
		NbrFluxValues = (NbrSites+2)/2;
	    }
	}
    }

  char* OutputName;
  char* LatticeName = Lattice->GeometryString();
  if ( (OutputName = Manager.GetString("output-file")) == NULL)
    {
      char *NewExtension = new char[20];
      
      if (Manager.GetBoolean("particle"))
	sprintf(NewExtension,"imp-p_s_%d.dat",ImpSublattice);
      else
	sprintf(NewExtension,"imp-h_s_%d.dat",ImpSublattice);
      OutputName = ReplaceExtensionToFileName(Manager.GetString("reference-state"),"vec",NewExtension);
      if (OutputName==0)
	{
	  OutputName = new char[512];
	  sprintf(OutputName,"%s.%s",Manager.GetString("reference-state"),NewExtension);
	}
      delete[] NewExtension;
    }

  cout << "Lattice geometry: "<<Lattice->GetLatticeLength(0)<<"x"<<Lattice->GetLatticeLength(1)<<endl;

    
  ParticleOnLattice* Space, *TargetSpace;
  if (NbrBosons==1)
    {
      Space = new SingleParticleOnLatticeGeneric(Lattice, NbrFluxQuanta, SolenoidX, SolenoidY);
    }
  else
    {
      if (HardCore)
	Space = new HardCoreBosonOnLatticeGeneric(NbrBosons, Lattice, NbrFluxQuanta, MemorySpace, SolenoidX, SolenoidY);
      else Space = new BosonOnLatticeGeneric(NbrBosons, Lattice, NbrFluxQuanta, MemorySpace, SolenoidX, SolenoidY);
    }
  if (TargetNbrBosons==1)
      {
      TargetSpace = new SingleParticleOnLatticeGeneric(Lattice, NbrFluxQuanta, SolenoidX, SolenoidY);
    }
  else
    {
      if (HardCore)
	TargetSpace = new HardCoreBosonOnLatticeGeneric(TargetNbrBosons, Lattice, NbrFluxQuanta, MemorySpace, SolenoidX, SolenoidY);
      else TargetSpace = new BosonOnLatticeGeneric(TargetNbrBosons, Lattice, NbrFluxQuanta, MemorySpace, SolenoidX, SolenoidY);
    }

  Space->SetTargetSpace(TargetSpace);

  char* StateFileName = Manager.GetString("reference-state");
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
  if (State.GetVectorDimension()!=Space->GetHilbertSpaceDimension())
    {
      cout << "error: vector and Hilbert-space have unequal dimensions"<<endl;
      return -1;
    }
  
  Architecture.GetArchitecture()->SetDimension(TargetSpace->GetHilbertSpaceDimension());

  AbstractQHEOnLatticeHamiltonian* Hamiltonian;
  if (Manager.GetString("external-two-body")==NULL)
    Hamiltonian = new ParticleOnLatticeGenericHamiltonian(TargetSpace, TargetNbrBosons, Lattice, NbrFluxQuanta, ContactU,
		       				  ReverseHopping, Architecture.GetArchitecture(), Memory, LoadPrecalculationFileName, Manager.GetDouble("cont-flux"), !Manager.GetBoolean("no-hermitian"), Manager.GetBoolean("hopping-only"));
  else
    Hamiltonian = new ParticleOnLatticeExternalHamiltonian(TargetSpace, TargetNbrBosons, NbrSites, Manager.GetString("one-body"),
							   Manager.GetString("external-two-body"), Architecture.GetArchitecture(),
							   Memory, LoadPrecalculationFileName, !Manager.GetBoolean("no-hermitian"));

  BackUpFile(OutputName);
  ofstream File(OutputName, ios::out);

  File << "#kx ky E"<< endl;
  for (int kx=0; kx<Lattice->GetLatticeLength(0); ++kx)
    for (int ky=0; ky<Lattice->GetLatticeLength(1); ++ky)
      {
	ComplexVector ImpurityState(TargetSpace->GetHilbertSpaceDimension());
	ComplexVector TmpState(TargetSpace->GetHilbertSpaceDimension());

	ParticleOnLatticeVacancyOperator VacancyOperator (Space, TargetSpace, Lattice, kx, ky, Manager.GetInteger("impurity-sublattice"),ParticleImpurity);

	Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	VectorOperatorMultiplyOperation Operation1(&VacancyOperator, &State, &ImpurityState);
	Operation1.ApplyOperation(Architecture.GetArchitecture());
	ImpurityState/=ImpurityState.Norm();

	Architecture.GetArchitecture()->SetDimension(TargetSpace->GetHilbertSpaceDimension());
	VectorHamiltonianMultiplyOperation Operation2 (Hamiltonian, &ImpurityState, &TmpState);
	Operation2.ApplyOperation(Architecture.GetArchitecture());
	Complex EnergyValue = ImpurityState*TmpState;
	double DisplayKx = 2.0/Lattice->GetLatticeLength(0)*M_PI*kx;	
	double DisplayKy = 2.0/Lattice->GetLatticeLength(1)*M_PI*ky;
	if (DisplayKx>M_PI) DisplayKx-=2.0*M_PI;
	if (DisplayKy>M_PI) DisplayKy-=2.0*M_PI;
	cout << kx << " " << ky << " " << DisplayKx <<" "<< DisplayKy << " " << EnergyValue.Re << endl;
	File << kx << " " << ky << " " << DisplayKx <<" "<< DisplayKy << " " << EnergyValue.Re << endl;
      }
  File.close();
  
  delete Hamiltonian;
  delete Space;
  delete TargetSpace;
  delete Lattice;
  delete [] LatticeName;
  delete [] OutputName;
  return 0;
}
