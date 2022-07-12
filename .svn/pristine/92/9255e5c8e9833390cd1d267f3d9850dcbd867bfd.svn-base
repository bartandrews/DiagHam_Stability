#include "HilbertSpace/BosonOnLatticeGeneric.h"
//#include "HilbertSpace/FermionOnLatticeGeneric.h"
#include "HilbertSpace/HardCoreBosonOnLatticeGeneric.h"
#include "HilbertSpace/SingleParticleOnLatticeGeneric.h"
#include "Hamiltonian/ParticleOnLatticeGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeExternalHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "MainTask/QHEOnLatticeMainTask.h"
#include "MainTask/GenericRealMainTask.h"

#include "Tools/FQHESpectrum/LatticePhases.h"
#include "Tools/FQHEWaveFunction/GutzwillerOnLatticeWaveFunction.h"

#include "LanczosAlgorithm/LanczosManager.h"

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


// store imaginary Hamiltonian into a complex matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix

ComplexMatrix& GetHamiltonianIm (AbstractHamiltonian *H, ComplexMatrix& M)
{
  ComplexVector TmpV1 (H->GetHilbertSpaceDimension(), true);
  ComplexVector TmpV2 (H->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < H->GetHilbertSpaceDimension(); i++)
    {
      TmpV1.Im(i) = 1.0;
      if (H->IsHermitian())
	H->HermitianLowLevelMultiply(TmpV1, TmpV2);
      else
	H->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = 0; j < H->GetHilbertSpaceDimension(); j++)
	{
	  M.SetMatrixElement(i, j, TmpV2[j]);
	}
      TmpV1.Im(i) = 0.0;
    }
  return M;  
}

// store imaginary Hamiltonian into an complex matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix

ComplexMatrix& GetHamiltonian (AbstractHamiltonian *H, ComplexMatrix& M)
{
  ComplexVector TmpV1 (H->GetHilbertSpaceDimension(), true);
  ComplexVector TmpV2 (H->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < H->GetHilbertSpaceDimension(); i++)
    {
      TmpV1.Re(i) = 1.0;
      if (H->IsHermitian())
	H->HermitianLowLevelMultiply(TmpV1, TmpV2);
      else
	H->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = 0; j < H->GetHilbertSpaceDimension(); j++)
	{
	  M.SetMatrixElement(i, j, TmpV2[j]);
	}
      TmpV1.Re(i) = 0.0;
    }
  return M;  
}


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHELatticeBosonsGeneric" , "0.01");  
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  LatticePhases::AddOptionGroup(&Manager);
  Architecture.AddOptionGroup(&Manager);
  
  LanczosManager Lanczos(true);  // functions in parallel to complex lattice main task...
  Lanczos.AddOptionGroup(&Manager);
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
  (*SystemGroup) += new BooleanOption  ('\n', "positive-hopping", "choose positive sign of hopping terms", false);
  (*SystemGroup) += new BooleanOption  ('\n', "hopping-only", "evaluate only energy of hopping terms, excluding local potentials", false);
  (*SystemGroup) += new BooleanOption  ('\n', "all-flux", "calculate all values of the flux to test symmetry under n_phi->1-n_phi", false);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hamiltonian-shift", "additional shift of Hamiltonian",0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "shift-dice", "shift to zero of half-filled dice-lattice");

  (*PrecalculationGroup) += new BooleanOption ('\n', "no-hermitian", "do not use hermitian symmetry of the hamiltonian");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption('\n', "optimize-condensate", "optimize a trial condensate wavefunction instead of diagonalizing");
  (*MiscGroup) += new SingleStringOption('\n', "init-parameters", "file with initial parameters");
  (*MiscGroup) += new SingleIntegerOption('\n', "trial-symmetry", "symmetry restrictions to trial state (1=odd/even)",0);
  (*MiscGroup) += new SingleDoubleOption('\n', "tolerance", "tolerance for variational parameters in condensate",1e-6);
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
  
  (*MiscGroup) += new BooleanOption  ('\n', "test-hamiltonian", "test hermiticity of Hamiltonian");
  (*MiscGroup) += new BooleanOption  ('\n', "force-complex", "always use complex vectors, even if Hamiltonian is real");
  (*MiscGroup) += new BooleanOption  ('\n', "get-hvalue", "show energy expectation value for eigenstates", false);
  (*MiscGroup) += new  BooleanOption ('\n',"show-basis", "show the basis of the Hilbert-space");
  (*MiscGroup) += new  BooleanOption ('\n',"show-hamiltonian", "show Hamiltonian matrix, and exit");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int NbrFluxQuanta = Manager.GetInteger("flux");  
  double SolenoidX=0.0, SolenoidY=0.0;
  {
    int tmpI;
    double *Fluxes=Manager.GetDoubles("solenoid-flux", tmpI);
    if (tmpI>0) SolenoidX=Fluxes[0];
    if (tmpI>1) SolenoidY=Fluxes[1];
    if (tmpI>0) delete [] Fluxes;
  }
  bool ReverseHopping = Manager.GetBoolean("positive-hopping");
  bool HardCore = Manager.GetBoolean("hard-core");
  double ContactU = Manager.GetDouble("contactU");
  double Shift = Manager.GetDouble("hamiltonian-shift");
  if (Manager.GetBoolean("shift-dice"))
    Shift = -NbrBosons*sqrt(6.0);
  if (HardCore) ContactU=0.0;
  if (ULONG_MAX>>20 < (unsigned long)Manager.GetInteger("memory"))
    cout << "Warning: integer overflow in memory request - you might want to use 64 bit code."<<endl;
  unsigned long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  bool FirstRun = true;

  if (Manager.GetString("energy-expectation") != 0 ) Memory = 0x0ul;

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
  char reverseHoppingString[4]="";
  char interactionStr[100]="";
  if ( (OutputName = Manager.GetString("output-file")) == NULL)
    {
      OutputName = new char [1024];
      if (Manager.GetString("external-two-body")==NULL)
	{
	  if (ReverseHopping)
	    sprintf(reverseHoppingString,"_rh");
	  if (HardCore)
	    {
	      sprintf(interactionStr,"_hardcore");
	    }
	  else
	    {
	      sprintf(interactionStr,"_u_%g",ContactU);
	    }
	  if ((SolenoidX!=0.0)||(SolenoidY!=0.0))
	    {
	      sprintf(interactionStr,"%s_s_%g_%g",interactionStr,SolenoidX,SolenoidY);
	    }
	  if (NbrFluxValues == 1)
	    {
	      if (Manager.GetDouble("cont-flux")!=0.0)
		sprintf (OutputName, "bosons_lattice_%s_n_%d%s%s_Q_%gs.dat", LatticeName, NbrBosons, interactionStr, reverseHoppingString, Manager.GetDouble("cont-flux"));
	      else
		sprintf (OutputName, "bosons_lattice_%s_n_%d%s%s_q_%d.dat", LatticeName, NbrBosons, interactionStr, reverseHoppingString, NbrFluxQuanta);
	    }
	  else
	    sprintf (OutputName, "bosons_lattice_%s_n_%d%s%s_q.dat", LatticeName, NbrBosons, interactionStr, reverseHoppingString);
	}
      else
	{
	  NbrFluxValues = 1;
	  NbrFluxQuanta = 3*NbrSites/2;
	  if (HardCore)
	    {
	      sprintf(interactionStr,"_hardcore");
	    }
	  if ((SolenoidX!=0.0)||(SolenoidY!=0.0))
	    {
	      sprintf(interactionStr,"%s_s_%g_%g",interactionStr,SolenoidX,SolenoidY);
	    }
	  sprintf (OutputName, "bosons_lattice_%s_n_%d_%s%s_q_%d.dat", LatticeName, NbrBosons, Manager.GetString("external-name"),
		   interactionStr, NbrFluxQuanta);
	}
    }
  ParticleOnLattice* Space;
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
  
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  
  AbstractQHEOnLatticeHamiltonian* Hamiltonian;
  if (Manager.GetString("external-two-body")==NULL)
    Hamiltonian = new ParticleOnLatticeGenericHamiltonian(Space, NbrBosons, Lattice, NbrFluxQuanta, ContactU,
							  ReverseHopping, Architecture.GetArchitecture(), Memory, LoadPrecalculationFileName, Manager.GetDouble("cont-flux"), !Manager.GetBoolean("no-hermitian"), Manager.GetBoolean("hopping-only"));
  else
    Hamiltonian = new ParticleOnLatticeExternalHamiltonian(Space, NbrBosons, NbrSites, Manager.GetString("one-body"),
							   Manager.GetString("external-two-body"), Architecture.GetArchitecture(),
							   Memory, LoadPrecalculationFileName, !Manager.GetBoolean("no-hermitian"));

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
	  if (State.GetVectorDimension()!=Space->GetHilbertSpaceDimension())
	    {
	      cout << "error: vector and Hilbert-space have unequal dimensions"<<endl;
	      return -1;
	    }
	  ComplexVector TmpState(Space->GetHilbertSpaceDimension());
	  VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  Complex EnergyValue = State*TmpState;
	  cout << "< Energy > = "<<EnergyValue<<endl;
	  return 0;
	}

//   // testing Hamiltonian:
  if (Manager.GetBoolean("test-hamiltonian"))    
    {
      if (Hamiltonian->GetHilbertSpaceDimension()>5000)
	{
	  cout << "Attention, debug mode of FQHELatticeBosonsGeneric run for large Hilbert-space"<<endl;
	}
      else
	{
	  ComplexMatrix HRe(Hamiltonian->GetHilbertSpaceDimension(),Hamiltonian->GetHilbertSpaceDimension());
	  ComplexMatrix HIm(Hamiltonian->GetHilbertSpaceDimension(),Hamiltonian->GetHilbertSpaceDimension());
	  GetHamiltonian(Hamiltonian,HRe);
	  GetHamiltonianIm(Hamiltonian,HIm);
	  Complex one, two, M_I(0.0,1.0);
	  for (int i=0; i<Hamiltonian->GetHilbertSpaceDimension(); ++i)
	    for (int j=0; j<Hamiltonian->GetHilbertSpaceDimension(); ++j)
	      {
		HRe.GetMatrixElement(i,j,one);
		HIm.GetMatrixElement(i,j,two);
		one*=M_I;
		if (Norm(one-two)>1e-10)
		  cout << "Discrepancy in "<<i<<", "<<j<<": "<<one << " vs " << two << endl;
	      }
	  cout << "HRe="<<endl<<HRe;
	  for (int i=0; i<Hamiltonian->GetHilbertSpaceDimension(); ++i)
	    {
	      HRe.GetMatrixElement(i,i,one);
	      if (Norm(one.Im)>1e-10)
		cout << "Matrix not hermitian in "<<i<<", "<<i<<": "<<one << endl;
	      for (int j=0; j<i; ++j)
		{
		  HRe.GetMatrixElement(i,j,one);
		  HRe.GetMatrixElement(j,i,two);
		  if (Norm(one-Conj(two))>1e-10)
		    cout << "Matrix not hermitian in "<<i<<", "<<j<<": "<<one << " vs " << two << endl;
		}
	    }
	  
	  ComplexVector TmpV1a (Hamiltonian->GetHilbertSpaceDimension(), true);
	  ComplexVector TmpV1b (Hamiltonian->GetHilbertSpaceDimension(), true);
	  ComplexVector TmpV2a (Hamiltonian->GetHilbertSpaceDimension(), true);
	  ComplexVector TmpV2b (Hamiltonian->GetHilbertSpaceDimension(), true);
	  for (int i = 0; i < Hamiltonian->GetHilbertSpaceDimension(); i++)
	    {
	      TmpV1a.Re(i) = (rand() - 32767) * 0.5;
	      TmpV1a.Im(i) = (rand() - 32767) * 0.5;
	    }
	  TmpV1a /= TmpV1a.Norm();
	  TmpV1b = TmpV1a*M_I;
	  Hamiltonian->LowLevelMultiply(TmpV1a, TmpV2a);
	  Hamiltonian->LowLevelMultiply(TmpV1b, TmpV2b);
	  for (int j=0; j<Hamiltonian->GetHilbertSpaceDimension(); ++j)
	    {
	      one = TmpV2a[j];
	      two = TmpV2b[j];
	      one = one*M_I;
	      if (Norm(one-two)>1e-10)
		cout << "Discrepancy in "<<j<<": "<<one << " vs " << two << endl;
	    }
	}
    }
  char *SubspaceLegend=new char[20];
  strcpy(SubspaceLegend,"NbrFlux");
  char* SubspaceStr=new char[5];
  for (int iter=0; iter<NbrFluxValues; ++iter, ++NbrFluxQuanta)
    {
      cout << "----------------------------------------------------------------" << endl;
      cout << "NbrFluxQuanta="<<NbrFluxQuanta<<endl;
      
      if (!FirstRun) Hamiltonian->SetNbrFluxQuanta(NbrFluxQuanta);
      
      sprintf(SubspaceStr,"%d",NbrFluxQuanta);
      char* EigenvectorName = 0;
      if ((Manager.GetBoolean("eigenstate")||(Manager.GetBoolean("optimize-condensate"))))
	{
	  EigenvectorName = new char [1024];
	  if ((NbrFluxValues == 1)&&(Manager.GetDouble("cont-flux")!=0.0))
	    sprintf (EigenvectorName, "bosons_lattice_%s_n_%d%s%s_Q_%gs", LatticeName, NbrBosons, interactionStr, reverseHoppingString, Manager.GetDouble("cont-flux"));
	  else
	    {
	      if (Manager.GetString("external-two-body")==NULL)
		sprintf (EigenvectorName, "bosons_lattice_%s_n_%d%s%s_q_%d", LatticeName, NbrBosons, interactionStr, reverseHoppingString, NbrFluxQuanta);
	      else
		sprintf (EigenvectorName, "bosons_lattice_%s_n_%d_%s%s_q_%d", LatticeName, NbrBosons, Manager.GetString("external-name"),
			 interactionStr, NbrFluxQuanta);
	    }
	}
      if (Manager.GetBoolean("optimize-condensate"))
	{	  
	  int Counter=0;
	  char* ParameterName = GetUniqueFileName(EigenvectorName,Counter,".cond.par");
	  char* WaveFunctionName = GetUniqueFileName(EigenvectorName,Counter,".cond.vec");
	  RealVector *InitialParameters = NULL;
	  if ((Manager.GetString("init-parameters")!=NULL))
	    {
	      InitialParameters = new RealVector;
	      if (InitialParameters->ReadVector(Manager.GetString("parameters"))==false)
		{
		  cout << "Could not read vector of initial parameters" <<Manager.GetString("parameters")<<endl;
		  exit(1);
		}
	    }
	  GutzwillerOnLatticeWaveFunction Condensate(NbrBosons, HardCore, Space, InitialParameters, Manager.GetInteger("trial-symmetry"));
	  Condensate.SetHamiltonian(Hamiltonian);
	  Condensate.SetArchitecture(Architecture.GetArchitecture());
	  Condensate.SetToRandomPhase();
	  int MaxEval = NbrSites*(NbrBosons+1)*2*Manager.GetInteger("nbr-iter");
	  double Energy=Condensate.Optimize(Manager.GetDouble("tolerance"), MaxEval);
	  Condensate.GetLastWaveFunction().WriteVector(WaveFunctionName);
	  Condensate.GetVariationalParameters().WriteVector(ParameterName);
	  cout << "Found condensate state with energy: "<<Energy<<endl<<WaveFunctionName<<endl;
	  ofstream File;
	  ifstream TestFile;
	  char *LogFileName = ReplaceExtensionToFileName(OutputName,"dat","cond.log");
	  TestFile.open(LogFileName, ios::in);
	  if (TestFile.is_open())
	    {
	      TestFile.close();
	      File.open(LogFileName, ios::app );
	    }
	  else
	    {
	      File.open(LogFileName, ios::out );
	      File << "#Index\tE_tot\tParameterFile"<<endl;
	    }
	  File << Counter << "\t" << Energy << "\t"<<ParameterName<<endl;
	  File.close();
	  if (InitialParameters!=NULL)
	    delete InitialParameters;
	  delete [] ParameterName;
	  delete [] LogFileName;
	  delete [] WaveFunctionName;
	}
      else
	{
	  if (Hamiltonian->IsComplex()||Manager.GetBoolean("force-complex"))
	    {
	      QHEOnLatticeMainTask Task (&Manager, Space,&Lanczos ,Hamiltonian, NbrFluxQuanta, Shift, OutputName, FirstRun, EigenvectorName);
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	    }
	  else
	    {
// 	      cout << "Attention: could simplify calculation, as this Hamiltonian is purely real!"<<endl;
// 	      QHEOnLatticeMainTask Task (&Manager, Space, Hamiltonian, NbrFluxQuanta, Shift, OutputName, FirstRun, EigenvectorName);
// 	      MainTaskOperation TaskOperation (&Task);
// 	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      GenericRealMainTask Task(&Manager, (AbstractHilbertSpace*)Space, &Lanczos, (AbstractHamiltonian*)Hamiltonian, SubspaceStr, SubspaceLegend,
				  Shift, OutputName, FirstRun, EigenvectorName, /* FakeComplex */ true);
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	    }
	}
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      if (FirstRun == true)
	FirstRun = false;
    }
  
  delete Hamiltonian;
  delete Space;
  delete Lattice;
  delete [] LatticeName;
  delete [] OutputName;
  delete [] SubspaceLegend;
  delete [] SubspaceStr;
  return 0;
}
