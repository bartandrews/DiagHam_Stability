#include "HilbertSpace/FermionOnLattice.h"
#include "Hamiltonian/ParticleOnLatticeDeltaHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "MainTask/QHEOnLatticeMainTask.h"

#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
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

  OptionManager Manager ("FQHELatticeFermions" , "0.01");  
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  QHEOnLatticeMainTask::AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 8);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice", 5);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice", 1);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice (-1=all)", -1);
  (*SystemGroup) += new SingleDoubleOption  ('u', "contactU", "prefactor U of the contact interaction (kinetic term ~ 1)", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('d', "deltaPotential", "Introduce a delta-potential at the origin", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('R', "randomPotential", "Introduce a random potential at all sites", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "positive-hopping", "choose positive sign of hopping terms", false);
  (*SystemGroup) += new BooleanOption  ('\n', "all-flux", "calculate all values of the flux to test symmetry under n_phi->1-n_phi", false);

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrFermions = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  int NbrSites = Lx*Ly;  
  bool ReverseHopping = Manager.GetBoolean("positive-hopping");
  double ContactU = Manager.GetDouble("contactU");
  double Delta = Manager.GetDouble("deltaPotential");
  double Random = Manager.GetDouble("randomPotential");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  bool FirstRun = true;

  if (Manager.GetString("energy-expectation") != 0 ) Memory = 0x0l;

  int NbrFluxValues = 1;
  if (NbrFluxQuanta == -1)
    {
      NbrFluxQuanta = 0;
      if (Manager.GetBoolean("all-flux"))
	NbrFluxValues = NbrSites+1;
      else
	NbrFluxValues = (NbrSites+2)/2;
    }

  char* OutputName;
  char reverseHoppingString[4]="";
  char deltaString[20]="";
  char interactionStr[20]="";
  if ( (OutputName = Manager.GetString("output-file")) == NULL)
    {
      OutputName = new char [256];      
      if (ReverseHopping)
	sprintf(reverseHoppingString,"rh_");
      if (Delta!=0.0)
	sprintf(deltaString,"d_%g_",Delta);
      if (Random!=0.0)
	sprintf(deltaString,"R_%g_",Random);
      sprintf(interactionStr,"_u_%g", ContactU);
      if (NbrFluxValues == 1)
	sprintf (OutputName, "fermions_lattice_n_%d_x_%d_y_%d%s_%s%sq_%d.dat", NbrFermions, Lx, Ly, interactionStr, reverseHoppingString, deltaString, NbrFluxQuanta);
      else
	sprintf (OutputName, "fermions_lattice_n_%d_x_%d_y_%d%s_%s%sq.dat", NbrFermions, Lx, Ly, interactionStr, reverseHoppingString, deltaString);
    }
  ParticleOnLattice* Space = new FermionOnLattice(NbrFermions, Lx, Ly, NbrFluxQuanta, MemorySpace);

  
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  
  AbstractQHEOnLatticeHamiltonian* Hamiltonian;
  Hamiltonian = new ParticleOnLatticeDeltaHamiltonian(Space, NbrFermions, Lx, Ly, NbrFluxQuanta, ContactU,
						      ReverseHopping, Delta, Random, Architecture.GetArchitecture(), Memory, LoadPrecalculationFileName);

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

  for (int iter=0; iter<NbrFluxValues; ++iter, ++NbrFluxQuanta)
    {
      cout << "----------------------------------------------------------------" << endl;
      cout << "NbrFluxQuanta="<<NbrFluxQuanta<<endl;
      
      if (!FirstRun) Hamiltonian->SetNbrFluxQuanta(NbrFluxQuanta);
  
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate"))	
	{
	  EigenvectorName = new char [64];
	  sprintf (EigenvectorName, "fermions_lattice_n_%d_x_%d_y_%d%s_%s%sq_%d", NbrFermions, Lx, Ly, interactionStr, reverseHoppingString, deltaString, NbrFluxQuanta);
	}
      QHEOnLatticeMainTask Task (&Manager, Space, Hamiltonian, NbrFluxQuanta, 0.0, OutputName, FirstRun, EigenvectorName);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      if (FirstRun == true)
	FirstRun = false;
    }
  
  delete Hamiltonian;
  delete Space;  

  return 0;
}
