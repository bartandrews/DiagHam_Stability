#include "HilbertSpace/BosonOnLatticeKy.h"
#include "HilbertSpace/HardCoreBosonOnLatticeKy.h"
#include "HilbertSpace/BosonOnLatticeKyLong.h"
#include "Hamiltonian/ParticleOnLatticeWithKyNBodyDeltaHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"


#include "MainTask/QHEOnLatticeMainTask.h"

#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"

#include <iostream>
#include <cstdlib>
#include <climits>
#include <cmath>
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

  OptionManager Manager ("FQHELatticeBosonsKySymNBodyHardcore" , "0.01");  
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
  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 8);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice", 5);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice", 1);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice (-1=all)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('k', "ky", "constraint of momentum in y-direction (-1=all)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-body", "number of body involved in the hardcore", 2);
	
  (*SystemGroup) += new SingleDoubleOption  ('u', "contactU", "prefactor U of the three-body contact interaction (kinetic term ~ 1)", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "2body", "coefficient of an additional 2-body contact interaction", 0.0);
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons");
  (*SystemGroup) += new SingleDoubleOption  ('R', "randomPotential", "Introduce a random potential at all sites", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "positive-hopping", "choose positive sign of hopping terms", false);
  (*SystemGroup) += new BooleanOption  ('\n', "all-flux", "calculate all values of the flux to test symmetry under n_phi->1-n_phi", false);
  (*SystemGroup) += new BooleanOption  ('\n', "all-ky", "calculate all values of the flux to test symmetry under ky->ky_max-ky", false);

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new BooleanOption  ('\n', "get-hvalue", "show energy expectation value for eigenstates", false);
  (*MiscGroup) += new  BooleanOption ('\n',"show-basis", "show the basis of the Hilbert-space");
  (*MiscGroup) += new  BooleanOption ('\n',"show-hamiltonian", "show Hamiltonian matrix, and exit");
#ifdef HAVE_ARPACK
  (*MiscGroup) += new  BooleanOption ('\n',"use-arpack","use ARPACK routines for Lanczos algorithm");
#endif
  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  int NbrSites = Lx*Ly;
  int NbrBody = Manager.GetInteger("nbr-body");
  bool ReverseHopping = Manager.GetBoolean("positive-hopping");
  bool HardCore = Manager.GetBoolean("hard-core");
  double ContactU = Manager.GetDouble("contactU");
  if (HardCore) ContactU=0.0;
  double Random = Manager.GetDouble("randomPotential");
  if (ULONG_MAX>>20 < (unsigned long)Manager.GetInteger("memory"))
    cout << "Warning: integer overflow in memory request - you might want to use 64 bit code."<<endl;
  unsigned long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  double TwoBodyU = Manager.GetDouble("2body");
  
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
  char randomString[20]="";
  char kyString[20]="";
  char interactionStr[20]="";

  if ( (OutputName = Manager.GetString("output-file")) == NULL)
    {
      OutputName = new char [256];      
      if (ReverseHopping)
	sprintf(reverseHoppingString,"_rh");
      if (Random!=0.0)
	sprintf(randomString,"_R_%g",Random);
      if (HardCore)
	sprintf(interactionStr,"_hardcore");
      else sprintf(interactionStr,"_u_%g", ContactU);
      if (Manager.GetInteger("ky")>=0)
	sprintf(kyString,"_k_%ld",Manager.GetInteger("ky"));
      else
	sprintf(kyString,"_k");
      if (NbrFluxValues == 1)
	sprintf (OutputName, "bosons_lattice_n_%d_x_%d_y_%d%s%s%s%s_q_%d.dat", NbrBosons, Lx, Ly, interactionStr, reverseHoppingString, randomString, kyString, NbrFluxQuanta);
      else
	sprintf (OutputName, "bosons_lattice_n_%d_x_%d_y_%d%s%s%s%s_q.dat", NbrBosons, Lx, Ly, interactionStr, reverseHoppingString, randomString, kyString);
    }
  ParticleOnLattice* Space=0;
  
  AbstractQHEOnLatticeHamiltonian* Hamiltonian=0;
  
  if (Manager.GetString("energy-expectation") != 0 )
    {
      if (HardCore)
	{
	  cout << "Hard-Core bosons not defined with translations!"<<endl;
	  exit(1);
	  //Space =new HardCoreBosonOnLatticeKy(NbrBosons, Lx, Ly, Manager.GetInteger("ky"), NbrFluxQuanta, MemorySpace);
	}
      else 
	{
	  if (NbrBosons + NbrSites - 1 <63)
	    {
	      Space = new BosonOnLatticeKy(NbrBosons, Lx, Ly, Manager.GetInteger("ky"), NbrFluxQuanta, MemorySpace);
	    }
	  else
	    {
	      if (NbrBosons + NbrSites - 1 <127)
		{
		  Space = new BosonOnLatticeKyLong(NbrBosons, Lx, Ly, Manager.GetInteger("ky"), NbrFluxQuanta, MemorySpace);
		}
	      else
		{
		  cout <<"The number of coding bits needded excesses 128 which is not supported"<<endl;
		  exit(1);
		}
	    }
	}
      
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      
      Hamiltonian = new ParticleOnLatticeWithKyNBodyDeltaHamiltonian(Space, NbrBosons, Lx, Ly, ((BosonOnLatticeKy*)Space)->GetMaximumKy(),  NbrFluxQuanta,NbrBody, TwoBodyU, ContactU, ReverseHopping, Random, Architecture.GetArchitecture(), Memory, LoadPrecalculationFileName);
      
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
  
  //   ComplexMatrix HRe(Hamiltonian->GetHilbertSpaceDimension(),Hamiltonian->GetHilbertSpaceDimension());
  //   ComplexMatrix HIm(Hamiltonian->GetHilbertSpaceDimension(),Hamiltonian->GetHilbertSpaceDimension());
//   GetHamiltonian(Hamiltonian,HRe);
//   GetHamiltonianIm(Hamiltonian,HIm);
//   Complex one, two, M_I(0.0,1.0);
//   for (int i=0; i<Hamiltonian->GetHilbertSpaceDimension(); ++i)
//     for (int j=0; j<Hamiltonian->GetHilbertSpaceDimension(); ++j)
//       {
// 	HRe.GetMatrixElement(i,j,one);
// 	HIm.GetMatrixElement(i,j,two);
// 	one= one*M_I;
// 	if (Norm(one-two)>1e-10)
// 	  cout << "Discrepancy in "<<i<<", "<<j<<": "<<one << " vs " << two << endl;
//       }
//   for (int i=0; i<Hamiltonian->GetHilbertSpaceDimension(); ++i)
//     for (int j=0; j<i; ++j)
//       {
// 	HRe.GetMatrixElement(i,j,one);
// 	HRe.GetMatrixElement(j,i,two);
// 	if (Norm(one-Conj(two))>1e-10)
// 	  cout << "Matrix not hermitian in "<<i<<", "<<j<<": "<<one << " vs " << two << endl;
//       }

//   ComplexVector TmpV1a (Hamiltonian->GetHilbertSpaceDimension(), true);
//   ComplexVector TmpV1b (Hamiltonian->GetHilbertSpaceDimension(), true);
//   ComplexVector TmpV2a (Hamiltonian->GetHilbertSpaceDimension(), true);
//   ComplexVector TmpV2b (Hamiltonian->GetHilbertSpaceDimension(), true);
//   for (int i = 0; i < Hamiltonian->GetHilbertSpaceDimension(); i++)
//     {
//       TmpV1a.Re(i) = (rand() - 32767) * 0.5;
//       TmpV1a.Im(i) = (rand() - 32767) * 0.5;
//     }
//   TmpV1a /= TmpV1a.Norm();
//   TmpV1b = TmpV1a*M_I;
//   Hamiltonian->LowLevelMultiply(TmpV1a, TmpV2a);
//   Hamiltonian->LowLevelMultiply(TmpV1b, TmpV2b);
//   for (int j=0; j<Hamiltonian->GetHilbertSpaceDimension(); ++j)
//       {
// 	one = TmpV2a[j];
// 	two = TmpV2b[j];
// 	one = one*M_I;
// 	if (Norm(one-two)>1e-10)
// 	  cout << "Discrepancy in "<<j<<": "<<one << " vs " << two << endl;
//       }  
  bool FirstRun=true;
  
  for (int iter=0; iter<NbrFluxValues; ++iter, ++NbrFluxQuanta)
    {
      cout << "================================================================" << endl;
      cout << "NbrFluxQuanta="<<NbrFluxQuanta<<endl;
      
      int Ky = Manager.GetInteger("ky");
      int MaxK = (Ky<0?1:Ky+1);
      int UpperLimit = MaxK;      
      for (int k=(Ky>=0?Ky:0); k<UpperLimit; ++k)
	{
	  if (Space!=0) 
	    delete Space;
	  if (Hamiltonian!=0) 
	    delete Hamiltonian;
	  
	  if (HardCore)
	    {
	      cout << "Hard-Core bosons not defined with translations!"<<endl;
	      exit(1);
	      //Space =new HardCoreBosonOnLatticeKy(NbrBosons, Lx, Ly, k, NbrFluxQuanta, MemorySpace);
	    }
	  else 
	    {
	      if (NbrBosons + NbrSites - 1 <63)
		{
		  Space = new BosonOnLatticeKy(NbrBosons, Lx, Ly, k, NbrFluxQuanta, MemorySpace);
		}
	      else
		{
		  if (NbrBosons + NbrSites - 1 <127)
		    {
		      Space = new BosonOnLatticeKyLong(NbrBosons, Lx, Ly, k, NbrFluxQuanta, MemorySpace);
		    }
		  else
		    {
		      cout <<"The number of coding bits needded excesses 128 which is not supported"<<endl;
		      exit(1);
		    }
		}
	    }
	  
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	  
	  if (Ky<0)
	    {
	      MaxK = ((BosonOnLatticeKy*)Space)->GetMaximumKy();
	      cout << "Maximum Ky="<<MaxK-1<<endl;
	      if (Manager.GetBoolean("all-ky"))
		{
		  UpperLimit = MaxK;
		}
	      else
		{
		  UpperLimit = (MaxK+2)/2;
		  cout << "Evaluating up to Ky="<<UpperLimit-1<<endl;
		}
	    }
	  
	  cout << "Momentum Ky="<<k<<": dim="<<Space->GetHilbertSpaceDimension()<<endl;
	  
	  Hamiltonian = new ParticleOnLatticeWithKyNBodyDeltaHamiltonian(Space, NbrBosons, Lx, Ly, ((BosonOnLatticeKy*)Space)->GetMaximumKy(),  NbrFluxQuanta,NbrBody, 0.0, ContactU, ReverseHopping, Random, Architecture.GetArchitecture(), Memory, LoadPrecalculationFileName);  
          
	  char* EigenvectorName = 0;
	  if (Manager.GetBoolean("eigenstate"))	
	    {
	      EigenvectorName = new char [100];	      
	      sprintf (EigenvectorName, "bosons_lattice_n_%d_x_%d_y_%d%s%s%s_k_%d_q_%d", NbrBosons, Lx, Ly, interactionStr, reverseHoppingString, randomString, k, NbrFluxQuanta);
	    }
	  QHEOnLatticeMainTask Task (&Manager, Space,&Lanczos, Hamiltonian, NbrFluxQuanta, 0.0, OutputName, FirstRun, EigenvectorName, k);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  FirstRun=false;
	  if (EigenvectorName != 0)
	    {
	      delete[] EigenvectorName;
	    }
	  cout << "----------------------------------------------------------------" << endl;
	}
    }
  cout << "================================================================" << endl;
  delete Hamiltonian;
  delete Space;  
  
  return 0;
}
