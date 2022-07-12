#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/HardCoreBosonOnLattice.h"
#include "HilbertSpace/SingleParticleOnLattice.h"
#include "Hamiltonian/ParticleOnLatticeDeltaHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeKapitMuellerHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeKapitMuellerMultiLayerHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/QHEOnLatticeMainTask.h"

#include "Tools/FQHEWaveFunction/GutzwillerOnLatticeWaveFunction.h"

#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
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

  OptionManager Manager ("FQHELatticeBosons" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-body", "number of body involed in the contact interaction", 2);
  (*SystemGroup) += new SingleDoubleOption  ('u', "contactU", "prefactor U of the contact interaction (kinetic term ~ 1)", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('w', "contactW", "prefactor W of an inter-layer contact interaction", 0.0);
  (*SystemGroup) += new MultipleDoubleOption  ('s', "solenoid-flux", "twist in periodic boundary conditions phi_x[,phi_y])",',');
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons");
  (*SystemGroup) += new SingleStringOption('l',"landau-axis","potential in the Landau gauge along axis","y");
  
  (*SystemGroup) += new SingleDoubleOption  ('D', "deltaPotential", "Introduce a delta-potential at the origin", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('R', "randomPotential", "Introduce a random potential at all sites", 0.0);
  (*SystemGroup) += new BooleanOption ('\n', "cylinder", "Consider a cylinder geometry (only for delta Hamiltonian)");
  (*SystemGroup) += new BooleanOption  ('\n', "positive-hopping", "choose positive sign of hopping terms", false);
  (*SystemGroup) += new BooleanOption  ('\n', "all-flux", "calculate all values of the flux to test symmetry under n_phi->1-n_phi", false);
  (*SystemGroup) += new SingleDoubleOption  ('K', "Kapit-Mueller", "Use the Kapit-Mueller hoppings with the argument being the maximum range", -1.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-layers", "number of layers for the Kapit-Mueller case", 1);
  (*SystemGroup) += new MultipleDoubleOption  ('\n', "branch-cuts", "branch cuts with 4 arguments (xi1,yi1,xj1,yj1) for each branch cut", ',', '_');
  (*SystemGroup) += new MultipleIntegerOption  ('\n', "branch-shift", "shifts for the individual branch-cuts s1,s2,... ", ',', '_');
  
  
  (*PrecalculationGroup) += new BooleanOption ('\n', "no-hermitian", "do not use hermitian symmetry of the hamiltonian");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);	
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption('\n', "optimize-condensate", "optimize a trial condensate wavefunction instead of diagonalizing");
  (*MiscGroup) += new SingleDoubleOption('\n', "tolerance", "tolerance for variational parameters in condensate",1e-6);
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
(*MiscGroup) += new BooleanOption  ('\n', "get-hvalue", "show energy expectation value for eigenstates", false);
  (*MiscGroup) += new  BooleanOption ('\n',"show-basis", "show the basis of the Hilbert-space");
  (*MiscGroup) += new  BooleanOption ('\n',"show-hamiltonian", "show Hamiltonian matrix, and exit");
#ifdef HAVE_ARPACK
  (*MiscGroup) += new  BooleanOption ('\n',"use-arpack","use ARPACK routines for Lanczos algorithm");
#endif
  (*MiscGroup) += new  BooleanOption ('\n',"no-single-particle-basis", "do not use the single-particle basis when -p 1 (for testing)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  int NbrSites = Lx*Ly;
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
  char LandauQuantization = Manager.GetString("landau-axis")[0];
  double ContactU = Manager.GetDouble("contactU");
  double ContactW = Manager.GetDouble("contactW"); // inter-layer interaction for Kapit-Mueller
  int NbrBody = Manager.GetInteger("nbr-body");
  if (HardCore) ContactU=0.0;
  double Delta = Manager.GetDouble("deltaPotential");
  double Random = Manager.GetDouble("randomPotential");
  if (ULONG_MAX>>20 < (unsigned long)Manager.GetInteger("memory"))
    cout << "Warning: integer overflow in memory request - you might want to use 64 bit code."<<endl;
  unsigned long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  bool FirstRun = true;

  int length;
  double *branchCuts = Manager.GetDoubles("branch-cuts",length);
  if (branchCuts!=NULL && length % 4 !=0)
    {
      std::cerr<<"error: need to have a multiple of four coordinates for --branch-cuts"<<endl;
      exit(1);
    }
  int NbrCuts = length / 4;
  int *branchShift = Manager.GetIntegers("branch-shift",length);
  if (branchShift!=NULL && length != NbrCuts)
    {
      std::cerr<<"error: number of arguments for --branch-shift needs to match number of branch cuts"<<endl;
      exit(1);
    }

  int NbrLayers = Manager.GetInteger("nbr-layers");

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
  char auxArguments[128]="";
  char deltaString[20]="";
  char interactionStr[100]="";
  int offset=0;
  if ( (OutputName = Manager.GetString("output-file")) == NULL)
    {
      OutputName = new char [256];      
      if (Manager.GetDouble("Kapit-Mueller")>0.0)
	{
	  offset+=sprintf(auxArguments,"KM_%g_", Manager.GetDouble("Kapit-Mueller"));
	  
	  if (NbrLayers>1)
	    offset+=sprintf(auxArguments+offset,"l_%d_", NbrLayers);

	  if (NbrCuts>0)
	    {
	      MultipleDoubleOption* option = (MultipleDoubleOption*)Manager["branch-cuts"];
	      char * coords = option->GetAsAString();
	      offset+=sprintf(auxArguments+offset,"cuts_%s_", coords);
	      delete [] coords;
	      if (branchShift!=0)
		{
		  bool nonTrivial=false;
		  for (int i=0; i<NbrCuts; ++i)
		    if (branchShift[i]!=1) nonTrivial=true;
		  
		  if (nonTrivial)
		    {
		      MultipleIntegerOption* option = (MultipleIntegerOption*)Manager["branch-cuts"];
		      
		      char * shift = option->GetAsAString();
		      offset+=sprintf(auxArguments+offset,"sh_%s_", shift);
		      delete [] shift;
		    }
		}
	    }
	}
      if (ReverseHopping)
	offset+=sprintf(auxArguments+offset,"rh_");
      if (Delta!=0.0)
	sprintf(deltaString,"d_%g_",Delta);
      if (Random!=0.0)
	sprintf(deltaString,"R_%g_",Random);
      if ((HardCore)&&(LandauQuantization=='x'))
	{
	  //sprintf(interactionStr,"_qx_hardcore");
	  // LandauQuantization not defined yet for hardcore bosons!
	  sprintf(interactionStr,"_hardcore");
	}
      else
	{
	  if (HardCore)
	    sprintf(interactionStr,"_hardcore");
	  else
	    {
	      if (LandauQuantization=='x')
		sprintf(interactionStr,"_qx_u_%g", ContactU);
	      else
		sprintf(interactionStr,"_u_%g", ContactU);	      
	    }
	  if (ContactW != 0.0)
	    {
	      sprintf(interactionStr,"%s_v%g", interactionStr, ContactW);
	    }
	}
      if ((SolenoidX!=0.0)||(SolenoidY!=0.0))
	    {
	      sprintf(interactionStr,"%s_s_%g_%g",interactionStr,SolenoidX,SolenoidY);
	    }
      int offset2;
      offset2 = sprintf(OutputName,"bosons_lattice_");
      if (Manager.GetBoolean("cylinder"))
	  offset2 += sprintf(OutputName+offset2,"cyl_");
      if (NbrFluxValues == 1)
	sprintf (OutputName+offset2, "n_%d_x_%d_y_%d%s_%s%sq_%d.dat", NbrBosons, Lx, Ly, interactionStr, auxArguments, deltaString, NbrFluxQuanta);
      else
	sprintf (OutputName+offset2, "n_%d_x_%d_y_%d%s_%s%sq.dat", NbrBosons, Lx, Ly, interactionStr, auxArguments, deltaString);
    }
  ParticleOnLattice* Space;
  if (NbrBosons>1)
    {
      if (HardCore)
	Space =new HardCoreBosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace, SolenoidX, SolenoidY, NbrLayers);
      else Space = new BosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace, SolenoidX, SolenoidY, LandauQuantization, NbrLayers);
    }
  else
    {
      if (!Manager.GetBoolean("no-single-particle-basis"))	
	Space = new SingleParticleOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace, SolenoidX, SolenoidY, LandauQuantization, NbrLayers);      
      else
       Space = new BosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace, SolenoidX, SolenoidY, LandauQuantization, NbrLayers);
    }
      
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    Memory = Architecture.GetArchitecture()->GetLocalMemory();
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  
  AbstractQHEOnLatticeHamiltonian* Hamiltonian;
  if (Manager.GetDouble("Kapit-Mueller")<=0.0)
    Hamiltonian = new ParticleOnLatticeDeltaHamiltonian(Space, NbrBosons, Lx, Ly, NbrFluxQuanta, ContactU, ReverseHopping, Delta,
							Random, Architecture.GetArchitecture(), NbrBody, Memory, LoadPrecalculationFileName,
							!Manager.GetBoolean("no-hermitian"), Manager.GetBoolean("cylinder"));  
  else
    {
      if (NbrLayers>1)
	Hamiltonian = new ParticleOnLatticeKapitMuellerMultiLayerHamiltonian(Space, NbrBosons, Lx, Ly, NbrFluxQuanta, ContactU, ContactW, ReverseHopping, Delta,
									     Random, Manager.GetDouble("Kapit-Mueller"), NbrCuts, branchCuts, branchShift, Architecture.GetArchitecture(), NbrBody, Memory, LoadPrecalculationFileName,
								 !Manager.GetBoolean("no-hermitian"));
      else
	Hamiltonian = new ParticleOnLatticeKapitMuellerHamiltonian(Space, NbrBosons, Lx, Ly, NbrFluxQuanta, ContactU, ReverseHopping, Delta,
								   Random, Manager.GetDouble("Kapit-Mueller"), Architecture.GetArchitecture(), NbrBody, Memory, LoadPrecalculationFileName,
								   !Manager.GetBoolean("no-hermitian"));
    }


	/*int NbrDiagonalInteractionFactors = NbrSites;
  double * DiagonalInteractionFactors = new double[NbrDiagonalInteractionFactors];
  int * DiagonalQValues=new int[NbrDiagonalInteractionFactors];
  for (int i=0; i<NbrDiagonalInteractionFactors; ++i)
	{
	  DiagonalQValues[i]=i;
	  DiagonalInteractionFactors[i]=ContactU;
	}
	
	for (int i = 0 ; i < Space->GetHilbertSpaceDimension(); i ++)
	{
		cout <<"i = "<<i<<" "<<Space->PrintState(cout, i)<<" "<<Space->ProdAdProdADiagonal(i,2, NbrDiagonalInteractionFactors,
							   DiagonalInteractionFactors, DiagonalQValues)<<" "<<Space->AdAdAADiagonal(i, NbrDiagonalInteractionFactors,
							 DiagonalInteractionFactors, DiagonalQValues) <<endl;

	}*/
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


  for (int iter=0; iter<NbrFluxValues; ++iter, ++NbrFluxQuanta)
    {
      cout << "----------------------------------------------------------------" << endl;
      cout << "NbrFluxQuanta="<<NbrFluxQuanta<<endl;
      
      if (!FirstRun) Hamiltonian->SetNbrFluxQuanta(NbrFluxQuanta);
  
      char* EigenvectorName = 0;
      if ((Manager.GetBoolean("eigenstate")||(Manager.GetBoolean("optimize-condensate"))))
	{
	  EigenvectorName = new char [64];
	  int offset;
	  offset = sprintf(EigenvectorName,"bosons_lattice_");
	  if (Manager.GetBoolean("cylinder"))
	      offset += sprintf(EigenvectorName+offset,"cyl_");
	  sprintf (EigenvectorName+offset, "n_%d_x_%d_y_%d%s_%s%sq_%d", NbrBosons, Lx, Ly, interactionStr, auxArguments, deltaString, NbrFluxQuanta);
	      }
      if (Manager.GetBoolean("optimize-condensate"))
	{
	  char *ParameterName = new char[strlen(EigenvectorName)+10];
	  sprintf(ParameterName,"%s.cond.par",EigenvectorName);
	  sprintf(EigenvectorName,"%s.cond.vec",EigenvectorName);
	  GutzwillerOnLatticeWaveFunction Condensate(NbrBosons, HardCore, Space);
	  Condensate.SetHamiltonian(Hamiltonian);
	  Condensate.SetArchitecture(Architecture.GetArchitecture());
	  Condensate.SetToRandomPhase();
	  int MaxEval = NbrSites*(NbrBosons+1)*2*Manager.GetInteger("nbr-iter");
	  double Energy=Condensate.Optimize(Manager.GetDouble("tolerance"), MaxEval);
	  Condensate.GetLastWaveFunction().WriteVector(EigenvectorName);
	  Condensate.GetVariationalParameters().WriteVector(ParameterName);
	  cout << "Found condensate state with energy: "<<Energy<<endl<<EigenvectorName<<endl;
	  delete [] ParameterName;
	}
      else
	{
	  QHEOnLatticeMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, NbrFluxQuanta, 0.0, OutputName, FirstRun, EigenvectorName);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
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
  return 0;
}
