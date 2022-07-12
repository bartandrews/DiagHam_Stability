#include "HilbertSpace/FermionOnLatticeGenericMomentumSpace.h"
#include "HilbertSpace/BosonOnLatticeGenericMomentumSpace.h"
#include "Hamiltonian/ParticleOnLatticeGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeExternalHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "MainTask/GenericComplexMainTask.h"

#include "LanczosAlgorithm/LanczosManager.h"

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

  OptionManager Manager ("FQHELatticeGenericMomentum" , "0.01");  
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  LanczosManager Lanczos(true);  
  Lanczos.AddOptionGroup(&Manager);

  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 8);
  (*SystemGroup) += new SingleIntegerOption  ('x', "length-x", "number of k-points in x-direction", 4);
  (*SystemGroup) += new SingleIntegerOption  ('y', "length-y", "number of k-points in y-direction", 4);
  (*SystemGroup) += new SingleIntegerOption  ('b', "nbr-bands", "number of bands", 1);
  (*SystemGroup) += new SingleIntegerOption  ('X', "x-momentum", "constraint on the total momentum in the x direction (negative if all)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('Y', "y-momentum", "constraint on the total momentum in the y direction (negative if all)", -1);

  (*SystemGroup) += new BooleanOption('f',"fermions","Assume fermionic statistics");
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons");
  
  (*SystemGroup) += new SingleStringOption  ('v', "external-two-body", "use definition of two-body interactions from a file");
  (*SystemGroup) += new SingleStringOption  ('V', "external-name", "descriptor of external interaction (if in use)");
  (*SystemGroup) += new SingleStringOption  ('u', "one-body", "use definition of one-body potential/hoppings from a file");
  (*SystemGroup) += new SingleStringOption  ('U', "potential-name", "descriptor of external potential");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hamiltonian-shift", "additional shift of Hamiltonian",0.0);

  (*SystemGroup) += new BooleanOption  ('\n', "all-points", "calculate all points in BZ", false);

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
  (*MiscGroup) += new SingleDoubleOption('\n', "tolerance", "tolerance for variational parameters in condensate",1e-6);
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
  (*MiscGroup) += new BooleanOption  ('\n', "test-hamiltonian", "test hermiticity of Hamiltonian");
  (*MiscGroup) += new BooleanOption  ('\n', "show-basis", "print basis states and quit");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrParticles = Manager.GetInteger("nbr-particles");
  int Nx = Manager.GetInteger("length-x");
  int Ny = Manager.GetInteger("length-y");
  int Kx = Manager.GetInteger("x-momentum");
  int Ky = Manager.GetInteger("y-momentum");
  int NbrBands = Manager.GetInteger("nbr-bands");
  int NbrStates = Nx*Ny*NbrBands;
  bool FermionicStatistics = Manager.GetBoolean("fermions");
  bool HardCore = Manager.GetBoolean("hard-core");
  double Shift = Manager.GetDouble("hamiltonian-shift");

  if (ULONG_MAX>>20 < (unsigned long)Manager.GetInteger("memory"))
    cout << "Warning: integer overflow in memory request - you might want to use 64 bit code."<<endl;

  unsigned long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  bool FirstRun = true;

  if (Manager.GetString("energy-expectation") != 0 ) Memory = 0x0ul;

  char* OutputName;
  char* OutputBase;
  if ( (OutputName = Manager.GetString("output-file")) == NULL)
    {
      OutputName = new char [1024];
      OutputBase = new char [1024];
      char *InteractionStr = new char[100];
      InteractionStr[0]='\0';
      char *StatisticsStr = new char[10];
      if (Manager.GetString("external-two-body")!=NULL)
	{
	  if (Manager.GetString("external-name")==NULL)
	    {
	      if (HardCore)
		sprintf(InteractionStr,"_ext_hardcore");
	      else
		sprintf(InteractionStr,"_ext");
	    }
	  else
	    strcpy(InteractionStr,Manager.GetString("external-name"));	  
	}
      if (Manager.GetString("one-body")!=NULL)
	{
	  if (Manager.GetString("potential-name")==NULL)
	    {
	      sprintf(InteractionStr,"%s_pot",InteractionStr);
	    }
	  else
	    sprintf(InteractionStr,"%s_%s",InteractionStr,Manager.GetString("potential-name"));	  
	}
      if (Manager.GetBoolean("fermions"))
	sprintf(StatisticsStr,"fermions");
      else
	sprintf(StatisticsStr,"bosons");
      sprintf(OutputBase,"%s_n_%d_x_%d_y_%d_b_%d%s.dat", StatisticsStr, NbrParticles, Nx, Ny, NbrBands, InteractionStr);
    }
  else
    {
      OutputBase = RemoveExtensionFromFileName(OutputName, ".dat");
      if (OutputBase==0)
	{
	  OutputBase = new char[strlen(OutputName)+1];
	  strcpy(OutputBase,OutputName);
	}
    }

  int XMaxMomentum = Nx-1;
    bool GenerateMomenta = false;
  if ((Kx < 0)||(Ky < 0))
    GenerateMomenta = true;
  if (Kx < 0)
    Kx = 0;
  else
    XMaxMomentum = Kx;
  int YMaxMomentum = (Ny - 1);
  if (Ky < 0)
    Ky = 0;
  else
    YMaxMomentum = Ky;

  int NbrMomenta;
  int *XMomenta;
  int *YMomenta;
  //  int *Multiplicities = NULL;
  //  int CenterX=0, CenterY=0;

  if (GenerateMomenta==false)
    {
      NbrMomenta=1;
      XMomenta = new int[1];
      YMomenta = new int[1];
      XMomenta[0]=Kx;
      YMomenta[0]=Ky;
    }
  else
    {
      // if (Manager.GetBoolean("all-points"))
	{
	  int Pos=0;
	  NbrMomenta = (XMaxMomentum-Kx+1)*(YMaxMomentum-Ky+1);
	  XMomenta = new int[NbrMomenta];
	  YMomenta = new int[NbrMomenta];
	  for (; Kx <= XMaxMomentum; ++Kx)
	    for (int Ky2 = Ky; Ky2<= YMaxMomentum; ++Ky2)
	      {
		XMomenta[Pos]=Kx;
		YMomenta[Pos]=Ky2;
		++Pos;
		//cout << "Pos="<<Pos<<endl;
	      }
	}
	//else {}
    }

  if (Manager.GetString("output-file")==NULL)
    {
      if (NbrMomenta>1)
	sprintf(OutputName,"%s.dat",OutputBase);
      else
	sprintf(OutputName,"%s_kx_%d_ky_%d.dat", OutputBase, XMomenta[0], YMomenta[0]);
    }
  
  char *SubspaceStr = new char[50];  
  char *SubspaceLegend= new char[10];
  sprintf(SubspaceLegend,"kx ky");

  for (int Pos=0; Pos<NbrMomenta; ++Pos)
    {
      sprintf(SubspaceStr,"%d %d", XMomenta[Pos], YMomenta[Pos]);

      ParticleOnLattice* Space = NULL;
      if (FermionicStatistics)
	{
	  Space = new FermionOnLatticeGenericMomentumSpace(NbrParticles, Nx, Ny, XMomenta[Pos], YMomenta[Pos], NbrBands, MemorySpace);
	}
      else
	{
	  Space = new BosonOnLatticeGenericMomentumSpace(NbrParticles, Nx, Ny, XMomenta[Pos], YMomenta[Pos], NbrBands, MemorySpace);
	}
      
      if (Manager.GetBoolean("show-basis"))
	for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i)
	  Space->PrintState(cout,i)<<endl;
      
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      
      AbstractQHEOnLatticeHamiltonian* Hamiltonian;
      Hamiltonian = new ParticleOnLatticeExternalHamiltonian(Space, NbrParticles, NbrStates, Manager.GetString("one-body"),
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
      cout << "----------------------------------------------------------------" << endl;
      cout << "K=["<<XMomenta[Pos]<<","<<YMomenta[Pos]<<"]"<<endl;
      
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate"))
	{
	  EigenvectorName = new char [1024];
	  sprintf(EigenvectorName,"%s_kx_%d_ky_%d", OutputBase, XMomenta[Pos], YMomenta[Pos]);
	}
      GenericComplexMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, SubspaceStr, SubspaceLegend,
				   Shift, OutputName, FirstRun, EigenvectorName);
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
  delete [] OutputBase;
  delete [] OutputName;
  delete [] SubspaceStr;
  delete [] SubspaceLegend;
  return 0;
}

