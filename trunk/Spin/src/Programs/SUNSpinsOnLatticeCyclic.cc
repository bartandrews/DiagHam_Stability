#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"

#include "HilbertSpace/GenericSUNSpinCollection.h"
#include "Hamiltonian/SUNSpinOnLatticeQuadraticHamiltonian.h"
#include "Tools/LatticeConnections.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "GeneralTools/FilenameTools.h"

#include "MainTask/GenericComplexMainTask.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <iostream>
#include <limits.h>

using std::cout;
using std::endl;

// recursively count number of cartan quantum numbers
// nbrSpins = number of spins
// level = remaining cartan quantum numbers
int CountCartanQuantumNumbers(int nbrSpins, int level, int maxPerLevel)
{
  //cout << "CountCartanQuantumNumbers("<< nbrSpins<<", "<< level<<", "<< maxPerLevel<<")"<<endl;
  int sum=0;
  if (level>2)
    {
      for (int N=(int)ceil((double)nbrSpins/level); N<=maxPerLevel; ++N)
	sum+=CountCartanQuantumNumbers(nbrSpins-N,level-1, N);
    }
  else // level = 2: only one degree of freedom left
    {
      for (int N=(int)ceil((double)nbrSpins/2.0); (N<=maxPerLevel) && (N<=nbrSpins); ++N)
	++sum;
    }
  return sum;
}

// recursively count number of cartan quantum numbers
// cartanQuantumNumbers = array where to store quantum numbers
// maxPos = maximum number of sets to generate
// levelN = degree of SUN spins
// nbrSpins = number of spins
// level = remaining cartan quantum numbers
// pos = place where to start storing things
// return = total number of quantum numbers
int GenerateCartanQuantumNumbers(int **cartanQuantumNumbers, int maxPos, int maxPerLevel, int levelN, int nbrSpins, int level, int pos)
{
  if (level>2)
    {      
      for (int N=(int)ceil((double)nbrSpins/level); (N<=maxPerLevel) && (pos<maxPos) ; ++N)
	{
	  int TmpPos = pos;
	  pos = GenerateCartanQuantumNumbers(cartanQuantumNumbers, maxPos, N, levelN, nbrSpins-N, level-1, TmpPos);
	  for (int i=TmpPos; i<pos; ++i)
	    cartanQuantumNumbers[i][levelN-level]=N;
	}
    }
  else // level = 2: only one degree of freedom left
    {
      for (int N=(int)ceil((double)nbrSpins/2.0); (N<=maxPerLevel) && (N<=nbrSpins)&& (pos<maxPos) ; ++N)
	{
	  cartanQuantumNumbers[pos]=new int[levelN];
	  cartanQuantumNumbers[pos][levelN-2]=N;
	  cartanQuantumNumbers[pos++][levelN-1]=nbrSpins-N;
	}
    }
  return pos;
}


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("SpinOnLatticeWithSUN" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");      
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  LatticeConnections::AddOptionGroup(&Manager);
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleIntegerOption  ('N', "level-n", "level of SU(N) symmetry group", 3);
  (*SystemGroup) += new MultipleIntegerOption  ('t', "cartan", "eigenvalues of the generators of the cartan algebra",',');
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-sectors", "number of cartan-sectors to be calculated", 1);
  (*SystemGroup) += new SingleDoubleOption  ('c', "cyclic", "prefactor of cyclic permutation operators around plaquettes", 1.0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout); 

  int LevelN=Manager.GetInteger("level-n");  
  unsigned long MemorySpace = ((unsigned long)Manager.GetInteger("fast-search")) << 20;
  if (ULONG_MAX>>20 < (unsigned long)Manager.GetInteger("memory"))
    cout << "Warning: integer overflow in memory request - you might want to use 64 bit code."<<endl;
  unsigned long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  char *LoadPrecalculationFileName = Manager.GetString("load-precalculation");

  int NbrSectors=Manager.GetInteger("nbr-sectors");
  int NbrCartan;
  int *CartanQuantumNumbers;
  int **AllCartanQuantumNumbers=NULL;
  double CyclicTerms = Manager.GetDouble("cyclic");
  CartanQuantumNumbers = Manager.GetIntegers("cartan",NbrCartan);  
  cout << "NbrCartan="<<NbrCartan<<", NbrSectors="<<NbrSectors <<endl;
  LatticeConnections *Lattice = new LatticeConnections();
  
  int NbrSpins=Lattice->GetNbrSites();
  cout << "NbrSpins="<<NbrSpins<<endl;
  if (NbrCartan>0)
    {
      if (NbrCartan<LevelN-1)
	{
	  cout << "Please select the desired subspace!"<<endl
	       << "Use option -t to indicate the first N-1 eigenvalues of the Cartan operators"<<endl;
	  exit(-1);
	}
      if (NbrCartan>LevelN-1)
	cout << "More values of Cartan operators than needed. Will ignore all beyond N-1"<<endl;  
      if (NbrSectors>1)
	{	  
	  cout << "Attention, evaluating a single sector only"<<endl;
	}
      else
	{
	  AllCartanQuantumNumbers = new int*[1];
	  AllCartanQuantumNumbers[0]=CartanQuantumNumbers;
	}
    }
  else
    {
      int MaxNbrSectors = CountCartanQuantumNumbers(NbrSpins, LevelN, NbrSpins);
      if (NbrSectors<0)
	NbrSectors=MaxNbrSectors;
      else if (NbrSectors>MaxNbrSectors)
	NbrSectors=MaxNbrSectors;
      AllCartanQuantumNumbers = new int*[NbrSectors];
      int TmpSectors = GenerateCartanQuantumNumbers(AllCartanQuantumNumbers, NbrSectors, NbrSpins, LevelN,
							     NbrSpins, LevelN, 0);
      CartanQuantumNumbers = AllCartanQuantumNumbers[0];
      if (TmpSectors!=NbrSectors)
	{
	  cout << "Problem with generation of quantum numbers"<<endl;
	}
    }

  bool FirstRun = true;

  double Shift=0.0;

  char *SubspaceStr = new char[50];  
  char *SubspaceLegend= new char[50];
  sprintf(SubspaceLegend,"C0");
  for (int n=1; n < LevelN; ++n)
    sprintf(SubspaceLegend,"%s C%d", SubspaceLegend, n);
  
  
  char *Geometry=Lattice->GeometryString();
  char *OutputFileName = new char[100];

  sprintf (OutputFileName, "spins_SU%d_%s_cyc_%g_n_%d_c", LevelN, Geometry, CyclicTerms, NbrSpins);

  // if there is a single set of cartan quantum-numbers -> add to filename
  if (NbrSectors==1)
    {
      for (int n=0; n < LevelN; ++n)
	sprintf(OutputFileName,"%s_%d", OutputFileName, CartanQuantumNumbers[n]);
      sprintf(OutputFileName,"%s.dat", OutputFileName);
    }
  else
    {
      sprintf(OutputFileName,"%s.dat", OutputFileName);
    }

  // [loop over cartan quantum numbers]
  for (int sec=0; sec<NbrSectors; ++sec) 
  {
    CartanQuantumNumbers=AllCartanQuantumNumbers[sec];
    sprintf(SubspaceStr,"%d", CartanQuantumNumbers[0]);
    for (int n=1; n < LevelN; ++n)
      sprintf(SubspaceStr,"%s %d", SubspaceStr, CartanQuantumNumbers[n]);
    
    GenericSUNSpinCollection *Space = new GenericSUNSpinCollection(LevelN, NbrSpins, CartanQuantumNumbers, MemorySpace);
    
    Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
    
    AbstractSUNSpinHamiltonian *Hamiltonian =
      new SUNSpinOnLatticeQuadraticHamiltonian(Space, Lattice, Architecture.GetArchitecture(),
					       Memory, LoadPrecalculationFileName, CyclicTerms);


    char* EigenvectorName = 0;
    if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
      {
	EigenvectorName = new char [64];
	sprintf (EigenvectorName, "spins_SU%d_%s_cyc_%g_n_%d_c", LevelN, Geometry, CyclicTerms, NbrSpins);
	for (int n=0; n < LevelN-1; ++n)
	  sprintf(EigenvectorName,"%s_%d", EigenvectorName, CartanQuantumNumbers[n]);
      }

    GenericComplexMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, SubspaceStr, SubspaceLegend,
				 Shift, OutputFileName, FirstRun, EigenvectorName);
    MainTaskOperation TaskOperation (&Task);
    TaskOperation.ApplyOperation(Architecture.GetArchitecture());

    if (EigenvectorName != 0)
      {
	delete[] EigenvectorName;
      }
    
    if (FirstRun == true)
      FirstRun = false;
    
    delete Space;
    delete Hamiltonian;
    delete [] CartanQuantumNumbers;    
  }
  
  delete [] Geometry;
  delete [] OutputFileName;
  delete [] SubspaceLegend;
  delete [] SubspaceStr;
  delete [] AllCartanQuantumNumbers;
  delete Lattice;
}
