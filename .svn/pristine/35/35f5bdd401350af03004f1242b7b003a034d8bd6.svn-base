#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "MathTools/BinomialCoefficients.h"

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

  // some running options and help
  OptionManager Manager ("FQHESphereFermionsTwoLandauLevels" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);
  ParticleOnSphereManager ParticleManager(true, false, 1);
  ParticleManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");

  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('\n', "level-difference", "Landau level index difference", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new SingleDoubleOption  ('e', "cyclotron-energy", "cyclotron energy in e^2/(epsilon l_b) unit ", 0.0);

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
      cout << "see man page for option syntax or type FQHESphereFermionsTwoLandauLevels -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  if (ULONG_MAX>>20 < (unsigned long) Manager.GetInteger("memory"))
    cout << "Warning: integer overflow in memory request - you might want to use 64 bit code."<<endl;
  unsigned long Memory = Manager.GetInteger("memory");
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool FirstRun = true;
  double* PseudoPotentials = 0;
  double* OneBodyPotentials = 0;
  int LandauLevelIndexDifference = Manager.GetInteger("level-difference");
  int LzMaxUp = LzMax + (2 * LandauLevelIndexDifference);
  int LzMaxDown = LzMax;
  double CyclotronEnergy = Manager.GetDouble("cyclotron-energy");


//   if (((SingleStringOption*) Manager["interaction-file"])->GetString() == 0)
//     {
//       cout << "an interaction file has to be provided" << endl;
//       return -1;
//     }
//   else
//     {
//       ConfigurationParser InteractionDefinition;
//       if (InteractionDefinition.Parse(((SingleStringOption*) Manager["interaction-file"])->GetString()) == false)
// 	{
// 	  InteractionDefinition.DumpErrors(cout) << endl;
// 	  return -1;
// 	}
//       int TmpNbrPseudoPotentials;
//       if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', PseudoPotentials, TmpNbrPseudoPotentials) == false)
// 	{
// 	  cout << "Weights is not defined or has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
// 	  return -1;
// 	}
//       cout << "LzMax= " << TmpNbrPseudoPotentials << " " << LzMax << endl;
//       if (TmpNbrPseudoPotentials != (LzMax +1))
// 	{
// 	  cout << "Invalid number of pseudo-potentials" << endl;
// 	  return -1;	  
// 	}
//     }

//   char* OutputNameLz = new char [256 + strlen(((SingleStringOption*) Manager["interaction-name"])->GetString())];
//   sprintf (OutputNameLz, "fermions_%s_n_%d_2s_%d_lz.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrParticles, LzMax);

  int TmpNbrParticleDown = (NbrParticles - LandauLevelIndexDifference) / 2;
  int Max = (((LzMax - TmpNbrParticleDown + 1) * TmpNbrParticleDown) +
	     ((LzMax + (2 * LandauLevelIndexDifference) - (NbrParticles - TmpNbrParticleDown) + 1) * (NbrParticles - TmpNbrParticleDown)));

  int  L = InitialLz;
  if (L < -Max)
    L = -Max;
  else
    if (L > Max)
      L = Max;
  if ((abs(Max) & 1) != (abs(InitialLz) & 1))
    L += 1;
  else
    {
      if (NbrLz > 0)
	{
	  Max = L + (2 * (NbrLz - 1));
	}
    }
  int TotalDim = 0;

  for (; L <= Max; L += 2)
    {

      ParticleOnSphere* Space = new FermionOnSphereTwoLandauLevels (NbrParticles, L, LzMaxUp, LzMaxDown);
//       for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
// 	{
// 	  cout << i << " : ";
// 	  Space->PrintState(cout, i) << endl;
// 	}

//       cout << L << " : " <<  Space->GetHilbertSpaceDimension() << endl;
//       if (L == 0)
// 	TotalDim += Space->GetHilbertSpaceDimension();
//       else
// 	TotalDim += (2 * Space->GetHilbertSpaceDimension());

      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      //      AbstractQHEOnSphereWithSpinHamiltonian* Hamiltonian = 0;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
 	Memory = Architecture.GetArchitecture()->GetLocalMemory();
//       if (OneBodyPotentials == 0)
// 	Hamiltonian = new ParticleOnSphereGenericHamiltonian(Space, NbrParticles, LzMax, PseudoPotentials,
// 							     ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(),
// 							     Architecture.GetArchitecture(), 
// 							     Memory, DiskCacheFlag,
// 							     LoadPrecalculationFileName);
//       else
// 	Hamiltonian = new ParticleOnSphereGenericHamiltonian(Space, NbrParticles, LzMax, PseudoPotentials, OneBodyPotentials,
// 							     ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(),
// 							     Architecture.GetArchitecture(), 
// 							     Memory, DiskCacheFlag,
// 							     LoadPrecalculationFileName);

//       double Shift = - 0.5 * ((double) (NbrParticles * NbrParticles)) / (0.5 * ((double) LzMax));
      
//       if (Manager.GetString("energy-expectation") != 0 )
// 	{
// 	  char* StateFileName = Manager.GetString("energy-expectation");
// 	  if (IsFile(StateFileName) == false)
// 	    {
// 	      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
// 	      return -1;           
// 	    }
// 	  RealVector State;
// 	  if (State.ReadVector(StateFileName) == false)
// 	    {
// 	      cout << "error while reading " << StateFileName << endl;
// 	      return -1;
// 	    }
// 	  if (State.GetVectorDimension()!=Space->GetHilbertSpaceDimension())
// 	    {
// 	      cout << "error: vector and Hilbert-space have unequal dimensions"<<endl;
// 	      return -1;
// 	    }
// 	  RealVector TmpState(Space->GetHilbertSpaceDimension());
// 	  VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
// 	  Operation.ApplyOperation(Architecture.GetArchitecture());
// 	  double EnergyValue = State*TmpState;
// 	  cout << "< Energy > = "<<EnergyValue<<endl;
// 	  cout << "< shifted energy > = "<<EnergyValue + Shift<<endl;
// 	  return 0;
// 	}
      
//       Hamiltonian->ShiftHamiltonian(Shift);
//       char* EigenvectorName = 0;
//       if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
// 	{
// 	  EigenvectorName = new char [64];
// 	  sprintf (EigenvectorName, "fermions_%s_n_%d_2s_%d_lz_%d", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrParticles, LzMax, L);
// 	}
      
//       QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
//       MainTaskOperation TaskOperation (&Task);
//       TaskOperation.ApplyOperation(Architecture.GetArchitecture());
//       if (EigenvectorName != 0)
// 	{
// 	  delete[] EigenvectorName;
// 	}

//       delete Hamiltonian;
//       if (FirstRun == true)
// 	FirstRun = false;
     }
//   cout << "total dim = " << TotalDim << endl;
//   BinomialCoefficients TmpCoef(2 * (LzMax + 1 + LandauLevelIndexDifference));
//   cout << TmpCoef(2 * (LzMax + 1 + LandauLevelIndexDifference), NbrParticles) << endl;
  return 0;
}
