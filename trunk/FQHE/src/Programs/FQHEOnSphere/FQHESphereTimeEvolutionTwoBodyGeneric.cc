#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereGenericHamiltonian.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"

#include "Operator/ParticleOnSphereDensityOperator.h"


#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"

#include "GeneralTools/ListIterator.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <cstring> 


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereTimeEvolutionTwoBodyGeneric" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  Manager += ToolsGroup;

  (*SystemGroup) += new SingleStringOption('\n', "initial-state", "name of the file containing the initial vector upon which e^{-iHt} acts");
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistics (default value))");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "complex", "initial vector is a complex vector");
  (*SystemGroup) += new BooleanOption  ('\n', "compute-energy", "compute the energy of each time-evolved vector");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "sma");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potentials");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  
  (*SystemGroup) += new SingleDoubleOption ('\n', "time-step", "time interval between two snap shots", 0.1);
  (*SystemGroup) += new SingleDoubleOption ('\n', "time-shift", "time shift to add in the vectors name, if initial state for current run is not at t = 0", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-steps", "number of points to evaluate", 1);
  
  (*SystemGroup) += new SingleDoubleOption ('\n', "precision", "convergence precision", 1.0e-14);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximal number of iterations", 100);
  
  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereTimeEvolutionTwoBodyGeneric -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz = Manager.GetInteger("total-lz");
  bool Statistics = true;
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();
  bool DiskCacheFlag = ((BooleanOption*) Manager["disk-cache"])->GetBoolean();
  
  double* PseudoPotentials0 = 0;
  double* PseudoPotentials;
  double* OneBodyPotentials = 0; 
    
  double TmpTime = Manager.GetDouble("time-step");
  int NbrTimeSteps = Manager.GetInteger("nbr-steps");
  double TimeShift = Manager.GetDouble("time-shift");
  
  if (Manager.GetString("initial-state") == 0)
  {
    cout << " Error: an initial state must be provided" << endl;
    return -1;
  }
  if (Manager.GetBoolean("boson") == true)
    {
      Statistics = false;
    }
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("initial-state"),
						   NbrParticles, LzMax, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("initial-state") << endl;
      return -1;
    }
  cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << LzMax << " Lz=" << TotalLz << " ";
  if (Statistics == false)
    {
      cout << "bosons";
    }
  else
    {
      cout << "fermions";
    }
  cout << endl;
  
  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      ConfigurationParser InteractionDefinition;
      if (InteractionDefinition.Parse(Manager.GetString("interaction-file")) == false)
	{
	  InteractionDefinition.DumpErrors(cout) << endl;
	  return -1;
	}
      int TmpNbrPseudoPotentials;
      if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', PseudoPotentials0, TmpNbrPseudoPotentials) == false)
	{
	  cout << "Weights is not defined or as a wrong value in " << Manager.GetString("interaction-file") << endl;
	  return -1;
	}
      if (TmpNbrPseudoPotentials > (LzMax +1))
	{
	  cout << "Invalid number of pseudo-potentials" << endl;
	  return -1;	  
	}
      else
	{
	  PseudoPotentials = new double[LzMax + 1];
	  for (int i = 0; i < TmpNbrPseudoPotentials; ++i)
	    PseudoPotentials[i] = PseudoPotentials0[i];
	  for (int i = TmpNbrPseudoPotentials; i < LzMax + 1; ++i)
	    PseudoPotentials[i] = 0.0;
	}
      if (InteractionDefinition.GetAsDoubleArray("Onebodypotentials", ' ', OneBodyPotentials, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials" << endl;
	      return -1;	  
	    }
	}
    }

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
 
  char* OutputNamePrefix = new char [512];
  char* NormName = new char[512];
  char* EnergyName;
  char* InteractionName = new char[strlen(Manager.GetString("interaction-name")) + 1];
  strcpy(InteractionName, Manager.GetString("interaction-name"));
  
  ParticleOnSphere* Space = 0;
  
  if (Statistics == false)
    {
#ifdef  __64_BITS__
      if ((LzMax + NbrParticles - 1) < 63)
#else
	if ((LzMax + NbrParticles - 1) < 31)	
#endif
	  {
	    Space = new BosonOnSphereShort(NbrParticles, TotalLz, LzMax);
	  }
	else
	  {
	    Space = new BosonOnSphere(NbrParticles, TotalLz, LzMax);
	  }
      sprintf (OutputNamePrefix, "bosons_sphere_%s_n_%d_2s_%d_t", Manager.GetString("interaction-name"), NbrParticles, LzMax);
      sprintf (NormName, "bosons_sphere_%s_n_%d_2s_%d_dt_%g_t0_%g_nbrsteps_%d_norm.dat", Manager.GetString("interaction-name"), NbrParticles, LzMax, TmpTime, TimeShift, NbrTimeSteps);
      if (Manager.GetBoolean("compute-energy"))
      {
	EnergyName = new char[512];
	sprintf (EnergyName, "bosons_sphere_%s_n_%d_2s_%d_dt_%g_t0_%g_nbrsteps_%d_energy.dat", InteractionName, NbrParticles, LzMax, TmpTime, TimeShift, NbrTimeSteps);
      }
      
    }
  else
    {
      Space = new FermionOnSphere (NbrParticles, TotalLz, LzMax);
      sprintf (OutputNamePrefix, "fermions_sphere_%s_n_%d_2s_%d_t", Manager.GetString("interaction-name"), NbrParticles, LzMax);
      sprintf (NormName, "fermions_sphere_%s_n_%d_2s_%d_dt_%g_t0_%g_nbrsteps_%d_norm.dat", Manager.GetString("interaction-name"), NbrParticles, LzMax, TmpTime, TimeShift, NbrTimeSteps);
      if (Manager.GetBoolean("compute-energy"))
      {
	EnergyName = new char[512];
	sprintf (EnergyName, "fermions_sphere_%s_n_%d_2s_%d_dt_%g_t0_%g_nbrsteps_%d_energy.dat", InteractionName, NbrParticles, LzMax, TmpTime, TimeShift, NbrTimeSteps);
      }
    }

  delete[] InteractionName;
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

  char* StateFileName = Manager.GetString("initial-state");
  if (IsFile(StateFileName) == false)
    {
      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }

  ComplexVector TmpInitialState (Space->GetHilbertSpaceDimension());
  if (Manager.GetBoolean("complex") == false)
  {
    RealVector InputState;
    if (InputState.ReadVector(StateFileName) == false)
    {
      cout << "error while reading " << StateFileName << endl;
      return -1;
    }
    if (InputState.GetVectorDimension() != Space->GetHilbertSpaceDimension())
    {
      cout << "error: vector and Hilbert-space have unequal dimensions " << InputState.GetVectorDimension() << " "<< Space->GetHilbertSpaceDimension() << endl;
      return -1;
    }
    TmpInitialState = InputState;
  }
  else
    {
    ComplexVector InputState;
    if (InputState.ReadVector(StateFileName) == false)
    {
      cout << "error while reading " << StateFileName << endl;
      return -1;
    }
    if (InputState.GetVectorDimension() != Space->GetHilbertSpaceDimension())
    {
      cout << "error: vector and Hilbert-space have unequal dimensions " << InputState.GetVectorDimension() << " "<< Space->GetHilbertSpaceDimension() << endl;
      return -1;
    }
    TmpInitialState = InputState;
  }
  
  
    
  cout << " Initial state Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;
  
  
  
  ParticleOnSphereGenericHamiltonian* Hamiltonian = 0;
  
   if (OneBodyPotentials == 0)
        Hamiltonian = new ParticleOnSphereGenericHamiltonian(Space, NbrParticles, LzMax, PseudoPotentials,
							   Manager.GetDouble("l2-factor"),
							   Architecture.GetArchitecture(), 
							   Memory, DiskCacheFlag,
							   LoadPrecalculationFileName);
   else
        Hamiltonian = new ParticleOnSphereGenericHamiltonian(Space, NbrParticles, LzMax, PseudoPotentials, OneBodyPotentials, 
							   Manager.GetDouble("l2-factor"),
							   Architecture.GetArchitecture(), 
							   Memory, DiskCacheFlag,
							   LoadPrecalculationFileName);
											 
//    double Shift = - 0.5 * ((double) (NbrParticles * NbrParticles)) / (0.5 * ((double) LzMax));
	
	

  ofstream File;
  File.open(NormName, ios::binary | ios::out);
  File.precision(14);
  File << "# t Norm dNorm" << endl;
  
  ofstream FileEnergy;
  if (Manager.GetBoolean("compute-energy"))
  {
    FileEnergy.open(EnergyName, ios::binary | ios::out);
    FileEnergy.precision(14);
    FileEnergy << "# t E "<< endl;
  }
  
  double Norm;
  int TmpExpansionOrder;
  ComplexVector TmpState (Space->GetHilbertSpaceDimension()) ;
  ComplexVector TmpState1 (Space->GetHilbertSpaceDimension()) ;
  Complex TmpCoefficient;
  for (int i = 0; i < NbrTimeSteps; ++i)
  {
    TmpState.Copy(TmpInitialState);
    Norm = TmpState.Norm();
    double TmpNorm = 1.0;
    TmpExpansionOrder = 0;
    TmpCoefficient = 1.0;
    cout << "Computing state " << (i + 1) << "/" << NbrTimeSteps << " at t = " << (TmpTime * i) << endl;
    while (((fabs(TmpNorm) > Manager.GetDouble("precision")) || (TmpExpansionOrder < 1)) && (TmpExpansionOrder <= Manager.GetInteger("iter-max")))
    {
      TmpExpansionOrder += 1;
      TmpCoefficient = -TmpCoefficient * TmpTime * Complex(0.0, 1.0) / ((double) TmpExpansionOrder);
      VectorHamiltonianMultiplyOperation Operation (Hamiltonian, (&TmpState), (&TmpState1));
      Operation.ApplyOperation(Architecture.GetArchitecture());
      TmpState.Copy(TmpState1);
      TmpNorm = sqrt(TmpCoefficient.Re*TmpCoefficient.Re + TmpCoefficient.Im*TmpCoefficient.Im) * TmpState.Norm();
      AddComplexLinearCombinationOperation Operation1 (&TmpInitialState, &TmpState1, 1, &TmpCoefficient);
      Operation1.ApplyOperation(Architecture.GetArchitecture());
      Norm = TmpInitialState.Norm();
      
      cout << "Norm = " << Norm << " +/- " << TmpNorm << " for step " << TmpExpansionOrder << endl;      
    }
    File << (TimeShift + (i + 1)*TmpTime) << " " << Norm << " " << TmpNorm << endl;
    cout << endl;  
    
    char* OutputName = new char [strlen(OutputNamePrefix)+ 16];
    sprintf (OutputName, "%s_%g_lz_%d.vec", OutputNamePrefix, TimeShift + (i + 1)*TmpTime, TotalLz);
    TmpInitialState.WriteVector(OutputName);
    delete[] OutputName;
    
    if (Manager.GetBoolean("compute-energy"))
    {
      VectorHamiltonianMultiplyOperation Operation (Hamiltonian, (&TmpInitialState), (&TmpState));
      Operation.ApplyOperation(Architecture.GetArchitecture());
      Complex Energy = (TmpInitialState) * (TmpState);
      FileEnergy << (TimeShift + (i + 1)*TmpTime) << " " << Energy.Re << endl;
      cout << "E = " << Energy.Re << endl;
    }
  }
  
  File.close();
  if (Manager.GetBoolean("compute-energy"))
  {
    FileEnergy.close();
    delete[] EnergyName;
  }
  delete Hamiltonian;
  delete[] OutputNamePrefix;
  delete[] NormName;
  
 
  return 0;
}
