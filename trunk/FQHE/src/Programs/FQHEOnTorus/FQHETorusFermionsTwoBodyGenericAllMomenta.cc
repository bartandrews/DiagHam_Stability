#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnTorus.h"

#include "Hamiltonian/ParticleOnTorusGenericHamiltonianAllMomenta.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/FQHEOnTorusMainTask.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Tools/FQHEFiles/FQHETorusPseudopotentialTools.h"

#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "QuantumNumber/AbstractQuantumNumber.h"

#include "GeneralTools/ListIterator.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

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


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHETorusFermionsTwoBodyGenericAllMomenta" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption ('l', "max-momentum", "maximum momentum for a single particle", 18);
  (*SystemGroup) += new SingleDoubleOption ('r', "ratio", "ratio between the two torus lengths", 1.0);
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new SingleStringOption ('\n', "confining-file", "file describing the confining potential");
  (*SystemGroup) += new BooleanOption  ('\n', "redundant-kymomenta", "Calculate all subspaces up to Ky  = MaxMomentum-1", false);
  (*SystemGroup) += new BooleanOption  ('\n', "mass-anisotropy", "use a mass anisotropy for the system");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "anisotropy", "value of the anisotropy parameter alpha (i.e. q_g^2 = alpha q_x^2 + q_y^2 / alpha)", 1.0);
  (*SystemGroup) += new SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
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
      cout << "see man page for option syntax or type FQHETorusFermionsTwoBodyGenericAllMomenta -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int NbrParticles = Manager.GetInteger("nbr-particles");
  int MaxMomentum = Manager.GetInteger("max-momentum");
  double XRatio = Manager.GetDouble("ratio");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  bool FirstRun = true;
  
  double* PseudoPotentials;
  int NbrPseudoPotentials = 0;
  double* OneBodyPotentials = 0;
  double** OffDiagonalOneBodyPotentials = 0;
  int MaximumMomentumTransfer = 0;

  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      if (FQHETorusGetPseudopotentials(Manager.GetString("interaction-file"), MaxMomentum, NbrPseudoPotentials, PseudoPotentials, OneBodyPotentials) == false)
	return -1;
    }
  if (Manager.GetString("confining-file") != 0)
    {
      MultiColumnASCIIFile ConfiningFile;
      if (ConfiningFile.Parse(Manager.GetString("confining-file")) == false)
	{
	  ConfiningFile.DumpErrors(cout);
	  return -1;
	}
      if ((ConfiningFile.GetNbrColumns() < 3) || (ConfiningFile.GetNbrLines() == 0))
	{
	  cout << "error, " << Manager.GetString("confining-file") << " has an invalid format" << endl;
	  return -1;
	}
      int* CreationIndices = ConfiningFile.GetAsIntegerArray(0);
      int* AnnihilationIndices = ConfiningFile.GetAsIntegerArray(1);
      double* TmpPotentials = ConfiningFile.GetAsDoubleArray(2);
      for (int i = 0; i < ConfiningFile.GetNbrLines(); ++i)
	{
	  if ((CreationIndices[i] < 0) || (AnnihilationIndices[i] < 0) || (CreationIndices[i] >= MaxMomentum) || (AnnihilationIndices[i] >= MaxMomentum))
	    {
	      cout << "invalid indices line " << (i+1) << endl;
	      return -1;
	    }
	  if (abs(CreationIndices[i] - AnnihilationIndices[i]) > MaximumMomentumTransfer)
	    {
	      MaximumMomentumTransfer = abs(CreationIndices[i] - AnnihilationIndices[i]);
	    }
	}
      OneBodyPotentials = new double[MaxMomentum + 1];
      for (int i = 0; i < MaxMomentum; ++i)
	{
	  OneBodyPotentials[i] = 0.0;
	}
      if (MaximumMomentumTransfer > 0)
	{
	  OffDiagonalOneBodyPotentials = new double*[MaxMomentum];
	  for (int i = 0; i < MaxMomentum; ++i)
	    {
	      OffDiagonalOneBodyPotentials[i] = new double[MaximumMomentumTransfer];
	      for (int j = 0; j < MaximumMomentumTransfer; ++j)	  
		{
		  OffDiagonalOneBodyPotentials[i][j] = 0.0;
		}
	    }	  
	}
      for (int i = 0; i < ConfiningFile.GetNbrLines(); ++i)
	{
	  int TmpMomentumTransfer = CreationIndices[i] - AnnihilationIndices[i];
	  if (TmpMomentumTransfer == 0)
	    {
	      OneBodyPotentials[CreationIndices[i]] = TmpPotentials[i];
	    }
	  else
	    {
	      if (TmpMomentumTransfer > 0)
		{
		  OffDiagonalOneBodyPotentials[CreationIndices[i]][TmpMomentumTransfer - 1] = TmpPotentials[i];
		}
	      else
		{
		  OffDiagonalOneBodyPotentials[CreationIndices[i]][-TmpMomentumTransfer - 1] = TmpPotentials[i];
		}
	    }
	}
    }

  char* OutputNamePrefix = new char [1024];
  if (Manager.GetBoolean("mass-anisotropy") == false)
    {
      sprintf (OutputNamePrefix, "fermions_torus_noky_%s_n_%d_2s_%d_ratio_%f", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum, XRatio);
    }
  else
    {
      sprintf (OutputNamePrefix, "fermions_torus_noky_%s_anisotropy_%f_n_%d_2s_%d_ratio_%f", Manager.GetString("interaction-name"), 
	       Manager.GetDouble("anisotropy"), NbrParticles, MaxMomentum, XRatio);
    }
  char* OutputNameLz = new char [strlen(OutputNamePrefix) + 8];
  sprintf (OutputNameLz, "%s.dat", OutputNamePrefix);
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  File.precision(14);

  cout << "----------------------------------------------------------------" << endl;
  cout << " Ratio = " << XRatio << endl;
  
  ParticleOnTorus* Space = new FermionOnTorus (NbrParticles, MaxMomentum);
  cout << " Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;
  
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    Memory = Architecture.GetArchitecture()->GetLocalMemory();
  
  AbstractQHEHamiltonian* Hamiltonian = 0;
  Hamiltonian = new ParticleOnTorusGenericHamiltonianAllMomenta (Space, NbrParticles, MaxMomentum, XRatio, NbrPseudoPotentials, PseudoPotentials,
								 OneBodyPotentials, 0.0, 0.0, Architecture.GetArchitecture(), Memory);
  double Shift = -10.0;
  Hamiltonian->ShiftHamiltonian(Shift);
  char* EigenvectorName = 0;
  if (Manager.GetBoolean("eigenstate") == true)	
    {
      EigenvectorName = new char [64 + strlen(OutputNamePrefix)];
      sprintf (EigenvectorName, "%s", OutputNamePrefix);
    }
  FQHEOnTorusMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, 0, Shift, OutputNameLz, FirstRun, EigenvectorName);
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
  File.close();

  return 0;
}
