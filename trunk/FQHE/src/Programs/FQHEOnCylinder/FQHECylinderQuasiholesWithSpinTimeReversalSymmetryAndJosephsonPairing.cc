#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairing.h"
#include "HilbertSpace/QuasiholeOnSphereWithSpin.h"
#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairingAllMomenta.h"

#include "Hamiltonian/ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing.h"
#include "Hamiltonian/ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"

#include "MainTask/GenericRealMainTask.h"
#include "MainTask/GenericComplexMainTask.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <cstring>
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
  OptionManager Manager ("FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing" , "0.01");
  OptionGroup* SystemGroup  = new OptionGroup("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  //  LanczosManager Lanczos(false);
  LanczosManager Lanczos(true);
  

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 8);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-ky", "if ky is conserved, indicate twice the inital momentum for the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-ky", "number of ky value to evaluate if ky is conserved", -1);
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "if positive, fix the total number of particles", -1);
  (*SystemGroup) += new BooleanOption ('\n', "all-fixednbrparticles", "do the calculation for all the fixed particle number sector from 0 to --nbr-particles");
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new SingleStringOption ('\n', "confining-file", "file describing the confining potential");
  (*SystemGroup) += new SingleStringOption ('\n', "superconducting-file", "file describing the superconducting order paarameter");
  (*SystemGroup) += new SingleDoubleOption ('\n', "charging-energy", "factor in front of the charging energy (i.e 1/(2C))", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "average-nbrparticles", "average number of particles", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "force-negativeky", "manually force to compute the negative ky sectors if ky is conserved");
  (*SystemGroup) += new SingleStringOption ('\n', "directory", "use a specific directory for the input data instead of the current one");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "aspect-ratio", "aspect ratio of the cylinder", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "cylinder-perimeter", "if non zero, fix the cylinder perimeter (in magnetic length unit) instead of the aspect ratio", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "use-bosons", "use bosons instead of fermions");

  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*SystemGroup) += new BooleanOption  ('\n', "compute-sparsity", "compute sparsity of Hamiltonian matrix");
  
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleStringOption  ('\n', "export-binhamiltonian", "export the hamiltonian as a binary file");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int LzMax = Manager.GetInteger("lzmax");
  int TotalSz = Manager.GetInteger("total-sz");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  int InitialLz = Manager.GetInteger("initial-ky");
  int NbrLz = Manager.GetInteger("nbr-ky");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool FirstRun = true;
  bool Statistics = !Manager.GetBoolean("use-bosons");
  
  int KValue = 1;
  int RValue = 2;

  double Ratio = Manager.GetDouble("aspect-ratio");
  double Perimeter = Manager.GetDouble("cylinder-perimeter");
  if (Perimeter > 0.0)
    {
      Ratio = (Perimeter * Perimeter) / (2.0 * M_PI * (LzMax + 1));
    }  

  char* FilePrefix = new char[512];

  if (Perimeter > 0.0)	
    {
      if (Statistics == true)
	{
	  sprintf (FilePrefix, "fermions_cylinder_perimeter_%.6f", Perimeter);
	}
      else
	{
	  sprintf (FilePrefix, "bosons_cylinder_perimeter_%.6f", Perimeter);
	}
    }
  else
    {
      if (Statistics == true)
	{
	  sprintf (FilePrefix, "fermions_cylinder_ratio_%.6f", Ratio);
	}
      else
	{
	  sprintf (FilePrefix, "bosons_cylinder_ratio_%.6f", Ratio);
	}      
    }

  double* OneBodyPotentialUpUp = 0;
  double* OneBodyPotentialDownDown = 0;
  double* OneBodyPotentialUpDown = 0;
  double* OneBodyPotentialPairing = 0;
  Complex** OffDiagonalOneBodyPotentialUpUp = 0;
  Complex** OffDiagonalOneBodyPotentialDownDown = 0;
  Complex* ComplexOneBodyPotentialPairing = 0;
  Complex** OffDiagonalComplexOneBodyPotentialPairing = 0;
  double** PseudoPotentials  = 0;
  bool PreserveKySymmetryFlag = false;
  int MaximumMomentumTransfer = 0;
 
  if ((Manager.GetString("interaction-file") == 0) && (Manager.GetString("confining-file") == 0))
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }

  if (Manager.GetString("interaction-file") != 0) 
    {
      if (FQHESphereSU2GetPseudopotentialsWithPairing(Manager.GetString("interaction-file"), LzMax, PseudoPotentials, 
						      OneBodyPotentialUpUp, OneBodyPotentialDownDown, OneBodyPotentialUpDown, OneBodyPotentialPairing) == false)
	{
	  return -1;
	}
      if (OneBodyPotentialUpDown != 0)
	{
	  cout << "warning, OneBodyPotentialUpDown is not supported" << endl;
	}
      if (PseudoPotentials != 0)
	{
	  cout << "warning, there should be no two-body pseudopotentials. Two-body interactions are implemented through the restriction to the quasihole Hilbert space" << endl;
	}
      if (OneBodyPotentialPairing != 0)
	{
	  ComplexOneBodyPotentialPairing = new Complex[LzMax + 1];
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      ComplexOneBodyPotentialPairing[i] = OneBodyPotentialPairing[i];
	    }
	}
    }
  else
    {
      MultiColumnASCIIFile ConfiningFile;
      if (ConfiningFile.Parse(Manager.GetString("confining-file")) == false)
	{
	  ConfiningFile.DumpErrors(cout);
	  return -1;
	}
      if ((ConfiningFile.GetNbrColumns() < 4) || (ConfiningFile.GetNbrLines() == 0))
	{
	  cout << "error, " << Manager.GetString("confining-file") << " has an invalid format" << endl;
	  return -1;
	}
      int* CreationIndices = ConfiningFile.GetAsIntegerArray(0);
      int* AnnihilationIndices = ConfiningFile.GetAsIntegerArray(1);
      double* TmpUpUpPotential = ConfiningFile.GetAsDoubleArray(2);
      double* TmpDownDownPotential = ConfiningFile.GetAsDoubleArray(3);
      double* TmpUpUpPhasePotential;
      if (ConfiningFile.GetNbrColumns() > 4)
	{
	  TmpUpUpPhasePotential = ConfiningFile.GetAsDoubleArray(4);
	}
      else
	{
	  TmpUpUpPhasePotential = new double[ConfiningFile.GetNbrLines()];
	  for (int i = 0; i < ConfiningFile.GetNbrLines(); ++i)
	    {
	      TmpUpUpPhasePotential[i] = 0.0;
	    }
	}
      double* TmpDownDownPhasePotential;
      if (ConfiningFile.GetNbrColumns() > 5)
	{
	  TmpDownDownPhasePotential = ConfiningFile.GetAsDoubleArray(5);
	}
      else
	{
	  TmpDownDownPhasePotential = new double[ConfiningFile.GetNbrLines()];
	  for (int i = 0; i < ConfiningFile.GetNbrLines(); ++i)
	    {
	      TmpDownDownPhasePotential[i] = 0.0;
	    }
	}
      for (int i = 0; i < ConfiningFile.GetNbrLines(); ++i)
	{
	  if ((CreationIndices[i] < 0) || (AnnihilationIndices[i] < 0) || (CreationIndices[i] > LzMax) || (AnnihilationIndices[i] > LzMax))
	    {
	      cout << "invalid indices line " << (i+1) << endl;
	      return -1;
	    }
	  if (abs(CreationIndices[i] - AnnihilationIndices[i]) > MaximumMomentumTransfer)
	    {
	      MaximumMomentumTransfer = abs(CreationIndices[i] - AnnihilationIndices[i]);
	    }
	}
      OneBodyPotentialUpUp = new double[LzMax + 1];
      OneBodyPotentialDownDown = new double[LzMax + 1];
      for (int i = 0; i <= LzMax; ++i)
	{
	  OneBodyPotentialUpUp[i] = 0.0;
	  OneBodyPotentialDownDown[i] = 0.0;	  
	}
      if (MaximumMomentumTransfer > 0)
	{
	  OffDiagonalOneBodyPotentialUpUp = new Complex*[LzMax + 1];
	  OffDiagonalOneBodyPotentialDownDown = new Complex*[LzMax + 1];
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      OffDiagonalOneBodyPotentialUpUp[i] = new Complex[MaximumMomentumTransfer];
	      OffDiagonalOneBodyPotentialDownDown[i] = new Complex[MaximumMomentumTransfer];
	      for (int j = 0; j < MaximumMomentumTransfer; ++j)	  
		{
		  OffDiagonalOneBodyPotentialUpUp[i][j] = 0.0;
		  OffDiagonalOneBodyPotentialDownDown[i][j] = 0.0;
		}
	    }	  
	}
      for (int i = 0; i < ConfiningFile.GetNbrLines(); ++i)
	{
	  int TmpMomentumTransfer = CreationIndices[i] - AnnihilationIndices[i];
	  if (TmpMomentumTransfer == 0)
	    {
	      OneBodyPotentialUpUp[CreationIndices[i]] = TmpUpUpPotential[i];
	      OneBodyPotentialDownDown[CreationIndices[i]] = TmpDownDownPotential[i];
	    }
	  else
	    {
	      if (TmpMomentumTransfer > 0)
		{
		  OffDiagonalOneBodyPotentialUpUp[CreationIndices[i]][TmpMomentumTransfer - 1] = (TmpUpUpPotential[i] * Phase(M_PI * TmpUpUpPhasePotential[i]));
		  OffDiagonalOneBodyPotentialDownDown[CreationIndices[i]][TmpMomentumTransfer - 1] = (TmpDownDownPotential[i] * Phase(M_PI * TmpDownDownPhasePotential[i]));
		}
	      else
		{
		  OffDiagonalOneBodyPotentialUpUp[CreationIndices[i]][-TmpMomentumTransfer - 1] = (TmpUpUpPotential[i] * Phase(-M_PI * TmpUpUpPhasePotential[i]));
		  OffDiagonalOneBodyPotentialDownDown[CreationIndices[i]][-TmpMomentumTransfer - 1] = (TmpDownDownPotential[i] * Phase(-M_PI * TmpDownDownPhasePotential[i]));
		}
	    }
	}
      if (Manager.GetString("superconducting-file") != 0)
	{
	  MultiColumnASCIIFile SuperconductingFile;
	  if (SuperconductingFile.Parse(Manager.GetString("superconducting-file")) == false)
	    {
	      SuperconductingFile.DumpErrors(cout);
	      return -1;
	    }
	  if ((SuperconductingFile.GetNbrColumns() < 3) || (SuperconductingFile.GetNbrLines() == 0))
	    {
	      cout << "error, " << Manager.GetString("confining-file") << " has an invalid format" << endl;
	      return -1;
	    }
	  int* CreationIndices = SuperconductingFile.GetAsIntegerArray(0);
	  int* AnnihilationIndices = SuperconductingFile.GetAsIntegerArray(1);
	  double* TmpPairingPotential = SuperconductingFile.GetAsDoubleArray(2);
	  double* TmpPairingPhasePotential;
	  if (SuperconductingFile.GetNbrColumns() > 3)
	    {
	      TmpPairingPhasePotential = SuperconductingFile.GetAsDoubleArray(3);
	    }
	  else
	    {
	      TmpPairingPhasePotential = new double[SuperconductingFile.GetNbrLines()];
	      for (int i = 0; i < SuperconductingFile.GetNbrLines(); ++i)
		{
		  TmpPairingPhasePotential[i] = 0.0;
		}
	    }
	  for (int i = 0; i < SuperconductingFile.GetNbrLines(); ++i)
	    {
	      if ((CreationIndices[i] < 0) || (AnnihilationIndices[i] < 0) || (CreationIndices[i] > LzMax) || (AnnihilationIndices[i] > LzMax))
		{
		  cout << "invalid indices line " << (i+1) << endl;
		  return -1;
		}
	      if (abs(CreationIndices[i] - AnnihilationIndices[i]) > MaximumMomentumTransfer)
		{
		  MaximumMomentumTransfer = abs(CreationIndices[i] - AnnihilationIndices[i]);
		}
	    }
	  ComplexOneBodyPotentialPairing = new Complex[LzMax + 1];
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      ComplexOneBodyPotentialPairing[i] = 0.0;
	    }
	  if (MaximumMomentumTransfer > 0)
	    {
	      OffDiagonalComplexOneBodyPotentialPairing = new Complex*[LzMax + 1];
	      for (int i = 0; i <= LzMax; ++i)
		{
		  OffDiagonalComplexOneBodyPotentialPairing[i] = new Complex[(2 * MaximumMomentumTransfer) + 1];
		  for (int j = 0; j < MaximumMomentumTransfer; ++j)	  
		    {
		      OffDiagonalComplexOneBodyPotentialPairing[i][j] = 0.0;
		    }
		}	  
	    }
	  for (int i = 0; i < SuperconductingFile.GetNbrLines(); ++i)
	    {
	      int TmpMomentumTransfer = CreationIndices[i] - AnnihilationIndices[i];
	      if (TmpMomentumTransfer == 0)
		{
		  ComplexOneBodyPotentialPairing[CreationIndices[i]] = TmpPairingPotential[i] * Phase(M_PI * TmpPairingPhasePotential[i]);
		}
	      else
		{
		  OffDiagonalComplexOneBodyPotentialPairing[CreationIndices[i]][MaximumMomentumTransfer + TmpMomentumTransfer] = TmpPairingPotential[i] * Phase(M_PI * TmpPairingPhasePotential[i]);
		}
	    }	  
	}
    }
  if (MaximumMomentumTransfer == 0)
    {
      PreserveKySymmetryFlag = true;
    }



  char* OutputNameLz = new char [256 + strlen(Manager.GetString("interaction-name"))];
  int TmpFixedNbrParticles = 0;
  if (Manager.GetInteger("nbr-particles") >= 0)
    TmpFixedNbrParticles = Manager.GetInteger("nbr-particles");
  if (PreserveKySymmetryFlag == true)
    {
      sprintf (OutputNameLz, "%s_su2_quasiholes_%s_cenergy_%.6f_n0_%.6f_pairing_n_%d_2s_%d_sz_%d_lz.dat", FilePrefix, Manager.GetString("interaction-name"), 
	       Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), TmpFixedNbrParticles, LzMax, TotalSz);
    }
  else
    {
      sprintf (OutputNameLz, "%s_su2_quasiholes_%s_cenergy_%.6f_n0_%.6f_pairing_n_%d_2s_%d_sz_%d.dat", FilePrefix, Manager.GetString("interaction-name"), 
	       Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), TmpFixedNbrParticles, LzMax, TotalSz);
    }

  int MinNbrParticles = abs(TotalSz);
  int MaxNbrParticles = (2 * (LzMax + 1)) - abs(TotalSz);
  if (Manager.GetInteger("nbr-particles") >= 0)
    {
      MinNbrParticles = Manager.GetInteger("nbr-particles");
      MaxNbrParticles = MinNbrParticles;
      if (OneBodyPotentialPairing != 0)
	{
	  cout << "discarding superconducting coupling" << endl;
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      OneBodyPotentialPairing[i] = 0.0;
	    }
	}
    }


  if (PreserveKySymmetryFlag == true)
    {
      long TotalDimension = 0l;
      int MaxL = 0;
      for (int TmpNbrParticles = MinNbrParticles; TmpNbrParticles <= MaxNbrParticles; TmpNbrParticles += 2)
	{
	  int NbrUp = (TmpNbrParticles + TotalSz) / 2;
	  int NbrDown = (TmpNbrParticles - TotalSz) / 2;
	  if ((NbrUp <= (LzMax + 1)) && (NbrUp >= 0) && (NbrDown <= (LzMax + 1)) && (NbrDown >= 0))
	    {
	      // warning, this should not work for k > 1
	      if (KValue > 1)
		cout << "please fix your code for k>1" << endl;
	      int MaxTotalLzUp = (LzMax * NbrUp) - ((KValue + RValue) * ((NbrUp - 1) * NbrUp) / KValue);
	      int MaxTotalLzDown = (LzMax * NbrDown) - ((KValue + RValue) * ((NbrDown - 1) * NbrDown) / KValue);
	      if ((MaxTotalLzUp + MaxTotalLzDown) > MaxL)
		{
		  MaxL = (MaxTotalLzUp + MaxTotalLzDown);
		}
	    }
	}

      int  L = 0;
      L = InitialLz;
      if ((abs(MaxL) & 1) != 0)
	L |= 1;
      else
	L &= ~0x1;
      
      if (NbrLz > 0)
	{
	  if (L + (2 * (NbrLz - 1)) < MaxL)
	    MaxL = L + (2 * (NbrLz - 1));
	}
      
      if (Manager.GetBoolean("force-negativeky"))
	L = -MaxL;
      
      
      for (; L <= MaxL; L += 2)
	{
	  int LocalNbrParticles = -1;
	  int MaxLocalNbrParticles = -1;
	  if (Manager.GetInteger("nbr-particles") >= 0)
	    {
	      LocalNbrParticles = Manager.GetInteger("nbr-particles");
	      MaxLocalNbrParticles = LocalNbrParticles;
	      if (Manager.GetBoolean("all-fixednbrparticles") == true)
		{
		  if ((TotalSz & 1) == 0)
		    LocalNbrParticles = 0;
		  else
		    LocalNbrParticles = 1;
		}
	    }
	  for (; LocalNbrParticles <= MaxLocalNbrParticles; LocalNbrParticles += 2)
	    {
	      QuasiholeOnSphereWithSpinAndPairing* Space = 0;
	      if (LocalNbrParticles >= 0)
		{
		  Space = new QuasiholeOnSphereWithSpin (KValue, RValue, L, LzMax, LocalNbrParticles, TotalSz, Manager.GetString("directory"), FilePrefix);
		}
	      else
		{
		  Space = new QuasiholeOnSphereWithSpinAndPairing (KValue, RValue, L, LzMax, TotalSz, Manager.GetString("directory"), FilePrefix);
		}
	      if (Space->GetLargeHilbertSpaceDimension() > 0l)
		{
		  TotalDimension += Space->GetLargeHilbertSpaceDimension();
		  cout << "Ky=" << L;
		  if (Manager.GetInteger("nbr-particles") >= 0)
		    cout << ", N=" << LocalNbrParticles;
		  cout << endl;
		  
		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
		  AbstractHamiltonian* Hamiltonian = 0;
		  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		    Memory = Architecture.GetArchitecture()->GetLocalMemory();
		  // 	      Hamiltonian = new ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing (Space, LzMax, OneBodyPotentialUpUp, 
		  // 													     OneBodyPotentialDownDown,
		  // 													     OneBodyPotentialPairing,
		  // 													     Manager.GetDouble("charging-energy"), 
		  // 													     Manager.GetDouble("average-nbrparticles"),
		  // 													     Architecture.GetArchitecture(), 
		  // 													     Memory);
		  Hamiltonian = new ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta (Space, LzMax, MaximumMomentumTransfer, 
															   OneBodyPotentialUpUp, OneBodyPotentialDownDown,
															   0, 0,
															   ComplexOneBodyPotentialPairing, 0,
															   Manager.GetDouble("charging-energy"), 
															   Manager.GetDouble("average-nbrparticles"),
															   Architecture.GetArchitecture(), 
															   Memory);
	      
		  double Shift = 0.0;
		  Hamiltonian->ShiftHamiltonian(Shift);
		  char* EigenvectorName = 0;
		  if ((Manager.GetBoolean("eigenstate") == true) || (Manager.GetBoolean("all-eigenstates") == true))
		    {
		      EigenvectorName = new char [256 + strlen(FilePrefix) + strlen(Manager.GetString("interaction-name"))];
		      if (LocalNbrParticles >= 0)
			{
			  sprintf (EigenvectorName, "%s_su2_quasiholes_%s_cenergy_%.6f_n0_%.6f_pairing_n_%d_2s_%d_sz_%d_lz_%d", 
				   FilePrefix, Manager.GetString("interaction-name"), 
				   Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), LocalNbrParticles, LzMax, TotalSz, L);
			}
		      else
			{
			  sprintf (EigenvectorName, "%s_su2_quasiholes_%s_cenergy_%.6f_n0_%.6f_pairing_n_0_2s_%d_sz_%d_lz_%d", 
				   FilePrefix, Manager.GetString("interaction-name"), 
				   Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), LzMax, TotalSz, L);
			}
		    }
		  
		  char* ContentPrefix = new char[256];
		  if (Manager.GetBoolean("all-fixednbrparticles") == true)
		    {
		      sprintf (ContentPrefix, "%d %d", LocalNbrParticles, L);
		    }
		  else
		    {
		      sprintf (ContentPrefix, "%d", L);
		    }
		  char* SubspaceLegend = new char[256];
		  if (Manager.GetBoolean("all-fixednbrparticles") == true)
		    {
		      sprintf (SubspaceLegend, "N ky");
		    }
		  else
		    {
		      sprintf (SubspaceLegend, "ky");
		    }
		  
		  GenericComplexMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, ContentPrefix, SubspaceLegend, Shift, OutputNameLz, FirstRun, EigenvectorName);
		  //	      GenericRealMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, ContentPrefix, SubspaceLegend, Shift, OutputNameLz, FirstRun, EigenvectorName);
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  if (EigenvectorName != 0)
		    {
		      delete[] EigenvectorName;
		    }
		  delete Hamiltonian;
		  if (FirstRun == true)
		    FirstRun = false;
		}
	      delete Space;
	    }
	}
      cout << "total Hilbert space dimension =" << TotalDimension << endl; 
    }
  else
    {
      int LocalNbrParticles = -1;
      int MaxLocalNbrParticles = -1;
      if (Manager.GetInteger("nbr-particles") >= 0)
	{
	  LocalNbrParticles = Manager.GetInteger("nbr-particles");
	  MaxLocalNbrParticles = LocalNbrParticles;
	  if (Manager.GetBoolean("all-fixednbrparticles") == true)
	    {
	      if ((TotalSz & 1) == 0)
		LocalNbrParticles = 0;
	      else
		LocalNbrParticles = 1;
	    }
	}
      for (; LocalNbrParticles <= MaxLocalNbrParticles; LocalNbrParticles += 2)
	{
	  QuasiholeOnSphereWithSpinAndPairing* Space = 0;
	  if (LocalNbrParticles >= 0)
	    {
	      //	      Space = new QuasiholeOnSphereWithSpinAllMomenta (KValue, RValue, LzMax, LocalNbrParticles, TotalSz, MaximumMomentumTransfer, 
	      // Manager.GetString("directory"), FilePrefix);
	    }
	  else
	    {
	      Space = new QuasiholeOnSphereWithSpinAndPairingAllMomenta (KValue, RValue, LzMax, TotalSz, MaximumMomentumTransfer, Manager.GetString("directory"), FilePrefix);
	    }
	  if (Space->GetLargeHilbertSpaceDimension() > 0l)
	    {
	      if (Manager.GetInteger("nbr-particles") >= 0)
		cout << ", N=" << LocalNbrParticles;
	      cout << endl;
	      
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	      AbstractHamiltonian* Hamiltonian = 0;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      cout << "MaximumMomentumTransfer = " << MaximumMomentumTransfer << " "  << OffDiagonalOneBodyPotentialUpUp << endl;
	      Hamiltonian = new ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta (Space, LzMax, MaximumMomentumTransfer, 
														       OneBodyPotentialUpUp, OneBodyPotentialDownDown,
														       OffDiagonalOneBodyPotentialUpUp, OffDiagonalOneBodyPotentialDownDown,
														       ComplexOneBodyPotentialPairing, OffDiagonalComplexOneBodyPotentialPairing,
														       Manager.GetDouble("charging-energy"), 
														       Manager.GetDouble("average-nbrparticles"),
														       Architecture.GetArchitecture(), 
														       Memory);
	      
	      double Shift = 0.0;
	      Hamiltonian->ShiftHamiltonian(Shift);
	      char* EigenvectorName = 0;
	      if ((Manager.GetBoolean("eigenstate") == true) || (Manager.GetBoolean("all-eigenstates") == true))
		{
		  EigenvectorName = new char [256 + strlen(FilePrefix) + strlen(Manager.GetString("interaction-name"))];
		  if (LocalNbrParticles >= 0)
		    {
		      sprintf (EigenvectorName, "%s_su2_quasiholes_%s_cenergy_%.6f_n0_%.6f_pairing_n_%d_2s_%d_sz_%d", 
			       FilePrefix, Manager.GetString("interaction-name"), 
			       Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), LocalNbrParticles, LzMax, TotalSz);
		    }
		  else
		    {
		      sprintf (EigenvectorName, "%s_su2_quasiholes_%s_cenergy_%.6f_n0_%.6f_pairing_n_0_2s_%d_sz_%d", 
			       FilePrefix, Manager.GetString("interaction-name"), 
			       Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), LzMax, TotalSz);
		    }
		}
	      
	      char* ContentPrefix = new char[256];
	      if (Manager.GetBoolean("all-fixednbrparticles") == true)
		{
		  sprintf (ContentPrefix, "%d", LocalNbrParticles);
		}
	      else
		{
		  sprintf (ContentPrefix, "");
		}
	      char* SubspaceLegend = new char[256];
	      if (Manager.GetBoolean("all-fixednbrparticles") == true)
		{
		  sprintf (SubspaceLegend, "N");
		}
	      else
		{
		  sprintf (SubspaceLegend, "");
		}
	      
	      GenericComplexMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, ContentPrefix, SubspaceLegend, Shift, OutputNameLz, FirstRun, EigenvectorName);
	      //	      GenericRealMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, ContentPrefix, SubspaceLegend, Shift, OutputNameLz, FirstRun, EigenvectorName);
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      if (EigenvectorName != 0)
		{
		  delete[] EigenvectorName;
		}
	      delete Hamiltonian;
	      if (FirstRun == true)
		FirstRun = false;
	    }
	  delete Space;
	}
    }
  return 0;
}
