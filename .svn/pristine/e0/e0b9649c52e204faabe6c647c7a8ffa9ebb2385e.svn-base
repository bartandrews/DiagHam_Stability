#include "Options/Options.h"

#include "Tools/FTITightBinding/TightBindingModelTwoOrbitalSquareLattice.h"
#include "Tools/FTITightBinding/TightBindingModelCylinderTwoOrbitalSquareLattice.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


// extract the correlation matrix for a given region out of the correlation matrix for a bigger region
//
// correlationMatrix = correlation matrix for the bigger region
// sourceNbrSitesX = number of sites along x for the bigger region
// sourceNbrSitesY = number of sites along y for the bigger region
// targetNbrSitesX = number of sites along x for the smaller region
// targetNbrSitesY = number of sites along y for the smaller region
// nbrBands = number of orbitals per unit cell
HermitianMatrix TIEntanglementEntropyExtractCorrelationMatrix(HermitianMatrix& correlationMatrix, int sourceNbrSitesX, int sourceNbrSitesY, 
							      int targetNbrSitesX, int targetNbrSitesY, int nbrBands);


int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("TIEntanglementEntropyRealSpacePartition" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of unit cells along the x direction", 4);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of unit cells along the y direction", 4);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "imag part of inter-orbital hopping between nearest neighbors along the x direction", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "inter-orbital hopping between nearest neighbors along the y direction", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t3", "intra-orbital hopping between nearest neighbors", 1.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "folding", "folding factor for the momenta along sigma_x and sigma_y", 1);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "use-cylinder", "use a cylinder geometry instead of a torus geometry (open boundary conditions along the x axis)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrsitex-a", "number of unit cells along the x direction for the part A", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrsitey-a", "number of unit cells along the y direction for the part A", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-nbrsitexa", "maximum number of unit cells along the x direction for the part A (equal to --nbrsitex-a if negative)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-nbrsiteya", "maximum number of unit cells along the y direction for the part A (equal to --nbrsitey-a if negative)", -1);
  (*SystemGroup) += new SingleStringOption  ('\n', "cuts", "provide the description of all cuts as a two column formatted text file");
  (*SystemGroup) += new BooleanOption ('\n', "use-random", "generate a random non-vacuum state");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "run-id", "add an additional run id to the file name when using the --use-random option", 0);  
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");  
  (*SystemGroup) += new BooleanOption ('\n', "use-approximation", "use a saddle appoximation to evaluate the entanglement entropy");
  (*SystemGroup) += new BooleanOption ('\n', "use-rational", "use rational number to overcome accuracy issues");
  (*SystemGroup) += new BooleanOption ('\n', "disable-ytranslation", "do not use the translation along y even if the cut preserves this symmetry");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the Chern number of the fully filled band (only available in singleparticle-spectrum mode)");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type TIEntanglementEntropyRealSpacePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  bool ShowTimeFlag = Manager.GetBoolean("show-time");

  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey"); 
  int NbrSitesXA = Manager.GetInteger("nbrsitex-a"); 
  int NbrSitesYA = Manager.GetInteger("nbrsitey-a"); 
  int NbrSitesA = NbrSitesXA * NbrSitesYA;
  int NbrCuts = 1;
  int* CutX = 0;
  int* CutY = 0;
  bool PreserveYTranslation = false;
  bool CylinderFlag = Manager.GetBoolean("use-cylinder");

  if (Manager.GetString("cuts") == 0)
    {
      if ((Manager.GetInteger("max-nbrsitexa") > NbrSitesXA) || (Manager.GetInteger("max-nbrsiteya") > NbrSitesYA))
	{
	  NbrCuts = Manager.GetInteger("max-nbrsitexa") - NbrSitesXA + 1;
	  CutX = new int [NbrCuts];
	  CutY = new int [NbrCuts];
	  for (int i = 0; i < NbrCuts; ++i)
	    {
	      CutX[i] = NbrSitesXA + i;
	      CutY[i] = NbrSitesXA + i;	      
	    }	  
	}
      else
	{
	  CutX = new int [NbrCuts];
	  CutY = new int [NbrCuts];
	  CutX[0] = NbrSitesXA;
	  CutY[0] = NbrSitesYA;
	  if ((NbrSitesYA == NbrSitesY) && (Manager.GetBoolean("disable-ytranslation") == false))
	    {
	      PreserveYTranslation = true;
	    }
	}
    }
  else
    {
      MultiColumnASCIIFile CutFile;
      if (CutFile.Parse(Manager.GetString("cuts")) == false)
	{
	  CutFile.DumpErrors(cout);
	  return -1;
	}
      NbrCuts = CutFile.GetNbrLines();
      CutX = CutFile.GetAsIntegerArray(0);
      CutY = CutFile.GetAsIntegerArray(1);
    }
  int MaxNbrSitesXA = CutX[NbrCuts - 1];
  int MaxNbrSitesYA = CutY[NbrCuts - 1];
  int MaxNbrSitesA = MaxNbrSitesXA * MaxNbrSitesYA;

  char* FilePrefix = new char [512];
  if (CylinderFlag == false)
    {
      sprintf(FilePrefix, "fermions_twoorbitals_n_%d_x_%d_y_%d_t1_%f_t2_%f_t3_%f_mus_%f_gx_%f_gy_%f", (NbrSitesX * NbrSitesY), NbrSitesX, NbrSitesY,
	      Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), Manager.GetDouble("mu-s"), 
	      Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
    }
  else
    {
      sprintf(FilePrefix, "fermions_cylinder_twoorbitals_n_%d_x_%d_y_%d_t1_%f_t2_%f_t3_%f_mus_%f_gx_%f_gy_%f", (NbrSitesX * NbrSitesY), NbrSitesX, NbrSitesY,
	      Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), Manager.GetDouble("mu-s"), 
	      Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
    }

  Abstract1DTightBindingModel* TightBindingModel;
  if (CylinderFlag == false)
    {
      TightBindingModel = new TightBindingModelTwoOrbitalSquareLattice (NbrSitesX, NbrSitesY, 
									Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), 
									Manager.GetInteger("folding"), Manager.GetDouble("mu-s"), 
									Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), 
									true);
      if (Manager.GetBoolean("singleparticle-chernnumber") == true)
	cout << "Chern number = " << ((Abstract2DTightBindingModel*) TightBindingModel)->ComputeChernNumber(0) << endl;
    }
  else
    {
      TightBindingModel = new TightBindingModelCylinderTwoOrbitalSquareLattice (NbrSitesY, NbrSitesX, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), 
										Manager.GetInteger("folding"), Manager.GetDouble("mu-s"), 
										Manager.GetDouble("gamma-x"), false,  Architecture.GetArchitecture(), true);
    }      
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      char* EigenvalueOutputFile = new char [strlen(FilePrefix) + 64];
      sprintf (EigenvalueOutputFile, "%s.dat", FilePrefix);
      TightBindingModel->WriteAsciiSpectrum(EigenvalueOutputFile);
    }
  double* TightBindingModelEnergies = 0;
  int* TightBindingModelLinearizedMomenta = 0;
  int* TightBindingModelBandIndices = 0;
  if (CylinderFlag == false)
    {
      TightBindingModel->GetEnergies(TightBindingModelEnergies, TightBindingModelLinearizedMomenta, 0);
      TightBindingModelBandIndices = new int [TightBindingModel->GetNbrStatePerBand()];
      for (int i = 0; i < TightBindingModel->GetNbrStatePerBand(); ++i)
	TightBindingModelBandIndices[i] = 0;
    }
  else
    {
      TightBindingModel->GetAllEnergies(TightBindingModelEnergies, TightBindingModelLinearizedMomenta, TightBindingModelBandIndices);
    }

  int* VacuumOneBodyLinearizedMomenta = 0; 
  int* VacuumOneBodyLinearizedBandIndices = 0;
  double VacuumTotalEnergy = 0.0;
  int VacuumXMomentum = 0;
  int VacuumYMomentum = 0;
  int TmpMomentumX;
  int TmpMomentumY;
  int VacuumNbrParticles = (TightBindingModel->GetNbrStatePerBand() * TightBindingModel->GetNbrBands()) / 2; 
  VacuumOneBodyLinearizedMomenta = new int[VacuumNbrParticles];
  VacuumOneBodyLinearizedBandIndices = new int[VacuumNbrParticles];
  for (int i = 0; i < VacuumNbrParticles; ++i)
    {
      TightBindingModel->GetLinearizedMomentumIndex(TightBindingModelLinearizedMomenta[i], TmpMomentumX, TmpMomentumY);
      VacuumXMomentum += TmpMomentumX;
      VacuumYMomentum += TmpMomentumY;
      VacuumOneBodyLinearizedMomenta[i] = TightBindingModel->GetLinearizedMomentumIndex(TmpMomentumX, TmpMomentumY);
      VacuumOneBodyLinearizedBandIndices[i] = TightBindingModelBandIndices[i];
      VacuumTotalEnergy += TightBindingModel->GetEnergy(VacuumOneBodyLinearizedBandIndices[i], VacuumOneBodyLinearizedMomenta[i]);
    }
  
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  
  char* EntropyFileName = new char [512];
  sprintf(EntropyFileName, "%s_xa_%d_ya_%d.ent", FilePrefix, MaxNbrSitesXA, MaxNbrSitesYA);
  ofstream File;
  File.open(EntropyFileName, ios::binary | ios::out);
  File.precision(14);

  char* OneBodyEntropyFileName = new char [512];
  sprintf(OneBodyEntropyFileName, "%s_xa_%d_ya_%d.onebody.ent", FilePrefix, MaxNbrSitesXA, MaxNbrSitesYA);
  ofstream OneBodyFile;
  OneBodyFile.open(OneBodyEntropyFileName, ios::binary | ios::out);
  OneBodyFile.precision(14);


  if (PreserveYTranslation == false)
    {
      if (ShowTimeFlag == true)
	{
	  gettimeofday (&(TotalStartingTime), 0);
	}
      HermitianMatrix EntanglementHamiltonian;
      if (CylinderFlag == false)
	{
	  EntanglementHamiltonian = ((Abstract2DTightBindingModel*) TightBindingModel)->EvaluateFullTwoPointCorrelationFunction(MaxNbrSitesXA, MaxNbrSitesYA, VacuumOneBodyLinearizedMomenta, VacuumNbrParticles, 0);
	}
      else
	{
	  EntanglementHamiltonian = TightBindingModel->EvaluateFullTwoPointCorrelationFunction(MaxNbrSitesXA, MaxNbrSitesYA, VacuumOneBodyLinearizedMomenta, VacuumOneBodyLinearizedBandIndices, VacuumNbrParticles);
	}
      if (ShowTimeFlag == true)
	{
	  gettimeofday (&(TotalEndingTime), 0);
	  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	  cout << "correlation matrix evaluated in " << Dt << "s" << endl;
	}
      for (int CurrentCutIndex = 0; CurrentCutIndex < NbrCuts ; ++CurrentCutIndex)
	{
	  NbrSitesXA = CutX[CurrentCutIndex];
	  NbrSitesYA = CutY[CurrentCutIndex];
	  cout << "computing entropy of a " << NbrSitesXA << "x" << NbrSitesYA << " patch" << endl;
	  int TotalNbrSitesA = 0;
	  if (CylinderFlag == false)
	    {
	      TotalNbrSitesA = NbrSitesXA * NbrSitesYA * TightBindingModel->GetNbrBands();
	    }
	  else
	    {
	      TotalNbrSitesA = NbrSitesXA * NbrSitesYA * 2;
	    }
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalStartingTime), 0);
	    }
	  RealDiagonalMatrix VacuumOneBodyEntanglementEnergies(TotalNbrSitesA, true);
	  HermitianMatrix TmpEntanglementHamiltonian;
	  if (CylinderFlag == false)
	    {
	      TmpEntanglementHamiltonian = TIEntanglementEntropyExtractCorrelationMatrix(EntanglementHamiltonian, MaxNbrSitesXA, MaxNbrSitesYA, NbrSitesXA, NbrSitesYA, 
											 TightBindingModel->GetNbrBands());
	    }
	  else
	    {
	      TmpEntanglementHamiltonian = TIEntanglementEntropyExtractCorrelationMatrix(EntanglementHamiltonian, MaxNbrSitesXA, MaxNbrSitesYA, NbrSitesXA, NbrSitesYA, 2);
	    }
#ifdef __LAPACK__
	  TmpEntanglementHamiltonian.LapackDiagonalize(VacuumOneBodyEntanglementEnergies);
#else
	  TmpEntanglementHamiltonian.Diagonalize(VacuumOneBodyEntanglementEnergies);
#endif
	  VacuumOneBodyEntanglementEnergies.SortMatrixUpOrder();
	  int MinOneBodyEntanglementEnergyIndex = 0;
	  int MaxOneBodyEntanglementEnergyIndex = VacuumOneBodyEntanglementEnergies.GetNbrRow() - 1;  
	  while ((MinOneBodyEntanglementEnergyIndex <= MaxOneBodyEntanglementEnergyIndex) && (VacuumOneBodyEntanglementEnergies[MinOneBodyEntanglementEnergyIndex] <= MACHINE_PRECISION))
	    ++MinOneBodyEntanglementEnergyIndex;
	  while ((MinOneBodyEntanglementEnergyIndex <= MaxOneBodyEntanglementEnergyIndex) && (VacuumOneBodyEntanglementEnergies[MaxOneBodyEntanglementEnergyIndex] >= (1.0 - MACHINE_PRECISION)))
	    --MaxOneBodyEntanglementEnergyIndex;
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalEndingTime), 0);
	      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "diagonalization done in " << Dt << "s" << endl;
	    }
	  
	  double EntanglementEntropy = 0.0;
	  double NonVacuumEntanglementEntropy = 0.0;
	  double NbrParticleFluctuation = 0.0;
	  int NbrVacuumOneBodyEntanglementTrimmedEnergies = MaxOneBodyEntanglementEnergyIndex - MinOneBodyEntanglementEnergyIndex + 1;
	  int NbrRejectedOneBodyEntropies = VacuumOneBodyEntanglementEnergies.GetNbrRow() - NbrVacuumOneBodyEntanglementTrimmedEnergies;
	  double* VacuumOneBodyEntanglementTrimmedEnergies = new double[NbrVacuumOneBodyEntanglementTrimmedEnergies];
	  for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies; ++i)
	    {
	      VacuumOneBodyEntanglementTrimmedEnergies[i] = VacuumOneBodyEntanglementEnergies[i + MinOneBodyEntanglementEnergyIndex];	
	      OneBodyFile << NbrSitesXA << " " << NbrSitesYA << " " << VacuumOneBodyEntanglementTrimmedEnergies[i] << endl;
	    }
	  double SumAlphaFactor = 0.0;      
	  for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies; ++i)
	    {
	      EntanglementEntropy -= VacuumOneBodyEntanglementTrimmedEnergies[i] * log (VacuumOneBodyEntanglementTrimmedEnergies[i]);
	      EntanglementEntropy -= (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i]) * log (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i]);
	      NbrParticleFluctuation += VacuumOneBodyEntanglementTrimmedEnergies[i] * (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i]);
	    }
	  NonVacuumEntanglementEntropy = EntanglementEntropy;
	  SumAlphaFactor = 1.0;
	  cout << "Normalization = " << SumAlphaFactor << endl;
	  cout << "Entanglement entropy = " << EntanglementEntropy << endl;
	  cout << "Nbr Rejected one-body entanglement energies = " << NbrRejectedOneBodyEntropies << " / " << TotalNbrSitesA << endl;
	  cout << "Fluctuation of the number of particles = " << NbrParticleFluctuation << endl;
	  File << NbrSitesXA << " " << NbrSitesYA << " " << EntanglementEntropy << " " << NonVacuumEntanglementEntropy << " " << NbrParticleFluctuation << endl;
	  delete[] VacuumOneBodyEntanglementTrimmedEnergies;
	}
    }
  else
    {
      NbrSitesXA = CutX[0];
      NbrSitesYA = CutY[0];
      cout << "computing entropy of a " << NbrSitesXA << "x" << NbrSitesYA << " patch" << endl;
      int TotalNbrSitesA = 0;
      if (CylinderFlag == false)
	{
	  TotalNbrSitesA = NbrSitesXA * TightBindingModel->GetNbrBands();
	}
      else
	{
	  TotalNbrSitesA = NbrSitesXA * 2;
	}
      double EntanglementEntropy = 0.0;
      double NbrParticleFluctuation = 0.0;
      double** VacuumOneBodyEntanglementTrimmedEnergies = new double*[NbrSitesY];
      int* NbrVacuumOneBodyEntanglementTrimmedEnergies = new int[NbrSitesY];
      int TotalNbrRejectedOneValues = 0;
      for (int TmpKy = 0; TmpKy < NbrSitesY; ++TmpKy)
	{
	  RealDiagonalMatrix VacuumOneBodyEntanglementEnergies(TotalNbrSitesA, true);
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalStartingTime), 0);
	    }
	  HermitianMatrix TmpEntanglementHamiltonian;
	  TmpEntanglementHamiltonian = TightBindingModel->EvaluateFullMixedTwoPointCorrelationFunctionWithK(MaxNbrSitesXA, TmpKy, VacuumOneBodyLinearizedMomenta, 
													    VacuumOneBodyLinearizedBandIndices, VacuumNbrParticles);
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalEndingTime), 0);
	      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "correlation matrix  at ky = " << TmpKy << " evaluated in " << Dt << "s" << endl;
	    }
#ifdef __LAPACK__
	  TmpEntanglementHamiltonian.LapackDiagonalize(VacuumOneBodyEntanglementEnergies);
#else
	  TmpEntanglementHamiltonian.Diagonalize(VacuumOneBodyEntanglementEnergies);
#endif
	  VacuumOneBodyEntanglementEnergies.SortMatrixUpOrder();
	  int MinOneBodyEntanglementEnergyIndex = 0;
	  int MaxOneBodyEntanglementEnergyIndex = VacuumOneBodyEntanglementEnergies.GetNbrRow() - 1;  
	  while ((MinOneBodyEntanglementEnergyIndex <= MaxOneBodyEntanglementEnergyIndex) && (VacuumOneBodyEntanglementEnergies[MinOneBodyEntanglementEnergyIndex] <= MACHINE_PRECISION))
	    ++MinOneBodyEntanglementEnergyIndex;
	  while ((MinOneBodyEntanglementEnergyIndex <= MaxOneBodyEntanglementEnergyIndex) && (VacuumOneBodyEntanglementEnergies[MaxOneBodyEntanglementEnergyIndex] >= (1.0 - MACHINE_PRECISION)))
	    {
	      --MaxOneBodyEntanglementEnergyIndex;
	      ++TotalNbrRejectedOneValues;
	    }
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalEndingTime), 0);
	      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "diagonalization done in " << Dt << "s" << endl;
	    }
	  
	  NbrVacuumOneBodyEntanglementTrimmedEnergies[TmpKy] = MaxOneBodyEntanglementEnergyIndex - MinOneBodyEntanglementEnergyIndex + 1;
	  int NbrRejectedOneBodyEntropies = VacuumOneBodyEntanglementEnergies.GetNbrRow() - NbrVacuumOneBodyEntanglementTrimmedEnergies[TmpKy];
	  VacuumOneBodyEntanglementTrimmedEnergies[TmpKy] = new double[NbrVacuumOneBodyEntanglementTrimmedEnergies[TmpKy]];
	  for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies[TmpKy]; ++i)
	    {
	      VacuumOneBodyEntanglementTrimmedEnergies[TmpKy][i] = VacuumOneBodyEntanglementEnergies[i + MinOneBodyEntanglementEnergyIndex];	
	      OneBodyFile << NbrSitesXA << " " << NbrSitesYA << " " << TmpKy << " " << VacuumOneBodyEntanglementTrimmedEnergies[TmpKy][i] << endl;
	    }
	  for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies[TmpKy]; ++i)
	    {
	      EntanglementEntropy -= VacuumOneBodyEntanglementTrimmedEnergies[TmpKy][i] * log (VacuumOneBodyEntanglementTrimmedEnergies[TmpKy][i]);
	      EntanglementEntropy -= (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[TmpKy][i]) * log (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[TmpKy][i]);
	      NbrParticleFluctuation += VacuumOneBodyEntanglementTrimmedEnergies[TmpKy][i] * (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[TmpKy][i]);
	    }
	  cout << "Nbr Rejected one-body entanglement energies = " << NbrRejectedOneBodyEntropies << " / " << TotalNbrSitesA << endl;
	}
      File << NbrSitesXA << " " << NbrSitesYA << " " << EntanglementEntropy << " " << NbrParticleFluctuation << endl;
      cout << "Entanglement entropy = " << EntanglementEntropy << endl;
      cout << "Fluctuation of the number of particles = " << NbrParticleFluctuation << endl;
//      int TotalNbrVacuumOneBodyEntanglementTrimmedEnergies = 0;
//      for (int TmpKy = 0; TmpKy < NbrSitesY; ++TmpKy)
//        {
// 	 TotalNbrVacuumOneBodyEntanglementTrimmedEnergies += NbrVacuumOneBodyEntanglementTrimmedEnergies[TmpKy];
//        }
//      double* TotalVacuumOneBodyEntanglementTrimmedEnergies = new double[TotalNbrVacuumOneBodyEntanglementTrimmedEnergies];
//      int* TotalVacuumOneBodyEntanglementTrimmedMomenta  = new int[TotalNbrVacuumOneBodyEntanglementTrimmedEnergies];
//      TotalNbrVacuumOneBodyEntanglementTrimmedEnergies = 0;
//      for (int TmpKy = 0; TmpKy < NbrSitesY; ++TmpKy)
//        {
// 	 for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies[TmpKy]; ++i)
// 	   {
// 	     TotalVacuumOneBodyEntanglementTrimmedEnergies[TotalNbrVacuumOneBodyEntanglementTrimmedEnergies] = VacuumOneBodyEntanglementTrimmedEnergies[TmpKy][i];
// 	     TotalVacuumOneBodyEntanglementTrimmedMomenta[TotalNbrVacuumOneBodyEntanglementTrimmedEnergies] = TmpKy;
// 	     ++TotalNbrVacuumOneBodyEntanglementTrimmedEnergies;
// 	   }
//        }
//      SortArrayDownOrdering<int>(TotalVacuumOneBodyEntanglementTrimmedEnergies, TotalVacuumOneBodyEntanglementTrimmedMomenta, TotalNbrVacuumOneBodyEntanglementTrimmedEnergies);
//      int NaturalNbrParticles = NbrSitesXA * NbrSitesYA;
//      int InitialNbrParticlesA = NaturalNbrParticles - TotalNbrRejectedOneValues + 1;
//      int NaturalMomentumA = 0;
// //      for (int i = 0; i < InitialNbrParticlesA; ++i)
// //        NaturalMomentumA  += TotalVacuumOneBodyEntanglementTrimmedMomenta[i];
//      bool TmpZeroFlag = false;
//      for (int i = InitialNbrParticlesA; (i < TotalNbrVacuumOneBodyEntanglementTrimmedEnergies) && (TmpZeroFlag == false); ++i)
//        {
// 	 double TmpLargestEigenvalue = 1.0;
// 	 int TmpLargestEigenvalueMomentum = -NaturalMomentumA;
// 	 for (int j = InitialNbrParticlesA; j <= i; ++j)
// 	   {
// 	     TmpLargestEigenvalue *= TotalVacuumOneBodyEntanglementTrimmedEnergies[j];
// 	     TmpLargestEigenvalueMomentum += TotalVacuumOneBodyEntanglementTrimmedMomenta[j];
// 	   }
// 	 for (int j = InitialNbrParticlesA; j <= i; ++j)
// 	   {
// 	     TmpLargestEigenvalue /= (1.0 - TotalVacuumOneBodyEntanglementTrimmedEnergies[j]);
// 	   }	 
// 	 cout << i << " " << TmpLargestEigenvalue << " " << (TmpLargestEigenvalueMomentum % NbrSitesY) << endl;
// // 	 remove % NbrSitesY?
// 	 if (TmpLargestEigenvalue < 1e-50)
// 	   {
// 	     TmpZeroFlag = true;
// 	   }
//        }
//      InitialNbrParticlesA = NaturalNbrParticles - TotalNbrRejectedOneValues - 1;
//      TmpZeroFlag = false;
//      for (int i = InitialNbrParticlesA; (i >= 0 ) && (TmpZeroFlag == false); --i)
//        {
// 	 double TmpLargestEigenvalue = 1.0;
// 	 int TmpLargestEigenvalueMomentum = -NaturalMomentumA;
// 	 for (int j = i; j <= InitialNbrParticlesA; ++j)
// 	   {
// 	     TmpLargestEigenvalue /= TotalVacuumOneBodyEntanglementTrimmedEnergies[j];
// 	     TmpLargestEigenvalueMomentum -= TotalVacuumOneBodyEntanglementTrimmedMomenta[j];
// 	   }
// 	 for (int j = i; j <= InitialNbrParticlesA; ++j)
// 	   {
// 	     TmpLargestEigenvalue *= (1.0 - TotalVacuumOneBodyEntanglementTrimmedEnergies[j]);
// 	   }	 
// 	 cout << i << " " << TmpLargestEigenvalue << " " << (TmpLargestEigenvalueMomentum % NbrSitesY) << endl;
// 	 if (TmpLargestEigenvalue < 1e-50)
// 	   {
// 	     TmpZeroFlag = true;
// 	   }
//        }
    }
  File.close();
  OneBodyFile.close();

  return 0;
}

// extract the correlation matrix for a given region out of the correlation matrix for a bigger region
//
// correlationMatrix = correlation matrix for the bigger region
// sourceNbrSitesX = number of sites along x for the bigger region
// sourceNbrSitesY = number of sites along y for the bigger region
// targetNbrSitesX = number of sites along x for the smaller region
// targetNbrSitesY = number of sites along y for the smaller region
// nbrBands = number of orbitals per unit cell

HermitianMatrix TIEntanglementEntropyExtractCorrelationMatrix(HermitianMatrix& correlationMatrix, int sourceNbrSitesX, int sourceNbrSitesY, 
							      int targetNbrSitesX, int targetNbrSitesY, int nbrBands)
{
  int TargetTotalNbrSites = targetNbrSitesX * targetNbrSitesY * nbrBands;
  HermitianMatrix TmpMatrix (TargetTotalNbrSites, true);
  Complex Tmp;
  for (int i = 0; i < TargetTotalNbrSites; ++i)
    {
      int TmpTargetX1 = i / (targetNbrSitesY * nbrBands);
      int TmpTargetY1 = (i % (targetNbrSitesY * nbrBands)) / nbrBands;
      int TmpTargetOrbital1 = i % nbrBands;
      int TmpSourceIndex1 = ((TmpTargetX1 * sourceNbrSitesY) + TmpTargetY1) * nbrBands + TmpTargetOrbital1;
      for (int j = i; j < TargetTotalNbrSites; ++j)
	{
	  int TmpTargetX2 = j / (targetNbrSitesY * nbrBands);
	  int TmpTargetY2 = (j % (targetNbrSitesY * nbrBands)) / nbrBands;
	  int TmpTargetOrbital2 = j % nbrBands;
	  int TmpSourceIndex2 = ((TmpTargetX2 * sourceNbrSitesY) + TmpTargetY2) * nbrBands + TmpTargetOrbital2;
	  correlationMatrix.GetMatrixElement(TmpSourceIndex1, TmpSourceIndex2, Tmp);
	  TmpMatrix.SetMatrixElement(i, j, Tmp);
	}
    }
  return TmpMatrix;
}

