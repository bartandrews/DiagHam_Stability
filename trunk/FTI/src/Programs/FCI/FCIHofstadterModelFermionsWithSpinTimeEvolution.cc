#include "Options/Options.h"

#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU4SpinMomentumSpace.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong.h"
#include "HilbertSpace/BosonOnLatticeRealSpace.h"
#include "HilbertSpace/BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong.h"


#include "Hamiltonian/ParticleOnLatticeHofstadterSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeHofstadterSingleBandGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeTwoBandHofstadterHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeFourBandHofstadterHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"
#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

// try to guess system information from file name
//
// filename = vector file name
// tauI = reference to the value of tau
// tfI = 
// nbrstepI = number of steps in previous simulation
// stepI = reference to the number of steps already done
// return value = true if no error occured

bool FTITimeEvolutionFindSystemInfoFromVectorFileName (char* filename, double& tauI, double& tfI, int& nbrStepI,int& stepI, double& uFinalI);

// create file containing the necessary energies for time-evolution with Hilbert-space truncation
//
// statisticPrefixEnergy
// kx = momentum along x
// ky = momentum along y
// sz = total magnetization
// szParity = parity under spin-inversion
// truncationIndex = number of low-energy states that are conserved for the time-evolution

char* FTITimeEvolutionCreateTruncatedEnergyFile (char* statisticPrefixEnergy, double gammaX, double gammaY, int kx, int ky, int sz, int szParity, int truncationIndex, double uFinal, double tau, double finalTime, int nbrSteps, bool createFile);

int main(int argc, char** argv)
{
  OptionManager Manager ("FCIHofstadterModelTimeEvolution" , "0.01");
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

  (*SystemGroup) += new SingleStringOption('\n', "initial-state", "name of the file containing the initial vector upon which e^{-iHt} acts");  
  (*SystemGroup) += new BooleanOption  ('\n', "complex", "initial vector is a complex vector");
  (*SystemGroup) += new BooleanOption  ('\n', "compute-energy", "compute the energy of each time-evolved vector");
  
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-final", "final value of the onsite potential strength", -10.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tau", "time where the final value of U is reached", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "final-time", "time where the code stops", 10.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-steps", "number of points to evaluate", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "save-period", "save time-evolved state every n steps, save only final step if negative", 1);
  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "truncate-basis", "if positive, use a truncated instantaneous basis with this number of states. Vector names are guessed from the name of the initial state. A path to the energy spectra and to the vectors has to be provided", 0);
  (*SystemGroup) += new SingleStringOption ('\n', "energy-path", "name of the folder where energy spectra are stored for truncated time-evolution");
  (*SystemGroup) += new SingleStringOption ('\n', "state-path", "name of the folder where vectors are stored for truncated time-evolution");
  (*SystemGroup) += new BooleanOption  ('\n', "preprocess-energy", "create energy file for truncated time-evolution");

  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  
  (*SystemGroup) += new BooleanOption  ('\n', "landau-x", "Use Landau gauge along the x-axis within unit cell");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrParticles; 
  int NbrCellX; 
  int NbrCellY;
  int UnitCellX; 
  int UnitCellY;
  int FluxPerCell;
  int MinNbrSinglets;
  int Sz;
  int Kx,Ky;
  int SzSymmetry;
  bool StatisticFlag;
  bool GutzwillerFlag;
  bool TranslationFlag;
  bool UsingConstraintMinNbrSinglets = false;
  char Axis ='y';
  bool SzSymmetryFlag;
  if (Manager.GetBoolean("landau-x"))
    Axis ='x';
  
  int TruncationIndex = Manager.GetInteger("truncate-basis");
  
  bool ResumeFlag = false;
  double TauI;
  double TfI;
  int NbrStepI;
  int StepI;
  double UFinalI;
  
  int NbrSteps = Manager.GetInteger("nbr-steps");
  int SavePeriod = Manager.GetInteger("save-period");
  double UFinal = Manager.GetDouble("u-final");
  double Tau = Manager.GetDouble("tau");
  double FinalTime = Manager.GetDouble("final-time");
  double TimeStep = FinalTime/ ((double)NbrSteps);
  double UPotential;
  
  double GammaX = Manager.GetDouble("gamma-x");
  double GammaY = Manager.GetDouble("gamma-y");
  int NbrStepsTau = (int) (NbrSteps * Tau / FinalTime) + 1;    
  
  
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  if (Manager.GetString("initial-state") == 0)
    {
      cout << " Error: an initial state must be provided" << endl;
      return -1;
    }
  char* StateFileName = Manager.GetString("initial-state");
  if (IsFile(StateFileName) == false)
    {
      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }
  
  if (FTIHofstadterdModelWith2DTranslationFindSystemInfoFromVectorFileName(StateFileName, NbrParticles, Kx, Ky, FluxPerCell , NbrCellX,  NbrCellY, UnitCellX, UnitCellY,  StatisticFlag, GutzwillerFlag, TranslationFlag) == false )
    {
      return -1;
    }

  if (FTIHofstadterModelWithSzFindSystemInfoFromVectorFileName(StateFileName, Sz, SzSymmetry, MinNbrSinglets, SzSymmetryFlag,UsingConstraintMinNbrSinglets) == false)
    {
      return -1;
    }

  if (FTITimeEvolutionFindSystemInfoFromVectorFileName(StateFileName, TauI, TfI, NbrStepI, StepI, UFinalI))
  {
    ResumeFlag = true;
    cout << "tau = " << Tau << ", tf = " << TfI << ", NbrSteps = " << NbrStepI << ", initial step = " << StepI << ", U final = " << UFinalI << endl;
    UFinal = UFinalI;
    Tau = TauI;
    NbrSteps = NbrStepI;
    FinalTime = TfI;
    TimeStep = FinalTime/ ((double)NbrSteps);
  }
    
  char*** EigenstateFile = 0;
  double** Energies = 0;
  if (TruncationIndex > 0)
  {
    EigenstateFile = new char** [NbrStepsTau];
    for (int i = 0; i < NbrStepsTau; ++i)
    {
      EigenstateFile[i] = new char* [TruncationIndex];
      for (int j = 0; j < TruncationIndex; ++j)
	EigenstateFile[i][j] = new char[strlen(Manager.GetString("state-path")) + 512];
    }
    
    char* StatisticPrefix = new char [128];
    Energies = new double*[TruncationIndex];
    
    int lenPrefixState = 0;
    if (MinNbrSinglets <= 0 )
      lenPrefixState += sprintf (StatisticPrefix + lenPrefixState, "fermions_su2");
    else
      lenPrefixState += sprintf (StatisticPrefix + lenPrefixState, "fermions_su2_minnbrsinglet_%d", MinNbrSinglets);
    
    lenPrefixState +=  sprintf (StatisticPrefix + lenPrefixState, "_realspace_hofstadter_X_%d_Y_%d_q_%d_n_%d_x_%d_y_%d", UnitCellX, UnitCellY, FluxPerCell, NbrParticles, NbrCellX, NbrCellY);
    
    
    char* StatisticPrefixEnergy = new char[strlen(Manager.GetString("energy-path")) + strlen(StatisticPrefix) + 2];
    sprintf(StatisticPrefixEnergy, "%s/%s", Manager.GetString("energy-path"), StatisticPrefix);
    char* EnergyFileName = FTITimeEvolutionCreateTruncatedEnergyFile(StatisticPrefixEnergy, GammaX, GammaY, Kx, Ky, Sz, SzSymmetry, TruncationIndex, UFinal, Tau, FinalTime, NbrSteps, Manager.GetBoolean("preprocess-energy"));
    delete[] StatisticPrefixEnergy;
    if (EnergyFileName == 0)
    {
      cout << "Energies could not all be retrieved for truncated time-evolution" << endl;
      return -1;
    }
    
    for(int i = 0 ; i < NbrStepsTau; i++)
    {
      double UPotential =  UFinal/ Tau * TimeStep * i;
      if (i == 0)
	UPotential = 0.0;
      for (int j = 0; j < TruncationIndex; ++j)
      {
	sprintf (EigenstateFile[i][j],"%s/%s_u_%g_gx_%g_gy_%g_sz_%d_szsym_%d_kx_%d_ky_%d.%d.vec", Manager.GetString("state-path"), StatisticPrefix, UPotential, GammaX, GammaY, Sz, SzSymmetry, Kx, Ky, j);
	if (IsFile(EigenstateFile[i][j]) == false)
	{
	  cout << "state " << EigenstateFile[i][j] << " does not exist or can't be opened" << endl;
	  return -1;           
	}
      }
    }
    
    delete[] StatisticPrefix;
    
    MultiColumnASCIIFile InputEnergies;
    if (InputEnergies.Parse(EnergyFileName) == false)
    {
      InputEnergies.DumpErrors(cout) << endl;
      return -1;
    }
    for (int j = 0; j < TruncationIndex; ++j)
      Energies[j] = InputEnergies.GetAsDoubleArray(j + 1);
    delete[] EnergyFileName;
  }
  
  
  
  Abstract2DTightBindingModel *TightBindingModel;
  
  TightBindingModel = new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), true, false);
  
  
  ParticleOnSphere* Space = 0;
  AbstractQHEHamiltonian* Hamiltonian = 0;
  
  int NbrSites = TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand();
  
  if (TranslationFlag)
    {
      if ((SzSymmetryFlag) && (Sz == 0))
	{				      
#ifdef __64_BITS__
	  if ( NbrSites < 31)
#else
	    if (NbrSites < 15)
#endif
	      {
		if (UsingConstraintMinNbrSinglets == false)
		  {
		    Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (NbrParticles, Sz,NbrSites ,(SzSymmetry == -1), Kx, NbrCellX, Ky,  NbrCellY );	  
		  }
		else
		  {
		    Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets (NbrParticles, MinNbrSinglets, Sz, NbrSites, (SzSymmetry == -1),   Kx, NbrCellX, Ky,  NbrCellY, 10000000);
		  } 
	      }
	    else
	      {
		if (UsingConstraintMinNbrSinglets == false)
		  {
		    Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong (NbrParticles, Sz, NbrSites, (SzSymmetry == -1),Kx, NbrCellX, Ky,  NbrCellY );	  
		  }
		else
		  {
		    Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong (NbrParticles, MinNbrSinglets, Sz, NbrSites , (SzSymmetry == -1),   Kx, NbrCellX, Ky,  NbrCellY);
		  }
	      }
	}
      else
	{
#ifdef __64_BITS__
	  if (NbrSites < 31)
#else
	    if ( NbrSites < 15)
#endif
	      {
		if (UsingConstraintMinNbrSinglets == false)
		  {
		    Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, Sz, NbrSites, Kx, NbrCellX, Ky,  NbrCellY);	  
		  }
		else
		  {
		    Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets (NbrParticles, MinNbrSinglets, Sz,  NbrSites , Kx, NbrCellX, Ky,  NbrCellY, 10000000ul);
		  }
	      }
	    else
	      {
		if (UsingConstraintMinNbrSinglets == false)
		  {
		    Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong (NbrParticles, Sz,  NbrSites, Kx, NbrCellX, Ky,  NbrCellY);	  
		  }
		else
		  {
		    Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong (NbrParticles, MinNbrSinglets, Sz,  NbrSites, Kx, NbrCellX, Ky,  NbrCellY, 10000000ul);
		  }
	      }
	}
    }
  else
    {
      if ((SzSymmetryFlag) && (Sz == 0))
	{
	  Space = new FermionOnLatticeWithSpinSzSymmetryRealSpace (NbrParticles, Sz,  NbrSites, (SzSymmetry == -1));	  
	}
      else
	{
	  Space = new FermionOnLatticeWithSpinRealSpace (NbrParticles, Sz, NbrSites,10000000);	  
	}
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
  

  RealSymmetricMatrix DensityDensityInteractionupup(NbrSites, true);
  RealSymmetricMatrix DensityDensityInteractiondowndown(NbrSites, true);
  RealSymmetricMatrix DensityDensityInteractionupdown(NbrSites, true);
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    Memory = Architecture.GetArchitecture()->GetLocalMemory();
      

  char* OutputNamePrefix = new char [512];
  if (ResumeFlag == false)
  {
    OutputNamePrefix =  RemoveExtensionFromFileName(StateFileName,".vec");
    StepI = -1;
  }
  else
  {
    char* ParameterStringPreviousSimulation = new char[512];
    sprintf(ParameterStringPreviousSimulation,"_uf_%f_tau_%f_tf_%f_nstep_%d.%d.vec",UFinal,TauI, TfI, NbrStepI, StepI);
    OutputNamePrefix = RemoveExtensionFromFileName(StateFileName, ParameterStringPreviousSimulation);
    delete[] ParameterStringPreviousSimulation;
  }

  char * ParameterString = new char [512];
  sprintf(ParameterString,"uf_%f_tau_%f_tf_%f_nstep_%d",UFinal,Tau,FinalTime,NbrSteps);
  char * EnergyNameFile = 0;
  char* NormFileName = 0;
  ofstream FileEnergy;
  ofstream FileNorm;
  if ((Manager.GetBoolean("compute-energy")) || (TruncationIndex > 0))
    {
      EnergyNameFile = new char[512];
      if (ResumeFlag == false)
      {
	if (TruncationIndex == 0)
	  sprintf (EnergyNameFile, "%s_%s_energy.dat", OutputNamePrefix,ParameterString);
	else
	  sprintf (EnergyNameFile, "%s_truncatedbasis_%d_%s_energy.dat", OutputNamePrefix, TruncationIndex, ParameterString);
      }
      else
	sprintf (EnergyNameFile, "%s_%s_energy.resume.dat", OutputNamePrefix,ParameterString);

      FileEnergy.open(EnergyNameFile, ios::out);
      FileEnergy.precision(14);
      FileEnergy << "# t E "<< endl;
    }
  if (TruncationIndex > 0)
  {
    NormFileName = new char[512];
    sprintf (NormFileName, "%s_truncatedbasis_%d_%s_norm.dat", OutputNamePrefix, TruncationIndex, ParameterString);
    FileNorm.open(NormFileName, ios::out);
    FileNorm.precision(14);
    FileNorm << "# t Norm "<< endl;
  }
  
  bool FinalHamiltonian = false;
  double FinalNorm = 1.0;
  for(int i = StepI + 1 ; i < NbrSteps; i++)
    {
      double UPotential =  UFinal/ Tau * TimeStep * i;
      if ( fabs(UPotential) > fabs(UFinal) ) 
	UPotential = UFinal ;
      
      for (int p = 0; p < NbrSites; ++p)
	{
	  DensityDensityInteractionupdown.SetMatrixElement(p, p, UPotential);
	}
      
      if ((FinalHamiltonian == false) && ((TruncationIndex == 0) || (Manager.GetBoolean("compute-energy")) || (i == NbrStepsTau - 1)))
      {
	cout << "build H" << endl;
	if (TranslationFlag)
	  {
	    Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian((ParticleOnSphereWithSpin  *)  Space, NbrParticles, NbrSites, Kx, NbrCellX, Ky,  NbrCellY,
											  TightBindingMatrix, TightBindingMatrix,
											  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
											  DensityDensityInteractionupdown, 
											  Architecture.GetArchitecture(), Memory);
	  }
	else
	  {
	    Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceHamiltonian((ParticleOnSphereWithSpin  *)  Space, NbrParticles, NbrSites,
									  TightBindingMatrix, TightBindingMatrix,
									  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
									  DensityDensityInteractionupdown, 
									  Architecture.GetArchitecture(), Memory);
	  }
      }

      double Norm;
      int TmpExpansionOrder;
      ComplexVector TmpState (Space->GetHilbertSpaceDimension()) ;
      ComplexVector TmpState1 (Space->GetHilbertSpaceDimension()) ;
      Complex TmpCoefficient;
      
      TmpState.Copy(TmpInitialState);
      Norm = TmpState.Norm();
      double TmpNorm = 1.0;
      cout << "Computing state " << (i + 1) << "/" <<  NbrSteps << " at t = " << (TimeStep * i) <<" with Interaction "<< UPotential  <<endl;
      
      if (TruncationIndex == 0)
      {
	TmpExpansionOrder = 0;
	TmpCoefficient = 1.0;
	while ( ((fabs(TmpNorm) > 1e-14 ) || (TmpExpansionOrder < 1)) && (TmpExpansionOrder <= Manager.GetInteger("iter-max")))
	  {
	    TmpExpansionOrder += 1;
	    TmpCoefficient = -TmpCoefficient * TimeStep * Complex(0.0, 1.0) / ((double) TmpExpansionOrder);
	    VectorHamiltonianMultiplyOperation Operation (Hamiltonian, (&TmpState), (&TmpState1));
	    Operation.ApplyOperation(Architecture.GetArchitecture());
	    TmpState.Copy(TmpState1);
	    TmpNorm = sqrt(TmpCoefficient.Re*TmpCoefficient.Re + TmpCoefficient.Im*TmpCoefficient.Im) * TmpState.Norm();
	    AddComplexLinearCombinationOperation Operation1 (&TmpInitialState, &TmpState1, 1, &TmpCoefficient);
	    Operation1.ApplyOperation(Architecture.GetArchitecture());
	    Norm = TmpInitialState.Norm();
	    cout << "Norm = " << Norm << " +/- " << TmpNorm << " for step " << TmpExpansionOrder << endl;      
	  } 
      }
      else
      {
	int TmpIndex = i;
	if (TmpIndex > NbrStepsTau - 1)
	  TmpIndex = NbrStepsTau - 1;
	TmpState1.ClearVector();
	Complex TmpOverlap;
	for (int j = 0; j < TruncationIndex; ++j)
	{
	  TmpState.ReadVector(EigenstateFile[TmpIndex][j]);
	  TmpOverlap = TmpState * TmpInitialState;
	  TmpNorm = sqrt(TmpOverlap.Re*TmpOverlap.Re + TmpOverlap.Im*TmpOverlap.Im);
	  TmpCoefficient = TmpOverlap * Phase(- TimeStep * Energies[j][TmpIndex]);
	  AddComplexLinearCombinationOperation Operation (&TmpState1, &TmpState, 1, &TmpCoefficient);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  Norm = TmpState1.Norm();
	  cout << "Norm = " << Norm << " +/- " << TmpNorm << " for vector " << j << endl;  
	}
	
	FinalNorm *= Norm; 
	FileNorm <<  (TimeStep * i) << " " << Norm << " " << FinalNorm << endl;
	TmpInitialState.Copy(TmpState1);
	TmpInitialState.Normalize();
      }

      if ((((i % SavePeriod) == 0) && (SavePeriod >= 0)) || ((SavePeriod < 0) && (i == NbrSteps - 1)))
      {
	char* OutputName = new char [strlen(OutputNamePrefix)+ strlen(ParameterString) + 32];
	if (TruncationIndex == 0)
	  sprintf (OutputName, "%s_%s.%d.vec",OutputNamePrefix,ParameterString,i);
	else
	  sprintf (OutputName, "%s_%s_truncatedbasis_%d.%d.vec",OutputNamePrefix,ParameterString, TruncationIndex, i);
	TmpInitialState.WriteVector(OutputName);      
	delete [] OutputName;
      }

      if ((Manager.GetBoolean("compute-energy")) || ((TruncationIndex > 0) && (i >= NbrStepsTau)))
	{
	  VectorHamiltonianMultiplyOperation Operation (Hamiltonian, (&TmpInitialState), (&TmpState));
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  Complex Energy = (TmpInitialState) * (TmpState);
	  FileEnergy <<  (TimeStep * i) << " " << Energy << endl;
	  cout << "E = " << Energy << endl;
	}
      if (UPotential == UFinal) 
	FinalHamiltonian = true;
      
      if (((FinalHamiltonian == false) && (Hamiltonian != 0)) || (i == NbrSteps - 1))
	delete  Hamiltonian;
    }  
  delete[] OutputNamePrefix;
  delete[]  ParameterString;
  delete TightBindingModel;
  delete Space;
  if (EnergyNameFile != 0)
    delete[] EnergyNameFile;
  if (NormFileName != 0)
    delete[] NormFileName;
  
  if (EigenstateFile != 0)
  {
    for(int i = StepI + 1 ; i < NbrStepsTau; i++)
    {
      for (int j = 0; j < TruncationIndex; ++j)
	delete[] EigenstateFile[i][j];
      delete[] EigenstateFile[i];
    }
    delete[] EigenstateFile;
  }
   
  if (Energies != 0)
  {
    for (int i = 0; i < TruncationIndex; ++i)
      delete[] Energies[i];
    delete[] Energies;
  }
  
  return 0;
}


// try to guess system information from file name
//
// filename = vector file name
// tauI = reference to the value of tau
// tfI = 
// nbrstepI = number of steps in previous simulation
// stepI = reference to the number of steps already done
// return value = true if no error occured

bool FTITimeEvolutionFindSystemInfoFromVectorFileName (char* filename, double& tauI, double& tfI, int& nbrStepI,int& stepI, double& uFinalI)
{
  char* StrTau;
  StrTau = strstr(filename, "_tau_");
  if (StrTau != 0)
    {
      StrTau += 5;
      int SizeString = 0;
      while ((StrTau[SizeString] != '\0') && (StrTau[SizeString] != '_') && (StrTau[SizeString] <= '9'))
	++SizeString;
      if ((StrTau[SizeString] == '_') && (SizeString != 0))
	{
          char TmpChar = StrTau[SizeString];
	  StrTau[SizeString] = '\0';
	  tauI = atoi(StrTau);
	  StrTau[SizeString] = TmpChar;
	  StrTau += SizeString;
	}
      else
	StrTau = 0;
    }
  if (StrTau == 0)
    {
      return false;            
    }
  StrTau = strstr(filename, "_tf_");
  if (StrTau != 0)
    {
      StrTau += 4;
      int SizeString = 0;
      while ((StrTau[SizeString] != '\0') && (StrTau[SizeString] != '_') && (StrTau[SizeString] <= '9'))
	++SizeString;
      if ((StrTau[SizeString] == '_') && (SizeString != 0))
	{
          char TmpChar = StrTau[SizeString];
	  StrTau[SizeString] = '\0';
	  tfI = atoi(StrTau);
// 	  cout << "tfi = " << tfI << endl;
	  StrTau[SizeString] = TmpChar;
	  StrTau += SizeString;
	}
      else
	StrTau = 0;
    }
  if (StrTau == 0)
    {
      return false;
    }
  StrTau = strstr(filename, "_nstep_");
  if (StrTau != 0)
    {
      StrTau += 7;
      int SizeString = 0;
      while ((StrTau[SizeString] != '\0') && (StrTau[SizeString] != '.') && (StrTau[SizeString] >= '0') 
	     && (StrTau[SizeString] <= '9'))
	++SizeString;
      if ((StrTau[SizeString] == '.') && (SizeString != 0))
	{
          char TmpChar = StrTau[SizeString];
	  StrTau[SizeString] = '\0';
	  nbrStepI = atoi(StrTau);
// 	  cout << "NbrStepI = " << nbrStepI << endl;
	  StrTau[SizeString] = TmpChar;
	  StrTau += SizeString;
	}
      else
	StrTau = 0;
      
      StrTau += 1;
      SizeString = 0;
      while ((StrTau[SizeString] != '\0') && (StrTau[SizeString] != '.') 
	     && (StrTau[SizeString] <= '9'))
	++SizeString;
      if ((StrTau[SizeString] == '.') && (SizeString != 0))
	{
          char TmpChar = StrTau[SizeString];
	  StrTau[SizeString] = '\0';
	  stepI = atoi(StrTau);
// 	  cout << "stepI  = " << stepI << endl;
	  StrTau[SizeString] = TmpChar;
	  StrTau += SizeString;
	}
      else
	StrTau = 0;
    }
  if (StrTau == 0)
    {
      return false;
    }
  
  
  StrTau = strstr(filename, "_uf_");
  if (StrTau != 0)
    {
      StrTau += 4;
      int SizeString = 0;
      while ((StrTau[SizeString] != '\0') && (StrTau[SizeString] != '_') && (StrTau[SizeString] <= '9'))
	++SizeString;
      if ((StrTau[SizeString] == '_') && (SizeString != 0))
	{
          char TmpChar = StrTau[SizeString];
	  StrTau[SizeString] = '\0';
	  uFinalI = atoi(StrTau);
// 	  cout << "tau = " << tauI << endl;
	  StrTau[SizeString] = TmpChar;
	  StrTau += SizeString;
	}
      else
	StrTau = 0;
    }
  
  return true;
}


// create file containing the necessary energies for time-evolution with Hilbert-space truncation
//
// statisticPrefixEnergy
// kx = momentum along x
// ky = momentum along y
// sz = total magnetization
// szParity = parity under spin-inversion
// truncationIndex = number of low-energy states that are conserved for the time-evolution
char* FTITimeEvolutionCreateTruncatedEnergyFile (char* statisticPrefixEnergy, double gammaX, double gammaY, int kx, int ky, int sz, int szParity, int truncationIndex, double uFinal, double tau, double finalTime, int nbrSteps, bool createFile)
{
  char* SpectrumFileName = new char[256 + strlen(statisticPrefixEnergy)];
  double uFactor = uFinal/ tau * finalTime / ((double) nbrSteps);
  sprintf(SpectrumFileName, "%s_gx_%g_gy_%g_sz_%d_szsym_%d_kx_%d_ky_%d_uf_%f_tau_%f_tf_%f_nstep_%d_truncatedbasis_%d.energies.dat", statisticPrefixEnergy, gammaX, gammaY, sz, szParity, kx, ky, uFinal, tau, finalTime, nbrSteps, truncationIndex);
  if (createFile == false)
    return SpectrumFileName;
  
  char* TmpFileName = new char[256 + strlen(statisticPrefixEnergy)];
  ofstream EnergyFile;
  EnergyFile.open(SpectrumFileName, ios::out);
  int lineIndex;
  int nbrEnergies;
  for(int i = 0 ; i < nbrSteps; i++)
    {
      double UPotential =  uFactor * i;
      if (i == 0)
	UPotential = 0.0;
      if ( fabs(UPotential) > fabs(uFinal) ) 
	UPotential = uFinal ;
      EnergyFile << UPotential << " ";
      sprintf(TmpFileName, "%s_u_%g_gx_%g_gy_%g.dat", statisticPrefixEnergy, UPotential, gammaX, gammaY);
      MultiColumnASCIIFile TmpFile;
      if (TmpFile.Parse(TmpFileName) == false)
      {
	TmpFile.DumpErrors(cout) << endl;
	return 0;
      }
      lineIndex = 0;
      while ((atoi(TmpFile(0, lineIndex)) != sz) || (atoi(TmpFile(1, lineIndex)) != szParity) || (atoi(TmpFile(2, lineIndex)) != kx) || (atoi(TmpFile(3, lineIndex)) != ky))
	++lineIndex;
     
      nbrEnergies = 0;
      for (int j = 0; j < truncationIndex;  ++j)
      {
	if ((atoi(TmpFile(0, lineIndex)) == sz) && (atoi(TmpFile(1, lineIndex)) == szParity) && (atoi(TmpFile(2, lineIndex)) == kx) && (atoi(TmpFile(3, lineIndex)) == ky))
	{
	  EnergyFile << TmpFile(4, lineIndex) << " ";
	  ++nbrEnergies;
	}
	++lineIndex;
      }
      EnergyFile << endl;    
      if (nbrEnergies < truncationIndex)
	return 0;
    }
    
  EnergyFile.close();
  delete[] TmpFileName;
  return SpectrumFileName;
}
