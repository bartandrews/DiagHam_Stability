#include "Options/Options.h"

#include "HilbertSpace/BosonOnLatticeRealSpace.h"
#include "HilbertSpace/BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong.h"


#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong.h"
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


#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceHamiltonian.h"
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

int main(int argc, char** argv)
{
  OptionManager Manager ("FCIHofstadterModelPeriodicDriving" , "0.01");
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
  (*SystemGroup) += new BooleanOption  ('\n', "compute-energy", "compute the energy of each time-evolved vector");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-hamiltonian", "number of hamiltonians to store (equal to number of samples if 0)", 1);
  
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "amplitude of on-site interaction (if hard core constraint not implemented)", 0.0);
  
  (*SystemGroup) += new SingleDoubleOption  ('\n', "omega", "frequency of drive", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "e-field", "amplitude of the driving field", 1.0);
  
  (*SystemGroup) += new BooleanOption  ('\n', "compute-FGR", "compute the absorption rate from Fermi's Golden rule");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the Chern number of the fully filled band (only available in singleparticle-spectrum mode)");
  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-periods", "number of periods to do the time-evolution with", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-samples", "number of points to evaluate per period", 10);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "save-period", "save time-evolved state every n periods, save only final step if negative, save nothing if 0", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "time-step", "time step for the time evolution", 0.1);
  
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
  int Kx,Ky;
  bool StatisticFlag;
  bool GutzwillerFlag;
  bool TranslationFlag;
  bool SpinFlag = false;
  bool SzSymmetryFlag = false;
  int Sz;
  int SzSymmetry;
  char Axis ='y';
  if (Manager.GetBoolean("landau-x"))
    Axis ='x';
  double UPotential = Manager.GetDouble("u-potential");
  int TruncationIndex = Manager.GetInteger("truncate-basis");
  
  double Omega = fabs(Manager.GetDouble("omega"));
  double Sign = Manager.GetDouble("omega") / Omega;
  double EField = Manager.GetDouble("e-field");
  
  int NbrSamples = Manager.GetInteger("nbr-samples");
  int NbrPeriods = Manager.GetInteger("nbr-periods");
  int SavePeriod = Manager.GetInteger("save-period");
  double TimePeriod = 2.0 * M_PI / ((double) Omega);
  double TimeStep = TimePeriod / ((double) NbrSamples);
  double TimeF;
  
  bool ResumeFlag = false;
  int NbrStoredHamiltonians = Manager.GetInteger("nbr-hamiltonian");
  if (NbrStoredHamiltonians == 0)
    NbrStoredHamiltonians = NbrSamples;
  if (NbrStoredHamiltonians < 0)
    NbrStoredHamiltonians = 0;
  
  if (NbrPeriods == 0)
  {
    cout << "Working with a time step not commensurate with the period of time-evolution" << endl;
    TimeStep = Manager.GetDouble("time-step");
    TimeF = TimeStep * NbrSamples;
    NbrStoredHamiltonians = 0;
  }
  
  int StepI;
  
  double GammaX = Manager.GetDouble("gamma-x");
  double GammaY = Manager.GetDouble("gamma-y");
  
  
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

  int MinNbrSinglets;
  bool UsingConstraintMinNbrSinglets;
  if (StatisticFlag)
  {
    if (FTIHofstadterModelWithSzFindSystemInfoFromVectorFileName(StateFileName, Sz, SzSymmetry, MinNbrSinglets, SzSymmetryFlag, UsingConstraintMinNbrSinglets) == false)
    {
      return -1;
    }
    SpinFlag = true;
    cout << Sz << " "  << SzSymmetryFlag << endl;
  }
    
    if (StatisticFlag)
      cout << "fermion" << endl;
    else
      cout << "boson" << endl;
    if (SpinFlag)
      cout << "spin" << endl;

    
  char*** EigenstateFile = 0;
  double** Energies = 0;
  if (TruncationIndex > 0)
  {
    cout << "truncated hilbert space time evolution is not implemented" << endl;
    return -1;
  }
  
  cout << NbrCellX << " " << NbrCellY << " " << UnitCellX << " " << UnitCellY << " " << FluxPerCell << " " << GammaX << " " << GammaY << endl;
  
  
  Abstract2DTightBindingModel *TightBindingModel;
  
  TightBindingModel = new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis, GammaX, GammaY, Architecture.GetArchitecture(), true, false);
  
  ParticleOnSphere* Space = 0;  
  int NbrSites = NbrCellX * NbrCellY * UnitCellX * UnitCellY;
  
  if (StatisticFlag)
  {
    if (SpinFlag == false)
    {
      if (TranslationFlag)
      {
	cout << "Case not implemented yet" << endl;
	return -1;
      }
      else
      {
	Space =  new FermionOnLatticeRealSpace (NbrParticles, NbrSites);
      }
    }
    else
    {
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
		Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (NbrParticles, Sz,NbrSites, (SzSymmetry == -1), Kx, NbrCellX, Ky,  NbrCellY);
	      }
	    else
	      {
		Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong (NbrParticles, Sz, NbrSites, (SzSymmetry == -1),Kx, NbrCellX, Ky,  NbrCellY );
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
		Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, Sz, NbrSites, Kx, NbrCellX, Ky,  NbrCellY);	  
	      }
	    else
	      {
		Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong (NbrParticles, Sz,  NbrSites, Kx, NbrCellX, Ky,  NbrCellY);	   
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
    }
  }
  else
  {
    if (TranslationFlag)
    {
      if (GutzwillerFlag)
      {
	if (UnitCellX*NbrCellX*  UnitCellY*NbrCellY <= 63 ) 
	  {
	    Space = new BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation(NbrParticles, UnitCellX*NbrCellX, UnitCellY*NbrCellY, Kx, NbrCellX, Ky,  NbrCellY);
	  }
	else
	  {
	    Space = new BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong(NbrParticles, UnitCellX*NbrCellX, UnitCellY*NbrCellY, Kx, NbrCellX, Ky,  NbrCellY);      
	  }
      }
      else
      {
	Space = new BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation(NbrParticles, UnitCellX*NbrCellX, UnitCellY*NbrCellY, Kx, NbrCellX, Ky,  NbrCellY); 
      }
    }
    else
    {
      if (GutzwillerFlag)
      {
	Space = new BosonOnLatticeGutzwillerProjectionRealSpace(NbrParticles, NbrSites);
      }
      else
      {
	Space = new BosonOnLatticeRealSpace(NbrParticles, NbrSites);
      }
    }
  }

  
  ComplexVector TmpInitialState;
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
  
  TmpInitialState.Copy(InputState);
  
  

  RealSymmetricMatrix DensityDensityInteraction(NbrSites, true);
  
  RealSymmetricMatrix DensityDensityInteractionupup(NbrSites, true);
  RealSymmetricMatrix DensityDensityInteractiondowndown(NbrSites, true);
  RealSymmetricMatrix DensityDensityInteractionupdown(NbrSites, true);
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    Memory = Architecture.GetArchitecture()->GetLocalMemory();
        
  if (GutzwillerFlag == false)
  {
    for (int x = 0; x <  NbrCellX; ++x)
    {
      for (int y = 0; y <  NbrCellY; ++y)
      {
	for (int k = 0; k < TightBindingModel->GetNbrBands(); ++k)
	{
	  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y, k),TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y, k), UPotential);
	}
      }
    }
  }
  
  char* OutputNamePrefix = new char [512];
  if (ResumeFlag == false)
  {
    OutputNamePrefix =  RemoveExtensionFromFileName(StateFileName,".vec");
    StepI = -1;
  }
  else
  {
    cout << "Resume option is not implemented" << endl;
//     char* ParameterStringPreviousSimulation = new char[512];
//     sprintf(ParameterStringPreviousSimulation,"E_%f_omega_%f_nbrSamples_%d_nbrPeriods_%d.%d.vec", EField, Manager.GetDouble("omega"), NbrSamples, NbrPeriods);
//     OutputNamePrefix = RemoveExtensionFromFileName(StateFileName, ParameterStringPreviousSimulation);
//     delete[] ParameterStringPreviousSimulation;
  }

  char * ParameterString = new char [512];
  if (NbrPeriods != 0)
    sprintf(ParameterString,"E_%f_omega_%f_nbrSamples_%d_nbrPeriods_%d", EField, Manager.GetDouble("omega"), NbrSamples, NbrPeriods);
  else
    sprintf(ParameterString,"E_%f_omega_%f_timeF_%f_timeStep_%f", EField, Manager.GetDouble("omega"), TimeF, TimeStep);
  char * EnergyNameFile = 0;
  char* EnergyNameFile0 = 0;
  char* NormFileName = 0;
  ofstream FileEnergy;
  ofstream FileEnergy0;
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
 
    EnergyNameFile0 = new char[512];
    sprintf (EnergyNameFile0, "%s_%s_H0_energy.dat", OutputNamePrefix, ParameterString);
    FileEnergy0.open(EnergyNameFile0, ios::out);
    FileEnergy0.precision(14);
    FileEnergy0 << "# t E "<< endl;
    
    NormFileName = new char[512];
    if (TruncationIndex > 0)
      sprintf (NormFileName, "%s_truncatedbasis_%d_%s_norm.dat", OutputNamePrefix, TruncationIndex, ParameterString);
    else
      sprintf (NormFileName, "%s_%s_overlap.dat", OutputNamePrefix,ParameterString);
    FileNorm.open(NormFileName, ios::out);
    FileNorm.precision(14);
    FileNorm << "# t Norm " << endl;
  
  
  double FinalNorm = 1.0;
  int NbrSteps = NbrSamples * NbrPeriods;
  if (NbrPeriods == 0)
    NbrSteps = NbrSamples;
  double TmpGammaX;
  double TmpGammaY;
  double t;
  double FluxDensity =  (((double) FluxPerCell)/( (double) (UnitCellX*UnitCellY)));
  double PhaseTranslationX = 2.0* M_PI * FluxDensity * UnitCellX;
  double PhaseTranslationY = 0.0;
  Complex TmpEnergy;
  
  delete TightBindingModel;
  HermitianMatrix* TightBindingMatrix = new HermitianMatrix [NbrSamples];
  AbstractQHEHamiltonian** Hamiltonian;
  if (NbrStoredHamiltonians == NbrSamples)
    Hamiltonian= new AbstractQHEHamiltonian* [NbrStoredHamiltonians];
  else
    Hamiltonian= new AbstractQHEHamiltonian* [NbrStoredHamiltonians + 1];
  
  for(int i = 0 ; i < NbrSamples; i++)
    {
      t = i * TimeStep;
      
      TmpGammaX = GammaX - EField * NbrCellX * sin(Omega * t) / (Omega * M_PI);
      TmpGammaY = GammaY + Sign * NbrCellY * EField * cos(Omega * t) / (Omega * M_PI);
      
      TightBindingModel = new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis, TmpGammaX, TmpGammaY, Architecture.GetArchitecture(), true, false);      
      TightBindingMatrix[i] = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
      delete TightBindingModel;
      
      if (i < NbrStoredHamiltonians)
      {
	cout << "build H for t = " << t << " + T" << endl;
	if (SpinFlag)
	{
	  Hamiltonian[i] = new ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian((ParticleOnSphereWithSpin  *)  Space, NbrParticles, NbrSites, Kx, NbrCellX, Ky,  NbrCellY,
											  TightBindingMatrix[i], TightBindingMatrix[i],
											  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
											  DensityDensityInteractionupdown, 
											  Architecture.GetArchitecture(), Memory);
	}
	else
	{
	  if (TranslationFlag)
	  {
	    cout << NbrParticles << " " << NbrSites << " " << Kx << " " << NbrCellX << " " << Ky << " " << NbrCellY << endl;
	    
	    Hamiltonian[i] = new ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian (Space, NbrParticles, NbrSites, Kx,  NbrCellX, Ky,  NbrCellY, PhaseTranslationX, PhaseTranslationY, TightBindingMatrix[i], DensityDensityInteraction, Architecture.GetArchitecture(), Memory);
	  }
	  else
	  {
	    Hamiltonian[i] = new ParticleOnLatticeRealSpaceHamiltonian (Space, NbrParticles, NbrSites, TightBindingMatrix[i], DensityDensityInteraction, Architecture.GetArchitecture(), Memory);
	  }
	}
      }    
    }  
   
   if (NbrStoredHamiltonians == NbrSamples)
    delete[] TightBindingMatrix;
  
  int TmpIndex;
  Complex TmpOverlap;
  for(int i = StepI + 1; i < NbrSteps; i++)
    {
      t = i * TimeStep;
      TmpIndex = i % NbrSamples;

      double Norm;
      int TmpExpansionOrder;
      ComplexVector TmpState (Space->GetHilbertSpaceDimension()) ;
      ComplexVector TmpState1 (Space->GetHilbertSpaceDimension()) ;
      Complex TmpCoefficient;
      
      TmpState.Copy(TmpInitialState);
      Norm = TmpState.Norm();
      double TmpNorm = 1.0;
      cout << "Computing state " << (i + 1) << "/" <<  NbrSteps << " at t = " << t << endl;
      
      if (TmpIndex >= NbrStoredHamiltonians)
      {
	if (SpinFlag)
	{
	  Hamiltonian[NbrStoredHamiltonians] = new ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian((ParticleOnSphereWithSpin  *)  Space, NbrParticles, NbrSites, Kx, NbrCellX, Ky,  NbrCellY,
											  TightBindingMatrix[TmpIndex], TightBindingMatrix[TmpIndex],
											  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
											  DensityDensityInteractionupdown, 
											  Architecture.GetArchitecture(), Memory);
	}
	else
	{
	  if (TranslationFlag)
	  {
	    Hamiltonian[NbrStoredHamiltonians] = new ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian (Space, NbrParticles, NbrSites, Kx,  NbrCellX, Ky,  NbrCellY, 	PhaseTranslationX, PhaseTranslationY, TightBindingMatrix[TmpIndex], DensityDensityInteraction, Architecture.GetArchitecture(), Memory);
	  }
	  else
	  {
	    Hamiltonian[NbrStoredHamiltonians] = new ParticleOnLatticeRealSpaceHamiltonian (Space, NbrParticles, NbrSites, TightBindingMatrix[TmpIndex], DensityDensityInteraction, Architecture.GetArchitecture(), Memory);
	  }
	}
	TmpIndex = NbrStoredHamiltonians;
      }
	  
      if (TruncationIndex == 0)
      {
	TmpExpansionOrder = 0;
	TmpCoefficient = 1.0;
	while ( ((fabs(TmpNorm) > 1e-14 ) || (TmpExpansionOrder < 1)) && (TmpExpansionOrder <= Manager.GetInteger("iter-max")))
	  {
	    TmpExpansionOrder += 1;
	    TmpCoefficient = -TmpCoefficient * TimeStep * Complex(0.0, 1.0) / ((double) TmpExpansionOrder);
	    VectorHamiltonianMultiplyOperation Operation (Hamiltonian[TmpIndex], (&TmpState), (&TmpState1));
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
	TmpState1.ClearVector();
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

      if ((SavePeriod != 0) && ((((i % SavePeriod) == 0) && (SavePeriod > 0)) || ((SavePeriod < 0) && (i >= NbrSteps - NbrSamples))))
      {
	char* OutputName = new char [strlen(OutputNamePrefix)+ strlen(ParameterString) + 32];
	if (TruncationIndex == 0)
	  sprintf (OutputName, "%s_%s.%d.vec",OutputNamePrefix,ParameterString,i);
	else
	  sprintf (OutputName, "%s_%s_truncatedbasis_%d.%d.vec",OutputNamePrefix,ParameterString, TruncationIndex, i);
	TmpInitialState.WriteVector(OutputName);      
	delete [] OutputName;
      }

      if (Manager.GetBoolean("compute-energy"))
	{
	  VectorHamiltonianMultiplyOperation Operation (Hamiltonian[TmpIndex], (&TmpInitialState), (&TmpState));
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  Complex Energy = (TmpInitialState) * (TmpState);
	  FileEnergy <<  (TimeStep * i) << " " << Energy.Re << endl;
	  cout << "E(t) = " << Energy << endl;
	  
	  VectorHamiltonianMultiplyOperation Operation1 (Hamiltonian[0], (&TmpInitialState), (&TmpState));
	  Operation1.ApplyOperation(Architecture.GetArchitecture());
	  TmpEnergy = TmpInitialState * TmpState;
	  FileEnergy0 << (TimeStep * i) << " " << (TmpEnergy.Re) << endl;
	  cout << "E_0 = " << (TmpEnergy.Re) << endl;
	  
	}
      
      TmpOverlap = InputState * TmpInitialState;
      cout << TmpOverlap << endl;
      FileNorm << (TimeStep * i) << " " << (sqrt(TmpOverlap.Re*TmpOverlap.Re + TmpOverlap.Im*TmpOverlap.Im)) << endl;
      
      
      
      if (TmpIndex == NbrStoredHamiltonians)
	delete  Hamiltonian[NbrStoredHamiltonians];
    }  
    
  
   for (int i = 0; i < NbrStoredHamiltonians; ++i)
     delete Hamiltonian[i];
    delete[] Hamiltonian;
  
  delete[] OutputNamePrefix;
  delete[]  ParameterString;
  delete Space;
  if (EnergyNameFile != 0)
    delete[] EnergyNameFile;
  if (NormFileName != 0)
    delete[] NormFileName;
  if (EnergyNameFile0 != 0)
    delete[] EnergyNameFile0;
    
  return 0;
}


