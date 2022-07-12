#include "Options/Options.h"

#include "HilbertSpace/BosonOnLatticeRealSpace.h"
#include "HilbertSpace/BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong.h"


#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
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

#include "Tools/FTITightBinding/TightBindingModelCheckerboardLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
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

bool FTITimeEvolutionFindSystemInfoFromVectorFileName (char* filename, double& eFieldI, double& omegaI, int& nbrSamplesI,int& nbrPeriodsI, int& stepI);

int main(int argc, char** argv)
{
  OptionManager Manager ("FCICheckerboardModelPeriodicDriving" , "0.01");
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
  
//   (*SystemGroup) += new SingleIntegerOption  ('\n', "nx1", "first coordinate of the first spanning vector of the tilted lattice", 0);
//   (*SystemGroup) += new SingleIntegerOption  ('\n', "ny1", "second coordinate of the first spanning vector of the tilted lattice", 0);
//   (*SystemGroup) += new SingleIntegerOption  ('\n', "nx2", "first coordinate of the second spanning vector of the tilted lattice", 0);
//   (*SystemGroup) += new SingleIntegerOption  ('\n', "ny2", "second coordinate of the second spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "real-offset", "second coordinate in real space of the second spanning vector of the real space lattice (0 if lattice is untilted)", 0);
  
  
    
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
  int NbrSiteX; 
  int NbrSiteY;
  int NbrSites;
  int Kx,Ky;
  bool StatisticFlag;
  bool GutzwillerFlag;
  bool TranslationFlag;
 
  
//   int nx1 = Manager.GetInteger("nx1");
//   int ny1 = Manager.GetInteger("ny1");
//   int nx2 = Manager.GetInteger("nx2");
//   int ny2 = Manager.GetInteger("ny2");
  int OffsetReal = Manager.GetInteger("real-offset");
  int nx1 = 0;
  int ny1 = 0;
  int nx2 = 0;
  int ny2 = 0;
  int Offset = 0;
  bool TiltedFlag = true;
  
  
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
  
  if (FTIHubbardModelFindSystemInfoFromVectorFileName(StateFileName, NbrParticles, NbrSites, StatisticFlag, GutzwillerFlag) == false)
  {
    cout << "error while retrieving system parameters from file name " << StateFileName << endl;
    return -1;
  
  }
  TranslationFlag = FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(StateFileName, NbrParticles, NbrSites, Kx, Ky, 
											     NbrSiteX, NbrSiteY, StatisticFlag, GutzwillerFlag);
    
	
  
  double GammaX = 0.0;
  double GammaY = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double tpp = 0.0;
  double muS = 0.0;
  double UPotential = 0.0;
  double VPotential = 0.0;
  
  if ( !FilenameDoubleSearch(t1, StateFileName, "_t1_") || (!FilenameDoubleSearch(t2, StateFileName, "_t2_")) || (!FilenameDoubleSearch(tpp, StateFileName, "_tpp_")))
    return -1;
  
  if (!FilenameDoubleSearch(UPotential, StateFileName, "_u_"))
    UPotential = 0.0;
  if (!FilenameDoubleSearch(VPotential, StateFileName, "_v_"))
    VPotential = 0.0;
   
  if (fabs(t2 - (1.0 - 0.5 * M_SQRT2)) < 1.0e-5)
    t2 = 1.0 - 0.5 * M_SQRT2;
  if (fabs(tpp - 0.5 * (M_SQRT2 - 1.0)) < 1.0e-5)
    tpp = 0.5 * (M_SQRT2 - 1.0);
  
  if (!FilenameDoubleSearch(muS, StateFileName, "_muS_"))
    muS = 0.0;
    
  if (StatisticFlag)
    cout << "fermion" << endl;
  else
  cout << "boson" << endl;
  
  if (!FilenameIntegerSearch(nx1, StateFileName, "_nx1_") || !FilenameIntegerSearch(ny1, StateFileName, "_ny1_") || !FilenameIntegerSearch(nx2, StateFileName, "_nx2_") || !FilenameIntegerSearch(ny2, StateFileName, "_ny2_"))
    TiltedFlag = false;
          
  char*** EigenstateFile = 0;
  double** Energies = 0;
  if (TruncationIndex > 0)
  {
    cout << "truncated hilbert space time evolution is not implemented" << endl;
    return -1;
  }
  
  cout << NbrParticles << " " << NbrSites << " " << NbrSiteX << " " << NbrSiteY << " " << Kx << " " << Ky << " " << GammaX << " " << GammaY << endl;
  
  cout << "t1 = " << t1 << " t2 = " << t2 << " tpp = " << tpp << " U = " << UPotential << " V = " << VPotential << endl;
  if (TiltedFlag)
    cout << "nx1 = " << nx1 << " ny1 = " << ny1 << " nx2 = " << nx2 << " ny2 = " << ny2 << " OffsetReal = " << OffsetReal << endl;
  Abstract2DTightBindingModel *TightBindingModel;
  if (!TiltedFlag)
    TightBindingModel = new TightBindingModelCheckerboardLattice (NbrSiteX, NbrSiteY, t1, t2, tpp, muS, GammaX, GammaY, Architecture.GetArchitecture(), true, false);
  else
    TightBindingModel = new TightBindingModelCheckerboardLattice (NbrSiteX, NbrSiteY, nx1, ny1, nx2, ny2, Offset, t1, t2, tpp, 
							     muS, GammaX, GammaY, Architecture.GetArchitecture(), OffsetReal, true);
    
  if (Manager.GetBoolean("singleparticle-chernnumber") == true)
  {
      cout << "Chern number = " << TightBindingModel->ComputeChernNumber(0) << endl;
      return 0;
  }
  
  ParticleOnSphere* Space = 0;  
  
  if (StatisticFlag)
  {
    if (TranslationFlag)
      Space = new FermionOnLatticeRealSpaceAnd2DTranslation (NbrParticles, NbrSites, Kx, NbrSiteX, Ky, NbrSiteY);
    else
      Space = new FermionOnLatticeRealSpace (NbrParticles, NbrSites);
  }
  else
  {
    if (TranslationFlag)
    {
      if (GutzwillerFlag)
      {
	if (NbrSites <= 63 ) 
	  {
	    
	    Space = new BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation(NbrParticles, NbrSites, Kx, NbrSiteX, Ky, NbrSiteY);
	  }
	else
	  {
	    cout << "Number of orbitals is larger than 64" << endl;
	    return -1;     
	  }
      }
      else
      {
	Space = new BosonOnLatticeRealSpaceAnd2DTranslation(NbrParticles, NbrSites, Kx, NbrSiteX, Ky, NbrSiteY);
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

  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    Memory = Architecture.GetArchitecture()->GetLocalMemory();
        
  if (GutzwillerFlag == false)
  {
    for (int x = 0; x <  NbrSiteX; ++x)
    {
      for (int y = 0; y <  NbrSiteY; ++y)
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
    cout << OutputNamePrefix << endl;
  }
  else
  {
    cout << "Resume option is not implemented" << endl;
//     char* ParameterStringPreviousSimulation = new char[512];
//     if (NbrPeriods != 0)
//       sprintf(ParameterStringPreviousSimulation,"E_%f_omega_%f_nbrSamples_%d_nbrPeriods_%d.%d.vec", EField, Manager.GetDouble("omega"), NbrSamples, NbrPeriods);
//     else
//       sprintf(ParameterStringPreviousSimulation,"E_%f_omega_%f_timeF_%f_timeStep_%f.%d.vec", EField, Manager.GetDouble("omega"), TimeF, TimeStep);
//     
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

      cout << EnergyNameFile << endl;
      FileEnergy.open(EnergyNameFile, ios::out);
      FileEnergy.precision(14);
      FileEnergy << "# t E "<< endl;
    }
 
    EnergyNameFile0 = new char[512];
    sprintf (EnergyNameFile0, "%s_%s_H0_energy.dat", OutputNamePrefix, ParameterString);
    FileEnergy0.open(EnergyNameFile0, ios::out);
    FileEnergy0.precision(14);
    FileEnergy0 << "# t E " << endl;
    
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
  Complex TmpEnergy;
  
  HermitianMatrix* TightBindingMatrix = new HermitianMatrix [NbrSamples];
  AbstractQHEHamiltonian** Hamiltonian;
  if (NbrStoredHamiltonians == NbrSamples)
    Hamiltonian= new AbstractQHEHamiltonian* [NbrStoredHamiltonians];
  else
    Hamiltonian= new AbstractQHEHamiltonian* [NbrStoredHamiltonians + 1];
  
  for(int i = 0 ; i < NbrSamples; i++)
    {
      t = i * TimeStep;
      
//       TmpGammaX = GammaX + Sign * EField * cos(2.0 * M_PI * Omega * t) / (Omega * M_PI);
//       TmpGammaY = GammaY - EField * sin(2.0 * M_PI * Omega * t) / (Omega * M_PI);
      
      TmpGammaX = GammaX - EField * NbrSiteX * sin(Omega * t) / (Omega * M_PI);
      TmpGammaY = GammaY + Sign * NbrSiteY * EField * cos(Omega * t) / (Omega * M_PI);
      cout << "(gx, gy) = " << TmpGammaX << " " << TmpGammaY << endl;
      
      if (!TiltedFlag)
	TightBindingModel = new TightBindingModelCheckerboardLattice (NbrSiteX, NbrSiteY, t1, t2, tpp, muS, TmpGammaX, TmpGammaY, Architecture.GetArchitecture(), true, false);
      else
	TightBindingModel = new TightBindingModelCheckerboardLattice (NbrSiteX, NbrSiteY, nx1, ny1, nx2, ny2, Offset, t1, t2, tpp, 
							     muS, TmpGammaX, TmpGammaY, Architecture.GetArchitecture(), OffsetReal, true);
      TightBindingMatrix[i] = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
      delete TightBindingModel;
      
      if (i < NbrStoredHamiltonians)
      {
	cout << "build H for t = " << t << " + T" << endl;

	  if (TranslationFlag)
	  {
	    cout << NbrParticles << " " << NbrSites << " " << Kx << " " << NbrSiteX << " " << Ky << " " << NbrSiteY << endl;
	    
	    Hamiltonian[i] = new ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian (Space, NbrParticles, NbrSites, Kx, NbrSiteX, Ky, NbrSiteY, TightBindingMatrix[i], DensityDensityInteraction, Architecture.GetArchitecture(), Memory);
	  }
	  else
	  {
	    Hamiltonian[i] = new ParticleOnLatticeRealSpaceHamiltonian (Space, NbrParticles, NbrSites, 
									       TightBindingMatrix[i], DensityDensityInteraction,
									       Architecture.GetArchitecture(), Memory);
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

	  if (TranslationFlag)
	  {
	    Hamiltonian[NbrStoredHamiltonians] = new ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian (Space, NbrParticles, NbrSites, 
											       Kx, NbrSiteX, Ky, NbrSiteY,
											       TightBindingMatrix[TmpIndex], DensityDensityInteraction,
											       Architecture.GetArchitecture(), Memory);
	  }
	  else
	  {
	    Hamiltonian[NbrStoredHamiltonians] = new ParticleOnLatticeRealSpaceHamiltonian (Space, NbrParticles, NbrSites, 
									       TightBindingMatrix[TmpIndex], DensityDensityInteraction,
									       Architecture.GetArchitecture(), Memory);
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
// 	  if (i == 0)
// 	  {
// 	    HermitianMatrix HamiltonianMatrix (Space->GetHilbertSpaceDimension(), true);
// 	    Hamiltonian[i]->GetHamiltonian(HamiltonianMatrix);
// 	    cout << i << " " <<  HamiltonianMatrix << endl;
// 	  }
// 	  
	  VectorHamiltonianMultiplyOperation Operation1 (Hamiltonian[0], (&TmpInitialState), (&TmpState));
	  Operation1.ApplyOperation(Architecture.GetArchitecture());
	  TmpEnergy = TmpInitialState * TmpState;
	  FileEnergy0 << (TimeStep * i) << " " << (TmpEnergy.Re) << endl;
	  cout << "E_0 = " << (TmpEnergy.Re) << endl;
	  
	}
      
      TmpOverlap = InputState * TmpInitialState;
      cout << TmpOverlap <<  " " << (sqrt(TmpOverlap.Re*TmpOverlap.Re + TmpOverlap.Im*TmpOverlap.Im)) << endl;
      FileNorm << (TimeStep * i) << " " << (sqrt(TmpOverlap.Re*TmpOverlap.Re + TmpOverlap.Im*TmpOverlap.Im)) << endl;
      
      
      
      if (TmpIndex == NbrStoredHamiltonians)
	delete  Hamiltonian[NbrStoredHamiltonians];
    }  
    
   
   FileNorm.close();
   if (Manager.GetBoolean("compute-energy"))
   {
    FileEnergy.close();
    FileEnergy0.close();
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


// try to guess system information from file name
//
// filename = vector file name
// tauI = reference to the value of tau
// tfI = 
// nbrstepI = number of steps in previous simulation
// stepI = reference to the number of steps already done
// return value = true if no error occured

bool FTITimeEvolutionFindSystemInfoFromVectorFileName (char* filename, double& eFieldI, double& omegaI, int& nbrSamplesI,int& nbrPeriodsI, int& stepI)
{
//   char* StrTau;
//   StrTau = strstr(filename, "_tau_");
//   if (StrTau != 0)
//     {
//       StrTau += 5;
//       int SizeString = 0;
//       while ((StrTau[SizeString] != '\0') && (StrTau[SizeString] != '_') && (StrTau[SizeString] <= '9'))
// 	++SizeString;
//       if ((StrTau[SizeString] == '_') && (SizeString != 0))
// 	{
//           char TmpChar = StrTau[SizeString];
// 	  StrTau[SizeString] = '\0';
// 	  tauI = atoi(StrTau);
// 	  StrTau[SizeString] = TmpChar;
// 	  StrTau += SizeString;
// 	}
//       else
// 	StrTau = 0;
//     }
//   if (StrTau == 0)
//     {
//       return false;            
//     }
//   StrTau = strstr(filename, "_tf_");
//   if (StrTau != 0)
//     {
//       StrTau += 4;
//       int SizeString = 0;
//       while ((StrTau[SizeString] != '\0') && (StrTau[SizeString] != '_') && (StrTau[SizeString] <= '9'))
// 	++SizeString;
//       if ((StrTau[SizeString] == '_') && (SizeString != 0))
// 	{
//           char TmpChar = StrTau[SizeString];
// 	  StrTau[SizeString] = '\0';
// 	  tfI = atoi(StrTau);
// // 	  cout << "tfi = " << tfI << endl;
// 	  StrTau[SizeString] = TmpChar;
// 	  StrTau += SizeString;
// 	}
//       else
// 	StrTau = 0;
//     }
//   if (StrTau == 0)
//     {
//       return false;
//     }
//   StrTau = strstr(filename, "_nstep_");
//   if (StrTau != 0)
//     {
//       StrTau += 7;
//       int SizeString = 0;
//       while ((StrTau[SizeString] != '\0') && (StrTau[SizeString] != '.') && (StrTau[SizeString] >= '0') 
// 	     && (StrTau[SizeString] <= '9'))
// 	++SizeString;
//       if ((StrTau[SizeString] == '.') && (SizeString != 0))
// 	{
//           char TmpChar = StrTau[SizeString];
// 	  StrTau[SizeString] = '\0';
// 	  nbrStepI = atoi(StrTau);
// // 	  cout << "NbrStepI = " << nbrStepI << endl;
// 	  StrTau[SizeString] = TmpChar;
// 	  StrTau += SizeString;
// 	}
//       else
// 	StrTau = 0;
//       
//       StrTau += 1;
//       SizeString = 0;
//       while ((StrTau[SizeString] != '\0') && (StrTau[SizeString] != '.') 
// 	     && (StrTau[SizeString] <= '9'))
// 	++SizeString;
//       if ((StrTau[SizeString] == '.') && (SizeString != 0))
// 	{
//           char TmpChar = StrTau[SizeString];
// 	  StrTau[SizeString] = '\0';
// 	  stepI = atoi(StrTau);
// // 	  cout << "stepI  = " << stepI << endl;
// 	  StrTau[SizeString] = TmpChar;
// 	  StrTau += SizeString;
// 	}
//       else
// 	StrTau = 0;
//     }
//   if (StrTau == 0)
//     {
//       return false;
//     }
//   
//   
//   StrTau = strstr(filename, "_uf_");
//   if (StrTau != 0)
//     {
//       StrTau += 4;
//       int SizeString = 0;
//       while ((StrTau[SizeString] != '\0') && (StrTau[SizeString] != '_') && (StrTau[SizeString] <= '9'))
// 	++SizeString;
//       if ((StrTau[SizeString] == '_') && (SizeString != 0))
// 	{
//           char TmpChar = StrTau[SizeString];
// 	  StrTau[SizeString] = '\0';
// 	  uFinalI = atoi(StrTau);
// // 	  cout << "tau = " << tauI << endl;
// 	  StrTau[SizeString] = TmpChar;
// 	  StrTau += SizeString;
// 	}
//       else
// 	StrTau = 0;
//     }
  
  return true;
}
