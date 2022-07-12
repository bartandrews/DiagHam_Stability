#include "Options/Options.h"

#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelSimpleSquareLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"

#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/FactorialCoefficient.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/LongRational.h"


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



// compute the z_0 value obtained from the saddle point approximation
//
// oneBodyEntanglementTrimmedEnergies= array containing the one-body entanglement energies
// nbrOneBodyEntanglementTrimmedEnergies = number of one-body entanglement energies
// nbrParticlesA = number of particles in the subsystem A
double GetZ0Value (double* oneBodyEntanglementTrimmedEnergies, int nbrOneBodyEntanglementTrimmedEnergies, 
		   int nbrParticlesA);

// compute the z_0 value obtained from the saddle point approximation for a generic Renyi entropy
//
// oneBodyEntanglementTrimmedEnergies= array containing the one-body entanglement energies
// nbrOneBodyEntanglementTrimmedEnergies = number of one-body entanglement energies
// nbrParticlesA = number of particles in the subsystem A
// renyiIndex = Renyi index
double GetZ0Value (double* oneBodyEntanglementTrimmedEnergies, int nbrOneBodyEntanglementTrimmedEnergies, 
		   int nbrParticlesA, double renyiIndex);

double GetZ0DefintionSum (double* oneBodyEntanglementTrimmedEnergies, int nbrOneBodyEntanglementTrimmedEnergies, 
			  int nbrParticlesA, double z0);

double GetZ0DefintionSum (double* oneBodyEntanglementTrimmedEnergies, int nbrOneBodyEntanglementTrimmedEnergies, 
			  int nbrParticlesA, double z0, double renyiIndex);


// extract the correlation matrix for a given region out of the correlation matrix for a bigger region
//
// correlationMatrix = correlation matrix for the bigger region
// sourceNbrSitesX = number of sites along x for the bigger region
// sourceNbrSitesY = number of sites along y for the bigger region
// targetNbrSitesX = number of sites along x for the smaller region
// targetNbrSitesY = number of sites along y for the smaller region
HermitianMatrix EtaPairaingEntanglementEntropyExtractCorrelationMatrix(HermitianMatrix& correlationMatrix, int sourceNbrSitesX, int sourceNbrSitesY, 
								       int targetNbrSitesX, int targetNbrSitesY);


// compute the contribution to the Renyi entropy for a given number of particles in the subregion A, using the exact evaluation of the partitions
//
// oneBodyEntanglementTrimmedEnergies= array containing the one-body entanglement energies
// nbrOneBodyEntanglementTrimmedEnergies = number of one-body entanglement energies
// nbrParticlesA = number of particles in the subsystem A
// currentOrbitalIndex = current orbital that is considered
// currentFactor = current factor for a single partition
// entropies = array that contains the Renyi entropies
// maxRenyiEntropy = maximum Renyi entropy that has to be evaluated
// alpha = reference total weight of the reduced density matrix for the sector with nbrParticlesA particles
void GetEntanglementEntropyPerNbrParticlesA(double* oneBodyEntanglementTrimmedEnergies, int nbrOneBodyEntanglementTrimmedEnergies, 
					    int nbrParticlesA, int currentOrbitalIndex, double currentFactor, double* entropies, int maxRenyiEntropy, double& alpha);


void ComputeThermalQuantities (double beta, double mu, int nbrStates, double* stateEnergies, double& thermalNbrParticles, double& thermalEnergy, double* thermalEntropy, int maxRenyiEntropy,
			       double& thermalNbrParticlesMuDerivative, double& thermalNbrParticlesBetaDerivative, 
			       double& thermalEnergyMuDerivative, double& thermalEnergyBetaDerivative);

// evaluate the contribution of the eta pairing to the Renyi entropies
// 
// nbrRenyiEntropies = number of Renyi Entropies to evaluate
// nbrSites = total number of sites 
// nbrPairs = number of eta pairing pairs
// vacuumNbrParticles = number of particles for the vacuum states 
// totalNbrSitesA = otal number of sites in the region A
// nbrParticlesA = number of particles in the region A
// useRational = true if rational numbers have to be used for intermediate calculations
// rationalCoefficient = reference on the temporary rational coefficient
// binomial = reference on the binomial coefficients
// return value = array containing the eta pairing contribution for each Renyi entropy
double* EvaluateEtaPairingContribution(int nbrRenyiEntropies, int nbrSites, int nbrPairs, int vacuumNbrParticles, int totalNbrSitesA, int nbrParticlesA,
				       bool useRational, LongRational& rationalCoefficient, BinomialCoefficients& binomial);



int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("HubbardSquareLatticeEtaPairingEntanglementEntropyRealSpacePartition" , "0.01");
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

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-pairs", "number of pairs", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of unit cells along the x direction", 4);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of unit cells along the y direction", 4);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrsitex-a", "number of unit cells along the x direction for the part A", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrsitey-a", "number of unit cells along the y direction for the part A", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-nbrsitexa", "maximum number of unit cells along the x direction for the part A (equal to --nbrsitex-a if negative)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-nbrsiteya", "maximum number of unit cells along the y direction for the part A (equal to --nbrsitey-a if negative)", -1);
  (*SystemGroup) += new SingleStringOption  ('\n', "cuts", "provide the description of all cuts as a two column formatted text file");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nearbyeta-x", "x distance of the broken pair when generating a nearby eta pairing state", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nearbyeta-y", "y distance of the broken pair when generating a nearby eta pairing state", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "use-nonvacuum", "apply the eta^+ operators to a non-vacuum state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-particles", "number of particles for the non-vacuum state", 0);
  (*SystemGroup) += new BooleanOption ('\n', "use-fermisea", "apply the eta^+ operators to the Fermi sea");
  (*SystemGroup) += new SingleStringOption  ('\n', "nonvacuum-file", "provide the description of the non-vacuum state as a two-column ASCII file");
  (*SystemGroup) += new BooleanOption ('\n', "use-random", "generate a random non-vacuum state");
  (*SystemGroup) += new SingleIntegerOption ('\n', "run-id", "add an additional run id to the file name when using the --use-random option", 0);  
  (*SystemGroup) += new BooleanOption ('\n', "randomize-nonvacuumfile", "randomize the non-vacuum state file by randomly changing occupied states while preserving the total energy");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "randomize-nbrmoves", "number of particles that should be moved before testing is the total energy has changed within the error bar", 10);
  (*SystemGroup) += new SingleDoubleOption ('\n', "randomize-energyerror", "acceptable error on the total energy when randomizing the non-vacuum state file by a given number of moves", 0.001);
  (*SystemGroup) += new SingleIntegerOption ('\n', "randomize-maxnbrgroupmoves", "maximum number of moves of --randomize-nbrmoves particles", 10);
  (*SystemGroup) += new BooleanOption  ('\n', "show-nonvacuum", "show the non-vacuum state in the momentum basis");
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");  
  (*SystemGroup) += new BooleanOption ('\n', "use-approximation", "use a saddle appoximation to evaluate the entanglement entropy");
  (*SystemGroup) += new BooleanOption ('\n', "use-rational", "use rational number to overcome accuracy issues");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-renyi", "maximum Renyi entropy to compute", 1);
  (*SystemGroup) += new BooleanOption ('\n', "test-thermalentropy", "check if the state satisfies the ETH comparing the entropies");
  (*SystemGroup) += new BooleanOption ('\n', "test-thermalcorrelation", "check if the state satisfies the ETH comparing the correlation matrices");
  (*SystemGroup) += new SingleDoubleOption ('\n', "thermal-beta", "inverse temperature value to use for thermal calculations", -0.0020268);
  (*SystemGroup) += new SingleDoubleOption ('\n', "thermal-mu", "chemical potential value to use for thermal calculations", 542.045);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardSquareLatticeEtaPairingEntanglementEntropyRealSpacePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  bool ShowTimeFlag = Manager.GetBoolean("show-time");
  int NbrPairs = Manager.GetInteger("nbr-pairs"); 

  int XMomentum = 0;
  int YMomentum = 0;
  int NbrSitesX = 0;
  int NbrSitesY = 0;
  int NbrSites = 0; 
  int TotalSz = 0;
  bool Statistics = true;
  int NbrRenyiEntropies = Manager.GetInteger("max-renyi"); 

  NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  NbrSitesY = Manager.GetInteger("nbr-sitey"); 
  NbrSites = NbrSitesX * NbrSitesY; 

  if ((NbrSitesX & 1) != 0)
    {
      cout << "error, eta pairing states require an even number of sites in the x direction" << endl;
      return 0;
    }
  if ((NbrSitesY & 1) != 0)
    {
      cout << "error, eta pairing states require an even number of sites in the y direction" << endl;
      return 0;
    }

  if ((NbrPairs & 1) != 0)
    {
      XMomentum = NbrSitesX >> 1;
      YMomentum = NbrSitesY >> 1;          
    }
  int NbrSitesXA = Manager.GetInteger("nbrsitex-a"); 
  int NbrSitesYA = Manager.GetInteger("nbrsitey-a"); 
  int NbrSitesA = NbrSitesXA * NbrSitesYA;
  int NbrCuts = 1;
  int* CutX = 0;
  int* CutY = 0;
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
 
  Abstract2DTightBindingModel* TightBindingModel;
  TightBindingModel = new TightBindingModelSimpleSquareLattice (NbrSitesX, NbrSitesY, 1.0, 0.0, 0.0, 0.0,
								Architecture.GetArchitecture(), true);

  double* TightBindingModelEnergies = 0;
  int* TightBindingModelLinearizedMomenta = 0;
  TightBindingModel->GetEnergies(TightBindingModelEnergies, TightBindingModelLinearizedMomenta, 0);

  int* VacuumOneBodyLinearizedMomenta = 0; 
  double VacuumTotalEnergy = 0.0;
  int VacuumXMomentum = 0;
  int VacuumYMomentum = 0;
  int TmpMomentumX;
  int TmpMomentumY;
  int VacuumNbrParticles = Manager.GetInteger("nbr-particles"); 
  int VacuumTotalSz = VacuumNbrParticles;


  if (Manager.GetBoolean("use-fermisea") == true)
    { 
      VacuumOneBodyLinearizedMomenta = new int[VacuumNbrParticles];
      for (int i = 0; i < VacuumNbrParticles; ++i)
	{
	  TightBindingModel->GetLinearizedMomentumIndex(TightBindingModelLinearizedMomenta[i], TmpMomentumX, TmpMomentumY);
	  VacuumXMomentum += TmpMomentumX;
	  VacuumYMomentum += TmpMomentumY;
	  VacuumOneBodyLinearizedMomenta[i] = TightBindingModel->GetLinearizedMomentumIndex(TmpMomentumX, TmpMomentumY);
	  VacuumTotalEnergy += TightBindingModel->GetEnergy(0, VacuumOneBodyLinearizedMomenta[i]);
	}
    }
  else
    {
     if (Manager.GetString("nonvacuum-file") != 0)
       {
	 MultiColumnASCIIFile NonVacuumFile;
	 if (NonVacuumFile.Parse(Manager.GetString("nonvacuum-file")) == false)
	   {
	     NonVacuumFile.DumpErrors(cout);
	     return -1;
	   }
	 VacuumNbrParticles = NonVacuumFile.GetNbrLines();
	 VacuumTotalSz = VacuumNbrParticles;
	 VacuumOneBodyLinearizedMomenta = new int[VacuumNbrParticles];
	 int* TmpXMomenta = NonVacuumFile.GetAsIntegerArray(0);
	 int* TmpYMomenta = NonVacuumFile.GetAsIntegerArray(1);
	 for (int i = 0; i < VacuumNbrParticles; ++i)
	   {
	     VacuumOneBodyLinearizedMomenta[i] = TightBindingModel->GetLinearizedMomentumIndex(TmpXMomenta[i], TmpYMomenta[i]);
	     VacuumTotalEnergy += TightBindingModel->GetEnergy(0, VacuumOneBodyLinearizedMomenta[i]);
	     VacuumXMomentum += TmpXMomenta[i];
	     VacuumYMomentum += TmpYMomenta[i];	     
	   }
	 if (Manager.GetBoolean("randomize-nonvacuumfile") == true)
	   {
	     int NbrDiscardedTightBindingModelEnergies = NbrSites - VacuumNbrParticles;
	     int* DiscardedLinearizedMomenta = new int[NbrDiscardedTightBindingModelEnergies];
	     NbrDiscardedTightBindingModelEnergies = 0;
	     int* ReferenceVacuumOneBodyLinearizedMomenta = new int [VacuumNbrParticles];
	     for (int i = 0; i < VacuumNbrParticles; ++i)
	       {
		 ReferenceVacuumOneBodyLinearizedMomenta[i] = VacuumOneBodyLinearizedMomenta[i];
	       }
	     for (int i = 0; i < NbrSites; ++i)
	       {
		 if (SearchInUnsortedArray(TightBindingModelLinearizedMomenta[i], VacuumOneBodyLinearizedMomenta, VacuumNbrParticles) == -1)
		   {
		     DiscardedLinearizedMomenta[NbrDiscardedTightBindingModelEnergies] = TightBindingModelLinearizedMomenta[i]; 
		     ++NbrDiscardedTightBindingModelEnergies;
		   }
	       }
	     AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (0);
	     RandomNumber->UseTimeSeed();

	     int NbrMoves = Manager.GetInteger("randomize-maxnbrgroupmoves");

	     long NbrRejetedMoves = 0;
	     long NbrAcceptedMoves = 0;
	     double DemonEnergy = 0.0;
	     while (NbrMoves > 0)
	       {
		 NbrRejetedMoves = 0;
		 NbrAcceptedMoves = 0;
		 for (int i = 0; i < VacuumNbrParticles; ++i)
		   {
		     int DestinationIndex = (int) (RandomNumber->GetRealRandomNumber() * ((double) NbrDiscardedTightBindingModelEnergies));
		     double TmpEnergyDifference = (TightBindingModel->GetEnergy(0, DiscardedLinearizedMomenta[DestinationIndex]) 
						   - TightBindingModel->GetEnergy(0, VacuumOneBodyLinearizedMomenta[i]));
		     if (TmpEnergyDifference < DemonEnergy)
		       {
			 DemonEnergy -= TmpEnergyDifference;
			 int TmpMomentum = VacuumOneBodyLinearizedMomenta[i];
			 VacuumOneBodyLinearizedMomenta[i] = DiscardedLinearizedMomenta[DestinationIndex];
			 DiscardedLinearizedMomenta[DestinationIndex] = TmpMomentum;
			 ++NbrAcceptedMoves;
		       }
		     else
		       {
			 ++NbrRejetedMoves;
		       }	
		   }
		 HermitianMatrix EntanglementHamiltonian = TightBindingModel->EvaluateFullTwoPointCorrelationFunction(MaxNbrSitesXA, MaxNbrSitesYA, VacuumOneBodyLinearizedMomenta, VacuumNbrParticles, 0);
		 Complex Tmp;
		 EntanglementHamiltonian.GetMatrixElement(1, 0, Tmp);
		 cout << Tmp.Re << " " << Tmp.Im << " " << DemonEnergy << endl;		 
// 		 for (int i = 0; i < VacuumNbrParticles; ++i)
// 		   {
// 		     VacuumOneBodyLinearizedMomenta[i] = ReferenceVacuumOneBodyLinearizedMomenta[i];
// 		   }
	       }
// 	     while (NbrAcceptedMoves < NbrMoves)
// 	       {
// 		 int SourceIndex = (int) (RandomNumber->GetRealRandomNumber() * ((double) VacuumNbrParticles));
// 		 int DestinationIndex = (int) (RandomNumber->GetRealRandomNumber() * ((double) NbrDiscardedTightBindingModelEnergies));
// 		 double TmpEnergyDifference = (TightBindingModel->GetEnergy(0, DiscardedLinearizedMomenta[DestinationIndex]) 
// 					       - TightBindingModel->GetEnergy(0, VacuumOneBodyLinearizedMomenta[SourceIndex]));
// 		 if (TmpEnergyDifference < DemonEnergy)
// 		   {
// 		     DemonEnergy -= TmpEnergyDifference;
//  		     int TmpMomentum = VacuumOneBodyLinearizedMomenta[SourceIndex];
//  		     VacuumOneBodyLinearizedMomenta[SourceIndex] = DiscardedLinearizedMomenta[DestinationIndex];
//  		     DiscardedLinearizedMomenta[DestinationIndex] = TmpMomentum;
// 		     ++NbrAcceptedMoves;
// 		     HermitianMatrix EntanglementHamiltonian = TightBindingModel->EvaluateFullTwoPointCorrelationFunction(MaxNbrSitesXA, MaxNbrSitesYA, VacuumOneBodyLinearizedMomenta, VacuumNbrParticles, 0);
// 		   }
// 		 else
// 		   {
// 		     ++NbrRejetedMoves;
// 		   }
// 	       }
// 	     NbrAcceptedMoves = 0l;
// 	     NbrRejetedMoves = 0l;
// 	     while (NbrAcceptedMoves < (NbrMoves * 4) )
// 	       {
// 		 int SourceIndex = (int) (RandomNumber->GetRealRandomNumber() * ((double) VacuumNbrParticles));
// 		 int DestinationIndex = (int) (RandomNumber->GetRealRandomNumber() * ((double) NbrDiscardedTightBindingModelEnergies));
// 		 double TmpEnergyDifference = (TightBindingModel->GetEnergy(0, DiscardedLinearizedMomenta[DestinationIndex]) 
// 					       - TightBindingModel->GetEnergy(0, VacuumOneBodyLinearizedMomenta[SourceIndex]));
// 		 if (TmpEnergyDifference < DemonEnergy)
// 		   {
// 		     DemonEnergy -= TmpEnergyDifference;
//  		     int TmpMomentum = VacuumOneBodyLinearizedMomenta[SourceIndex];
//  		     VacuumOneBodyLinearizedMomenta[SourceIndex] = DiscardedLinearizedMomenta[DestinationIndex];
//  		     DiscardedLinearizedMomenta[DestinationIndex] = TmpMomentum;
// 		     ++NbrAcceptedMoves;
// 		     HermitianMatrix EntanglementHamiltonian = TightBindingModel->EvaluateFullTwoPointCorrelationFunction(MaxNbrSitesXA, MaxNbrSitesYA, VacuumOneBodyLinearizedMomenta, VacuumNbrParticles, 0);
// 		     Complex Tmp;
// 		     EntanglementHamiltonian.GetMatrixElement(0, 1, Tmp);
// 		     cout << Tmp.Re << " " << Tmp.Im << " " << DemonEnergy << " " << TmpEnergyDifference << endl;
// 		   }
// 		 else
// 		   {
// 		     ++NbrRejetedMoves;
// 		   }
//	       }
	     
//	     double VacuumEnergyError = Manager.GetDouble("randomize-energyerror");
//	     int NbrMoves = Manager.GetInteger("randomize-nbrmoves");
// 	     int* SourceStates = new int [NbrMoves];
// 	     int* DestinationStates = new int [NbrMoves];
// 	     double TmpEnergyDifference = 0.0;
// 	     int MaxNbrGroupMoves = (int) (RandomNumber->GetRealRandomNumber() * ((double) Manager.GetInteger("randomize-maxnbrgroupmoves")));
// 	     if (MaxNbrGroupMoves <= 0)
// 	       MaxNbrGroupMoves = 1;
// 	     long NbrRejetedMoves = 0;
// 	     long NbrAcceptedMoves = 0;
// 	     while (MaxNbrGroupMoves > 0)
// 	       {
// 		 for (int i = 0; i < NbrMoves; ++i)
// 		   {
// 		     SourceStates[i] = (int) (RandomNumber->GetRealRandomNumber() * ((double) VacuumNbrParticles));
// 		     DestinationStates[i] = (int) (RandomNumber->GetRealRandomNumber() * ((double) NbrDiscardedTightBindingModelEnergies));
// 		   }
// 		 double TmpCurrentEnergyDifference = 0.0;
// 		 for (int i = 0; i < NbrMoves; ++i)
// 		   {
// 		     TmpCurrentEnergyDifference += (TightBindingModel->GetEnergy(0, DiscardedLinearizedMomenta[DestinationStates[i]]) 
// 						    - TightBindingModel->GetEnergy(0, VacuumOneBodyLinearizedMomenta[SourceStates[i]]));
// 		     int TmpMomentum = VacuumOneBodyLinearizedMomenta[SourceStates[i]];
// 		     VacuumOneBodyLinearizedMomenta[SourceStates[i]] = DiscardedLinearizedMomenta[DestinationStates[i]];
// 		     DiscardedLinearizedMomenta[DestinationStates[i]] = TmpMomentum;
// 		   }
// 		 if (fabs(TmpEnergyDifference + TmpCurrentEnergyDifference) < VacuumEnergyError)
// 		   {
// 		     TmpEnergyDifference += TmpCurrentEnergyDifference;
// 		     --MaxNbrGroupMoves;
// 		     VacuumXMomentum = 0;
// 		     VacuumYMomentum = 0;
// 		     VacuumTotalEnergy = 0.0;
// 		     int TmpXMomentum;
// 		     int TmpYMomentum;
// 		     for (int i = 0; i < VacuumNbrParticles; ++i)
// 		       {
// 			 TightBindingModel->GetLinearizedMomentumIndex(VacuumOneBodyLinearizedMomenta[i], TmpXMomentum, TmpYMomentum);
// 			 VacuumTotalEnergy += TightBindingModel->GetEnergy(0, VacuumOneBodyLinearizedMomenta[i]);
// 			 VacuumXMomentum += TmpXMomentum;
// 			 VacuumYMomentum += TmpYMomentum;	     
// 		       }
// 		     ++NbrAcceptedMoves;
// 		   }
// 		 else
// 		   {
// 		     for (int i = NbrMoves - 1; i >= 0; --i)
// 		       {
// 			 int TmpMomentum = VacuumOneBodyLinearizedMomenta[SourceStates[i]];
// 			 VacuumOneBodyLinearizedMomenta[SourceStates[i]] = DiscardedLinearizedMomenta[DestinationStates[i]];
// 			 DiscardedLinearizedMomenta[DestinationStates[i]] = TmpMomentum;
// 		       }
// 		     ++NbrRejetedMoves;
// 		   }
// 	       }
	     int NbrNewStates = 0;
	     SortArrayUpOrdering(ReferenceVacuumOneBodyLinearizedMomenta, VacuumNbrParticles);
	     for (int i = 0; i < VacuumNbrParticles; ++i)
	       {
		 if (SearchInArray(VacuumOneBodyLinearizedMomenta[i], ReferenceVacuumOneBodyLinearizedMomenta, VacuumNbrParticles) == -1)
		   {
		     ++NbrNewStates;
		   }
	       }
//	     cout << "difference of energy for the randomized states = " << TmpEnergyDifference << endl;
	     cout << "number of accepted moves = " << NbrAcceptedMoves << endl;
	     cout << "number of rejected moves = " << NbrRejetedMoves << endl;
	     cout << "number of updated occupied one-body states = " << NbrNewStates << endl;
	     delete[] ReferenceVacuumOneBodyLinearizedMomenta;
	     delete[] DiscardedLinearizedMomenta;
       }
       }
     else
       {
	 if (Manager.GetBoolean("use-random") == true)
	   {
	     VacuumOneBodyLinearizedMomenta = new int[VacuumNbrParticles];
	     AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (0);
	     RandomNumber->UseTimeSeed();
	     for (int i = 0; i < VacuumNbrParticles; ++i)
	       {
		 VacuumOneBodyLinearizedMomenta[i] = TightBindingModelLinearizedMomenta[i];
	       }
	     for (int i = VacuumNbrParticles; i < NbrSites; ++i)
	       {
		 int Tmp = (int) (RandomNumber->GetRealRandomNumber() * (((double) i) + 0.0001));
		 if (Tmp < VacuumNbrParticles)
		   VacuumOneBodyLinearizedMomenta[Tmp] = TightBindingModelLinearizedMomenta[i];
	       }
	     for (int i = 0; i < VacuumNbrParticles; ++i)
	       {
		 VacuumTotalEnergy += TightBindingModel->GetEnergy(0, VacuumOneBodyLinearizedMomenta[i]);
		 TightBindingModel->GetLinearizedMomentumIndex(VacuumOneBodyLinearizedMomenta[i], TmpMomentumX, TmpMomentumY);
		 VacuumXMomentum += TmpMomentumX;
		 VacuumYMomentum += TmpMomentumY;	     
	       }
	   }
       }
    }

  if (Manager.GetBoolean("show-nonvacuum") == true)
    {
      cout << "using states : " << endl;
      for (int i = 0; i < VacuumNbrParticles; ++i)
	{
	  TightBindingModel->GetLinearizedMomentumIndex(VacuumOneBodyLinearizedMomenta[i], TmpMomentumX, TmpMomentumY);
	  cout << "(" << TmpMomentumX << ", " << TmpMomentumY << ")" << endl;
	}
    }
  
  XMomentum += VacuumXMomentum;
  YMomentum += VacuumYMomentum;      
  VacuumXMomentum %= NbrSitesX;
  VacuumYMomentum %= NbrSitesY;
  XMomentum %= NbrSitesX;
  YMomentum %= NbrSitesY;

  int NbrParticles = VacuumNbrParticles; 
  NbrParticles += 2 * NbrPairs;
  TotalSz += VacuumTotalSz;
  char* StatisticPrefix = new char [64];
  sprintf (StatisticPrefix, "fermions_hubbard");
  char* FilePrefix = new char [256];
  if ((Manager.GetInteger("nearbyeta-x") == 0) && (Manager.GetInteger("nearbyeta-y") == 0))
    {
      sprintf (FilePrefix, "%s_square_etapairing_nbrpairs_%ld_x_%d_y_%d_n_%d_ns_%d", StatisticPrefix, Manager.GetInteger("nbr-pairs"), NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
    }
  else
    {
      sprintf (FilePrefix, "%s_square_nearbyetapairing_nbrpairs_%ld_alphax_%ld_alphay_%ld_x_%d_y_%d_n_%d_ns_%d", StatisticPrefix, Manager.GetInteger("nbr-pairs"), Manager.GetInteger("nearbyeta-x"), 
	       Manager.GetInteger("nearbyeta-y"), NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
    }

  if (Manager.GetBoolean("test-thermalentropy") == true)
    {
      double MinBeta = -0.002;
      double BetaStep = 0.0001;
      int NbrBetaSteps = 20;
      double MinMu = 500;//TightBindingModelEnergies[0] - fabs(TightBindingModelEnergies[0] * 2);
      double MaxMu = TightBindingModelEnergies[NbrSites - 1] + fabs(TightBindingModelEnergies[NbrSites - 1] * 2);
      int NbrMuSteps = 5;
      double MuStep = 20;//(MaxMu - MinMu) / ((double) NbrMuSteps);
      double* EnergyExponentials = new double[NbrSites];
      double** ThermalEnergies = new double*[NbrBetaSteps];
      double** ThermalNbrParticules = new double*[NbrBetaSteps];
      double** ThermalEntropy = new double*[NbrBetaSteps];
      for (int j = 0; j < NbrBetaSteps; ++j)
	{
	  ThermalEnergies[j] = new double[NbrMuSteps];
	  ThermalNbrParticules[j] = new double[NbrMuSteps];
	  ThermalEntropy[j] = new double[NbrMuSteps];
	}
      double MinErrorNbrParticules = (double) VacuumNbrParticles;
      double MinErrorEnergy = ((double) VacuumNbrParticles) * TightBindingModelEnergies[NbrSites - 1];
      double CurrentBeta = MinBeta;
      for (int j = 0; j < NbrBetaSteps; ++j)
	{
	  for (int k = 0 ; k < NbrSites; ++k)
	    {
	      EnergyExponentials[k] = exp(CurrentBeta * TightBindingModelEnergies[k]);
	    }
	  double CurrentMu = MinMu;
	  for (int i = 0; i < NbrMuSteps; ++i)
	    {
	      double MuFactor = exp(-CurrentBeta * CurrentMu);
	      double TmpEnergy = 0.0;
	      double TmpNbrParticules = 0.0;
	      double TmpThermalEntropy = 0.0;
	      for (int k = 0 ; k < NbrSites; ++k)
		{
		  double Tmp = 1.0 / (1.0 + (EnergyExponentials[k] * MuFactor));
		  double Tmp2 = EnergyExponentials[k] * MuFactor * Tmp;
		  TmpNbrParticules += Tmp;
		  TmpEnergy += Tmp * TightBindingModelEnergies[k];
		  TmpThermalEntropy -= (Tmp * log (Tmp)) + (Tmp2 * log (Tmp2));
		}
	      ThermalEnergies[j][i] = TmpEnergy;
	      ThermalNbrParticules[j][i] = TmpNbrParticules;
	      ThermalEntropy[j][i] = TmpThermalEntropy;
	      if (abs(TmpNbrParticules - VacuumNbrParticles) < MinErrorNbrParticules)
		{
		  MinErrorNbrParticules = abs(TmpNbrParticules - VacuumNbrParticles);
		}
	      if (fabs(ThermalEnergies[j][i] - VacuumTotalEnergy) < MinErrorEnergy)
		{
		  MinErrorEnergy = fabs(ThermalEnergies[j][i] - VacuumTotalEnergy);
		}
//  	      if (ThermalEnergies[j][i] > 0.0)
// 		cout << CurrentBeta << " " << CurrentMu << " " << ThermalNbrParticules[j][i] << " " << ThermalEnergies[j][i] << " " << ThermalEntropy[j][i] << endl;
	      CurrentMu += MuStep;
	    }
	  CurrentBeta += BetaStep;
	}
      cout << "min energy error = " << MinErrorEnergy << endl;
      cout << "min nbr particle error = " << MinErrorNbrParticules << endl;
      CurrentBeta = MinBeta;
      MinErrorNbrParticules += ((double) VacuumNbrParticles) * 0.05;
      cout << "best thermal nbr of particle approximations : " << endl;
      int NbrAcceptedValues = 0;
      MinErrorEnergy = ((double) VacuumNbrParticles) * TightBindingModelEnergies[NbrSites - 1];
      for (int j = 0; j < NbrBetaSteps; ++j)
	{
	  double CurrentMu = MinMu;
	  for (int i = 0; i < NbrMuSteps; ++i)
	    {
	      if (abs(ThermalNbrParticules[j][i] - VacuumNbrParticles) < MinErrorNbrParticules)
		{
 		  if (fabs(ThermalEnergies[j][i] - VacuumTotalEnergy) < MinErrorEnergy)
 		    {
 		      MinErrorEnergy = fabs(ThermalEnergies[j][i] - VacuumTotalEnergy);
 		    }
		  //		  cout << CurrentBeta << " " << CurrentMu << " " << ThermalNbrParticules[j][i] << " " << ThermalEnergies[j][i] << " " << ThermalEntropy[j][i] << endl;
		}
	      CurrentMu += MuStep;
	    }
	  CurrentBeta += BetaStep;
	}
      
      CurrentBeta = MinBeta;
      MinErrorEnergy += fabs(((double) VacuumTotalEnergy) * 0.1);
      cout << "best thermal energy approximations : " << endl;
      for (int j = 0; j < NbrBetaSteps; ++j)
	{
	  double CurrentMu = MinMu;
	  for (int i = 0; i < NbrMuSteps; ++i)
	    {
	      if ((abs(ThermalNbrParticules[j][i] - VacuumNbrParticles) < MinErrorNbrParticules) && (fabs(ThermalEnergies[j][i] - VacuumTotalEnergy) < MinErrorEnergy))
		//	      if (fabs(ThermalEnergies[j][i] - VacuumTotalEnergy) < MinErrorEnergy)
		{
		  cout << CurrentBeta << " " << CurrentMu << " " << ThermalNbrParticules[j][i] << " " << ThermalEnergies[j][i] << " " << ThermalEntropy[j][i] << endl;
		  ++NbrAcceptedValues;
		}
	      CurrentMu += MuStep;
	    }
	  CurrentBeta += BetaStep;
	}
      cout << "nbr of best thermal energy approximations = " << NbrAcceptedValues << endl;
    }
  HermitianMatrix ThermalCorrelationMatrix;
  if (Manager.GetBoolean("test-thermalcorrelation") == true)
    {
      double Mu = Manager.GetDouble("thermal-mu");
      double Beta = Manager.GetDouble("thermal-beta");
      double ThermalNbrParticles;
      double ThermalEnergy;
      double* ThermalEntropy = new double[NbrRenyiEntropies];
      double ThermalNbrParticlesMuDerivative;
      double ThermalNbrParticlesBetaDerivative;
      double ThermalEnergyMuDerivative;
      double ThermalEnergyBetaDerivative;
      ComputeThermalQuantities (Beta, Mu, NbrSites, TightBindingModelEnergies, 
				ThermalNbrParticles, ThermalEnergy, ThermalEntropy, NbrRenyiEntropies,
				ThermalNbrParticlesMuDerivative, ThermalNbrParticlesBetaDerivative, 
				ThermalEnergyMuDerivative, ThermalEnergyBetaDerivative);
      cout << "thermal nbr of particles = " << ThermalNbrParticles << endl;
      cout << "thermal energy = " << ThermalEnergy << endl;
      cout << "thermal entropy = " << ThermalEntropy[0] << endl;
      cout << "thermal entropy per volume = " << (ThermalEntropy[0] / (((double) NbrSitesX) * ((double) NbrSitesY)))<< endl;
      if (NbrRenyiEntropies > 1)
	{
	  for (int i = 1; i < NbrRenyiEntropies; ++i)
	    {
	      cout << "thermal Renyi entropy (n=" << (i + 1) << ") = " << ThermalEntropy[i] << endl;
	      cout << "thermal Renyi entropy (n=" << (i + 1) << ") per volume = " << (ThermalEntropy[i] / (((double) NbrSitesX) * ((double) NbrSitesY)))<< endl;
	    }
	}
      double* FermiFactors = new double[NbrSites];
      double Tmp = 0.0;
      for (int k = 0 ; k < NbrSites; ++k)
	{
	  FermiFactors[k] = 1.0 / (1.0 + exp(Beta * (TightBindingModelEnergies[k] - Mu)));
	  Tmp += FermiFactors[k];
	}
      Tmp /= ((double) NbrSites);
      int TmpNbrSitesA = NbrSitesXA * NbrSitesYA;
      int TmpMomentumX;
      int TmpMomentumY;
      ThermalCorrelationMatrix = HermitianMatrix(TmpNbrSitesA, true);
      Complex** PhaseFactors = new Complex*[2 * NbrSitesXA + 1];
      for (int i = 0; i <= (2 * NbrSitesXA); ++i)
	{
	  PhaseFactors[i] = new Complex[2 * NbrSitesYA + 1];
	  double TmpXFactor = 2.0 * M_PI * ((double) (i - NbrSitesXA)) / ((double) NbrSitesX);
	  for (int j = 0; j <= (2 * NbrSitesYA); ++j)
	    {
	      double TmpYFactor = 2.0 * M_PI * ((double) (j - NbrSitesYA)) / ((double) NbrSitesY);
	      Complex Tmp2 = 0.0;
	      for (int k = 0 ; k < NbrSites; ++k)
		{
		  TightBindingModel->GetLinearizedMomentumIndex(TightBindingModelLinearizedMomenta[k], TmpMomentumX, TmpMomentumY);
		  Tmp2 += Phase((((double) TmpMomentumX) * TmpXFactor) + (((double) TmpMomentumY) * TmpYFactor)) * FermiFactors[k];
		}
	      PhaseFactors[i][j] = Tmp2 / ((double) NbrSites);
	    }
	}
      for (int i = 0 ; i < TmpNbrSitesA; ++i)
	{
	  ThermalCorrelationMatrix.SetMatrixElement(i, i, Tmp);
	  int TmpXA1 = i / NbrSitesYA;
	  int TmpYA1 = i % NbrSitesYA;
	  for (int j = i + 1 ; j < TmpNbrSitesA; ++j)
	    {
	      int TmpXA2 = j / NbrSitesYA;
	      int TmpYA2 = j % NbrSitesYA;
	      ThermalCorrelationMatrix.SetMatrixElement(i, j, PhaseFactors[(TmpXA1 - TmpXA2) + NbrSitesXA][(TmpYA1 - TmpYA2) + NbrSitesYA]);
	    }
	}
      delete[] FermiFactors;
      for (int i = 0; i <= (2 * NbrSitesXA); ++i)
	{
	  delete[] PhaseFactors[i];
	}
      delete[] PhaseFactors;
      char* TmpMatrixFileName = new char[512];
      if ((Manager.GetBoolean("use-random") == true) || (Manager.GetBoolean("randomize-nonvacuumfile") == true))
	{
	  sprintf(TmpMatrixFileName, "%s_sz_%d_xa_%d_ya_%d_runid_%ld.thermal.mat.txt", FilePrefix, TotalSz, 
		  NbrSitesXA, NbrSitesYA, Manager.GetInteger("run-id"));
	}
      else
	{
	  sprintf(TmpMatrixFileName, "%s_sz_%d_xa_%d_ya_%d.thermal.mat.txt", FilePrefix, TotalSz, 
		  NbrSitesXA, NbrSitesYA);
	}
      ThermalCorrelationMatrix.WriteAsciiMatrix(TmpMatrixFileName, true);
//       RealDiagonalMatrix ThermalCorrelationMatrixEigenvalues(ThermalCorrelationMatrix.GetNbrRow(), true);
// #ifdef __LAPACK__
//       ThermalCorrelationMatrix.LapackDiagonalize(ThermalCorrelationMatrixEigenvalues);
// #else
//       ThermalCorrelationMatrix.Diagonalize(ThermalCorrelationMatrixEigenvalues);
// #endif
//       for (int i = 0; i < ThermalCorrelationMatrixEigenvalues.GetNbrRow(); ++i)
// 	cout << ThermalCorrelationMatrixEigenvalues[i] << endl;
    }

  char* EntropyFileName = new char [512];
  if ((Manager.GetBoolean("use-random") == true) || (Manager.GetBoolean("randomize-nonvacuumfile") == true))
    {
      sprintf(EntropyFileName, "%s_sz_%d_runid_%ld_xa_%d_ya_%d.ent", FilePrefix, TotalSz, Manager.GetInteger("run-id"), 
	      MaxNbrSitesXA, MaxNbrSitesYA);
    }
  else
    {
      sprintf(EntropyFileName, "%s_sz_%d_xa_%d_ya_%d.ent", FilePrefix, TotalSz, MaxNbrSitesXA, MaxNbrSitesYA);
    }

  char* OneBodyEntropyFileName = new char [512];
  if ((Manager.GetBoolean("use-random") == true) || (Manager.GetBoolean("randomize-nonvacuumfile") == true))
    {
      sprintf(OneBodyEntropyFileName, "%s_sz_%d_runid_%ld_xa_%d_ya_%d.onebody.ent", FilePrefix, TotalSz, Manager.GetInteger("run-id"), 
	      MaxNbrSitesXA, MaxNbrSitesYA);
    }
  else
    {
      sprintf(OneBodyEntropyFileName, "%s_sz_%d_xa_%d_ya_%d.onebody.ent", FilePrefix, TotalSz, MaxNbrSitesXA, MaxNbrSitesYA);
    }

  ofstream File;
  File.open(EntropyFileName, ios::binary | ios::out);
  File.precision(14);

  ofstream OneBodyFile;
  OneBodyFile.open(OneBodyEntropyFileName, ios::binary | ios::out);
  OneBodyFile.precision(14);

  if (Manager.GetBoolean("use-nonvacuum") == false)
    {
      FactorialCoefficient TmpCoefficient;
      for (; NbrSitesXA <= MaxNbrSitesXA; ++NbrSitesXA)
	{
	  if ((Manager.GetInteger("max-nbrsitexa") > 0) && (Manager.GetInteger("max-nbrsiteya") > 0))
	    NbrSitesYA = NbrSitesXA;
	  cout << "computing entropy of a " << NbrSitesXA << "x" << NbrSitesYA << " patch" << endl;
	  int TotalNbrSitesA = NbrSitesXA * NbrSitesYA;
	  int MaxSumIndex = TotalNbrSitesA;
	  if (MaxSumIndex > NbrPairs)
	    MaxSumIndex = NbrPairs;
	  double Tmp = 0.0;
	  double* EntanglementEntropies = new double [NbrRenyiEntropies];
	  for (int k = 0; k < NbrRenyiEntropies; ++k)
	    EntanglementEntropies[k] = 0.0;
	  for (int j = 0; j <= MaxSumIndex; ++j)
	    {
	      TmpCoefficient.SetToOne();
	      TmpCoefficient.BinomialMultiply(TotalNbrSitesA, j);
	      TmpCoefficient.BinomialMultiply(NbrSites - TotalNbrSitesA, NbrPairs - j);
	      TmpCoefficient.BinomialDivide(NbrSites, NbrPairs);
	      double Tmp2 = TmpCoefficient.GetNumericalValue();
	      cout << TotalNbrSitesA << " " << NbrSites << " " << j << " " << NbrPairs << " : " << Tmp2 << endl;
	      EntanglementEntropies[0] -= Tmp2 * log (Tmp2);
	      for (int k = 1; k < NbrRenyiEntropies; ++k)
		EntanglementEntropies[k] += pow(Tmp2, (double) (k + 1));
	    }
	  File << NbrSitesXA << " " << NbrSitesYA << " " << EntanglementEntropies[0];
	  for (int k = 1; k < NbrRenyiEntropies; ++k)
	    {
	      EntanglementEntropies[k] = -log(EntanglementEntropies[k]) / ((double) k);	  
	      File << " " << EntanglementEntropies[k];
	    }
	  File << endl;
	}
      File.close();
      return 0;
    }


  if (Manager.GetString("nonvacuum-file") == 0)
    {
      char* NonVacuumFileName = new char [512];
      if ((Manager.GetBoolean("use-random") == true) || (Manager.GetBoolean("randomize-nonvacuumfile") == true))
	{
	  sprintf(NonVacuumFileName, "%s_sz_%d_runid_%ld_nonvacuum.dat", FilePrefix, TotalSz, Manager.GetInteger("run-id"));
	}
      else
	{
	  sprintf(NonVacuumFileName, "%s_sz_%d_nonvacuum.dat", FilePrefix, TotalSz);
	}
      ofstream NonVacuumFile;
      NonVacuumFile.open(NonVacuumFileName, ios::binary | ios::out);
      NonVacuumFile.precision(14);
      NonVacuumFile << "# Non-vacuum state total momentum along x = " << VacuumXMomentum << endl;
      NonVacuumFile << "# Non-vacuum state total momentum along y = " << VacuumYMomentum << endl;
      NonVacuumFile << "# Non-vacuum state total energy = " << VacuumTotalEnergy << endl;
      for (int i = 0; i < VacuumNbrParticles; ++i)
	{
	  TightBindingModel->GetLinearizedMomentumIndex(VacuumOneBodyLinearizedMomenta[i], TmpMomentumX, TmpMomentumY);
	  NonVacuumFile << TmpMomentumX << " " << TmpMomentumY << endl;
	}
      NonVacuumFile.close();
    }

  File << "# Non-vacuum state total momentum along x = " << VacuumXMomentum << endl;
  File << "# Non-vacuum state total momentum along y = " << VacuumYMomentum << endl;
  File << "# Non-vacuum state total energy = " << VacuumTotalEnergy << endl;
  cout << "Non-vacuum state total momentum along x = " << VacuumXMomentum << endl;
  cout << "Non-vacuum state total momentum along y = " << VacuumYMomentum << endl;
  cout << "Non-vacuum state total energy = " << VacuumTotalEnergy << endl;

  timeval TotalStartingTime;
  timeval TotalEndingTime;
  if (ShowTimeFlag == true)
    {
      gettimeofday (&(TotalStartingTime), 0);
    }
  HermitianMatrix EntanglementHamiltonian = TightBindingModel->EvaluateFullTwoPointCorrelationFunction(MaxNbrSitesXA, MaxNbrSitesYA, VacuumOneBodyLinearizedMomenta, VacuumNbrParticles, 0);
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
      int TotalNbrSitesA = NbrSitesXA * NbrSitesYA;
      if (ShowTimeFlag == true)
	{
	  gettimeofday (&(TotalStartingTime), 0);
	}
      RealDiagonalMatrix VacuumOneBodyEntanglementEnergies(TotalNbrSitesA, true);
      HermitianMatrix TmpEntanglementHamiltonian = EtaPairaingEntanglementEntropyExtractCorrelationMatrix(EntanglementHamiltonian, MaxNbrSitesXA, MaxNbrSitesYA, NbrSitesXA, NbrSitesYA);
       if (Manager.GetBoolean("test-thermalcorrelation") == true)
 	{
	  char* TmpMatrixFileName = new char[512];
	  if ((Manager.GetBoolean("use-random") == true) || (Manager.GetBoolean("randomize-nonvacuumfile") == true))
	    {
	      sprintf(TmpMatrixFileName, "%s_sz_%d_xa_%d_ya_%d_runid_%ld.corr.mat.txt", FilePrefix, TotalSz, 
		      NbrSitesXA, NbrSitesYA, Manager.GetInteger("run-id"));
	    }
	  else
	    {
	      sprintf(TmpMatrixFileName, "%s_sz_%d_xa_%d_ya_%d.corr.mat.txt", FilePrefix, TotalSz, 
		      NbrSitesXA, NbrSitesYA);
	    }
 	  TmpEntanglementHamiltonian.WriteAsciiMatrix(TmpMatrixFileName, true);
// 	  ComplexMatrix TmpMatrix1(TmpEntanglementHamiltonian);
// 	  ComplexMatrix TmpMatrix2(ThermalCorrelationMatrix);
// 	  ComplexMatrix TmpMatrix3 = HermitianMultiply(TmpMatrix1, TmpMatrix2);
// 	  ComplexMatrix TmpMatrix4 = HermitianMultiply(TmpMatrix1, TmpMatrix1);
// 	  ComplexMatrix TmpMatrix5 = HermitianMultiply(TmpMatrix2, TmpMatrix2);
// 	  cout << "scalar = " << (TmpMatrix3.ComplexTr() / (sqrt(TmpMatrix4.Tr() * TmpMatrix5.Tr()))) << endl;
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
      
      double* EntanglementEntropies = new double[NbrRenyiEntropies];
      double* NonVacuumEntanglementEntropies = new double[NbrRenyiEntropies];
      for (int i = 0; i < NbrRenyiEntropies; ++i)
	{
	  EntanglementEntropies[i] = 0.0;
	  NonVacuumEntanglementEntropies[i] = 0.0;
	}
      int NbrVacuumOneBodyEntanglementTrimmedEnergies = MaxOneBodyEntanglementEnergyIndex - MinOneBodyEntanglementEnergyIndex + 1;
      int NbrRejectedOneBodyEntropies = VacuumOneBodyEntanglementEnergies.GetNbrRow() - NbrVacuumOneBodyEntanglementTrimmedEnergies;
      double* VacuumOneBodyEntanglementTrimmedEnergies = new double[NbrVacuumOneBodyEntanglementTrimmedEnergies];
      for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies; ++i)
	{
	  VacuumOneBodyEntanglementTrimmedEnergies[i] = VacuumOneBodyEntanglementEnergies[i + MinOneBodyEntanglementEnergyIndex];	
	  OneBodyFile << NbrSitesXA << " " << NbrSitesYA << " " << VacuumOneBodyEntanglementTrimmedEnergies[i] << endl;
	}
      double* SumAlphaFactors = new double[NbrRenyiEntropies];      
      for (int i = 0; i < NbrRenyiEntropies; ++i)
	{
	  SumAlphaFactors[i] = 0.0;
	}
      if (NbrPairs == 0)
	{
	  for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies; ++i)
	    {
	      EntanglementEntropies[0] -= VacuumOneBodyEntanglementTrimmedEnergies[i] * log (VacuumOneBodyEntanglementTrimmedEnergies[i]);
	      EntanglementEntropies[0] -= (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i]) * log (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i]);
	    }
	  NonVacuumEntanglementEntropies[0] = EntanglementEntropies[0];	
	  for (int j = 1; j < NbrRenyiEntropies; ++j)
	    {
	      for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies; ++i)
		{
		  EntanglementEntropies[j] -= 1.0 / ((double) j) * log(powl(VacuumOneBodyEntanglementTrimmedEnergies[i], (double) (j + 1)) + powl(1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i], (double) (j + 1)));
		}
	      NonVacuumEntanglementEntropies[j] = EntanglementEntropies[j];
	    }
	  for (int i = 0; i < NbrRenyiEntropies; ++i)
	    {
	      SumAlphaFactors[i] = 1.0;
	    }
	}
      else
	{
	  int MaxNbrParticlesA = VacuumNbrParticles;
	  if (MaxNbrParticlesA > TotalNbrSitesA)
	    MaxNbrParticlesA = TotalNbrSitesA;
	  LongRational TmpCoefficient;
	  int MaxBinomial = NbrSites;
	  if (Manager.GetBoolean("use-rational") == true)
	    {
	      MaxBinomial = 2;
	    }
	  BinomialCoefficients  Binomial(MaxBinomial);
	  double CurrentEntanglementEntropyContribution = 1.0;


	  if (Manager.GetBoolean("use-approximation") == false)
	    {
	      int OptimalNbrParticlesA = (TotalNbrSitesA * VacuumNbrParticles) / NbrSites;
	      for (int TmpNbrParticlesA = OptimalNbrParticlesA; ((TmpNbrParticlesA <= MaxNbrParticlesA) && 
								 ((EntanglementEntropies[0] + CurrentEntanglementEntropyContribution) != EntanglementEntropies[0])); ++TmpNbrParticlesA)
		{
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalStartingTime), 0);
		    }
		  double AlphaFactor = 0.0;
		  double* TmpEntanglementEntropies = new double[NbrRenyiEntropies];
		  for (int i = 0; i < NbrRenyiEntropies; ++i)
		    TmpEntanglementEntropies[i] = 0.0;
		  GetEntanglementEntropyPerNbrParticlesA(VacuumOneBodyEntanglementTrimmedEnergies, NbrVacuumOneBodyEntanglementTrimmedEnergies, TmpNbrParticlesA, 
							 0, 1.0, TmpEntanglementEntropies, NbrRenyiEntropies, AlphaFactor);
		  int MaxSumIndex = TotalNbrSitesA - TmpNbrParticlesA;
		  if (MaxSumIndex > NbrPairs)
		    MaxSumIndex = NbrPairs;
		  double* Tmp = EvaluateEtaPairingContribution(NbrRenyiEntropies, NbrSites, NbrPairs, VacuumNbrParticles, TotalNbrSitesA, TmpNbrParticlesA,
							       Manager.GetBoolean("use-rational"), TmpCoefficient, Binomial);
		  CurrentEntanglementEntropyContribution = (AlphaFactor * Tmp[0]) + TmpEntanglementEntropies[0];
		  EntanglementEntropies[0] += CurrentEntanglementEntropyContribution;
		  NonVacuumEntanglementEntropies[0] += TmpEntanglementEntropies[0];
		  for (int j = 1; j < NbrRenyiEntropies; ++j)
		    {
		      double TmpCurrentEntanglementEntropyContribution = TmpEntanglementEntropies[j] * Tmp[j];
		      EntanglementEntropies[j] += TmpCurrentEntanglementEntropyContribution;
		      NonVacuumEntanglementEntropies[j] += TmpEntanglementEntropies[j];
		    }
		  delete[] Tmp;
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalEndingTime), 0);
		      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		      cout << TmpNbrParticlesA << " : " << EntanglementEntropies[0] << " " << CurrentEntanglementEntropyContribution << " " << AlphaFactor << " " << "(" << Dt << "s)" << endl;
		    }
		  else
		    {
		      cout << TmpNbrParticlesA << " : " << EntanglementEntropies[0] << " " << CurrentEntanglementEntropyContribution << " " << AlphaFactor << endl;
		    }	  
		  if (NbrRenyiEntropies > 1)
		    {
		      cout << "Renyi alphas : ";
		      for (int j = 1; j < NbrRenyiEntropies; ++j)
			{
			  cout << " " << TmpEntanglementEntropies[j];
			}
		      cout << endl;	      
		    }
		  delete[] TmpEntanglementEntropies;
		}
	      CurrentEntanglementEntropyContribution = EntanglementEntropies[0];
	      for (int TmpNbrParticlesA = OptimalNbrParticlesA - 1; ((TmpNbrParticlesA >= 0) && 
								     ((EntanglementEntropies[0] + CurrentEntanglementEntropyContribution) != EntanglementEntropies[0])); --TmpNbrParticlesA)
		{
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalStartingTime), 0);
		    }
		  double AlphaFactor = 0.0;
		  double* TmpEntanglementEntropies = new double[NbrRenyiEntropies];
		  for (int i = 0; i < NbrRenyiEntropies; ++i)
		    {
		      TmpEntanglementEntropies[i] = 0.0;
		    }
		  GetEntanglementEntropyPerNbrParticlesA(VacuumOneBodyEntanglementTrimmedEnergies, NbrVacuumOneBodyEntanglementTrimmedEnergies, TmpNbrParticlesA, 
							 0, 1.0, TmpEntanglementEntropies, NbrRenyiEntropies, AlphaFactor);
		  int MaxSumIndex = TotalNbrSitesA - TmpNbrParticlesA;
		  if (MaxSumIndex > NbrPairs)
		    MaxSumIndex = NbrPairs;
		  double* Tmp = EvaluateEtaPairingContribution(NbrRenyiEntropies, NbrSites, NbrPairs, VacuumNbrParticles, TotalNbrSitesA, TmpNbrParticlesA,
							       Manager.GetBoolean("use-rational"), TmpCoefficient, Binomial);
		  CurrentEntanglementEntropyContribution = (AlphaFactor * Tmp[0]) + TmpEntanglementEntropies[0];
		  EntanglementEntropies[0] += CurrentEntanglementEntropyContribution;
		  NonVacuumEntanglementEntropies[0] += TmpEntanglementEntropies[0];
		  for (int j = 1; j < NbrRenyiEntropies; ++j)
		    {
		      double TmpCurrentEntanglementEntropyContribution = TmpEntanglementEntropies[j] * Tmp[j];
		      EntanglementEntropies[j] += TmpCurrentEntanglementEntropyContribution;
		      NonVacuumEntanglementEntropies[j] += TmpEntanglementEntropies[j];
		    }
		  delete[] Tmp;
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalEndingTime), 0);
		      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		      cout << TmpNbrParticlesA << " : " << EntanglementEntropies[0] << " " << CurrentEntanglementEntropyContribution << " " << AlphaFactor << " " << "(" << Dt << "s)" << endl;
		    }
		  else
		    {
		      cout << TmpNbrParticlesA << " : " << EntanglementEntropies[0] << " " << CurrentEntanglementEntropyContribution << " " << AlphaFactor << endl;
		    }	  
		  if (NbrRenyiEntropies > 1)
		    {
		      cout << "Renyi alphas : ";
		      for (int j = 1; j < NbrRenyiEntropies; ++j)
			{
			  cout << " " << TmpEntanglementEntropies[j];
			}
		      cout << endl;	      
		    }
		  delete[] TmpEntanglementEntropies;
		}
	    }
	  else
	    {	  

	      // use the approximated formula

	      int* OptimalNbrParticlesAPerRenyiEntropy = new int [NbrRenyiEntropies];
	      OptimalNbrParticlesAPerRenyiEntropy[0] = (TotalNbrSitesA * VacuumNbrParticles) / NbrSites;
	      cout << "optimal number of particles for the Von Neumann entropy  = " << OptimalNbrParticlesAPerRenyiEntropy[0] << endl;
	      for (int j = 1; j < NbrRenyiEntropies; ++j)
		{
		  double Tmp = 0.0;
		  double Sigma = 0.0;
		  double TmpRenyiIndex = (double) (j + 1);
		  double TmpEntanglementEntropy = 0.0;
		  for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies; ++i)
		    {
		      double Tmp2 = (powl(VacuumOneBodyEntanglementTrimmedEnergies[i], TmpRenyiIndex) / 
				     (powl(VacuumOneBodyEntanglementTrimmedEnergies[i], TmpRenyiIndex) + 
				      powl(1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i], TmpRenyiIndex)));
		      Tmp += Tmp2;
		      Sigma += Tmp2 * Tmp2 * powl((1.0 / VacuumOneBodyEntanglementTrimmedEnergies[i]) - 1.0, TmpRenyiIndex);
		      TmpEntanglementEntropy -= (log(powl(VacuumOneBodyEntanglementTrimmedEnergies[i], TmpRenyiIndex) + 
						     powl(1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i], TmpRenyiIndex)));
		    }
		  TmpEntanglementEntropy /= (double) j;
		  OptimalNbrParticlesAPerRenyiEntropy[j] = (int) Tmp;
		  cout << "optimal number of particles for the Renyi entropy " << (j + 1) << "  = " << OptimalNbrParticlesAPerRenyiEntropy[j] << " (" << Tmp << ")" << ", sigma=" << Sigma << " " << TmpEntanglementEntropy << endl;		  
		}

	      // Von Neumann entropy

	      double AlphaFactor = 1.0;
	      int OptimalNbrParticlesA = (TotalNbrSitesA * VacuumNbrParticles) / NbrSites;
	      for (int TmpNbrParticlesA = OptimalNbrParticlesA; ((TmpNbrParticlesA <= MaxNbrParticlesA) && 
								 ((SumAlphaFactors[0] + AlphaFactor) != SumAlphaFactors[0])); ++TmpNbrParticlesA)
		{
		  double TmpZ0 = GetZ0Value(VacuumOneBodyEntanglementTrimmedEnergies, NbrVacuumOneBodyEntanglementTrimmedEnergies, TmpNbrParticlesA);
		  double LogAlphaFactor = (-((double) (TmpNbrParticlesA + 1)) * log(TmpZ0)) - (0.5 * log (2.0 * M_PI));
		  double TmpSum = ((double) (TmpNbrParticlesA + 1)) / (TmpZ0 * TmpZ0);
		  for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies; ++i)
		    {
		      LogAlphaFactor += log((1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i] + 
					     (TmpZ0 * VacuumOneBodyEntanglementTrimmedEnergies[i])));
		      TmpSum -= ((VacuumOneBodyEntanglementTrimmedEnergies[i] * VacuumOneBodyEntanglementTrimmedEnergies[i]) 
				 / ((1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i] + (TmpZ0 * VacuumOneBodyEntanglementTrimmedEnergies[i])) 
				    * (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i] + (TmpZ0 * VacuumOneBodyEntanglementTrimmedEnergies[i]))));
		    }
		  AlphaFactor = exp(LogAlphaFactor);
		  AlphaFactor /= sqrt(TmpSum);
		  SumAlphaFactors[0] += AlphaFactor;
		}
	      AlphaFactor = SumAlphaFactors[0];
	      for (int TmpNbrParticlesA = OptimalNbrParticlesA - 1; ((TmpNbrParticlesA >= 0) && 
								     ((SumAlphaFactors[0] + AlphaFactor) != SumAlphaFactors[0])); --TmpNbrParticlesA)
		{
		  double TmpZ0 = GetZ0Value(VacuumOneBodyEntanglementTrimmedEnergies, NbrVacuumOneBodyEntanglementTrimmedEnergies, TmpNbrParticlesA);
		  double TmpSum = ((double) (TmpNbrParticlesA + 1)) / (TmpZ0 * TmpZ0);
		  double LogAlphaFactor = (-((double) (TmpNbrParticlesA + 1)) * log(TmpZ0)) - (0.5 * log (2.0 * M_PI));
		  for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies; ++i)
		    {
		      LogAlphaFactor += log(1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i] + 
					    (TmpZ0 * VacuumOneBodyEntanglementTrimmedEnergies[i]));
		      TmpSum -= ((VacuumOneBodyEntanglementTrimmedEnergies[i] * VacuumOneBodyEntanglementTrimmedEnergies[i]) 
				 / ((1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i] + (TmpZ0 * VacuumOneBodyEntanglementTrimmedEnergies[i])) 
				    * (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i] + (TmpZ0 * VacuumOneBodyEntanglementTrimmedEnergies[i]))));
		    }
		  AlphaFactor = exp(LogAlphaFactor);
		  AlphaFactor /= sqrt(TmpSum);
		  SumAlphaFactors[0] += AlphaFactor;
		}
	      for (int TmpNbrParticlesA = OptimalNbrParticlesA; ((TmpNbrParticlesA <= MaxNbrParticlesA) && 
								 ((EntanglementEntropies[0] + CurrentEntanglementEntropyContribution) != EntanglementEntropies[0])); ++TmpNbrParticlesA)
		{
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalStartingTime), 0);
		    }
		  double AlphaFactor = 0.0;
		  double TmpEntanglementEntropy = 0.0;
		  double TmpZ0 = GetZ0Value(VacuumOneBodyEntanglementTrimmedEnergies, NbrVacuumOneBodyEntanglementTrimmedEnergies, TmpNbrParticlesA);
		  AlphaFactor = pow (TmpZ0, -((double) (TmpNbrParticlesA + 1))) / sqrt (2.0 * M_PI);
		  double TmpSum = ((double) (TmpNbrParticlesA + 1)) / (TmpZ0 * TmpZ0);
		  double LogAlphaFactor = (-((double) (TmpNbrParticlesA + 1)) * log(TmpZ0)) - (0.5 * log (2.0 * M_PI));
		  for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies; ++i)
		    {
		      LogAlphaFactor += log((1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i] + 
					     (TmpZ0 * VacuumOneBodyEntanglementTrimmedEnergies[i])));
		      TmpSum -= ((VacuumOneBodyEntanglementTrimmedEnergies[i] * VacuumOneBodyEntanglementTrimmedEnergies[i]) 
				 / ((1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i] + (TmpZ0 * VacuumOneBodyEntanglementTrimmedEnergies[i])) 
				    * (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i] + (TmpZ0 * VacuumOneBodyEntanglementTrimmedEnergies[i]))));
		    }
		  AlphaFactor = exp(LogAlphaFactor);
		  AlphaFactor /= sqrt(TmpSum) * SumAlphaFactors[0];
		  int MaxSumIndex = TotalNbrSitesA - TmpNbrParticlesA;
		  if (MaxSumIndex > NbrPairs)
		    MaxSumIndex = NbrPairs;
		  double* Tmp = EvaluateEtaPairingContribution(1, NbrSites, NbrPairs, VacuumNbrParticles, TotalNbrSitesA, TmpNbrParticlesA,
							       Manager.GetBoolean("use-rational"), TmpCoefficient, Binomial);
		  CurrentEntanglementEntropyContribution = (AlphaFactor * Tmp[0]) + TmpEntanglementEntropy;
		  EntanglementEntropies[0] += CurrentEntanglementEntropyContribution;
		  NonVacuumEntanglementEntropies[0] += TmpEntanglementEntropy;
		  delete[] Tmp;
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalEndingTime), 0);
		      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		      cout << TmpNbrParticlesA << " : " << EntanglementEntropies[0] << " " << CurrentEntanglementEntropyContribution << " " << AlphaFactor << " " << "(" << Dt << "s)" << endl;
		    }
		  else
		    {
		      cout << TmpNbrParticlesA << " : " << EntanglementEntropies[0] << " " << CurrentEntanglementEntropyContribution << " " << AlphaFactor << endl;
		    }	  
		}
	      CurrentEntanglementEntropyContribution = EntanglementEntropies[0];
	      for (int TmpNbrParticlesA = OptimalNbrParticlesA - 1; ((TmpNbrParticlesA >= 0) && 
								     ((EntanglementEntropies[0] + CurrentEntanglementEntropyContribution) != EntanglementEntropies[0])); --TmpNbrParticlesA)
		{
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalStartingTime), 0);
		    }
		  double AlphaFactor = 0.0;
		  double TmpEntanglementEntropy = 0.0;
		  double TmpZ0 = GetZ0Value(VacuumOneBodyEntanglementTrimmedEnergies, NbrVacuumOneBodyEntanglementTrimmedEnergies, TmpNbrParticlesA);
		  AlphaFactor = pow (TmpZ0, -((double) (TmpNbrParticlesA + 1))) / sqrt (2.0 * M_PI);
		  double TmpSum = ((double) (TmpNbrParticlesA + 1)) / (TmpZ0 * TmpZ0);
		  double LogAlphaFactor = (-((double) (TmpNbrParticlesA + 1)) * log(TmpZ0)) - (0.5 * log (2.0 * M_PI));
		  for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies; ++i)
		    {
		      LogAlphaFactor += log((1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i] + 
					     (TmpZ0 * VacuumOneBodyEntanglementTrimmedEnergies[i])));
		      TmpSum -= ((VacuumOneBodyEntanglementTrimmedEnergies[i] * VacuumOneBodyEntanglementTrimmedEnergies[i]) 
				 / ((1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i] + (TmpZ0 * VacuumOneBodyEntanglementTrimmedEnergies[i])) 
				    * (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i] + (TmpZ0 * VacuumOneBodyEntanglementTrimmedEnergies[i]))));
		    }
		  AlphaFactor = exp(LogAlphaFactor);
		  AlphaFactor /= sqrt(TmpSum) * SumAlphaFactors[0];
		
		  int MaxSumIndex = TotalNbrSitesA - TmpNbrParticlesA;
		  if (MaxSumIndex > NbrPairs)
		    MaxSumIndex = NbrPairs;
		  double* Tmp = EvaluateEtaPairingContribution(1, NbrSites, NbrPairs, VacuumNbrParticles, TotalNbrSitesA, TmpNbrParticlesA,
							       Manager.GetBoolean("use-rational"), TmpCoefficient, Binomial);
		  CurrentEntanglementEntropyContribution = (AlphaFactor * Tmp[0]) + TmpEntanglementEntropy;
		  EntanglementEntropies[0] += CurrentEntanglementEntropyContribution;
		  NonVacuumEntanglementEntropies[0] += TmpEntanglementEntropy;
		  delete[] Tmp;
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalEndingTime), 0);
		      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		      cout << TmpNbrParticlesA << " : " << EntanglementEntropies[0] << " " << CurrentEntanglementEntropyContribution << " " << AlphaFactor << " " << "(" << Dt << "s)" << endl;
		    }
		  else
		    {
		      cout << TmpNbrParticlesA << " : " << EntanglementEntropies[0] << " " << CurrentEntanglementEntropyContribution << " " << AlphaFactor << endl;
		    }	  
		}


	      // Renyi entropies
	      for (int RenyiIndex = 2; RenyiIndex <= NbrRenyiEntropies; ++RenyiIndex)
		{
		  cout << "Warning, --use-approximation is not working for the Renyi entropies" << endl;
		}
	    }
	}
      NonVacuumEntanglementEntropies[0] = 0.0;
      for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies; ++i)
	{
	  NonVacuumEntanglementEntropies[0] -= VacuumOneBodyEntanglementTrimmedEnergies[i] * log (VacuumOneBodyEntanglementTrimmedEnergies[i]);
	  NonVacuumEntanglementEntropies[0] -= (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i]) * log (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i]);
	}
      EntanglementEntropies[0] += NonVacuumEntanglementEntropies[0];
      if (NbrPairs != 0)
	{
	  for (int j = 1; j < NbrRenyiEntropies; ++j)
	    {
	      EntanglementEntropies[j] = -1.0 / ((double) j) * log(EntanglementEntropies[j]);
	      NonVacuumEntanglementEntropies[j] = -1.0 / ((double) j) * log(NonVacuumEntanglementEntropies[j]);
	    }
	}
      cout << "Normalization = " << SumAlphaFactors[0] << endl;
      cout << "Entanglement entropy = " << EntanglementEntropies[0] << endl;
      for (int i = 1; i < NbrRenyiEntropies; ++i)
	cout << "Renyi entropy " << (i + 1) << " = " << EntanglementEntropies[i] << " " << NonVacuumEntanglementEntropies[i] << endl;
      cout << "Nbr Rejected one-body entanglement energies = " << NbrRejectedOneBodyEntropies << " / " << TotalNbrSitesA << endl;
      File << NbrSitesXA << " " << NbrSitesYA << " " << EntanglementEntropies[0] << " " << NonVacuumEntanglementEntropies[0];
      if (NbrRenyiEntropies > 1)
	{
	  for (int i = 1; i < NbrRenyiEntropies; ++i)
	    File << " " << EntanglementEntropies[i] << " " << NonVacuumEntanglementEntropies[i];
	}
      File << endl;
      delete[] VacuumOneBodyEntanglementTrimmedEnergies;
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

HermitianMatrix EtaPairaingEntanglementEntropyExtractCorrelationMatrix(HermitianMatrix& correlationMatrix, int sourceNbrSitesX, int sourceNbrSitesY, 
								       int targetNbrSitesX, int targetNbrSitesY)
{
  int TargetTotalNbrSites = targetNbrSitesX * targetNbrSitesY;
  HermitianMatrix TmpMatrix (TargetTotalNbrSites, true);
  Complex Tmp;
  for (int i = 0; i < TargetTotalNbrSites; ++i)
    {
      int TmpTargetX1 = i / targetNbrSitesY;
      int TmpTargetY1 = i % targetNbrSitesY;
      int TmpSourceIndex1 = (TmpTargetX1 * sourceNbrSitesY) + TmpTargetY1;
      for (int j = i; j < TargetTotalNbrSites; ++j)
	{
	  int TmpTargetX2 = j / targetNbrSitesY;
	  int TmpTargetY2 = j % targetNbrSitesY;
	  int TmpSourceIndex2 = (TmpTargetX2 * sourceNbrSitesY) + TmpTargetY2;
	  correlationMatrix.GetMatrixElement(TmpSourceIndex1, TmpSourceIndex2, Tmp);
	  TmpMatrix.SetMatrixElement(i, j, Tmp);
	}
    }
  return TmpMatrix;
}


// compute the contribution to the Renyi entropy for a given number of particles in the subregion A, using the exact evaluation of the partitions
//
// oneBodyEntanglementTrimmedEnergies= array containing the one-body entanglement energies
// nbrOneBodyEntanglementTrimmedEnergies = number of one-body entanglement energies
// nbrParticlesA = number of particles in the subsystem A
// currentOrbitalIndex = current orbital that is considered
// currentFactor = current factor for a single partition
// entropies = array that contains the Renyi entropies
// maxRenyiEntropy = maximum Renyi entropy that has to be evaluated
// alpha = reference total weight of the reduced density matrix for the sector with nbrParticlesA particles

void GetEntanglementEntropyPerNbrParticlesA(double* oneBodyEntanglementTrimmedEnergies, int nbrOneBodyEntanglementTrimmedEnergies, 
					    int nbrParticlesA, int currentOrbitalIndex, double currentFactor, double* entropies, int maxRenyiEntropy, double& alpha)
{
  if (currentOrbitalIndex > nbrOneBodyEntanglementTrimmedEnergies)
    return;
  if (nbrParticlesA == 0)
    {      
      for (; currentOrbitalIndex < nbrOneBodyEntanglementTrimmedEnergies; ++currentOrbitalIndex)
	{
	  currentFactor *= (1.0 - oneBodyEntanglementTrimmedEnergies[currentOrbitalIndex]);
	}
      entropies[0] -= currentFactor * log(currentFactor);
      for (int i = 1; i < maxRenyiEntropy; ++i)
	{
	  entropies[i] += pow(currentFactor, (double) (i + 1));
	}
      alpha += currentFactor;
      return;
    }
  if (currentOrbitalIndex == nbrOneBodyEntanglementTrimmedEnergies)
    return;
  GetEntanglementEntropyPerNbrParticlesA(oneBodyEntanglementTrimmedEnergies, nbrOneBodyEntanglementTrimmedEnergies, nbrParticlesA, currentOrbitalIndex + 1, currentFactor * (1.0 - oneBodyEntanglementTrimmedEnergies[currentOrbitalIndex]), entropies, maxRenyiEntropy, alpha);
  GetEntanglementEntropyPerNbrParticlesA(oneBodyEntanglementTrimmedEnergies, nbrOneBodyEntanglementTrimmedEnergies, nbrParticlesA - 1, currentOrbitalIndex + 1, currentFactor * oneBodyEntanglementTrimmedEnergies[currentOrbitalIndex], entropies, maxRenyiEntropy, alpha);
  return;
}

// compute the z_0 value obtained from the saddle point approximation
//
// oneBodyEntanglementTrimmedEnergies= array containing the one-body entanglement energies
// nbrOneBodyEntanglementTrimmedEnergies = number of one-body entanglement energies
// nbrParticlesA = number of particles in the subsystem A

double GetZ0Value (double* oneBodyEntanglementTrimmedEnergies, int nbrOneBodyEntanglementTrimmedEnergies, int nbrParticlesA)
{
  double MinZ0 = 0.0;
  double MaxZ0 = 10000.0;
  while (GetZ0DefintionSum(oneBodyEntanglementTrimmedEnergies, nbrOneBodyEntanglementTrimmedEnergies, nbrParticlesA, MaxZ0) < 0.0)
    {
      MaxZ0 *= 10.0;
    }
  while (fabs(MinZ0 - MaxZ0) > (MACHINE_PRECISION * fabs(MaxZ0)))
    {
      double TmpZ0 = 0.5 * (MaxZ0 + MinZ0);
      double TmpSum =  GetZ0DefintionSum(oneBodyEntanglementTrimmedEnergies, nbrOneBodyEntanglementTrimmedEnergies, nbrParticlesA, TmpZ0);
      if (TmpSum > 0.0)
	MaxZ0 = TmpZ0;
      else
	MinZ0 = TmpZ0;
    }
  return MaxZ0;
}

// compute the z_0 value obtained from the saddle point approximation for a generic Renyi entropy
//
// oneBodyEntanglementTrimmedEnergies= array containing the one-body entanglement energies
// nbrOneBodyEntanglementTrimmedEnergies = number of one-body entanglement energies
// nbrParticlesA = number of particles in the subsystem A
// renyiIndex = Renyi index

double GetZ0Value (double* oneBodyEntanglementTrimmedEnergies, int nbrOneBodyEntanglementTrimmedEnergies, int nbrParticlesA, double renyiIndex)
{
  double MinZ0 = 0.0;
  double MaxZ0 = 10000.0;
  while (GetZ0DefintionSum(oneBodyEntanglementTrimmedEnergies, nbrOneBodyEntanglementTrimmedEnergies, nbrParticlesA, MaxZ0, renyiIndex) < 0.0)
    {
      MaxZ0 *= 10.0;
    }
  while (fabs(MinZ0 - MaxZ0) > (MACHINE_PRECISION * fabs(MaxZ0)))
    {
      double TmpZ0 = 0.5 * (MaxZ0 + MinZ0);
      double TmpSum =  GetZ0DefintionSum(oneBodyEntanglementTrimmedEnergies, nbrOneBodyEntanglementTrimmedEnergies, nbrParticlesA, TmpZ0, renyiIndex);
      if (TmpSum > 0.0)
	MaxZ0 = TmpZ0;
      else
	MinZ0 = TmpZ0;
    }
  return MaxZ0;
}

double GetZ0DefintionSum (double* oneBodyEntanglementTrimmedEnergies, int nbrOneBodyEntanglementTrimmedEnergies, 
			  int nbrParticlesA, double z0)
{
  double TmpSum = 0.0;
  for (int i = 0; i < nbrOneBodyEntanglementTrimmedEnergies; ++i)
    {
      TmpSum += oneBodyEntanglementTrimmedEnergies[i]  / ((1.0 - oneBodyEntanglementTrimmedEnergies[i]) + (z0 * oneBodyEntanglementTrimmedEnergies[i]));
    }
  TmpSum *= z0;
  return (TmpSum - ((double) (nbrParticlesA + 1)));
}

double GetZ0DefintionSum (double* oneBodyEntanglementTrimmedEnergies, int nbrOneBodyEntanglementTrimmedEnergies, 
			  int nbrParticlesA, double z0, double renyiIndex)
{
  double TmpSum = 0.0;
  for (int i = 0; i < nbrOneBodyEntanglementTrimmedEnergies; ++i)
    {
      TmpSum += pow(oneBodyEntanglementTrimmedEnergies[i], renyiIndex)  / (pow(1.0 - oneBodyEntanglementTrimmedEnergies[i], renyiIndex) 
									   + (z0 * pow(oneBodyEntanglementTrimmedEnergies[i], renyiIndex)));
    }
  TmpSum *= z0;
  return (TmpSum - ((double) (nbrParticlesA + 1)));
}

void ComputeThermalQuantities (double beta, double mu, int nbrStates, double* stateEnergies, double& thermalNbrParticles, double& thermalEnergy, double* thermalEntropy, int maxRenyiEntropy, 
			       double& thermalNbrParticlesMuDerivative, double& thermalNbrParticlesBetaDerivative, 
			       double& thermalEnergyMuDerivative, double& thermalEnergyBetaDerivative)
{
  thermalEnergy = 0.0;
  thermalNbrParticles = 0.0;
  thermalEnergyMuDerivative = 0.0;
  thermalEnergyBetaDerivative = 0.0;
  thermalNbrParticlesMuDerivative = 0.0;
  thermalNbrParticlesBetaDerivative = 0.0;
  for (int i = 0; i < maxRenyiEntropy; ++i)    
    thermalEntropy[i] = 0.0;
  for (int k = 0 ; k < nbrStates; ++k)
    {
      double TmpExp = exp (beta * (stateEnergies[k] - mu));
      double Tmp = 1.0 / (1.0 + TmpExp);
      thermalNbrParticles += Tmp;
      thermalEnergy += Tmp * stateEnergies[k];
      thermalEntropy[0] -= (Tmp * log (Tmp)) + ((1.0 - Tmp) * log (1.0 - Tmp));
      for (int i = 1; i < maxRenyiEntropy; ++i)    
	thermalEntropy[i] += -1.0 / ((double) i) * log(powl(Tmp, (double) (i + 1)) + powl(1.0 - Tmp, (double) (i + 1)));
      Tmp *= Tmp;
      Tmp *= TmpExp;
      thermalEnergyMuDerivative += Tmp * beta * stateEnergies[k];
      thermalEnergyBetaDerivative -= Tmp * stateEnergies[k] * (stateEnergies[k] - mu);
      thermalNbrParticlesMuDerivative += Tmp * beta;
      thermalNbrParticlesBetaDerivative -= Tmp * (stateEnergies[k] - mu);
    }
}

// evaluate the contribution of the eta pairing to the Renyi entropies
// 
// nbrRenyiEntropies = number of Renyi Entropies to evaluate
// nbrSites = total number of sites 
// nbrPairs = number of eta pairing pairs
// vacuumNbrParticles = number of particles for the vacuum states 
// totalNbrSitesA = otal number of sites in the region A
// nbrParticlesA = number of particles in the region A
// useRational = true if rational numbers have to be used for intermediate calculations
// rationalCoefficient = reference on the temporary rational coefficient
// binomial = reference on the binomial coefficients
// return value = array containing the eta pairing contribution for each Renyi entropy

double* EvaluateEtaPairingContribution(int nbrRenyiEntropies, int nbrSites, int nbrPairs, int vacuumNbrParticles, int totalNbrSitesA, int nbrParticlesA,
				       bool useRational, LongRational& rationalCoefficient, BinomialCoefficients& binomial)
{
  double* Tmp = new double[nbrRenyiEntropies];
  for (int j = 0; j < nbrRenyiEntropies; ++j)
    {
      Tmp[j] = 0.0;
    }
  int MaxSumIndex = totalNbrSitesA;
  if (MaxSumIndex > nbrPairs)
    MaxSumIndex = nbrPairs;
  for (int j = 0; j <= MaxSumIndex; ++j)
    {
      double Tmp2;
      if (useRational == true)
	{
	  rationalCoefficient.SetToOne();
	  rationalCoefficient.BinomialMultiply(totalNbrSitesA - nbrParticlesA, j);
	  rationalCoefficient.BinomialMultiply(nbrSites - totalNbrSitesA - vacuumNbrParticles + nbrParticlesA, nbrPairs - j);
	  rationalCoefficient.BinomialDivide(nbrSites - vacuumNbrParticles, nbrPairs);
	  Tmp2 = rationalCoefficient.GetNumericalValue();
//	  cout << Tmp2 << endl;
	}
      else
	{
	  Tmp2 = ((binomial.GetNumericalCoefficient(totalNbrSitesA - nbrParticlesA, j) / binomial.GetNumericalCoefficient(nbrSites - vacuumNbrParticles, nbrPairs))
		  * binomial.GetNumericalCoefficient(nbrSites - totalNbrSitesA - vacuumNbrParticles + nbrParticlesA, nbrPairs - j));
	}
      if (Tmp2 > 0.0)
	{
	  Tmp[0] -= Tmp2 * log(Tmp2);
	  for (int k = 1; k < nbrRenyiEntropies; ++k)
	    {
	      Tmp[k] += pow(Tmp2, (double) (k + 1));
	    }
	}
    }
  return Tmp;
}
