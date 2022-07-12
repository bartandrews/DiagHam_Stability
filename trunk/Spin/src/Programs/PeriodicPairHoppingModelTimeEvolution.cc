#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslations.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainWithTranslations.h"
#include "HilbertSpace/PairHoppingGenericPAsSpin1ChainWithTranslations.h"
#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetry.h"
#include "HilbertSpace/PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry.h"

#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslationsLong.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainWithTranslationsLong.h"
#include "HilbertSpace/PairHoppingGenericPAsSpin1ChainWithTranslationsLong.h"
#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong.h"
#include "HilbertSpace/PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetryLong.h"

#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Hamiltonian/PairHoppingHamiltonianWithTranslations.h"
#include "Hamiltonian/PairHoppingRealHamiltonianWithTranslations.h"

#include "Matrix/ComplexMatrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"
#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"

#include "config.h"


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("PeriodicPairHoppingModelTimeEvolution" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new  SingleIntegerOption ('p', "p-value", "value that defines the filling factor p/(2p+1)", 1);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-spin", "number of spins", 10);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "momentum", "momentum sector", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "inversion-symmetry", "set the inversion symmetry sector (+1 or -1)", 1);
  (*SystemGroup) += new  SingleDoubleOption ('t', "time-step", "time step", 0.1);
  (*SystemGroup) += new  SingleIntegerOption ('T', "nbr-timesteps", "number of time steps", 100);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "max-nbriterations", "maximum number of iterations when performing the time evolution", 100);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "norm-accuracy", "accuracy required when computing the state norm at a given time t", 1e-13); 
  (*SystemGroup) += new  SingleDoubleOption ('\n', "v1-0", "strength of the first electrostatic coupling", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "v2-0", "strength of the second electrostatic coupling", 0.0);
  (*OutputGroup) += new BooleanOption  ('\n', "use-root", "us the root state instead of the z2 state as the initial state");
  (*OutputGroup) += new BooleanOption  ('\n', "store-states", "store all the states generated at each time step in a single binary matrix file");
  (*OutputGroup) += new BooleanOption  ('\n', "store-initialstate", "store all the initial state");
  (*OutputGroup) += new BooleanOption  ('\n', "disable-entropy", "do not evaluate the entanglement entropy  at each time step");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PeriodicPairHoppingModelTimeEvolution -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int PValue = Manager.GetInteger("p-value");
  int NbrSpins = Manager.GetInteger("nbr-spin");
  if ((NbrSpins % PValue) != 0)
    {
      cout << "--nbr-spin (here " << NbrSpins << ") should be a multiple of --p-value (here " << PValue << ")" << endl;
      return -1;
    }
  int SubsystemSize = (((NbrSpins / PValue) / 2)) * PValue;

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  int NbrTimeSteps =  Manager.GetInteger("nbr-timesteps");
  double TimeStep = Manager.GetDouble("time-step");
  int MaxNbrIterations = Manager.GetInteger("max-nbriterations");
  double NormConvergence = Manager.GetDouble("norm-accuracy");
  double V10 = Manager.GetDouble("v1-0");
  double V20 = Manager.GetDouble("v2-0");
  bool FirstRun = true;
  int Momentum = Manager.GetInteger("momentum");
  if  (!((Momentum == 0) || ((((NbrSpins / PValue) & 1) == 0) && (Momentum == ((NbrSpins / PValue) >> 1)))))
    {
      cout << "only k=0 and pi momenta are available" << endl;
      return -1;
    }
  int InversionSymmetrySector = Manager.GetInteger("inversion-symmetry");

  
  char* OutputFileName = new char [512];
  if ((V10 != 0.0) || (V20 != 0.0))
    {
      if (Manager.GetBoolean("use-root") == true)
	{
	  sprintf (OutputFileName, "spin_1_periodicpairhopping_p_%d_timeevol_root_v10_%.6f_v20_%.6f_n_%d_invsym_%d_k_%d", PValue, V10, V20, NbrSpins, InversionSymmetrySector, Momentum);
	}
      else
	{
	  sprintf (OutputFileName, "spin_1_periodicpairhopping_p_%d_timeevol_z2_v10_%.6f_v20_%.6f_n_%d_invsym_%d_k_%d", PValue, V10, V20, NbrSpins, InversionSymmetrySector, Momentum);
	}
    }
  else
    {
      if (Manager.GetBoolean("use-root") == true)
	{
	  sprintf (OutputFileName, "spin_1_periodicpairhopping_p_%d_timeevol_root_n_%d_invsym_%d_k_%d", PValue, NbrSpins, InversionSymmetrySector, Momentum);
	}
      else
	{
	  sprintf (OutputFileName, "spin_1_periodicpairhopping_p_%d_timeevol_z2_n_%d_invsym_%d_k_%d", PValue, NbrSpins, InversionSymmetrySector, Momentum);
	}
    }
  
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 64];
  if (Manager.GetBoolean("disable-entropy") == false)
    {
      sprintf (FullOutputFileName, "%s.ent.fidelity", OutputFileName);
    }
  else
    {
      sprintf (FullOutputFileName, "%s.fidelity", OutputFileName);
    }
  char* StateOutputFileName = new char [strlen(OutputFileName)+ 64];
  sprintf (StateOutputFileName, "%s.eigenvec.mat", OutputFileName);  
  ofstream File;
  File.open(FullOutputFileName, ios::binary | ios::out);
  File.precision(14);
  File << "# periodic pair hopping model with p=" << PValue << " in spin 1 language with " << NbrSpins << " sites K=" << Momentum << " InvSym=" << InversionSymmetrySector << endl;
  File << "# time fidelity S_A Tr(rho_A) nbr_iter state_norm computation_time" << endl;
  AbstractSpinChainWithTranslations* Chain = 0;
  if (NbrSpins <= 32)
    {
      switch (PValue)
	{
	case 1 :
	  Chain = new PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry (NbrSpins, Momentum, InversionSymmetrySector);
	  break;
	case 2 :
	  Chain = new PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetry (NbrSpins, Momentum, InversionSymmetrySector);
	  break;
	default :
	  {
	    Chain = new PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry (NbrSpins, Momentum, InversionSymmetrySector, PValue);
	  }
	}
    }
  else
    {
      switch (PValue)
	{
	case 1 :
	  Chain = new PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong (NbrSpins, Momentum, InversionSymmetrySector);
	  break;
	case 2 :
	  Chain = new PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong (NbrSpins, Momentum, InversionSymmetrySector);
	  break;
	default :
	  {
	    Chain = new PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetryLong (NbrSpins, Momentum, InversionSymmetrySector, PValue);
	    break;
	  }
	}
    }
  cout << "Hilbert space dimension = " <<  Chain->GetHilbertSpaceDimension() << endl;
  
  if (Chain->GetHilbertSpaceDimension() > 0)
    {
      Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());
      PairHoppingHamiltonianWithTranslations* Hamiltonian = 0;
      Hamiltonian = new PairHoppingHamiltonianWithTranslations(Chain, NbrSpins, PValue, V10, V20, Architecture.GetArchitecture(), Memory);
      ComplexVector InitialState (Chain->GetLargeHilbertSpaceDimension(), true);
      timeval StartingTime;
      timeval EndingTime;  
      int InitialStateIndex = 0;
      if (NbrSpins <= 32)
	{
	  unsigned long InitialConfiguration = 0x0ul;
	  if (Manager.GetBoolean("use-root") == true)
	    {
	      for (int i = 0; i < NbrSpins; ++i)
		{
		  InitialConfiguration |= (0x2ul) << (i << 1);
		}
	    }
	  else
	    {
	      for (int i = 0; i < (NbrSpins / PValue); i += 2)
		{
		  for (int j = 0; j < PValue; ++j)
		    {
		      InitialConfiguration |= (0x3ul) << (((i * PValue) + j) << 1);
		    }
		}
	    }
	  InitialStateIndex = ((PairHoppingP1AsSpin1ChainWithTranslations*) Chain)->FindStateIndex(InitialConfiguration);
	}
      else
 	{
	  ULONGLONG InitialConfiguration = ((ULONGLONG) 0x0ul);
	  if (Manager.GetBoolean("use-root") == true)
	    {
	      for (int i = 0; i < NbrSpins; ++i)
		{
		  InitialConfiguration |= ((ULONGLONG) 0x2ul) << (i << 1);
		}
	    }
	  else
	    {
	      for (int i = 0; i < (NbrSpins / PValue); i += 2)
		{
		  for (int j = 0; j < PValue; ++j)
		    {
		      InitialConfiguration |= ((ULONGLONG) 0x3ul) << (((i * PValue) + j) << 1);
		    }
		}
	    }
	  InitialStateIndex = ((PairHoppingP1AsSpin1ChainWithTranslationsLong*) Chain)->FindStateIndex(InitialConfiguration);
	}
      if (InitialStateIndex == Chain->GetHilbertSpaceDimension())
	{
	  cout << "error while retrieving the initial configuration index" << endl;
	  return -1;       
	}
      InitialState[InitialStateIndex] = 1.0;
      if (Manager.GetBoolean("store-initialstate") == true)
	{
	  char* InitialStateOutputFileName = new char [strlen(OutputFileName)+ 64];
	  sprintf (InitialStateOutputFileName, "%s.0.vec", OutputFileName);
	  InitialState.WriteVector(InitialStateOutputFileName);
	  delete[] InitialStateOutputFileName;	  
	}      
      ComplexMatrix AllStates;
      if (Manager.GetBoolean("store-states") == true)
	{
	  AllStates = ComplexMatrix (Chain->GetHilbertSpaceDimension(), NbrTimeSteps + 1);
	  AllStates[0].Copy(InitialState);
	}

      ComplexVector TmpState1 (Chain->GetLargeHilbertSpaceDimension(), true);
      ComplexVector TmpState2 (Chain->GetLargeHilbertSpaceDimension(), true);
      ComplexVector TmpTotalState (Chain->GetLargeHilbertSpaceDimension());
      TmpTotalState.Copy(InitialState);
      double CurrentTime = 0.0;
      ComplexMatrix PartialEntanglementMatrix;
      gettimeofday (&(StartingTime), 0);
      double* TmpValues;
      double ReducedDensityMatrixTrace = 0.0;
      double EntanglementEntropy = 0.0;
      if (Manager.GetBoolean("disable-entropy") == false)
	{
	  for (int TmpSector = 0; TmpSector < ((PValue + 1) * (PValue + 1)); ++TmpSector)
	    {
	      PartialEntanglementMatrix = Chain->EvaluatePartialEntanglementMatrix(SubsystemSize, TmpSector, TmpTotalState);
	      cout << "entanglement matrix size = " << PartialEntanglementMatrix.GetNbrRow() << "x" << PartialEntanglementMatrix.GetNbrColumn() << endl;
	      int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
	      if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
		{
		  TmpDimension = PartialEntanglementMatrix.GetNbrRow();
		}
	      TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
	      for (int i = 0; i < TmpDimension; ++i)
		{
		  TmpValues[i] *= TmpValues[i];
		  ReducedDensityMatrixTrace += TmpValues[i];
		  if (TmpValues[i] > 0.0)
		    {
		      EntanglementEntropy -=  TmpValues[i] * log(TmpValues[i]);
		    }
		}
	      delete[] TmpValues;
	    }
	}
      gettimeofday (&(EndingTime), 0);
      double DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			   ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
      cout << CurrentTime << " " << Norm(TmpTotalState * InitialState) << " " <<  EntanglementEntropy << " " << ReducedDensityMatrixTrace << " " << 0 <<  " " <<  TmpTotalState.Norm() << " " << DeltaTime << endl;
      File << CurrentTime << " " << Norm(TmpTotalState * InitialState) << " " <<  EntanglementEntropy << " " << ReducedDensityMatrixTrace << " " << 0 <<  " " <<   TmpTotalState.Norm() << " " << DeltaTime << endl;
      CurrentTime += TimeStep;
      for (int Index = 1; Index <= NbrTimeSteps; ++Index)
	{
	  double TmpNorm = 0.0;
	  TmpState1.Copy(TmpTotalState);
	  //TmpTotalState.ClearVector();
	  int NbrIterations;
	  gettimeofday (&(StartingTime), 0);
	  for (NbrIterations = 1; (NbrIterations < MaxNbrIterations) && (fabs(TmpNorm - 1.0) > NormConvergence); ++NbrIterations)
	    {
	      VectorHamiltonianMultiplyOperation TmpOperation (Hamiltonian, &TmpState1, &TmpState2);
	      TmpOperation.ApplyOperation(Architecture.GetArchitecture());
	      TmpState2 *= Complex (0, TimeStep / ((double) NbrIterations));
	      TmpTotalState += TmpState2;
	      ComplexVector TmpState3 = TmpState2;
	      TmpState2 = TmpState1;
	      TmpState1 = TmpState3;
	      TmpNorm = TmpTotalState.Norm();
	      //cout << CurrentTime << " " << NbrIterations << " " << TmpNorm << " " <<  NbrIterations <<  " " << << endl;
	    }
	  ReducedDensityMatrixTrace = 0.0;
	  EntanglementEntropy = 0.0;
	  if (Manager.GetBoolean("disable-entropy") == false)
	    {
	      for (int TmpSector = 0; TmpSector < ((PValue + 1) * (PValue + 1)); ++TmpSector)
		{
		  PartialEntanglementMatrix = Chain->EvaluatePartialEntanglementMatrix(SubsystemSize, TmpSector, TmpTotalState);
		  int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
		  if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
		    {
		      TmpDimension = PartialEntanglementMatrix.GetNbrRow();
		    }
		  TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
		  for (int i = 0; i < TmpDimension; ++i)
		    {
		      TmpValues[i] *= TmpValues[i];
		      ReducedDensityMatrixTrace += TmpValues[i];
		      if (TmpValues[i] > 0.0)
			{
			  EntanglementEntropy -=  TmpValues[i] * log(TmpValues[i]);
			}
		    }
		  delete[] TmpValues;
		}
	    }
          gettimeofday (&(EndingTime), 0);
          double DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
                               ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	  cout << CurrentTime << " " << SqrNorm(TmpTotalState * InitialState) << " " <<  EntanglementEntropy << " " << ReducedDensityMatrixTrace << " " << NbrIterations <<  " " << TmpNorm << " " << DeltaTime << endl;
	  File << CurrentTime << " " << SqrNorm(TmpTotalState * InitialState) << " " <<  EntanglementEntropy << " " << ReducedDensityMatrixTrace << " " << NbrIterations <<  " " << TmpNorm << " " << DeltaTime << endl;
	  CurrentTime += TimeStep;
	  TmpTotalState /= TmpTotalState.Norm();
	  if (Manager.GetBoolean("store-states") == true)
	    {
	      AllStates[Index].Copy(TmpTotalState);
	    }
	}
      if (Manager.GetBoolean("store-states") == true)
	{
	  AllStates.WriteMatrix(StateOutputFileName);
	}
      delete Hamiltonian;
    }  
  File.close();
  return 0;
}
