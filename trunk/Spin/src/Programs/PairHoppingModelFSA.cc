#include "HilbertSpace/PairHoppingP1AsSpin1Chain.h"
#include "HilbertSpace/PairHoppingP2AsSpin1Chain.h"
#include "HilbertSpace/PairHoppingGenericPAsSpin1Chain.h"

#include "HilbertSpace/PairHoppingP1AsSpin1ChainLong.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainLong.h"
#include "HilbertSpace/PairHoppingGenericPAsSpin1ChainLong.h"

#include "Hamiltonian/PairHoppingHamiltonian.h"

#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Operator/PairHoppingHPlusOperator.h"

#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

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
  OptionManager Manager ("PairHoppingModelFSA" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new  SingleIntegerOption ('p', "p-value", "value that defines the filling factor p/(2p+1)", 1);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-spin", "number of spins", 10);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PairHoppingModelFSA -h" << endl;
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

  char* OutputFileName = new char [512];
  sprintf (OutputFileName, "spin_1_pairhopping_p_%d_n_%d", PValue, NbrSpins);


  char* FullOutputFileName = new char [strlen(OutputFileName)+ 64];
  sprintf (FullOutputFileName, "%s.fsa", OutputFileName);

  ofstream File;
  File.open(FullOutputFileName, ios::binary | ios::out);
  File.precision(14);
  File << "# periodic pair hopping model FSA with p=" << PValue << " in spin 1 language with " << NbrSpins << " sites" << endl;
  File << "# FSA_energy" << endl;
  AbstractSpinChain* Chain = 0;
  if (NbrSpins <= 32)
    {
      switch (PValue)
	{
	case 1 :
	  Chain = new PairHoppingP1AsSpin1Chain (NbrSpins, true, 1000000);
	  break;
	case 2 :
	  Chain = new PairHoppingP2AsSpin1Chain (NbrSpins, true, 1000000);
	  break;
	default :
	  {
	    Chain = new PairHoppingGenericPAsSpin1Chain (NbrSpins, PValue, true, 1000000);
	  }
	}
    }
  else
    {
      switch (PValue)
	{
	case 1 :
	  Chain = new PairHoppingP1AsSpin1ChainLong (NbrSpins, true, 1000000);
	  break;
	case 2 :
	  Chain = new PairHoppingP2AsSpin1ChainLong (NbrSpins, true, 1000000);
	  break;
	default :
	  {
	    Chain = new PairHoppingGenericPAsSpin1ChainLong (NbrSpins, PValue, true, 1000000);
	  }
	}
    }
  cout << "Hilbert space dimension = " <<  Chain->GetHilbertSpaceDimension() << endl;
  return 0;
  
  if (Chain->GetHilbertSpaceDimension() > 0)
    {
      Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());
      int NbrUnitCells = (NbrSpins / PValue);
      //     int NbrStates = NbrUnitCells  + 1;
      int NbrStates = (NbrSpins *  PValue) + 1;
      
      RealVector* HPlusNStates = new RealVector[NbrStates];
      HPlusNStates[0] = RealVector(Chain->GetHilbertSpaceDimension(), true);
      if (NbrSpins <= 32)
	{
	  int InitialStateIndex = 0;
	  unsigned long InitialConfiguration = 0x0ul;
	  for (int i = 0; i < (NbrSpins / PValue); i += 2)
	    {
	      for (int j = 0; j < PValue; ++j)
	       	{
		   InitialConfiguration |= (0x3ul) << (((i * PValue) + j) << 1);
		}
	    }
	  InitialStateIndex = ((PairHoppingP1AsSpin1Chain*) Chain)->FindStateIndex(InitialConfiguration);
	  ((PairHoppingP1AsSpin1Chain*) Chain)->PrintState(cout, InitialStateIndex) << endl;
	  HPlusNStates[0][InitialStateIndex] = 1.0;// / M_SQRT2;
	  InitialConfiguration = 0x0ul;
	  for (int i = 1; i < (NbrSpins / PValue); i += 2)
	    {
	      for (int j = 0; j < PValue; ++j)
	       	{
		   InitialConfiguration |= (0x3ul) << (((i * PValue) + j) << 1);
		}
	    }
	  InitialStateIndex = ((PairHoppingP1AsSpin1Chain*) Chain)->FindStateIndex(InitialConfiguration);
	  ((PairHoppingP1AsSpin1Chain*) Chain)->PrintState(cout, InitialStateIndex) << endl;
	  // HPlusNStates[0][InitialStateIndex] = 1.0;// / M_SQRT2;
	}
      else
	{
	  int InitialStateIndex = 0;
	  ULONGLONG InitialConfiguration = ((ULONGLONG) 0x0ul);
	  for (int i = 0; i < (NbrSpins / PValue); i += 2)
	    {
	      for (int j = 0; j < PValue; ++j)
		{
		  InitialConfiguration |= ((ULONGLONG) 0x3ul) << (((i * PValue) + j) << 1);
		}
	    }
	  InitialStateIndex = ((PairHoppingP1AsSpin1ChainLong*) Chain)->FindStateIndex(InitialConfiguration);
	  HPlusNStates[0][InitialStateIndex] = 1.0;// / M_SQRT2;
	  InitialConfiguration <<= (PValue << 1);
	  InitialStateIndex = ((PairHoppingP1AsSpin1ChainLong*) Chain)->FindStateIndex(InitialConfiguration);
	  //	  HPlusNStates[0][InitialStateIndex] = 1.0; //1.0 / M_SQRT2;
	}
      PairHoppingHPlusOperator Operator (Chain, NbrSpins, PValue, true);
      RealSymmetricMatrix TmpHamiltonian (NbrStates, true);
      RealDiagonalMatrix TmpDiag (NbrStates);
      for (int i = 1; i < NbrStates; ++i)
	{
	  HPlusNStates[i] = RealVector(Chain->GetHilbertSpaceDimension(), true);
	  VectorOperatorMultiplyOperation TmpOperation(&Operator, &(HPlusNStates[i - 1]), &(HPlusNStates[i]));
	  TmpOperation.ApplyOperation(Architecture.GetArchitecture());
	  double TmpNorm = HPlusNStates[i].Norm();
	  cout << "beta_" << i << " : " << TmpNorm << endl;
 	  TmpHamiltonian.SetMatrixElement(i - 1, i, TmpNorm);
	  HPlusNStates[i] /= TmpNorm;
	}
     TmpHamiltonian.LapackDiagonalize(TmpDiag);
      File.precision(14);
      for (int i = 0; i < NbrStates; ++i)
	{
	  File << TmpDiag[i] << endl;
	}
    }  
  File.close();
  return 0;
}
