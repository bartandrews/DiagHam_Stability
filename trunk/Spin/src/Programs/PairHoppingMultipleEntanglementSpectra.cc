#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "HilbertSpace/PairHoppingP1AsSpin1Chain.h"
#include "HilbertSpace/PairHoppingP2AsSpin1Chain.h"
#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslations.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainWithTranslations.h"
#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetry.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Architecture/ArchitectureOperation/SpinChainMultipleEntanglementSpectrumOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "Options/Options.h"

#include "GeneralTools/Endian.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;



int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("PairHoppingMultipleEntanglementSpectra" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
 
  (*SystemGroup) += new  SingleIntegerOption ('p', "p-value", "value that defines the filling factor p/(2p+1)", 1);
  (*SystemGroup) += new SingleStringOption  ('\n', "multiple-states", "provide as a matrix a series of states whose entanglement spectrum/entropy have to be computed");  
  (*SystemGroup) += new SingleIntegerOption ('\n', "min-multiplestates", "index of the first state to consider in the matrix provided by the --multiple-groundstate option", 0);  
  (*SystemGroup) += new SingleIntegerOption ('\n', "max-multiplestates", "index of the last state to consider in the matrix provided by the --multiple-groundstate option (negative if this is the last available state)", -1);  
  (*SystemGroup) += new BooleanOption  ('c', "complex", "consider complex wave function");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "la", "subsystem size (negative if half of the system has to be considered)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");  
  (*OutputGroup) += new BooleanOption ('\n', "show-entropies", "show the entangelement entropy and trace of the reduced density matrix for each state");
  (*OutputGroup) += new BooleanOption ('\n', "export-entropies", "export in a text file the entangelement entropy and trace of the reduced density matrix for each state");
  (*OutputGroup) += new BooleanOption ('\n', "counting", "count the number of non zero eigenvalues for each entanglement spectrum");
  (*OutputGroup) += new SingleDoubleOption ('\n', "counting-error", "error below which an eigenvalue is consideredto be eual to zero", 1e-28);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PairHoppingMultipleEntanglementSpectra -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("multiple-states") == 0)
    {
      cout << "error, an eigenstate file should be provided. See man page for option syntax or type PairHoppingMultipleEntanglementSpectra -h" << endl;
      return -1;
    }

  int PValue = Manager.GetInteger("p-value");
  int NbrSpins = 0;
  int SpinValue = 0;
  int SubsystemSize = Manager.GetInteger("la");
  bool InversionFlag = false;
  bool Momentum1DFlag = false;
  int XMomentum = 0;
  int XPeriodicity = 0;
  int InversionSector = 0;
  int SzSymmetrySector = 0;
  int* SubsystemSites = 0;
  bool ShowTimeFlag = Manager.GetBoolean("show-time");

  char* StrPValue = strstr(Manager.GetString("multiple-states"), "_p_");
  if (StrPValue != 0)
    {
      StrPValue += 3;
      int SizeString = 0;
      while ((StrPValue[SizeString] != '\0') && (StrPValue[SizeString] >= '0') && (StrPValue[SizeString] <= '9'))
	++SizeString;
      if (SizeString != 0)
	{
	  char Tmp = StrPValue[SizeString];
	  StrPValue[SizeString] = '\0';
	  PValue = atoi(StrPValue);
	  StrPValue[SizeString] = Tmp;
	  StrPValue += SizeString;
	}
      else
	StrPValue = 0;
    }
  if (StrPValue == 0)
    {
      cout << "can't guess p value from file name " << Manager.GetString("multiple-states") << endl;
    }

  if (SpinAllSzFindSystemInfoFromVectorFileName(Manager.GetString("multiple-states"), NbrSpins, SpinValue, 
						XMomentum, InversionSector, SzSymmetrySector) == false)
    {
      if (SpinAllSzFindSystemInfoFromVectorFileName(Manager.GetString("multiple-states"), NbrSpins, SpinValue, InversionSector, SzSymmetrySector) == false)
	{
	  if (SpinFindSystemInfoFromFileName(Manager.GetString("multiple-states"), NbrSpins, SpinValue) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << Manager.GetString("multiple-states") << endl;
	      return -1;
	    }
	}
    }
  else
    {
      Momentum1DFlag = true;
    }

  cout << "N=" << NbrSpins << " p=" << PValue;
  if (InversionSector != 0)
    cout << " IP=" << InversionSector;
  cout << endl;

  if (SubsystemSize < 0)
    {
      SubsystemSize = NbrSpins / 2;
    }



  char* TmpExtenstion = new char[32];
  sprintf(TmpExtenstion, "_la_%d.full.ent", SubsystemSize);
  char* DensityMatrixFileName = ReplaceExtensionToFileName(Manager.GetString("multiple-states"), ".mat", TmpExtenstion);
  sprintf(TmpExtenstion, "_la_%d.ent", SubsystemSize);
  char* EntropyFileName = ReplaceExtensionToFileName(Manager.GetString("multiple-states"), ".mat", TmpExtenstion);

  int MinIndex = Manager.GetInteger("min-multiplestates");
  int MaxIndex = Manager.GetInteger("max-multiplestates");
  RealMatrix RealEigenstates;
  ComplexMatrix ComplexEigenstates;

  if (Manager.GetBoolean("complex") == false)
    {
      if (RealEigenstates.ReadMatrix(Manager.GetString("multiple-states")) == false)
	{
	  cout << "can't read " << Manager.GetString("multiple-states") << endl;
	  return -1;
	}
      if (MaxIndex < 0)
	MaxIndex = RealEigenstates.GetNbrColumn() - 1;
    }
  else
    {
      if (ComplexEigenstates.ReadMatrix(Manager.GetString("multiple-states")) == false)
	{
	  cout << "can't read " << Manager.GetString("multiple-states") << endl;
	  return -1;
	}
      if (MaxIndex < 0)
	MaxIndex = ComplexEigenstates.GetNbrColumn() - 1;
    }
  
  int NbrStates = (MaxIndex - MinIndex + 1);

  AbstractSpinChain* Space;

  if (Momentum1DFlag == false)
    {
      switch (PValue)
	{
	case 1 :
	  Space = new PairHoppingP1AsSpin1Chain (NbrSpins, false, 1000000);
	  break;
	case 2 :
	  Space = new PairHoppingP2AsSpin1Chain (NbrSpins, false, 1000000);
	  break;
	default :
	  {
	    cout << "p value > 2 are not available" << endl;
	    return -1;
	  }
	}	  
     }
  else
    {
      switch (PValue)
	{
	case 1 :
	  {
	    if (InversionSector != 0)
	      {
		Space = new PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry (NbrSpins, XMomentum, InversionSector);
	      }
	    else
	      {
		Space = new PairHoppingP1AsSpin1ChainWithTranslations (NbrSpins, XMomentum);
	      }
	  }
	  break;
	case 2 :
	  {
	    if (InversionSector != 0)
	      {
		Space = new PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetry (NbrSpins, XMomentum, InversionSector);
	      }
	    else
	      {
		Space = new PairHoppingP2AsSpin1ChainWithTranslations (NbrSpins, XMomentum);
	      }
	  }
	  break;
	default :
	  {
	    cout << "p value > 2 are not available" << endl;
	    return -1;
	  }
	}
    }      
 
  int MaxSzA = (PValue + 1) * (PValue + 1) - 1;
  int MinSzA = 0;
  double*** EntanglementSpectra = new double**[MaxSzA + 1];
  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; ++TmpSzA)
    {
      EntanglementSpectra[TmpSzA] = new double*[NbrStates];
    }
  int* EntanglementSpectrumDimension = new int[MaxSzA + 1];

  timeval TotalStartingTime;
  timeval TotalEndingTime;
  if (Manager.GetBoolean("complex") == false)
    {
      for (int TmpSzA = 0; TmpSzA <= MaxSzA; ++TmpSzA)
	{
	  cout << "processing subsytem size " << SubsystemSize << "  sector=(" << (TmpSzA / (PValue + 1)) << "," << (TmpSzA % (PValue + 1)) << ")" << endl;
	  SpinChainMultipleEntanglementSpectrumOperation TmpOperation(Space, RealEigenstates, MinIndex, MaxIndex, SubsystemSize, TmpSzA, SubsystemSites);
	  TmpOperation.ApplyOperation(Architecture.GetArchitecture());
	  EntanglementSpectrumDimension[TmpSzA] = TmpOperation.GetEntanglementSpectrumDimension();	      
	  EntanglementSpectra[TmpSzA] = TmpOperation.GetEntanglementSpectra();
	}
    }
  else
    {
      for (int TmpSzA = 0; TmpSzA <= MaxSzA; ++TmpSzA)
	{
	  cout << "processing subsytem size " << SubsystemSize << "  sector=(" << (TmpSzA / (PValue + 1)) << "," << (TmpSzA % (PValue + 1)) << ")" << endl;
	  SpinChainMultipleEntanglementSpectrumOperation TmpOperation(Space, ComplexEigenstates, MinIndex, MaxIndex, SubsystemSize, TmpSzA, SubsystemSites);
	  TmpOperation.ApplyOperation(Architecture.GetArchitecture());
	  EntanglementSpectrumDimension[TmpSzA] = TmpOperation.GetEntanglementSpectrumDimension();	      
	  EntanglementSpectra[TmpSzA] = TmpOperation.GetEntanglementSpectra();
	}
    }
  int TotalDimensionPerState = 0;
  int NbrSzASectors = 0;
  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; ++TmpSzA)
    {
      TotalDimensionPerState += EntanglementSpectrumDimension[TmpSzA];
      if (EntanglementSpectrumDimension[TmpSzA] > 0)
	{
	  NbrSzASectors++;
	}
    }
  ofstream File;
  File.open(DensityMatrixFileName, ios::binary | ios::out);
  WriteLittleEndian(File, NbrStates);  
  WriteLittleEndian(File, NbrSzASectors);  
  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; ++TmpSzA)
    {
      if (EntanglementSpectrumDimension[TmpSzA] > 0)
	{
	  WriteLittleEndian(File, TmpSzA);  
	}
    }
  
  if (Manager.GetBoolean("show-entropies"))
    {
      cout << "# state_index entropy trace" << endl;
    }
  if (Manager.GetBoolean("export-entropies"))
    {
      ofstream FileEntropy;
      FileEntropy.open(EntropyFileName, ios::binary | ios::out);
      FileEntropy << "# Index S_A Tr_A(rho_A)" << endl;
      FileEntropy.close();
    }
  for (int Index = MinIndex; Index <= MaxIndex; ++Index)
    {
      double TmpEntanglementEntropy = 0.0;
      double TmpTrace = 0.0;
      for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; ++TmpSzA)
	{
	  if (EntanglementSpectrumDimension[TmpSzA] > 0)
	    {
	      WriteLittleEndian(File, EntanglementSpectrumDimension[TmpSzA]);
	      WriteBlockLittleEndian(File, EntanglementSpectra[TmpSzA][Index - MinIndex], EntanglementSpectrumDimension[TmpSzA]);
	      for (int i = 0; i < EntanglementSpectrumDimension[TmpSzA]; ++i)
		{
		  double Tmp = EntanglementSpectra[TmpSzA][Index - MinIndex][i];
		  if (Tmp > 0.0)
		    TmpEntanglementEntropy -= Tmp * log(Tmp);
		  TmpTrace +=  Tmp;
		}
	    }
	}
      if (Manager.GetBoolean("show-entropies"))
	{
	  cout << Index << " " << TmpEntanglementEntropy << " " << TmpTrace << endl;
	}
      if (Manager.GetBoolean("export-entropies"))
	{
	  ofstream FileEntropy;
	  FileEntropy.open(EntropyFileName, ios::binary | ios::out | ios::app);
	  FileEntropy.precision(14);
	  FileEntropy << Index << " " << TmpEntanglementEntropy << " " << TmpTrace << endl;
	  FileEntropy.close();
	}
    }
  File.close();

  if (Manager.GetBoolean("counting"))
    {
      char* DensityMatrixCountingFileName = new char[strlen(DensityMatrixFileName) + 32];
      sprintf (DensityMatrixCountingFileName, "%s.counting", DensityMatrixFileName);
      ofstream File;
      File.open(DensityMatrixCountingFileName, ios::binary | ios::out);
      File << "# index total_counting total_nbr";
      for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; ++TmpSzA)
	{
	  File << " counting(Sza=" << TmpSzA << ")";
	}
      File << endl;
      int* TmpCounting = new int [((MaxSzA - MinSzA) / 2) + 1];
      double CountingError = Manager.GetDouble("counting-error");
      for (int Index = MinIndex; Index <= MaxIndex; ++Index)
	{
	  long TotalNbr = 0;
	  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; ++TmpSzA)
	    {
	      TmpCounting[TmpSzA] = 0;
	      for (int i = 0; i < EntanglementSpectrumDimension[TmpSzA]; ++i)
		{
		  ++TotalNbr;
		  if (EntanglementSpectra[TmpSzA][Index - MinIndex][i] > CountingError)
		    ++TmpCounting[TmpSzA];
		}
	    }		  
	  int TmpTotalCounting = 0;
	  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; ++TmpSzA)
	    {
	      TmpTotalCounting += TmpCounting[TmpSzA];
	    }
	  File << Index << " " << TmpTotalCounting << " " << TotalNbr;
	  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; ++TmpSzA)
	    {
	      File << " " << TmpCounting[TmpSzA];
	    }	  
	  File << endl;
	}
      File.close();
    }
  return 0;
}


