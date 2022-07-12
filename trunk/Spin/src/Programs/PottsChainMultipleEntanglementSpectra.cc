#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Hamiltonian/SpinChainHamiltonian.h"

#include "HilbertSpace/Potts3Chain.h"
#include "HilbertSpace/Potts3ChainWithTranslations.h"
#include "HilbertSpace/Potts3ChainWithTranslationsAndInversion.h"


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
  OptionManager Manager ("PottsChainMultipleEntanglementSpectra" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption ('n', "zn", "Zn value of the Potts model", 3);
  (*SystemGroup) += new SingleStringOption  ('\n', "multiple-states", "provide as a matrix a series of states whose entanglement spectrum/entropy have to be computed");  
  (*SystemGroup) += new SingleIntegerOption ('\n', "min-multiplestates", "index of the first state to consider in the matrix provided by the --multiple-groundstate option", 0);  
  (*SystemGroup) += new SingleIntegerOption ('\n', "max-multiplestates", "index of the last state to consider in the matrix provided by the --multiple-groundstate option (negative if this is the last available state)", -1);  
  (*SystemGroup) += new BooleanOption  ('c', "complex", "consider complex wave function");
  (*SystemGroup) += new SingleStringOption  ('\n', "generic-cut", "provide a list of sites that define the subsystem instead of the default cut");  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "la", "subsystem size (negative if half of the system has to be considered)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");  
  (*OutputGroup) += new BooleanOption ('\n', "show-entropies", "show the entangelement entropy and trace of the reduced density matrix for each state");
  (*OutputGroup) += new BooleanOption ('\n', "export-entropies", "export in a text file the entangelement entropy and trace of the reduced density matrix for each state");
  (*OutputGroup) += new BooleanOption ('\n', "counting", "count the number of non zero eigenvalues for each entanglement spectrum");
  (*OutputGroup) += new SingleDoubleOption ('\n', "counting-error", "error below which an eigenvalue is consideredto be eual to zero", 1e-28);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PottsChainMultipleEntanglementSpectra -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("multiple-states") == 0)
    {
      cout << "error, an eigenstate file should be provided. See man page for option syntax or type PottsChainMultipleEntanglementSpectra -h" << endl;
      return -1;
    }

  int ZnValue = Manager.GetInteger("zn");
  int NbrSpins = 0;
  int TotalQ = 0;
  int SubsystemSize = Manager.GetInteger("la");
  bool SzFlag = true;
  bool Momentum1DFlag = false;
  bool InversionFlag = false;  
  int XMomentum = 0;
  int XPeriodicity = 0;
  int InversionSector = 0;
  bool GenericCutFlag = false;
  int* SubsystemSites = 0;
  bool ShowTimeFlag = Manager.GetBoolean("show-time");

  if (PottsFindSystemInfoFromVectorFileName(Manager.GetString("multiple-states"), NbrSpins, ZnValue, TotalQ, 
					    XMomentum, InversionSector) == false)
    {
      if (PottsFindSystemInfoFromVectorFileName(Manager.GetString("multiple-states"), NbrSpins, ZnValue, TotalQ) == false)
	{
	  SzFlag = false;
	  //	  if (SpinFindSystemInfoFromFileName(Manager.GetString("multiple-states"), NbrSpins, ZnValue) == false)
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

  if (SzFlag == true)
    cout << "N=" << NbrSpins << " Q=" <<  TotalQ << " Zn=" << ZnValue;
  else
    cout << "N=" << NbrSpins << " Zn=" << ZnValue;
  if (Momentum1DFlag == true)
    {
      cout << " k=" << XMomentum;
    }
  if (InversionSector != 0)
    cout << " I=" << InversionSector;
  cout << endl;


  if (SubsystemSize < 0)
    {
      SubsystemSize = NbrSpins / 2;
    }

  // if (Manager.GetString("generic-cut") != 0)
  //   {
  //     GenericCutFlag = true;
  //     MultiColumnASCIIFile GenericCutFile;
  //     if (GenericCutFile.Parse(Manager.GetString("generic-cut")) == false)
  // 	{
  // 	  GenericCutFile.DumpErrors(cout);
  // 	  return -1;
  // 	}
  //     SubsystemSize = GenericCutFile.GetNbrLines();
  //     if (GenericCutFile.GetNbrColumns() == 1)
  // 	{
  // 	  SubsystemSites = GenericCutFile.GetAsIntegerArray(0);
  // 	}
  //     else
  // 	{
  // 	  SubsystemSites =  new int [SubsystemSize];
  // 	  int* TmpXPositions = GenericCutFile.GetAsIntegerArray(0);
  // 	  int* TmpYPositions = GenericCutFile.GetAsIntegerArray(1);
  // 	  for (int i = 0; i < SubsystemSize; ++i)
  // 	    {
  // 	      SubsystemSites[i] = PottsChainMultipleEntanglementSpectraGetLinearizedIndex(TmpXPositions[i], XPeriodicity, TmpYPositions[i], YPeriodicity);
  // 	    }
  // 	}
  //   }
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

  if (SzFlag == true)
    {
      if (Momentum1DFlag == false)
	{
	  switch (ZnValue)
	    {
	    case 3:
	      {
		if (InversionSector != 0)
		  {
		    //		    Space = new Potts3Chain (NbrSpins, InversionSector, TotalQ, 1000000);
		  }
		else
		  {
		    Space = new Potts3Chain (NbrSpins, TotalQ, 1000000);
		  }
	      }
	      break;
	    default :
	      {
		  cout << "Potts " << ZnValue << " are not available" << endl;
		return -1;
	      }
	    }
	}
      else
	{
	  switch (ZnValue)
	    {
	    case 3:
	      {
		if (InversionSector != 0)
		  {
		    Space = new Potts3ChainWithTranslationsAndInversion (NbrSpins, TotalQ, XMomentum, InversionSector, 1000000);
		  }
		else
		  {
		    Space = new Potts3ChainWithTranslations (NbrSpins, TotalQ, XMomentum, 1000000);
		  }
	      }
	      break;
	    default :
	      {
		cout << "Potts " << ZnValue << " are not available" << endl;
		return -1;
	      }
	    }
	}
    }
  else
    {
      switch (ZnValue)
	{
	default :
	  {
	    cout << "Potts " << ZnValue << " are not available" << endl;
	    return -1;
	  }
	}
    }
 
  int MaxQA = ZnValue - 1;
  int MinQA = 0;
  int MaxQB = ZnValue - 1;
  int MinQB = 0;
  if (SzFlag == false)
    {
      MaxQA = 0;
      MinQA = 0;
      MaxQB = 0;
      MinQB =0;
      TotalQ = 0;
    }
  double*** EntanglementSpectra = new double**[(MaxQA - MinQA) + 1];
  for (int TmpQA = MinQA;  TmpQA <= MaxQA; ++TmpQA)
    {
      EntanglementSpectra[TmpQA - MinQA] = new double*[NbrStates];
    }
  int* EntanglementSpectrumDimension = new int[(MaxQA - MinQA) + 1];

  timeval TotalStartingTime;
  timeval TotalEndingTime;
  if (Manager.GetBoolean("complex") == false)
    {
      for (int TmpQA = MinQA; TmpQA <= MaxQA; ++TmpQA)
	{
	  int QB = TotalQ - TmpQA;
	  if (QB < 0)
	    QB += ZnValue;
	  if ((QB <= MaxQB) && (QB >= MinQB))
	    {
	      if (SzFlag == true)
		{
		  cout << "processing subsytem size " << SubsystemSize << " QA=" << TmpQA << endl;
		}
	      else
		{
		  cout << "processing subsytem size " << SubsystemSize << endl;
		}

	      SpinChainMultipleEntanglementSpectrumOperation TmpOperation(Space, RealEigenstates, MinIndex, MaxIndex, SubsystemSize, TmpQA, SubsystemSites);
	      TmpOperation.ApplyOperation(Architecture.GetArchitecture());
	      EntanglementSpectrumDimension[TmpQA - MinQA] = TmpOperation.GetEntanglementSpectrumDimension();	      
	      EntanglementSpectra[TmpQA - MinQA] = TmpOperation.GetEntanglementSpectra();
	    }
	  else
	    {
	      EntanglementSpectrumDimension[TmpQA - MinQA] = 0;
	    }  
	}
    }
  else
    {
      for (int TmpQA = MinQA; TmpQA <= MaxQA; ++TmpQA)
	{
	  int QB = TotalQ - TmpQA;
	  if (QB < 0)
	    QB += ZnValue;
	  if ((QB <= MaxQB) && (QB >= MinQB))
	    {
	      if (SzFlag == true)
		{
		  cout << "processing subsytem size " << SubsystemSize << " QA=" << TmpQA << endl;
		}
	      else
		{
		  cout << "processing subsytem size " << SubsystemSize << endl;
		}
	      SpinChainMultipleEntanglementSpectrumOperation TmpOperation(Space, ComplexEigenstates, MinIndex, MaxIndex, SubsystemSize, TmpQA, SubsystemSites);
	      TmpOperation.ApplyOperation(Architecture.GetArchitecture());
	      EntanglementSpectrumDimension[TmpQA - MinQA] = TmpOperation.GetEntanglementSpectrumDimension();	      
	      EntanglementSpectra[TmpQA - MinQA] = TmpOperation.GetEntanglementSpectra();
	    }
	  else
	    {
	      EntanglementSpectrumDimension[TmpQA - MinQA] = 0;
	    }  
	}
    }
  
  int TotalDimensionPerState = 0;
  int NbrSzASectors = 0;
  for (int TmpQA = MinQA;  TmpQA <= MaxQA; ++TmpQA)
    {
      TotalDimensionPerState += EntanglementSpectrumDimension[TmpQA - MinQA];
      if (EntanglementSpectrumDimension[TmpQA - MinQA] > 0)
	{
	  NbrSzASectors++;
	}
    }
  ofstream File;
  File.open(DensityMatrixFileName, ios::binary | ios::out);
  WriteLittleEndian(File, NbrStates);  
  WriteLittleEndian(File, NbrSzASectors);  
  for (int TmpQA = MinQA;  TmpQA <= MaxQA; ++TmpQA)
    {
      if (EntanglementSpectrumDimension[TmpQA - MinQA] > 0)
	{
	  WriteLittleEndian(File, TmpQA);  
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
      for (int TmpQA = MinQA;  TmpQA <= MaxQA; ++TmpQA)
	{
	  if (EntanglementSpectrumDimension[TmpQA - MinQA] > 0)
	    {
	      WriteLittleEndian(File, EntanglementSpectrumDimension[TmpQA - MinQA]);
	      WriteBlockLittleEndian(File, EntanglementSpectra[TmpQA - MinQA][Index - MinIndex], EntanglementSpectrumDimension[TmpQA - MinQA]);
	      for (int i = 0; i < EntanglementSpectrumDimension[TmpQA - MinQA]; ++i)
		{
		  double Tmp = EntanglementSpectra[TmpQA - MinQA][Index - MinIndex][i];
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
      for (int TmpQA = MinQA;  TmpQA <= MaxQA; ++TmpQA)
	{
	  File << " counting(Qa=" << TmpQA << ")";
	}
      File << endl;
      int* TmpCounting = new int [(MaxQA - MinQA) + 1];
      double CountingError = Manager.GetDouble("counting-error");
      for (int Index = MinIndex; Index <= MaxIndex; ++Index)
	{
	  long TotalNbr = 0;
	  for (int TmpQA = MinQA;  TmpQA <= MaxQA; ++TmpQA)
	    {
	      TmpCounting[TmpQA - MinQA] = 0;
	      for (int i = 0; i < EntanglementSpectrumDimension[TmpQA - MinQA]; ++i)
		{
		  ++TotalNbr;
		  if (EntanglementSpectra[TmpQA - MinQA][Index - MinIndex][i] > CountingError)
		    ++TmpCounting[TmpQA - MinQA];
		}
	    }		  
	  int TmpTotalCounting = 0;
	  for (int TmpQA = MinQA;  TmpQA <= MaxQA; ++TmpQA)
	    {
	      TmpTotalCounting += TmpCounting[TmpQA - MinQA];
	    }
	  File << Index << " " << TmpTotalCounting << " " << TotalNbr;
	  for (int TmpQA = MinQA;  TmpQA <= MaxQA; ++TmpQA)
	    {
	      File << " " << TmpCounting[TmpQA - MinQA];
	    }	  
	  File << endl;
	}
      File.close();
    }
  return 0;
}
