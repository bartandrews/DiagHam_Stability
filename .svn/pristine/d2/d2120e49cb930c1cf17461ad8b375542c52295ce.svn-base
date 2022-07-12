#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Hamiltonian/SpinChainHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin2Chain.h"
#include "HilbertSpace/Spin1_2ChainFull.h"
#include "HilbertSpace/Spin1_2ChainFullAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainFullInversionAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzInversionSymmetries.h"
#include "HilbertSpace/Spin1_2ChainLong.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsLong.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzSymmetryLong.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndInversionSymmetryLong.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzInversionSymmetriesLong.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzInversionSymmetries.h"
#include "HilbertSpace/Spin2ChainWithTranslations.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndSzInversionSymmetries.h"
#include "HilbertSpace/Spin1ChainWithSzSymmetry.h"
#include "HilbertSpace/Spin1ChainWithInversionSymmetry.h"
#include "HilbertSpace/Spin1ChainWithSzInversionSymmetries.h"


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



// get a linearized position index from the 2d coordinates
//
// xPosition = position along the x direction
// xPeriodicity = system size in the x direction
// yPosition = position along the y direction
// yPeriodicity = system size in the y direction
// return value = linearized index
int SpinChainMultipleEntanglementSpectraGetLinearizedIndex(int xPosition, int xPeriodicity, int yPosition, int yPeriodicity);


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("SpinChainMultipleEntanglementSpectra" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 1);
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
      cout << "see man page for option syntax or type SpinChainMultipleEntanglementSpectra -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("multiple-states") == 0)
    {
      cout << "error, an eigenstate file should be provided. See man page for option syntax or type SpinChainMultipleEntanglementSpectra -h" << endl;
      return -1;
    }

  int SpinValue = 0;
  int NbrSpins = 0;
  int TotalSz = 0;
  int SubsystemSize = Manager.GetInteger("la");
  bool SzFlag = true;
  bool Momentum2DFlag = false;
  bool Momentum1DFlag = false;
  bool InversionFlag = false;  
  int XMomentum = 0;
  int XPeriodicity = 0;
  int YMomentum = 0;
  int YPeriodicity = 0;
  int InversionSector = 0;
  int SzSymmetrySector = 0;
  bool GenericCutFlag = false;
  int* SubsystemSites = 0;
  bool ShowTimeFlag = Manager.GetBoolean("show-time");

  if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("multiple-states"), NbrSpins, TotalSz, SpinValue, XMomentum, XPeriodicity,
							    YMomentum, YPeriodicity) == false)
    {
      if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("multiple-states"), NbrSpins, SpinValue, XMomentum, XPeriodicity, 
								YMomentum, YPeriodicity) == false)
	{
	  if (SpinFindSystemInfoFromVectorFileName(Manager.GetString("multiple-states"), NbrSpins, TotalSz, SpinValue, 
						   XMomentum, InversionSector, SzSymmetrySector) == false)
	    {
	      if (SpinFindSystemInfoFromVectorFileName(Manager.GetString("multiple-states"), NbrSpins, TotalSz, SpinValue, InversionSector, SzSymmetrySector) == false)
		{
		  if (SpinFindSystemInfoFromVectorFileName(Manager.GetString("multiple-states"), NbrSpins, TotalSz, SpinValue) == false)
		    {
		      SzFlag = false;
		      if (SpinFindSystemInfoFromFileName(Manager.GetString("multiple-states"), NbrSpins, SpinValue) == false)
			{
			  cout << "error while retrieving system parameters from file name " << Manager.GetString("multiple-states") << endl;
			  return -1;
			}
		    }
		}
	    }
	  else
	    {
	      Momentum1DFlag = true;
	    }
	}
      else
	{
	  Momentum2DFlag = true;
	  SzFlag = false;
	  InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(Manager.GetString("multiple-states"), NbrSpins, SpinValue, XMomentum, XPeriodicity, 
											 YMomentum, YPeriodicity, InversionSector);
	}
    }
  else
    {
      InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(Manager.GetString("multiple-states"), NbrSpins, TotalSz, SpinValue, XMomentum, XPeriodicity,
										     YMomentum, YPeriodicity, InversionSector);
      SzFlag = true;
      Momentum2DFlag = true;
    }
  if (Momentum2DFlag == false)
    {
      if (SzFlag == true)
	cout << "N=" << NbrSpins << " Sz=" <<  TotalSz << " 2s=" << SpinValue;
      else
	cout << "N=" << NbrSpins << " 2s=" << SpinValue;
      if (InversionSector != 0)
	cout << " I=" << InversionSector;
      if (SzSymmetrySector != 0)
	cout << " P=" << SzSymmetrySector;
      cout << endl;
    }
  else
    {
       if (SzFlag == true)
	 {
	   if ((InversionFlag == false) || (InversionSector == 0))
	     {
	       cout << "N=" << NbrSpins << "=" << XPeriodicity << "x" << YPeriodicity << " Sz=" <<  TotalSz << " 2s=" << SpinValue << endl;
	     }
	   else
	     {
	       cout << "N=" << NbrSpins << "=" << XPeriodicity << "x" << YPeriodicity << " Sz=" <<  TotalSz << " 2s=" << SpinValue << " I=" << InversionSector << endl;
	     }
	 }
       else
	 {
	   if ((InversionFlag == false) || (InversionSector == 0))
	     {
	       cout << "N=" << NbrSpins << "=" << XPeriodicity << "x" << YPeriodicity << " 2s=" << SpinValue << endl;
	     }
	   else
	     {
	       cout << "N=" << NbrSpins << "=" << XPeriodicity << "x" << YPeriodicity << " 2s=" << SpinValue << " I=" << InversionSector << endl;
	     }
	 }
   }
  if (SubsystemSize < 0)
    {
      SubsystemSize = NbrSpins / 2;
    }

  if (Manager.GetString("generic-cut") != 0)
    {
      GenericCutFlag = true;
      MultiColumnASCIIFile GenericCutFile;
      if (GenericCutFile.Parse(Manager.GetString("generic-cut")) == false)
	{
	  GenericCutFile.DumpErrors(cout);
	  return -1;
	}
      SubsystemSize = GenericCutFile.GetNbrLines();
      if (GenericCutFile.GetNbrColumns() == 1)
	{
	  SubsystemSites = GenericCutFile.GetAsIntegerArray(0);
	}
      else
	{
	  SubsystemSites =  new int [SubsystemSize];
	  int* TmpXPositions = GenericCutFile.GetAsIntegerArray(0);
	  int* TmpYPositions = GenericCutFile.GetAsIntegerArray(1);
	  for (int i = 0; i < SubsystemSize; ++i)
	    {
	      SubsystemSites[i] = SpinChainMultipleEntanglementSpectraGetLinearizedIndex(TmpXPositions[i], XPeriodicity, TmpYPositions[i], YPeriodicity);
	    }
	}
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

  if ((Manager.GetBoolean("complex") == false) && (Momentum2DFlag == false))
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
	  switch (SpinValue)
	    {
	    case 1 :
	      {
		if (NbrSpins < 64)
		  {
		    Space = new Spin1_2Chain (NbrSpins, TotalSz, 1000000);
		  }
		else
		  {
		    Space = new Spin1_2ChainLong (NbrSpins, TotalSz, 1000000);
		  }
	      }
	      break;
	    case 2 :
	      {
		if (InversionSector != 0)
		  {
		    if (SzSymmetrySector != 0)
		      {
			Space = new Spin1ChainWithSzInversionSymmetries (NbrSpins, InversionSector, SzSymmetrySector, TotalSz, 1000000);
		      }
		    else
		      {
			Space = new Spin1ChainWithInversionSymmetry (NbrSpins, InversionSector, TotalSz, 1000000);
		      }
		  }
		else
		  {
		    if (SzSymmetrySector != 0)
		      {
			Space = new Spin1ChainWithSzSymmetry (NbrSpins, SzSymmetrySector, TotalSz, 1000000);
		      }
		    else
		      {
			Space = new Spin1Chain (NbrSpins, TotalSz, 1000000);
		      }
		  }
	      }
	      break;
	    case 4 :
	      {
		Space = new Spin2Chain (NbrSpins, TotalSz, 1000000);
	      }
	      break;
	    default :
	      {
		if ((SpinValue & 1) == 0)
		  cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		else 
		  cout << "spin " << SpinValue << "/2 are not available" << endl;
		return -1;
	      }
	    }
	}
      else
	{
	  switch (SpinValue)
	    {
	    case 1 :
	      {
		if (InversionSector != 0)
		  {
		    if (SzSymmetrySector != 0)
		      {
			if (NbrSpins < 64)
			  {
			    Space = new Spin1_2ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, XMomentum, 1, InversionSector, SzSymmetrySector, TotalSz, 1000000, 1000000);
			  }
			else
			  {
			    Space = new Spin1_2ChainWithTranslationsAndSzInversionSymmetriesLong (NbrSpins, XMomentum, 1, InversionSector, SzSymmetrySector, TotalSz, 1000000, 1000000);
			  }
		      }
		    else
		      {
			if (NbrSpins < 64)
			  {
			    Space = new Spin1_2ChainWithTranslationsAndInversionSymmetry (NbrSpins, XMomentum, 1, InversionSector, TotalSz, 1000000, 1000000);
			  }
			else
			  {
			    Space = new Spin1_2ChainWithTranslationsAndInversionSymmetryLong (NbrSpins, XMomentum, 1, InversionSector, TotalSz, 1000000, 1000000);
			  }
		      }
		  }
		else
		  {
		    if (SzSymmetrySector != 0)
		      {
			if (NbrSpins < 64)
			  {
			    Space = new Spin1_2ChainWithTranslationsAndSzSymmetry (NbrSpins, XMomentum, 1, SzSymmetrySector, TotalSz, 1000000, 1000000);
			  }
			else
			  {
			    Space = new Spin1_2ChainWithTranslationsAndSzSymmetryLong (NbrSpins, XMomentum, 1, SzSymmetrySector, TotalSz, 1000000, 1000000);
			  }
		      }
		    else
		      {
			if (NbrSpins < 64)
			  {
			    Space = new Spin1_2ChainWithTranslations (NbrSpins, XMomentum, 1, TotalSz, 1000000, 1000000);
			  }
			else
			  {
			    Space = new Spin1_2ChainWithTranslationsLong (NbrSpins, XMomentum, 1, TotalSz, 1000000, 1000000);
			  }
		      }
		  }
	      }
	      break;
	    case 2 :
	      {
		if (InversionSector != 0)
		  {
		    if (SzSymmetrySector != 0)
		      {
			Space = new Spin1ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, XMomentum, InversionSector, SzSymmetrySector, TotalSz);
		      }
		    else
		      {
			Space = new Spin1ChainWithTranslationsAndInversionSymmetry (NbrSpins, XMomentum, InversionSector, TotalSz);
		      }
		  }
		else
		  {
		    if (SzSymmetrySector != 0)
		      {
			Space = new Spin1ChainWithTranslationsAndSzSymmetry (NbrSpins, XMomentum, SzSymmetrySector, TotalSz);
		      }
		    else
		      {
			Space = new Spin1ChainWithTranslations (NbrSpins, XMomentum, TotalSz);
		      }
		  }
	      }
	      break;
	    case 4 :
	      {
		if (InversionSector != 0)
		  {
		    if (SzSymmetrySector != 0)
		      {
			Space = new Spin2ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, XMomentum, InversionSector, SzSymmetrySector, TotalSz);
		      }
		    else
		      {
			Space = new Spin2ChainWithTranslationsAndInversionSymmetry (NbrSpins, XMomentum, InversionSector, TotalSz);
		      }
		  }
		else
		  {
		    if (SzSymmetrySector != 0)
		      {
			Space = new Spin2ChainWithTranslationsAndSzSymmetry (NbrSpins, XMomentum, SzSymmetrySector, TotalSz);
		      }
		    else
		      {
			Space = new Spin2ChainWithTranslations (NbrSpins, XMomentum, TotalSz);
		      }
		  }
	      }
	      break;
	    default :
	      {
		if ((SpinValue & 1) == 0)
		  cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		else 
		  cout << "spin " << SpinValue << "/2 are not available" << endl;
		return -1;
	      }
	    }
	}
    }
  else
    {
      switch (SpinValue)
	{
	case 1 :
	  Space = new Spin1_2ChainFull (NbrSpins);
	  break;
	default :
	  {
	    if ((SpinValue & 1) == 0)
	      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
	    else 
	  cout << "spin " << SpinValue << "/2 are not available" << endl;
	    return -1;
	  }
	}
    }
 
  if (Momentum2DFlag == true)
    {
      AbstractSpinChain* TmpSpace;
      
      if ((InversionFlag == false) || (InversionSector == 0))
	{
	  if (SzFlag == true)
	    {
	      switch (SpinValue)
		{
		default :
		  {
		    if ((SpinValue & 1) == 0)
		      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		    else 
		      cout << "spin " << SpinValue << "/2 are not available" << endl;
		    return -1;
		  }
		}
	    }
	  else
	    {
	      switch (SpinValue)
		{
		case 1 :
		  TmpSpace = new Spin1_2ChainFullAnd2DTranslation (XMomentum, XPeriodicity, YMomentum, YPeriodicity);
		  break;
		default :
		  {
		    if ((SpinValue & 1) == 0)
		      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		    else 
		      cout << "spin " << SpinValue << "/2 are not available" << endl;
		    return -1;
		  }
		}
	    }
	}
      else
	{
	  if (SzFlag == true)
	    {
	      switch (SpinValue)
		{
		default :
		  {
		    if ((SpinValue & 1) == 0)
		      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		    else 
		      cout << "spin " << SpinValue << "/2 are not available" << endl;
		    return -1;
		  }
		}
	    }
	  else
	    {
	      switch (SpinValue)
		{
		case 1 :
		  TmpSpace = new Spin1_2ChainFullInversionAnd2DTranslation (InversionSector, XMomentum, XPeriodicity, YMomentum, YPeriodicity);
		  break;
		default :
		  {
		    if ((SpinValue & 1) == 0)
		      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		    else 
		      cout << "spin " << SpinValue << "/2 are not available" << endl;
		    return -1;
		  }
		}
	    }
	}
      ComplexMatrix TmpComplexEigenstates(Space->GetHilbertSpaceDimension(), NbrStates);
      for (int i = MinIndex; i <= MaxIndex; ++i)
	{
	  TmpComplexEigenstates[i - MinIndex] = TmpSpace->ConvertFromKxKyBasis(ComplexEigenstates[i], Space);
	}
      ComplexEigenstates = TmpComplexEigenstates;
      MinIndex = 0;
      MaxIndex = NbrStates - 1;
      delete TmpSpace;
    }
  
  int MaxSzA = (SubsystemSize * SpinValue);
  int MinSzA = -MaxSzA;
  int MaxSzB = ((NbrSpins - SubsystemSize) * SpinValue);
  int MinSzB = -MaxSzB;
  if (SzFlag == false)
    {
      MaxSzA = 0;
      MinSzA = 0;
      MaxSzB = 0;
      MinSzB =0;
      TotalSz = 0;
    }
  double*** EntanglementSpectra = new double**[((MaxSzA - MinSzA) / 2) + 1];
  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; TmpSzA += 2)
    {
      EntanglementSpectra[(TmpSzA - MinSzA) / 2] = new double*[NbrStates];
    }
  int* EntanglementSpectrumDimension = new int[((MaxSzA - MinSzA) / 2) + 1];

  timeval TotalStartingTime;
  timeval TotalEndingTime;
  if ((Manager.GetBoolean("complex") == false) && (Momentum2DFlag == false))
    {
      for (int TmpSzA = MinSzA; TmpSzA <= MaxSzA; TmpSzA += 2)
	{
	  int SzB = TotalSz - TmpSzA;
	  if ((SzB <= MaxSzB) && (SzB >= MinSzB))
	    {
	      if (SzFlag == true)
		{
		  cout << "processing subsytem size " << SubsystemSize << " SzA=" << TmpSzA << endl;
		}
	      else
		{
		  cout << "processing subsytem size " << SubsystemSize << endl;
		}

	      SpinChainMultipleEntanglementSpectrumOperation TmpOperation(Space, RealEigenstates, MinIndex, MaxIndex, SubsystemSize, TmpSzA, SubsystemSites);
	      TmpOperation.ApplyOperation(Architecture.GetArchitecture());
	      EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] = TmpOperation.GetEntanglementSpectrumDimension();	      
	      EntanglementSpectra[(TmpSzA - MinSzA) / 2] = TmpOperation.GetEntanglementSpectra();
	    }
	  else
	    {
	      EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] = 0;
	    }  
	}
    }
  else
    {
      for (int TmpSzA = MinSzA; TmpSzA <= MaxSzA; TmpSzA += 2)
	{
	  int SzB = TotalSz - TmpSzA;
	  if ((SzB <= MaxSzB) && (SzB >= MinSzB))
	    {
	      if (SzFlag == true)
		{
		  cout << "processing subsytem size " << SubsystemSize << " SzA=" << TmpSzA << endl;
		}
	      else
		{
		  cout << "processing subsytem size " << SubsystemSize << endl;
		}
	      SpinChainMultipleEntanglementSpectrumOperation TmpOperation(Space, ComplexEigenstates, MinIndex, MaxIndex, SubsystemSize, TmpSzA, SubsystemSites);
	      TmpOperation.ApplyOperation(Architecture.GetArchitecture());
	      EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] = TmpOperation.GetEntanglementSpectrumDimension();	      
	      EntanglementSpectra[(TmpSzA - MinSzA) / 2] = TmpOperation.GetEntanglementSpectra();
	    }
	  else
	    {
	      EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] = 0;
	    }  
	}
    }
  
  int TotalDimensionPerState = 0;
  int NbrSzASectors = 0;
  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; TmpSzA += 2)
    {
      TotalDimensionPerState += EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2];
      if (EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] > 0)
	{
	  NbrSzASectors++;
	}
    }
  ofstream File;
  File.open(DensityMatrixFileName, ios::binary | ios::out);
  WriteLittleEndian(File, NbrStates);  
  WriteLittleEndian(File, NbrSzASectors);  
  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; TmpSzA += 2)
    {
      if (EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] > 0)
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
      for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; TmpSzA += 2)
	{
	  if (EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] > 0)
	    {
	      WriteLittleEndian(File, EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2]);
	      WriteBlockLittleEndian(File, EntanglementSpectra[(TmpSzA - MinSzA) / 2][Index - MinIndex], EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2]);
	      for (int i = 0; i < EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2]; ++i)
		{
		  double Tmp = EntanglementSpectra[(TmpSzA - MinSzA) / 2][Index - MinIndex][i];
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
      for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; TmpSzA += 2)
	{
	  File << " counting(Sza=" << TmpSzA << ")";
	}
      File << endl;
      int* TmpCounting = new int [((MaxSzA - MinSzA) / 2) + 1];
      double CountingError = Manager.GetDouble("counting-error");
      for (int Index = MinIndex; Index <= MaxIndex; ++Index)
	{
	  long TotalNbr = 0;
	  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; TmpSzA += 2)
	    {
	      TmpCounting[(TmpSzA - MinSzA) / 2] = 0;
	      for (int i = 0; i < EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2]; ++i)
		{
		  ++TotalNbr;
		  if (EntanglementSpectra[(TmpSzA - MinSzA) / 2][Index - MinIndex][i] > CountingError)
		    ++TmpCounting[(TmpSzA - MinSzA) / 2];
		}
	    }		  
	  int TmpTotalCounting = 0;
	  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; TmpSzA += 2)
	    {
	      TmpTotalCounting += TmpCounting[(TmpSzA - MinSzA) / 2];
	    }
	  File << Index << " " << TmpTotalCounting << " " << TotalNbr;
	  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; TmpSzA += 2)
	    {
	      File << " " << TmpCounting[(TmpSzA - MinSzA) / 2];
	    }	  
	  File << endl;
	}
      File.close();
    }
  return 0;
}


// get a linearized position index from the 2d coordinates
//
// xPosition = position along the x direction
// xPeriodicity = system size in the x direction
// yPosition = position along the y direction
// yPeriodicity = system size in the y direction
// return value = linearized index

inline int SpinChainMultipleEntanglementSpectraGetLinearizedIndex(int xPosition, int xPeriodicity, int yPosition, int yPeriodicity)
{
  return ((xPosition * yPeriodicity) + yPosition);
  
}

