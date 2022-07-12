#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Operator/SpinS2Operator.h"
#include "Operator/SpinWith1DTranslationS2Operator.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainNew.h"
#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin2Chain.h"
#include "HilbertSpace/Spin1_2ChainFull.h"
#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzInversionSymmetries.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzInversionSymmetries.h"
#include "HilbertSpace/Spin1_2ChainFullAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainFullInversionAnd2DTranslation.h"
#include "HilbertSpace/Spin2ChainWithTranslations.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndSzInversionSymmetries.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "Options/Options.h"

#include "GeneralTools/Endian.h"
#include "GeneralTools/StringTools.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;




int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("SpinChainComputeS2" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "state whose total S^2 has to be computed");  
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerate-states", "a single column formatted ASCCI file that contains a list of states whose total S^2 have to be computed");  
  (*SystemGroup) += new BooleanOption  ('c', "complex", "consider complex wave function");
  (*SystemGroup) += new BooleanOption  ('\n', "force-real", "assume real eigenstates even for momentum eigenstates");
  (*OutputGroup) += new BooleanOption  ('\n', "compute-eigenstates", "compute the S^2 eigenstates when using the --degenerate-states option");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinChainComputeS2 -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if ((Manager.GetString("input-state") == 0) && (Manager.GetString("degenerate-states") == 0))
    {
      cout << "error, an eigenstate file should be provided. See man page for option syntax or type SpinChainComputeS2 -h" << endl;
      return -1;
    }

  int SpinValue = 0;
  int NbrSpins = 0;
  int TotalSz = 0;
  bool SzFlag = true;
  bool Momentum2DFlag = false;
  bool Momentum1DFlag = false;
  bool SzSymmetryFlag = false;  
  bool InversionFlag = false;  
  int XMomentum = 0;
  int XPeriodicity = 0;
  int YMomentum = 0;
  int YPeriodicity = 0;
  int InversionSector = 0;
  int SzSymmetrySector = 0;
  
  int NbrStates = 1;
  char** InputStateNames = 0;
  if (Manager.GetString("degenerate-states") == 0)
    {      
      InputStateNames = new char*[NbrStates];
      InputStateNames[0] = new char[strlen(Manager.GetString("input-state")) + 1];
      strcpy (InputStateNames[0], Manager.GetString("input-state"));
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      NbrStates = DegenerateFile.GetNbrLines();
      InputStateNames = new char*[NbrStates];
      for (int i = 0; i < NbrStates; ++i)
	{
	  InputStateNames[i] = new char[strlen(DegenerateFile(0, i)) + 1];
	  strcpy (InputStateNames[i], DegenerateFile(0, i));
	}
    }
  if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(InputStateNames[0], NbrSpins, TotalSz, SpinValue, XMomentum, XPeriodicity,
							    YMomentum, YPeriodicity) == false)
    {
      if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(InputStateNames[0], NbrSpins, SpinValue, XMomentum, XPeriodicity, 
								YMomentum, YPeriodicity) == false)
	{
	  if (SpinFindSystemInfoFromVectorFileName(InputStateNames[0], NbrSpins, TotalSz, SpinValue, XMomentum, InversionSector, SzSymmetrySector) == false)
	    {
	      if (SpinFindSystemInfoFromVectorFileName(InputStateNames[0], NbrSpins, TotalSz, SpinValue, XMomentum) == false)
		{
		  if (SpinFindSystemInfoFromVectorFileName(InputStateNames[0], NbrSpins, TotalSz, SpinValue) == false)
		    {
		      SzFlag = false;
		      if (SpinFindSystemInfoFromFileName(InputStateNames[0], NbrSpins, SpinValue) == false)
			{
			  cout << "error while retrieving system parameters from file name " << InputStateNames[0] << endl;
			  return -1;
			}
		    }
		}
	      else
		{
		  XPeriodicity = NbrSpins;
		  Momentum1DFlag = true;	      
		}
	    }
	  else
	    {
	      XPeriodicity = NbrSpins;
	      Momentum1DFlag = true;	      
	      if (InversionSector != 0)
		InversionFlag = true;
	      if (SzSymmetrySector != 0)
		SzSymmetryFlag = true;
	    }
	}
      else
	{
	  Momentum2DFlag = true;
	  SzFlag = false;
	  InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(InputStateNames[0], NbrSpins, SpinValue, XMomentum, XPeriodicity, 
											 YMomentum, YPeriodicity, InversionSector);
	}
    }
  else
    {
      InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(InputStateNames[0], NbrSpins, TotalSz, SpinValue, XMomentum, XPeriodicity,
										     YMomentum, YPeriodicity, InversionSector);
      SzFlag = true;
      Momentum2DFlag = true;
    }
  if ((Momentum2DFlag == false) && (Momentum1DFlag == false))
    {
      if (SzFlag == true)
	cout << "N=" << NbrSpins << " Sz=" <<  TotalSz << " 2s=" << SpinValue << endl;
      else
	cout << "N=" << NbrSpins << " 2s=" << SpinValue << endl;
    }
  else
    {
      if (Momentum1DFlag == false)
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
      else
	{
	  if (SzFlag == true)
	    cout << "N=" << NbrSpins << " Sz=" <<  TotalSz << " 2s=" << SpinValue << " kx=" << XMomentum << endl;
	  else
	    cout << "N=" << NbrSpins << " 2s=" << SpinValue << " kx=" << XMomentum  << endl;
	}
   }

  if ((Momentum2DFlag == false) && (Momentum1DFlag == false))
    {
      AbstractSpinChain* Space;
      
      if (SzFlag == true)
	{
	  switch (SpinValue)
	    {
	    case 1 :
	      Space = new Spin1_2ChainNew (NbrSpins, TotalSz, 1000000);
	      break;
	    case 2 :
	      Space = new Spin1Chain (NbrSpins, TotalSz, 1000000);
	      break;
	    case 4 :
	      Space = new Spin2Chain (NbrSpins, TotalSz, 1000000);
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
      SpinS2Operator TmpOperator(Space, NbrSpins);
      if (Manager.GetBoolean("complex") == false)
	{
	  RealVector TmpState;
	  if (TmpState.ReadVector (InputStateNames[0]) == false)
	    {
	      cout << "can't open vector file " << InputStateNames[0] << endl;
	      return -1;      
	    }
	  
	  Complex TmpS2 = TmpOperator.MatrixElement(TmpState, TmpState);
	  cout << "S^2 = " << TmpS2.Re << endl; 
	  cout << "2S = " << (sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0) << " " << round(sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0) << endl; 
	}
      else
	{
	  ComplexVector TmpState;
	  if (TmpState.ReadVector (InputStateNames[0]) == false)
	    {
	      cout << "can't open vector file " << InputStateNames[0] << endl;
	      return -1;      
	    }
	  
	  Complex TmpS2 = TmpOperator.MatrixElement(TmpState, TmpState);
	  cout << "<S^2>=" << TmpS2.Re << " <S>=" << (0.5 * (sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0)) << endl;
	  cout << "round(<2S>)=" << round(sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0) << endl; 
	}
    }
  else
    {      
      if (Momentum1DFlag == false)
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
	}
      else
	{
	  AbstractSpinChainWithTranslations* Space = 0;
	  if (SzFlag == true)
	    {
	      switch (SpinValue)
		{
		case 1 :
		  {
		    if (InversionFlag == true)
		      {
			if (SzSymmetryFlag == true)
			  {
			    Space = new Spin1_2ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, XMomentum, 1, InversionSector, SzSymmetrySector, TotalSz, 1000000, 1000000);
			  }
			else
			  {
			    Space = new Spin1_2ChainWithTranslationsAndInversionSymmetry (NbrSpins, XMomentum, 1, InversionSector, TotalSz, 1000000, 1000000);
			  }
		      }
		    else
		      {
			if (SzSymmetryFlag == true)
			  {
			    Space = new Spin1_2ChainWithTranslationsAndSzSymmetry (NbrSpins, XMomentum, 1, SzSymmetrySector, TotalSz, 1000000, 1000000);
			  }
			else
			  {
			    Space = new Spin1_2ChainWithTranslations (NbrSpins, XMomentum, 1, TotalSz, 1000000, 1000000);
			  }
		      }
		  }
		  break;
		case 2 :
		  {
		    if (InversionFlag == true)
		      {
			if (SzSymmetryFlag == true)
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
			if (SzSymmetryFlag == true)
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
		    if (InversionFlag == true)
		      {
			if (SzSymmetryFlag == true)
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
			if (SzSymmetryFlag == true)
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
	  else
	    {
	    }

	  if (Manager.GetBoolean("force-real") == true)
	    {
	      SpinWith1DTranslationS2Operator TmpOperator(Space, NbrSpins);
	      
	      if (Manager.GetString("degenerate-states") == 0)
		{
		  RealVector TmpState;
		  if (TmpState.ReadVector (InputStateNames[0]) == false)
		    {
		      cout << "can't open vector file " << InputStateNames[0] << endl;
		      return -1;      
		    }
		  
		  Complex TmpS2 = TmpOperator.MatrixElement(TmpState, TmpState);
		  cout << "<S^2>=" << TmpS2.Re << " <S>=" << (0.5 * (sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0)) << endl;
		  cout << "round(<2S>)=" << round(sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0) << endl; 
		}
	      else
		{
		  MultiColumnASCIIFile DegenerateFile;
		  RealVector* InputStates = new RealVector[NbrStates];
		  if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
		    {
		      DegenerateFile.DumpErrors(cout);
		      return -1;
		    }
		  if (InputStates[0].ReadVector (DegenerateFile(0, 0)) == false)
		    {
		      cout << "can't open vector file " << DegenerateFile(0, 0) << endl;
		      return -1;      
		    }	  
		  for (int i = 1; i < NbrStates; ++i)
		    {
		      if (InputStates[i].ReadVector (DegenerateFile(0, i)) == false)
			{
			  cout << "can't open vector file " << DegenerateFile(0, i) << endl;
			  return -1;      
			}	  
		      if (InputStates[0].GetVectorDimension() != InputStates[i].GetVectorDimension())
			{
			  cout << "error, " << DegenerateFile(0, 0) << " and " <<  DegenerateFile(0, i) << "don't have the same  dimension (" << InputStates[0].GetVectorDimension() << " and " << InputStates[i].GetVectorDimension()<< ")" << endl;
			  return -1;
			}
		    }
		  RealSymmetricMatrix S2Matrix(NbrStates, true);
		  for (int i = 0; i < NbrStates; ++i)
		    {
		      for (int j = i; j < NbrStates; ++j)
			{
			  Complex TmpS2 = TmpOperator.MatrixElement(InputStates[i], InputStates[j]);
			  S2Matrix.SetMatrixElement(j, i, TmpS2.Re);
			}
		    }
		  RealDiagonalMatrix TmpS2Eigenvalues(NbrStates);
		  RealMatrix TmpBasis (NbrStates, NbrStates);
		  TmpBasis.SetToIdentity();
		  S2Matrix.LapackDiagonalize(TmpS2Eigenvalues, TmpBasis);
		  for (int i = 0; i < NbrStates; ++i)
		    {
		      double TmpS2 = TmpS2Eigenvalues[i];
		      cout << "<S^2>=" << TmpS2 << " <S>=" << (0.5 * (sqrt((4.0 * TmpS2) + 1.0) - 1.0)) << endl;
		      cout << "round(<2S>)=" <<  round(sqrt((4.0 * TmpS2) + 1.0) - 1.0) << endl; 
		    }
		  if (Manager.GetBoolean("compute-eigenstates"))
		    {
		      RealMatrix TmpMatrix (InputStates, NbrStates);
		      TmpMatrix.Multiply(TmpBasis);
		      for (int i = 0; i < NbrStates; ++i)
			{
			  char* SzString = new char [16];
			  sprintf (SzString, "sz_%d", TotalSz);
			  char* SzS2String = new char [32];
			  sprintf (SzS2String, "sz_%d_s_%d", TotalSz, (int) (round(sqrt((4.0 * TmpS2Eigenvalues[i]) + 1.0) - 1.0)));
			  char* TmpOutputName = ReplaceString(InputStateNames[i], SzString, SzS2String);
			  cout << "writing " << TmpOutputName << endl;		      
			  if (TmpMatrix[i].WriteVector(TmpOutputName) == false)
			    {
			      cout << "error, can't write " << TmpOutputName << endl;
			    }
			  delete[] TmpOutputName;
			}
		    }
		  
		}
	    }
	  else
	    {
	      SpinWith1DTranslationS2Operator TmpOperator(Space, NbrSpins);
	      
	      if (Manager.GetString("degenerate-states") == 0)
		{
		  ComplexVector TmpState;
		  if (TmpState.ReadVector (InputStateNames[0]) == false)
		    {
		      cout << "can't open vector file " << InputStateNames[0] << endl;
		      return -1;      
		    }
		  
		  Complex TmpS2 = TmpOperator.MatrixElement(TmpState, TmpState);
		  cout << "<S^2>=" << TmpS2.Re << " <S>=" << (0.5 * (sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0)) << endl;
		  cout << "round(<2S>)=" << round(sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0) << endl; 
		}
	      else
		{
		  MultiColumnASCIIFile DegenerateFile;
		  ComplexVector* InputStates = new ComplexVector[NbrStates];
		  if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
		    {
		      DegenerateFile.DumpErrors(cout);
		      return -1;
		    }
		  if (InputStates[0].ReadVector (DegenerateFile(0, 0)) == false)
		    {
		      cout << "can't open vector file " << DegenerateFile(0, 0) << endl;
		      return -1;      
		    }	  
		  for (int i = 1; i < NbrStates; ++i)
		    {
		      if (InputStates[i].ReadVector (DegenerateFile(0, i)) == false)
			{
			  cout << "can't open vector file " << DegenerateFile(0, i) << endl;
			  return -1;      
			}	  
		      if (InputStates[0].GetVectorDimension() != InputStates[i].GetVectorDimension())
			{
			  cout << "error, " << DegenerateFile(0, 0) << " and " <<  DegenerateFile(0, i) << "don't have the same  dimension (" << InputStates[0].GetVectorDimension() << " and " << InputStates[i].GetVectorDimension()<< ")" << endl;
			  return -1;
			}
		    }
		  HermitianMatrix S2Matrix(NbrStates, true);
		  for (int i = 0; i < NbrStates; ++i)
		    {
		      for (int j = i; j < NbrStates; ++j)
			{
		      Complex TmpS2 = TmpOperator.MatrixElement(InputStates[i], InputStates[j]);
		      S2Matrix.SetMatrixElement(j, i, TmpS2);
			}
		    }
		  RealDiagonalMatrix TmpS2Eigenvalues(NbrStates);
		  ComplexMatrix TmpBasis (NbrStates, NbrStates);
		  TmpBasis.SetToIdentity();
		  S2Matrix.LapackDiagonalize(TmpS2Eigenvalues, TmpBasis);
		  for (int i = 0; i < NbrStates; ++i)
		    {
		      double TmpS2 = TmpS2Eigenvalues[i];
		      cout << "<S^2>=" << TmpS2 << " <S>=" << (0.5 * (sqrt((4.0 * TmpS2) + 1.0) - 1.0)) << endl;
		      cout << "round(<2S>)=" <<  round(sqrt((4.0 * TmpS2) + 1.0) - 1.0) << endl; 
		    }
		  if (Manager.GetBoolean("compute-eigenstates"))
		    {
		      ComplexMatrix TmpMatrix (InputStates, NbrStates);
		      TmpMatrix.Multiply(TmpBasis);
		      for (int i = 0; i < NbrStates; ++i)
			{
			  char* SzString = new char [16];
			  sprintf (SzString, "sz_%d", TotalSz);
			  char* SzS2String = new char [32];
			  sprintf (SzS2String, "sz_%d_s_%d", TotalSz, (int) (round(sqrt((4.0 * TmpS2Eigenvalues[i]) + 1.0) - 1.0)));
			  char* TmpOutputName = ReplaceString(InputStateNames[i], SzString, SzS2String);
			  cout << "writing " << TmpOutputName << endl;		      
			  if (TmpMatrix[i].WriteVector(TmpOutputName) == false)
			    {
			      cout << "error, can't write " << TmpOutputName << endl;
			    }
			  delete[] TmpOutputName;
			}
		    }
		  
		}
	    }
	}
    }
  
  return 0;
}

