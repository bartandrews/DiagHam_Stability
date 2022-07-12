#include "Vector/ComplexVector.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin1_2ChainFull.h"
#include "HilbertSpace/Spin1_2ChainFullAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainFullInversionAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainNew.h"
#include "HilbertSpace/Spin1_2ChainNewAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainNewSzSymmetryAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzInversionSymmetries.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzInversionSymmetries.h"
#include "HilbertSpace/Spin2ChainWithTranslations.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndSzInversionSymmetries.h"



#include "Options/Options.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <cstring> 

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("SpinSystemConvertFromTranslationInvariantBasis" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption('i', "input-state", "name of the file containing the state defined in the translation invariant basis");
  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinSystemConvertFromTranslationInvariantBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSites = 0;
  int* XMomentum = 0;
  int* YMomentum = 0;
  int* InversionSector = 0;
  int* SzValue = 0;
  int* SzParitySector = 0;
  bool TotalSpinConservedFlag = true;
  bool InversionFlag = false;
  bool SzSymmetryFlag = false;
  int XPeriodicity = 0;
  int YPeriodicity = 0;
  bool Statistics = true;
  int NbrInputStates = 0;
  int SpinValue = 1;

  if (Manager.GetString("input-state") == 0)
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type SpinSystemConvertFromTranslationInvariantBasis -h" << endl;
      return -1;
    }

  if (Manager.GetString("input-state") != 0)
    {
      NbrInputStates = 1;
      XMomentum = new int[1];
      YMomentum = new int[1];
      InversionSector = new int[1];
      InversionSector[0] = 0;
      SzValue = new int [1];
      SzParitySector = new int[1];
      SzParitySector[0] = 0;
      if (IsFile(Manager.GetString("input-state")) == false)
	{
	  cout << "can't open file " << Manager.GetString("input-state") << endl;
	  return -1;
	}      
      if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SzValue[0], SpinValue, XMomentum[0], XPeriodicity,
								YMomentum[0], YPeriodicity) == false)
	{
	  TotalSpinConservedFlag = true;
	  if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SpinValue, XMomentum[0], XPeriodicity, 
								    YMomentum[0], YPeriodicity) == false)
	    {
	      YPeriodicity = 0;
	      if (SpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SzValue[0], SpinValue, XMomentum[0], InversionSector[0], SzParitySector[0]) == false)
		{
		  if (SpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SzValue[0], SpinValue, XMomentum[0]) == false)
		    {
		      if (SpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SzValue[0], SpinValue) == false)
			{
			  TotalSpinConservedFlag = false;
			  if (SpinFindSystemInfoFromFileName(Manager.GetString("input-state"), NbrSites, SpinValue) == false)
			    {
			      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state") << endl;
			      return -1;
			    }
			}
		    }
		  else
		    {
		      XPeriodicity = NbrSites;
		    }
		}
	      else
		{
		  XPeriodicity = NbrSites;
		  if (InversionSector[0] != 0)
		    InversionFlag = true;
		  if (SzParitySector[0] != 0)
		    SzSymmetryFlag = true;
		}
	    }	  
	  else
	    {
	      InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SpinValue, XMomentum[0], XPeriodicity, 
											     YMomentum[0], YPeriodicity, InversionSector[0]);
	      if (InversionFlag == false)
		{
		  cout << "error while retrieving system parameters from file name " <<Manager.GetString("input-state")  << endl;
		  return -1;
		}
	    }
	}
      else
	{
	  InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SzValue[0], SpinValue, XMomentum[0], XPeriodicity,
											 YMomentum[0], YPeriodicity, InversionSector[0]);
	  if (SzValue[0] == 0)
	    {
	      SpinFindSystemInfoFromVectorFileName((Manager.GetString("input-state")), NbrSites, SzValue[0], SpinValue, InversionSector[0], SzParitySector[0]);
	      cout << NbrSites << " " << SzValue[0] << " " << SpinValue << " " << XMomentum[0] << " " << XPeriodicity << " " << YMomentum[0] << " " <<  YPeriodicity << " " << SzParitySector[0] << endl;
	      if (SzParitySector[0] != 0)
		SzSymmetryFlag = true;
	    }
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      NbrInputStates = DegenerateFile.GetNbrLines();
      XMomentum = new int [NbrInputStates];
      YMomentum = new int [NbrInputStates];
      InversionSector = new int[NbrInputStates];
      SzValue = new int [NbrInputStates];
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  InversionSector[i] = 0;
	}
      
    }
  
  if (YPeriodicity == 0)
    {
      cout << "Convert from kx basis to real space basis" << endl;
      for (int i = 0; i < NbrInputStates; ++i)
	cout << " Nbr sites=" << NbrSites << " Kx=" << XMomentum[i] << endl;
    }
  else
    {
      cout << "Convert from (kx,ky) basis to real space basis" << endl;
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  if ((InversionFlag == false) || (InversionSector[i] == 0))
	    {
	      cout << " Nbr sites=" << NbrSites << "=" << XPeriodicity << "x" << YPeriodicity << " Kx=" << XMomentum[i] << " Ky=" << YMomentum[i] << endl;
	    }
	  else
	    {
	      cout << " Nbr sites=" << NbrSites << "=" << XPeriodicity << "x" << YPeriodicity << " Kx=" << XMomentum[i] << " Ky=" << YMomentum[i] << " I=" << InversionSector[i] << endl;
	    }
	}
    }

  ComplexVector* InputStates = new ComplexVector[NbrInputStates];
  char** InputStateNames = new char*[NbrInputStates];
  if (Manager.GetString("input-state") != 0)
    {
      InputStateNames[0] = new char [strlen(Manager.GetString("input-state")) + 1];
      strcpy (InputStateNames[0], Manager.GetString("input-state"));
      if (InputStates[0].ReadVector (Manager.GetString("input-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("input-state") << endl;
	  return -1;      
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
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
      InputStateNames[0] = new char [strlen(DegenerateFile(0, 0)) + 1];
      strcpy (InputStateNames[0], DegenerateFile(0, 0));
      for (int i = 1; i < NbrInputStates; ++i)
	{
	  if (InputStates[i].ReadVector (DegenerateFile(0, i)) == false)
	    {
	      cout << "can't open vector file " << DegenerateFile(0, i) << endl;
	      return -1;      
	    }	  
	  if (InputStates[0].GetVectorDimension() != InputStates[i].GetVectorDimension())
	    {
	      cout << "error, " << DegenerateFile(0, 0) << " and " <<  DegenerateFile(0, i) << " don't have the same  dimension (" << InputStates[0].GetVectorDimension() << " and " << InputStates[i].GetVectorDimension()<< ")" << endl;
	      return -1;
	    }
	  InputStateNames[i] = new char [strlen(DegenerateFile(0, i)) + 1];
	  strcpy (InputStateNames[i], DegenerateFile(0, i));
	}
    }


  AbstractSpinChain** InputSpace = new AbstractSpinChain* [NbrInputStates];
  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (YPeriodicity == 0)
	{
	  if (TotalSpinConservedFlag == false)
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
		  {
		    if (InversionFlag == true)
		      {
			if (SzSymmetryFlag == true)
			  {
			    InputSpace[i] = new Spin1_2ChainWithTranslationsAndSzInversionSymmetries (NbrSites, XMomentum[i], 1, InversionSector[i], SzParitySector[i], SzValue[i], 1000000, 1000000);
			  }
			else
			  {
			    InputSpace[i] = new Spin1_2ChainWithTranslationsAndInversionSymmetry (NbrSites, XMomentum[i], 1, InversionSector[i], SzValue[i], 1000000, 1000000);
			  }
		      }
		    else
		      {
			if (SzSymmetryFlag == true)
			  {
			    InputSpace[i] = new Spin1_2ChainWithTranslationsAndSzSymmetry (NbrSites, XMomentum[i], 1, SzParitySector[i], SzValue[i], 1000000, 1000000);
			  }
			else
			  {
			    InputSpace[i] = new Spin1_2ChainWithTranslations (NbrSites, XMomentum[i], 1, SzValue[i], 1000000, 1000000);
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
			    InputSpace[i] = new Spin1ChainWithTranslationsAndSzInversionSymmetries (NbrSites, XMomentum[i], InversionSector[i], SzParitySector[i], SzValue[i]);
			  }
			else
			  {
			    InputSpace[i] = new Spin1ChainWithTranslationsAndInversionSymmetry (NbrSites, XMomentum[i], InversionSector[i], SzValue[i]);
			  }
		      }
		    else
		      {
			if (SzSymmetryFlag == true)
			  {
			    InputSpace[i] = new Spin1ChainWithTranslationsAndSzSymmetry (NbrSites, XMomentum[i], SzParitySector[i], SzValue[i]);
			  }
			else
			  {
			    InputSpace[i] = new Spin1ChainWithTranslations (NbrSites, XMomentum[i], SzValue[i]);
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
			    InputSpace[i] = new Spin2ChainWithTranslationsAndSzInversionSymmetries (NbrSites, XMomentum[i], InversionSector[i], SzParitySector[i], SzValue[i]);
			  }
			else
			  {
			    InputSpace[i] = new Spin2ChainWithTranslationsAndInversionSymmetry (NbrSites, XMomentum[i], InversionSector[i], SzValue[i]);
			  }
		      }
		    else
		      {
			if (SzSymmetryFlag == true)
			  {
			    InputSpace[i] = new Spin2ChainWithTranslationsAndSzSymmetry (NbrSites, XMomentum[i], SzParitySector[i], SzValue[i]);
			  }
			else
			  {
			    InputSpace[i] = new Spin2ChainWithTranslations (NbrSites, XMomentum[i], SzValue[i]);
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
	  if (TotalSpinConservedFlag == false)
	    {
	      if ((InversionFlag == false) || (InversionSector[i] == 0))
		{
		  switch (SpinValue)
		    {
		    case 1 :
		      InputSpace[i] = new Spin1_2ChainFullAnd2DTranslation (XMomentum[i], XPeriodicity, YMomentum[i], YPeriodicity);
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
		      InputSpace[i] = new Spin1_2ChainFullInversionAnd2DTranslation (InversionSector[i], XMomentum[i], XPeriodicity, YMomentum[i], YPeriodicity);
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
		 {
		   if (SzSymmetryFlag == false)
		     InputSpace[i] = new Spin1_2ChainNewAnd2DTranslation (NbrSites, SzValue[i], XMomentum[i], XPeriodicity, YMomentum[i], YPeriodicity);
		   else
		   {
		     cout << "Create HilbertSpace" << endl;
		     InputSpace[i] = new Spin1_2ChainNewSzSymmetryAnd2DTranslation (NbrSites, SzValue[i], (1 - SzParitySector[0])/2, XMomentum[i], XPeriodicity, YMomentum[i], YPeriodicity);
		   }
		  break;
		 }
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
  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (InputSpace[i]->GetHilbertSpaceDimension() != InputStates[i].GetVectorDimension())
	{
	  cout << "error, dimension mismatch for vector " << i << " (" << InputStates[i].GetVectorDimension() << ", should be " << InputSpace[i]->GetHilbertSpaceDimension() << ")" << endl;
	  return -1;
	}
    }
  
  AbstractSpinChain** OutputSpace = new AbstractSpinChain* [NbrInputStates];
  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (TotalSpinConservedFlag == false)
	{
	  switch (SpinValue)
	    {
	    case 1 :
	      OutputSpace[i] = new Spin1_2ChainFull (NbrSites);
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
	      OutputSpace[i] = new Spin1_2ChainNew (NbrSites, SzValue[i], 1000000);
	      break;
	    case 2 :
	      OutputSpace[i] = new Spin1Chain (NbrSites, SzValue[i], 1000000);
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

  char* XPeriodicityString = new char[256];
  char* XMomentumString = new char[256];
  
  
  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (YPeriodicity == 0)
	{
	  sprintf (XPeriodicityString, "_x_%d", XPeriodicity);
	  if (strstr(InputStateNames[i], "_kx_") == 0)
	    {
	      if (SzSymmetryFlag == false)
		{
		  sprintf (XMomentumString, "_k_%d", XMomentum[i]);
		}
	      else
		{
		  sprintf (XMomentumString, "_szsym_%d_k_%d", SzParitySector[i], XMomentum[i]);
		}
	    }
	  else
	    {
	      if (SzSymmetryFlag == false)
		{
		  sprintf (XMomentumString, "_kx_%d", XMomentum[i]);
		}
	      else
		{
		  sprintf (XMomentumString, "_szsym_%d_kx_%d", SzParitySector[i], XMomentum[i]);
		}
	    }
	}
      else
	{
	  sprintf (XPeriodicityString, "_x_%d_y_%d", XPeriodicity, YPeriodicity);
	  if (SzSymmetryFlag == false)
	    sprintf (XMomentumString, "_kx_%d_ky_%d", XMomentum[i], YMomentum[i]);
	  else
	    sprintf (XMomentumString, "_szsym_%d_kx_%d_ky_%d", SzParitySector[i], XMomentum[i], YMomentum[i]);
	}
      
      char* VectorOutputName = ReplaceString(InputStateNames[i], XMomentumString, "");
      if (VectorOutputName == 0)
	{
	  cout << "can't replace pattern " << XMomentumString << " in file name " << InputStateNames[i] << endl;
	  return -1;
	}
      ComplexVector TmpVector;
      TmpVector = InputSpace[i]->ConvertFromKxKyBasis(InputStates[i], OutputSpace[i]);
      char* Extension = new char[256];
      sprintf (Extension, "%d.vec", i);
      VectorOutputName = ReplaceExtensionToFileName(VectorOutputName, "vec", Extension);
      if (TmpVector.WriteVector(VectorOutputName) == false)
	{
	  cout << "error, can't write vector " << VectorOutputName << endl;
	}
      delete[] Extension;
      delete[] VectorOutputName;
    }
}
