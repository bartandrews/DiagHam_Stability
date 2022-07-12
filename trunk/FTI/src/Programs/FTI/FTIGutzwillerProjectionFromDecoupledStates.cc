#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"
#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"
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

int main(int argc, char** argv)
{
  OptionManager Manager ("FTIGutzwillerProjectionFromDecoupledStates" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('1', "input-file1", "name of the first decoulped state");
  (*SystemGroup) += new SingleStringOption  ('2', "input-file2", "name of the second decoulped state");  
  (*SystemGroup) += new SingleStringOption  ('\n', "file-list", "two column file describing a list to states that have to be projected");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-weightcutoff", "do not discard the projected state has its square norm below machine accuracy");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (while happen .x.vec at the end of each stored vector)");
  (*OutputGroup) += new BooleanOption ('\n', "disable-gutzwillerbasis", "do not express the projected wave function in the Gutzwiller reduced Hilbert space but the full SU(2) basis");
  (*OutputGroup) += new BooleanOption ('\n', "unnormalized", "store the Gutzwiller projected states without normalizing them");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIGutzwillerProjectionFromDecoupledStates -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
    
  int NbrSpaces = 1;
  ComplexVector* InputStates = 0;
  char** InputStateFiles = 0;
  int* NbrParticles = 0;
  int* TotalKx = 0;
  int* TotalKy = 0;
  int NbrSites = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  bool Statistics = true;
  int* TotalSpin = 0;
  bool TwoDTranslationFlag = false;
  bool TotalSpinConservedFlag = false;
  bool SU2SpinFlag = false;
  bool GutzwillerFlag = false;

  if (((Manager.GetString("input-file1") == 0) || (Manager.GetString("input-file2") == 0)) && (Manager.GetString("file-list") == 0))
    {
      cout << "error, an input state file should be provided. See man page for option syntax or type FTIGutzwillerProjectionFromDecoupledStates -h" << endl;
      return -1;
    }
  if ((Manager.GetString("input-file1") != 0) && 
      (IsFile(Manager.GetString("input-file1")) == false))
    {
      cout << "can't open file " << Manager.GetString("input-file1") << endl;
      return -1;
    }
  if ((Manager.GetString("input-file2") != 0) && 
      (IsFile(Manager.GetString("input-file2")) == false))
    {
      cout << "can't open file " << Manager.GetString("input-file2") << endl;
      return -1;
    }
  if ((Manager.GetString("file-list") != 0) && 
      (IsFile(Manager.GetString("file-list")) == false))
    {
      cout << "can't open file " << Manager.GetString("file-list") << endl;
      return -1;
    }

  if (Manager.GetString("file-list") == 0)
    {
      NbrSpaces = 2;
      InputStateFiles = new char* [2];
      NbrParticles = new int[2];
      TotalKx = new int[2];
      TotalKy = new int[2];
      TotalSpin = new int[2];
      InputStateFiles[0] = new char [strlen(Manager.GetString("input-file1")) + 1];
      InputStateFiles[1] = new char [strlen(Manager.GetString("input-file2")) + 1];
      strcpy (InputStateFiles[0], Manager.GetString("input-file1"));
      strcpy (InputStateFiles[1], Manager.GetString("input-file2"));
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("file-list")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
       NbrSpaces = 2 * DegeneratedFile.GetNbrLines();
       InputStateFiles = new char* [NbrSpaces];
       NbrParticles = new int[NbrSpaces];
       TotalKx = new int[NbrSpaces];
       TotalKy = new int[NbrSpaces];
       TotalSpin = new int[NbrSpaces];
       for (int i = 0; i < DegeneratedFile.GetNbrLines(); i++)
	 {
	   InputStateFiles[2 * i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (InputStateFiles[2 * i], DegeneratedFile(0, i));		   
	   InputStateFiles[2 * i + 1] = new char [strlen(DegeneratedFile(1, i)) + 1];
	   strcpy (InputStateFiles[2 * i + 1], DegeneratedFile(1, i));		   
	 }
    }

  int TmpNbrParticles = 0;
  if (FTIHubbardModelFindSystemInfoFromVectorFileName(InputStateFiles[0], TmpNbrParticles, NbrSites, Statistics, GutzwillerFlag) == false)
    {
      cout << "error while retrieving system parameters from file name " << InputStateFiles[0] << endl;
      return -1;      
    }
  TwoDTranslationFlag = FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(InputStateFiles[0],
											 TmpNbrParticles, NbrSites, TotalKx[0], TotalKy[0], NbrSiteX, NbrSiteY, Statistics, GutzwillerFlag);
  TotalSpinConservedFlag = FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(InputStateFiles[0], TmpNbrParticles, NbrSites, TotalSpin[0], Statistics, GutzwillerFlag);
  if (TotalSpinConservedFlag == true)
    {
      SU2SpinFlag = true;
    }

  if (TwoDTranslationFlag == true)
    { 
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  NbrParticles[i] = 0;
	  TotalKx[i] = 0;
	  TotalKy[i] = 0;
	  TotalSpin[i] = 0;
	  if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(InputStateFiles[i],
									       NbrParticles[i], NbrSites, TotalKx[i], TotalKy[i], NbrSiteX, NbrSiteY, Statistics, GutzwillerFlag) == false)
	    {
	      cout << "error while retrieving 2D translation parameters from file name " << InputStateFiles[i] << endl;
	      return -1;
	    }
	  if (TotalSpinConservedFlag == true)
	    {
	      if (FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(InputStateFiles[i], NbrParticles[i], NbrSites, TotalSpin[i], Statistics, GutzwillerFlag) == false)
		{
		  
		  cout << "error while retrieving Sz value from file name " << InputStateFiles[i] << endl;
		  return -1;
		}
	    }
	} 
    }
  else
    {
       for (int i = 0; i < NbrSpaces; ++i)
	{
	  NbrParticles[i] = 0;
	  if (FTIHubbardModelFindSystemInfoFromVectorFileName(InputStateFiles[i], NbrParticles[i], NbrSites, Statistics, GutzwillerFlag) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << InputStateFiles[0] << endl;
	      return -1;      
	    }
	}
   }
  
  int TotalNbrSites = NbrSiteX * NbrSiteY;
  if (TwoDTranslationFlag == false)
    TotalNbrSites = 1;
  int* NbrInputStatePerMomentumSector = new int[TotalNbrSites];
  ComplexVector** InputStatePerMomentumSector = new ComplexVector*[TotalNbrSites];
  double** CoefficientPerMomentumSector = new double*[TotalNbrSites];
  for (int i = 0; i < TotalNbrSites; ++i)
    {
      NbrInputStatePerMomentumSector[i] = 0;
      InputStatePerMomentumSector[i] = 0;
    }

  InputStates = new ComplexVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (InputStates[i].ReadVector (InputStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << InputStateFiles[i] << endl;
	  return -1;      
	}
    }

  int MaxNbrSpaces;
  if (TwoDTranslationFlag == false)
    MaxNbrSpaces = 1;
  else
    MaxNbrSpaces = NbrSiteX * NbrSiteY;
  ParticleOnSphere** Spaces = new ParticleOnSphere*[MaxNbrSpaces];
  ParticleOnSphereWithSpin** TargetSpaces = new ParticleOnSphereWithSpin*[MaxNbrSpaces];
  for (int i = 0; i < MaxNbrSpaces; ++i)
    {
      Spaces[i] = 0;
      TargetSpaces[i] = 0;      
    }
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (InputStates[i].ReadVector (InputStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << InputStateFiles[i] << endl;
	  return -1;      
	}
      int TmpIndex;
      if(TwoDTranslationFlag == true)
	TmpIndex = (TotalKx[i] * NbrSiteY) + TotalKy[i];
      else
	TmpIndex = 0;
      NbrInputStatePerMomentumSector[TmpIndex]++; 
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      int TmpIndex;
      if (TwoDTranslationFlag == false)
	TmpIndex = 0;
      else
	TmpIndex = TotalKx[i] * NbrSiteY + TotalKy[i];
      if (Spaces[TmpIndex] == 0)
	{
	  if (Statistics == true)
	    {
	      if (TwoDTranslationFlag == false)
		{
		  cout << "building Hilbert spinless space N=" << NbrParticles[i] << " Ns=" << NbrSites << endl;
		  Spaces[TmpIndex] = new FermionOnLatticeRealSpace (NbrParticles[i], NbrSites);
		}
	      else
		{
		  cout << "building Hilbert spinless space N=" << NbrParticles[i] << " Ns=" << NbrSites << " kx=" << TotalKx[i] << "/" << NbrSiteX << " ky=" << TotalKy[i] << "/" << NbrSiteY << endl;
		  Spaces[TmpIndex] = new FermionOnLatticeRealSpaceAnd2DTranslation (NbrParticles[i], NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
		}
	    }
	}
      if (Spaces[TmpIndex]->GetLargeHilbertSpaceDimension() != InputStates[i].GetLargeVectorDimension())
	{
	  cout << Spaces[TmpIndex]->GetLargeHilbertSpaceDimension() << " " << InputStates[i].GetLargeVectorDimension() << endl;
	  cout << "dimension mismatch between Hilbert space and ground state" << endl;
	  return 0;
	}
    }

  
  for (int i = 0; i < NbrSpaces; i += 2)
    {
      int TmpIndex;
      if (TwoDTranslationFlag == false)
	TmpIndex = 0;
      else
	TmpIndex = ((TotalKx[i] + TotalKx[i + 1]) % NbrSiteX) * NbrSiteY + ((TotalKy[i] + TotalKy[i + 1]) % NbrSiteY);
      if (TargetSpaces[TmpIndex] == 0)
	{
	  if (Statistics == true)
	    {
	      if (TwoDTranslationFlag == false)
		{
		  if (!Manager.GetBoolean("disable-gutzwillerbasis"))
		    {
		      cout << "building Hilbert spinful Gutzwiller space N=" << (NbrParticles[i] + NbrParticles[i + 1]) << " Ns=" << NbrSites << " 2Sz=" << (NbrParticles[i] - NbrParticles[i + 1]) << endl;
		      TargetSpaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles[i] + NbrParticles[i + 1], NbrParticles[i] - NbrParticles[i + 1], NbrSites);
		    }
		  else
		    {
		      cout << "building Hilbert spinful space N=" << (NbrParticles[i] + NbrParticles[i + 1]) << " Ns=" << NbrSites << " 2Sz=" << (NbrParticles[i] - NbrParticles[i + 1]) << " kx=" << ((TotalKx[i] + TotalKx[i + 1]) % NbrSiteX) << " ky=" << ((TotalKy[i] + TotalKy[i + 1]) % NbrSiteY) << endl;
		      TargetSpaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpace (NbrParticles[i] + NbrParticles[i + 1], NbrParticles[i] - NbrParticles[i + 1], NbrSites);
		    }
		}
	      else
		{
		  
		  cout << "building Hilbert spinless space N=" << NbrParticles[i] << " Ns=" << NbrSites << " kx=" << TotalKx[i] << "/" << NbrSiteX << " ky=" << TotalKy[i] << "/" << NbrSiteY << endl;
		  if (!Manager.GetBoolean("disable-gutzwillerbasis"))
		    {
		      TargetSpaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles[i] + NbrParticles[i + 1],
															     NbrParticles[i] - NbrParticles[i + 1], NbrSites,
															     (TotalKx[i] + TotalKx[i + 1]) % NbrSiteX, NbrSiteX,
															     (TotalKy[i] + TotalKy[i + 1]) % NbrSiteY, NbrSiteY);
		    }
		  else
		    {
		      TargetSpaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles[i] + NbrParticles[i + 1],
												      NbrParticles[i] - NbrParticles[i + 1], NbrSites,
												      (TotalKx[i] + TotalKx[i + 1]) % NbrSiteX, NbrSiteX,
												      (TotalKy[i] + TotalKy[i + 1]) % NbrSiteY, NbrSiteY);
		    }
		}
	    }
	  else
	    {
	      cout << " Bosonic statistics not implemented" << endl;
	    }
	}
    }

  char* TmpGutzwillerOutputName = new char [256];
  if (Manager.GetBoolean("disable-gutzwillerbasis") == true)
    {
      if (Manager.GetBoolean("unnormalized") == false)
	{
	  sprintf(TmpGutzwillerOutputName, "realspace_gutzwillerprojected");
	}
      else
	{
	  sprintf(TmpGutzwillerOutputName, "realspace_unnormalized_gutzwillerprojected");
	}
    }
  else
    {
      if (Manager.GetBoolean("unnormalized") == false)
	{
	  sprintf(TmpGutzwillerOutputName, "realspace_gutzwiller_gutzwillerprojected");
 	}
      else
	{
	  sprintf(TmpGutzwillerOutputName, "realspace_gutzwiller_unnormalized_gutzwillerprojected");
 	}
   }
	      
  for (int i = 0; i < NbrSpaces; i += 2)
    {
      int TmpIndex1;
      int TmpIndex2;
      int TmpTargetIndex;
      if (TwoDTranslationFlag == false)
	{
	  TmpIndex1 = 0;
	  TmpIndex2 = 0;
	  TmpTargetIndex = 0;
	}
      else
	{
	  TmpIndex1 = TotalKx[i] * NbrSiteY + TotalKy[i];
	  TmpIndex2 = TotalKx[i + 1] * NbrSiteY + TotalKy[i + 1];
	  TmpTargetIndex = ((TotalKx[i] + TotalKx[i + 1]) % NbrSiteX) * NbrSiteY + ((TotalKy[i] + TotalKy[i + 1]) % NbrSiteY);
	}
      cout << "projecting " << InputStateFiles[i] << " and " << InputStateFiles[i + 1] << endl;
      ComplexVector TmpVector;
      if (TargetSpaces[TmpTargetIndex] != 0)
	{
	  TmpVector = TargetSpaces[TmpTargetIndex]->ForgeSU2FromU1 (InputStates[i], Spaces[TmpIndex1], InputStates[i + 1], Spaces[TmpIndex2]);
	}
      double TmpWeight = TmpVector.SqrNorm();
      cout << "   weight of the projected state = " << TmpWeight << endl;
      if ((TmpWeight > MACHINE_PRECISION) || (Manager.GetBoolean("disable-weightcutoff")))
	{
	  if (Manager.GetBoolean("unnormalized") == false)
	    {
	      TmpVector /= sqrt(TmpWeight);
	    }
	  char* TmpOutputName;
	  if (Manager.GetString("output-file") != 0)
	    {
	      TmpOutputName = new char[strlen(Manager.GetString("output-file"))+ 24];
	      sprintf(TmpOutputName, "%s.%d.vec", Manager.GetString("output-file"), i);
	    }
	  else
	    {
	      char* TmpString = strstr (InputStateFiles[i], "_realspace_");
	      if (TmpString == 0)
		{
		  cout << "no occurence of _realspace_ was find in " << InputStateFiles[i] << " file name, cannot build output file name" << endl;
		  return 0;
		}
	      TmpOutputName =  new char[strlen(InputStateFiles[i])+ strlen(TmpGutzwillerOutputName) + 256];
	      char Tmp = (*TmpString);
	      TmpString[0] = '\0';
	      sprintf(TmpOutputName, "%s_%s_%s", InputStateFiles[i], TmpGutzwillerOutputName, TmpString + 11);
	      TmpString[0] = Tmp;
	      
	      TmpString = strstr (TmpOutputName, "_n_");
	      Tmp = (*TmpString);
	      if (TmpString == 0)
		{
		  cout << "no occurence of _n_ was find in " << TmpOutputName << " file name, cannot build output file name" << endl;
		  return 0;
		}
	      char* TmpString2 = strstr (TmpOutputName, "_ns_");
	      if (TmpString2 == 0)
		{
		  cout << "no occurence of _ns_ was find in " << TmpOutputName << " file name, cannot build output file name" << endl;
		  return 0;
		}
	      TmpString[0] = '\0';
	      char* TmpOutputName2 =  new char[strlen(TmpOutputName) + strlen(TmpString2) + 24];
	      sprintf(TmpOutputName2, "%s_n_%d%s", TmpOutputName, (NbrParticles[i] + NbrParticles[i + 1]), TmpString2);
	      TmpString[0] = Tmp;
	      delete[]  TmpOutputName;
	      TmpOutputName = TmpOutputName2;

	      int TmpStringSize = strlen(TmpOutputName);
	      int TmpDotPosition = TmpStringSize;
	      for (int i = TmpStringSize; i >= 0; --i)
		{
		  if (TmpOutputName[i] == '.')
		    {		     
		      if (TmpDotPosition == TmpStringSize)
			{
			  TmpDotPosition = i;
			}
		      else
			{
			  TmpDotPosition = i;
			  i = -1;
			}
		    }
		}
	      char* TmpKString = strstr(TmpOutputName, "_kx_");
	      Tmp = TmpKString[0];
	      TmpKString[0] = '\0';
	      TmpOutputName2 =  new char[strlen(TmpOutputName) + 256];
	      sprintf(TmpOutputName2, "%s_kx_%d_ky_%d_sz_%d.%s", TmpOutputName, ((TotalKx[i] + TotalKx[i + 1]) % NbrSiteX),
		      ((TotalKy[i] + TotalKy[i + 1]) % NbrSiteY),
		      (NbrParticles[i] - NbrParticles[i + 1]), TmpOutputName + TmpDotPosition + 1);
	      TmpOutputName[TmpDotPosition] = Tmp;
	      delete[]  TmpOutputName;
	      TmpOutputName = TmpOutputName2;
	    }
	  if (TmpVector.WriteVector(TmpOutputName) == false)
	    {
	      cout << "can't write " << TmpOutputName << endl;
	      return 0;
	    }
	}
    }

  for (int i = 0; i < MaxNbrSpaces; ++i)
    {
      if (Spaces[i] != 0)
	{
	  delete Spaces[i];
	}
    }
  return 0;
}
