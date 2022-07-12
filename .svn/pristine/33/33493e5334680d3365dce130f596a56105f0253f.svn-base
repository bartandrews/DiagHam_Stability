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
  OptionManager Manager ("FTIGutzwillerProjection" , "0.01");
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

  (*SystemGroup) += new SingleStringOption  ('i', "input-file", "name of the file on which the Gutzwiller projection has to be applied");
  (*SystemGroup) += new SingleStringOption  ('\n', "file-list", "single column file describing a list to states that have to be projected");
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "particles have a SU(2) spin");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (while happen .x.vec at the end of each stored vector)");
  (*OutputGroup) += new BooleanOption ('\n', "gutzwiller-basis", "express the projected wave function in th Gutzwiller reduced Hilbert space");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIGutzwillerProjection -h" << endl;
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
  int* TotalKx = 0;
  int* TotalKy = 0;
  int NbrParticles = 0;
  int NbrSites = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  bool Statistics = true;
  int* TotalSpin = 0;
  bool TwoDTranslationFlag = false;
  bool TotalSpinConservedFlag = false;
  bool SU2SpinFlag = Manager.GetBoolean("su2-spin");
  bool GutzwillerFlag = false;

  if ((Manager.GetString("input-file") == 0) && (Manager.GetString("file-list") == 0))
    {
      cout << "error, an input state file should be provided. See man page for option syntax or type FTIGutzwillerProjection -h" << endl;
      return -1;
    }
  if ((Manager.GetString("input-file") != 0) && 
      (IsFile(Manager.GetString("input-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("input-file") << endl;
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
      InputStateFiles = new char* [1];
      TotalKx = new int[1];
      TotalKy = new int[1];
      TotalSpin = new int[1];
      InputStateFiles[0] = new char [strlen(Manager.GetString("input-file")) + 1];
      strcpy (InputStateFiles[0], Manager.GetString("input-file"));
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("file-list")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
       NbrSpaces = DegeneratedFile.GetNbrLines();
       InputStateFiles = new char* [NbrSpaces];
       TotalKx = new int[NbrSpaces];
       TotalKy = new int[NbrSpaces];
       TotalSpin = new int[NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   InputStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (InputStateFiles[i], DegeneratedFile(0, i));		   
	 }
    }

  if (FTIHubbardModelFindSystemInfoFromVectorFileName(InputStateFiles[0], NbrParticles, NbrSites, Statistics, GutzwillerFlag) == false)
    {
      cout << "error while retrieving system parameters from file name " << InputStateFiles[0] << endl;
      return -1;
      
    }
  TwoDTranslationFlag = FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(InputStateFiles[0],
											 NbrParticles, NbrSites, TotalKx[0], TotalKy[0], NbrSiteX, NbrSiteY, Statistics, GutzwillerFlag);
  TotalSpinConservedFlag = FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(InputStateFiles[0], NbrParticles, NbrSites, TotalSpin[0], Statistics, GutzwillerFlag);
  if (TotalSpinConservedFlag == true)
    {
      SU2SpinFlag = true;
    }

  if (TwoDTranslationFlag == true)
    { 
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  TotalKx[i] = 0;
	  TotalKy[i] = 0;
	  TotalSpin[i] = 0;
	  if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(InputStateFiles[i],
									       NbrParticles, NbrSites, TotalKx[i], TotalKy[i], NbrSiteX, NbrSiteY, Statistics, GutzwillerFlag) == false)
	    {
	      cout << "error while retrieving 2D translation parameters from file name " << InputStateFiles[i] << endl;
	      return -1;
	    }
	  if (TotalSpinConservedFlag == true)
	    {
	      if (FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(InputStateFiles[i], NbrParticles, NbrSites, TotalSpin[i], Statistics, GutzwillerFlag) == false)
		{
		  
		  cout << "error while retrieving Sz value from file name " << InputStateFiles[i] << endl;
		  return -1;
		}
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
  ParticleOnSphere** TargetSpaces = new ParticleOnSphere*[MaxNbrSpaces];
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
	      if (SU2SpinFlag == false)
		{
		  if (TwoDTranslationFlag == false)
		    {
		      Spaces[TmpIndex] = new FermionOnLatticeRealSpace (NbrParticles, NbrSites);
		    }
		  else
		    {
		      Spaces[TmpIndex] = new FermionOnLatticeRealSpaceAnd2DTranslation (NbrParticles, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
		    }
		}
	      else
		{
		  if (GutzwillerFlag == false)
		    {
		      if (TotalSpinConservedFlag == false)
			{
			  if (TwoDTranslationFlag == false)
			    { 
			      Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpace (NbrParticles, NbrSites);
			      if (Manager.GetBoolean("gutzwiller-basis"))
				{
				  TargetSpaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
				}
			    }
			  else
			    {
			      Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
			      if (Manager.GetBoolean("gutzwiller-basis"))
				{
				  TargetSpaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
				}
			    }
			}
		      else
			{
			  if (TwoDTranslationFlag == false)
			  {
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpace (NbrParticles, TotalSpin[i], NbrSites, 10000000);
			      if (Manager.GetBoolean("gutzwiller-basis"))
				{
				  TargetSpaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, TotalSpin[i], NbrSites, 10000000);
				}
			  }
			  else
			  {
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, TotalSpin[i], NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY, 10000000);
			    if (Manager.GetBoolean("gutzwiller-basis"))
				{
				  TargetSpaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, TotalSpin[i], NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
				}
			  }
			}
		    }
		  else
		    {
		      if (TotalSpinConservedFlag == false)
			{
			  if (TwoDTranslationFlag == false)
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
			  else
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
			}
		      else
			{
			  if (TwoDTranslationFlag == false)
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, TotalSpin[i], NbrSites, 10000000);
			  else
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, TotalSpin[i], NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY, 10000000);
			}
		      
		    }
		}
	    }
	  else
	    {
	      cout << " Bosonic statistics not implemented" << endl;
	    }
	}
      
      if (Spaces[TmpIndex]->GetLargeHilbertSpaceDimension() != InputStates[i].GetLargeVectorDimension())
	{
	  cout << Spaces[TmpIndex]->GetLargeHilbertSpaceDimension() << " " << InputStates[i].GetLargeVectorDimension() << endl;
	  cout << "dimension mismatch between Hilbert space and ground state" << endl;
	  return 0;
	}
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      int TmpIndex;
      if (TwoDTranslationFlag == false)
	TmpIndex = 0;
      else
	TmpIndex = TotalKx[i] * NbrSiteY + TotalKy[i];
      cout << "projecting " << InputStateFiles[i] << endl;
      ComplexVector TmpVector;
      if (TargetSpaces[TmpIndex] != 0)
	TmpVector = TargetSpaces[TmpIndex]->GutzwillerProjection(InputStates[i], Spaces[TmpIndex]);
      else
	TmpVector = Spaces[TmpIndex]->GutzwillerProjection(InputStates[i], Spaces[TmpIndex]);
      double TmpWeight = TmpVector.SqrNorm();
      cout << "   weight of the projected state = " << TmpWeight << endl;
      if (TmpWeight > MACHINE_PRECISION)
	{
	  TmpVector /= sqrt(TmpWeight);
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
	      TmpOutputName =  new char[strlen(InputStateFiles[i])+ 24];
	      char Tmp = (*TmpString);
	      TmpString[0] = '\0';
	      if (Manager.GetBoolean("gutzwiller-basis") == false)
		sprintf(TmpOutputName, "%s_gutzwillerprojected_%s", InputStateFiles[i], TmpString + 11);
	      else
		sprintf(TmpOutputName, "%s_gutzwiller_gutzwillerprojected_%s", InputStateFiles[i], TmpString + 11);
	      TmpString[0] = Tmp;
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
