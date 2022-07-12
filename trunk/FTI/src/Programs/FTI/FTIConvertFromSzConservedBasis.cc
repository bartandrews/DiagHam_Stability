#include "Vector/ComplexVector.h"

#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnTorusWithSpinAndMagneticTranslations.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "Options/Options.h"

#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"
#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

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
  OptionManager Manager ("FTIConvertFromSzConservedBasis" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption('i', "input-state", "name of the file containing the state whose Kx momentum has to be computed");
  (*SystemGroup) += new SingleStringOption('\n', "degenerate-states", "name of the file containing a list of states (override input-state)");
  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIConvertFromSzConservedBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int NbrSites = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  int XMomentum = 0;
  int YMomentum = 0;
  bool TwoDTranslationFlag;
  int SzValue = 0;
  bool Statistics = true;
  bool GutzwillerFlag = false;
  int NbrInputStates = 0;

  if ((Manager.GetString("input-state") == 0) && (Manager.GetString("degenerate-states") == 0))
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type FTIConvertFromSzConservedBasis -h" << endl;
      return -1;
    }

  if (Manager.GetString("input-state") != 0)
    {
      NbrInputStates = 1;
      if (IsFile(Manager.GetString("input-state")) == false)
	{
	  cout << "can't open file " << Manager.GetString("input-state") << endl;
	  return -1;
	}
      TwoDTranslationFlag = FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSiteX, NbrSiteY, XMomentum, YMomentum, SzValue, Statistics);
      
      if (FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSiteX, NbrSiteY, SzValue, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " <<Manager.GetString("input-state")  << endl;
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
      NbrInputStates = DegenerateFile.GetNbrLines();
      TwoDTranslationFlag = FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), NbrParticles, NbrSiteX, NbrSiteY, XMomentum, YMomentum, SzValue, Statistics);
      if (FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), NbrParticles, NbrSiteX, NbrSiteY, SzValue, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << DegenerateFile(0, 0) << endl;
	  return -1;
	}
    }
  NbrSites = 2*NbrSiteX * NbrSiteY;
  if (TwoDTranslationFlag == false)
  {
    cout << "Convert from Sz basis to full basis" << endl;
    cout << "Nbr particles=" << NbrParticles << " Nbr sites=" << NbrSites << " Sz = " << SzValue << endl;
  }
  else
  {
    cout << "Convert from Sz, kx, ky basis to full basis" << endl;
    cout << "Nbr particles=" << NbrParticles << " Nbr sites=" << NbrSites << " kx = " << XMomentum << " ky = " << YMomentum << " Sz = " << SzValue << endl;
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
	      cout << "error, " << DegenerateFile(0, 0) << " and " <<  DegenerateFile(0, i) << "don't have the same  dimension (" << InputStates[0].GetVectorDimension() << " and " << InputStates[i].GetVectorDimension()<< ")" << endl;
	      return -1;
	    }
	  InputStateNames[i] = new char [strlen(DegenerateFile(0, i)) + 1];
	  strcpy (InputStateNames[i], DegenerateFile(0, i));
	}
    }


  ParticleOnSphereWithSpin* InputSpace = 0;
//   FermionOnLatticeWithSpinRealSpace* InputSpace = 0;
//   FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* TwoDInputSpace = 0;
  
  if (Statistics == true)
    {
      if (TwoDTranslationFlag == false)
      {
	if (GutzwillerFlag == false)
	    InputSpace = new FermionOnLatticeWithSpinRealSpace (NbrParticles, SzValue, NbrSites, 10000000);
	else
	    InputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, SzValue, NbrSites, 10000000);
      }
      else
      {
	if (GutzwillerFlag == false)
	    InputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, SzValue, NbrSites, XMomentum, NbrSiteX, YMomentum, NbrSiteY, 10000000);
	else
	    InputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, SzValue, NbrSites, XMomentum, NbrSiteX, YMomentum, NbrSiteY, 10000000);
      }
    }
  else
    {
      cout << "not available for bosons" << endl;
      return -1;
    }
//   if (TwoDTranslationFlag == false)
  {
    if (InputSpace->GetHilbertSpaceDimension() != InputStates[0].GetVectorDimension())
      {
	cout << "error, " << Manager.GetString("input-state")  << " has a wrong dimension (" << InputStates[0].GetVectorDimension() << ", should be " << InputSpace->GetHilbertSpaceDimension() << ")" << endl;
	return -1;
      }
  }
//   else
//   {
//     if (GutzwillerFlag == false)
// 	  TwoDInputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, SzValue, NbrSites, XMomentum, NbrSiteX, YMomentum, NbrSiteY, 10000000);
//       else
// 	  TwoDInputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, SzValue, NbrSites, XMomentum, NbrSiteX, YMomentum, NbrSiteY, 10000000);
//     if (TwoDInputSpace->GetHilbertSpaceDimension() != InputStates[0].GetVectorDimension())
//       {
// 	cout << "error, " << Manager.GetString("input-state")  << " has a wrong dimension (" << InputStates[0].GetVectorDimension() << ", should be " << TwoDInputSpace->GetHilbertSpaceDimension() << ")" << endl;
// 	return -1;
//       }
//   }
  

//   FermionOnSphereWithSpin* OutputSpace = 0;
  ParticleOnSphereWithSpin* OutputSpace = 0;
  if (Statistics == true)
    {
      if (TwoDTranslationFlag == false)
      {
	if (GutzwillerFlag == false)
	    OutputSpace = new FermionOnLatticeWithSpinRealSpace (NbrParticles, NbrSites);
	else
	    OutputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
      }
      else
      {
	if (GutzwillerFlag == false)
	    OutputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum, NbrSiteX, YMomentum, NbrSiteY);
	else
	    OutputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles,NbrSites, XMomentum, NbrSiteX, YMomentum, NbrSiteY);	
      }
    }
  else
    {
      cout << "not available for bosons" << endl;
      return -1;
    }

  char* SzValueString = new char[256];
  sprintf (SzValueString, "_sz_%d", SzValue);
  
  for (int i = 0; i < NbrInputStates; ++i)
    {
      char* VectorOutputName = ReplaceString(InputStateNames[i], SzValueString, "");
      ComplexVector TmpVector;
      if (TwoDTranslationFlag == false)
      {
      }
      else
      {
	TmpVector = ((FermionOnTorusWithSpinAndMagneticTranslations*) InputSpace)->ConvertFromNbodyBasis (InputStates[i], *((FermionOnTorusWithSpinAndMagneticTranslations*)OutputSpace));
	cout << "dimension = " << TmpVector.GetVectorDimension() <<  " " << ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) OutputSpace)->GetHilbertSpaceDimension() << endl;
// 	TmpVector = ((FermionOnSphereWithSpin*) InputSpace)->ConvertFromNbodyBasis (InputStates[i], *OutputSpace);
      }
//       {
// 	ComplexVector TmpVector1 = TwoDInputSpace->ConvertFromKxKyBasis(InputStates[i], InputSpace);
// 	TmpVector = InputSpace->ConvertFromNbodyBasis (TmpVector1, *OutputSpace);
//       }
      if (TmpVector.WriteVector(VectorOutputName) == false)
	{
	  cout << "error, can't write vector " << VectorOutputName << endl;
	}
      delete[] VectorOutputName;
    }
}
