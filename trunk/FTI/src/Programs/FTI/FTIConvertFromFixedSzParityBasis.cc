#include "Vector/ComplexVector.h"

#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnTorusWithSpinAndMagneticTranslations.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include "Operator/ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator.h"

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
  OptionManager Manager ("FTIConvertFromFixedSzParityBasis" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption ('\n', "force-symmetrizesector", "when symmetrizing the input state (instead of unsymmetrizing it), force the symmetry sector (can be either +1 or -1, 0 of autodetection should be used)", 0);
  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIConvertFromFixedSzParityBasis -h" << endl;
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
  bool TwoDTranslationFlag = false;
  int SzValue = 0;
  int SzSymmetrySector = 0;
  bool Statistics = true;
  bool GutzwillerFlag = false;
  int NbrInputStates = 0;

  if ((Manager.GetString("input-state") == 0) && (Manager.GetString("degenerate-states") == 0))
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type FTIConvertFromFixedSzParityBasis -h" << endl;
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
      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites, SzValue, SzSymmetrySector, XMomentum, YMomentum, NbrSiteX, NbrSiteY, Statistics, GutzwillerFlag) == false)
	{
	  if (FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites, SzValue, SzSymmetrySector, Statistics, GutzwillerFlag) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state") << endl;
	      return -1;
	    }
	}
      else
	{
	  TwoDTranslationFlag = true;
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
      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), NbrParticles, NbrSites, SzValue, SzSymmetrySector, XMomentum, YMomentum, NbrSiteX, NbrSiteY, Statistics, GutzwillerFlag) == false)
	{
	  if (FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), NbrParticles, NbrSites, SzValue, SzSymmetrySector, Statistics, GutzwillerFlag) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << DegenerateFile(0, 0) << endl;
	      return -1;
	    }
	}
      else
	{
	  TwoDTranslationFlag = true;
	}
    }
  if (TwoDTranslationFlag == false)
    {
      cout << "Convert from fixed Sz parity basis to full basis" << endl;
      cout << "Nbr particles=" << NbrParticles << " Nbr sites=" << NbrSites << " Sz = " << SzValue << endl;
    }
  else
    {
      cout << "Convert from fixed Sz parity, kx, ky basis to full basis" << endl;
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


  ParticleOnSphereWithSpin* OutputSpace = 0;
  if (Statistics == true)
    {
      if (TwoDTranslationFlag == false)
	{
	  if (GutzwillerFlag == false)
	    {
	      OutputSpace = new FermionOnLatticeWithSpinRealSpace (NbrParticles, SzValue, NbrSites);
	    }
	  else
	    {
	      OutputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, SzValue, NbrSites);
	    }
	}
      else
	{
	  if (GutzwillerFlag == false)
	    {
	      OutputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, SzValue, NbrSites, XMomentum, NbrSiteX, YMomentum, NbrSiteY);
	    }
	  else
	    {
	      OutputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, SzValue, NbrSites, XMomentum, NbrSiteX, YMomentum, NbrSiteY);	
	    }
	}
    }
  else
    {
      cout << "not available for bosons" << endl;
      return -1;
    }

  int InputSzSymmetrySector = SzSymmetrySector;
  if (InputSzSymmetrySector == 0)
    {
      if (Manager.GetInteger("force-symmetrizesector") == 0)
	{
	  if (OutputSpace->GetHilbertSpaceDimension() != InputStates[0].GetVectorDimension())
	    {
	      cout << "error, " << Manager.GetString("input-state")  << " has a wrong dimension (" << InputStates[0].GetVectorDimension() << ", should be " << OutputSpace->GetHilbertSpaceDimension() << ")" << endl;
	      return -1;
	    }
	  Architecture.GetArchitecture()->SetDimension(OutputSpace->GetHilbertSpaceDimension());

 	  ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator TmpOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) OutputSpace);
	  OperatorMatrixElementOperation Operation(&TmpOperator, InputStates[0], InputStates[0]);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
//	  Complex TmpElement = TmpOperator.MatrixElement(InputStates[0], InputStates[0]);
	  Complex TmpElement = Operation.GetScalar();
	  cout << "Sz parity of " <<  Manager.GetString("input-state") << " = " << TmpElement << endl;
	  if (TmpElement.Re > 0.0)
	    {
	      InputSzSymmetrySector = 1;
	    }
	  else
	    {
	      InputSzSymmetrySector = -1;
	    }
	}
      else
	{
	  InputSzSymmetrySector = Manager.GetInteger("force-symmetrizesector");
	}
    }
  else
    {
    }
 
  ParticleOnSphereWithSpin* InputSpace = 0;
//   FermionOnLatticeWithSpinRealSpace* InputSpace = 0;
//   FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* TwoDInputSpace = 0;
  if (Statistics == true)
    {
      if (TwoDTranslationFlag == false)
	{
	  if (GutzwillerFlag == false)
	    {
	      InputSpace = new FermionOnLatticeWithSpinSzSymmetryRealSpace (NbrParticles, SzValue, NbrSites, (InputSzSymmetrySector == -1), 10000000);
	    }
	  else
	    {
	      InputSpace = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace (NbrParticles, SzValue, NbrSites, (InputSzSymmetrySector == -1), 10000000);
	    }
	}
      else
	{
	  if (GutzwillerFlag == false)
	    {
	      InputSpace = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (NbrParticles, SzValue, NbrSites, (InputSzSymmetrySector == -1),
											    XMomentum, NbrSiteX, YMomentum, NbrSiteY, 10000000);
	    }
	  else
	    {
	      InputSpace = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, SzValue, NbrSites, (InputSzSymmetrySector == -1),
														   XMomentum, NbrSiteX, YMomentum, NbrSiteY, 10000000);
	    }
	}
    }
  else
    {
      cout << "not available for bosons" << endl;
      return -1;
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
  


//   if (TwoDTranslationFlag == false)
  {
    if (((SzSymmetrySector != 0) && (InputSpace->GetHilbertSpaceDimension() != InputStates[0].GetVectorDimension()))
	|| ((SzSymmetrySector == 0) && (OutputSpace->GetHilbertSpaceDimension() != InputStates[0].GetVectorDimension())))
      {
	cout << "error, " << Manager.GetString("input-state")  << " has a wrong dimension (" << InputStates[0].GetVectorDimension() << ", should be " << InputSpace->GetHilbertSpaceDimension() << ")" << endl;
	return -1;
      }
  }

  if (SzSymmetrySector == 0)
    {
      char* SzSymmetryValueString = new char[256];
      sprintf (SzSymmetryValueString, "sz_0_szsym_%d_", InputSzSymmetrySector);

      for (int i = 0; i < NbrInputStates; ++i)
	{
	  char* VectorOutputName = ReplaceString(InputStateNames[i], "sz_0_", SzSymmetryValueString);
	  ComplexVector TmpVector;
	  if (TwoDTranslationFlag == false)
	    {
	    }
	  else
	    {
	      TmpVector = ((FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation*) InputSpace)->ConvertFromNbodyBasis (InputStates[i], ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*)OutputSpace));
	      cout << "dimension = " << TmpVector.GetVectorDimension() <<  " " << InputSpace->GetHilbertSpaceDimension() << endl;
	    }
	  if (TmpVector.WriteVector(VectorOutputName) == false)
	    {
	      cout << "error, can't write vector " << VectorOutputName << endl;
	    }
	  delete[] VectorOutputName;
	}
    }
  else
    {
      char* SzSymmetryValueString = new char[256];
      sprintf (SzSymmetryValueString, "_szsym_%d", SzSymmetrySector);
      
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  char* VectorOutputName = ReplaceString(InputStateNames[i], SzSymmetryValueString, "");
	  ComplexVector TmpVector;
	  if (TwoDTranslationFlag == false)
	    {
	    }
	  else
	    {
	      TmpVector = ((FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation*) InputSpace)->ConvertToNbodyBasis (InputStates[i], ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*)OutputSpace));
	      cout << "dimension = " << TmpVector.GetVectorDimension() <<  " " << OutputSpace->GetHilbertSpaceDimension() << endl;
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
}
