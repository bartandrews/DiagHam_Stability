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
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong.h"

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
  OptionManager Manager ("FTIMinNbrSingletProjection" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption ('\n', "min-nbrsinglets", "mininum number of on-site singlet that defines the projected space", 0);
  (*SystemGroup) += new BooleanOption ('\n', "only-weights", "do not save the projected states, only compute their weights in the unprojected basis", 0);  
  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIMinNbrSingletProjection -h" << endl;
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
  int InputMinNbrSinglets = 0;

  if ((Manager.GetString("input-state") == 0) && (Manager.GetString("degenerate-states") == 0))
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type FTIMinNbrSingletProjection -h" << endl;
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
      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites, SzValue, SzSymmetrySector, InputMinNbrSinglets,
									   XMomentum, YMomentum, NbrSiteX, NbrSiteY, Statistics, GutzwillerFlag) == false)
	{
	  if (FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites, SzValue, SzSymmetrySector,
								    Statistics, GutzwillerFlag) == false)
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
      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), NbrParticles, NbrSites, SzValue, SzSymmetrySector, InputMinNbrSinglets,
									   XMomentum, YMomentum, NbrSiteX, NbrSiteY, Statistics, GutzwillerFlag) == false)
	{
	  if (FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), NbrParticles, NbrSites, SzValue, SzSymmetrySector,
								    Statistics, GutzwillerFlag) == false)
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


  int OutputMinNbrSinglets = Manager.GetInteger("min-nbrsinglets");
  if ((OutputMinNbrSinglets * 2) > NbrParticles)
    {
      cout << "error, the minimal number of on-site singlets cannot be higher than " << (NbrParticles / 2) << endl;
      return -1;
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


  char* MinNbrSingletsString = new char[256];
  char* NbrParticlesString = new char[256];
  sprintf (MinNbrSingletsString, "_minnbrsinglet_%d_n_%d_", OutputMinNbrSinglets, NbrParticles);
  if (InputMinNbrSinglets != 0)
    {
      sprintf (NbrParticlesString, "_minnbrsinglet_%d_n_%d_", InputMinNbrSinglets, NbrParticles);
    }
  else
    {
      sprintf (NbrParticlesString, "_n_%d_", NbrParticles);
    }

  //  ParticleOnSphereWithSpin* InputSpace = 0;
  //  FermionOnLatticeWithSpinRealSpace* InputSpace = 0;
  if (Statistics == true)
    {
      if (TwoDTranslationFlag == false)
	{
// 	  if (GutzwillerFlag == false)
// 	    {
// 	      if (SzSymmetrySector == 0)
// 		{
// 		  InputSpace = new FermionOnLatticeWithSpinRealSpace (NbrParticles, SzValue, NbrSites, 10000000);
// 		}
// 	      else
// 		{
// 		  InputSpace = new FermionOnLatticeWithSpinSzSymmetryRealSpace (NbrParticles, SzValue, NbrSites, (SzSymmetrySector == -1), 10000000);
// 		}
// 	    }
// 	  else
// 	    {
// 	      if (SzSymmetrySector == 0)
// 		{
// 		  InputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, SzValue, NbrSites, 10000000);
// 		}
// 	      else
// 		{
// 		  InputSpace = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace (NbrParticles, SzValue, NbrSites, (SzSymmetrySector == -1), 10000000);
// 		}
// 	    }
	}
      else
	{
	  ParticleOnTorusWithSpinAndMagneticTranslations* InputSpace = 0;
	  ParticleOnTorusWithSpinAndMagneticTranslations* OutputSpace  = 0;
// 	  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* InputSpace = 0;
// 	  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* OutputSpace = 0;
	  if (GutzwillerFlag == false)
	    {
	      if (SzSymmetrySector == 0)
		{
		  if (InputMinNbrSinglets == 0)
		    {
#ifdef __64_BITS__
		      if (NbrSites < 31)
#else
		      if (NbrSites < 15)
#endif
			{
			  InputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, SzValue, NbrSites, 
											      XMomentum, NbrSiteX, YMomentum, NbrSiteY);
			}
		      else
			{
			  InputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong (NbrParticles, SzValue, NbrSites, 
												  XMomentum, NbrSiteX, YMomentum, NbrSiteY);
			}
		    }
		  else
		    {
#ifdef __64_BITS__
		      if (NbrSites < 31)
#else
		      if (NbrSites < 15)
#endif
			{
			  InputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets (NbrParticles, InputMinNbrSinglets, SzValue, NbrSites, 
													    XMomentum, NbrSiteX, YMomentum, NbrSiteY);
			}
		      else
			{
			  InputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong (NbrParticles, InputMinNbrSinglets, SzValue, NbrSites, 
														XMomentum, NbrSiteX, YMomentum, NbrSiteY);
			}
		    }
#ifdef __64_BITS__
		  if (NbrSites < 31)
#else
		    if (NbrSites < 15)
#endif
		      {
			OutputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets (NbrParticles, OutputMinNbrSinglets,
													   SzValue, NbrSites, XMomentum, NbrSiteX, YMomentum, NbrSiteY);
		      }
		    else
		      {
			OutputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong (NbrParticles, OutputMinNbrSinglets,
													       SzValue, NbrSites, XMomentum, NbrSiteX, YMomentum, NbrSiteY);
		      }
		}
	      else
		{
		  if (InputMinNbrSinglets == 0)
		    {
#ifdef __64_BITS__
		  if (NbrSites < 31)
#else
		    if (NbrSites < 15)
#endif
		      {
			InputSpace = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (NbrParticles, SzValue, NbrSites, (SzSymmetrySector == -1),
												      XMomentum, NbrSiteX, YMomentum, NbrSiteY);
		      }
		    else
		      {
			InputSpace = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong (NbrParticles, SzValue, NbrSites, (SzSymmetrySector == -1),
													  XMomentum, NbrSiteX, YMomentum, NbrSiteY);
		      }
		    }
		  else
		    {
#ifdef __64_BITS__
		  if (NbrSites < 31)
#else
		    if (NbrSites < 15)
#endif
		      {
			InputSpace = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets (NbrParticles, InputMinNbrSinglets, SzValue, 
														    NbrSites, (SzSymmetrySector == -1),
														    XMomentum, NbrSiteX, YMomentum, NbrSiteY);
		      }
		    else
		      {
			InputSpace = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong (NbrParticles, InputMinNbrSinglets, SzValue, 
															NbrSites, (SzSymmetrySector == -1),
															XMomentum, NbrSiteX, YMomentum, NbrSiteY);
		      }

		    }
#ifdef __64_BITS__
		  if (NbrSites < 31)
#else
		    if (NbrSites < 15)
#endif
		      {
			OutputSpace = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets (NbrParticles, OutputMinNbrSinglets, SzValue, NbrSites, (SzSymmetrySector == -1),
														     XMomentum, NbrSiteX, YMomentum, NbrSiteY);
		      }
		    else
		      {
			OutputSpace = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong (NbrParticles, OutputMinNbrSinglets, SzValue, NbrSites, (SzSymmetrySector == -1),
															 XMomentum, NbrSiteX, YMomentum, NbrSiteY);
		      }
		}
	    }
	  else
	    {
	      if (SzSymmetrySector == 0)
		{
// 		  InputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, SzValue, NbrSites,
// 													     XMomentum, NbrSiteX, YMomentum, NbrSiteY, 10000000);
		}
	      else
		{
// 		  InputSpace = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, SzValue, NbrSites, (SzSymmetrySector == -1),
// 														       XMomentum, NbrSiteX, YMomentum, NbrSiteY, 10000000);
		}
	    }
	  for (int i = 0; i < NbrInputStates; ++i)
	    {
	      if (InputStates[i].GetVectorDimension() != InputSpace->GetHilbertSpaceDimension())
		{
		  cout << "error, " << InputStateNames[i] << " does not have the correct dimension (is " << InputStates[i].GetVectorDimension()
		       << ", should be " << InputSpace->GetHilbertSpaceDimension() << " )" << endl;
		    return -1;
		}
	      ComplexVector TmpVector;
#ifdef __64_BITS__
	      if (NbrSites < 31)
#else
		if (NbrSites < 15)
#endif
		  {
		    TmpVector = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) OutputSpace)->FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::ConvertToNbodyBasis(InputStates[i], 
															(FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) InputSpace);
		  }
		else
		  {
		    TmpVector = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong*) OutputSpace)->FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong::ConvertToNbodyBasis(InputStates[i], 
															    (FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong*) InputSpace);
		  }
	      double TmpSqrNorm = TmpVector.SqrNorm(); 
	      cout << InputStateNames[i] << " : " << (TmpVector.SqrNorm()) << endl;
	      if (Manager.GetBoolean("only-weights") == false)
		{
		  if (TmpSqrNorm != 0.0)
		    {
		      char* VectorOutputName = ReplaceString(InputStateNames[i], NbrParticlesString, MinNbrSingletsString);
		      TmpVector /= sqrt(TmpSqrNorm);
		      if (TmpVector.WriteVector(VectorOutputName) == false)
			{
			  cout << "error, can't write " << VectorOutputName << endl;
			  return -1;
			}
		    }
		}
	    }
	}
    }
  else
    {
      cout << "not available for bosons" << endl;
      return -1;
    }
  

  return 0;
}
