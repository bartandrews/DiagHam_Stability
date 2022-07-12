#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include "Operator/ParticleOnSphereWithSpinSuperconductorOrderParameterOperator.h"
#include "Operator/ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator.h"


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


// Diagionalize the superconductor order parameter matrix
//
// orderParameter = reference to the superconductor order parameter matrix
// output = reference on the output stream
void HubbardSuperconductorOrderParameterMatrixDiagonalize(ComplexMatrix& orderParameter, ostream& output);

// compute the Fourier transform of the superconductor order parameter
//
// nbrUnitCellX = number of unit cells in the x direction
// nbrUnitCellY = number of unit cells in the y direction
// nbrOrbitals = number of orbitals per unit cells
// spinIndex1 = spin index of the first creation operator
// spinIndex2 = spin index of the second creation operator
// orderParameter = reference to the superconductor order parameter matrix
// output = reference on the output stream
void HubbardSuperconductorOrderParameterFourrierTransform(int nbrUnitCellX, int nbrUnitCellY, int nbrOrbitals, int spinIndex1, int spinIndex2,
							  ComplexMatrix& orderParameter, ostream& output);


int main(int argc, char** argv)
{
  OptionManager Manager ("HubbardSuperconductorOrderParameter" , "0.01");
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

  (*SystemGroup) += new SingleStringOption  ('\n', "left-state", "name of the file corresponding to the state |Psi_L> (the order parameter being <Psi_L|c^+c^+|Psi_R>");
  (*SystemGroup) += new SingleStringOption  ('\n', "right-state", "name of the file corresponding to the state |Psi_R> (the order parameter being <Psi_L|c^+c^+|Psi_R>");   (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerate-leftstates", "single column file describing a set of degenerate left states");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerate-rightstates", "single column file describing a set of degenerate right states");
  (*SystemGroup) += new BooleanOption ('\n', "only-cc", "compute only the parameters c^+_sigma c^+_sigma' instead of their linear combinations");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardSuperconductorOrderParameter -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int LeftNbrParticles = 0;
  int LeftNbrSites = 0;
  bool LeftStatistics = true;
  bool LeftGutzwillerFlag = false;
  int NbrLeftStates = 0;
  int LeftMomentumFlag = false;
  int LeftKxMomentum = 0;
  int LeftXPeriodicity = 0;
  int LeftKyMomentum = 0;
  int LeftYPeriodicity = 0;
  if ((Manager.GetString("left-state") == 0) && (Manager.GetString("degenerate-leftstates") == 0))
    {
      cout << "error, a left state file should be provided. See man page for option syntax or type HubbardSuperconductorOrderParameter -h" << endl;
      return -1;
    }
  if (Manager.GetString("left-state") != 0)
    {
      NbrLeftStates = 1;
      if (IsFile(Manager.GetString("left-state")) == false)
	{
	  cout << "can't open file " << Manager.GetString("left-state") << endl;
	  return -1;
	}
      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("left-state"), LeftNbrParticles, LeftNbrSites, 
									   LeftKxMomentum, LeftKyMomentum, LeftXPeriodicity, LeftYPeriodicity, 
									   LeftStatistics, LeftGutzwillerFlag) == false)
	{
	  if (FTIHubbardModelFindSystemInfoFromVectorFileName(Manager.GetString("left-state"), LeftNbrParticles, LeftNbrSites, LeftStatistics, LeftGutzwillerFlag) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << Manager.GetString("left-state") << endl;
	      return -1;
	    }
	}
      else
	{
	  LeftMomentumFlag = true;
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-leftstates")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      NbrLeftStates = DegenerateFile.GetNbrLines();
      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), LeftNbrParticles, LeftNbrSites, 
									   LeftKxMomentum, LeftKyMomentum, LeftXPeriodicity, LeftYPeriodicity, 
									   LeftStatistics, LeftGutzwillerFlag) == false)
	{
	  if (FTIHubbardModelFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), LeftNbrParticles, LeftNbrSites, LeftStatistics, LeftGutzwillerFlag) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << DegenerateFile(0, 0) << endl;
	      return -1;
	    }
	}
      else
	{
	  LeftMomentumFlag = true;
	}
    }

  int RightNbrParticles = 0;
  int RightNbrSites = 0;
  bool RightStatistics = true;
  bool RightGutzwillerFlag = false;
  int NbrRightStates = 0;
  int RightMomentumFlag = false;
  int RightKxMomentum = 0;
  int RightXPeriodicity = 0;
  int RightKyMomentum = 0;
  int RightYPeriodicity = 0;
  if ((Manager.GetString("right-state") == 0) && (Manager.GetString("degenerate-rightstates") == 0))
    {
      cout << "error, a right state file should be provided. See man page for option syntax or type HubbardSuperconductorOrderParameter -h" << endl;
      return -1;
    }
  if (Manager.GetString("right-state") != 0)
    {
      NbrRightStates = 1;
      if (IsFile(Manager.GetString("right-state")) == false)
	{
	  cout << "can't open file " << Manager.GetString("right-state") << endl;
	  return -1;
	}
      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("right-state"), RightNbrParticles, RightNbrSites, 
									   RightKxMomentum, RightKyMomentum, RightXPeriodicity, RightYPeriodicity, 
									   RightStatistics, RightGutzwillerFlag) == false)
	{
	  if (FTIHubbardModelFindSystemInfoFromVectorFileName(Manager.GetString("right-state"), RightNbrParticles, RightNbrSites, RightStatistics, RightGutzwillerFlag) == false)
	    {
	      cout << "error while retrieving system parameters from file name " <<Manager.GetString("right-state")  << endl;
	      return -1;
	    }
	}
      else
	{
	  RightMomentumFlag = true;
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-rightstates")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      NbrRightStates = DegenerateFile.GetNbrLines();
      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), RightNbrParticles, RightNbrSites, 
									   RightKxMomentum, RightKyMomentum, RightXPeriodicity, RightYPeriodicity, 
									   RightStatistics, RightGutzwillerFlag) == false)
	{
	  if (FTIHubbardModelFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), RightNbrParticles, RightNbrSites, RightStatistics, RightGutzwillerFlag) == false)
	    {
	      cout << "error while retrieving system parameters from file name " <<  DegenerateFile(0, 0) << endl;
	      return -1;
	    }
	}
      else
	{
	  RightMomentumFlag = true;
	}
    }

  if (RightNbrSites != LeftNbrSites)
    {
      cout << "error, left and right states don't have the same number of sites" << endl;
      return -1;      
    }

  if (LeftNbrParticles != (RightNbrParticles + 2))
    {
      cout << "error, left and right states don't have the proper number of particles" << endl;
      return -1;      
    }
  if (RightMomentumFlag != LeftMomentumFlag)
    {
      cout << "error, left and right states should be both written as momentum eigenstate" << endl;
      return -1;      
    }
  if ((RightMomentumFlag == true) && ((RightXPeriodicity != LeftXPeriodicity) || (RightYPeriodicity != LeftYPeriodicity)))
    {
      cout << "error, left and right states are not defined on the same lattice (x-periodicity = " << LeftXPeriodicity << " vs " <<  RightXPeriodicity 
	   << ", y-periodicity = " << LeftYPeriodicity << " vs " <<  RightYPeriodicity << ")" << endl;
      return -1;      
    }
      
  ComplexVector* LeftStates = new ComplexVector[NbrLeftStates];
  if (Manager.GetString("left-state") != 0)
    {
      if (LeftStates[0].ReadVector (Manager.GetString("left-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("left-state") << endl;
	  return -1;      
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-leftstates")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      if (LeftStates[0].ReadVector (DegenerateFile(0, 0)) == false)
	{
	  cout << "can't open vector file " << DegenerateFile(0, 0) << endl;
	  return -1;      
	}	  
      for (int i = 1; i < NbrLeftStates; ++i)
	{
	  if (LeftStates[i].ReadVector (DegenerateFile(0, i)) == false)
	    {
	      cout << "can't open vector file " << DegenerateFile(0, i) << endl;
	      return -1;      
	    }	  
	  if (LeftStates[0].GetVectorDimension() != LeftStates[i].GetVectorDimension())
	    {
	      cout << "error, " << DegenerateFile(0, 0) << " and " <<  DegenerateFile(0, i) << "don't have the same  dimension (" << LeftStates[0].GetVectorDimension() << " and " << LeftStates[i].GetVectorDimension()<< ")" << endl;
	      return -1;
	    }
	}
    }
  ParticleOnSphereWithSpin* LeftSpace = 0;
  if (LeftMomentumFlag == false)
    {
      if (LeftStatistics == true)
	{
	  if (LeftGutzwillerFlag == false)
	    LeftSpace = new FermionOnLatticeWithSpinRealSpace (LeftNbrParticles, LeftNbrSites);
	  else
	    LeftSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (LeftNbrParticles, LeftNbrSites);
	}
      else
	{
	  cout << "not available for bosons" << endl;
	  return -1;
	}
    }
  else
    {
      if (LeftStatistics == true)
	{
	  if (LeftGutzwillerFlag == false)
	    LeftSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (LeftNbrParticles, LeftNbrSites, LeftKxMomentum, LeftXPeriodicity, LeftKyMomentum, LeftYPeriodicity);
	  else
	    LeftSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (LeftNbrParticles, LeftNbrSites, LeftKxMomentum, LeftXPeriodicity, LeftKyMomentum, LeftYPeriodicity);
	}
      else
	{
	  cout << "not available for bosons" << endl;
	  return -1;
	}
    }

  if (LeftSpace->GetHilbertSpaceDimension() != LeftStates[0].GetVectorDimension())
    {
      cout << "error, " << Manager.GetString("left-state")  << " has a wrong dimension (" <<LeftStates[0].GetVectorDimension() << ", should be " << LeftSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  

  ComplexVector* RightStates = new ComplexVector[NbrRightStates];
  if (Manager.GetString("right-state") != 0)
    {
      if (RightStates[0].ReadVector (Manager.GetString("right-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("right-state") << endl;
	  return -1;      
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-rightstates")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      if (RightStates[0].ReadVector (DegenerateFile(0, 0)) == false)
	{
	  cout << "can't open vector file " << DegenerateFile(0, 0) << endl;
	  return -1;      
	}	  
      for (int i = 1; i < NbrRightStates; ++i)
	{
	  if (RightStates[i].ReadVector (DegenerateFile(0, i)) == false)
	    {
	      cout << "can't open vector file " << DegenerateFile(0, i) << endl;
	      return -1;      
	    }	  
	  if (RightStates[0].GetVectorDimension() != RightStates[i].GetVectorDimension())
	    {
	      cout << "error, " << DegenerateFile(0, 0) << " and " <<  DegenerateFile(0, i) << "don't have the same  dimension (" << RightStates[0].GetVectorDimension() << " and " << RightStates[i].GetVectorDimension()<< ")" << endl;
	      return -1;
	    }
	}
    }

  ParticleOnSphereWithSpin* RightSpace = 0;
  if (RightMomentumFlag == false)
    {
      if (RightStatistics == true)
	{
	  if (RightGutzwillerFlag == false)
	    RightSpace = new FermionOnLatticeWithSpinRealSpace (RightNbrParticles, RightNbrSites);
	  else
	    RightSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (RightNbrParticles, RightNbrSites);
	}
      else
	{
	  cout << "not available for bosons" << endl;
	  return -1;
	}
    }
  else
    {
      if (RightStatistics == true)
	{
	  if (RightGutzwillerFlag == false)
	    RightSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (RightNbrParticles, RightNbrSites, 
										RightKxMomentum, RightXPeriodicity, RightKyMomentum, RightYPeriodicity);
	  else
	    RightSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (RightNbrParticles, RightNbrSites, 
												       RightKxMomentum, RightXPeriodicity, RightKyMomentum, RightYPeriodicity);
	}
      else
	{
	  cout << "not available for bosons" << endl;
	  return -1;
	}
    }
  if (RightSpace->GetHilbertSpaceDimension() != RightStates[0].GetVectorDimension())
    {
      cout << "error, " << Manager.GetString("right-state")  << " has a wrong dimension (" << RightStates[0].GetVectorDimension() << ", should be " << RightSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }

  RightSpace->SetTargetSpace(LeftSpace);
    
  ofstream File;
  char* OutputFileName;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      File.open(OutputFileName, ios::binary | ios::out);
    }
  else
    {
      if (Manager.GetString("left-state") != 0)
	{
	  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("left-state"), "vec", "orderparam.dat");
	  if (OutputFileName == 0)
	    {
	      cout << "no vec extension was find in " << Manager.GetString("left-state") << " file name" << endl;
	      return 0;
	    }
	  File.open(OutputFileName, ios::binary | ios::out);
	}
      else
	{
	  MultiColumnASCIIFile DegenerateFile;
	  if (DegenerateFile.Parse(Manager.GetString("degenerate-leftstates")) == false)
	    {
	      DegenerateFile.DumpErrors(cout);
	      return -1;
	    }
	  OutputFileName = ReplaceExtensionToFileName(DegenerateFile(0, 0), "vec", "orderparam.dat");
	  if (OutputFileName == 0)
	    {
	      cout << "no vec extension was find in " << DegenerateFile(0, 0) << " file name" << endl;
	      return 0;
	    }
	  File.open(OutputFileName, ios::binary | ios::out);
	}
    }
  File.precision(14);
  cout.precision(14);

  ofstream FileFourierTransform;
  if ((RightMomentumFlag == true) && (LeftMomentumFlag == true))
    {
      char* TmpFileName = ReplaceExtensionToFileName(OutputFileName, "dat", "fourier.dat");
      if (TmpFileName == 0)
	{
	  cout << "no dat extension was find in " << OutputFileName << " file name" << endl;
	  return 0;
	}
      FileFourierTransform.open(TmpFileName, ios::binary | ios::out);
      FileFourierTransform.precision(14);
    }

  if (Manager.GetBoolean("only-cc") == true)
    {
      if ((NbrLeftStates == 1) && (NbrRightStates == 1))
	{
	  File << "# <Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R> with sigma,sigma' = 0 (down) or 1 (up)" << endl
	       << "# i j sigma sigma' |<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>|^2 |<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>| Arg(<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>)" << endl;
	}
      else
	{
	  File << "# <Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R> with sigma,sigma' = 0 (down) or 1 (up)" << endl
	       << "# i j sigma sigma' |<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>|^2" << endl;
	}
      double NormalizationFactor = 1.0 / sqrt((double) (NbrLeftStates * NbrRightStates));
      for (int i = 0; i < RightNbrSites; ++i)
	{
	  int Index = i;
	  if (RightStatistics == true)
	    {
	      if ((RightGutzwillerFlag == false) && (LeftGutzwillerFlag == false))
		{
		  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorDownUpDiag = 0;
		  if (RightMomentumFlag == false)
		    OperatorDownUpDiag = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator (RightSpace, i, 0, i, 1);
		  else
		    OperatorDownUpDiag = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 0, i, 1);
		  ComplexMatrix TmpMatrix (NbrLeftStates, NbrRightStates);
		  for (int k = 0; k < NbrLeftStates; ++k)
		    for (int l = 0; l < NbrRightStates; ++l)
		      {
			OperatorMatrixElementOperation OperationDownUp(OperatorDownUpDiag, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
			OperationDownUp.ApplyOperation(Architecture.GetArchitecture());
			TmpMatrix[l][k] = OperationDownUp.GetScalar() * NormalizationFactor;
		      }
		  File << i << " " << i << " 0 1";
		  HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
		  File << endl;
		  delete OperatorDownUpDiag;
		  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorUpDownDiag = 0;
		  if (RightMomentumFlag == false)
		    OperatorUpDownDiag = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator (RightSpace, i, 1, i, 0);
		  else
		    OperatorUpDownDiag = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, i, 0);
		  for (int k = 0; k < NbrLeftStates; ++k)
		    for (int l = 0; l < NbrRightStates; ++l)
		      {
			OperatorMatrixElementOperation OperationUpDown(OperatorUpDownDiag, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
			OperationUpDown.ApplyOperation(Architecture.GetArchitecture());
			TmpMatrix[l][k] = OperationUpDown.GetScalar() * NormalizationFactor;
		      }
		  File << i << " " << i << " 1 0";
		  HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
		  File << endl;
		  delete OperatorUpDownDiag;
		}
	      ++Index;
	    }
	  for (int j = Index; j < RightNbrSites; ++j)
	    {
	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorDownDown = 0;
	      if (RightMomentumFlag == false)
		OperatorDownDown = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 0, j, 0);
	      else
		OperatorDownDown = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 0, j, 0);
	      ComplexMatrix TmpMatrix (NbrLeftStates, NbrRightStates);
	      for (int k = 0; k < NbrLeftStates; ++k)
		for (int l = 0; l < NbrRightStates; ++l)
		  {
		    OperatorMatrixElementOperation OperationDownDown(OperatorDownDown, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
		    OperationDownDown.ApplyOperation(Architecture.GetArchitecture());
		    TmpMatrix[l][k] = OperationDownDown.GetScalar() * NormalizationFactor;
		  }
	      File << i << " " << j << " 0 0";
	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
	      File << endl;
	      delete OperatorDownDown;
	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorDownUp = 0;
	      if (RightMomentumFlag == false)
		OperatorDownUp = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 0, j, 1);
	      else
		OperatorDownUp = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 0, j, 1);
	      for (int k = 0; k < NbrLeftStates; ++k)
		for (int l = 0; l < NbrRightStates; ++l)
		  {
		    OperatorMatrixElementOperation OperationDownUp(OperatorDownUp, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
		    OperationDownUp.ApplyOperation(Architecture.GetArchitecture());
		    TmpMatrix[l][k] = OperationDownUp.GetScalar() * NormalizationFactor;
		  }
	      File << i << " " << j << " 0 1";
	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
	      File << endl;
	      delete OperatorDownUp;
	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorUpDown = 0;
	      if (RightMomentumFlag == false)
		OperatorUpDown = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 1, j, 0);
	      else
		OperatorUpDown = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, j, 0);
	      for (int k = 0; k < NbrLeftStates; ++k)
		for (int l = 0; l < NbrRightStates; ++l)
		  {
		    OperatorMatrixElementOperation OperationUpDown(OperatorUpDown, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
		    OperationUpDown.ApplyOperation(Architecture.GetArchitecture());
		    TmpMatrix[l][k] = OperationUpDown.GetScalar() * NormalizationFactor;
		  }
	      File << i << " " << j << " 1 0";
	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
	      File << endl;
	      delete OperatorUpDown;

	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorUpUp = 0;
	      if (RightMomentumFlag == false)
		OperatorUpUp = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 1, j, 1);
	      else
		OperatorUpUp = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, j, 1);
	      for (int k = 0; k < NbrLeftStates; ++k)
		for (int l = 0; l < NbrRightStates; ++l)
		  {
		    OperatorMatrixElementOperation OperationUpUp(OperatorUpUp, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
		    OperationUpUp.ApplyOperation(Architecture.GetArchitecture());
		    TmpMatrix[l][k] = OperationUpUp.GetScalar() * NormalizationFactor;
		  }
	      File << i << " " << j << " 1 1";
	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
	      File << endl;
	      delete OperatorUpUp;
	    }
	}
    }
  else
    {
      File << "# <Psi_L| c^+_{i,sigma} c^+_{j,sigma'} +/- c^+_{i,sigma} c^+_{j,sigma'}|Psi_R> with sigma,sigma' = 0 (down) or 1 (up)" << endl
	   << "# for each case is given the (norm)^2, the norm and the argument" << endl;
      File << "# i j ";
      if ((NbrLeftStates == 1) && (NbrRightStates == 1))
	{
	  File << "<Psi_L| (c^+_{i,up} c^+_{j,up} + c^+_{i,down} c^+_{j,down} |Psi_R> <Psi_L| (c^+_{i,up} c^+_{j,up} - c^+_{i,down} c^+_{j,down} |Psi_R> <Psi_L| (c^+_{i,up} c^+_{j,down} + c^+_{i,up} c^+_{j,down} |Psi_R> <Psi_L| (c^+_{i,down} c^+_{j,up} - c^+_{i,down} c^+_{j,up} |Psi_R>" << endl;
	}
      else
	{
	  File << "|<Psi_L| (c^+_{i,up} c^+_{j,up} + c^+_{i,down} c^+_{j,down} |Psi_R> |^2 |<Psi_L| (c^+_{i,up} c^+_{j,up} - c^+_{i,down} c^+_{j,down} |Psi_R> |^2 |<Psi_L| (c^+_{i,up} c^+_{j,down} + c^+_{i,down} c^+_{j,up} |Psi_R> |^2 |<Psi_L| (c^+_{i,up} c^+_{j,down} - c^+_{i,down} c^+_{j,up} |Psi_R> |^2" << endl;
	}
      double NormalizationFactor = 1.0 / sqrt((double) (NbrLeftStates * NbrRightStates));
      for (int i = 0; i < RightNbrSites; ++i)
	{
	  int Index = i;
	  if (RightStatistics == true)
	    {
	      if ((RightGutzwillerFlag == false) && (LeftGutzwillerFlag == false))
		{
		  ComplexMatrix TmpMatrix (NbrLeftStates, NbrRightStates);
		  File << i << " " << i << " (0,0) 0 0 (0,0) 0 0";

		  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorDownUpDiag = 0;
		  if (RightMomentumFlag == false)
		    OperatorDownUpDiag = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 0, i, 1, 1, 0, 1.0);
		  else
		    OperatorDownUpDiag = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 0, i, 1, 1, 0, 1.0);
		  for (int k = 0; k < NbrLeftStates; ++k)
		    for (int l = 0; l < NbrRightStates; ++l)
		      {
			OperatorMatrixElementOperation OperationDownUp(OperatorDownUpDiag, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
			OperationDownUp.ApplyOperation(Architecture.GetArchitecture());
			TmpMatrix[l][k] = OperationDownUp.GetScalar() * NormalizationFactor;
		      }
		  HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
		  delete OperatorDownUpDiag;

		  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorUpDownDiag = 0;
		  if (RightMomentumFlag == false)
		    OperatorUpDownDiag = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 1, i, 0, 1, 0, 1.0);
		  else
		    OperatorUpDownDiag = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, i, 0, 1, 0, 1.0);
		  for (int k = 0; k < NbrLeftStates; ++k)
		    for (int l = 0; l < NbrRightStates; ++l)
		      {
			OperatorMatrixElementOperation OperationUpDown (OperatorUpDownDiag, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
			OperationUpDown.ApplyOperation(Architecture.GetArchitecture());
			TmpMatrix[l][k] = OperationUpDown.GetScalar() * NormalizationFactor;
		      }
		  HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
		  delete OperatorUpDownDiag;
		  File << endl;
		}
	      ++Index;
	    }
	  for (int j = Index; j < RightNbrSites; ++j)
	    {
	      ComplexMatrix TmpMatrix (NbrLeftStates, NbrRightStates);
	      File << i << " " << j << " ";

	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorUpUp = 0;
	      if (RightMomentumFlag == false)
		OperatorUpUp = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 1, j, 1, 0, 0, 1.0);
	      else
		OperatorUpUp = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, j, 1, 0, 0, 1.0);
	      for (int k = 0; k < NbrLeftStates; ++k)
		for (int l = 0; l < NbrRightStates; ++l)
		  {
		    OperatorMatrixElementOperation OperationUpUp(OperatorUpUp, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
		    OperationUpUp.ApplyOperation(Architecture.GetArchitecture());
		    TmpMatrix[l][k] = OperationUpUp.GetScalar() * NormalizationFactor;
		  }
	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
	      delete OperatorUpUp;

	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorDownDown = 0;
	      if (RightMomentumFlag == false)
		OperatorDownDown = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 1, j, 1, 0, 0, -1.0);
	      else
		OperatorDownDown = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, j, 1, 0, 0, -1.0);
	      for (int k = 0; k < NbrLeftStates; ++k)
		for (int l = 0; l < NbrRightStates; ++l)
		  {
		    OperatorMatrixElementOperation OperationDownDown(OperatorDownDown, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
		    OperationDownDown.ApplyOperation(Architecture.GetArchitecture());
		    TmpMatrix[l][k] = OperationDownDown.GetScalar() * NormalizationFactor;
		  }
	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
	      delete OperatorDownDown;

	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorDownUp = 0;
	      if (RightMomentumFlag == false)
		OperatorDownUp = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 1, j, 0, 0, 1, 1.0);
	      else
		OperatorDownUp = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, j, 0, 0, 1, 1.0);
	      for (int k = 0; k < NbrLeftStates; ++k)
		for (int l = 0; l < NbrRightStates; ++l)
		  {
		    OperatorMatrixElementOperation OperationDownUp(OperatorDownUp, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
		    OperationDownUp.ApplyOperation(Architecture.GetArchitecture());
		    TmpMatrix[l][k] = OperationDownUp.GetScalar() * NormalizationFactor;
		  }
	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
	      delete OperatorDownUp;

	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorUpDown = 0;
	      if (RightMomentumFlag == false)
		OperatorUpDown = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 1, j, 0, 0, 1, -1.0);
	      else
		OperatorUpDown = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, j, 0, 0, 1, -1.0);
	      for (int k = 0; k < NbrLeftStates; ++k)
		for (int l = 0; l < NbrRightStates; ++l)
		  {
		    OperatorMatrixElementOperation OperationUpDown(OperatorUpDown, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
		    OperationUpDown.ApplyOperation(Architecture.GetArchitecture());
		    TmpMatrix[l][k] = OperationUpDown.GetScalar() * NormalizationFactor;
		  }
	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
	      delete OperatorUpDown;
	      File << endl;
	    }
	}
    }

  File.close();
  if ((RightMomentumFlag == true) && (LeftMomentumFlag == true))
    {
      FileFourierTransform.close();
    }
  return 0;
}

// Diagionalize the superconductor order parameter matrix
//
// orderParameter = reference to the superconductor order parameter matrix
// output = reference on the output stream

void HubbardSuperconductorOrderParameterMatrixDiagonalize(ComplexMatrix& orderParameter, ostream& output)
{
  if ((orderParameter.GetNbrRow() > 1) || (orderParameter.GetNbrColumn() > 1))
    {
      ComplexMatrix TmpConjugate;
      TmpConjugate.Copy(orderParameter);
      TmpConjugate.HermitianTranspose();
      if (orderParameter.GetNbrRow() >= orderParameter.GetNbrColumn())
	{      
	  HermitianMatrix TmpMatrix (TmpConjugate * orderParameter);
	  RealDiagonalMatrix TmpDiag (TmpMatrix.GetNbrRow());
#ifdef __LAPACK__
	  TmpMatrix.LapackDiagonalize(TmpDiag);
#else
	  TmpMatrix.Diagonalize(TmpDiag);
#endif		  
	  TmpDiag.SortMatrixDownOrder();
	  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
	    output << " " << TmpDiag[i];
	}
      else
	{
	  HermitianMatrix TmpMatrix2 (orderParameter * TmpConjugate);
	  RealDiagonalMatrix TmpDiag2 (TmpMatrix2.GetNbrRow());
#ifdef __LAPACK__
	  TmpMatrix2.LapackDiagonalize(TmpDiag2);
#else
	  TmpMatrix2.Diagonalize(TmpDiag2);
#endif		  
	  TmpDiag2.SortMatrixDownOrder();
	  for (int i = 0; i < TmpDiag2.GetNbrRow(); ++i)
	    output << " " << TmpDiag2[i];
	}
    }
  else
    {
      output << " " << SqrNorm(orderParameter[0][0]) << " " << Norm(orderParameter[0][0]) << " " << Arg(orderParameter[0][0]);
    }
}


// compute the Fourier transform of the superconductor order parameter
//
// nbrUnitCellX = number of unit cells in the x direction
// nbrUnitCellY = number of unit cells in the y direction
// nbrOrbitals = number of orbitals per unit cells
// spinIndex1 = spin index of the first creation operator
// spinIndex2 = spin index of the second creation operator
// orderParameter = reference to the superconductor order parameter matrix
// output = reference on the output stream

void HubbardSuperconductorOrderParameterFourrierTransform(int nbrUnitCellX, int nbrUnitCellY, int nbrOrbitals, int spinIndex1, int spinIndex2,
							  ComplexMatrix& orderParameter, ostream& output)
{
  int TotalNbrUnitCells = nbrUnitCellX * nbrUnitCellY;
  Complex** TmpPhases = new Complex*[nbrUnitCellX];
  for (int Kx1 = 0; Kx1 < nbrUnitCellX; ++Kx1)
    {
      for (int Ky1 = 0; Ky1 < nbrUnitCellY; ++Ky1)
	{
	  ComplexMatrix TmpMatrix1;
	  TmpMatrix1.Copy(orderParameter);
	  for (int x = 0; x < nbrUnitCellX; ++x)
	    {
	      for (int y = 0; y < nbrUnitCellY; ++y)
		{
		  Complex TmpPhase = Phase (2.0 * M_PI * ((((double) (Kx1 * x)) / ((double) nbrUnitCellX)) + (((double) (Ky1 * y)) / ((double) nbrUnitCellY))));
		  for (int i = 0; i < nbrOrbitals; ++i)
		    {
		      int TmpIndex = i + ((x * nbrUnitCellY) + y) * nbrOrbitals;
		      for (int j = 0; j < TmpMatrix1.GetNbrRow(); ++j)
			{
			  TmpMatrix1[TmpIndex][j] *= TmpPhase;
			}
		    }
		}
	    }
	  for (int Kx2 = 0; Kx2 < nbrUnitCellX; ++Kx2)
	    {
	      for (int Ky2 = 0; Ky2 < nbrUnitCellY; ++Ky2)
		{
		  ComplexMatrix TmpMatrix2;
		  TmpMatrix2.Copy(TmpMatrix1);
		  for (int x = 0; x < nbrUnitCellX; ++x)
		    {
		      for (int y = 0; y < nbrUnitCellY; ++y)
			{
			  Complex TmpPhase = Phase (2.0 * M_PI * ((((double) (Kx1 * x)) / ((double) nbrUnitCellX)) + (((double) (Ky1 * y)) / ((double) nbrUnitCellY))));
			  for (int i = 0; i < nbrOrbitals; ++i)
			    {
			      int TmpIndex = i + ((x * nbrUnitCellY) + y) * nbrOrbitals;
			      for (int j = 0; j < TmpMatrix2.GetNbrColumn(); ++j)
				{
				  TmpMatrix2[j][TmpIndex] *= TmpPhase;
				}
			    }
			}
		    }
		  for (int i = 0; i < nbrOrbitals; ++i)
		    {
		      for (int j = 0; j < nbrOrbitals; ++j)
			{
			  Complex Tmp = 0.0;
			  for (int k = i; k < TmpMatrix2.GetNbrRow(); k += nbrOrbitals)
			    {
			      for (int l = j; l < TmpMatrix2.GetNbrColumn(); l += nbrOrbitals)
				{				  
				  Tmp += TmpMatrix2[k][l];
				}
			    }
			  output << Kx1 << " " << Ky1 << " " << Kx2 << " " << Ky2 << " " << spinIndex1 << " " << spinIndex2 << " " << i << " " << j << " " << Tmp.Re << " " << Tmp.Im << endl;
			}
		    }
		}	      
	    }
	}
    }
}
