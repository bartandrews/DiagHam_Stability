#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"


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
  OptionManager Manager ("FTIRotateOneBodyBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-file", "name of the file corresponding to the state to be projected");
  (*SystemGroup) += new SingleStringOption  ('t', "transformation-file", "name of the ASCIII file provinding the transformation matrices");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name to store the projected state");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIRotateOneBodyBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int TotalKx = 0;
  int TotalKy = 0;
  int NbrParticles = 0;
  int NbrSitesX = 0;
  int NbrSitesY = 0;
  bool Statistics = true;
  int TotalSpin = 0;

  
  if ( Manager.GetString("input-file") == 0)
    {
      cout << "error, a input state file should be provided. See man page for option syntax or type  FTIRotateOneBodyBasis -h" << endl;
      return -1;
    }
  if ((Manager.GetString("input-file") != 0) &&  (IsFile(Manager.GetString("input-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("input-file") << endl;
      return -1;
    }
  char* OutputFileName = 0;
  if (Manager.GetString("output-file") == 0)
    {
      OutputFileName = ReplaceString(Manager.GetString("input-file"), "twoband", "rotated_twoband");
      if (OutputFileName == 0)
	{
	  cout << "can't guess output name from " << Manager.GetString("input-file") << "(should contain twoband)" << endl;
	  return 0;
	}
    }
  else
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(Manager.GetString("input-file"), NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " <<  Manager.GetString("input-file") << endl;
      return -1;
    }
  
  
  ComplexVector InitialState;
  if (InitialState.ReadVectorTest(Manager.GetString("input-file")) == false)
    {
      RealVector RealInitialState;
      if (RealInitialState.ReadVector(Manager.GetString("input-file")) == false)
	{
	  cout << "can't open vector file " <<  Manager.GetString("input-file") << endl;
	  return -1;      
	}
      InitialState = ComplexVector(RealInitialState);
    }
  else
    {  
      if (InitialState.ReadVector(Manager.GetString("input-file")) == false)
	{
	  cout << "can't open vector file " <<  Manager.GetString("input-file") << endl;
	  return -1;      
	}
    }

  Abstract2DTightBindingModel* TightBindingModel;
  double* DummyChemicalPotentials = new double[2];
  DummyChemicalPotentials[0] = 0.0;
  DummyChemicalPotentials[1] = 0.0;
  TightBindingModel = new TightBindingModel2DAtomicLimitLattice (NbrSitesX, NbrSitesY, 2, DummyChemicalPotentials,
								 0.0, 0.0, Architecture.GetArchitecture(), true);

  MultiColumnASCIIFile TransformationMatrixFile;
  if (TransformationMatrixFile.Parse(Manager.GetString("transformation-file")) == false)
    {
      TransformationMatrixFile.DumpErrors(cout) << endl;
      return 0;
    }
  if (TransformationMatrixFile.GetNbrLines() == 0)
    {
      cout << Manager.GetString("transformation-file") << " is an empty file" << endl;
      return 0;
    }
  if (TransformationMatrixFile.GetNbrColumns() < 6)
    {
      cout << Manager.GetString("transformation-file") << " has a wrong number of columns (has "
	   << TransformationMatrixFile.GetNbrColumns() << ", should be at least 13)" << endl;
      return 0;
    }
  int NbrTransformationMatrices = TransformationMatrixFile.GetNbrLines();
  if (NbrTransformationMatrices != (NbrSitesX * NbrSitesY))
    {
       cout << Manager.GetString("transformation-file") << " has a wrong number of lins (has "
	    << TransformationMatrixFile.GetNbrColumns() << ", should be " << (NbrSitesX * NbrSitesY) << ")" << endl;
       return 0;
    }

  int* TmpKx = TransformationMatrixFile.GetAsIntegerArray(0);
  int* TmpKy = TransformationMatrixFile.GetAsIntegerArray(1);
  int* TmpLinearizedK = new int[NbrTransformationMatrices];
  Complex* TmpMatrix00 = TransformationMatrixFile.GetAsComplexArray(2);
  Complex* TmpMatrix10 = TransformationMatrixFile.GetAsComplexArray(3);
  Complex* TmpMatrix01 = TransformationMatrixFile.GetAsComplexArray(4);
  Complex* TmpMatrix11 = TransformationMatrixFile.GetAsComplexArray(5);
  ComplexMatrix* TransformationMatrices = new ComplexMatrix[NbrTransformationMatrices];
  for (int i = 0; i < NbrTransformationMatrices; ++i)
    {
      TmpLinearizedK[i] = TightBindingModel->GetLinearizedMomentumIndex(TmpKx[i], TmpKy[i]);
      TransformationMatrices[TmpLinearizedK[i]] = ComplexMatrix(2, 2, true);
      TransformationMatrices[TmpLinearizedK[i]].SetMatrixElement(0 , 0, TmpMatrix00[TmpLinearizedK[i]]);
      TransformationMatrices[TmpLinearizedK[i]].SetMatrixElement(0 , 1, TmpMatrix01[TmpLinearizedK[i]]);
      TransformationMatrices[TmpLinearizedK[i]].SetMatrixElement(1 , 0, TmpMatrix10[TmpLinearizedK[i]]);
      TransformationMatrices[TmpLinearizedK[i]].SetMatrixElement(1 , 1, TmpMatrix11[TmpLinearizedK[i]]);
    }

	 
  ParticleOnSphereWithSpin* Space = 0;
  if (Statistics == true)
    {
      if ((NbrSitesX * NbrSitesY) <= 32)
	{
	  Space = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
	}
      else
	{
	  Space = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
	}
    }
  else
    {
      Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, TotalKx, TotalKy);
    }
  
  
  ComplexVector FinalState (Space->GetHilbertSpaceDimension(), true);
  Space->TransformOneBodyBasis(InitialState, FinalState, TransformationMatrices, 0l, Space->GetLargeHilbertSpaceDimension());
  FinalState.WriteVector(OutputFileName);

  cout << "Norm of the rotated many-body state = " << FinalState.Norm() << endl;
  cout << "Overlap between the orginal and the rotated states = "<< (InitialState * FinalState) << endl;
  
  delete Space;
  delete TightBindingModel;
  
  return 0;
}
