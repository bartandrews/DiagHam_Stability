#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"


#include "Operator/ParticleOnSphereDensityOperator.h"
#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "FunctionBasis/ParticleOnCylinderFunctionBasis.h"
#include "Hamiltonian/ParticleOnCylinderStructureFactor.h"
#include "Hamiltonian/ParticleOnCylinderPairAmplitude.h"
#include "Hamiltonian/ParticleOnCylinderDensityDensity.h"

#include "Tools/FQHEFiles/FQHEOnCylinderFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "Options/OptionManager.h"
#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHECylinderDensity" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('\n', "input-state", "name of the file containing the state to convert");
  (*SystemGroup) += new SingleStringOption('\n', "degenerate-states", "name of the file containing a list of states (override input-state)");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state of the output state");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-points", "number of points along the cylinder axis", 100);
  (*OutputGroup) += new SingleDoubleOption  ('\n', "offset", "additional length along the cylinder axis on each side of the [-Lx/2,Lx/2] region where the density should be computed", 5.0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with rhorho extension");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECylinderDensity -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int KyMax = 0;
  int TotalKy = 0;
  bool Statistics = true;  
  int NbrPoints = Manager.GetInteger("nbr-points");
  double Ratio = 0.0;
  double Perimeter = 0.0;
  int LandauLevel = 0;
  int NbrInputStates = 0;
  char** InputStateNames = 0;

  if ((Manager.GetString("input-state") == 0) && (Manager.GetString("degenerate-states") == 0) )
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type FQHECylinderDensity -h" << endl;
      return -1;
    }

  if (Manager.GetString("input-state") != 0)
    {
      if (FQHEOnCylinderFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, KyMax, TotalKy,Statistics, Ratio, Perimeter) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state")  << endl;
	  return -1;
	}
      NbrInputStates = 1;
      InputStateNames = new char*[NbrInputStates];
      InputStateNames[0] = new char [strlen(Manager.GetString("input-state")) + 1];
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
      NbrInputStates = DegenerateFile.GetNbrLines();
      InputStateNames = new char*[NbrInputStates];
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  InputStateNames[i] = new char [strlen(DegenerateFile(0, i)) + 1];
	  strcpy (InputStateNames[i], DegenerateFile(0, i));
	}
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  if (FQHEOnCylinderFindSystemInfoFromVectorFileName(InputStateNames[i], NbrParticles, KyMax, TotalKy, Statistics, Ratio, Perimeter) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << InputStateNames[i] << endl;
	      return -1;
	    }
	}
    }

  if ((Perimeter == 0.0) & (Ratio == 0.0))
    {
      cout << "error, neither the aspect ratio or the cylinder perimeter is defined" << endl;
      return -1;
    }
  if (Perimeter == 0.0)
    {
      Perimeter = sqrt(2.0 * M_PI * (KyMax + 1) * Ratio);
    }
  else
    {
      Ratio = (Perimeter * Perimeter) / (2.0 * M_PI * (KyMax + 1));
    }
  double CylinderLength = Perimeter / Ratio;

  ParticleOnSphere* Space = 0;
  if (Statistics == true)
    {
#ifdef __64_BITS__
      if (TotalKy <= 62)
#else
	if (TotalKy <= 30)
#endif
	  {
	    if (Manager.GetBoolean("haldane") == true)
	      {
		int* ReferenceState = 0;
		if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, KyMax, ReferenceState) == false)
		  return -1;
		if (Manager.GetString("load-hilbert") != 0)
		  Space = new FermionOnSphereHaldaneBasis(Manager.GetString("load-hilbert"));	  
		else
		  Space = new FermionOnSphereHaldaneBasis(NbrParticles, TotalKy, KyMax, ReferenceState);
	      }
	    else
	      {
		Space = new FermionOnSphere(NbrParticles, TotalKy, KyMax);
	      }
	  }      
	else
#ifdef __128_BIT_LONGLONG__
	  if (TotalKy <= 126)
#else
	    if (TotalKy <= 62)
#endif
	      Space = new FermionOnSphereLong(NbrParticles, TotalKy, KyMax);
	    else
	      Space = new FermionOnSphereUnlimited(NbrParticles, TotalKy, KyMax);
    }

  cout << " Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;

  
  RealVector* InputStates = new RealVector[NbrInputStates];
  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (InputStates[i].ReadVector (InputStateNames[i]) == false)
	{
	  cout << "can't open vector file " << InputStateNames[i] << endl;
	  return -1;      
		}	  
      if (Space->GetHilbertSpaceDimension() != InputStates[i].GetVectorDimension())
	{
	  cout << "error, " << InputStateNames[i] << " does not have the correct dimension (is "
	       << InputStates[i].GetVectorDimension() << ", should be " << Space->GetHilbertSpaceDimension() << ")" << endl;
	  return -1;
	}
    }

  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = ReplaceExtensionToFileName(InputStateNames[0], "vec", "rho.dat");
    }


  ofstream File;
  File.precision(14);
  File.open(OutputFileName, ios::binary | ios::out);
  File << "# Length = " << CylinderLength << ", Perimeter = " << Perimeter  << ", Aspect ratio = " << Ratio << endl;
  HermitianMatrix* OccupationMatrixElements = new HermitianMatrix[KyMax + 1];
  
      
  for (int i = 0; i <= KyMax; ++i)
    {
      OccupationMatrixElements[i] = HermitianMatrix(NbrInputStates, true);
      ParticleOnSphereDensityOperator TmpOperator (Space, i);
      RealVector TmpState(Space->GetHilbertSpaceDimension(), true);
      for (int j = 0; j < NbrInputStates; ++j)
	{
	  VectorOperatorMultiplyOperation Operation (&TmpOperator, &InputStates[j], &TmpState);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  for (int k = j; k <  NbrInputStates; ++k)
	    {
	      Complex Tmp = InputStates[k] * TmpState;
	      OccupationMatrixElements[i].SetMatrixElement(k, j, Tmp);
	      File << "# " << i << " " << j << " " << k << " " << Tmp << endl;
	    }
	}
    }

  ParticleOnCylinderFunctionBasis Basis (KyMax, LandauLevel, Ratio);
  HermitianMatrix TmpMatrix (NbrInputStates, true);
  HermitianMatrix TmpMatrix2 (NbrInputStates, true);
  RealDiagonalMatrix TmpEigenvalues(NbrInputStates);
  RealDiagonalMatrix TmpEigenvalues2(NbrInputStates);
  double Offset = Manager.GetDouble("offset");
  double XPosition = -(0.5 * CylinderLength) - Offset;
  double XStep = -2.0 * XPosition / ((double) (NbrPoints + 1));
  double TmpPrefactor1 =  1.0 / sqrt(M_PI);
  double TmpPrefactor2 =  2.0 * M_PI / Perimeter;
  double TmpShift = 0.5 * ((double) KyMax) * TmpPrefactor2;
  double* TotalCharges = new double[NbrInputStates];
  for (int i = 0; i < NbrInputStates; ++i)
    {   
      TotalCharges[i] = 0.0;
    }
  for (int TmpX = 0; TmpX <= NbrPoints; ++TmpX)
    { 
      TmpMatrix.ClearMatrix();
      TmpMatrix2.ClearMatrix();
      for (int i = 0; i <= KyMax; ++i)
	{
	  double TmpFactor = TmpPrefactor1 * exp(-(XPosition + TmpShift - (TmpPrefactor2 * ((double) i))) * (XPosition + TmpShift - (TmpPrefactor2 * ((double) i))));
	  TmpMatrix.AddLinearCombination(TmpFactor, OccupationMatrixElements[i]);
	  TmpFactor = 0.5 * (1.0 + erf(XPosition + TmpShift - (TmpPrefactor2 * ((double) i))));
	  TmpMatrix2.AddLinearCombination(TmpFactor, OccupationMatrixElements[i]);
	}
      TmpMatrix.LapackDiagonalize(TmpEigenvalues);
      TmpMatrix2.LapackDiagonalize(TmpEigenvalues2);
      File << XPosition;
      for (int i = 0; i < NbrInputStates; ++i)
	{   
	  File << " " << TmpEigenvalues[i] << " " << TmpEigenvalues2[i];
	  TotalCharges[i] +=  TmpEigenvalues[i] * XStep;
	}      
      File << endl;
      XPosition += XStep;
    }

  cout << "total charge : " << endl;
  for (int i = 0; i < NbrInputStates; ++i)
    {   
      cout << TotalCharges[i] << " ";
    }
  cout << endl;
  File.close();

  return 0;
}

