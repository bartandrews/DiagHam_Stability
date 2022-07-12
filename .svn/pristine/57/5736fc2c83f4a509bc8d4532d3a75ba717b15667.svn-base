#include "Vector/RealVector.h"
#include "Matrix/ComplexMatrix.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Operator/ParticleOnSphereDensityOperator.h"
#include "FunctionBasis/ParticleOnDiskFunctionBasis.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


// get the Hilbert space and the vector state form the input file name
//
// inputState = input file name
// nbrParticles = reference on the number of particles 
// nbrSitesX = reference on the number of site along the x direction
// nbrSitesY = reference on the number of site along the y direction
// kxMomentum = reference on the momentum along the x direction
// kyMomentum = reference on the momentum along the y direction
// statistics = reference on the statistic flag
// space = reference on the pointer to the Hilbert space
// state = reference on the state vector
// spaces = array where the Hilbert spaces are stored
// recomputeHilbert = true if the Hilbert space should not be stored
// return value = true if no error occured
bool FCIImpuritiesGetHilbertSpace(char* inputState, int& nbrParticles, int& nbrSitesX, int& nbrSitesY,
				  int& kxMomentum, int& kyMomentum, bool& statistics,
				  ParticleOnSphere*& space, ComplexVector& state,
				  ParticleOnSphere** spaces, bool recomputeHilbert);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FCIImpurities" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PlotOptionGroup = new OptionGroup ("plot options");  
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += PlotOptionGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\n', "input-states", "ASCII column formatted file that describes the low energy states of the system without impurities");
  (*SystemGroup) += new SingleStringOption  ('\n', "impurities", "ASCII column formatted file that gives the location and strength of each impurity");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a file");
  (*SystemGroup) += new BooleanOption('\n', "recompute-hilbert", "do not store the Hilbert spaces (decreasing the memory consumption)");

  (*PrecalculationGroup) += new SingleStringOption  ('\n', "use-precomputed", "use precomputed matrix elements to perform the calculation");
  (*PlotOptionGroup) += new SingleStringOption ('\n', "output", "output file name (default output name replace the .vec extension of the input file with .impurities.dat)", 0);
  (*PlotOptionGroup) += new BooleanOption ('\n', "binary-output", "export the eignestates as binary vectors instead of text files");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCIImpurities -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("input-states") == 0)
    {
      cout << "FCIImpurities requires an input state" << endl;
      return -1;
    }

  int NbrParticles = 0;
  int NbrSitesX = 0; 
  int NbrSitesY = 0; 
  bool Statistics = false;
  char* OutputName = Manager.GetString("output");
  int BandIndex = 0;

  if (Manager.GetString("import-onebody") == 0)
    {
      cout << "error, a tight binding model has to be provided"  << endl;
      return -1;
    }

  MultiColumnASCIIFile InputVectors;
  if (InputVectors.Parse(Manager.GetString("input-states")) == false)
    {
      InputVectors.DumpErrors(cout) << endl;
      return -1;
    }

  if (Manager.GetString("impurities") == 0)
    {
      cout << "error, impurity locations have to be provided"  << endl;
      return -1;
    }
  MultiColumnASCIIFile Impurities;
  if (Impurities.Parse(Manager.GetString("impurities")) == false)
    {
      Impurities.DumpErrors(cout) << endl;
      return -1;
    }

  int NbrImpurities = Impurities.GetNbrLines();
  if (NbrImpurities == 0)
    {
      cout << "no impurity location provided in " << Manager.GetString("impurities") << endl;
      return -1;
    }
  if (Impurities.GetNbrColumns() < 3)
    {
      cout << "wrong number of columns in " << Manager.GetString("impurities") << endl;
      return -1;
    }
  int* ImpurityXPositions = Impurities.GetAsIntegerArray(0);
  int* ImpurityYPositions = Impurities.GetAsIntegerArray(1);
  int* ImpurityOrbitals = Impurities.GetAsIntegerArray(2);
  double* ImpurityStrengths = 0;
  if (Impurities.GetNbrColumns() > 3)
    {
      ImpurityStrengths = Impurities.GetAsDoubleArray(3);;
    }
  else
    {
      ImpurityStrengths = new double [NbrImpurities];
      for (int i = 0; i < NbrImpurities; ++i)
	ImpurityStrengths[i] = 1.0;
    }
  Generic2DTightBindingModel TightBindingModel(Manager.GetString("import-onebody")); 
  int CurrentKx = 0;
  int CurrentKy = 0;
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(InputVectors(0, 0), NbrParticles, NbrSitesX, NbrSitesY, CurrentKx, CurrentKy, Statistics) == false)
    {
      return -1;      
    }
  for (int i = 1; i < InputVectors.GetNbrLines(); ++i)
    {
      int TmpNbrParticles = 0;
      int TmpNbrSitesX = 0;
      int TmpNbrSitesY = 0;
      int TmpCurrentKx = 0;
      int TmpCurrentKy = 0;
      if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(InputVectors(0, i), TmpNbrParticles, TmpNbrSitesX, TmpNbrSitesY, TmpCurrentKx, TmpCurrentKy, Statistics) == false)
	{
	  return -1;      
	}
      if (NbrParticles != TmpNbrParticles)
	{
	  cout << InputVectors(0, i) << " and " << InputVectors(0, 0) << " have different numbers of particles" << endl;
	  return -1;
	}
      if (TmpNbrSitesX != NbrSitesX)
	{
	  cout << InputVectors(0, i) << " and " << InputVectors(0, 0) << " have different number of sites in the x direction" << endl;
	  return -1;
	}
      if (TmpNbrSitesY != NbrSitesY)
	{
	  cout << InputVectors(0, i) << " and " << InputVectors(0, 0) << " have different number of sites in the y direction" << endl;
	  return -1;
	}
    }
  
  int ForceMaxMomentum = NbrSitesX * NbrSitesY - 1;      
  ComplexMatrix** RawPrecalculatedValues = new ComplexMatrix*[InputVectors.GetNbrLines()];
  for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
    {
      RawPrecalculatedValues[i] = new ComplexMatrix[InputVectors.GetNbrLines()];
      for (int j = 0; j < InputVectors.GetNbrLines(); ++j)
	RawPrecalculatedValues[i][j] = ComplexMatrix(ForceMaxMomentum + 1, ForceMaxMomentum + 1, true);
    }
  
  if (Manager.GetString("use-precomputed") == 0)
    {	  
      ParticleOnSphere** Spaces = new ParticleOnSphere* [NbrSitesX * NbrSitesY];
      for (int m = 0; m <= ForceMaxMomentum; ++m)
	Spaces[m] = 0;
      for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
	{
	  ParticleOnSphere* LeftSpace = 0;
	  ComplexVector LeftState;
	  int LeftKxMomentum = 0;
	  int LeftKyMomentum = 0;
	  if (FCIImpuritiesGetHilbertSpace(InputVectors(0, i), NbrParticles, NbrSitesX, NbrSitesY, 
					   LeftKxMomentum, LeftKyMomentum, Statistics, 
					   LeftSpace, LeftState, Spaces, Manager.GetBoolean("recompute-hilbert")) == false)
	    return -1;
	  for (int j = 0; j < InputVectors.GetNbrLines(); ++j)
	    {
	      ParticleOnSphere* RightSpace = 0;
	      ComplexVector RightState;
	      int RightKxMomentum = 0;
	      int RightKyMomentum = 0;
	      if (FCIImpuritiesGetHilbertSpace(InputVectors(0, j), NbrParticles, NbrSitesX, NbrSitesY, 
					       RightKxMomentum, RightKyMomentum, Statistics, 
					       RightSpace, RightState, Spaces, Manager.GetBoolean("recompute-hilbert")) == false)
		return -1;
	      RightSpace->SetTargetSpace(LeftSpace);
	      for (int m = 0; m <= ForceMaxMomentum; ++m)
		{
		  int TotalKx1;
		  int TotalKy1;
		  TightBindingModel.GetLinearizedMomentumIndex(m, TotalKx1, TotalKy1);
		  for (int n = 0; n <= ForceMaxMomentum; ++n)
		    {
		      int TotalKx2;
		      int TotalKy2;
		      TightBindingModel.GetLinearizedMomentumIndex(n, TotalKx2, TotalKy2);
		      int TmpRightKxMomentum = RightKxMomentum - TotalKx2;
		      if (TmpRightKxMomentum < 0)
			TmpRightKxMomentum += NbrSitesX;
		      int TmpLeftKxMomentum = LeftKxMomentum - TotalKx1;
		      if (TmpLeftKxMomentum < 0)
			TmpLeftKxMomentum += NbrSitesX;
		      int TmpRightKyMomentum = RightKyMomentum - TotalKy2;
		      if (TmpRightKyMomentum < 0)
			TmpRightKyMomentum += NbrSitesY;
		      int TmpLeftKyMomentum = LeftKyMomentum - TotalKy1;
		      if (TmpLeftKyMomentum < 0)
			TmpLeftKyMomentum += NbrSitesY;
		      if ((TmpRightKxMomentum == TmpLeftKxMomentum) && (TmpRightKyMomentum == TmpLeftKyMomentum))
			{
			  ParticleOnSphereDensityOperator Operator (RightSpace, m, n);
			  OperatorMatrixElementOperation Operation(&Operator, LeftState, RightState, RightState.GetLargeVectorDimension());
			  Operation.ApplyOperation(Architecture.GetArchitecture());
			  RawPrecalculatedValues[i][j].AddToMatrixElement(m, n, Operation.GetScalar());
			}
		    }
		}
	      if (Manager.GetBoolean("recompute-hilbert") == true)
		delete RightSpace;
	    }
	}
    }
  else
    {
      MultiColumnASCIIFile PrecomputedElements;
      if (PrecomputedElements.Parse(Manager.GetString("use-precomputed")) == false)
	{
	  PrecomputedElements.DumpErrors(cout) << endl;
	  return -1;
	}
      int* StateIndexLeft = PrecomputedElements.GetAsIntegerArray(0);
      int* StateIndexRight = PrecomputedElements.GetAsIntegerArray(1);
      int* OperatorIndexLeft = PrecomputedElements.GetAsIntegerArray(2);
      int* OperatorIndexRight = PrecomputedElements.GetAsIntegerArray(3);	  
      Complex* OneBodyCoefficients = PrecomputedElements.GetAsComplexArray(4);	  
      for (int i = 0; i < PrecomputedElements.GetNbrLines(); ++i)
	{
	  RawPrecalculatedValues[StateIndexLeft[i]][StateIndexRight[i]].AddToMatrixElement(OperatorIndexLeft[i], OperatorIndexRight[i], OneBodyCoefficients[i]);
	}
    }



  if (OutputName == 0)
    OutputName = ReplaceExtensionToFileName(Manager.GetString("input-states"), "dat", "coefficients.dat");
  if (Manager.GetString("use-precomputed") == 0)
    {	  
      ofstream File;
      File.precision(14);
      File.open(OutputName, ios::binary | ios::out);
      File << "# density coefficients for " << Manager.GetString("input-states") << endl;
      File << "#" << endl << "# state_index_left state_index_right m  n  <left|c^+_m c_n|right>" << endl;
      for (int m = 0; m < InputVectors.GetNbrLines(); ++m)
	for (int n = 0; n < InputVectors.GetNbrLines(); ++n)
	  for (int i = 0; i <= ForceMaxMomentum; ++i)
	    for (int j = 0; j <= ForceMaxMomentum; ++j)
	      {
		if ((RawPrecalculatedValues[m][n][j][i].Re != 0.0) || (RawPrecalculatedValues[m][n][j][i].Im != 0.0))
		  File << m << " " << n << " " << i << " " << j << " " << RawPrecalculatedValues[m][n][j][i] << endl;
	      }
      File.close();
    }

  double Normalization = 1.0 /  ((double) (NbrSitesX * NbrSitesY));  
  int TmpHilbertSpaceDimension = InputVectors.GetNbrLines();
  HermitianMatrix HRep (TmpHilbertSpaceDimension, true);
  for (int m = 0; m < TmpHilbertSpaceDimension; ++m)
    for (int n = m; n < TmpHilbertSpaceDimension; ++n)
      {
	Complex TmpElement = 0.0;
	for (int i = 0; i <= ForceMaxMomentum; ++i)
	  {
	    int TotalKx1;
	    int TotalKy1;
	    TightBindingModel.GetLinearizedMomentumIndex(i, TotalKx1, TotalKy1);
	    for (int j = 0; j <= ForceMaxMomentum; ++j)
	      {
		int TotalKx2;
		int TotalKy2;
		TightBindingModel.GetLinearizedMomentumIndex(j, TotalKx2, TotalKy2);
		for (int p = 0; p < NbrImpurities; ++p)
		  {
		    TmpElement += (Phase(-2.0 * M_PI * ((double) ((TotalKx1 - TotalKx2) * ImpurityXPositions[p])) / ((double) NbrSitesX) 
					 - 2.0 * M_PI * ((double) ((TotalKy1 - TotalKy2) * ImpurityYPositions[p])) / ((double) NbrSitesY)) 
				   * Conj(TightBindingModel.GetOneBodyMatrix(i)[BandIndex][ImpurityOrbitals[p]]) * (TightBindingModel.GetOneBodyMatrix(j)[BandIndex][ImpurityOrbitals[p]])
				   * Normalization * RawPrecalculatedValues[m][n][i][j]);
		  }
	      }
	  }
	HRep.SetMatrixElement(m, n, TmpElement);
      }


  
  RealDiagonalMatrix TmpDiag (TmpHilbertSpaceDimension);
  ComplexMatrix TmpEigenvector (TmpHilbertSpaceDimension, TmpHilbertSpaceDimension);	      
  if (TmpHilbertSpaceDimension > 1)
    {
      TmpEigenvector.SetToIdentity();
#ifdef __LAPACK__
      HRep.LapackDiagonalize(TmpDiag, TmpEigenvector);
#else
      HRep.Diagonalize(TmpDiag, TmpEigenvector);
#endif
    }
  else
    {
      TmpEigenvector[0][0] = 1.0;
      TmpDiag[0] = HRep(0, 0);
    }

  char* TruncatedOutputNameSpectrum = ReplaceExtensionToFileName(OutputName, ".dat", "");
  char* TmpOutputNamePrefix = new char [strlen(OutputName) + strlen(Manager.GetString("impurities")) + 128];
  sprintf (TmpOutputNamePrefix, "%s.%s", TruncatedOutputNameSpectrum, Manager.GetString("impurities"));
  char* OutputNamePrefix = ReplaceExtensionToFileName(TmpOutputNamePrefix, ".dat", "");
  char* OutputNameSpectrum = AddExtensionToFileName(OutputNamePrefix, "dat");
  ofstream File2;
  File2.precision(14);
  File2.open(OutputNameSpectrum, ios::binary | ios::out);
  for (int i = 0; i < TmpHilbertSpaceDimension; ++i)
    {
      File2 << TmpDiag[i] << endl;
      char* VectorSuffix = new char [128];
      if (Manager.GetBoolean("binary-output"))
	{
	  sprintf(VectorSuffix, "%d.vec", i);
	  char* OutputNameVector = AddExtensionToFileName(OutputNamePrefix, VectorSuffix);
	  if (TmpEigenvector[i].WriteVector(OutputNameVector) == false)
	    {
	      cout << "can't write vector " << OutputNameVector << endl;
	    }
	  delete[] OutputNameVector;
	}
      else
	{
	  sprintf (VectorSuffix, "%d.vec.txt", i);
	  char* OutputNameVector = AddExtensionToFileName(OutputNamePrefix, VectorSuffix);
	  ofstream File3;
	  File3.precision(14);
	  File3.open(OutputNameVector, ios::binary | ios::out);
	  File3 << "# vector coefficient" << endl;
	  for (int j = 0; j < InputVectors.GetNbrLines(); ++j)
	    {
	      File3 << InputVectors(0, j) << " " << Conj(TmpEigenvector[i][j]) << endl;
	    }
	  File3.close();	  
	  delete[] OutputNameVector;
	}
      delete[] VectorSuffix;
    }
  File2.close();

}

// get the Hilbert space and the vector state form the input file name
//
// inputState = input file name
// nbrParticles = reference on the number of particles 
// nbrSitesX = reference on the number of site along the x direction
// nbrSiteY = reference on the number of site along the y direction
// kxMomentum = reference on the momentum along the x direction
// kyMomentum = reference on the momentum along the y direction
// statistics = reference on the statistic flag
// space = reference on the pointer to the Hilbert space
// state = reference on the state vector
// spaces = array where the Hilbert spaces are stored
// recomputeHilbert = true if the Hilbert space should not be stored
// return value = true if no error occured

bool FCIImpuritiesGetHilbertSpace(char* inputState, int& nbrParticles, int& nbrSitesX, int& nbrSitesY, int& kxMomentum, int& kyMomentum, bool& statistics,
				  ParticleOnSphere*& space, ComplexVector& state, ParticleOnSphere** spaces, bool recomputeHilbert)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(inputState, nbrParticles, nbrSitesX, nbrSitesY, kxMomentum, kyMomentum, statistics) == false)
    {
      return false;      
    }
  if (state.ReadVector (inputState) == false)
    {
      cout << "can't open vector file " << inputState << endl;
      return false;      
    }

  if (spaces[kxMomentum + nbrSitesX * kyMomentum] == 0)
    {
      if (statistics == true)
	{
#ifdef __64_BITS__
	  if ((nbrSitesX * nbrSitesY) <= 63)
#else
	    if ((nbrSitesX * nbrSitesY) <= 31)
#endif
	      {
		space = new FermionOnSquareLatticeMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, kxMomentum, kyMomentum);
	      }
	    else
	      {
		space = new FermionOnSquareLatticeMomentumSpaceLong (nbrParticles, nbrSitesX, nbrSitesY, kxMomentum, kyMomentum);
	      }
	}
      else
	{
	  space = new BosonOnSquareLatticeMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, kxMomentum, kyMomentum);
	}
      if (recomputeHilbert == false)
	spaces[kxMomentum + nbrSitesX * kyMomentum] = space;
    }
  else
    {
      space = spaces[kxMomentum + nbrSitesX * kyMomentum];
    }
  if (space->GetLargeHilbertSpaceDimension() != state.GetLargeVectorDimension())
    {
      cout << "dimension mismatch between the state (" << state.GetLargeVectorDimension() << ") and the Hilbert space (" << space->GetLargeHilbertSpaceDimension() << ")" << endl;
      return false;
    }
  return true;
}
