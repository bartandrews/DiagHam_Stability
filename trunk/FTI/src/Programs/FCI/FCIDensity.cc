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
bool FCIDensityGetHilbertSpace(char* inputState, int& nbrParticles, int& nbrSitesX, int& nbrSitesY,
			       int& kxMomentum, int& kyMomentum, bool& statistics,
			       ParticleOnSphere*& space, ComplexVector& state, ParticleOnSphere** spaces, bool recomputeHilbert);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FCIDensity" , "0.01");
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

  (*SystemGroup) += new SingleStringOption  ('\n', "input-states", "use a file to describe the state as a linear combination");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (overriding the one found in the vector file name if greater than 0)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a file");
  (*SystemGroup) += new BooleanOption  ('\n', "coefficients-only", "only compute the one body coefficients that are requested to evaluate the density profile", false);
  (*SystemGroup) += new BooleanOption('\n', "recompute-hilbert", "do not store the Hilbert spaces (decreasing the memory consumption)");

  (*PrecalculationGroup) += new SingleStringOption  ('\n', "use-precomputed", "use precomputed matrix elements to do the plot");
  (*PlotOptionGroup) += new SingleStringOption ('\n', "output", "output file name (default output name replace the .vec extension of the input file with .rho.dat)", 0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCIDensity -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("input-states") == 0)
    {
      cout << "FCIDensity requires an input state" << endl;
      return -1;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey"); 
  bool Statistics = false;
  bool CoefficientOnlyFlag = Manager.GetBoolean("coefficients-only");
  char* OutputName = Manager.GetString("output");
  int BandIndex = 0;

  if ((CoefficientOnlyFlag == false) && (Manager.GetString("import-onebody") == 0))
    {
      cout << "error, a tight binding model as to be provided"  << endl;
      return -1;
    }

  MultiColumnASCIIFile InputVectors;
  if (InputVectors.Parse(Manager.GetString("input-states")) == false)
    {
      InputVectors.DumpErrors(cout) << endl;
      return -1;
    }

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
  Complex* Coefficients = 0;
  if (InputVectors(1, 0) != 0)
    {
      Coefficients = InputVectors.GetAsComplexArray(1);
    }
  else
    {
      if (CoefficientOnlyFlag == true)
	{
	  Coefficients = new Complex [InputVectors.GetNbrLines()];
	  for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
	    {
	      Coefficients[i] = 1.0;
	    }
	}
      else
	{
	  cout << "no coefficients defined in " << Manager.GetString("input-states") << endl;
	  return -1; 
	}
    }
  
  int ForceMaxMomentum = NbrSitesX * NbrSitesY - 1;      
  ComplexMatrix PrecalculatedValues(ForceMaxMomentum + 1, ForceMaxMomentum + 1, true);
  ComplexMatrix** RawPrecalculatedValues = new ComplexMatrix*[InputVectors.GetNbrLines()];
  for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
    {
      RawPrecalculatedValues[i] = new ComplexMatrix[InputVectors.GetNbrLines()];
      for (int j = 0; j < InputVectors.GetNbrLines(); ++j)
	RawPrecalculatedValues[i][j] = ComplexMatrix(ForceMaxMomentum + 1, ForceMaxMomentum + 1, true);
    }
  
  Generic2DTightBindingModel TightBindingModel(Manager.GetString("import-onebody")); 
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
	  if (FCIDensityGetHilbertSpace(InputVectors(0, i), NbrParticles, NbrSitesX, NbrSitesY, 
					LeftKxMomentum, LeftKyMomentum, Statistics, 
					LeftSpace, LeftState, Spaces, Manager.GetBoolean("recompute-hilbert")) == false)
	    return -1;
	  for (int j = 0; j < InputVectors.GetNbrLines(); ++j)
	    {
	      ParticleOnSphere* RightSpace = 0;
	      ComplexVector RightState;
	      int RightKxMomentum = 0;
	      int RightKyMomentum = 0;
	      if (FCIDensityGetHilbertSpace(InputVectors(0, j), NbrParticles, NbrSitesX, NbrSitesY, 
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
      char* TmpOutputName = ReplaceExtensionToFileName(Manager.GetString("input-states"), "dat", "coefficients.dat");
      if (TmpOutputName == 0)
	TmpOutputName = ReplaceExtensionToFileName(Manager.GetString("input-states"), "txt", "coefficients.dat");
      ofstream File;
      File.precision(14);
      File.open(TmpOutputName, ios::binary | ios::out);
      File << "# density coefficients for " << Manager.GetString("input-states") << endl;
      File << "#" << endl << "# state_index_left state_index_right m  n  <left|c^+_m c_n|right> " << endl;
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

  for (int m = 0; m <= ForceMaxMomentum; ++m)
    {
      for (int n = 0; n <= ForceMaxMomentum; ++n)
	{
	  for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
	    {
	      for (int j = 0; j < InputVectors.GetNbrLines(); ++j)
		{
		  Complex Tmp = RawPrecalculatedValues[i][j][n][m];
		  Tmp *= (Conj(Coefficients[i]) * Coefficients[j]);
		  PrecalculatedValues.AddToMatrixElement(m, n, Tmp);
		}
	    }
	}
    }

  int ImpurityXPositions = 1;
  int ImpurityYPositions = 0;
  int ImpurityOrbitals = 1;
  Complex TmpElement = 0.0;
  double Normalization = 1.0 / (((double) NbrSitesX) * ((double) NbrSitesY));      
  for (int m = 0; m < InputVectors.GetNbrLines(); ++m)
    for (int n = 0; n < InputVectors.GetNbrLines(); ++n)
      {
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
		cout << TotalKx1 << " " << TotalKx2 << " " << NbrSitesX << " | " << TotalKy1 << " " << TotalKy2 << " " << NbrSitesY << endl;
		TmpElement += (Phase(-2.0 * M_PI * ((double) ((TotalKx1 - TotalKx2) * ImpurityXPositions)) / ((double) NbrSitesX) 
				     - 2.0 * M_PI * ((double) ((TotalKy1 - TotalKy2) * ImpurityYPositions)) / ((double) NbrSitesY)) 
			       * Conj(TightBindingModel.GetOneBodyMatrix(i)[BandIndex][ImpurityOrbitals]) * (TightBindingModel.GetOneBodyMatrix(j)[BandIndex][ImpurityOrbitals])
			       * Normalization * RawPrecalculatedValues[m][n][i][j]) * (Conj(Coefficients[m]) * Coefficients[n]);
	      }
	  }
      }
  cout << "energy = " << TmpElement << endl;



  if (OutputName == 0)
    OutputName = ReplaceExtensionToFileName(Manager.GetString("input-states"), "dat", "rho.dat");
  if (OutputName == 0)
    OutputName = ReplaceExtensionToFileName(Manager.GetString("input-states"), "txt", "rho.dat");
  ofstream File;
  File.precision(14);
  File.open(OutputName, ios::binary | ios::out);
  File << "# density coefficients for " << Manager.GetString("input-states") << endl;
  File << "#" << endl << "# m  n  <psi| c^+_m c_n |psi>" << endl;
  for (int i = 0; i <= ForceMaxMomentum; ++i)
    for (int j = 0; j <= ForceMaxMomentum; ++j)
      File << "# " << i << " " << j << " " << PrecalculatedValues[i][j] << endl;
  File << "#" << endl;
  File << "# unit_cell_x unit_cell_y orbital_index density " << endl;
  Complex Sum = 0.0;
  //  double Normalization = 1.0 / (((double) NbrSitesX) * ((double) NbrSitesY));      
  for (int i = 0; i < NbrSitesX; ++i)
    {
      for (int j = 0; j < NbrSitesY; ++j)
	{
	  for (int alpha = 0; alpha < TightBindingModel.GetNbrBands(); ++alpha)
	    {
	      File << i << " " << j <<  " " << alpha;
	      Complex Density = 0.0;
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
		      Density += Phase(-2.0 * M_PI * ((double) ((TotalKx1 - TotalKx2) * i)) / ((double) NbrSitesX) 
				       - 2.0 * M_PI * ((double) ((TotalKy1 - TotalKy2) * j)) / ((double) NbrSitesY)) * Normalization * Conj(TightBindingModel.GetOneBodyMatrix(m)[BandIndex][alpha]) * (TightBindingModel.GetOneBodyMatrix(n)[BandIndex][alpha]) * PrecalculatedValues[n][m];
		    }
		}
	      File << " " << Density.Re << " " << Density << endl;
	      Sum += Density;
	    }
	  File << endl;
	}
      File << endl;
    }
  cout << "Sum = " << Sum << endl;

  File.close();
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

bool FCIDensityGetHilbertSpace(char* inputState, int& nbrParticles, int& nbrSitesX, int& nbrSitesY,
			       int& kxMomentum, int& kyMomentum, bool& statistics,
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
