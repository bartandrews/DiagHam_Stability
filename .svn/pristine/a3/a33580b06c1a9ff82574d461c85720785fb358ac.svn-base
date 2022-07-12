#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"

#include "FunctionBasis/ParticleOnChernInsulatorSingleBandFunctionBasis.h"
#include "FunctionBasis/ParticleOnCheckerboardLatticeFunctionBasis.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

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

int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FCICorrelation" , "0.01");
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

  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the vector file describing the state whose density has to be plotted");
  (*SystemGroup) += new BooleanOption  ('\n', "density", "plot density instead of density-density correlation", false);
  (*SystemGroup) += new BooleanOption  ('\n', "k-space", "compute the density/correlation in momentum space", false);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", 1.0 - 0.5 * M_SQRT2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tpp", "second next nearest neighbor hoping amplitude", 0.5 * (M_SQRT2 - 1.0));
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);

  (*PlotOptionGroup) += new SingleStringOption ('\n', "output", "output file name (default output name replace the .vec extension of the input file with .rho or .rhorho)", 0);
  (*PlotOptionGroup) += new SingleIntegerOption ('\n', "nbr-samplesx", "number of samples along the x direction", 100, true, 10);
  (*PlotOptionGroup) += new SingleIntegerOption ('\n', "nbr-samplesy", "number of samples along the y direction", 100, true, 10);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCICorrelation -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("state") == 0)
    {
      cout << "FCICorrelation requires an input state" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state")) == false)
    {
      cout << "can't find vector file " << Manager.GetString("state") << endl;
      return -1;      
    }

  int NbrParticles = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  int MomentumX = 0;
  int MomentumY = 0;
  double Mass = 0.0;
  int NbrSamplesX = Manager.GetInteger("nbr-samplesx");
  int NbrSamplesY = Manager.GetInteger("nbr-samplesy");
  bool DensityFlag = Manager.GetBoolean("density");
  bool Statistics = true;

  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(Manager.GetString("state"),
							  NbrParticles, NbrSiteX, NbrSiteY, MomentumX, MomentumY, Mass, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state") << endl;
      return -1;
    }
  cout << "N=" << NbrParticles << " Nx=" << NbrSiteX << " Ny=" << NbrSiteY << " kx=" << MomentumX << " ky=" << MomentumY << " m=" << Mass << endl; 
  ParticleOnSphere* Space = 0;
  if (Statistics == true)
    Space = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, MomentumX, MomentumY);
  else
    Space = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, MomentumX, MomentumY);
  ComplexVector ComplexState;
  if (ComplexState.ReadVector (Manager.GetString("state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("state") << endl;
      return -1;      
    }
  Complex* PrecalculatedValues = 0;
  int* PrecalculatedIndices = 0;
  int NbrPrecalculatedValues = 0;
  if (DensityFlag == false)
    {
      for (int kx1 =0; kx1 < NbrSiteX; ++kx1)
	for (int kx2 =0; kx2 < NbrSiteX; ++kx2)
	  for (int kx3 =0; kx3 < NbrSiteX; ++kx3)
	    for (int kx4 =0; kx4 < NbrSiteX; ++kx4)
	      {
		if (((kx1 + kx2 - kx3 - kx4) % NbrSiteX) == 0)
		  {
		    for (int ky1 = 0; ky1 < NbrSiteY; ++ky1)
		      for (int ky2 = 0; ky2 < NbrSiteY; ++ky2)
			for (int ky3 = 0; ky3 < NbrSiteY; ++ky3)
			  for (int ky4 = 0; ky4 < NbrSiteY; ++ky4)
			    {
			      if (((ky1 + ky2 - ky3 - ky4) % NbrSiteY) == 0)
				{
				  ++NbrPrecalculatedValues;
				}
			    }
			}
	      }
      PrecalculatedValues = new Complex [NbrPrecalculatedValues];
      PrecalculatedIndices = new int [4 * NbrPrecalculatedValues];
      NbrPrecalculatedValues = 0; 
      for (int kx1 =0; kx1 < NbrSiteX; ++kx1)
	for (int kx2 =0; kx2 < NbrSiteX; ++kx2)
	  for (int kx3 =0; kx3 < NbrSiteX; ++kx3)
	    for (int kx4 =0; kx4 < NbrSiteX; ++kx4)
	      {
		if (((kx1 + kx2 - kx3 - kx4) % NbrSiteX) == 0)
		  {
		    for (int ky1 = 0; ky1 < NbrSiteY; ++ky1)
		      for (int ky2 = 0; ky2 < NbrSiteY; ++ky2)
			for (int ky3 = 0; ky3 < NbrSiteY; ++ky3)
			  for (int ky4 = 0; ky4 < NbrSiteY; ++ky4)
			    {
			      if (((ky1 + ky2 - ky3 - ky4) % NbrSiteY) == 0)
				{
				  int Index1 = (kx1 * NbrSiteY) + ky1;
				  int Index2 = (kx2 * NbrSiteY) + ky2;
				  int Index3 = (kx3 * NbrSiteY) + ky3;
				  int Index4 = (kx4 * NbrSiteY) + ky4;
				  ParticleOnSphereDensityDensityOperator Operator (Space, Index1, Index2, Index3, Index4);
				  PrecalculatedValues[NbrPrecalculatedValues] = Operator.MatrixElement(ComplexState, ComplexState);
				  PrecalculatedIndices[(NbrPrecalculatedValues << 2)] = Index1;
				  PrecalculatedIndices[(NbrPrecalculatedValues << 2) + 1] = Index2;
				  PrecalculatedIndices[(NbrPrecalculatedValues << 2) + 2] = Index3;
				  PrecalculatedIndices[(NbrPrecalculatedValues << 2) + 3] = Index4;
				  ++NbrPrecalculatedValues;
				}
			    }
		  }
	      }
    }
  else
    {
      NbrPrecalculatedValues = NbrSiteX * NbrSiteY;
      PrecalculatedValues = new Complex [NbrPrecalculatedValues];
      for (int kx =0; kx < NbrSiteX; ++kx)
	for (int ky = 0; ky < NbrSiteY; ++ky)
	  {
	    int Index = (kx * NbrSiteY) + ky;
	    ParticleOnSphereDensityOperator Operator (Space, Index);	    
	    PrecalculatedValues[Index] = Operator.MatrixElement(ComplexState, ComplexState);
	  }
    }  
  delete Space;
  ofstream File;
  File.precision(14);
  double XStep = ((double) NbrSiteX) / ((double) NbrSamplesX);
  double YStep = ((double) NbrSiteY) / ((double) NbrSamplesY);
  RealVector Position(2, true);
  if (Manager.GetString("output") != 0)
    File.open(Manager.GetString("output"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = 0;
      if (DensityFlag == false)
	{
	  TmpFileName = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", "rhorho");
	}
      else
	{
	  TmpFileName = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", "rho");
	}
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << Manager.GetString("state") << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }

  if (Manager.GetBoolean("k-space") == true)
    {
      if (DensityFlag == true)
	{
	  File << "# kx ky n(kx,ky)" << endl;
	  for (int kx =0; kx < NbrSiteX; ++kx)
	    for (int ky = 0; ky < NbrSiteY; ++ky)
	      {
		int Index = (kx * NbrSiteY) + ky;
		File << kx << " " << ky << " " << PrecalculatedValues[Index].Re << endl;;
	      }
	  File.close();
	}
      return 0;
    }
  
  ParticleOnChernInsulatorSingleBandFunctionBasis Basis(NbrSiteX, NbrSiteY, Mass);
  int TotalNbrSites = NbrSiteX * NbrSiteY;
  Complex* Coefficients = new Complex[TotalNbrSites];
  Complex* Coefficients2 = new Complex[TotalNbrSites];
  for (int i = 0; i < TotalNbrSites; ++i)
    {
      Basis.GetFunctionValue(Position, Coefficients[i], i);      
    }
  for (int i = 0; i < NbrSamplesX; ++i)
    {
      Position[1] = 0.0;
      for (int j = 0; j < NbrSamplesY; ++j)
	{
	  Complex TmpValue = 0.0;
	  if (DensityFlag == false)
	    {
	      for (int i = 0; i < TotalNbrSites; ++i)
		{
		  Basis.GetFunctionValue(Position, Coefficients2[i], i);      
		}
	      for (int i = 0; i < NbrPrecalculatedValues; ++i) 
		{
		  TmpValue -= (PrecalculatedValues[i] * Conj(Coefficients[PrecalculatedIndices[(i << 2)]]) *
			       Coefficients[PrecalculatedIndices[(i << 2) + 2]] * Conj(Coefficients2[PrecalculatedIndices[(i << 2) + 1]]) 
			       * Coefficients2[PrecalculatedIndices[(i << 2) + 3]]);
		}
	    }
	  else
	    {	      
	      for (int i = 0; i < NbrPrecalculatedValues; ++i) 
		{
		  Complex TmpValue2;
		  Basis.GetFunctionValue(Position, TmpValue2, i);
		  TmpValue += PrecalculatedValues[i] * SqrNorm(TmpValue2);
		}
	    }
	  File << Position[0] << " " << Position[1] << " " << TmpValue.Re << endl;
	  Position[1] += YStep;
	}
      Position[0] += XStep;
      File << endl;
    }
  File.close();
  delete[] Coefficients;
  delete[] Coefficients2;
  return 0;
}
