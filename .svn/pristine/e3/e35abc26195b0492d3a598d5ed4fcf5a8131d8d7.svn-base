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
  OptionManager Manager ("FQHECheckerboardLatticeModelForgeState" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\n', "input-state", "file describing the state real space occupation");
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", 1.0 - 0.5 * M_SQRT2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tpp", "second next nearest neighbor hoping amplitude", 0.5 * (M_SQRT2 - 1.0));
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECheckerboardLatticeModelForgeState -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("input-state") == 0)
    {
      cout << "FQHECheckerboardLatticeModelForgeState requires an input state" << endl;
      return -1;
    }

  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey"); 
  double NNHoping = Manager.GetDouble("t1");
  double NextNNHoping = Manager.GetDouble("t2");
  double SecondNextNNHoping = Manager.GetDouble("tpp");
  double MuS = Manager.GetDouble("mu-s");

  bool Statistics = true;
  
  MultiColumnASCIIFile InputPositions;
  if (InputPositions.Parse(Manager.GetString("input-state")) == false)
    {
      InputPositions.DumpErrors(cout) << endl;
      return -1;
    }
  if (InputPositions.GetNbrLines() == 0)
    {
      cout << "no position provided in " << Manager.GetString("input-state") << endl;
      return -1;
    }
  int NbrParticles = InputPositions.GetNbrLines();
  int* SublatticeIndex =  InputPositions.GetAsIntegerArray(2);
  int* LatticeXPosition = InputPositions.GetAsIntegerArray(0);
  int* LatticeYPosition = InputPositions.GetAsIntegerArray(1);

  ComplexMatrix* OneBodyBasis = new ComplexMatrix [NbrSitesX * NbrSitesY];
  double GammaX = 0.0;
  double GammaY = 0.0;
  for (int kx = 0; kx < NbrSitesX; ++kx)
    for (int ky = 0; ky < NbrSitesY; ++ky)
      {
	int Index = (kx * NbrSitesY) + ky;
	Complex B1 = 4.0 * NNHoping * Complex (cos (1.0 * M_PI * (((double) kx) + GammaX) / ((double) NbrSitesX)) * cos (1.0 * M_PI * (((double) ky) + GammaY) / ((double) NbrSitesY)) * cos(M_PI * 0.25), 
						     sin (1.0 * M_PI * (((double) kx) + GammaX) / ((double) NbrSitesX)) * sin (1.0 * M_PI * (((double) ky) + GammaY) / ((double) NbrSitesY)) * sin(M_PI * 0.25));
	double d1 = 4.0 * SecondNextNNHoping * cos (2.0 * M_PI * (((double) kx) + GammaX) / ((double) NbrSitesX)) * cos (2.0 * M_PI * (((double) ky) + GammaY) / ((double) NbrSitesY));
	double d3 =  MuS + (2.0 * NextNNHoping * (cos (2.0 * M_PI * (((double) kx) + GammaX) / ((double) NbrSitesX))
							      - cos (2.0 * M_PI * (((double) ky) + GammaY) / ((double) NbrSitesY))));
	HermitianMatrix TmpOneBobyHamiltonian(2, true);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 0, d1 + d3);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 1, B1);
	TmpOneBobyHamiltonian.SetMatrixElement(1, 1, d1 - d3);
	ComplexMatrix TmpMatrix(2, 2, true);
	TmpMatrix[0][0] = 1.0;
	TmpMatrix[1][1] = 1.0;
	RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	TmpOneBobyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBobyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	OneBodyBasis[Index] = TmpMatrix;	
      }


  char* EigenstateListOutputFile = new char [512];
  if (Manager.GetDouble("mu-s") == 0.0)
    sprintf (EigenstateListOutputFile, "fermions_singleband_checkerboardlattice_twoparticles_n_%d_x_%d_y_%d_t1_%f_t2_%f_gx_%f_gy_%f.linear", NbrParticles, NbrSitesX, NbrSitesY, 
	     Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
  else
    sprintf (EigenstateListOutputFile, "fermions_singleband_checkerboardlattice_twoparticles_n_%d_x_%d_y_%d_t1_%f_t2_%f_gx_%f_gy_%f_mus_%f.linear", NbrParticles, NbrSitesX, NbrSitesY, 
	     Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
  ofstream File;
  File.precision(14);
  File.open(EigenstateListOutputFile, ios::binary | ios::out);
  File << "# vector_name coefficient" << endl;

  double Normalization = 1.0 / sqrt (((double) NbrSitesX) * ((double) NbrSitesY));
  double Sum = 0.0;
  for (int i = 0; i < NbrSitesX; ++i)
    {
      for (int j = 0; j < NbrSitesY; ++j)
	{
	  FermionOnSquareLatticeMomentumSpace* Space = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSitesX, NbrSitesY, i, j);
	  ComplexVector State (Space->GetHilbertSpaceDimension(), true);
	  int* MomentumArray = new int [2 * NbrParticles];
	  if (NbrParticles == 1)
	    {
	      MomentumArray[0] = i;
	      MomentumArray[1] = j;
	      int Index = Space->FindStateIndexFromArray(MomentumArray);
	      State[Index] += Phase(2.0 * M_PI * (((double) i) * ((0.5 * ((double) SublatticeIndex[0])) + (double) LatticeXPosition[0])) / ((double) NbrSitesX) 
				    +  2.0 * M_PI * (((double) j) * ((0.5 * ((double) SublatticeIndex[0])) + (double) LatticeYPosition[0])) / ((double) NbrSitesY)) * Normalization * (OneBodyBasis[(i * NbrSitesY) + j ][0][SublatticeIndex[0]]);
	    }
	  else
	    {
	      for (int kx1 = 0; kx1 < NbrSitesX; ++kx1)
		for (int kx2 = 0; kx2 < NbrSitesX; ++kx2)
		  {
		    if ((((kx1 + kx2) % NbrSitesX) == i) && (kx1 != kx2))
		      {
			for (int ky1 = 0; ky1 < NbrSitesY; ++ky1)
			  for (int ky2 = 0; ky2 < NbrSitesY; ++ky2)
			    {
			      if ((((ky1 + ky2) % NbrSitesY) == j) && (ky1 != ky2))
				{
				  MomentumArray[0] = kx1;
				  MomentumArray[1] = ky1;
				  MomentumArray[2] = kx2;
				  MomentumArray[3] = ky2;
				  int Index = Space->FindStateIndexFromArray(MomentumArray);
				  State[Index] += (Phase(2.0 * M_PI * (((double) kx1) * ((0.5 * ((double) SublatticeIndex[0])) + (double) LatticeXPosition[0])) / ((double) NbrSitesX) 
							 +  2.0 * M_PI * (((double) ky1) * ((0.5 * ((double) SublatticeIndex[0])) + (double) LatticeYPosition[0])) / ((double) NbrSitesY)) * Normalization * (OneBodyBasis[(kx1 * NbrSitesY) + ky1][0][SublatticeIndex[0]]))
				    * (Phase(2.0 * M_PI * (((double) kx2) * ((0.5 * ((double) SublatticeIndex[1])) + (double) LatticeXPosition[1])) / ((double) NbrSitesX) 
							 +  2.0 * M_PI * (((double) ky2) * ((0.5 * ((double) SublatticeIndex[1])) + (double) LatticeYPosition[1])) / ((double) NbrSitesY)) * Normalization * (OneBodyBasis[(kx2 * NbrSitesY) + ky2][0][SublatticeIndex[1]]));
				  
				}  
			    }
		      }
		  }
	    }
	  double Norm = State.Norm();
	  if (Norm > MACHINE_PRECISION)
	    {
	      Sum += Norm * Norm;
	      State /= Norm;
	      char* EigenstateOutputFile = new char [512];
	      if (Manager.GetDouble("mu-s") == 0.0)
		sprintf (EigenstateOutputFile, "fermions_singleband_checkerboardlattice_twoparticles_n_%d_x_%d_y_%d_t1_%f_t2_%f_gx_%f_gy_%f_kx_%d_ky_%d.0.vec", NbrParticles, NbrSitesX, NbrSitesY, 
			 Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j);
	      else
		sprintf (EigenstateOutputFile, "fermions_singleband_checkerboardlattice_twoparticles_n_%d_x_%d_y_%d_t1_%f_t2_%f_gx_%f_gy_%f_mus_%f_kx_%d_ky_%d.0.vec", NbrParticles, NbrSitesX, NbrSitesY, 
			 Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"), i, j);
	      State.WriteVector(EigenstateOutputFile);
	      Complex TmpComponent = Norm;
	      File << EigenstateOutputFile << " " << TmpComponent << endl;
	      delete[] EigenstateOutputFile;
	    }
	  delete Space;
	  delete[] MomentumArray;
	}
    }
  cout << "normalization = " << sqrt(Sum) << endl;

  File.close();
  delete[] EigenstateListOutputFile;

  return 0;
}
