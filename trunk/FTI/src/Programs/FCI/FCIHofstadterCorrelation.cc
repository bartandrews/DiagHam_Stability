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

#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"


#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <cassert>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;
//new namespaces (added by ba340)
using std::setw;

int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FCIHofstadterCorrelation" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('x', "reference-x", "x-coordinate or reference site", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "reference-y", "y-coordinate or reference site", 0);
  
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
  int NbrCellX = 0;
  int NbrCellY = 0;
  int MomentumX = 0;
  int MomentumY = 0;
  int UnitCellX=0;
  int UnitCellY=0;     
  int FluxPerCell=0;   
  char Axis='y';         
  double Interaction=0; 
  double GammaX=0;
  double GammaY=0;
  bool EmbeddingFlag=false;
  bool Hardcore=false;	
  int NbrState=0; 
  int NbrSamplesX = Manager.GetInteger("nbr-samplesx");
  int NbrSamplesY = Manager.GetInteger("nbr-samplesy");
  bool DensityFlag = Manager.GetBoolean("density");
  bool Statistics = true;
  int NbrBands = 1;
  bool EnlargeCell = false;
  double MuPotential = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  
  
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName_Hofstadter(Manager.GetString("state"), NbrParticles, NbrCellX, NbrCellY, Interaction, FluxPerCell, NbrState, Statistics, Hardcore, EmbeddingFlag, Axis, GammaX, GammaY, MomentumX, MomentumY, UnitCellX, UnitCellY, EnlargeCell, MuPotential, NbrBands) == false)
    {
      if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName_Hofstadter(Manager.GetString("state"), NbrParticles, NbrCellX, NbrCellY, Interaction, FluxPerCell, NbrState, Statistics, Hardcore, EmbeddingFlag, Axis, GammaX, GammaY, MomentumX, MomentumY, UnitCellX, UnitCellY) == false)
      {
	cout << "error while retrieving system parameters from file name " << Manager.GetString("state") << endl;
	return -1;
      }
    }
   
  int FullUnitCellX = UnitCellX;
  if (EnlargeCell)
  {
    FullUnitCellX = 2 * UnitCellX;
    cout << "Warning: autodetect of t1 and t2 not implemented" << endl;
  }

  cout << setw(20) << std::left << "Statistics" << setw(20) << std::left << Statistics << endl;
  cout << setw(20) << std::left << "UnitCellX" << setw(20) << std::left << UnitCellX << endl;
  cout << setw(20) << std::left << "UnitCellY" << setw(20) << std::left << UnitCellY << endl;
  cout << setw(20) << std::left << "FluxPerCell" << setw(20) << std::left << FluxPerCell << endl;
  cout << setw(20) << std::left << "Axis" << setw(20) << std::left << Axis << endl;
  cout << setw(20) << std::left << "NbrParticles" << setw(20) << std::left << NbrParticles << endl;
  cout << setw(20) << std::left << "NbrCellX" << setw(20) << std::left << NbrCellX << endl;
  cout << setw(20) << std::left << "NbrCellY" << setw(20) << std::left << NbrCellY << endl;
  cout << setw(20) << std::left << "Hardcore" << setw(20) << std::left << Hardcore << endl;
  cout << setw(20) << std::left << "Interaction" << setw(20) << std::left << Interaction << endl;
  cout << setw(20) << std::left << "GammaX" << setw(20) << std::left << GammaX << endl;
  cout << setw(20) << std::left << "GammaY" << setw(20) << std::left << GammaY << endl;
  cout << setw(20) << std::left << "EmbeddingFlag" << setw(20) << std::left << EmbeddingFlag << endl;
  cout << setw(20) << std::left << "MomentumX" << setw(20) << std::left << MomentumX << endl;
  cout << setw(20) << std::left << "MomentumY" << setw(20) << std::left << MomentumY << endl;
  cout << setw(20) << std::left << "NbrState" << setw(20) << std::left << NbrState << endl;
  cout << setw(20) << std::left << "NbrBands" << setw(20) << std::left << NbrBands << endl;
  cout << setw(20) << std::left << "Mus" << setw(20) << std::left << MuPotential << endl;

  ParticleOnSphere* Space = 0;
  if (Statistics == true)
    Space = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrCellX, NbrCellY, MomentumX, MomentumY);
  else
  {
    switch(NbrBands)
    {
      case 1:
	Space = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrCellX, NbrCellY, MomentumX, MomentumY);
	break;
      case 2:
	Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace(NbrParticles, NbrCellX, NbrCellY, MomentumX, MomentumY);
	break;
    }
  }
  ComplexVector ComplexState;
  if (ComplexState.ReadVector (Manager.GetString("state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("state") << endl;
      return -1;      
    }
  Complex* PrecalculatedValues_rho = 0;
  Complex* PrecalculatedValues_rhorho = 0;
  int* PrecalculatedIndices = 0;
  int NbrPrecalculatedValues = 0;

  // use tight binding model to provide function basis and mappings of momenta
  TightBindingModelHofstadterSquare *tightBindingModel;
  if (EnlargeCell == false)
    tightBindingModel = new TightBindingModelHofstadterSquare (NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis, GammaX, GammaY, Architecture.GetArchitecture(), true, EmbeddingFlag);
  else
    tightBindingModel = new TightBindingModelHofstadterSquare (NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, MuPotential, FullUnitCellX, t1, t2, Axis, GammaX, GammaY, Architecture.GetArchitecture(), true, EmbeddingFlag);

  if (DensityFlag == false)
    {
      for (int kx1 =0; kx1 < NbrCellX; ++kx1)
	for (int kx2 =0; kx2 < NbrCellX; ++kx2)
	  for (int kx3 =0; kx3 < NbrCellX; ++kx3)
	    for (int kx4 =0; kx4 < NbrCellX; ++kx4)
	      {
		if (((kx1 + kx2 - kx3 - kx4) % NbrCellX) == 0)
		  {
		    for (int ky1 = 0; ky1 < NbrCellY; ++ky1)
		      for (int ky2 = 0; ky2 < NbrCellY; ++ky2)
			for (int ky3 = 0; ky3 < NbrCellY; ++ky3)
			  for (int ky4 = 0; ky4 < NbrCellY; ++ky4)
			    {
			      if (((ky1 + ky2 - ky3 - ky4) % NbrCellY) == 0)
				{
				  ++NbrPrecalculatedValues;
				}
			    }
		  }
	      }
      PrecalculatedValues_rhorho = new Complex [NbrPrecalculatedValues];
      PrecalculatedIndices = new int [4 * NbrPrecalculatedValues];
      NbrPrecalculatedValues = 0; 
      for (int kx1 =0; kx1 < NbrCellX; ++kx1)
	for (int kx2 =0; kx2 < NbrCellX; ++kx2)
	  for (int kx3 =0; kx3 < NbrCellX; ++kx3)
	    for (int kx4 =0; kx4 < NbrCellX; ++kx4)
	      {
		if (((kx1 + kx2 - kx3 - kx4) % NbrCellX) == 0)
		  {
		    for (int ky1 = 0; ky1 < NbrCellY; ++ky1)
		      for (int ky2 = 0; ky2 < NbrCellY; ++ky2)
			for (int ky3 = 0; ky3 < NbrCellY; ++ky3)
			  for (int ky4 = 0; ky4 < NbrCellY; ++ky4)
			    {
			      if (((ky1 + ky2 - ky3 - ky4) % NbrCellY) == 0)
				{
				  int Index1 = tightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
				  int Index2 = tightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
				  int Index3 = tightBindingModel->GetLinearizedMomentumIndex(kx3, ky3);
				  int Index4 = tightBindingModel->GetLinearizedMomentumIndex(kx4, ky4);
				  ParticleOnSphereDensityDensityOperator Operator (Space, Index1, Index2, Index3, Index4);
				  PrecalculatedValues_rhorho[NbrPrecalculatedValues] = Operator.MatrixElement(ComplexState, ComplexState);
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
  // else
  //   {
  //     NbrPrecalculatedValues = NbrCellX * NbrCellY;
  //     PrecalculatedValues_rhorho = new Complex [NbrPrecalculatedValues];
  //     for (int kx =0; kx < NbrCellX; ++kx)
  // 	for (int ky = 0; ky < NbrCellY; ++ky)
  // 	  {
  // 	    int Index = (kx * NbrCellY) + ky;
  // 	    ParticleOnSphereDensityOperator Operator (Space, Index);	    
  // 	    PrecalculatedValues_rho[Index] = Operator.MatrixElement(ComplexState, ComplexState);
  // 	  }
  //   }
   
  // PrecalculatedValues_rho = new Complex [NbrPrecalculatedValues];
  PrecalculatedValues_rho = new Complex [NbrCellX * NbrCellY];
  for (int kx =0; kx < NbrCellX; ++kx)
    for (int ky = 0; ky < NbrCellY; ++ky)
      {
	int Index = tightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	ParticleOnSphereDensityOperator Operator (Space, Index);	    
	PrecalculatedValues_rho[Index] = Operator.MatrixElement(ComplexState, ComplexState);
      }
      
  delete Space;
  ofstream File;
  File.precision(14);
  double XStep = ((double) NbrCellX) / ((double) NbrSamplesX);
  double YStep = ((double) NbrCellY) / ((double) NbrSamplesY);
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
	  for (int kx =0; kx < NbrCellX; ++kx)
	    for (int ky = 0; ky < NbrCellY; ++ky)
	      {
		int Index2 = (kx * NbrCellY) + ky;
		int Index = tightBindingModel->GetLinearizedMomentumIndex(kx, ky);
		assert(Index == Index2);
		File << kx << " " << ky << " " << PrecalculatedValues_rho[Index].Re << endl;
	      }
	  File.close();
	}
      return 0;
    }

  //
  int TotalNbrMomenta = NbrCellX * NbrCellY;
  int NbrSublattices = FullUnitCellX * UnitCellY;
  Complex* Coefficients = new Complex[TotalNbrMomenta];
  Complex* Coefficients2 = new Complex[TotalNbrMomenta];
  Position[0] = Manager.GetInteger("reference-x");
  Position[1] = Manager.GetInteger("reference-y");
  for (int i = 0; i < TotalNbrMomenta; ++i){
    tightBindingModel->GetFunctionValue(Coefficients[i], Position, i, 0); 
  }
    
  //routine to calculate the two-particle correlation function of the form <psi|n_i n_0|psi> - <psi|n_i|psi><psi|n_0|psi>, with <psi|n_0|psi>(<psi|n_0|psi>-1) at the origin
  //
  double Normalisation = (1.0/(double)(NbrCellX*NbrCellY)); //normalisation factor for the <c^+ c> term i.e. 1/N_c
  double Normalisation2 = (1.0/(double)(NbrCellX*NbrCellY*NbrCellX*NbrCellY)); //normalisation factor for the <c^+ c^+ c c> term i.e. 1/N_c^2

  double DensityPrefactor = ( ( double ) ( TotalNbrMomenta*NbrSublattices ) / ( double ) NbrParticles );

  RealVector NumTranslations(2);
  int subX, subY;
  double *Correlations = new double[NbrCellX*NbrCellY*NbrSublattices];
  
  for (int Rjx = 0; Rjx < NbrCellX; ++Rjx) //loop over MUCs
    {
      NumTranslations[0]=Rjx;
      for (int Rjy = 0; Rjy < NbrCellY; ++Rjy)
	{
	  NumTranslations[1]=Rjy;
	  for (int alphaJ=0; alphaJ < NbrSublattices; ++alphaJ) // sublattice for r_j
	    {
	      tightBindingModel->GetSitePosition(Position, NumTranslations, alphaJ);
              
              Complex TmpValue = 0.0; //<psi|n_i n_0|psi>
              
              for (int i = 0; i < TotalNbrMomenta; ++i)
	      {
		tightBindingModel->GetFunctionValue(Coefficients2[i], Rjx, Rjy, alphaJ, i, 0);
	      }
	      
	      if (DensityFlag == false)
		{
		  for (int i = 0; i < NbrPrecalculatedValues; ++i)
		    {
		      TmpValue += DensityPrefactor * Normalisation2*(PrecalculatedValues_rhorho[i] 
						  * Conj(Coefficients[PrecalculatedIndices[(i << 2)]])
						  * Coefficients[PrecalculatedIndices[(i << 2) + 3]] 
						  * Conj(Coefficients2[PrecalculatedIndices[(i << 2) + 1]])
						  * Coefficients2[PrecalculatedIndices[(i << 2) + 2]]); //calculate <psi|n_i n_0|psi>
		    }
		  TmpValue-=( ( double ) NbrParticles / ( double ) ( TotalNbrMomenta*NbrSublattices ) );
		}
	      else
		{
		  for (int i = 0; i < TotalNbrMomenta; ++i) 
		    {
		      TmpValue += PrecalculatedValues_rho[i] * SqrNorm(Coefficients2[i]);
		    }
		  TmpValue *= Normalisation;
		}
              cout << Position[0] << " " << Position[1] << " " << TmpValue.Re << endl;
	      Correlations[tightBindingModel->GetRealSpaceTightBindingLinearizedIndex(Rjx, Rjy, alphaJ)] = TmpValue.Re;
	      assert (fabs(TmpValue.Im)<1e-12);
            }
	}
    }
  // output in format suitable for Gnuplot:
  int tx, ty, subl;
  for (int x=0; x< NbrCellX*FullUnitCellX; ++x)
    {
      Position[0]=x;
      for (int y=0; y< NbrCellY*UnitCellY; ++y)
	{
	  Position[1]=y;
	  tightBindingModel->PositionToLatticeCoordinates(Position, tx, ty, subl);
	  File << Position[0] << " " << Position[1] << " " << Correlations[tightBindingModel->GetRealSpaceTightBindingLinearizedIndex(tx, ty, subl)] << endl;
	}
      File << endl;
    }
  
  File << endl;
  File.close();
  if (PrecalculatedValues_rhorho!=0) delete [] PrecalculatedValues_rhorho;
  delete[] PrecalculatedValues_rho;
  delete[] Coefficients;
  delete[] Coefficients2;
  return 0;
}
