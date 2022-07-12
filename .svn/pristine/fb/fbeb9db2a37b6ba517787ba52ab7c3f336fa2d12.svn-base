#include "Options/Options.h"
#include "Options/SingleStringOption.h"

#include "Tools/FTITightBinding/TightBindingModelKagomeLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"
#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"
#include "Tools/FTITightBinding/TightBindingModelCheckerboardLattice.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"


#include "Vector/ComplexVector.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/StringTools.h"

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
  OptionManager Manager ("FCIChernNumberFluxInsertion" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  
  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-cellx", "number of unit cells along the x direction", 5);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-celly", "number of unit cells along the y direction", 1);
  
  (*SystemGroup) += new SingleIntegerOption  ('X', "unit-cellx", "number of sites in unit cell along the x direction (for Hofstadter model)", 1);
  (*SystemGroup) += new SingleIntegerOption  ('Y', "unit-celly", "number of sites in unit cell along the y direction (for Hofstadter model)", 7);
  
  (*SystemGroup) += new BooleanOption  ('\n', "landau-x", "Use Landau gauge along the x-axis within unit cell");
  (*SystemGroup) += new BooleanOption  ('\n', "embedding", "compute the band structure with the embedding");
  (*SystemGroup) += new BooleanOption  ('\n', "checkerboard", "use checkerboard tightbiding model for both spins");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", 1.0 - 0.5 * M_SQRT2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tpp", "second next nearest neighbor hoping amplitude", 0.5 * (M_SQRT2 - 1.0));
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux-per-cell", "number of flux quanta per unit cell", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "band-index", "band-index", 0);
  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-gammax", "number of points in the discretization of gamma_x for the computation of Berry curvature", 10);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-gammay", "number of points in the discretization of gamma_y for the computation of Berry curvature", 10);
  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "singleparticle-precision", "length of floating point representations used for single-particle diagonalization, in bits", 64);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCIOFLNOrbitalTriangularLatticeModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
    
    
  
  bool CheckerboardFlag = Manager.GetBoolean("checkerboard");
  bool HofstadterFlag = false;
  if (CheckerboardFlag == false)
    HofstadterFlag = true;
  int NbrSitesX = Manager.GetInteger("nbr-cellx"); 
  int NbrSitesY = Manager.GetInteger("nbr-celly");
  int NbrStatesPerBand;
  
  int UnitCellX = Manager.GetInteger("unit-cellx"); 
  int UnitCellY = Manager.GetInteger("unit-celly");
  
  int FluxPerCell = Manager.GetInteger("flux-per-cell");
  
  int NbrBands = UnitCellX * UnitCellY / FluxPerCell;
  if (CheckerboardFlag)
    NbrBands = 2;
  
  char Axis ='y';
  
  if (Manager.GetBoolean("landau-x"))
    Axis ='x';
  
  int NbrPointX = Manager.GetInteger("nbr-gammax");
  int NbrPointY = Manager.GetInteger("nbr-gammay");
  int NbrPoint = 4;
  
  int IncNbrPointX =  NbrPointX+3;
  int IncNbrPointY =  NbrPointY+3;
  int BandIndex =  Manager.GetInteger("band-index");
  
  double TrueGammaX;
  double TrueGammaY;
  int TmpGammaIndex;
  
  
  bool ExportOneBody = true; 
  bool EmbeddingFlag = Manager.GetBoolean("embedding");
  int Precision = Manager.GetInteger("singleparticle-precision");
  
//   char* StatisticPrefix = new char [16];
     
  Abstract2DTightBindingModel *TightBindingModel;
  ComplexMatrix** OneBodyBasis = new ComplexMatrix*[NbrPoint];
  NbrStatesPerBand = NbrSitesX * NbrSitesY;
  double*** OneBodyEnergy = new double**[NbrBands];
  for (int i = 0; i < NbrBands; ++i)
    OneBodyEnergy[i] = new double*[NbrPoint];
  
  for (int GammaX = 0; GammaX < 2; ++GammaX)
    {
      for (int GammaY = 0; GammaY < 2; ++GammaY)
      {
	TmpGammaIndex = GammaX + GammaY * 2;
	TrueGammaX = ((double) GammaX) / ( (double) NbrPointX);
	TrueGammaY = ((double) GammaY) / ( (double) NbrPointY);
	
	for (int i = 0; i < NbrBands; ++i)
	  OneBodyEnergy[i][TmpGammaIndex] = new double[NbrStatesPerBand];
	
	if (HofstadterFlag)
	{
	  TightBindingModel = new TightBindingModelHofstadterSquare(NbrSitesX, NbrSitesY, UnitCellX, UnitCellY, FluxPerCell, Axis, TrueGammaX, TrueGammaY, Architecture.GetArchitecture(), ExportOneBody, EmbeddingFlag, Precision);
	  
	  
	}
	
	if (CheckerboardFlag)
	{
	  cout << "Checkerboard model" << endl;
	  TightBindingModel = new TightBindingModelCheckerboardLattice (NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), 
							     Manager.GetDouble("mu-s"), TrueGammaX, TrueGammaY, Architecture.GetArchitecture(), ExportOneBody);
	}

	OneBodyBasis[TmpGammaIndex] = new ComplexMatrix[NbrStatesPerBand];
	for (int TmpMomentumIndex= 0; TmpMomentumIndex < NbrStatesPerBand; ++ TmpMomentumIndex)
	{
	  OneBodyBasis[TmpGammaIndex][TmpMomentumIndex] = TightBindingModel->GetOneBodyMatrix(TmpMomentumIndex);	
	  for (int i = 0; i < NbrBands; ++i)
	    OneBodyEnergy[i][TmpGammaIndex][TmpMomentumIndex] = TightBindingModel->GetEnergy(i, TmpMomentumIndex);
	}
      }
    }
  
  
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);
  Complex TmpChernNumber = 0.0;
  Complex TmpChernNumberNoLog = 0.0;
  Complex TmpChernNumberHamiltonian = 0.0;
  Complex Tmp[4];
  Complex Tmp1[4];
  Complex Tmp2[4];
  
  Complex TmpX;
  Complex TmpY;
  
  double delta_2 = 4.0 * M_PI * M_PI / ((double) NbrPointX * NbrPointY * NbrSitesX * NbrSitesY);
  double GlobalFactor = 2.0 * M_PI / (delta_2 * NbrSitesX * NbrSitesY);
  
  for (int Kx = 0; Kx < NbrSitesX; ++Kx)
  {
   for (int Ky = 0; Ky < NbrSitesY; ++Ky)
   {
	  int TmpMomentumIndex = TightBindingModel->GetLinearizedMomentumIndex(Kx, Ky);
	  
	  ComplexMatrix& LocalBasis = OneBodyBasis[0][TmpMomentumIndex];
	  ComplexMatrix& LocalBasisIncX = OneBodyBasis[1][TmpMomentumIndex];
	  ComplexMatrix& LocalBasisIncY = OneBodyBasis[2][TmpMomentumIndex];
	  ComplexMatrix& LocalBasisIncXY = OneBodyBasis[3][TmpMomentumIndex];
      
  
	  Tmp[0] = 0.0;
	  Tmp[1] = 0.0;
	  Tmp[2] = 0.0;
	  Tmp[3] = 0.0;
	  
	  for (int i = 0; i < NbrBands; ++i)
	  {
	    Tmp[0] += Conj(LocalBasis[BandIndex][i]) * LocalBasisIncX[BandIndex][i];
	    Tmp[1] += Conj(LocalBasisIncX[BandIndex][i]) * LocalBasisIncXY[BandIndex][i];
	    Tmp[2] +=  Conj(LocalBasisIncXY[BandIndex][i]) * LocalBasisIncY[BandIndex][i];
	    Tmp[3] += Conj(LocalBasisIncY[BandIndex][i]) * LocalBasis[BandIndex][i];
	  }
	  
	  
	  
	  for (int TmpBandIndex1 = 0; TmpBandIndex1 < NbrBands; ++TmpBandIndex1)
	  {
	    if (TmpBandIndex1 != BandIndex)
	    {
	      double EnergyFactor = (OneBodyEnergy[BandIndex][0][TmpMomentumIndex] - OneBodyEnergy[TmpBandIndex1][0][TmpMomentumIndex]);
	      EnergyFactor *= EnergyFactor;
	      TmpX = 0.0;
	      TmpY = 0.0;
	      for (int TmpBandIndex2 = 0; TmpBandIndex2 < NbrBands; ++TmpBandIndex2)
	      {
		Tmp1[0] = 0.0;
		Tmp1[1] = 0.0;
		Tmp1[2] = 0.0;
		Tmp2[0] = 0.0;
		Tmp2[1] = 0.0;
		Tmp2[2] = 0.0;

		for (int i = 0; i < NbrBands; ++i)
		{
		  Tmp1[0] += Conj(LocalBasis[TmpBandIndex1][i]) * LocalBasisIncX[TmpBandIndex2][i];
		  Tmp1[2] += Conj(LocalBasis[TmpBandIndex1][i]) * LocalBasis[TmpBandIndex2][i];
		  Tmp1[1] += Conj(LocalBasis[TmpBandIndex1][i]) * LocalBasisIncY[TmpBandIndex2][i];

		  Tmp2[0] += Conj(LocalBasisIncX[TmpBandIndex2][i]) * LocalBasis[BandIndex][i];
		  Tmp2[2] += Conj(LocalBasis[TmpBandIndex2][i]) * LocalBasis[BandIndex][i];
		  Tmp2[1] += Conj(LocalBasisIncY[TmpBandIndex2][i]) * LocalBasis[BandIndex][i];
		}
		
// 		cout << Tmp1[0] << " " << Tmp2[0] << " " << Tmp1[1] << " " << Tmp1[1] << " " << Tmp1[2] << " " << Tmp1[2] << endl;
		TmpX += OneBodyEnergy[TmpBandIndex2][1][TmpMomentumIndex] * Tmp1[0] * Tmp2[0];
		TmpX -= OneBodyEnergy[TmpBandIndex2][0][TmpMomentumIndex] * Tmp1[2] * Tmp2[2];
		TmpY += OneBodyEnergy[TmpBandIndex2][2][TmpMomentumIndex] * Tmp1[1] * Tmp2[1];
		TmpY -= OneBodyEnergy[TmpBandIndex2][0][TmpMomentumIndex] * Tmp1[2] * Tmp2[2];
		
	      }
	      
	      TmpChernNumberHamiltonian += (Conj(TmpX) * TmpY / EnergyFactor);
// 	      cout << TmpX << " " << TmpY << " " << TmpChernNumberHamiltonian << endl;
	    }
	  }

	  TmpChernNumber += ln(Tmp[0] * Tmp[1] * Tmp[2] * Tmp[3]);
	  TmpChernNumberNoLog += (Tmp[0] * Tmp[1] * Tmp[2] * Tmp[3]);
	  
	  cout << Kx << " " << Ky << " "  << (ln(Tmp[0] * Tmp[1] * Tmp[2] * Tmp[3]) * GlobalFactor) << endl;
	
	}
      }
  TmpChernNumber *= GlobalFactor;
  TmpChernNumberNoLog *= GlobalFactor;
  TmpChernNumberHamiltonian *= 2.0*GlobalFactor;
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
  cout << "Chern number computed in  " << Dt << "s" << endl;
  
  cout << (TmpChernNumber.Im) << " " << (TmpChernNumberNoLog.Im)  << " " << (TmpChernNumberHamiltonian.Im) << endl;
  return 0;
  
  
  
  
}

