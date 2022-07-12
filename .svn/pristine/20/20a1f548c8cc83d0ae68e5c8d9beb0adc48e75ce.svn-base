#include "Options/Options.h"

#include "Tools/FTITightBinding/TightBindingModelKagomeLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"
#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"

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
  OptionManager Manager ("FCIHofstadterModelChernNumberFluxInsertion" , "0.01");
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
  
  (*SystemGroup) += new SingleIntegerOption  ('X', "unit-cellx", "number of sites in unit cell along the x direction", 1);
  (*SystemGroup) += new SingleIntegerOption  ('Y', "unit-celly", "number of sites in unit cell along the y direction", 7);
  
  (*SystemGroup) += new BooleanOption  ('\n', "landau-x", "Use Landau gauge along the x-axis within unit cell");
  (*SystemGroup) += new BooleanOption  ('\n', "embedding", "compute the band structure with the embedding");
  
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux-per-cell", "number of flux quanta per unit cell", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "band-index", "band-index", 0);
  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-gammax", "number of points in the discretization of gamma_x for the computation of Berry curvature", 10);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-gammay", "number of points in the discretization of gamma_y for the computation of Berry curvature", 10);
  (*SystemGroup) += new BooleanOption  ('\n', "allowed-momenta", "use flux insertion to compute Berry curvature, but compute Chern number based on allowed momenta only");
  
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
    
  int NbrCellX = Manager.GetInteger("nbr-cellx"); 
  int NbrCellY = Manager.GetInteger("nbr-celly");
  int NbrStatesPerBand = NbrCellX * NbrCellY;
  
  int UnitCellX = Manager.GetInteger("unit-cellx"); 
  int UnitCellY = Manager.GetInteger("unit-celly");
  
  int FluxPerCell = Manager.GetInteger("flux-per-cell");
  
  int NbrBands = UnitCellX * UnitCellY / FluxPerCell;
  
  char Axis ='y';
  
  if (Manager.GetBoolean("landau-x"))
    Axis ='x';
  
  int NbrPointX = Manager.GetInteger("nbr-gammax");
  int NbrPointY = Manager.GetInteger("nbr-gammay");
  int NbrPoint = NbrPointX * NbrPointY;
  bool AllowedMomentaOnly = Manager.GetBoolean("allowed-momenta");
  
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
  
  
  for (int GammaX = 0; GammaX < NbrPointX; ++GammaX)
    {
      for (int GammaY = 0; GammaY < NbrPointY; ++GammaY)
      {
	TmpGammaIndex = GammaX + GammaY * NbrPointX;
	TrueGammaX =   ( ( (double) GammaX) / ( (double) NbrPointX));
	TrueGammaY =   ( ( (double) GammaY) / ( (double) NbrPointY));
	TightBindingModel = new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis, TrueGammaX, TrueGammaY, Architecture.GetArchitecture(), ExportOneBody, EmbeddingFlag, Precision);
	
	OneBodyBasis[TmpGammaIndex] = new ComplexMatrix[NbrStatesPerBand];
	for (int TmpMomentumIndex= 0; TmpMomentumIndex < NbrStatesPerBand; ++ TmpMomentumIndex)
	  OneBodyBasis[TmpGammaIndex][TmpMomentumIndex] = TightBindingModel->GetOneBodyMatrix(TmpMomentumIndex);	
      }
    }
  
  
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);
  Complex TmpChernNumber = 0.0;
  Complex Tmp1[4];
  Complex Tmp2[8];
  
  int IncKX;
  int DecKX;
  int IncKY;
  int DecKY;
  
  for (int Kx = 0; Kx < NbrCellX; ++Kx)
  {
    for (int GammaX = 0; GammaX < NbrPointX; ++GammaX)
    {
      for (int Ky = 0; Ky < NbrCellY; ++Ky)
      {
	for (int GammaY = 0; GammaY < NbrPointY; ++GammaY)
	{
	  TrueGammaX =   ( ( (double) GammaX) / ( (double) NbrPointX));
	  TrueGammaY =   ( ( (double) GammaY) / ( (double) NbrPointY));
	  TmpGammaIndex = GammaX + GammaY * NbrPointX;
	  int TmpMomentumIndex = TightBindingModel->GetLinearizedMomentumIndex(Kx, Ky);
	  
	  IncKX = Kx;
	  if (GammaX == NbrPointX - 1)
	    IncKX = (Kx + 1) % NbrCellX;
	  IncKY = Ky;
	  if (GammaY == NbrPointY - 1)
	    IncKY = (Ky + 1) % NbrCellY;
	  
	  DecKX = Kx;
	  if (GammaX == 0)
	  {
	    if (Kx > 0)
	      DecKX = Kx - 1;
	    else
	      DecKX = NbrCellX - 1;
	  }
	    
	  DecKY = Ky;
	  if (GammaY == 0)
	  {
	    if (Ky > 0)
	      DecKY = Ky - 1;
	    else
	      DecKY = NbrCellY - 1;
	  }
	  
	  int LinearizedMomentumIndex = TightBindingModel->GetLinearizedMomentumIndex(Kx, Ky);
	  int LinearizedMomentumIndexIncX = TightBindingModel->GetLinearizedMomentumIndex(IncKX, Ky);
	  int LinearizedMomentumIndexDecX = TightBindingModel->GetLinearizedMomentumIndex(DecKX, Ky);
	  int LinearizedMomentumIndexIncY = TightBindingModel->GetLinearizedMomentumIndex(Kx, IncKY);
	  int LinearizedMomentumIndexDecY = TightBindingModel->GetLinearizedMomentumIndex(Kx, DecKY);
	  
	  int TmpGammaIndexIncX = ((GammaX + 1) % NbrPointX) + GammaY * NbrPointX;
	  int TmpGammaIndexDecX = (GammaX + NbrPointX - 1) % NbrPointX + GammaY * NbrPointX;
	  int TmpGammaIndexIncY = GammaX + NbrPointX * ((GammaY + 1) % NbrPointY);
	  int TmpGammaIndexDecY = GammaX + NbrPointX * ((GammaY + NbrPointY - 1) % NbrPointY);

	  ComplexMatrix& LocalBasis = OneBodyBasis[TmpGammaIndex][LinearizedMomentumIndex];
	  ComplexMatrix& LocalBasisIncX = OneBodyBasis[TmpGammaIndexIncX][LinearizedMomentumIndexIncX];
	  ComplexMatrix& LocalBasisDecX = OneBodyBasis[TmpGammaIndexDecX][LinearizedMomentumIndexDecX];
	  ComplexMatrix& LocalBasisIncY = OneBodyBasis[TmpGammaIndexIncY][LinearizedMomentumIndexIncY];
	  ComplexMatrix& LocalBasisDecY = OneBodyBasis[TmpGammaIndexDecY][LinearizedMomentumIndexDecY];  
	  Tmp1[0] = 0.0;
	  Tmp1[1] = 0.0;
	  Tmp1[2] = 0.0;
	  Tmp1[3] = 0.0;

	  Tmp2[0] = 0.0;
	  Tmp2[1] = 0.0;
	  Tmp2[2] = 0.0;
	  Tmp2[3] = 0.0;
	  Tmp2[4] = 0.0;
	  Tmp2[5] = 0.0;
	  Tmp2[6] = 0.0;
	  Tmp2[7] = 0.0;

	  for (int i = 0; i < NbrBands; ++i)
	  {
	    Tmp1[0] += LocalBasis[BandIndex][i] * Conj(LocalBasisIncX[BandIndex][i]);
	    Tmp1[1] += LocalBasis[BandIndex][i] * Conj(LocalBasisDecX[BandIndex][i]);
	    Tmp1[2] += LocalBasis[BandIndex][i] * Conj(LocalBasisIncY[BandIndex][i]);
	    Tmp1[3] += LocalBasis[BandIndex][i] * Conj(LocalBasisDecY[BandIndex][i]);

	    Tmp2[0] += Conj(LocalBasisIncX[BandIndex][i]) * LocalBasisIncY[BandIndex][i];
	    Tmp2[1] += Conj(LocalBasisDecX[BandIndex][i]) * LocalBasisIncY[BandIndex][i];
	    Tmp2[2] += Conj(LocalBasisIncX[BandIndex][i]) * LocalBasisDecY[BandIndex][i];
	    Tmp2[3] += Conj(LocalBasisDecX[BandIndex][i]) * LocalBasisDecY[BandIndex][i];
	    Tmp2[4] += Conj(LocalBasisIncY[BandIndex][i]) * LocalBasisIncX[BandIndex][i];
	    Tmp2[5] += Conj(LocalBasisDecY[BandIndex][i]) * LocalBasisIncX[BandIndex][i];
	    Tmp2[6] += Conj(LocalBasisIncY[BandIndex][i]) * LocalBasisDecX[BandIndex][i];
	    Tmp2[7] += Conj(LocalBasisDecY[BandIndex][i]) * LocalBasisDecX[BandIndex][i];
	  }

	if ((AllowedMomentaOnly == false) or ((GammaX == 0) and (GammaY == 0)))
	{
	  TmpChernNumber += (Tmp1[2] * Conj(Tmp1[0]) * Tmp2[0]);
	  TmpChernNumber -= (Tmp1[2] * Conj(Tmp1[1]) * Tmp2[1]);
	  TmpChernNumber -= (Tmp1[3] * Conj(Tmp1[0]) * Tmp2[2]);
	  TmpChernNumber += (Tmp1[3] * Conj(Tmp1[1]) * Tmp2[3]);
	  
	  TmpChernNumber -= (Tmp1[0] * Conj(Tmp1[2]) * Tmp2[4]);
	  TmpChernNumber += (Tmp1[0] * Conj(Tmp1[3]) * Tmp2[5]);
	  TmpChernNumber += (Tmp1[1] * Conj(Tmp1[2]) * Tmp2[6]);
	  TmpChernNumber -= (Tmp1[1] * Conj(Tmp1[3]) * Tmp2[7]);
	}
	}
      }
    }
  }
  TmpChernNumber /= 8.0 * M_PI;
  if (AllowedMomentaOnly)
    TmpChernNumber *= (NbrPointX * NbrPointY);
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
  cout << "Chern number computed in  " << Dt << "s" << endl;
  
  cout << (TmpChernNumber.Im) << endl;
  return 0;
  
  
  
  
}

