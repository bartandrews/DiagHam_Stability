#include "Options/Options.h"

#include "Tools/FTITightBinding/TightBindingModelKagomeLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"
#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"
#include "Tools/FTITightBinding/TightBindingModelHaldaneHoneycombLattice.h"


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
  OptionManager Manager ("FCIHofstadterModelComputeMatrixElementsFermiGoldenRule" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  
  ArchitectureManager Architecture;
//   LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
//   Lanczos.AddOptionGroup(&Manager);
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
  
  
   (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t3", "next to next nearest neighbor hoping amplitude", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "phi", "Haldane phase on nnn hopping (multiples of pi)", 1.0/3.0);
  (*SystemGroup) += new BooleanOption  ('\n', "phase-in-pi", "Haldane phase on nnn hopping given in multiples of pi");  
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  
  
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
  
  double TrueGammaX;
  double TrueGammaY;
  int TmpGammaIndex;
  
  
  bool ExportOneBody = true; 
  bool EmbeddingFlag = Manager.GetBoolean("embedding");
  int Precision = Manager.GetInteger("singleparticle-precision");
  
  
  char* FileName = new char [128];
  sprintf(FileName, "transition_elements_hofstadter_X_%d_Y_%d_q_%d_x_%d_y_%d_nbrgx_%d_nbrgy_%d.dat", UnitCellX, UnitCellY, FluxPerCell, NbrCellX, NbrCellY, NbrPointX, NbrPointY);
  ofstream File;
  File.open(FileName, ios::out);
  File.precision(14);
  File << "# omega band kx ky V+ V+/omega V- V-/omega"<< endl;
     
  Abstract2DTightBindingModel *TightBindingModel;
  ComplexMatrix** OneBodyBasis = new ComplexMatrix*[NbrPoint];
  double*** OneBodyEnergy = new double**[NbrBands];
  for (int i = 0; i < NbrBands; ++i)
    OneBodyEnergy[i] = new double*[NbrPoint];
  
  
  for (int GammaX = 0; GammaX < NbrPointX; ++GammaX)
    {
      for (int GammaY = 0; GammaY < NbrPointY; ++GammaY)
      {
	TmpGammaIndex = GammaX + GammaY * NbrPointX;
	TrueGammaX =   ( ( (double) GammaX) / ( (double) NbrPointX));
	TrueGammaY =   ( ( (double) GammaY) / ( (double) NbrPointY));
	TightBindingModel = new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis, TrueGammaX, TrueGammaY, Architecture.GetArchitecture(), ExportOneBody, EmbeddingFlag, Precision);
	
	
// 	double HaldanePhi;
// 
//   if (Manager.GetBoolean("phase-in-pi"))
//     HaldanePhi = M_PI*Manager.GetDouble("phi");
//   else
//     HaldanePhi = Manager.GetDouble("phi");
// 	TightBindingModel = new TightBindingModelHaldaneHoneycombLattice (NbrCellX, NbrCellY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), 
// 								  HaldanePhi, Manager.GetDouble("mu-s"), TrueGammaX, TrueGammaY, Architecture.GetArchitecture(), ExportOneBody);
// 	
// 	
	
	OneBodyBasis[TmpGammaIndex] = new ComplexMatrix[NbrStatesPerBand];
	for (int i = 0; i < NbrBands; ++i)
	  OneBodyEnergy[i][TmpGammaIndex] = new double[NbrStatesPerBand];
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
  Complex Tmp1[4];
  Complex Tmp2[8];
  
  int IncKX;
  int DecKX;
  int IncKY;
  int DecKY;
  
  int GammaX = 0;
  int GammaY = 0;
  
  double Omega;
  
  Complex VgeElement;
  for (int TmpBandIndex1 = 1; TmpBandIndex1 < NbrBands; ++TmpBandIndex1)
  {
    for (int Kx = 0; Kx < NbrCellX; ++Kx)
    {
      for (int Ky = 0; Ky < NbrCellY; ++Ky)
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
	 
	 
	int LinearizedMomentumIndex = TightBindingModel->GetLinearizedMomentumIndex(Kx, Ky);
	int LinearizedMomentumIndexIncX = TightBindingModel->GetLinearizedMomentumIndex(IncKX, Ky);
	int LinearizedMomentumIndexIncY = TightBindingModel->GetLinearizedMomentumIndex(Kx, IncKY);
	     
	int TmpGammaIndexIncX = ((GammaX + 1) % NbrPointX) + GammaY * NbrPointX;
	int TmpGammaIndexIncY = GammaX + NbrPointX * ((GammaY + 1) % NbrPointY);
	ComplexMatrix& LocalBasis = OneBodyBasis[TmpGammaIndex][LinearizedMomentumIndex];
	ComplexMatrix& LocalBasisIncX = OneBodyBasis[TmpGammaIndexIncX][LinearizedMomentumIndexIncX];
	ComplexMatrix& LocalBasisIncY = OneBodyBasis[TmpGammaIndexIncY][LinearizedMomentumIndexIncY];
	
	Omega = OneBodyEnergy[TmpBandIndex1][TmpGammaIndex][LinearizedMomentumIndex] - OneBodyEnergy[0][TmpGammaIndex][LinearizedMomentumIndex];
	File << Omega << " " << TmpBandIndex1 << " " << Kx << " " << Ky << " " ;
	
	for (int Tmp = 0; Tmp < 2; ++Tmp)
	{
	  int TmpSign = Tmp * 2 - 1;
	  VgeElement = 0.0;
	  for (int TmpBandIndex2 = 0; TmpBandIndex2 < NbrBands; ++TmpBandIndex2)
	  {
	
// 	  for (int GammaX = 0; GammaX < NbrPointX; ++GammaX)
// 	  {
	    
// 	      for (int GammaY = 0; GammaY < NbrPointY; ++GammaY)
// 	      {
		
	  
		Tmp1[0] = 0.0;
		Tmp1[1] = 0.0;
		Tmp1[2] = 0.0;

		Tmp2[0] = 0.0;
		Tmp2[1] = 0.0;
		Tmp2[2] = 0.0;

		for (int i = 0; i < NbrBands; ++i)
		{
		  Tmp1[0] += Conj(LocalBasis[TmpBandIndex1][i]) * LocalBasisIncX[TmpBandIndex2][i];
		  Tmp1[1] += Conj(LocalBasis[TmpBandIndex1][i]) * LocalBasis[TmpBandIndex2][i];
		  Tmp1[2] += Conj(LocalBasis[TmpBandIndex1][i]) * LocalBasisIncY[TmpBandIndex2][i];

		  Tmp2[0] += Conj(LocalBasisIncX[TmpBandIndex2][i]) * LocalBasis[0][i];
		  Tmp2[1] += Conj(LocalBasis[TmpBandIndex2][i]) * LocalBasis[0][i];
		  Tmp2[2] += Conj(LocalBasisIncY[TmpBandIndex2][i]) * LocalBasis[0][i];
		}

		if ((AllowedMomentaOnly == false) or ((GammaX == 0) and (GammaY == 0)))
		{
		  VgeElement += (OneBodyEnergy[TmpBandIndex2][TmpGammaIndexIncX][LinearizedMomentumIndexIncX] * Tmp1[0] * Tmp2[0] - OneBodyEnergy[TmpBandIndex2][TmpGammaIndex][LinearizedMomentumIndex] * Tmp1[1] * Tmp2[1]) * NbrPointX * NbrCellX * Complex (0.0, -1.0) / (2.0 * M_PI);
		  VgeElement += (OneBodyEnergy[TmpBandIndex2][TmpGammaIndexIncY][LinearizedMomentumIndexIncY] * Tmp1[2] * Tmp2[2] - OneBodyEnergy[TmpBandIndex2][TmpGammaIndex][LinearizedMomentumIndex] * Tmp1[1] * Tmp2[1]) * NbrPointY * NbrCellY * (-TmpSign)  / (2.0 * M_PI);
		}
// 	      }
// 	    }
	  }
//   if (AllowedMomentaOnly)
//     VgeElement *= (NbrPointX * NbrPointY);
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));    
  cout << "Chern number computed in  " << Dt << "s" << endl;
  
  cout << (VgeElement) << endl;
      double V = 2.0 * M_PI * (VgeElement.Re*VgeElement.Re + VgeElement.Im*VgeElement.Im);
      
      File << V << " " << (V / (Omega*Omega)) << " " ;
	}
	
      File << endl;
      }
    }
  }
  File.close();
  return 0;
  
  
  
  
}

