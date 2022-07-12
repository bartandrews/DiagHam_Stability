#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "HilbertSpace/Spin1_2ChainNewAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainNewSzSymmetryAnd2DTranslation.h"

#include "Operator/AbstractOperator.h"
#include "Operator/SpinWith2DTranslationSpinSpinCorrelationOperator.h"
#include "Operator/SpinWith2DTranslationBondBondCorrelationOperator.h"


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


// get a linearized position index from the 2d coordinates
//
// xPosition = position along the x direction
// yPosition = position along the y direction
// return value = linearized index

int GetLinearizedIndex(int xPosition, int yPosition, int atomicIndex, int nbrSpinX, int nbrSpinY);

int main(int argc, char** argv)
{
  OptionManager Manager ("SpinKagomeComputeSpinSpinCorrelations" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
  Manager += PrecalculationGroup;

  
  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "name of the file corresponding to the state |Psi_R> (the order parameter being <Psi_L|c^+c^+|Psi_R>");   
  (*SystemGroup) += new BooleanOption  ('\n', "bond", "calculate bond-bond correlations");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-x", "only evalute a given x position (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-y", "only evalute a given y position sector (negative if all ky sectors have to be computed)", -1); 
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Sz symmetry)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Sz symmetry)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinKagomeComputeSpinSpinCorrelations -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSites = 0;
  int XMomentum = 0;
  int YMomentum = 0;
  int InversionSector = 0;
  int SzValue = 0;
  int SzParitySector = 0;
  bool TotalSpinConservedFlag = true;
  bool InversionFlag = false;
  bool SzSymmetryFlag = false;
  int XPeriodicity = 0;
  int YPeriodicity = 0;
  bool Statistics = true;
  int NbrInputStates = 0;
  int SpinValue = 1;
  int Offset = 0;
  
  bool BondFlag = Manager.GetBoolean("bond");
  
  if (Manager.GetString("input-state") == 0)
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type SpinSystemConvertFromTranslationInvariantBasis -h" << endl;
      return -1;
    }

  else
    {
      InversionSector = 0;
      SzParitySector = 0;
      if (IsFile(Manager.GetString("input-state")) == false)
	{
	  cout << "can't open file " << Manager.GetString("input-state") << endl;
	  return -1;
	}   
	
      SpinFindSystemInfoFromVectorFileName (Manager.GetString("input-state"), NbrSites, SzValue, SpinValue, XMomentum, InversionSector, SzParitySector, Offset);
      
//       (char* filename, int& nbrSpins, int& sz, int& spin, int& momentum, int& inversion, int& szSymmetry, int& offset)
      cout << Offset << endl;
      if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SzValue, SpinValue, XMomentum, XPeriodicity,
								YMomentum, YPeriodicity) == false)
	{
	  TotalSpinConservedFlag = false;
	  if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SpinValue, XMomentum, XPeriodicity, 
								    YMomentum, YPeriodicity) == false)
	    {
	      cout << "error while retrieving system parameters from file name " <<Manager.GetString("input-state")  << endl;
	      return -1;
	    }
	  InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SpinValue, XMomentum, XPeriodicity, 
											 YMomentum, YPeriodicity, InversionSector);
	  if (InversionFlag == false)
	    cout << "error" << endl;
	  else
	    cout << "OK" << " " << InversionSector << endl;
	}
      else
	{
	  InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SzValue, SpinValue, XMomentum, XPeriodicity,
											 YMomentum, YPeriodicity, InversionSector);
	  if (SzValue == 0)
	  {
	    SpinFindSystemInfoFromVectorFileName((Manager.GetString("input-state")), NbrSites, SzValue, SpinValue, InversionSector, SzParitySector);
	    cout << NbrSites << " " << SzValue << " " << SpinValue << " " << XMomentum << " " << XPeriodicity << " " << YMomentum << " " <<  YPeriodicity << " " << SzParitySector << endl;
	    if (SzParitySector != 0)
	      SzSymmetryFlag = true;
	  }
	}
    }
   
      
  ComplexVector State;
  if (Manager.GetString("input-state") != 0)
    {
      if (State.ReadVector (Manager.GetString("input-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("input-state") << endl;
	  return -1;      
	}
    }

  Spin1_2ChainNewAnd2DTranslation* Space = 0;
      
   if (SzSymmetryFlag == false)
    Space = new Spin1_2ChainNewAnd2DTranslation (NbrSites, SzValue, XMomentum, XPeriodicity, YMomentum, YPeriodicity);
   else
   {
      cout << "Create HilbertSpace" << endl;
      if (Manager.GetString("load-hilbert") != 0)
	{
	  Space = new Spin1_2ChainNewSzSymmetryAnd2DTranslation(Manager.GetString("load-hilbert"));
	}
      else
	Space = new Spin1_2ChainNewSzSymmetryAnd2DTranslation (NbrSites, SzValue, (1 - SzParitySector)/2, XMomentum, XPeriodicity, YMomentum, YPeriodicity);
   }

  

  if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
    {
      cout << "error, " << Manager.GetString("input-state")  << " has a wrong dimension (" << State.GetVectorDimension() << ", should be " << Space->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  
  if (Manager.GetString("save-hilbert") != 0)
  {
    Space->WriteHilbertSpace(Manager.GetString("save-hilbert"));
    return 0;
  }
 
  ofstream File;
  ofstream File1;
  char* OutputFileName;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      File.open(OutputFileName, ios::binary | ios::out);
    }
  else
    {
      if (Manager.GetString("input-state") != 0)
	{
	  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("input-state"), "vec", "neighboring_spin_correlation.dat");
	  if (OutputFileName == 0)
	    {
	      cout << "no vec extension was find in " << Manager.GetString("input-state") << " file name" << endl;
	      return 0;
	    }
	  File.open(OutputFileName, ios::binary | ios::out);
	  
	  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("input-state"), "vec", "spinspin.dat");
	  if (OutputFileName == 0)
	    {
	      cout << "no vec extension was find in " << Manager.GetString("input-state") << " file name" << endl;
	      return 0;
	    }
	  File1.open(OutputFileName, ios::binary | ios::out);
	}
    }
  File.precision(14);
  File1.precision(14);
  cout.precision(14);
  
  cout << "NbrSites = " << NbrSites << " " << " XPeriodicity = " << XPeriodicity << " YPeriodicity  = " << YPeriodicity << " Sz = " << SzValue << " SzParitySector = " << SzParitySector << " XMomentum = " << XMomentum << " YMomentum = " << YMomentum << endl;
  
  AbstractOperator* Operator;
  Complex* NeighborSpinSpinCorrelation = new Complex [6];
  int** TmpIndices = new int* [6];
  for (int i = 0; i < 6; ++i)
  {
    TmpIndices[i] = new int[2];
    NeighborSpinSpinCorrelation[i] = 0.0;
  }
  
  int MinX = 0;
  int MaxX = XPeriodicity - 1;
  int MinY = 0;
  int MaxY = YPeriodicity - 1;
  
  if (Manager.GetInteger("only-x") != -1)
  {
    MinX = Manager.GetInteger("only-x");
    MaxX = Manager.GetInteger("only-x");
  }
  if (Manager.GetInteger("only-y") != -1)
  {
    MinY = Manager.GetInteger("only-y");
    MaxY = Manager.GetInteger("only-y");
  }
  
  File << "# E_{AB} E_{BA} E_{BC} E_{CB} E_{AC} E_{CA}" << endl;
  cout << "# E_{AB} E_{BA} E_{BC} E_{CB} E_{AC} E_{CA}" << endl;
  for (int i = 0; i < XPeriodicity; ++i)
  {
    for (int j = 0; j < YPeriodicity; ++j)
    {
      TmpIndices[0][0] = GetLinearizedIndex(i, j, 0, XPeriodicity, YPeriodicity);
      TmpIndices[0][1] = GetLinearizedIndex(i, j, 1, XPeriodicity, YPeriodicity);
      TmpIndices[1][0] = GetLinearizedIndex(i, j, 0, XPeriodicity, YPeriodicity);
      TmpIndices[1][1] = GetLinearizedIndex(i - 1, j, 1, XPeriodicity, YPeriodicity);
      TmpIndices[2][0] = GetLinearizedIndex(i, j, 1, XPeriodicity, YPeriodicity);
      TmpIndices[2][1] = GetLinearizedIndex(i, j, 2, XPeriodicity, YPeriodicity);
      TmpIndices[3][0] = GetLinearizedIndex(i - Offset, j - 1, 2, XPeriodicity, YPeriodicity);
      TmpIndices[3][1] = GetLinearizedIndex(i - 1, j, 1, XPeriodicity, YPeriodicity);
      TmpIndices[4][0] = GetLinearizedIndex(i, j, 0, XPeriodicity, YPeriodicity);
      TmpIndices[4][1] = GetLinearizedIndex(i, j, 2, XPeriodicity, YPeriodicity);
      TmpIndices[5][0] = GetLinearizedIndex(i, j, 0, XPeriodicity, YPeriodicity);
      TmpIndices[5][1] = GetLinearizedIndex(i - Offset, j - 1, 2, XPeriodicity, YPeriodicity);
  
//       cout << i << " " << j << " " ;
//       File << i << " " << j << " " ;
      for (int l = 0; l < 6; ++l)
      {
	Operator = new SpinWith2DTranslationSpinSpinCorrelationOperator(Space, XMomentum, XPeriodicity, YMomentum, YPeriodicity, TmpIndices[l][0], TmpIndices[l][1]);
	OperatorMatrixElementOperation Operation(Operator, State, State, State.GetVectorDimension());
	Operation.ApplyOperation(Architecture.GetArchitecture());
	NeighborSpinSpinCorrelation[l] += Operation.GetScalar();
// 	cout << NeighborSpinSpinCorrelation[l] << " " ;
	delete Operator;
      }
//       cout << endl;
      
    }
  }
  
  Complex NematicOP = (-NeighborSpinSpinCorrelation[0] + NeighborSpinSpinCorrelation[1] + NeighborSpinSpinCorrelation[2] * Phase( - M_PI / 3.0 ) + NeighborSpinSpinCorrelation[3] *  Phase( 2.0 * M_PI / 3.0 ) + NeighborSpinSpinCorrelation[4]  * Phase( M_PI / 3.0 ) + NeighborSpinSpinCorrelation[5]  * Phase( -2.0 * M_PI / 3.0 )) / (XPeriodicity * YPeriodicity);
  for (int l = 0; l < 6; ++l)
  {
    NeighborSpinSpinCorrelation[l] = NeighborSpinSpinCorrelation[l] / (XPeriodicity * YPeriodicity);
    File << NeighborSpinSpinCorrelation[l]  << " " ;
    cout << NeighborSpinSpinCorrelation[l] << " " ;
  }
  File << endl;
  cout << endl;
  File << "C3 OP = " << sqrt(NematicOP.Re*NematicOP.Re + NematicOP.Im*NematicOP.Im) << endl;
  File.close();
  cout << "C3 OP = " << sqrt(NematicOP.Re*NematicOP.Re + NematicOP.Im*NematicOP.Im) << endl;

  
  Complex* SpinSpinCorrelation = new Complex [3];
  for (int i = 0; i < 3; ++i)
  {
    SpinSpinCorrelation[i] = 0.0;
  }
  int* TmpIndex1 = new int[3];
  for (int i = 0; i < XPeriodicity; ++i)
  {
    for (int j = 0; j < YPeriodicity; ++j)
    {
      for (int l = 0; l < 3; ++l)
	SpinSpinCorrelation[l] = 0.0;
      for (int nx = 0; nx < XPeriodicity; ++nx)
      {
	for (int ny = 0; ny < YPeriodicity; ++ny)
	{
	  int TmpIndex = GetLinearizedIndex(nx, ny, 0, XPeriodicity, YPeriodicity);
	  TmpIndex1[0] = GetLinearizedIndex(nx + i, ny + j, 0, XPeriodicity, YPeriodicity);
	  TmpIndex1[1] = GetLinearizedIndex(nx + i, ny + j, 1, XPeriodicity, YPeriodicity);
	  TmpIndex1[2] = GetLinearizedIndex(nx + i, ny + j, 2, XPeriodicity, YPeriodicity);
	 
	  for (int l = 0; l < 3; ++l)
	  {
	    Operator = new SpinWith2DTranslationSpinSpinCorrelationOperator(Space, XMomentum, XPeriodicity, YMomentum, YPeriodicity, TmpIndex, TmpIndex1[l]);
	    OperatorMatrixElementOperation Operation(Operator, State, State, State.GetVectorDimension());
	    Operation.ApplyOperation(Architecture.GetArchitecture());
	    SpinSpinCorrelation[l] += Operation.GetScalar();
// 	cout << NeighborSpinSpinCorrelation[l] << " " ;
	    delete Operator;
	  }
//       cout << endl;
	}
      }
      for (int l = 0; l < 3; ++l)
	File1 << i << " " << j << " " << l << " " << (SpinSpinCorrelation[l].Re/(XPeriodicity * YPeriodicity)) << endl;
    }
  }
  delete[] TmpIndex1;
  delete[] SpinSpinCorrelation;
  
  if (BondFlag == true)
  {
    OutputFileName = ReplaceExtensionToFileName(Manager.GetString("input-state"), "vec", "bondbond.dat");
    if (OutputFileName == 0)
    {
      cout << "no vec extension was find in " << Manager.GetString("input-state") << " file name" << endl;
      return 0;
    }
    File.open(OutputFileName, ios::binary | ios::out);
    File.precision(14);
    cout.precision(14);
  
    File << "# l_0 i j l_{i,j} <S_{0,a} S_{0,b} S_{i,j}S_{i,j,l}}> <S_{0,b} S_{0,c} S_{i,j}S_{i,j,l}}> <S_{0,a} S_{0,c} S_{i,j}S_{i,j,l}}>" << endl;
    int** TmpIndex2 = new int*[6];
    for (int l = 0; l < 6; ++l)
      TmpIndex2[l] = new int[2];
    int* TmpIndex1b = new int[2];
    Complex** BondBondCorrelations = new Complex*[6];
    for (int l = 0; l < 6; ++l)
      BondBondCorrelations[l] = new Complex[2];
    for (int i = MinX; i <= MaxX; ++i)
    {
      for (int j = MinY; j <= MaxY; ++j)
      {
	for (int l = 0; l < 6; ++l)
	  for (int m = 0; m < 2; ++m)
	    BondBondCorrelations[l][m] = 0.0;
	for (int nx = 0; nx < XPeriodicity; ++nx)
	{
	  for (int ny = 0; ny < YPeriodicity; ++ny)
	  {
	    int TmpIndex1a = GetLinearizedIndex(nx, ny, 0, XPeriodicity, YPeriodicity);
	    TmpIndex1b[0] = GetLinearizedIndex(nx, ny, 1, XPeriodicity, YPeriodicity);
	    TmpIndex1b[1] = GetLinearizedIndex(nx - 1, ny, 1, XPeriodicity, YPeriodicity);
	  
	    TmpIndex2[0][0] = GetLinearizedIndex(nx + i, ny + j, 0, XPeriodicity, YPeriodicity);
	    TmpIndex2[1][0] = GetLinearizedIndex(nx + i, ny + j, 0, XPeriodicity, YPeriodicity);
	    TmpIndex2[2][0] = GetLinearizedIndex(nx + i, ny + j, 1, XPeriodicity, YPeriodicity);
	    TmpIndex2[3][0] = GetLinearizedIndex(nx + i - Offset, ny + j - 1, 2, XPeriodicity, YPeriodicity);
	    TmpIndex2[4][0] = GetLinearizedIndex(nx + i, ny + j, 0, XPeriodicity, YPeriodicity);
	    TmpIndex2[5][0] = GetLinearizedIndex(nx + i, ny + j, 0, XPeriodicity, YPeriodicity);
	    
	    TmpIndex2[0][1] = GetLinearizedIndex(nx + i, ny + j, 1, XPeriodicity, YPeriodicity);
	    TmpIndex2[1][1] = GetLinearizedIndex(nx + i - 1, ny + j, 1, XPeriodicity, YPeriodicity);
	    TmpIndex2[2][1] = GetLinearizedIndex(nx + i, ny + j, 2, XPeriodicity, YPeriodicity);
	    TmpIndex2[3][1] = GetLinearizedIndex(nx + i - 1, ny + j, 1, XPeriodicity, YPeriodicity);
	    TmpIndex2[4][1] = GetLinearizedIndex(nx + i, ny + j, 2, XPeriodicity, YPeriodicity);
	    TmpIndex2[5][1] = GetLinearizedIndex(nx + i - Offset, ny + j - 1, 2, XPeriodicity, YPeriodicity);
	    
	     
	    for (int l = 0; l < 6; ++l)
	    {
	      for (int m = 0; m < 2; ++m)
	      {
		Operator = new SpinWith2DTranslationBondBondCorrelationOperator(Space, XMomentum, XPeriodicity, YMomentum, YPeriodicity, TmpIndex1a, TmpIndex1b[m], TmpIndex2[l][0], TmpIndex2[l][1]);
		OperatorMatrixElementOperation Operation(Operator, State, State, State.GetVectorDimension());
		Operation.ApplyOperation(Architecture.GetArchitecture());
		BondBondCorrelations[l][m] += Operation.GetScalar();
		delete Operator;
	      }
	    }
	  }
	}
	for (int l = 0; l < 6; ++l)
	{
// 	  File << i << " " << j << " " << l << " " << ((BondBondCorrelations[l][0].Re )/ (XPeriodicity * YPeriodicity) ) << " " << ((BondBondCorrelations[l][1].Re )/ (XPeriodicity * YPeriodicity)) << endl;
	  cout << (BondBondCorrelations[l][0].Re/ (XPeriodicity * YPeriodicity)) << " " << (NeighborSpinSpinCorrelation[0].Re) << " " << (NeighborSpinSpinCorrelation[l].Re) << " " << (NeighborSpinSpinCorrelation[0].Re * NeighborSpinSpinCorrelation[l].Re) << endl;
	  File << i << " " << j << " " << l << " " << (BondBondCorrelations[l][0].Re/ (XPeriodicity * YPeriodicity) - NeighborSpinSpinCorrelation[0].Re * NeighborSpinSpinCorrelation[l].Re) << " " << (BondBondCorrelations[l][1].Re/ (XPeriodicity * YPeriodicity)  - NeighborSpinSpinCorrelation[1].Re * NeighborSpinSpinCorrelation[l].Re) << endl;
	}
      }
    }

    File.close();
    
  for (int l = 0; l < 6; ++l)
    delete[] BondBondCorrelations[l];
  delete[] BondBondCorrelations;
  }
  
  
  return 0;
}


// get a linearized position index from the 2d coordinates
//
// xPosition = position along the x direction
// yPosition = position along the y direction
// return value = linearized index

int GetLinearizedIndex(int xPosition, int yPosition, int atomicIndex , int nbrSpinX, int nbrSpinY)
{
  if (xPosition < 0)
    xPosition += nbrSpinX;
  if (xPosition >= nbrSpinX)
    xPosition -= nbrSpinX;
  if (yPosition < 0)
    yPosition += nbrSpinY;
  if (yPosition >= nbrSpinY)
    yPosition -= nbrSpinY;
  return (3 * ((xPosition * nbrSpinY) + yPosition) + atomicIndex);
}