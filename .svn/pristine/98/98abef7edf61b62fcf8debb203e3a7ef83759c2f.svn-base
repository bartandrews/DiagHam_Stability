#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/FermionOnSphereMPSWrapper.h"
#include "HilbertSpace/FermionOnCylinderMPSWrapper.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "Tools/FQHEMPS/FQHEMPSMatrixManager.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Hamiltonian/TensorProductSparseMatrixHamiltonian.h"
#include "Hamiltonian/TensorProductSparseMatrixSelectedBlockHamiltonian.h"

#include "LanczosAlgorithm/BasicArnoldiAlgorithm.h"
#include "LanczosAlgorithm/BasicArnoldiAlgorithmWithDiskStorage.h"

#include "Matrix/SparseRealMatrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "MainTask/FQHEMPSEMatrixMainTask.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Options/Options.h"

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
  cout.precision(14); 
  
  OptionManager Manager ("FQHETorusMPSOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  FQHEMPSMatrixManager MPSMatrixManager (true, true);

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  OptionGroup* ArnoldiGroup  = new OptionGroup ("Arnoldi options");
  Architecture.AddOptionGroup(&Manager);
  Manager += ArnoldiGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-fluxquanta", "set the total number of flux quanta and deduce the number of particles", 0);
  (*SystemGroup) += new SingleIntegerOption ('\n', "left-topologicalsector", "set the topological sector of the left state", 0);
  (*SystemGroup) += new SingleIntegerOption ('\n', "right-topologicalsector", "set the topological sector of the right state", 0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");  

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusMPSOverlap -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFluxQuanta = 0;
  if (Manager.GetInteger("nbr-fluxquanta") <= 0)
    {
      cout << "invalid number of flux quanta" << endl;
      return -1;
    }
  NbrFluxQuanta = Manager.GetInteger("nbr-fluxquanta");

  AbstractFQHEMPSMatrix* LeftMPSMatrix = MPSMatrixManager.GetLeftMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 
  AbstractFQHEMPSMatrix* RightMPSMatrix = MPSMatrixManager.GetRightMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 
  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  int NbrBMatrices = LeftMPSMatrix->GetNbrMatrices();
  cout << "handling " << NbrBMatrices << " B matrices" << endl;
  int NbrParticles = LeftMPSMatrix->GetMatrixNaturalNbrParticles(NbrFluxQuanta, true);
  cout << "Nbr of particles = " << NbrParticles << " " << ", Nbr of flux quanta=" << NbrFluxQuanta << endl;
  SparseRealMatrix* LeftBMatrices = LeftMPSMatrix->GetMatrices();
  SparseRealMatrix* ConjugateLeftBMatrices = new SparseRealMatrix[NbrBMatrices];
  SparseRealMatrix* RightBMatrices = RightMPSMatrix->GetMatrices();
  SparseRealMatrix* ConjugateRightBMatrices = new SparseRealMatrix[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      ConjugateLeftBMatrices[i] = LeftBMatrices[i];//.Transpose();
      ConjugateRightBMatrices[i] = RightBMatrices[i];//.Transpose();
    } 
  SparseRealMatrix LeftStringMatrix;
  SparseRealMatrix RightStringMatrix;
  if (Manager.GetBoolean("boson") == true)
    {
      LeftStringMatrix = LeftMPSMatrix->GetTorusStringMatrix(0);
      RightStringMatrix = RightMPSMatrix->GetTorusStringMatrix(0);
    }
  else
    {
      LeftStringMatrix = LeftMPSMatrix->GetTorusStringMatrix(NbrParticles);
      RightStringMatrix = RightMPSMatrix->GetTorusStringMatrix(NbrParticles);
    }


  cout << "B matrix size = " << LeftBMatrices[0].GetNbrRow() << "x" << LeftBMatrices[0].GetNbrColumn() << endl;
  
  int NbrMPSLeftSumIndices = 0;
  int* MPSLeftSumIndices = LeftMPSMatrix->GetTopologicalSectorIndices(Manager.GetInteger("left-topologicalsector"), NbrMPSLeftSumIndices);
  int NbrMPSRightSumIndices = 0;
  int* MPSRightSumIndices = RightMPSMatrix->GetTopologicalSectorIndices(Manager.GetInteger("right-topologicalsector"), NbrMPSRightSumIndices);

  SparseRealMatrix TransferMatrixLeftLeft = TensorProduct(LeftBMatrices[0], ConjugateLeftBMatrices[0]) + TensorProduct(LeftBMatrices[1], ConjugateLeftBMatrices[1]);
  SparseRealMatrix TransferMatrixRightRight = TensorProduct(RightBMatrices[0], ConjugateRightBMatrices[0]) + TensorProduct(RightBMatrices[1], ConjugateRightBMatrices[1]);
  SparseRealMatrix TransferMatrixLeftRight = TensorProduct(LeftBMatrices[0], ConjugateRightBMatrices[0]) + TensorProduct(LeftBMatrices[1], ConjugateRightBMatrices[1]);
  
//   SparseRealMatrix LeftNormMatrix = TensorProduct(LeftStringMatrix, LeftStringMatrix);
//   SparseRealMatrix RightNormMatrix = TensorProduct(RightStringMatrix, RightStringMatrix);
//   SparseRealMatrix ScalarProductMatrix = TensorProduct(LeftStringMatrix, RightStringMatrix);

  int TmpNbrRow = LeftBMatrices[0].GetNbrRow() * LeftBMatrices[0].GetNbrRow();
  int* TmpNbrElementPerRow = new int [TmpNbrRow];
  for (int i = 0; i < TmpNbrRow; ++i)
    TmpNbrElementPerRow[i] = 0;
  for (int i = 0; i < NbrMPSLeftSumIndices; ++i)
    {
      for (int j = 0; j < NbrMPSLeftSumIndices; ++j)
	{
	  TmpNbrElementPerRow[MPSLeftSumIndices[i] * LeftBMatrices[0].GetNbrRow() + MPSLeftSumIndices[j]]++;
	}
    }
  SparseRealMatrix LeftNormMatrix (TmpNbrRow, LeftBMatrices[0].GetNbrColumn() * LeftBMatrices[0].GetNbrColumn(), TmpNbrElementPerRow);
  double TmpElement = 1.0;
  int TmpIndex;
  for (int i = 0; i < NbrMPSLeftSumIndices; ++i)
    {
      for (int j = 0; j < NbrMPSLeftSumIndices; ++j)
	{
	  TmpIndex = MPSLeftSumIndices[i] * LeftBMatrices[0].GetNbrRow() + MPSLeftSumIndices[j];
	  LeftNormMatrix.SetMatrixElement(TmpIndex, TmpIndex, TmpElement);
	}
    }  
  delete[] TmpNbrElementPerRow;

  TmpNbrRow = RightBMatrices[0].GetNbrRow() * RightBMatrices[0].GetNbrRow();
  TmpNbrElementPerRow = new int [TmpNbrRow];
  for (int i = 0; i < TmpNbrRow; ++i)
    TmpNbrElementPerRow[i] = 0;
  for (int i = 0; i < NbrMPSRightSumIndices; ++i)
    {
      for (int j = 0; j < NbrMPSRightSumIndices; ++j)
	{
	  TmpNbrElementPerRow[MPSRightSumIndices[i] * RightBMatrices[0].GetNbrRow() + MPSRightSumIndices[j]]++;
	}
    }
  SparseRealMatrix RightNormMatrix (TmpNbrRow, RightBMatrices[0].GetNbrColumn() * RightBMatrices[0].GetNbrColumn(), TmpNbrElementPerRow);
  for (int i = 0; i < NbrMPSRightSumIndices; ++i)
    {
      for (int j = 0; j < NbrMPSRightSumIndices; ++j)
	{
	  TmpIndex = MPSRightSumIndices[i] * RightBMatrices[0].GetNbrRow() + MPSRightSumIndices[j];
	  RightNormMatrix.SetMatrixElement(TmpIndex, TmpIndex, TmpElement);
	}
    }  
  delete[] TmpNbrElementPerRow;

  TmpNbrRow = LeftBMatrices[0].GetNbrRow() * RightBMatrices[0].GetNbrRow();
  TmpNbrElementPerRow = new int [TmpNbrRow];
  for (int i = 0; i < TmpNbrRow; ++i)
    TmpNbrElementPerRow[i] = 0;
  for (int i = 0; i < NbrMPSLeftSumIndices; ++i)
    {
      for (int j = 0; j < NbrMPSRightSumIndices; ++j)
	{
	  TmpNbrElementPerRow[MPSLeftSumIndices[i] * RightBMatrices[0].GetNbrRow() + MPSRightSumIndices[j]]++;
	}
    }
  SparseRealMatrix ScalarProductMatrix (TmpNbrRow, LeftBMatrices[0].GetNbrColumn() * RightBMatrices[0].GetNbrColumn(), TmpNbrElementPerRow);
  for (int i = 0; i < NbrMPSLeftSumIndices; ++i)
    {
      for (int j = 0; j < NbrMPSRightSumIndices; ++j)
	{
	  TmpIndex = MPSLeftSumIndices[i] * RightBMatrices[0].GetNbrRow() + MPSRightSumIndices[j];
	  ScalarProductMatrix.SetMatrixElement(TmpIndex, TmpIndex, TmpElement);
	}
    }  
  delete[] TmpNbrElementPerRow;

  LeftNormMatrix.Multiply(TensorProduct(LeftStringMatrix, LeftStringMatrix));
  RightNormMatrix.Multiply(TensorProduct(RightStringMatrix, RightStringMatrix));
  ScalarProductMatrix.Multiply(TensorProduct(LeftStringMatrix, RightStringMatrix));

  for (int i = 0; i < NbrFluxQuanta; ++i)
    {
      LeftNormMatrix.Multiply(TransferMatrixLeftLeft);
      RightNormMatrix.Multiply(TransferMatrixRightRight);
      ScalarProductMatrix.Multiply(TransferMatrixLeftRight);
    }

//   LeftNormMatrix.PrintNonZero(cout) << endl;
//   RightNormMatrix.PrintNonZero(cout) << endl;

//  double TmpElement;
  double LeftNorm = 0.0;
  for (int i = 0; i < NbrMPSLeftSumIndices; ++i)
    {
      for (int j = 0; j < NbrMPSLeftSumIndices; ++j)
	{
	  TmpIndex = MPSLeftSumIndices[i] * LeftBMatrices[0].GetNbrRow() + MPSLeftSumIndices[j];
	  LeftNormMatrix.GetMatrixElement(TmpIndex, TmpIndex, TmpElement);
	  LeftNorm += TmpElement;	  
	}
   }

  double RightNorm = 0.0;
  for (int i = 0; i < NbrMPSRightSumIndices; ++i)
    {
      for (int j = 0; j < NbrMPSRightSumIndices; ++j)
	{
	  TmpIndex = MPSRightSumIndices[i] * RightBMatrices[0].GetNbrRow() + MPSRightSumIndices[j];
	  RightNormMatrix.GetMatrixElement(TmpIndex, TmpIndex, TmpElement);
	  RightNorm += TmpElement;	  
	}
   }

  double ScalarProduct = 0.0;
  for (int i = 0; i < NbrMPSLeftSumIndices; ++i)
    {
      for (int j = 0; j < NbrMPSRightSumIndices; ++j)
	{
	  TmpIndex = MPSLeftSumIndices[i] * RightBMatrices[0].GetNbrRow() + MPSRightSumIndices[j];
	  ScalarProductMatrix.GetMatrixElement(TmpIndex, TmpIndex, TmpElement);
	  ScalarProduct += TmpElement;	  
	}
   }

  cout << "left square norm = " << LeftNorm << endl;
  cout << "right square norm = " << RightNorm << endl;
  cout << "unnormalized scalar product = " << ScalarProduct << endl;
  cout << "scalar product = " << (ScalarProduct / sqrt (LeftNorm * RightNorm)) << endl;
  
  return 0;
}

