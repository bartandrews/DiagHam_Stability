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

  OptionManager Manager ("FQHETorusMPSGeometricEntanglement" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  FQHEMPSMatrixManager MPSMatrixManager (false, true);

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  OptionGroup* ArnoldiGroup  = new OptionGroup ("Arnoldi options");
  Architecture.AddOptionGroup(&Manager);
  Manager += ArnoldiGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-fluxquanta", "set the total number of flux quanta and deduce the number of particles", 0);
  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-fusedblock", "number of blocks to be considered simultaneously", 1);
  (*SystemGroup) += new SingleIntegerOption ('\n', "topologicalsector", "set the topological sector of state", 0);
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
  int SizeBlock =  Manager.GetInteger("nbr-fusedblock");

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 

  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  int NbrParticles = MPSMatrix->GetMatrixNaturalNbrParticles(NbrFluxQuanta, true);
  int NbrBlock =  NbrParticles/SizeBlock;
  if ( (NbrParticles  % NbrBlock) != 0)
    {
      cout << "invalid number of flux quanta" << endl;
      return -1;
    }

  int NbrBMatrices = MPSMatrix->GetNbrMatrices();

  int NbrOrbitals = MPSMatrix->GetNbrOrbitals();
  int NbrStatesPerOrbital = MPSMatrix->GetMaximumOccupation() + 1;
  int NbrStatesPerBlock =  1;
  cout << "MPSMatrix->GetMaximumOccupation()  = "<< MPSMatrix->GetMaximumOccupation()<<endl;
  cout <<"NbrOrbitals = "<< NbrOrbitals<<endl;
  int DimensionPhysicalHilbertSpace = 1;
  
  for (int i = 0; i < NbrOrbitals * SizeBlock; i++)
     DimensionPhysicalHilbertSpace *= NbrStatesPerOrbital;
  
  for (int i = 0; i < NbrOrbitals ; i++)
    NbrStatesPerBlock *= NbrStatesPerOrbital;
  cout << "handling " << NbrBMatrices << " B matrices" << endl;
  cout << "Nbr of particles = " << NbrParticles << " " << ", Nbr of flux quanta=" << NbrFluxQuanta << endl;
  cout <<"Size block = "<< SizeBlock<< " , Nbr blocks = " << NbrBlock<<endl;
  cout <<"Dimension of the physical Hilbert space = " << DimensionPhysicalHilbertSpace <<endl;

  SparseRealMatrix* BMatrices = MPSMatrix->GetMatrices();
  SparseRealMatrix* ConjugateBMatrices = new SparseRealMatrix[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      ConjugateBMatrices[i] = BMatrices[i];//.Transpose();
    } 
  SparseRealMatrix StringMatrix;
  if (Manager.GetBoolean("boson") == true)
    {
      StringMatrix = MPSMatrix->GetTorusStringMatrix(0);
    }
  else
    {
      StringMatrix = MPSMatrix->GetTorusStringMatrix(NbrParticles);
    }
  

  
  int NbrMPSSumIndices = 0;
  int* MPSSumIndices = MPSMatrix->GetTopologicalSectorIndices(Manager.GetInteger("topologicalsector"), NbrMPSSumIndices);
  SparseRealMatrix TransferMatrix = TensorProduct(BMatrices[0], ConjugateBMatrices[0]);
  for(int i= 1; i < NbrStatesPerOrbital; i++)
    TransferMatrix =  TransferMatrix + TensorProduct(BMatrices[i], ConjugateBMatrices[i]);

  int TmpNbrRow = BMatrices[0].GetNbrRow() * BMatrices[0].GetNbrRow();
  int* TmpNbrElementPerRow = new int [TmpNbrRow];
  for (int i = 0; i < TmpNbrRow; ++i)
    TmpNbrElementPerRow[i] = 0;
  for (int i = 0; i < NbrMPSSumIndices; ++i)
    {
      for (int j = 0; j < NbrMPSSumIndices; ++j)
	{
	  TmpNbrElementPerRow[MPSSumIndices[i] * BMatrices[0].GetNbrRow() + MPSSumIndices[j]]++;
	}
    }
  SparseRealMatrix NormMatrix (TmpNbrRow, BMatrices[0].GetNbrColumn() * BMatrices[0].GetNbrColumn(), TmpNbrElementPerRow);
  double TmpElement = 1.0;
  int TmpIndex;
  for (int i = 0; i < NbrMPSSumIndices; ++i)
    {
      for (int j = 0; j < NbrMPSSumIndices; ++j)
	{
	  TmpIndex = MPSSumIndices[i] * BMatrices[0].GetNbrRow() + MPSSumIndices[j];
	  NormMatrix.SetMatrixElement(TmpIndex, TmpIndex, TmpElement);
	}
    }  
  delete[] TmpNbrElementPerRow;
  
//   cout << "NbrMPSSumIndices = " << NbrMPSSumIndices << endl;
//   for (int i = 0; i < NbrMPSSumIndices; ++i)
//     {
//       cout << i << " " << MPSSumIndices[i] << endl;
//     }
  
//  NormMatrix.PrintNonZero(cout) << endl;
  NormMatrix.Multiply(TensorProduct(StringMatrix, StringMatrix));
  for (int i = 0; i < (NbrFluxQuanta / NbrOrbitals); ++i)
    {
      NormMatrix.Multiply(TransferMatrix);
    }
//  TransferMatrix.PrintNonZero(cout) << endl;
//  NormMatrix.PrintNonZero(cout) << endl;
  
  double Norm = 0.0;
  for (int i = 0; i < NbrMPSSumIndices; ++i)
    {
      for (int j = 0; j < NbrMPSSumIndices; ++j)
	{
	  TmpIndex = MPSSumIndices[i] * BMatrices[0].GetNbrRow() + MPSSumIndices[j];
	  NormMatrix.GetMatrixElement(TmpIndex, TmpIndex, TmpElement);
	  Norm += TmpElement;	  
	}
    }
  
  double Normalisation = 1.0/sqrt(Norm);
  
  cout << "B matrix size = " << BMatrices[0].GetNbrRow() << "x" << BMatrices[0].GetNbrColumn() << endl;
  cout <<"Norm = " <<Norm<<endl;
  unsigned long * ArrayPhysicalIndice = MPSMatrix->GetPhysicalIndices();
  
  SparseRealMatrix* FusedBMatrices = new SparseRealMatrix [DimensionPhysicalHilbertSpace];
  
  int TmpI;
  for(int i =0 ; i <DimensionPhysicalHilbertSpace; i++)
    { 
      TmpI = i;			
      int Index = SearchInArray( (unsigned long)( TmpI %   NbrStatesPerBlock) , ArrayPhysicalIndice,  NbrBMatrices);
      if (Index <0)
	{
	  FusedBMatrices[i] = SparseRealMatrix(BMatrices[0].GetNbrRow(),BMatrices[0].GetNbrColumn());
	}
      else
	{
	  FusedBMatrices[i].Copy(BMatrices[Index]);
	}
      TmpI /= NbrStatesPerBlock;
      for(int p = 1; p < SizeBlock ; p++)
	{
	  int Index = SearchInArray( (unsigned long)( TmpI %   NbrStatesPerBlock) , ArrayPhysicalIndice,  NbrBMatrices);
	  if (Index <0)
	    {
	      FusedBMatrices[i].ClearMatrix ();
	    }
	  else
	    {
	      FusedBMatrices[i].Multiply(BMatrices[Index]);
	    }
	  TmpI /=NbrStatesPerBlock;
	}
    }
  
  
  RealVector CoefficientVector (DimensionPhysicalHilbertSpace,true);
  
  for (int i = 0; i <DimensionPhysicalHilbertSpace; i++)
    {
      CoefficientVector[i] = ((double) rand() / (RAND_MAX) - 0.5);
    }
  CoefficientVector.Normalize();
  
  double lamba = 0.0;
  double Newlamba = 100.0;
  
  cout <<"Start algorithm"<<endl;
  while ( fabs(lamba - Newlamba) > 1e-8  ) 
    {
      SparseRealMatrix TmpMatrix =  FusedBMatrices[0] * CoefficientVector[0];
      for (int i = 1; i<DimensionPhysicalHilbertSpace; i++)
	{
	  TmpMatrix = TmpMatrix + (FusedBMatrices[i] * CoefficientVector[i]);
	}
      
      SparseRealMatrix TmpMatrix2;
      TmpMatrix2.Copy(TmpMatrix);
      for (int i = 1; i< NbrBlock - 1; i++)
	{
	  TmpMatrix2.Multiply(TmpMatrix);
	}
      
      for (int i = 0; i < DimensionPhysicalHilbertSpace; i++)
	{
	  SparseRealMatrix TmpMatrix3 = Multiply(FusedBMatrices[i],TmpMatrix2);
	  CoefficientVector[i] = TmpMatrix3.PartialTr(MPSSumIndices, NbrMPSSumIndices) * Normalisation;
	}
      lamba = Newlamba ;
      Newlamba = CoefficientVector.Norm();
      CoefficientVector.Normalize();
      
      cout <<Newlamba << " "<<lamba<<endl;
    }
  
  
  SparseRealMatrix TmpMatrix =  FusedBMatrices[0] * CoefficientVector[0];
  
  for (int i = 1; i<DimensionPhysicalHilbertSpace; i++)
    {
      TmpMatrix = TmpMatrix + (FusedBMatrices[i] * CoefficientVector[i]);
    }
  

  SparseRealMatrix TmpMatrix2;
  TmpMatrix2.Copy(TmpMatrix);
  for (int i = 1; i < NbrBlock ; i++)
    {
      TmpMatrix2.Multiply(TmpMatrix);
    }
  
  double FinalResult = TmpMatrix2.PartialTr(MPSSumIndices, NbrMPSSumIndices) * Normalisation;
  
  cout <<" FinalResult = "<<FinalResult <<endl;
  
  for (int i = 0; i < DimensionPhysicalHilbertSpace;  i++)
    cout << i <<" " <<CoefficientVector[i]<<endl;
  
  return 0;
}

