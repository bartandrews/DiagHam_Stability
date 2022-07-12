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
#include "GeneralTools/MultiColumnASCIIFile.h"

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
  
  OptionManager Manager ("FQHECylinderMPSGeometricEntanglement" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  
  ArchitectureManager Architecture;
   FQHEMPSMatrixManager MPSMatrixManager (false, false);

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
  (*SystemGroup) += new BooleanOption  ('\n',"use-padding","use-padding");
  (*SystemGroup) += new SingleStringOption ('\n', "use-productstate", "use a prefined product state defined through a single column formatted ASCII file, instead of finding the optimal one");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");  


  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECylinderMPSGeometricEntanglement -h" << endl;
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

  int NbrParticles = MPSMatrix->GetMatrixNaturalNbrParticles(NbrFluxQuanta, Manager.GetBoolean("use-padding"));

  int NbrOrbitals = MPSMatrix->GetNbrOrbitals();
  int NbrBlock =   ((NbrFluxQuanta  + 1) / NbrOrbitals);
  cout << "MPSMatrix->GetMaximumOccupation()  = "<< MPSMatrix->GetMaximumOccupation()<<endl;
  cout <<"NbrOrbitals = "<< NbrOrbitals<<endl;
  cout <<" NbrBlock = " << NbrBlock<<endl;
  cout <<"NbrParticles = "<<NbrParticles<<endl;

  if ( ( NbrBlock   % NbrParticles) != 0)
    {
      cout << "invalid number of flux quanta" << endl;
      return -1;
    }

  int NbrBMatrices = MPSMatrix->GetNbrMatrices();

  int NbrStatesPerOrbital = MPSMatrix->GetMaximumOccupation() + 1;
  int NbrStatesPerBlock =  1;

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
  int MPSRowIndex = 0;
  int MPSColumnIndex = 0;
  MPSMatrix->GetMatrixBoundaryIndices(MPSRowIndex, MPSColumnIndex, Manager.GetBoolean("use-padding"));
  
  RealVector TmpVectorEMatrix (BMatrices[0].GetNbrRow() * BMatrices[0].GetNbrRow(), true);
  RealVector TmpVectorEMatrix2 (BMatrices[0].GetNbrRow() * BMatrices[0].GetNbrRow(), true);
  double* Coefficients = new double[NbrBMatrices];
  for(int i= 0; i < NbrBMatrices; i++)
    Coefficients[i] = 1.0;
  TmpVectorEMatrix[MPSRowIndex + (BMatrices[0].GetNbrRow() * MPSRowIndex)] = 1.0;
  TensorProductSparseMatrixHamiltonian* EHamiltonian = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, BMatrices, BMatrices, Coefficients,
												Architecture.GetArchitecture()); 
  for (int p = 0; p <  NbrBlock ;p++)
    {
      EHamiltonian->LowLevelMultiply(TmpVectorEMatrix, TmpVectorEMatrix2);
      RealVector TmpVectorEMatrix3 = TmpVectorEMatrix;
      TmpVectorEMatrix = TmpVectorEMatrix2;
      TmpVectorEMatrix2 = TmpVectorEMatrix3;
    }

  double Norm = TmpVectorEMatrix[MPSColumnIndex + (BMatrices[0].GetNbrColumn() * MPSColumnIndex)];
  double Normalisation = 1.0/ sqrt(Norm);


  cout << "B matrix size = " << BMatrices[0].GetNbrRow() << "x" << BMatrices[0].GetNbrColumn() << endl;
  unsigned long * ArrayPhysicalIndice = MPSMatrix->GetPhysicalIndices();
  
  SparseRealMatrix* FusedBMatrices = new SparseRealMatrix [DimensionPhysicalHilbertSpace];

  int TmpI;
  for(int i =0 ; i < DimensionPhysicalHilbertSpace; i++)
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
  if (Manager.GetString("use-productstate") == 0)
    {      
      for (int i = 0; i <DimensionPhysicalHilbertSpace; i++)
	{
	  CoefficientVector[i] = ((double) rand() / (RAND_MAX) - 0.5);
	}
      CoefficientVector.Normalize();

      cout <<"Start algorithm"<<endl;
      
      double PreviousOverlap = 1.0;
      double CurrentOverlap = 0.0;
      double Lambda = 1.0;
      double TotalCoefficientVectorSqrNorm =  pow(CoefficientVector.SqrNorm(), (double) NbrBlock);
      long Iteration = 0;
      double LogPrevious = -log (PreviousOverlap);
      double LogCurrent = -log (CurrentOverlap);
      RealVector TmpVector(FusedBMatrices[0].GetNbrRow(), true);
      RealVector TmpVector2(FusedBMatrices[0].GetNbrRow(), true);
      
      while(fabs(LogPrevious - LogCurrent ) > 1e-14)
	{
	  TmpVector.ClearVector();
	  SparseRealMatrix TmpMatrix = FusedBMatrices[0] * CoefficientVector[0];
	  SparseRealMatrix* TmpMatrix2 = new SparseRealMatrix[NbrBlock];
	  for (int i = 1; i < DimensionPhysicalHilbertSpace; i++)
	    {
	      TmpMatrix = TmpMatrix + (FusedBMatrices[i] * CoefficientVector[i]);
	    }

	  
	  TmpMatrix2[0].Copy(TmpMatrix);
	  
	  for (int i = 1; i< NbrBlock; i++)
	    {
	      TmpMatrix2[i] = Multiply(TmpMatrix2[i-1], TmpMatrix);
	    }
	  PreviousOverlap = CurrentOverlap / TotalCoefficientVectorSqrNorm;
	  LogPrevious = -log(PreviousOverlap) / ((double) NbrBlock) ;
	  cout << Iteration << " : " << PreviousOverlap << " " << (-log(PreviousOverlap) / ((double) NbrBlock)) << endl;
	  TmpMatrix2[NbrBlock-1].GetMatrixElement(MPSRowIndex,MPSColumnIndex, CurrentOverlap);
	  CurrentOverlap *= Normalisation;
	  
	  RealVector Derivative(DimensionPhysicalHilbertSpace, true);
	  for (int i = 0; i < DimensionPhysicalHilbertSpace; i++)
	    {
	      SparseRealMatrix TmpMatrix3 = Multiply(FusedBMatrices[i], TmpMatrix2[NbrBlock - 1 - 1]);
	      double Tmp;
	      TmpMatrix3.GetMatrixElement(MPSRowIndex,MPSColumnIndex ,Tmp);
	      Derivative[i] = Tmp;
	      for(int p = 1; p < (NbrBlock - 1);p++)
		{
		  SparseRealMatrix TmpMatrix4 = Multiply(Multiply(TmpMatrix2[p - 1], FusedBMatrices[i]), TmpMatrix2[(NbrBlock - 1) - p - 1]);
		  TmpMatrix4.GetMatrixElement(MPSRowIndex,MPSColumnIndex,Tmp);
		  Derivative[i] += Tmp;
		}
	      SparseRealMatrix TmpMatrix5 = Multiply(TmpMatrix2[NbrBlock - 1 - 1], FusedBMatrices[i]);
	      TmpMatrix5.GetMatrixElement(MPSRowIndex,MPSColumnIndex ,Tmp);
	      Derivative[i] += Tmp;
	    }
	  
	  Derivative *= 2.0 / CurrentOverlap;
	  for (int i = 0; i < DimensionPhysicalHilbertSpace; i++)
	    {
	      Derivative[i] -= 2.0 * NbrBlock * CoefficientVector[i] / CoefficientVector.SqrNorm();  
	    }
	  double Epsilon = 0.01;
	  for (int i = 0; i < DimensionPhysicalHilbertSpace; i++)
	    {
	      CoefficientVector[i] += Epsilon * Derivative[i];
	    }
	  
	  CoefficientVector.Normalize();
	  CurrentOverlap *= CurrentOverlap;
	  LogCurrent = -log (CurrentOverlap / TotalCoefficientVectorSqrNorm) / ((double) NbrBlock);
	  ++Iteration;
	  delete[] TmpMatrix2;
	}
      
      
      
      CoefficientVector.Normalize();
      for (int i = 0; i < DimensionPhysicalHilbertSpace;  i++)
	cout << i << " " <<CoefficientVector[i]<<endl;
      
    }
  else
    {
      MultiColumnASCIIFile ProductStateFile;
      if (ProductStateFile.Parse(Manager.GetString("use-productstate")) == false)
	{
	  ProductStateFile.DumpErrors(cout);
	  return -1;
	}
      if (ProductStateFile.GetNbrLines() != DimensionPhysicalHilbertSpace)
	{
	  cout << "error, " << Manager.GetString("use-productstate") << " does not have the proper number of components (" << DimensionPhysicalHilbertSpace << " vs " << ProductStateFile.GetNbrLines() << ")" << endl;
	  return -1;
	}
      CoefficientVector = RealVector(ProductStateFile.GetAsDoubleArray(0), DimensionPhysicalHilbertSpace);
      CoefficientVector.Normalize();

      for (int i = 0; i < DimensionPhysicalHilbertSpace; i++)
	{
	  FusedBMatrices[i] *= CoefficientVector[i];
	}

      RealVector TmpVector(FusedBMatrices[0].GetNbrRow(), true);
      RealVector TmpVector2(FusedBMatrices[0].GetNbrRow(), true);
      
      TmpVector[MPSRowIndex] = 1.0;

      
      for (int i = 0; i < NbrBlock; i++)
	{
	  FusedBMatrices[0].RightMultiply(TmpVector, TmpVector2);
	  for (int j = 1; j < DimensionPhysicalHilbertSpace; ++j)
	    {
	      FusedBMatrices[j].RightAddMultiply(TmpVector, TmpVector2);
	    }
	  RealVector TmpVector3 = TmpVector2;
	  TmpVector2 = TmpVector;
	  TmpVector = TmpVector3;
	}
      double Overlap = (TmpVector[MPSColumnIndex] * Normalisation);
      Overlap *= Overlap;
      cout << "Geometrical entropy = " << (-log(Overlap)  / ((double) NbrBlock)) << endl;
      cout << "Overlap = " << Overlap << endl;
    }


  
  return 0;
}

