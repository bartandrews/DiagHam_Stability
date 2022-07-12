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

#include "Tools/FQHEMPS/FQHEMPODensityOperator.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Hamiltonian/TensorProductSparseMatrixHamiltonian.h"
#include "Hamiltonian/TensorProductSparseMatrixSelectedBlockHamiltonian.h"
#include "Hamiltonian/TripleTensorProductSparseMatrixHamiltonian.h"

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

#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"

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


// compute the spectrum of a single egde
// 
// mPSMatrix = pointer tothe MPS matrix
// overlapMatrix= reference on the overlap matrix
// hamiltonianMatrix = reference on the hamiltonian matrix for the edge
// pSector = P sector that should be considered
// cFTSector = CFT sector that should be considered
// qSector = Q sector that should be considered
// eigenvalueError = relative error on the eigenvalues below which an eigenvalue is considered to be equal to zero
RealDiagonalMatrix FQHEMPSEvaluateSingleEdgeSpectrum(AbstractFQHEMPSMatrix* mPSMatrix, RealMatrix& overlapMatrix, RealMatrix& hamiltonianMatrix,
						     int pSector, int cFTSector, int qSector, double eigenvalueError);


int main(int argc, char** argv)
{
  cout.precision(14); 
  
  OptionManager Manager ("FQHECylinderMPSEdgeSpectrum" , "0.01");
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

  (*SystemGroup) += new BooleanOption ('\n', "infinite-cylinder", "evaluate the entanglement spectrum on the infinite cylinder");  
  (*SystemGroup) += new SingleStringOption  ('\n', "left-eigenstate", "file containing the transfer matrix left eigenstate");
  (*SystemGroup) += new SingleStringOption  ('\n', "right-eigenstate", "file containing the transfer matrix right eigenstate");
  (*SystemGroup) += new BooleanOption  ('\n',"use-padding","use-padding");
  (*SystemGroup) += new SingleStringOption ('\n', "left-interaction", "file describing the confining potential of the left hand side of the cylinder (optional)");
  (*SystemGroup) += new SingleStringOption ('\n', "right-interaction", "file describing the confining potential of the right hand side of the cylinder");
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");  


  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECylinderMPSEdgeSpectrum -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  double Error = 1e-13;
  int LeftNbrFluxQuanta = -1;
  double* LeftOneBodyPotential = 0;
  int RightNbrFluxQuanta = 0;
  double* RightOneBodyPotential = 0;
  if (Manager.GetString("right-interaction") == 0)
    {
      cout << "error, the right interaction description should be provided" << endl;
      return -1;
    }
  if (FQHESphereGetOneBodyPotentials(Manager.GetString("right-interaction"), RightNbrFluxQuanta, RightOneBodyPotential) == false)
    {
      return -1;
    }
  if (Manager.GetString("left-interaction") != 0)
    {
      if (FQHESphereGetOneBodyPotentials(Manager.GetString("left-interaction"), LeftNbrFluxQuanta, LeftOneBodyPotential) == false)
	{
	  return -1;
	}
    }
  int NbrFluxQuanta = (LeftNbrFluxQuanta + RightNbrFluxQuanta) + 1;

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 
  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }
  int NbrBMatrices = MPSMatrix->GetNbrMatrices();
  SparseRealMatrix* BMatrices = MPSMatrix->GetMatrices();
  double Normalisation = 1.0;

  if (Manager.GetBoolean("infinite-cylinder") == false)
    {
      double* TotalOneBodyPotentials = new double[NbrFluxQuanta + 1];
      if (Manager.GetString("left-interaction") != 0)
	{
	  for (int i = 0; i <= LeftNbrFluxQuanta; ++i)
	    {
	      TotalOneBodyPotentials[i] = LeftOneBodyPotential[i];
	    }
	}
      for (int i = 0; i <= RightNbrFluxQuanta; ++i)
	{
	  TotalOneBodyPotentials[i + LeftNbrFluxQuanta + 1] = RightOneBodyPotential[i];
	}

      cout << "Number of flux quanta = " << NbrFluxQuanta << endl;
      int MPSRowIndex = 0;
      int MPSColumnIndex = 0;
      MPSMatrix->GetMatrixBoundaryIndices(MPSRowIndex, MPSColumnIndex, Manager.GetBoolean("use-padding"));
      
      RealVector TmpVectorEMatrix (BMatrices[0].GetNbrRow() * BMatrices[0].GetNbrRow(), true);
      RealVector TmpVectorEMatrix2 (BMatrices[0].GetNbrRow() * BMatrices[0].GetNbrRow(), true);
      double* Coefficients = new double[NbrBMatrices];
      for(int i= 0; i < NbrBMatrices; i++)
	Coefficients[i] = 1.0;
      TmpVectorEMatrix[MPSRowIndex + (BMatrices[0].GetNbrRow() * MPSRowIndex)] = 1.0;

//       SparseRealMatrix* TmpMPOMatrices = new SparseRealMatrix [NbrBMatrices];
//       for (int i = 0; i < NbrBMatrices; ++i)
// 	{
// 	  TmpMPOMatrices[i] = SparseRealMatrix(1, 1);
// 	  TmpMPOMatrices[i].SetMatrixElement(0, 0, 1.0);
// 	}

      TensorProductSparseMatrixHamiltonian* TransferMatrix = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, BMatrices, BMatrices, Coefficients,
												      Architecture.GetArchitecture()); 
//       TripleTensorProductSparseMatrixHamiltonian* TransferMatrix = new TripleTensorProductSparseMatrixHamiltonian(NbrBMatrices, BMatrices, TmpMPOMatrices, BMatrices, Coefficients,
// 														  Architecture.GetArchitecture()); 
      for (int i = 0; i <= NbrFluxQuanta; ++i)
	{
	  TransferMatrix->LowLevelMultiply(TmpVectorEMatrix, TmpVectorEMatrix2);
//	  cout << TmpVectorEMatrix2 << "---------" << endl;
	  RealVector TmpVectorEMatrix3 = TmpVectorEMatrix;
	  TmpVectorEMatrix = TmpVectorEMatrix2;
	  TmpVectorEMatrix2 = TmpVectorEMatrix3;
	}      
      Normalisation = 1.0 / TmpVectorEMatrix[MPSColumnIndex + (BMatrices[0].GetNbrColumn() * MPSColumnIndex)];

      cout << "Normalisation = " << Normalisation << endl;
      
//       for (int i = 0; i < NbrBMatrices; ++i)
//  	{
//  	  BMatrices[i] = SparseRealMatrix(1, 1);
//  	  BMatrices[i].SetMatrixElement(0, 0, 1.0);
//  	}
//       MPSRowIndex = 0;
//       MPSColumnIndex = 0;

      FQHEMPODensityOperator TmpMPO (MPSMatrix->GetMaximumOccupation(), 1.0);
      int MPORowIndex = 0;
      int MPOColumnIndex = 0;
      TmpMPO.GetMatrixBoundaryIndices(MPORowIndex, MPOColumnIndex);
      SparseRealMatrix* TmpMPOMatrices = TmpMPO.GetMatrices();
      cout << "MPSRowIndex=" << MPSRowIndex << " MPORowIndex=" << MPORowIndex << " MPSColumnIndex=" << MPSColumnIndex << " MPOColumnIndex=" << MPOColumnIndex << endl;
      TmpVectorEMatrix = RealVector(BMatrices[0].GetNbrRow() * BMatrices[0].GetNbrRow() * TmpMPO.GetBondDimension(), true);
      TmpVectorEMatrix2 = RealVector(BMatrices[0].GetNbrRow() * BMatrices[0].GetNbrRow() * TmpMPO.GetBondDimension(), true);
      TmpVectorEMatrix[MPSColumnIndex + (((TmpMPOMatrices[0].GetNbrColumn() * MPSColumnIndex) + MPOColumnIndex) * BMatrices[0].GetNbrColumn())] = 1.0;
//       for (int i = 0; i < NbrBMatrices; ++i)
// 	{
// 	  cout << TmpMPOMatrices[i] << endl;
// 	}
      for (int i = 0; i <= NbrFluxQuanta; ++i)
	{
	  TmpMPO.SetPrefactor(TotalOneBodyPotentials[i]);
	  TripleTensorProductSparseMatrixHamiltonian* MPOTransferMatrix = new TripleTensorProductSparseMatrixHamiltonian(NbrBMatrices, BMatrices, TmpMPOMatrices, BMatrices, Coefficients,
															 Architecture.GetArchitecture()); 
	  MPOTransferMatrix->LowLevelMultiply(TmpVectorEMatrix, TmpVectorEMatrix2);
 	  RealVector TmpVectorEMatrix3 = TmpVectorEMatrix;
 	  TmpVectorEMatrix = TmpVectorEMatrix2;
 	  TmpVectorEMatrix2 = TmpVectorEMatrix3;	  
 	}
      double TmpNbrParticles = TmpVectorEMatrix[MPSRowIndex + (((TmpMPOMatrices[0].GetNbrRow() * MPSRowIndex) + MPORowIndex) * BMatrices[0].GetNbrRow())];
      cout << (TmpNbrParticles * Normalisation) << endl;
    }
  else
    {      
      ComplexVector TransferMatrixRightEigenstate;
      ComplexVector TransferMatrixLeftEigenstate;
      if (Manager.GetString("left-eigenstate") == 0)
	{
	  cout << "error, the transfer matrix left eigenstate should be provided" << endl;
	  return -1;
	}
      if (TransferMatrixLeftEigenstate.ReadVector(Manager.GetString("left-eigenstate")) == false)
	{
	  cout << "can't read " << Manager.GetString("left-eigenstate") << endl;
	  return 0;
	}            
      if (Manager.GetString("left-interaction") != 0)
	{
	  if (Manager.GetString("right-eigenstate") == 0)
	    { 
	      cout << "error, the transfer matrix right eigenstate should be provided" << endl;
	      return -1;
	    }
	  if (TransferMatrixRightEigenstate.ReadVector(Manager.GetString("right-eigenstate")) == false)
	    {
	      cout << "can't read " << Manager.GetString("right-eigenstate") << endl;
	      return 0;
	    }            
	}
      FQHEMPODensityOperator TmpMPO (MPSMatrix->GetMaximumOccupation(), 1.0);
      int MPORowIndex = 0;
      int MPOColumnIndex = 0;
      TmpMPO.GetMatrixBoundaryIndices(MPORowIndex, MPOColumnIndex);
      SparseRealMatrix* TmpMPOMatrices = TmpMPO.GetMatrices();
      double* Coefficients = new double[NbrBMatrices];
      for(int i= 0; i < NbrBMatrices; i++)
	Coefficients[i] = 1.0;
      TensorProductSparseMatrixHamiltonian* MPSTransferMatrix = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, BMatrices,BMatrices, Coefficients,
													 Architecture.GetArchitecture()); 
      TripleTensorProductSparseMatrixHamiltonian* MPSMPOTransferMatrix = new TripleTensorProductSparseMatrixHamiltonian(NbrBMatrices, BMatrices, TmpMPOMatrices, BMatrices, Coefficients,
															Architecture.GetArchitecture()); 



      ComplexVector TmpVector1 (MPSMatrix->GetBondDimension() * MPSMatrix->GetBondDimension(), true);
      ComplexVector TmpVector2 (MPSMatrix->GetBondDimension() * MPSMatrix->GetBondDimension(), true);
      MPSTransferMatrix->LowLevelMultiply(TransferMatrixLeftEigenstate, TmpVector1);
      MPSTransferMatrix->LowLevelMultiply(TmpVector1, TmpVector2);
      MPSTransferMatrix->LowLevelMultiply(TmpVector2, TmpVector1);
      double MPSNorm = pow(TmpVector1.Norm(), 1.0 / 3.0);
      MPSNorm = 1.0 / sqrt(MPSNorm);
      for(int i= 0; i < NbrBMatrices; i++)
	{
	  BMatrices[i] *= MPSNorm;
	}

      TmpVector1 = ComplexVector (MPSMatrix->GetBondDimension() * MPSMatrix->GetBondDimension() * TmpMPO.GetBondDimension(), true);
      TmpVector2 = ComplexVector (MPSMatrix->GetBondDimension() * MPSMatrix->GetBondDimension() * TmpMPO.GetBondDimension(), true);							
      int TmpLeftIndex;
      int TmpRightIndex;
      int TmpLinearizedIndex;
      for (int i = 0 ; i < TransferMatrixLeftEigenstate.GetVectorDimension(); ++i)
	{
	  MPSTransferMatrix->GetIndicesFromLinearizedIndex(i, TmpLeftIndex, TmpRightIndex);
	  TmpLinearizedIndex = MPSMPOTransferMatrix->GetLinearizedIndex(TmpLeftIndex, MPOColumnIndex, TmpRightIndex);
	  TmpVector1[TmpLinearizedIndex] = TransferMatrixLeftEigenstate[i];
	}
      if (Manager.GetString("left-interaction") == 0)
	{
	  for (int i = 0; i <= RightNbrFluxQuanta; ++i)
	    {
 	      TmpMPO.SetPrefactor(RightOneBodyPotential[i]);
 	      MPSMPOTransferMatrix->LowLevelMultiply(TmpVector1, TmpVector2);
 	      ComplexVector TmpVector3 = TmpVector1;
 	      TmpVector1 = TmpVector2;
 	      TmpVector2 = TmpVector3;	      
	    }
	  
	  TmpVector2 = ComplexVector(MPSMatrix->GetBondDimension() * MPSMatrix->GetBondDimension(), true);
	  for (int i = 0 ; i < TransferMatrixLeftEigenstate.GetVectorDimension(); ++i)
	    {
	      MPSTransferMatrix->GetIndicesFromLinearizedIndex(i, TmpLeftIndex, TmpRightIndex);
	      TmpLinearizedIndex = MPSMPOTransferMatrix->GetLinearizedIndex(TmpLeftIndex, MPORowIndex, TmpRightIndex);
	      TmpVector2[i] = TmpVector1[TmpLinearizedIndex];
	    }
	  cout << (TmpVector2 * TmpVector2) << endl;

	  int TmpDimension =  MPSMatrix->GetBondDimension();
	  RealMatrix OverlapMatrix(TmpDimension, TmpDimension, true);
	  RealMatrix HamiltonianMatrix(TmpDimension, TmpDimension, true);
 	  for (int i = 0; i < TmpDimension; ++i)
	    for (int j = 0; j < TmpDimension; ++j)
	      {
		OverlapMatrix.SetMatrixElement(i, j,  TransferMatrixLeftEigenstate[i * TmpDimension + j].Re); 
		HamiltonianMatrix.SetMatrixElement(i, j, TmpVector2[i * TmpDimension + j].Re); 
	      }
	  
	  double**** EdgeSpectrum = new double***[MPSMatrix->GetTruncationLevel() + 1];
	  int*** EdgeSpectrumDimension = new int**[MPSMatrix->GetTruncationLevel() + 1];
	  for (int CurrentPLevel = 0; CurrentPLevel <= MPSMatrix->GetTruncationLevel(); ++CurrentPLevel)
	    {
	      EdgeSpectrumDimension[CurrentPLevel] = new int*[MPSMatrix->GetNbrCFTSectors()];
	      EdgeSpectrum[CurrentPLevel] = new double**[MPSMatrix->GetNbrCFTSectors()];	      	      
	      for (int CurrentCFTSector = 0; CurrentCFTSector < MPSMatrix->GetNbrCFTSectors(); ++CurrentCFTSector)
		{
		  int LocalMinQValue;
		  int LocalMaxQValue;
		  MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, LocalMinQValue, LocalMaxQValue);
		  if (LocalMinQValue <=  LocalMaxQValue)
		    {
		      EdgeSpectrumDimension[CurrentPLevel][CurrentCFTSector] = new int[LocalMaxQValue - LocalMinQValue + 1];
		      EdgeSpectrum[CurrentPLevel][CurrentCFTSector] = new double*[LocalMaxQValue - LocalMinQValue + 1];	      
		    }
		  else
		    {
		      EdgeSpectrum[CurrentPLevel][CurrentCFTSector] = 0;
		    }
		  for (int LocalQValue =  LocalMinQValue; LocalQValue <= LocalMaxQValue; ++LocalQValue)
		    {
		      cout << "computing sector P=" << CurrentPLevel<< " CFT=" << CurrentCFTSector << " Q=" << LocalQValue << endl;
		      RealDiagonalMatrix TmpEdgeSpectrum = FQHEMPSEvaluateSingleEdgeSpectrum(MPSMatrix, OverlapMatrix, HamiltonianMatrix, CurrentPLevel, CurrentCFTSector, LocalQValue, Error);		      
		      if (TmpEdgeSpectrum.GetNbrRow() > 0)
			{
			  for (int i = 0 ; i < TmpEdgeSpectrum.GetNbrRow(); ++i)
			    {
			      if (TmpEdgeSpectrum[i] > 0.0)
				EdgeSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]++;
			    }
			  if (EdgeSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] > 0)
			    {
			      EdgeSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = new double[EdgeSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]];
			      EdgeSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = 0;
			      for (int i = 0 ; i < TmpEdgeSpectrum.GetNbrRow(); ++i)
				{
				  if (TmpEdgeSpectrum[i] > 0.0)
				    {
				      EdgeSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][EdgeSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]] = TmpEdgeSpectrum[i];
				      EdgeSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]++;
				    }
				}
			      SortArrayDownOrdering<double>(EdgeSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue],
							    EdgeSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]);
			    }
			}
		      else
			{
			  EdgeSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = 0;
			}		      
		    }
		}
	    }
	  for (int CurrentPLevel = 0; CurrentPLevel <= MPSMatrix->GetTruncationLevel(); ++CurrentPLevel)
	    {
	      for (int CurrentCFTSector = 0; CurrentCFTSector < MPSMatrix->GetNbrCFTSectors(); ++CurrentCFTSector)
		{
		  int LocalMinQValue;
		  int LocalMaxQValue;
		  MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, LocalMinQValue, LocalMaxQValue);
		  for (int LocalQValue =  LocalMinQValue; LocalQValue <= LocalMaxQValue; ++LocalQValue)
		    {
		      for (int i = 0; i < EdgeSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]; ++i)
			{
			  cout << CurrentCFTSector  << " " << LocalQValue << " " 
			       << CurrentPLevel << " "
			       <<  EdgeSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i]  << endl;
			}		      
		    }
		}
	    }

	}
    }
  
  return 0;
}

// compute the spectrum of a single egde
// 
// mPSMatrix = pointer tothe MPS matrix
// overlapMatrix= reference on the overlap matrix
// hamiltonianMatrix = reference on the hamiltonian matrix for the edge
// pSector = P sector that should be considered
// cFTSector = CFT sector that should be considered
// qSector = Q sector that should be considered
// eigenvalueError = relative error on the eigenvalues below which an eigenvalue is considered to be equal to zero
// eigenstateFileName = if non-zero, save the eigenstate of the reduced density matrix using eigenstateFileName as a prefix

RealDiagonalMatrix FQHEMPSEvaluateSingleEdgeSpectrum(AbstractFQHEMPSMatrix* mPSMatrix, RealMatrix& overlapMatrix, RealMatrix& hamiltonianMatrix,
						     int pSector, int cFTSector, int qSector, double eigenvalueError)
{
  RealMatrix OverlapMatrix = mPSMatrix->ExtractBlock(overlapMatrix, pSector, cFTSector,qSector, pSector, cFTSector, qSector);
  RealDiagonalMatrix TmpEnergySpectrum;
  if ((OverlapMatrix.GetNbrRow() > 0) && (OverlapMatrix.ComputeNbrNonZeroMatrixElements() > 0l))
    {
      cout << "diagonalizing the overlap" << endl;
      RealSymmetricMatrix SymOverlapMatrix (OverlapMatrix);
      RealMatrix TmpBasis(SymOverlapMatrix.GetNbrRow(), SymOverlapMatrix.GetNbrRow());
      TmpBasis.SetToIdentity();
      RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
      SymOverlapMatrix.LapackDiagonalize(TmpDiag, TmpBasis);
#else
      SymOverlapMatrix.Diagonalize(TmpDiag, TmpBasis);
#endif
      double LocalEigenvalueError = 0.0;
      for (int i = 0; i < TmpDiag.GetNbrColumn(); ++i)
	if (TmpDiag(i, i) > LocalEigenvalueError)
	  LocalEigenvalueError = TmpDiag(i, i);
      cout << "LocalEigenvalueError=" << LocalEigenvalueError << endl;
      LocalEigenvalueError = eigenvalueError;
      int NbrZeroEigenvalues = 0;
      for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
	{
	  if (TmpDiag(i, i) < LocalEigenvalueError)
	    {
	      ++NbrZeroEigenvalues;	    
	    }
	}
      cout << "nbr non zero eigenvalues = " << (TmpDiag.GetNbrRow() - NbrZeroEigenvalues) << " (full dim = " << TmpDiag.GetNbrRow() << ")" << endl;
      
      if (NbrZeroEigenvalues < SymOverlapMatrix.GetNbrRow())
	{
 	  RealMatrix TruncatedBasis (TmpDiag.GetNbrRow(), TmpDiag.GetNbrRow() -  NbrZeroEigenvalues, true);
 	  NbrZeroEigenvalues = 0;
 	  for (int i = 0; i < TmpBasis.GetNbrColumn(); ++i)
 	    {
 	      if (TmpDiag(i, i) > LocalEigenvalueError)
 		{
 		  TruncatedBasis[NbrZeroEigenvalues].Copy(TmpBasis[i]);
 		  TruncatedBasis[NbrZeroEigenvalues] *= sqrt(TmpDiag(i, i));
 		  ++NbrZeroEigenvalues;
 		}
 	    }

	  RealMatrix HamiltonianMatrix = mPSMatrix->ExtractBlock(hamiltonianMatrix, pSector, cFTSector,qSector, pSector, cFTSector, qSector);
	  RealSymmetricMatrix SymHamiltonianMatrix (HamiltonianMatrix);
	  RealSymmetricMatrix* TruncatedHamiltonianMatrix = (RealSymmetricMatrix*) SymHamiltonianMatrix.Conjugate(TruncatedBasis);
#ifdef __LAPACK__
	  TruncatedHamiltonianMatrix->LapackDiagonalize(TmpEnergySpectrum);
#else
	  TruncatedHamiltonianMatrix->Diagonalize(TmpEnergySpectrum);
#endif
	}
    }
  return TmpEnergySpectrum;
}
