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

#include "MainTask/FQHEMPSEMatrixMainTask.h"

#include "GeneralTools/ArrayTools.h"

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
  
  OptionManager Manager ("FQHESphereMPSEntanglementSpectrum" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  FQHEMPSMatrixManager MPSMatrixManager;

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  OptionGroup* ArnoldiGroup  = new OptionGroup ("Arnoldi options");
  Architecture.AddOptionGroup(&Manager);
  Manager += ArnoldiGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "file that describes the root configuration");
  (*SystemGroup) += new BooleanOption  ('\n', "use-padding", "root partitions use the extra zero padding");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "la", "number of orbitals in subsystem A", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "na", "number of particles in subsystem A", 0);
  (*SystemGroup) += new BooleanOption ('\n', "all-na", "print all charge sectors");
  (*SystemGroup) += new BooleanOption ('\n', "infinite-cylinder", "evaluate the entanglement spectrum on the infinite cylinder");
  (*SystemGroup) += new BooleanOption ('\n', "use-singlestate", "use a single real eigenstate of the E matrix  when evaluating the infinite entanglement spectrum");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "memory", "amount of memory that can used for precalculations (in Mb)", 500);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "ematrix-memory", "amount of memory that can used for precalculations of the E matrix (in Mb)", 500);
  (*OutputGroup) += new SingleStringOption  ('o', "output-file", "output file name");
  (*OutputGroup) += new BooleanOption ('\n', "suppress-output", "minimize the amount of output information");

  (*OutputGroup) += new BooleanOption ('n', "normalize-sphere", "express the MPS in the normalized sphere basis");
  (*ArnoldiGroup) += new SingleIntegerOption  ('\n', "full-diag", 
					       "maximum Hilbert space dimension for which full diagonalization is applied", 1000);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "disk", "enable disk storage for the Arnoldi algorithm", false);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "resume", "resume from disk datas", false);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Arnoldi iteration", false); 
  (*ArnoldiGroup) += new  SingleIntegerOption ('\n', "arnoldi-memory", "amount of memory when using the Arnoldi algorithm (in Mb)", 500); 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");  

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereMPSEntanglementSpectrum -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0; 
  int NbrFluxQuanta = 0;
  int TotalLz = 0;
  int* ReferenceState = 0;
  if ((Manager.GetBoolean("infinite-cylinder") == false) && (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false))
    return -1;

  int EntCut = Manager.GetInteger("la");
  int Na = Manager.GetInteger("na");

  bool MinimizeOutput = Manager.GetBoolean("suppress-output");

  bool CylinderFlag = Manager.GetBoolean("normalize-cylinder");

  int LandauLevel = 0;

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 
  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  int NbrBMatrices = MPSMatrix->GetNbrMatrices();
  SparseRealMatrix* BMatrices = MPSMatrix->GetMatrices();
  SparseRealMatrix* ConjugateBMatrices = new SparseRealMatrix[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    ConjugateBMatrices[i] = BMatrices[i].Transpose();

  cout << "B matrix size = " << BMatrices[0].GetNbrRow() << "x" << BMatrices[0].GetNbrColumn() << endl;
  
  int MPSRowIndex = 0;
  int MPSColumnIndex = 0;
  int NbrEigenstates = 0;
  int MinQ;
  int MaxQ;
  MPSMatrix->GetChargeIndexRange(0, MinQ, MaxQ);
  MPSMatrix->GetMatrixBoundaryIndices(MPSRowIndex, MPSColumnIndex, Manager.GetBoolean("use-padding"));
  NbrEigenstates = MPSMatrix->GetTransferMatrixLargestEigenvalueDegeneracy();

  ofstream File;
  File.precision(14);
  
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = new char [512];
      if (CylinderFlag == true)
	{
	  if (Manager.GetBoolean("infinite-cylinder"))
	    {
	      sprintf(TmpFileName, "fermions_infinite_cylinder_%s_perimeter_%f_plevel_%ld_n_0_2s_0_lz_0.0.full.ent", MPSMatrix->GetName(),
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"));
	    }
	  else
	    {
	      sprintf(TmpFileName, "fermions_cylinder_%s_plevel_%ld_n_%d_2s_%d_lz_%d.0.full.ent", MPSMatrix->GetName(),
		      Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz);
	    }
	}
      else
	{
	  if (Manager.GetBoolean("infinite-cylinder"))
	    {
	      sprintf(TmpFileName, "fermions_infinite_%s_plevel_%ld_n_0_2s_0_lz_0.0.full.ent", MPSMatrix->GetName(),
		      Manager.GetInteger("p-truncation"));
	    }
	  else
	    {
	      sprintf(TmpFileName, "fermions_%s_plevel_%ld_n_%d_2s_%d_lz_%d.0.full.ent", MPSMatrix->GetName(),
		      Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz);
	    }
	}
      File.open(TmpFileName, ios::binary | ios::out);     
   }


  if (Manager.GetBoolean("infinite-cylinder"))
    {
      double Error = 1e-13;

      SparseRealMatrix* SparseBMatrices = MPSMatrix->GetMatrices();
      SparseRealMatrix* SparseTransposeBMatrices = new SparseRealMatrix[NbrBMatrices];
      double* Coefficients = new double[NbrBMatrices];
      for (int i = 0; i < NbrBMatrices; ++i)
	{
	  Coefficients[i] = 1.0;
	  SparseTransposeBMatrices[i] = SparseBMatrices[i].Transpose();
	}
      
      int MinQValue = 0;
      int MaxQValue = 0;
      int TmpBMatrixDimension = SparseBMatrices[0].GetNbrRow();

      long EffectiveDimension = 0l;
      for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	{
	  MPSMatrix->GetChargeIndexRange(PLevel, MinQValue, MaxQValue);
	  for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
	    {
	      long Tmp = MPSMatrix->GetBondIndexRange(PLevel, QValue);
	      EffectiveDimension += Tmp * Tmp;
	    }
	}

      cout << "number of eigenvalues for the transfer matrix = " << NbrEigenstates << endl;
      cout << "computing effective transfer matrix indices " << endl;
      long** BlockIndexProductTable = new long* [TmpBMatrixDimension];
      int* BlockIndexProductTableNbrElements = new int [TmpBMatrixDimension];
      int* BlockIndexProductTableShift = new int [TmpBMatrixDimension];
      for (long i = 0; i < TmpBMatrixDimension; ++i)
	{
	  BlockIndexProductTableNbrElements[i] = 0;
	  BlockIndexProductTableShift[i] = -1;	  
	}

      long* EffectiveBlockIndices = new long [EffectiveDimension];
      EffectiveDimension = 0l;
      for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	{
	  MPSMatrix->GetChargeIndexRange(PLevel, MinQValue, MaxQValue);
	  long Tmp = MPSMatrix->GetBondIndexRange(PLevel, MaxQValue);
	  for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
	    {
	      for (int i = 0; i < Tmp; ++i)
		{
		  long Tmp2 = ((long) MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(i, PLevel, QValue));
		  BlockIndexProductTableNbrElements[Tmp2] = Tmp;
		  BlockIndexProductTableShift[Tmp2] = EffectiveDimension;
		  BlockIndexProductTable[Tmp2] = new long[Tmp];
		  long* TmpBlockIndexProductTable = BlockIndexProductTable[Tmp2];
		  Tmp2 *= TmpBMatrixDimension;
		  for (int j = 0; j < Tmp; ++j)
		    {
		      TmpBlockIndexProductTable[j] = Tmp2 + MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(j, PLevel, QValue);
		      EffectiveBlockIndices[EffectiveDimension] = TmpBlockIndexProductTable[j];		      
		      ++EffectiveDimension;
		    }
		}
	    }
	}
      for (long i = 1l; i < EffectiveDimension; ++i)
	if (EffectiveBlockIndices[i] < EffectiveBlockIndices[i - 1l])
	  cout << "error, unsorted indices" << endl;
      //      SortArrayUpOrdering(EffectiveBlockIndices, EffectiveDimension);
      cout << "E matrix effective dimension = " << EffectiveDimension << "( vs " << (SparseBMatrices[0].GetNbrRow() * SparseBMatrices[0].GetNbrRow()) << ")" << endl;
      
      Architecture.GetArchitecture()->SetDimension(EffectiveDimension);
      TensorProductSparseMatrixSelectedBlockHamiltonian* ETransposeHamiltonian = new TensorProductSparseMatrixSelectedBlockHamiltonian(NbrBMatrices, SparseBMatrices, SparseBMatrices, Coefficients, 
																       EffectiveDimension, EffectiveBlockIndices, 
																       BlockIndexProductTable, BlockIndexProductTableNbrElements, BlockIndexProductTableShift, 
 																       Architecture.GetArchitecture(), Manager.GetInteger("ematrix-memory") << 20);
      cout << "computing left eigenstates : " << endl;
      FQHEMPSEMatrixMainTask TaskLeft(&Manager, ETransposeHamiltonian, NbrEigenstates, true, true, 1e-10);
      MainTaskOperation TaskOperationLeft (&TaskLeft);
      TaskOperationLeft.ApplyOperation(Architecture.GetArchitecture());
      ComplexVector* LeftEigenstates = TaskLeft.GetEigenstates();
      Complex* LeftEigenvalues = TaskLeft.GetEigenvalues();
      delete ETransposeHamiltonian;

      Architecture.GetArchitecture()->SetDimension(EffectiveDimension);
      TensorProductSparseMatrixSelectedBlockHamiltonian* EHamiltonian = new TensorProductSparseMatrixSelectedBlockHamiltonian(NbrBMatrices, SparseTransposeBMatrices, SparseTransposeBMatrices, Coefficients, 
															      EffectiveDimension, EffectiveBlockIndices, 
 															      BlockIndexProductTable, BlockIndexProductTableNbrElements, BlockIndexProductTableShift, 
 															      Architecture.GetArchitecture(), Manager.GetInteger("ematrix-memory") << 20);
      cout << "computing right eigenstates : " << endl;
      FQHEMPSEMatrixMainTask TaskRight(&Manager, EHamiltonian, NbrEigenstates, true, false, 1e-10);
      MainTaskOperation TaskOperationRight (&TaskRight);
      TaskOperationRight.ApplyOperation(Architecture.GetArchitecture());
      ComplexVector* RightEigenstates = TaskRight.GetEigenstates();
      Complex* RightEigenvalues = TaskRight.GetEigenvalues();
      delete EHamiltonian;

      cout << "eigenvalues : " << endl;
      for (int i = 0; i < NbrEigenstates; ++i)
	cout << LeftEigenvalues[i] << " " << RightEigenvalues[i] 
	     << "   moduli : " << Norm(LeftEigenvalues[i]) << " " << Norm(RightEigenvalues[i]) <<endl;
      cout << endl;
      
      cout << "checking scalar products between left and right E eigenstates:" << endl;
      for (int i = 0; i < NbrEigenstates; ++i)
	for (int j = 0; j < NbrEigenstates; ++j)
	  {
	    Complex Test = 0.0;
	    for (int k = 0; k < LeftEigenstates[i].GetVectorDimension(); ++k)
	      Test += LeftEigenstates[i][k] * RightEigenstates[j][k];
	    cout << "< " << i << " | " << j << " > = " << EuclidianScalarProduct(LeftEigenstates[i], RightEigenstates[j]) << " " << Test << endl;
	  }

      File << "# la na lz shifted_lz lambda -log(lambda)" << endl;
      MPSMatrix->GetChargeIndexRange(0, MinQValue, MaxQValue);

      int NbrEigenstatesLinearSuperposition = NbrEigenstates;
      if (Manager.GetBoolean("use-singlestate"))
	NbrEigenstatesLinearSuperposition = 1;
      Complex* TmpLeftFactors = new Complex [NbrEigenstates];
      Complex* TmpRightFactors = new Complex [NbrEigenstates];
      int ReducedBoundaryIndex = SearchInArray<long>((((long) MPSRowIndex) * TmpBMatrixDimension) + MPSRowIndex, 
						    EffectiveBlockIndices, EffectiveDimension);
      for (int i = 0; i < NbrEigenstates; ++i)
	TmpLeftFactors[i] = RightEigenstates[i][ReducedBoundaryIndex] / EuclidianScalarProduct(LeftEigenstates[i], RightEigenstates[i]);
      ReducedBoundaryIndex = SearchInArray<long>((((long) MPSColumnIndex) * TmpBMatrixDimension) + MPSColumnIndex, 
						EffectiveBlockIndices, EffectiveDimension);
      for (int i = 0; i < NbrEigenstates; ++i)
	TmpRightFactors[i] = LeftEigenstates[i][ReducedBoundaryIndex] / EuclidianScalarProduct(LeftEigenstates[i], RightEigenstates[i]);

      double LeftEigenvalueError = 0.0;
      double RightEigenvalueError = 0.0;
      LeftEigenvalueError = Error;
      RightEigenvalueError = Error;

      double TotalTraceThoA = 0;

      double*** EntanglementSpectrum = new double**[MaxQValue - MinQValue + 1];
      int** EntanglementSpectrumDimension = new int*[MaxQValue - MinQValue + 1];
      
      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
	{
	  EntanglementSpectrum[QValue - MinQValue]  = new double* [Manager.GetInteger("p-truncation") + 1];
	  EntanglementSpectrumDimension[QValue - MinQValue]  = new int [Manager.GetInteger("p-truncation") + 1];
	  for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	    {
	      int LocalMinQValue;
	      int LocalMaxQValue;
	      MPSMatrix->GetChargeIndexRange(PLevel, LocalMinQValue, LocalMaxQValue);
	      if ((QValue >= LocalMinQValue) && (QValue <= LocalMaxQValue))
		{
		  int IndexRange = MPSMatrix->GetBondIndexRange(PLevel, QValue);
		  if (IndexRange >= 0)
		    {
		      EntanglementSpectrumDimension[QValue - MinQValue][PLevel] = 0;
		      EntanglementSpectrum[QValue - MinQValue][PLevel] = new double[IndexRange];
		    }	      
		}
	    }
	}


      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
	{
	  for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	    {
	      int LocalMinQValue;
	      int LocalMaxQValue;
	      MPSMatrix->GetChargeIndexRange(PLevel, LocalMinQValue, LocalMaxQValue);
	      if ((QValue >= LocalMinQValue) && (QValue <= LocalMaxQValue))
		{
		  int IndexRange = MPSMatrix->GetBondIndexRange(PLevel, QValue);
		  if (IndexRange >= 0)
		    {
		      cout << "----------------------------------------------------" << endl;
		      cout << "P=" << PLevel << " " << " Q=" << QValue << endl;
		      ComplexMatrix LeftMDaggerM (IndexRange, IndexRange);		    
		      for (int i = 0; i < IndexRange; ++i)
			for (int j = 0; j < IndexRange; ++j)
			  {
			    LeftMDaggerM[j][i] = 0.0;
			    for (int k = 0; k < NbrEigenstatesLinearSuperposition; ++k)
			      {
				int TmpIndex = SearchInArray<long>((((long) MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(j, PLevel, QValue)) * TmpBMatrixDimension) 
								   + MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(i, PLevel, QValue), 
								   EffectiveBlockIndices, EffectiveDimension);
				LeftMDaggerM[j][i] += LeftEigenstates[k][TmpIndex] * TmpLeftFactors[k];
			      }
			  }
		      
		      
		      ComplexMatrix RightMDaggerM (IndexRange, IndexRange);
		      for (int i = 0; i < IndexRange; ++i)
			for (int j = 0; j < IndexRange; ++j)
			  {
			    RightMDaggerM[j][i] = 0.0;
			    for (int k = 0; k < NbrEigenstatesLinearSuperposition; ++k)
			      {
				int TmpIndex = SearchInArray<long>((((long) MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(j, PLevel, QValue)) * TmpBMatrixDimension) 
								   + MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(i, PLevel, QValue), 
								   EffectiveBlockIndices, EffectiveDimension);
				RightMDaggerM[j][i] += RightEigenstates[k][TmpIndex] * TmpRightFactors[k];
			      }
			  }
		      
		      
		      RealSymmetricMatrix SymLeftMDaggerM (LeftMDaggerM);
		      RealMatrix TmpLeftBasis(SymLeftMDaggerM.GetNbrRow(), SymLeftMDaggerM.GetNbrRow());
		      TmpLeftBasis.SetToIdentity();
		      RealDiagonalMatrix TmpLeftDiag;
		      // 		  if (PLevel == 3)		    
		      // 		    cout << SymLeftMDaggerM << endl;
#ifdef __LAPACK__
		      SymLeftMDaggerM.LapackDiagonalize(TmpLeftDiag, TmpLeftBasis);
#else
		      SymLeftMDaggerM.Diagonalize(TmpLeftDiag, TmpLeftBasis);
#endif
		      int NbrZeroLeftEigenvalues = 0;
		      cout << "scalar product matrix for the left state : " << endl;
		      for (int i = 0; i < TmpLeftDiag.GetNbrRow(); ++i)
			{
			  if (fabs(TmpLeftDiag(i, i)) < LeftEigenvalueError)
			    {
			      ++NbrZeroLeftEigenvalues;	    
			    }
			  //		      cout << TmpLeftDiag(i, i) << endl;
			}
		      cout << "nbr non zero eigenvalues = " << (TmpLeftDiag.GetNbrRow() - NbrZeroLeftEigenvalues) << " (full dim = " << TmpLeftDiag.GetNbrRow() << ")" << endl;
		      
		      RealSymmetricMatrix SymRightMDaggerM (RightMDaggerM);
		      RealMatrix TmpRightBasis(SymRightMDaggerM.GetNbrRow(), SymRightMDaggerM.GetNbrRow());
		      RealDiagonalMatrix TmpRightDiag;
		      TmpRightBasis.SetToIdentity();
		      // 		  if (PLevel == 3)		    
		      // 		    cout << SymRightMDaggerM << endl;
#ifdef __LAPACK__
		      SymRightMDaggerM.LapackDiagonalize(TmpRightDiag, TmpRightBasis);
#else
		      SymRightMDaggerM.Diagonalize(TmpRightDiag, TmpRightBasis);
#endif
		      int NbrZeroRightEigenvalues = 0;
		      cout << "scalar product matrix for the right state : " << endl;
		      for (int i = 0; i < TmpRightDiag.GetNbrRow(); ++i)
			{
			  if (fabs(TmpRightDiag(i, i)) < RightEigenvalueError)
			    {
			      ++NbrZeroRightEigenvalues;	    
			    }
			  //		      cout << TmpRightDiag(i, i) << endl;
			}
		      cout << "nbr non zero eigenvalues = " << (TmpRightDiag.GetNbrRow() - NbrZeroRightEigenvalues) << " (full dim = " << TmpRightDiag.GetNbrRow() << ")"  << endl;
		      if (NbrZeroRightEigenvalues != NbrZeroLeftEigenvalues)
			{
			  cout << "left and right scalar product matrices have different ranks" << endl;
			}
		      if ((NbrZeroLeftEigenvalues < SymLeftMDaggerM.GetNbrRow()) && (NbrZeroRightEigenvalues < SymRightMDaggerM.GetNbrRow()))
			{
			  RealMatrix TruncatedLeftBasis (TmpLeftDiag.GetNbrRow(), TmpLeftDiag.GetNbrRow() -  NbrZeroLeftEigenvalues, true);
			  NbrZeroLeftEigenvalues = 0;
			  for (int i = 0; i < TmpLeftBasis.GetNbrColumn(); ++i)
			    {
			      if (fabs(TmpLeftDiag(i, i)) > LeftEigenvalueError)
				{
				  TruncatedLeftBasis[NbrZeroLeftEigenvalues].Copy(TmpLeftBasis[i]);
				  TruncatedLeftBasis[NbrZeroLeftEigenvalues] *= sqrt(TmpLeftDiag(i, i));
				  ++NbrZeroLeftEigenvalues;
				}
			    }
			  
			  RealMatrix TruncatedRightBasis (TmpRightDiag.GetNbrRow(), TmpRightDiag.GetNbrRow() -  NbrZeroRightEigenvalues, true);
			  NbrZeroRightEigenvalues = 0;
			  for (int i = 0; i < TmpRightBasis.GetNbrColumn(); ++i)
			    {
			      if (fabs(TmpRightDiag(i, i)) > RightEigenvalueError)
				{
				  TruncatedRightBasis[NbrZeroRightEigenvalues].Copy(TmpRightBasis[i]);
				  TruncatedRightBasis[NbrZeroRightEigenvalues] *= sqrt(TmpRightDiag(i, i));
				  ++NbrZeroRightEigenvalues;
				}
			    }
		      
			  
			  RealMatrix TranposedTruncatedRightBasis = TruncatedRightBasis.DuplicateAndTranspose();
			  TruncatedRightBasis.Multiply(TranposedTruncatedRightBasis);		      
			  RealMatrix TranposedTruncatedLeftBasis = TruncatedLeftBasis.DuplicateAndTranspose();
			  TranposedTruncatedLeftBasis.Multiply(TruncatedRightBasis);
			  TranposedTruncatedLeftBasis.Multiply(TruncatedLeftBasis);
			  
			  RealSymmetricMatrix ReducedDensityMatrix ((Matrix&) TranposedTruncatedLeftBasis);
			  
			  RealDiagonalMatrix TmpRhoADiag;
			  if (ReducedDensityMatrix.IsDiagonal() == true)
			    {
			      TmpRhoADiag = RealDiagonalMatrix(TranposedTruncatedLeftBasis);
			    }
			  else
			    {
#ifdef __LAPACK__
			      ReducedDensityMatrix.LapackDiagonalize(TmpRhoADiag);
#else
			      ReducedDensityMatrix.Diagonalize(TmpRhoADiag);
#endif
			    }
			  int NbrNonZeroEigenvalues = 0;
			  double Sum = 0.0;
			  for (int i = 0; i < TmpRhoADiag.GetNbrColumn(); ++i)
			    {
			      if (TmpRhoADiag[i] > 0.0)
				{
				  cout << "xi = " << TmpRhoADiag[i] << endl;
				  EntanglementSpectrum[QValue - MinQValue][PLevel][EntanglementSpectrumDimension[QValue - MinQValue][PLevel]] = TmpRhoADiag[i];
				  ++EntanglementSpectrumDimension[QValue - MinQValue][PLevel];
				  Sum += TmpRhoADiag(i, i);
			      ++NbrNonZeroEigenvalues;
				}
			    }
			  cout << "P=" << PLevel << " " << " Q=" << QValue << " NbrStates=" << NbrNonZeroEigenvalues << " Tr(rho_A)=" << Sum << endl;
			  TotalTraceThoA += Sum;
			}
		    }
		}
	    }
	}

      double EntanglementEntropy = 0.0;
      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
	{
	  for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	    {
	      int LocalMinQValue;
	      int LocalMaxQValue;
	      MPSMatrix->GetChargeIndexRange(PLevel, LocalMinQValue, LocalMaxQValue);
	      if ((QValue >= LocalMinQValue) && (QValue <= LocalMaxQValue))
		{
		  for (int i = 0; i < EntanglementSpectrumDimension[QValue - MinQValue][PLevel]; ++i)
		    {
		      File << "0 " << QValue << " " << PLevel << " " << PLevel << " " 
			   <<  (EntanglementSpectrum[QValue - MinQValue][PLevel][i] / TotalTraceThoA)  
			   <<  " " << (-log(EntanglementSpectrum[QValue - MinQValue][PLevel][i] / TotalTraceThoA)) <<endl;
		      EntanglementEntropy -= (log(EntanglementSpectrum[QValue - MinQValue][PLevel][i] / TotalTraceThoA)
					      * EntanglementSpectrum[QValue - MinQValue][PLevel][i] / TotalTraceThoA);
		    }
		}
	    }
	}
      cout << "S_A=" << EntanglementEntropy << endl;
      cout << "Tr(rho_A)=" << TotalTraceThoA << endl;
      delete[] TmpLeftFactors;
      delete[] TmpRightFactors;
      File.close();
      return 0;
    }

  int MatDim = BMatrices[0].GetNbrRow();
  int LambdaMax = Manager.GetInteger("p-truncation");
  int LaughlinIndex = Manager.GetInteger("laughlin-index");

  cout << "B matrix size = " << MatDim << "x" << MatDim << endl;

  SparseRealMatrix* FinalBMatrices = new SparseRealMatrix[NbrBMatrices];
  for (int j = 0; j < NbrBMatrices; ++j)
    FinalBMatrices[j] = SparseRealMatrix (MatDim, MatDim);

  for(int i = 0; i < MatDim; i++)
   {
    double Tmp;
    for (int j = 0; j < NbrBMatrices; ++j)
      {
	BMatrices[j].GetMatrixElement(i, MPSColumnIndex, Tmp);
	if (Tmp != 0.0)
	  FinalBMatrices[j].SetMatrixElement(i, MPSColumnIndex, Tmp);
      }
   }

  cout<<"Done preparing B matrices and the vectors at 0 and Nphi orbital"<<endl;

  double CutOff = 1e-14;
  double TmpTrace;

  cout<<"Proceed to calculate overlap matrix (full space dimension, and it will be stored) "<<endl;

  long MaxTmpMatrixElements = (((long) BMatrices[0].GetNbrRow()) * 
				((long) BMatrices[0].GetNbrRow() / 1l));
  if (MaxTmpMatrixElements  > (1l << 28))
    MaxTmpMatrixElements  = 1l << 28; 
  cout << "Requested memory for sparse matrix multiplications = " << ((MaxTmpMatrixElements * (2l * sizeof(double) + sizeof(int))) >> 20) << "Mb" << endl;
  double* TmpMatrixElements = new double [MaxTmpMatrixElements];
  int* TmpColumnIndices = new int [MaxTmpMatrixElements];
  double * TmpElements = new double [BMatrices[0].GetNbrRow()];

  double* NormalizationCoefficients = new double[NbrFluxQuanta + 1];
  BinomialCoefficients Binomial(NbrFluxQuanta);
  for (int i = 0; i <= NbrFluxQuanta; ++i)
    {      
      NormalizationCoefficients[i] = ((double) (NbrFluxQuanta + 1)) / Binomial.GetNumericalCoefficient(NbrFluxQuanta, i);
    }

  SparseRealMatrix OverlapMatrix (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());
  OverlapMatrix.SetMatrixElement(MPSRowIndex, MPSRowIndex, 1.0);
  SparseRealMatrix TmpMatrix1;
  SparseRealMatrix TmpMatrix2;
  SparseRealMatrix TmpMatrix3;
  for (int i = 0; i < EntCut; i++)
    {
      TmpMatrix1 = Multiply(ConjugateBMatrices[0], OverlapMatrix, 
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix2 = Multiply(TmpMatrix1, BMatrices[0],  
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix1 = Multiply(ConjugateBMatrices[1], OverlapMatrix, 
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix3 = Multiply(TmpMatrix1, BMatrices[1],  
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 

      if (CylinderFlag == false)
	TmpMatrix3 *= NormalizationCoefficients[i];
      
      OverlapMatrix = TmpMatrix2 + TmpMatrix3;
      
      TmpTrace = OverlapMatrix.Tr();
      if (fabs(TmpTrace) > CutOff)
        OverlapMatrix /= TmpTrace;
      else
        {
          cout << "Warning: trying to normalize with 0! " << endl;
          exit(1);
        } 
      
    }


  cout<<"Compute density matrix in nonorthogonal basis (full space dimension, will be stored)"<<endl;

  SparseRealMatrix RhoA;

  TmpMatrix2 = Multiply(FinalBMatrices[0], FinalBMatrices[0].Transpose(), 
 			TmpMatrixElements, TmpColumnIndices, TmpElements); 

  if (FinalBMatrices[1].ComputeNbrNonZeroMatrixElements() != 0)
   {
     TmpMatrix3 = Multiply(FinalBMatrices[1], FinalBMatrices[1].Transpose(), 
			   TmpMatrixElements, TmpColumnIndices, TmpElements); 
      if (CylinderFlag == false)
	TmpMatrix3 *= NormalizationCoefficients[NbrFluxQuanta];
     RhoA = TmpMatrix2 + TmpMatrix3;
   }
  else
    {
      RhoA = TmpMatrix2;
    }

  for (int i = (NbrFluxQuanta - 1); i >= EntCut; i--)
    {
      TmpMatrix1 = Multiply(BMatrices[0], RhoA, 
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix2 = Multiply(TmpMatrix1, ConjugateBMatrices[0],  
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix1 = Multiply(BMatrices[1], RhoA, 
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix3 = Multiply(TmpMatrix1, ConjugateBMatrices[1],  
			    TmpMatrixElements, TmpColumnIndices, TmpElements);   
      if (CylinderFlag == false)
	TmpMatrix3 *= NormalizationCoefficients[i];
      RhoA = TmpMatrix2 + TmpMatrix3;

      TmpTrace = RhoA.Tr();
      if (fabs(TmpTrace) > CutOff)
        RhoA /= TmpTrace;
      else
        {
          cout << "Warning: trying to normalize with 0! " << endl;
          exit(1);
        } 
    }

  //Free up some space that is no longer needed (this needs to be done in a cleaner way)

  delete[] BMatrices; 
  delete[] ConjugateBMatrices; 
  delete[] FinalBMatrices;
  delete[] TmpElements;

  cout<<"Proceed to calculate ES per momentum P sector (all N sectors)"<<endl;

  double TraceRho = 0.0;
  double* RhoEigenvalues = new double [MatDim];
  int* RhoQSector = new int [MatDim];
  int* RhoPSector = new int [MatDim];

  for (int i =0 ; i < MatDim; i++)
   {
     RhoEigenvalues[i] = 0.0; 
     RhoQSector[i] = 0; 
     RhoPSector[i] = 0;
   }

  long NbrEigenvalues = 0l;

  if (MinimizeOutput == false)
     OverlapMatrix.PrintNonZero(cout) << endl;

 int MinQValue;
 int MaxQValue;

  for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
    {
      MinQValue = 0;
      MaxQValue = 0;
      //cout << "warning, code not compatible with multi Q range" << endl;
      MPSMatrix->GetChargeIndexRange(PLevel, MinQValue, MaxQValue);
      cout << "|P|= " << PLevel << " MinQValue = " << MinQValue << " MaxQValue= " << MaxQValue << endl;
      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
        {

	  SparseRealMatrix TmpOverlapBlock = MPSMatrix->ExtractBlock(OverlapMatrix, PLevel, QValue, PLevel, QValue);
	  SparseRealMatrix RhoABlock = MPSMatrix->ExtractBlock(RhoA, PLevel, QValue, PLevel, QValue);

	  if ((TmpOverlapBlock.ComputeNbrNonZeroMatrixElements() != 0) && (RhoABlock.ComputeNbrNonZeroMatrixElements()))
	    {
              cout<<"------------sector P = "<<PLevel<<" Q = "<< QValue << "---------------"<<endl;

 	      //cout << "QValue=" << QValue << "  PLevel=" << PLevel << " : "<< endl;
	      // 	      TmpOverlapBlock.PrintNonZero(cout) << endl;
	      int TmpSectorDim = TmpOverlapBlock.GetNbrRow();
	      
	      RealSymmetricMatrix HRep (TmpOverlapBlock);
	      RealDiagonalMatrix TmpDiag (TmpSectorDim);
	      RealMatrix Q(TmpSectorDim, TmpSectorDim);
#ifdef __LAPACK__
	      HRep.LapackDiagonalize(TmpDiag, Q);
		      
#else
	      HRep.Diagonalize(TmpDiag, Q);
#endif
	      
	      int NbrNonZeroVectors = 0;
	      for (int j = 0; j < TmpSectorDim; ++j)
		{
		  cout << TmpDiag[j] << " ";
		  if (fabs(TmpDiag[j]) > CutOff)
		    {
		      ++NbrNonZeroVectors;
                      if (TmpDiag[j] < 0)
                        cout << "******* Negative j = " << j << " " << TmpDiag[j] << " ******* "<< endl;
		    }
		}

	      if (NbrNonZeroVectors > 0)
		{
		  RealMatrix BlockBasisLeftMatrix (TmpSectorDim, NbrNonZeroVectors);
		  RealMatrix BlockBasisRightMatrix (TmpSectorDim, NbrNonZeroVectors);
		  NbrNonZeroVectors = 0;
		  for (int j = 0; j < TmpSectorDim; ++j)
		    {
		      if (fabs(TmpDiag[j]) > CutOff)
			{  			
			  BlockBasisLeftMatrix[NbrNonZeroVectors].Copy(Q[j]);
			  if (TmpDiag[j] > 0.0)
			    BlockBasisLeftMatrix[NbrNonZeroVectors] /= sqrt(fabs(TmpDiag[j]));
			  else
			    BlockBasisLeftMatrix[NbrNonZeroVectors] /= -sqrt(fabs(TmpDiag[j]));
			  BlockBasisRightMatrix[NbrNonZeroVectors].Copy(Q[j]);
			  BlockBasisRightMatrix[NbrNonZeroVectors] /= sqrt(fabs(TmpDiag[j]));
			  ++NbrNonZeroVectors;
			}  
		    }
	      
		  
		  //cout<<"-------------------- Start computing rho in the new basis --------------------"<<endl;

		  TmpElements = new double [TmpOverlapBlock.GetNbrRow()];
		  SparseRealMatrix TmpMat = Conjugate(TmpOverlapBlock, RhoABlock, TmpOverlapBlock.Transpose(),
						      TmpMatrixElements, TmpColumnIndices, TmpElements);
		  delete[] TmpElements;
		  RealSymmetricMatrix HRepRho = TmpMat.Conjugate(BlockBasisLeftMatrix, BlockBasisRightMatrix);

		  RealDiagonalMatrix TmpDiagRho (NbrNonZeroVectors);
#ifdef __LAPACK__
		  HRepRho.LapackDiagonalize(TmpDiagRho);
#else
		  HRepRho.Diagonalize(TmpDiagRho);
#endif
		  for (int j = 0; j < NbrNonZeroVectors; ++j)
		    {
		      TraceRho += TmpDiagRho[j];
		      RhoEigenvalues[NbrEigenvalues] = TmpDiagRho[j];
		      RhoQSector[NbrEigenvalues] = QValue;
		      RhoPSector[NbrEigenvalues] = PLevel; 
		      NbrEigenvalues++;
		    }
		}
	    }
	}
    }
  
  cout<<"Trace rho = "<<TraceRho<<endl;
  if (TraceRho > 1e-200)
    for (int i = 0; i < MatDim; ++i)
      RhoEigenvalues[i] /= TraceRho;
  else
    {
      cout << "Very small trace... exiting!" << endl;
      exit(1);
    }  

  //Sort eigenvalues according to P and store reordering map
  int *ReorderingMap = new int [MatDim]; 
  for (int k = 0; k < MatDim; ++k) 
     ReorderingMap[k] = k;
  SortArrayUpOrdering(RhoPSector, ReorderingMap, MatDim);


  int p = 0;
  int q = 0;
  MPSMatrix->GetFillingFactor(p, q);
  cout << "Filling factor nu = p/q = "<<p<<"/"<<q<<endl;
  int gcd = FindGCD(p, q);
  if (gcd > 1)
    {
      p /= gcd;
      q /= gcd;
      cout << "Filling factor nu* = p/q = " << p << "/" << q << endl;
    }

  int TmpNa;
  File << "# l_a    Na    Lz    lambda" << endl;
  for (int i = 0; i < MatDim; ++i)
   { 
    if (((fabs(RhoEigenvalues[ReorderingMap[i]]) > CutOff))) 
      {
        TmpNa = RhoQSector[ReorderingMap[i]];

        MPSMatrix->ComputeGlobalChargeIndexRange(RhoPSector[i], MinQValue, MaxQValue);

        if ((p == 2) && (q == 4)) //Moore-Read
          TmpNa -= (MaxQValue - 1)/2;
        else
          TmpNa -= MaxQValue/2;

        if ((p==1) && (Manager.GetBoolean("use-padding") == false))
           TmpNa -= (q-1)/2;
        else
           if (Manager.GetBoolean("use-padding") == false)
              TmpNa -= p;

        double ParticlesInA =  (p * EntCut - TmpNa)/((double)q); 

        if (Manager.GetBoolean("all-na"))
          {
            cout<< "Na= " << (int)ParticlesInA << " P= " << RhoPSector[i] << " " << RhoEigenvalues[ReorderingMap[i]] << endl;  
            File << EntCut << " " << (int)ParticlesInA << " " << RhoPSector[i] << " " << RhoEigenvalues[ReorderingMap[i]] << endl;
          }
        else if ((int)ParticlesInA == Na)
          {
            cout<< "Na= " << (int)ParticlesInA << " P= " << RhoPSector[i] << " " << RhoEigenvalues[ReorderingMap[i]] << endl;  
            File << EntCut << " " << (int)ParticlesInA << " " << RhoPSector[i] << " " << RhoEigenvalues[ReorderingMap[i]] << endl;
          }
      }
   }
  delete[] TmpMatrixElements;
  delete[] TmpColumnIndices;
  delete[] RhoEigenvalues;
  delete[] RhoPSector;
  delete[] RhoQSector;
  delete[] ReorderingMap;

  File.close();
 
  return 0;
}

