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


// compute the entanglement spectrum from the overlap matrices of the left and right parts
// 
// mPSMatrix = pointer tothe MPS matrix
// leftOverlapMatrix= reference on the left overlap matrix
// rightOverlapMatrix= reference on the right overlap matrix
// normalizationMatrix = optional normalization matrix (for PES and RES)
// leftPSector = P sector that has to be selected for the left part
// leftCFTSector = CFT sector that has to be selected for the left part
// leftQSector = Q sector that has to be selected for the left part
// rightPSector = P sector that has to be selected for the right part
// rightCFTSector = CFT sector that has to be selected for the right part
// rightQSector = Q sector that has to be selected for the right part
// eigenvalueError = relative error on the eigenvalues below which an eigenvalue is considered to be equal to zero
// eigenstateFileName = if non-zero, save the eigenstate of the reduced density matrix using eigenstateFileName as a prefix
RealDiagonalMatrix FQHEMPSEvaluatePartialEntanglementSpectrum(AbstractFQHEMPSMatrix* mPSMatrix, SparseRealMatrix& leftOverlapMatrix, SparseRealMatrix& rightOverlapMatrix, 
							      SparseRealMatrix& normalizationMatrix,
							      int leftPSector, int leftCFTSector, int leftQSector, int rightPSector, int rightCFTSector, int rightQSector, 
							      double eigenvalueError, char* eigenstateFileName = 0);

// convert the orbital weights into coefficients for the transfer matrices
//
// mPSMatrix = pointer to the MPS matrices
// nbrBMatrices = number of B matrices
// weightOrbitals = array of orbital weights
// coefficients = array storing the coefficients for the transfer matrices
// orbitalIndex = reference on the current orbital index (will be updated by this function)
// reverseFlag = read the weight from the last orbital
void ConvertWeightsToCoefficients(AbstractFQHEMPSMatrix* mPSMatrix, int nbrBMatrices, double* weightOrbitals, double* coefficients, int& orbitalIndex, bool reverseFlag);



  
int main(int argc, char** argv)
{
  cout.precision(14); 
  
  OptionManager Manager ("FQHESphereMPSEntanglementSpectrumParticlePartition" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "na", "number of particles in subsystem A", 0);
  (*SystemGroup) += new BooleanOption ('\n', "all-na", "print all charge sectors");
  (*SystemGroup) += new BooleanOption ('\n', "infinite-cylinder", "evaluate the entanglement spectrum on the infinite cylinder");
  (*SystemGroup) += new BooleanOption  ('\n', "realspace-cut", "use real space partition instead of particle partition");
  (*SystemGroup) += new SingleStringOption  ('\n', "realspace-partition", "geometrical weights that define the real spac partition");
  (*SystemGroup) += new BooleanOption  ('\n', "realspace-bgcharge", "use adaptive background charge to improve accuracy");
  (*SystemGroup) += new BooleanOption ('\n', "use-singlestate", "use a single real eigenstate of the E matrix  when evaluating the infinite entanglement spectrum");
  (*SystemGroup) += new BooleanOption ('\n', "orbital-es", "compute the orbital entanglement spectrum");
  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-orbitals", "number of orbitals for the A part (i negative, use (N_phi+1) / 2)", -1);
  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-fluxquanta", "set the total number of flux quanta and deduce the root partition instead of using the reference-file", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "left-eigenstate", "file containing the transfer matrix left eigenstate");
  (*SystemGroup) += new SingleStringOption  ('\n', "right-eigenstate", "file containing the transfer matrix right eigenstate");
  (*SystemGroup) += new BooleanOption ('\n', "boundary-lefteigenstate", "use the left eigenstate initialized from the boundary conditions instead of the one provided by --left-eigenstate");
  (*SystemGroup) += new BooleanOption ('\n', "boundary-righteigenstate", "use the right eigenstate initialized from the boundary conditions instead of the one provided by --right-eigenstate");
  (*SystemGroup) += new BooleanOption  ('\n', "diagonal-block", "transfer matrix eigenstates are computed only from the block diagonal in P, CFT sector and Q (override autodetect from eigenvector file names)");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "memory", "amount of memory that can used for precalculations (in Mb)", 0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "ematrix-memory", "amount of memory that can used for precalculations of the E matrix (in Mb)", 500);
  (*OutputGroup) += new SingleStringOption  ('o', "output-file", "output file name");

  (*OutputGroup) += new BooleanOption ('n', "normalize-sphere", "express the MPS in the normalized sphere basis");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstates", "compute the reduced density matrix eigenstates");
  (*OutputGroup) += new SingleIntegerOption ('\n', "eigenstates-qsector", 
					     "compute the reduced density matrix eigenstates only in a given Q sector (-1 if all Q sectors have to be computed)", -1);
  (*OutputGroup) += new SingleIntegerOption ('\n', "eigenstates-psector", 
					     "compute the reduced density matrix eigenstates only in a given P sector (-1 if all P sectors have to be computed)", -1);
  (*OutputGroup) += new SingleIntegerOption ('\n', "eigenstates-cftsector", 
					     "compute the reduced density matrix eigenstates only in a given CFT sector (-1 if all CFT sectors have to be computed)", -1);
  (*ArnoldiGroup) += new SingleIntegerOption  ('\n', "full-diag", 
					       "maximum Hilbert space dimension for which full diagonalization is applied", 1000);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "disk", "enable disk storage for the Arnoldi algorithm", false);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "resume", "resume from disk datas", false);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Arnoldi iteration", false); 
  (*ArnoldiGroup) += new  SingleIntegerOption ('\n', "arnoldi-memory", "amount of memory when using the Arnoldi algorithm (in Mb)", 500); 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");  

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereMPSEntanglementSpectrumParticlePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  bool ShowTimeFlag = true;
  int NbrParticles = 0; 
  int NbrFluxQuanta = 0;
  int TotalLz = 0;
  int* ReferenceState = 0;
  if ((Manager.GetBoolean("infinite-cylinder") == false) && (Manager.GetInteger("nbr-fluxquanta") <= 0) && 
      (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false))
    {
      return -1;
    }
  if (Manager.GetInteger("nbr-fluxquanta") > 0)
    {
      NbrFluxQuanta = Manager.GetInteger("nbr-fluxquanta");
    }

  int Na = Manager.GetInteger("na");

  bool CylinderFlag = Manager.GetBoolean("normalize-cylinder");

  int LandauLevel = 0;

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 
  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  int NbrBMatrices = MPSMatrix->GetNbrMatrices();
  cout << "handling " << NbrBMatrices << " B matrices" << endl;
  SparseRealMatrix* BMatrices = MPSMatrix->GetMatrices();
  SparseRealMatrix* ConjugateBMatrices = new SparseRealMatrix[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      ConjugateBMatrices[i] = BMatrices[i].Transpose();
    }

  cout << "B matrix size = " << BMatrices[0].GetNbrRow() << "x" << BMatrices[0].GetNbrColumn() << endl;
  
  int MPSRowIndex = 0;
  int MPSColumnIndex = 0;
  int NbrEigenstates = 0;
  int MinQ;
  int MaxQ;
  MPSMatrix->GetChargeIndexRange(0, MinQ, MaxQ);
  MPSMatrix->GetMatrixBoundaryIndices(MPSRowIndex, MPSColumnIndex, Manager.GetBoolean("use-padding"));
  cout << "using boundary conditions : row index = " << MPSRowIndex << " column index = " << MPSColumnIndex << endl;
  NbrParticles = MPSMatrix->GetMatrixNaturalNbrParticles(NbrFluxQuanta, Manager.GetBoolean("use-padding"));
  NbrEigenstates = MPSMatrix->GetTransferMatrixLargestEigenvalueDegeneracy();
  int PLevel = MPSMatrix->GetTruncationLevel();
  int NbrCFTSectors = MPSMatrix->GetNbrCFTSectors();
  ofstream File;
  File.precision(14);
  ofstream File2;
  File2.precision(14);
  
  char* Extension = new char[256];  
  if (Manager.GetBoolean("orbital-es"))
    {
      sprintf (Extension, "ent");
    }
  else
    {
      if (Manager.GetBoolean("realspace-cut"))
	{
	  ConfigurationParser RealSpaceWeights;
	  if (RealSpaceWeights.Parse(Manager.GetString("realspace-partition")) == false)
	    {
	      RealSpaceWeights.DumpErrors(cout) << endl;
	      return -1;
	    }
	  double* TmpSquareWeights = 0;
	  int TmpNbrOrbitals = 0;
	  if (RealSpaceWeights.GetAsDoubleArray("OrbitalSquareWeights", ' ', TmpSquareWeights, TmpNbrOrbitals) == false)
	    {
	      cout << "OrbitalSquareWeights is not defined or as a wrong value" << endl;
	      return -1;
	    }
	  if (Manager.GetBoolean("realspace-bgcharge") == true)
	    {
	      if ((TmpNbrOrbitals & 1) == 1)
		{
		  cout << "The number of weights in " << Manager.GetString("realspace-partition") << " should be even when using --realspace-bgcharge" << endl;
		  return -1;
		}
	      sprintf (Extension, "norb_%d.acc.rsent", TmpNbrOrbitals);
	    }
	  else
	    {
	      sprintf (Extension, "norb_%d.rsent", TmpNbrOrbitals);
	    }
	}
      else
	{
	  sprintf (Extension, "parent");
	}
    }
  char* EigenstateOutputNamePrefix = 0;
  if (Manager.GetString("output-file") != 0)
    {      
      File2.open(Manager.GetString("output-file"), ios::binary | ios::out);
      char* TmpFileName2 = new char [strlen(Manager.GetString("output-file")) + 16];
      sprintf (TmpFileName2, "%s.full.parent", Manager.GetString("output-file"));
      File.open(TmpFileName2, ios::binary | ios::out);     
    }
  else
    {
      char* StatisticPrefix = new char[16];
      char* TruncationName = new char[32];
      if (Manager.GetBoolean("boson") == true)
	{
	  sprintf(StatisticPrefix, "bosons");
	  sprintf(TruncationName, "plevel_%ld_maxocc_%ld", Manager.GetInteger("p-truncation"), Manager.GetInteger("boson-truncation"));
	}
      else
	{
	  sprintf(StatisticPrefix, "fermions");
	  sprintf(TruncationName, "plevel_%ld", Manager.GetInteger("p-truncation"));
	}
      char* TmpFileName = new char [512];
      char* TmpFileName2 = new char [512];
      EigenstateOutputNamePrefix = new char [512];
      char* StateName = new char [256];      
      strcpy(StateName, MPSMatrix->GetName());
      if (CylinderFlag == true)
	{
	  if (Manager.GetBoolean("infinite-cylinder"))
	    {
	      sprintf(TmpFileName, "%s_infinite_cylinder_%s_perimeter_%f_%s_n_0_2s_0_lz_0.0.full.%s", StatisticPrefix, StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), TruncationName, Extension);
	      sprintf(TmpFileName2, "%s_infinite_cylinder_%s_perimeter_%f_%s_n_0_2s_0_lz_0.0.%s", StatisticPrefix, StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), TruncationName, Extension);
	      sprintf(EigenstateOutputNamePrefix, "%s_infinite_cylinder_%s_perimeter_%f_%s_%s", StatisticPrefix, StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), TruncationName, Extension);
	    }
	  else
	    {
	      sprintf(TmpFileName, "%s_cylinder_%s_perimeter_%f_%s_n_%d_2s_%d_lz_%d.0.full.%s", StatisticPrefix, StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), TruncationName, NbrParticles, NbrFluxQuanta, TotalLz, Extension);
	      sprintf(TmpFileName2, "%s_cylinder_%s_perimeter_%f_%s_n_%d_2s_%d_lz_%d.0.%s", StatisticPrefix, StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), TruncationName, NbrParticles, NbrFluxQuanta, TotalLz, Extension);
	      sprintf(EigenstateOutputNamePrefix, "%s_cylinder_%s_perimeter_%f_%s_n_%d_2s_%d_lz_%d.0.%s", StatisticPrefix, StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), TruncationName, NbrParticles, NbrFluxQuanta, TotalLz, Extension);
	    }
	}
      else
	{
	  if (Manager.GetBoolean("infinite-cylinder"))
	    {
	      sprintf(TmpFileName, "%s_infinite_%s_%s_n_0_2s_0_lz_0.0.full.%s", StatisticPrefix, StateName,
		      TruncationName, Extension);
	      sprintf(TmpFileName2, "%s_infinite_%s_%s_n_0_2s_0_lz_0.0.%s", StatisticPrefix, StateName,
		      TruncationName, Extension);
	      sprintf(EigenstateOutputNamePrefix, "%s_infinite_%s_%s_%s", StatisticPrefix, StateName,
		      TruncationName, Extension);
	    }
	  else
	    {
	      sprintf(TmpFileName, "%s_%s_%s_n_%d_2s_%d_lz_%d.0.full.%s", StatisticPrefix, StateName,
		      TruncationName, NbrParticles, NbrFluxQuanta, TotalLz, Extension);
	      sprintf(TmpFileName2, "%s_%s_%s_n_%d_2s_%d_lz_%d.0.%s", StatisticPrefix, StateName,
		      TruncationName, NbrParticles, NbrFluxQuanta, TotalLz, Extension);
	      sprintf(EigenstateOutputNamePrefix, "%s_%s_%s_n_%d_2s_%d_lz_%d.0.%s", StatisticPrefix, StateName,
		      TruncationName, NbrParticles, NbrFluxQuanta, TotalLz, Extension);
	    }
	}
      File.open(TmpFileName, ios::binary | ios::out);     
      File2.open(TmpFileName2, ios::binary | ios::out);     
   }
  if ((Manager.GetBoolean("orbital-es")) || (Manager.GetBoolean("realspace-cut")))
    {
      File << "# x    Q    P    lambda    -ln(lambda)" << endl;
    }
  else
    {
      File << "#  N    Lz    lambda" << endl;
    }

  if (Manager.GetBoolean("infinite-cylinder"))
    {
      if ((Manager.GetString("left-eigenstate") == 0) && (Manager.GetBoolean("boundary-lefteigenstate") == false))
	{
	  cout << "The transfer matrix left eigenstate has to be provided in infinite-cylinder mode" << endl;
	  return 0;
	}
      if ((Manager.GetString("right-eigenstate") == 0) && (Manager.GetBoolean("boundary-righteigenstate") == false))
	{
	  cout << "The transfer matrix right eigenstate has to be provided in infinite-cylinder mode" << endl;
	  return 0;
	}
      ComplexVector LeftEigenstate;
      ComplexVector RightEigenstate;
      if (Manager.GetBoolean("boundary-lefteigenstate") == false)
	{
	  if (LeftEigenstate.ReadVector(Manager.GetString("left-eigenstate")) == false)
	    {
	      cout << "can't read " << Manager.GetString("left-eigenstate") << endl;
	      return 0;
	    }      
	}
      else
	{
	  LeftEigenstate = ComplexVector(BMatrices[0].GetNbrRow() * BMatrices[0].GetNbrRow(), true);
	  LeftEigenstate[(BMatrices[0].GetNbrRow() * MPSColumnIndex) + MPSColumnIndex] = 1.0;
	}
      if (Manager.GetBoolean("boundary-righteigenstate") == false)
	{
	  if (RightEigenstate.ReadVector(Manager.GetString("right-eigenstate")) == false)
	    {
	      cout << "can't read " << Manager.GetString("right-eigenstate") << endl;
	      return 0;
	    }
	}
      else
	{
	  RightEigenstate = ComplexVector(BMatrices[0].GetNbrRow() * BMatrices[0].GetNbrRow(), true);
	  RightEigenstate[(BMatrices[0].GetNbrRow() * MPSRowIndex) + MPSRowIndex] = 1.0;
	}
      int MaxNbrFluxQuantaA = 0;
      int MaxNbrFluxQuantaB = 0;
      double* WeightAOrbitals = 0;
      double* WeightBOrbitals = 0;
      if (Manager.GetBoolean("orbital-es") == false)
	{
	  if (Manager.GetBoolean("realspace-cut") == true)
	    {
	      ConfigurationParser RealSpaceWeights;
	      if (RealSpaceWeights.Parse(Manager.GetString("realspace-partition")) == false)
		{
		  RealSpaceWeights.DumpErrors(cout) << endl;
		  return -1;
		}
	      double* TmpSquareWeights = 0;
	      int TmpNbrOrbitals = 0;
	      if (RealSpaceWeights.GetAsDoubleArray("OrbitalSquareWeights", ' ', TmpSquareWeights, TmpNbrOrbitals) == false)
		{
		  cout << "OrbitalSquareWeights is not defined or as a wrong value" << endl;
		  return -1;
		}
	      MaxNbrFluxQuantaA = TmpNbrOrbitals - 1;
	      MaxNbrFluxQuantaB = TmpNbrOrbitals - 1;
	      WeightAOrbitals = new double[MaxNbrFluxQuantaA + 1];
	      WeightBOrbitals = new double[MaxNbrFluxQuantaB + 1];
	      for (int i = 0; i <= MaxNbrFluxQuantaA ; ++i)
		{
		  WeightAOrbitals[i]= sqrt(TmpSquareWeights[i]);
		  WeightBOrbitals[i]= sqrt(1.0 - TmpSquareWeights[MaxNbrFluxQuantaA - i]);
		}
	    }
	  else
	    {
	      MaxNbrFluxQuantaA = 2 * PLevel;
	      MaxNbrFluxQuantaB = 2 * PLevel;
	      WeightAOrbitals = new double[MaxNbrFluxQuantaA + 1];
	      WeightBOrbitals = new double[MaxNbrFluxQuantaB + 1];
	      for (int i = 0; i < MaxNbrFluxQuantaA ; ++i)
		{
		  WeightAOrbitals[i]= M_SQRT1_2;
		  WeightBOrbitals[i]= M_SQRT1_2;
		}
	    }
	}
      
      int TmpDimension = BMatrices[0].GetNbrRow();
      int QSectorShift = (MaxNbrFluxQuantaA + 1) / MPSMatrix->GetNbrOrbitals();
      if (Manager.GetBoolean("orbital-es") == true)
	{
	  QSectorShift = 0;
	}
      cout << "QSectorShift = " << QSectorShift << endl;
      Complex Factor = EuclidianScalarProduct(LeftEigenstate, RightEigenstate);
      double InvFactor = 1.0 / Factor.Re;
      RealMatrix TmpFullLeftOverlapMatrix(TmpDimension, TmpDimension, true);
      RealMatrix TmpFullRightOverlapMatrix(TmpDimension, TmpDimension, true);
      SparseRealMatrix NormalizationMatrix(BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());
      NormalizationMatrix.SetToIdentity();

      for (int i = 0; i < QSectorShift; ++i)
	{
	  NormalizationMatrix.Multiply(BMatrices[0]);
	}
      if (Manager.GetBoolean("orbital-es") == false)
	{
	  // infinite cylinder real space ES with 
	  if ((Manager.GetBoolean("diagonal-block") == true) || 
	      ((Manager.GetString("left-eigenstate") != 0) && (Manager.GetString("right-eigenstate") != 0)
	       && (strstr(Manager.GetString("left-eigenstate"), "_diagblock_") != 0) 
	       && (strstr(Manager.GetString("right-eigenstate"), "_diagblock_") != 0)))
	    {

	      // diagonal block E matrix eigenstates
	      ComplexVector TmpLeftEigenstate (TmpDimension * TmpDimension, true);
	      ComplexVector TmpRightEigenstate (TmpDimension * TmpDimension, true);
	      int TmpBlockDimension = 0;
	      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
		{
		  for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
		    {
		      int MinQValue = 0;
		      int MaxQValue = 0;
		      MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
		      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
			{
			  int TmpLocalDimension = MPSMatrix->GetBondIndexRange(CurrentPLevel, QValue, CurrentCFTSector);
			  TmpBlockDimension += TmpLocalDimension * TmpLocalDimension;
			}
		    }
		}
	      if (LeftEigenstate.GetVectorDimension() != TmpBlockDimension)
		{
		  cout << "error, left eigenstate does not have the expected dimension (" << LeftEigenstate.GetVectorDimension() << " vs " << TmpBlockDimension << ") " << endl;
		  return 0;
		}
	      if (RightEigenstate.GetVectorDimension() != TmpBlockDimension)
		{
		  cout << "error, right eigenstate does not have the expected dimension (" << RightEigenstate.GetVectorDimension() << " vs " << TmpBlockDimension << ") " << endl;
		  return 0;
		}
	      int TmpBlockPosition = 0;
	      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
		{
		  for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
		    {
		      int MinQValue = 0;
		      int MaxQValue = 0;
		      MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
		      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
			{
			  int TmpLocalDimension = MPSMatrix->GetBondIndexRange(CurrentPLevel, QValue, CurrentCFTSector);
			  for (int i = 0; i < TmpLocalDimension; ++i)
			    {
			      int Tmp1 = MPSMatrix->GetBondIndexWithFixedChargePLevelCFTSector(i, CurrentPLevel, QValue, CurrentCFTSector);
			      for (int j = 0; j < TmpLocalDimension; ++j)
				{
				  int Tmp2 = MPSMatrix->GetBondIndexWithFixedChargePLevelCFTSector(j, CurrentPLevel, QValue, CurrentCFTSector);
				  TmpLeftEigenstate[(Tmp1 * TmpDimension) + Tmp2] = LeftEigenstate[TmpBlockPosition]; 
				  TmpRightEigenstate[(Tmp1 * TmpDimension) + Tmp2] = RightEigenstate[TmpBlockPosition]; 
				  ++TmpBlockPosition;
				}
			    }
			}
		    }
		}
	      
	      LeftEigenstate = TmpLeftEigenstate;
	      RightEigenstate = TmpRightEigenstate;
	    }
	  else
	    {
	      // standard E matrix eigenstates	      
	      if (LeftEigenstate.GetVectorDimension() != (TmpDimension * TmpDimension))
		{
		  cout << "error, left eigenstate does not have the expected dimension (" << LeftEigenstate.GetVectorDimension() << " vs " << TmpDimension << ") " << endl;
		  return 0;
		}
	      if (RightEigenstate.GetVectorDimension() != (TmpDimension * TmpDimension))
		{
		  cout << "error, right eigenstate does not have the expected dimension (" << RightEigenstate.GetVectorDimension() << " vs " << TmpDimension << ") " << endl;
		  return 0;
		}
	    }
	  if (Manager.GetBoolean("realspace-bgcharge") == true)
	    {
	      if (((MaxNbrFluxQuantaA + 1) % (2* MPSMatrix->GetNbrOrbitals())) != 0)
		{
		  cout << "The number of weights should be a multiple of 2 * " << MPSMatrix->GetNbrOrbitals() << endl;
		  return -1;
		}
	    }
	  else
	    {
	      if (((MaxNbrFluxQuantaA + 1) % MPSMatrix->GetNbrOrbitals()) != 0)
		{
		  cout << "The number of weights should be a multiple of " << MPSMatrix->GetNbrOrbitals() << endl;
		  return -1;
		}
	    }
	  ComplexVector TmpEigenstate (RightEigenstate.GetVectorDimension());
	  double* Coefficients = new double[NbrBMatrices];
	  int TmpOrbitalIndex = 0;
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalStartingTime), 0);
	    }
	  int NbrEMatrixEvolution = 0;
	  if (Manager.GetBoolean("realspace-bgcharge") == true)
	    {
	      NbrEMatrixEvolution = (MaxNbrFluxQuantaB + 1) / (2 * MPSMatrix->GetNbrOrbitals());
	      cout << "evolving right eigenstate with " << NbrEMatrixEvolution << " transfer matrices and background charge (covering " << ((MaxNbrFluxQuantaB + 1) / 2) << " orbitals)" << endl;	      
	    }
	  else
	    {
	      NbrEMatrixEvolution = (MaxNbrFluxQuantaB + 1) / MPSMatrix->GetNbrOrbitals();
	      cout << "evolving right eigenstate with " << NbrEMatrixEvolution << " transfer matrices (covering " << (MaxNbrFluxQuantaB + 1)<< " orbitals)" << endl;
	    }
	  for (int i = 0; i < NbrEMatrixEvolution; ++i)
	    {		  
	      ConvertWeightsToCoefficients(MPSMatrix, NbrBMatrices, WeightBOrbitals, Coefficients, TmpOrbitalIndex, false);
	      TensorProductSparseMatrixHamiltonian* ETransposeHamiltonian = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, BMatrices, BMatrices, Coefficients,
														     Architecture.GetArchitecture()); 
	      ETransposeHamiltonian->LowLevelMultiply(RightEigenstate, TmpEigenstate);
	      ComplexVector TmpVector = RightEigenstate;		  
	      RightEigenstate = TmpEigenstate;
	      TmpEigenstate = TmpVector;
	      delete ETransposeHamiltonian;
	    }
	  if (Manager.GetBoolean("realspace-bgcharge") == true)
	    {
	      cout << "evolving right eigenstate with " << NbrEMatrixEvolution << " transfer matrices and without background charge (covering " << ((MaxNbrFluxQuantaB + 1) / 2) << " orbitals)" << endl;	      
	      MPSMatrix->ComputeSiteDependentMatrices(0, ((MaxNbrFluxQuantaB + 1) / 2) - 1);
	      // SparseRealMatrix** SiteDependentMatrices;
	      // MPSMatrix->GetSiteDependentMatrices(SiteDependentMatrices);
	      for (int i = 0; i < NbrEMatrixEvolution; ++i)
		{		  
		  ConvertWeightsToCoefficients(MPSMatrix, NbrBMatrices, WeightBOrbitals, Coefficients, TmpOrbitalIndex, false);
		  TensorProductSparseMatrixHamiltonian* ETransposeHamiltonian = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, BMatrices, BMatrices, Coefficients, Architecture.GetArchitecture()); 
		  ETransposeHamiltonian->LowLevelMultiply(RightEigenstate, TmpEigenstate);
		  ComplexVector TmpVector = RightEigenstate;		  
		  RightEigenstate = TmpEigenstate;
		  TmpEigenstate = TmpVector;
		  delete ETransposeHamiltonian;
		}
	    }
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalEndingTime), 0);
	      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "evolution done in " << Dt << "s" << endl;
	    }

	  NbrEMatrixEvolution = 0;
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalStartingTime), 0);
	    }
	  if (Manager.GetBoolean("realspace-bgcharge") == true)
	    {
	      NbrEMatrixEvolution = ((MaxNbrFluxQuantaA + 1) / (2 * MPSMatrix->GetNbrOrbitals()));
	      cout << "evolving left eigenstate with " << NbrEMatrixEvolution <<  " transfer matrices and background charge(covering " << ((MaxNbrFluxQuantaA + 1) / 2) << " orbitals)" << endl;
	    }
	  else
	    {
	      NbrEMatrixEvolution = ((MaxNbrFluxQuantaA + 1) / MPSMatrix->GetNbrOrbitals());
	      cout << "evolving left eigenstate with " << NbrEMatrixEvolution <<  " transfer matrices (covering " << (MaxNbrFluxQuantaA + 1) << " orbitals)" << endl;
	    }
	  TmpOrbitalIndex = 0;
	  for (int i = 0; i < NbrEMatrixEvolution; ++i)
	    {
	      ConvertWeightsToCoefficients(MPSMatrix, NbrBMatrices, WeightAOrbitals, Coefficients, TmpOrbitalIndex, true);
	      TensorProductSparseMatrixHamiltonian* EHamiltonian = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, ConjugateBMatrices, ConjugateBMatrices, 
													    Coefficients, Architecture.GetArchitecture()); 
	      EHamiltonian->LowLevelMultiply(LeftEigenstate, TmpEigenstate);
	      ComplexVector TmpVector = LeftEigenstate;		  
	      LeftEigenstate = TmpEigenstate;
	      TmpEigenstate = TmpVector;
	      delete EHamiltonian;
	    }
	  if (Manager.GetBoolean("realspace-bgcharge") == true)
	    {
	      cout << "evolving left eigenstate with " << NbrEMatrixEvolution <<  " transfer matrices and without background charge(covering " << ((MaxNbrFluxQuantaA + 1) / 2) << " orbitals)" << endl;
	      MPSMatrix->ComputeSiteDependentMatrices(-1, -((MaxNbrFluxQuantaA + 1) / 2));
	      for (int i = 0; i < NbrEMatrixEvolution; ++i)
		{
		  ConvertWeightsToCoefficients(MPSMatrix, NbrBMatrices, WeightAOrbitals, Coefficients, TmpOrbitalIndex, true);
		  TensorProductSparseMatrixHamiltonian* EHamiltonian = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, ConjugateBMatrices, ConjugateBMatrices, Coefficients, Architecture.GetArchitecture()); 
		  EHamiltonian->LowLevelMultiply(LeftEigenstate, TmpEigenstate);
		  ComplexVector TmpVector = LeftEigenstate;		  
		  LeftEigenstate = TmpEigenstate;
		  TmpEigenstate = TmpVector;
		  delete EHamiltonian;
		}
	    }	  
	  delete[] Coefficients;

// 	  LeftEigenstate.PrintNonZero(cout) << endl;
// 	  cout << LeftEigenstate << endl;
// 	  RightEigenstate.PrintNonZero(cout) << endl;
	  for (int i = 0; i < TmpDimension; ++i)
	    for (int j = 0; j < TmpDimension; ++j)
	      {
		TmpFullLeftOverlapMatrix.SetMatrixElement(i, j, InvFactor * LeftEigenstate[i * TmpDimension + j].Re); 
		TmpFullRightOverlapMatrix.SetMatrixElement(i, j, InvFactor * RightEigenstate[i * TmpDimension + j].Re); 
	      }
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalEndingTime), 0);
	      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "evolution done in " << Dt << "s" << endl;
	    }
	}
      else
	{
	  if ((Manager.GetBoolean("diagonal-block") == true) || 
	      ((Manager.GetString("left-eigenstate") != 0) && (Manager.GetString("right-eigenstate") != 0)
	       && (strstr(Manager.GetString("left-eigenstate"), "_diagblock_") != 0) 
	       && (strstr(Manager.GetString("right-eigenstate"), "_diagblock_") != 0)))
	    {
	      // infinite cylinder orbital ES with diagonal block E matrix eigenstates
	      int TmpBlockDimension = 0;
	      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
		{
		  for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
		    {
		      int MinQValue = 0;
		      int MaxQValue = 0;
		      MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
		      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
			{
			  int TmpLocalDimension = MPSMatrix->GetBondIndexRange(CurrentPLevel, QValue, CurrentCFTSector);
			  TmpBlockDimension += TmpLocalDimension * TmpLocalDimension;
			}
		    }
		}
	      if (LeftEigenstate.GetVectorDimension() != TmpBlockDimension)
		{
		  cout << "error, left eigenstate does not have the expected dimension (" << LeftEigenstate.GetVectorDimension() << " vs " << TmpBlockDimension << ") " << endl;
		  return 0;
		}
	      if (RightEigenstate.GetVectorDimension() != TmpBlockDimension)
		{
		  cout << "error, right eigenstate does not have the expected dimension (" << RightEigenstate.GetVectorDimension() << " vs " << TmpBlockDimension << ") " << endl;
		  return 0;
		}
	      
	      int TmpBlockPosition = 0;
	      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
		{
		  for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
		    {
		      int MinQValue = 0;
		      int MaxQValue = 0;
		      MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
		      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
			{
			  int TmpLocalDimension = MPSMatrix->GetBondIndexRange(CurrentPLevel, QValue, CurrentCFTSector);
			  for (int i = 0; i < TmpLocalDimension; ++i)
			    {
			      int Tmp1 = MPSMatrix->GetBondIndexWithFixedChargePLevelCFTSector(i, CurrentPLevel, QValue, CurrentCFTSector);
			      for (int j = 0; j < TmpLocalDimension; ++j)
				{
				  int Tmp2 = MPSMatrix->GetBondIndexWithFixedChargePLevelCFTSector(j, CurrentPLevel, QValue, CurrentCFTSector);
				  TmpFullLeftOverlapMatrix.SetMatrixElement(Tmp1, Tmp2, InvFactor * LeftEigenstate[TmpBlockPosition].Re); 
				  TmpFullRightOverlapMatrix.SetMatrixElement(Tmp1, Tmp2, InvFactor * RightEigenstate[TmpBlockPosition].Re); 
				  ++TmpBlockPosition;
				}
			    }
			}
		    }
		}
	    }
	  else
	    {
	      for (int i = 0; i < TmpDimension; ++i)
		for (int j = 0; j < TmpDimension; ++j)
		  {
		    TmpFullLeftOverlapMatrix.SetMatrixElement(i, j, InvFactor * LeftEigenstate[i * TmpDimension + j].Re); 
		    TmpFullRightOverlapMatrix.SetMatrixElement(i, j, InvFactor * RightEigenstate[i * TmpDimension + j].Re); 
		  }
	    }
	}
      //       cout << TmpFullLeftOverlapMatrix.Tr() << " " << TmpFullRightOverlapMatrix.Tr() << endl;
      //       TmpFullLeftOverlapMatrix.PrintNonZero(cout) << endl;
      //       TmpFullRightOverlapMatrix.PrintNonZero(cout) << endl;
      if (TmpFullLeftOverlapMatrix.Tr() < 0.0)
	TmpFullLeftOverlapMatrix *= -1.0;
      if (TmpFullRightOverlapMatrix.Tr() < 0.0)
	TmpFullRightOverlapMatrix *= -1.0;
      if (TmpFullLeftOverlapMatrix.IsSymmetric())
	{
	  cout << "FullLeftOverlapMatrix is symmetric" << endl;
	}
      else
	{
	  cout << "error : FullLeftOverlapMatrix is not symmetric" << endl;
	  
	}
      if (TmpFullRightOverlapMatrix.IsSymmetric())
	{
	  cout << "FullRightOverlapMatrix is symmetric" << endl;
	}
      else
	{
	  cout << "error : FullRightOverlapMatrix is symmetric" << endl;
	}
      SparseRealMatrix FullLeftOverlapMatrix (TmpFullLeftOverlapMatrix);
      SparseRealMatrix FullRightOverlapMatrix (TmpFullRightOverlapMatrix);

      
      double Error = 1e-13;
      double LeftEigenvalueError = 0.0;
      double RightEigenvalueError = 0.0;
      LeftEigenvalueError = Error;
      RightEigenvalueError = Error;
      double TotalTraceThoA = 0;

      double**** EntanglementSpectrum = new double***[PLevel + 1];
      int*** EntanglementSpectrumDimension = new int**[PLevel + 1];
      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
	{
	  EntanglementSpectrumDimension[CurrentPLevel] = new int*[NbrCFTSectors];
	  EntanglementSpectrum[CurrentPLevel] = new double**[NbrCFTSectors];	      	      
	  for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
	    {
	      int LocalMinQValue;
	      int LocalMaxQValue;
	      MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, LocalMinQValue, LocalMaxQValue);
	      if (LocalMinQValue <=  LocalMaxQValue)
		{
		  EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector] = new int[LocalMaxQValue - LocalMinQValue + 1];
		  EntanglementSpectrum[CurrentPLevel][CurrentCFTSector] = new double*[LocalMaxQValue - LocalMinQValue + 1];	      
		}
	      else
		{
		  EntanglementSpectrum[CurrentPLevel][CurrentCFTSector] = 0;
		}
	      for (int LocalQValue =  LocalMinQValue; LocalQValue <= LocalMaxQValue; ++LocalQValue)
		{
		  cout << "computing sector P=" << CurrentPLevel<< " CFT=" << CurrentCFTSector << " Q=" << LocalQValue << endl;
		  int RightLocalQValue = LocalQValue - QSectorShift;
		  EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = 0;
		  if ((RightLocalQValue <= LocalMaxQValue) && (RightLocalQValue  >= LocalMinQValue))
		    {
		      RealDiagonalMatrix TmpRhoADiag;
		      if ((Manager.GetBoolean("density-eigenstates") == true) && 
			  ((Manager.GetInteger("eigenstates-qsector") < 0) || (Manager.GetInteger("eigenstates-qsector") == LocalQValue)) &&
			  ((Manager.GetInteger("eigenstates-psector") < 0) || (Manager.GetInteger("eigenstates-psector") == CurrentPLevel)) &&
			  ((Manager.GetInteger("eigenstates-cftsector") < 0) || (Manager.GetInteger("eigenstates-cftsector") == CurrentCFTSector)))
			{
			  TmpRhoADiag = FQHEMPSEvaluatePartialEntanglementSpectrum(MPSMatrix, FullLeftOverlapMatrix, FullRightOverlapMatrix, 
										   NormalizationMatrix,
										   CurrentPLevel, CurrentCFTSector, LocalQValue, 
										   CurrentPLevel, CurrentCFTSector, RightLocalQValue, Error, EigenstateOutputNamePrefix);
			}
		      else
			{
			  TmpRhoADiag = FQHEMPSEvaluatePartialEntanglementSpectrum(MPSMatrix, FullLeftOverlapMatrix, FullRightOverlapMatrix, 
										   NormalizationMatrix,
										   CurrentPLevel, CurrentCFTSector, LocalQValue, 
										   CurrentPLevel, CurrentCFTSector, RightLocalQValue, Error);
			}
		      if (TmpRhoADiag.GetNbrRow() > 0)
			{
			  for (int i = 0 ; i < TmpRhoADiag.GetNbrRow(); ++i)
			    {
			      if (TmpRhoADiag[i] > 0.0)
				EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]++;
			    }
			  if (EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] > 0)
			    {
			      EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = new double[EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]];
			      EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = 0;
			      for (int i = 0 ; i < TmpRhoADiag.GetNbrRow(); ++i)
				{
				  if (TmpRhoADiag[i] > 0.0)
				    {
				      EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]] = TmpRhoADiag[i];
				      EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]++;
				      TotalTraceThoA += TmpRhoADiag[i];
				    }
				}
			      SortArrayDownOrdering<double>(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue],
							    EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]);
			    }
			}
		    }
		}
	    }
	}
	  
      
      int GlobalMinQValue = 100000;
      int GlobalMaxQValue = 0;
      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
	{
	  for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
	    {
	      int MinQValue = 0;
	      int MaxQValue = 0;
	      MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
	      if (MinQValue < GlobalMinQValue)
		GlobalMinQValue = MinQValue;
	      if (MaxQValue > GlobalMaxQValue)
		GlobalMaxQValue = MaxQValue;		  
	    }
	}
      double EntanglementEntropy = 0.0;
      for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
	{
	  for (int LocalQValue =  GlobalMinQValue; LocalQValue <= GlobalMaxQValue; ++LocalQValue)
	    {
	      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
		{
		  int LocalMinQValue;
		  int LocalMaxQValue;
		  MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, LocalMinQValue, LocalMaxQValue);
		  if ((LocalQValue >= LocalMinQValue) && (LocalQValue <= LocalMaxQValue))
		    {
		      for (int i = 0; i < EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]; ++i)
			{
			  File << CurrentCFTSector  << " " << LocalQValue << " " 
			       << CurrentPLevel << " "
			       <<  (EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)  
			       <<  " " << (-log(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)) << endl;
			  EntanglementEntropy -= (log(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)
						  * EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA);
			}
		    }
		}
	    }
	}
      cout << "S_A=" << EntanglementEntropy << endl;
      File2 << "inf" << " " << EntanglementEntropy << " 1" << endl;

      File.close();
      File2.close();
      
      cout << "Tr(rho_A)=" << TotalTraceThoA << endl;
      return 0;
    }
  
  
  SparseRealMatrix FullLeftOverlapMatrix (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());
  SparseRealMatrix FullRightOverlapMatrix (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());

  long MaxTmpMatrixElements = (((long) BMatrices[0].GetNbrRow()) * 
				((long) BMatrices[0].GetNbrRow() / 1l));
  long MaxMemoryTmpMatrixElements = Manager.GetInteger("memory");
  MaxMemoryTmpMatrixElements <<= 20;
  MaxMemoryTmpMatrixElements /= (2l * sizeof(double) + sizeof(int));
  if ((MaxMemoryTmpMatrixElements != 0l) && (MaxTmpMatrixElements > MaxMemoryTmpMatrixElements))
    MaxTmpMatrixElements = MaxMemoryTmpMatrixElements;
  cout << "Requested memory for sparse matrix multiplications = " << ((MaxTmpMatrixElements * (2l * sizeof(double) + sizeof(int))) >> 20) << "Mb" << endl;
  double* TmpMatrixElements = new double [MaxTmpMatrixElements];
  int* TmpColumnIndices = new int [MaxTmpMatrixElements];
  double* TmpElements = new double [BMatrices[0].GetNbrRow()];

  int MaxNbrFluxQuantaA = 0;
  int MaxNbrFluxQuantaB = 0;
  double* WeightAOrbitals = 0;
  double* WeightBOrbitals = 0;
  if (Manager.GetBoolean("realspace-cut") == true)
    {
      ConfigurationParser RealSpaceWeights;
      if (RealSpaceWeights.Parse(Manager.GetString("realspace-partition")) == false)
	{
	  RealSpaceWeights.DumpErrors(cout) << endl;
	  return -1;
	}
      double* TmpSquareWeights = 0;
      int TmpNbrOrbitals = 0;
      if (RealSpaceWeights.GetAsDoubleArray("OrbitalSquareWeights", ' ', TmpSquareWeights, TmpNbrOrbitals) == false)
	{
	  cout << "OrbitalSquareWeights is not defined or as a wrong value" << endl;
	  return -1;
	}
      if (TmpNbrOrbitals > (NbrFluxQuanta + 1))
	{
	  cout << "error, the number of weights (" << TmpNbrOrbitals << ") cannot exceed the number of orbitals (" << (NbrFluxQuanta + 1) << ")" << endl;
	  return -1;
	}
      int NbrAOrbitals = (NbrFluxQuanta + 1 + TmpNbrOrbitals) / 2;
      int NbrBOrbitals = (NbrFluxQuanta + 1 + TmpNbrOrbitals) / 2;
      WeightAOrbitals = new double [NbrAOrbitals];
      for (int i = 0; i < (NbrAOrbitals - TmpNbrOrbitals); ++i)
	WeightAOrbitals[i] = 1.0;
      for (int i = NbrAOrbitals - TmpNbrOrbitals; i < NbrAOrbitals; ++i)
	{
	  WeightAOrbitals[i] = sqrt(TmpSquareWeights[i - NbrAOrbitals + TmpNbrOrbitals]);
	}
      WeightBOrbitals = new double [NbrBOrbitals];
      for (int i = 0; i < TmpNbrOrbitals; ++i)
	{
	  WeightBOrbitals[i] = sqrt(1.0 - TmpSquareWeights[i]);
	}
      for (int i = TmpNbrOrbitals; i < NbrBOrbitals; ++i)
	{
	  WeightBOrbitals[i] = 1.0;
	}     
      MaxNbrFluxQuantaA = NbrAOrbitals - 1;
      MaxNbrFluxQuantaB = NbrBOrbitals - 1;
    }
  else
    {
      if (Manager.GetBoolean("orbital-es") == true)
	{
	  if ((Manager.GetInteger("nbr-orbitals") < 0) || (Manager.GetInteger("nbr-orbitals") > (NbrFluxQuanta + 1)))
	    {
	      MaxNbrFluxQuantaA = ((NbrFluxQuanta + 1) / 2) - 1;
	    }
	  else
	    {
	      MaxNbrFluxQuantaA = Manager.GetInteger("nbr-orbitals") - 1;
	    }
	  MaxNbrFluxQuantaB = NbrFluxQuanta - MaxNbrFluxQuantaA - 1;
	}
      else
	{
	  MaxNbrFluxQuantaA = ((NbrFluxQuanta + 1) / 2) - 1 + PLevel;
	  if (MaxNbrFluxQuantaA > NbrFluxQuanta)
	    {
	      MaxNbrFluxQuantaA = NbrFluxQuanta;
	    }
	  MaxNbrFluxQuantaB = MaxNbrFluxQuantaA;
	  int NbrAOrbitals = MaxNbrFluxQuantaA + 1;
	  int NbrBOrbitals = MaxNbrFluxQuantaB + 1;
	  WeightAOrbitals = new double [NbrAOrbitals];
	  WeightBOrbitals = new double [NbrBOrbitals];
	  for (int i = 0; i < (NbrFluxQuanta - MaxNbrFluxQuantaB); ++i)
	    {
	      WeightAOrbitals[i] = 1.0;
	    }
	  for (int i = NbrFluxQuanta - MaxNbrFluxQuantaB; i < NbrAOrbitals; ++i)
	    {	  
	      WeightAOrbitals[i] = M_SQRT1_2;
	    }	  
	  for (int i = 0; i < ((NbrAOrbitals + NbrBOrbitals - NbrFluxQuanta - 1)); ++i)
	    {	  
	      WeightBOrbitals[i] = M_SQRT1_2;
	    }	  
	  for (int i = ((NbrAOrbitals + NbrBOrbitals - NbrFluxQuanta - 1)); i < NbrBOrbitals; ++i)
	    {
	      WeightBOrbitals[i] = 1.0;
	    }     
	}
    }

  cout << "keeping " << (MaxNbrFluxQuantaA + 1) << " orbitals in A" << endl; 
  FullLeftOverlapMatrix.SetMatrixElement(MPSRowIndex, MPSRowIndex, 1.0);  
  for (int i = 0; i <= MaxNbrFluxQuantaA; ++i)
    {
      SparseRealMatrix* TmpMatrices = new SparseRealMatrix[NbrBMatrices];
      for (int j = 0; j < NbrBMatrices; ++j)
	{
	  TmpMatrices[j] = Conjugate(&(ConjugateBMatrices[j]), &FullLeftOverlapMatrix, &(BMatrices[j]), 
				     TmpMatrixElements, TmpColumnIndices, MaxTmpMatrixElements, Architecture.GetArchitecture()); 
	  if ((WeightAOrbitals != 0) && (j > 0))
	    {
	      double TmpWeight = 1.0;
	      for (int k = 1; k <= j; ++k)
		TmpWeight *= WeightAOrbitals[i] * WeightAOrbitals[i];
	      TmpMatrices[j] *= TmpWeight;
	    }
	}
      FullLeftOverlapMatrix = TmpMatrices[0];
      for (int j = 1; j < NbrBMatrices; ++j)
	{
	  FullLeftOverlapMatrix = FullLeftOverlapMatrix + TmpMatrices[j];
	}
      delete[] TmpMatrices;
    }

  FullRightOverlapMatrix.SetMatrixElement(MPSColumnIndex, MPSColumnIndex, 1.0);  
  cout << "keeping " << (MaxNbrFluxQuantaB + 1) << " orbitals in B" << endl; 
  for (int i = 0; i <= MaxNbrFluxQuantaB; ++i)
    {
      SparseRealMatrix* TmpMatrices = new SparseRealMatrix[NbrBMatrices];
      for (int j = 0; j < NbrBMatrices; ++j)
	{
	  TmpMatrices[j] = Conjugate(&(BMatrices[j]), &FullRightOverlapMatrix, &(ConjugateBMatrices[j]), 
				     TmpMatrixElements, TmpColumnIndices, MaxTmpMatrixElements, Architecture.GetArchitecture()); 
	  if ((WeightBOrbitals != 0) && (j > 0))
	    {
	      double TmpWeight = 1.0;
	      for (int k = 1; k <= j; ++k)
		TmpWeight *= WeightBOrbitals[MaxNbrFluxQuantaB - i] * WeightBOrbitals[MaxNbrFluxQuantaB - i];
	      TmpMatrices[j] *= TmpWeight;
	    }	    
	}
      FullRightOverlapMatrix = TmpMatrices[0];
      for (int j = 1; j < NbrBMatrices; ++j)
	{
	  FullRightOverlapMatrix = FullRightOverlapMatrix + TmpMatrices[j];
	}
      delete[] TmpMatrices;
    }

  SparseRealMatrix NormalizationMatrix(BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());
  NormalizationMatrix.SetToIdentity();
  int QSectorShift = MaxNbrFluxQuantaA + 1 + MaxNbrFluxQuantaB + 1 - NbrFluxQuanta - 1;
  for (int i = 0; i < QSectorShift; ++i)
    {
      NormalizationMatrix.Multiply(BMatrices[0]);
    }
  double Error = 1e-13;
  double LeftEigenvalueError = 0.0;
  double RightEigenvalueError = 0.0;
  LeftEigenvalueError = Error;
  RightEigenvalueError = Error;
  double TotalTraceThoA = 0;
  int MinQValue;
  int MaxQValue;
  MPSMatrix->GetChargeIndexRange(0, MinQValue, MaxQValue);
  for (int CurrentPLevel = 1; CurrentPLevel <= PLevel; ++CurrentPLevel)
    {
      int TmpMinQValue;
      int TmpMaxQValue;
      MPSMatrix->GetChargeIndexRange(CurrentPLevel, TmpMinQValue, TmpMaxQValue);
    }


  double**** EntanglementSpectrum = new double***[PLevel + 1];
  int*** EntanglementSpectrumDimension = new int**[PLevel + 1];
  for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
    {
      EntanglementSpectrumDimension[CurrentPLevel] = new int*[NbrCFTSectors];
      EntanglementSpectrum[CurrentPLevel] = new double**[NbrCFTSectors];	      	      
      for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
	{
	  int LocalMinQValue;
	  int LocalMaxQValue;
	  MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, LocalMinQValue, LocalMaxQValue);
	  if (LocalMinQValue <=  LocalMaxQValue)
	    {
	      EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector] = new int[LocalMaxQValue - LocalMinQValue + 1];
	      EntanglementSpectrum[CurrentPLevel][CurrentCFTSector] = new double*[LocalMaxQValue - LocalMinQValue + 1];	      
	    }
	  else
	    {
	      EntanglementSpectrum[CurrentPLevel][CurrentCFTSector] = 0;
	    }
	  for (int LocalQValue =  LocalMinQValue; LocalQValue <= LocalMaxQValue; ++LocalQValue)
	    {
	      cout << "computing sector P=" << CurrentPLevel<< " CFT=" << CurrentCFTSector << " Q=" << LocalQValue << endl;
	      int RightLocalQValue = LocalQValue - QSectorShift;
	      EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = 0;
	      if ((RightLocalQValue <= LocalMaxQValue) && (RightLocalQValue  >= LocalMinQValue))
		{
		  RealDiagonalMatrix TmpRhoADiag = FQHEMPSEvaluatePartialEntanglementSpectrum(MPSMatrix, FullLeftOverlapMatrix, FullRightOverlapMatrix, NormalizationMatrix,
											      CurrentPLevel, CurrentCFTSector, LocalQValue, 
											      CurrentPLevel, CurrentCFTSector, RightLocalQValue, Error);
		  if (TmpRhoADiag.GetNbrRow() > 0)
		    {
		      for (int i = 0 ; i < TmpRhoADiag.GetNbrRow(); ++i)
			{
			  if (TmpRhoADiag[i] > 0.0)
			    EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]++;
			}
		      if (EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] > 0)
			{
			  EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = new double[EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]];
			  EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = 0;
			  for (int i = 0 ; i < TmpRhoADiag.GetNbrRow(); ++i)
			    {
			      if (TmpRhoADiag[i] > 0.0)
				{
				  EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]] = TmpRhoADiag[i];
				  EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]++;
				  TotalTraceThoA += TmpRhoADiag[i];
				}
			    }
			  SortArrayDownOrdering<double>(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue],
							EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]);
			}
		    }
		}
	    }
	}
    }
  
  
  int GlobalMinQValue = 100000;
  int GlobalMaxQValue = 0;
  for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
    {
      for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
	{
	  int MinQValue = 0;
	  int MaxQValue = 0;
	  MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
	  if (MinQValue < GlobalMinQValue)
	    GlobalMinQValue = MinQValue;
	  if (MaxQValue > GlobalMaxQValue)
	    GlobalMaxQValue = MaxQValue;		  
	}
    }
  double EntanglementEntropy = 0.0;
  for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
    {
      for (int LocalQValue =  GlobalMinQValue; LocalQValue <= GlobalMaxQValue; ++LocalQValue)
	{
	  for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
	    {
	      int LocalMinQValue;
	      int LocalMaxQValue;
	      MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, LocalMinQValue, LocalMaxQValue);
	      if ((LocalQValue >= LocalMinQValue) && (LocalQValue <= LocalMaxQValue))
		{
		  for (int i = 0; i < EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]; ++i)
		    {
		      File << CurrentCFTSector  << " " << LocalQValue << " " 
			   << CurrentPLevel << " "
			   <<  (EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)  
			   <<  " " << (-log(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)) << endl;
		      EntanglementEntropy -= (log(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)
					      * EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA);
		    }
		}
	    }
	}
    }
  cout << "S_A=" << EntanglementEntropy << endl;
  File2 << (MaxNbrFluxQuantaA + 1) << " " << EntanglementEntropy << " 1" << endl;
  File.close();
  File2.close();
  cout << "Tr(rho_A)=" << TotalTraceThoA << endl;

  return 0;
}


// compute the entanglement spectrum from the overlap matrices of the left and right parts
// 
// mPSMatrix = pointer tothe MPS matrix
// leftOverlapMatrix= reference on the left overlap matrix
// rightOverlapMatrix= reference on the right overlap matrix
// normalizationMatrix = optional normalization matrix (for PES and RES)
// leftPSector = P sector that has to be selected for the left part
// leftCFTSector = CFT sector that has to be selected for the left part
// leftQSector = Q sector that has to be selected for the left part
// rightPSector = P sector that has to be selected for the right part
// rightCFTSector = CFT sector that has to be selected for the right part
// rightQSector = Q sector that has to be selected for the right part
// eigenvalueError = relative error on the eigenvalues below which an eigenvalue is considered to be equal to zero
// eigenstateFileName = if non-zero, save the eigenstate of the reduced density matrix using eigenstateFileName as a prefix

RealDiagonalMatrix FQHEMPSEvaluatePartialEntanglementSpectrum(AbstractFQHEMPSMatrix* mPSMatrix, SparseRealMatrix& leftOverlapMatrix, SparseRealMatrix& rightOverlapMatrix, 
							      SparseRealMatrix& normalizationMatrix,
							      int leftPSector, int leftCFTSector, int leftQSector, 
							      int rightPSector, int rightCFTSector, int rightQSector, double eigenvalueError, char* eigenstateFileName)
{
  SparseRealMatrix RightOverlapMatrix = mPSMatrix->ExtractBlock(rightOverlapMatrix, rightPSector, rightCFTSector, rightQSector,
								rightPSector, rightCFTSector, rightQSector);
  SparseRealMatrix LeftOverlapMatrix = mPSMatrix->ExtractBlock(leftOverlapMatrix, leftPSector, leftCFTSector, leftQSector, 
							       leftPSector, leftCFTSector, leftQSector);
  RealDiagonalMatrix TmpRhoADiag;
  SparseRealMatrix NormalizationMatrix;
  if ((LeftOverlapMatrix.GetNbrRow() > 0) && (RightOverlapMatrix.GetNbrRow() > 0) && 
      (LeftOverlapMatrix.ComputeNbrNonZeroMatrixElements() > 0l) && (RightOverlapMatrix.ComputeNbrNonZeroMatrixElements() > 0l))
    {
      SparseRealMatrix NormalizationMatrix;
      if (normalizationMatrix.GetNbrRow() > 0)
	{
	  NormalizationMatrix = mPSMatrix->ExtractBlock(normalizationMatrix, rightPSector, rightCFTSector, rightQSector, leftPSector, leftCFTSector, leftQSector);
	  for (int i = 0; i < NormalizationMatrix.GetNbrRow(); ++i)
	    {
	      double Tmp;
	      NormalizationMatrix.GetMatrixElement(i, i, Tmp);
	      NormalizationMatrix.SetMatrixElement(i, i, 1.0 / Tmp);
	    }
	}
      cout << "scalar product matrix for the left part : " << endl;
      RealSymmetricMatrix SymLeftOverlapMatrix (LeftOverlapMatrix);
      RealMatrix TmpLeftBasis(SymLeftOverlapMatrix.GetNbrRow(), SymLeftOverlapMatrix.GetNbrRow());
      TmpLeftBasis.SetToIdentity();
      RealDiagonalMatrix TmpLeftDiag;
#ifdef __LAPACK__
      SymLeftOverlapMatrix.LapackDiagonalize(TmpLeftDiag, TmpLeftBasis);
#else
      SymLeftOverlapMatrix.Diagonalize(TmpLeftDiag, TmpLeftBasis);
#endif
      double LocalLeftEigenvalueError = 0.0;
      for (int i = 0; i < TmpLeftDiag.GetNbrColumn(); ++i)
	if (TmpLeftDiag(i, i) > LocalLeftEigenvalueError)
	  LocalLeftEigenvalueError = TmpLeftDiag(i, i);
      LocalLeftEigenvalueError *= eigenvalueError;
      int NbrZeroLeftEigenvalues = 0;
      for (int i = 0; i < TmpLeftDiag.GetNbrRow(); ++i)
	{
	  if (TmpLeftDiag(i, i) < LocalLeftEigenvalueError)
	    {
	      ++NbrZeroLeftEigenvalues;	    
	    }
	}
      cout << "nbr non zero eigenvalues = " << (TmpLeftDiag.GetNbrRow() - NbrZeroLeftEigenvalues) << " (full dim = " << TmpLeftDiag.GetNbrRow() << ")" << endl;
      
      RealSymmetricMatrix SymRightOverlapMatrix (RightOverlapMatrix);
      RealMatrix TmpRightBasis(SymRightOverlapMatrix.GetNbrRow(), SymRightOverlapMatrix.GetNbrRow());
      RealDiagonalMatrix TmpRightDiag;
      TmpRightBasis.SetToIdentity();
#ifdef __LAPACK__
      SymRightOverlapMatrix.LapackDiagonalize(TmpRightDiag, TmpRightBasis);
#else
      SymRightOverlapMatrix.Diagonalize(TmpRightDiag, TmpRightBasis);
#endif
      int NbrZeroRightEigenvalues = 0;
      cout << "scalar product matrix for the right part : " << endl;
      double LocalRightEigenvalueError = 0.0;
      for (int i = 0; i < TmpRightDiag.GetNbrColumn(); ++i)
	if (TmpRightDiag(i, i) > LocalRightEigenvalueError)
	  LocalRightEigenvalueError = TmpRightDiag(i, i);
      LocalRightEigenvalueError *= eigenvalueError;
      for (int i = 0; i < TmpRightDiag.GetNbrRow(); ++i)
	{
	  if (TmpRightDiag(i, i) < LocalRightEigenvalueError)
	    {
	      ++NbrZeroRightEigenvalues;	    
	    }
	}
      cout << "nbr non zero eigenvalues = " << (TmpRightDiag.GetNbrRow() - NbrZeroRightEigenvalues) << " (full dim = " << TmpRightDiag.GetNbrRow() << ")"  << endl;
      if ((NbrZeroLeftEigenvalues < SymLeftOverlapMatrix.GetNbrRow()) && (NbrZeroRightEigenvalues < SymRightOverlapMatrix.GetNbrRow()))
	{
	  RealMatrix TruncatedLeftBasis (TmpLeftDiag.GetNbrRow(), TmpLeftDiag.GetNbrRow() -  NbrZeroLeftEigenvalues, true);
	  NbrZeroLeftEigenvalues = 0;
	  for (int i = 0; i < TmpLeftBasis.GetNbrColumn(); ++i)
	    {
	      if (TmpLeftDiag(i, i) > LocalLeftEigenvalueError)
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
	      if (TmpRightDiag(i, i) > LocalRightEigenvalueError)
		{
		  TruncatedRightBasis[NbrZeroRightEigenvalues].Copy(TmpRightBasis[i]);
		  TruncatedRightBasis[NbrZeroRightEigenvalues] *= sqrt(TmpRightDiag(i, i));
		  ++NbrZeroRightEigenvalues;
		}
	    }
	  
	  RealMatrix TmpEntanglementMatrix = TruncatedLeftBasis.DuplicateAndTranspose();
 	  if (NormalizationMatrix.GetNbrRow() > 0)
 	    {
 	      RealMatrix TruncatedNormalizationMatrix(NormalizationMatrix);
 	      TmpEntanglementMatrix.Multiply(TruncatedNormalizationMatrix);
 	    }
	  TmpEntanglementMatrix.Multiply(TruncatedRightBasis);
#ifdef __LAPACK__

	  if (eigenstateFileName == 0)
	    {
	      double* TmpValues = TmpEntanglementMatrix.SingularValueDecomposition();
	      int TmpDimension = TmpEntanglementMatrix.GetNbrColumn();
	      if (TmpDimension > TmpEntanglementMatrix.GetNbrRow())
		{
		  TmpDimension = TmpEntanglementMatrix.GetNbrRow();
		}
	      for (int i = 0; i < TmpDimension; ++i)
		{
		  TmpValues[i] *= TmpValues[i];
		}
	      TmpRhoADiag = RealDiagonalMatrix(TmpValues, TmpDimension);
	    }
	  else
	    {
	      RealMatrix TmpU;
	      RealMatrix TmpV;
	      double* TmpValues = TmpEntanglementMatrix.SingularValueDecomposition(TmpU, TmpV);
	      int TmpDimension = TmpEntanglementMatrix.GetNbrColumn();
	      if (TmpDimension > TmpEntanglementMatrix.GetNbrRow())
		{
		  TmpDimension = TmpEntanglementMatrix.GetNbrRow();
		}
	      for (int i = 0; i < TmpDimension; ++i)
		{
		  TmpValues[i] *= TmpValues[i];
		}
	      TmpRhoADiag = RealDiagonalMatrix(TmpValues, TmpDimension);
	      char* OutputFileName = new char[strlen(eigenstateFileName) + 256];
	      for (int i = 0; i < TmpV.GetNbrColumn(); ++i)
		{
		  sprintf (OutputFileName, "%s_p_%d_q_%d_x_%d.%d.vec", eigenstateFileName, leftPSector, leftQSector, leftCFTSector, i);
		  TmpV[i].WriteVector(OutputFileName);
		}
	      delete[] OutputFileName;
	    }
#endif		
	}
    }
  return TmpRhoADiag;
}

// convert the orbital weights into coefficients for the transfer matrices
//
// mPSMatrix = pointer to the MPS matrices
// nbrBMatrices = number of B matrices
// weightOrbitals = array of orbital weights
// coefficients = array storing the coefficients for the transfer matrices
// orbitalIndex = reference on the current orbital index (will be updated by this function)
// reverseFlag = read the weight from the last orbital

void ConvertWeightsToCoefficients(AbstractFQHEMPSMatrix* mPSMatrix, int nbrBMatrices, double* weightOrbitals, double* coefficients, int& orbitalIndex, bool reverseFlag)
{
  unsigned long* TmpPhysicalIndex = new unsigned long[mPSMatrix->GetNbrOrbitals()];
  for (int j = 0; j < nbrBMatrices; ++j)
    {
      double Tmp = 1.0;
      mPSMatrix->GetPhysicalIndex(j, TmpPhysicalIndex);
      for (int k = 0; k <  mPSMatrix->GetNbrOrbitals(); ++k)
	{
	  for (int l = 1; l <= TmpPhysicalIndex[k]; ++l)
	    {
	      if (reverseFlag == true)
		{
		  Tmp *= weightOrbitals[(orbitalIndex + (mPSMatrix->GetNbrOrbitals() - 1 - k))] * weightOrbitals[(orbitalIndex + (mPSMatrix->GetNbrOrbitals() - 1 - k))];
		}
	      else
		{
		  Tmp *= weightOrbitals[(orbitalIndex + k)] * weightOrbitals[(orbitalIndex + k)];
		}
	    }
	}
      coefficients[j] = Tmp;
    }
  delete[] TmpPhysicalIndex;
  orbitalIndex += mPSMatrix->GetNbrOrbitals();
}
