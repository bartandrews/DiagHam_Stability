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


int main(int argc, char** argv)
{
  cout.precision(14); 
  
  OptionManager Manager ("FQHETorusMPSEntanglementSpectrum" , "0.01");
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

  (*SystemGroup) += new SingleIntegerOption ('\n', "topological-sector", "set the topological sector", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "na", "number of particles in subsystem A", 0);
  (*SystemGroup) += new BooleanOption ('\n', "all-na", "print all charge sectors");
  (*SystemGroup) += new BooleanOption  ('\n', "realspace-cut", "use real space partition instead of particle partition");
  (*SystemGroup) += new SingleStringOption  ('\n', "realspace-partition", "geometrical weights that define the real spac partition");
  (*SystemGroup) += new BooleanOption ('\n', "orbital-es", "compute the orbital entanglement spectrum");
  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-orbitals", "number of orbitals for the A part (i negative, use (N_phi+1) / 2)", -1);
  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-fluxquanta", "set the total number of flux quanta and deduce the root partition instead of using the reference-file", 0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "memory", "amount of memory that can used for precalculations (in Mb)", 0);
  (*OutputGroup) += new SingleStringOption  ('o', "output-file", "output file name");
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
      cout << "see man page for option syntax or type FQHETorusMPSEntanglementSpectrum -h" << endl;
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
  if (Manager.GetInteger("nbr-fluxquanta") <= 0)
    {
      cout << "invalid number of flux quanta" << endl;
      return -1;
    }
  NbrFluxQuanta = Manager.GetInteger("nbr-fluxquanta");

  int Na = Manager.GetInteger("na");

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
  NbrParticles = MPSMatrix->GetMatrixNaturalNbrParticles(NbrFluxQuanta, true);
  NbrEigenstates = MPSMatrix->GetTransferMatrixLargestEigenvalueDegeneracy();
  int PLevel = MPSMatrix->GetTruncationLevel();
  int NbrCFTSectors = MPSMatrix->GetNbrCFTSectors();
  double AspectRatio = Manager.GetDouble("aspect-ratio");

  ofstream File;
  File.precision(14);
  ofstream File2;
  File2.precision(14);
  
  char* Extension = new char[8];  
  if (Manager.GetBoolean("orbital-es"))
    {
      sprintf (Extension, "ent");
    }
  else
    {
      sprintf (Extension, "parent");
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
      if (Manager.GetDouble("angle") == 0.0)
	{
	  sprintf(TmpFileName, "%s_torus_kysym_mps_%s_%s_n_%d_2s_%d_ratio_%.6f_0.full.%s", StatisticPrefix, TruncationName, StateName,
		  NbrParticles, NbrFluxQuanta, AspectRatio, Extension);
	  sprintf(TmpFileName2, "%s_torus_kysym_mps_%s_%s_n_%d_2s_%d_ratio_%.6f_0.%s", StatisticPrefix, TruncationName, StateName,
		  NbrParticles, NbrFluxQuanta, AspectRatio, Extension);
	  sprintf(EigenstateOutputNamePrefix, "%s_torus_kysym_mps_%s_%s_n_%d_2s_%d_ratio_%.6f_0.%s", StatisticPrefix, TruncationName, StateName,
		  NbrParticles, NbrFluxQuanta, AspectRatio, Extension);
	}
      else
	{
	  sprintf(TmpFileName, "%s_torus_kysym_mps_%s_%s_n_%d_2s_%d_ratio_%.6f_angle_%.6f_0.full.%s", StatisticPrefix, TruncationName, StateName,
		  NbrParticles, NbrFluxQuanta, AspectRatio, (M_PI * Manager.GetDouble("angle")), Extension);
	  sprintf(TmpFileName2, "%s_torus_kysym_mps_%s_%s_n_%d_2s_%d_ratio_%.6f_angle_%.6f_0.%s", StatisticPrefix, TruncationName, StateName,
		  NbrParticles, NbrFluxQuanta, AspectRatio, (M_PI * Manager.GetDouble("angle")), Extension);
	  sprintf(EigenstateOutputNamePrefix, "%s_torus_kysym_mps_%s_%s_n_%d_2s_%d_ratio_%.6f_angle_%.6f_0.%s", StatisticPrefix, TruncationName, StateName,
		  NbrParticles, NbrFluxQuanta, AspectRatio, (M_PI * Manager.GetDouble("angle")), Extension);
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

  double**** EntanglementSpectrum = new double***[PLevel + 1];
  int*** EntanglementSpectrumDimension = new int**[PLevel + 1];
  double TotalTraceThoA = 0.0;
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
      WeightAOrbitals = new double [NbrAOrbitals];
      for (int i = 0; i < (NbrAOrbitals - TmpNbrOrbitals); ++i)
	WeightAOrbitals[i] = 1.0;
      for (int i = NbrAOrbitals - TmpNbrOrbitals; i < NbrAOrbitals; ++i)
	{
	  WeightAOrbitals[i] = 1.0;//sqrt(TmpSquareWeights[i - NbrAOrbitals + TmpNbrOrbitals]);
	}
      int NbrBOrbitals = (NbrFluxQuanta + 1 + TmpNbrOrbitals) / 2;
      WeightBOrbitals = new double [NbrBOrbitals];
      for (int i = 0; i < TmpNbrOrbitals; ++i)
	{
	  WeightBOrbitals[i] = 1.0;//sqrt(1.0 - TmpSquareWeights[i]);
	    }
      for (int i = TmpNbrOrbitals; i < NbrBOrbitals; ++i)
	{
	  WeightBOrbitals[i] = 1.0;
	}     
      NbrAOrbitals = 14;
      NbrBOrbitals = 14;
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


  if (Manager.GetDouble("angle") == 0.0)
    {

      SparseRealMatrix StringMatrix;
      if (Manager.GetBoolean("boson") == true)
	{
	  StringMatrix = MPSMatrix->GetTorusStringMatrix(0);
	}
      else
	{
	  StringMatrix = MPSMatrix->GetTorusStringMatrix(NbrParticles);
	}

      SparseRealMatrix TransferMatrix = TensorProduct(BMatrices[0], BMatrices[0]) + TensorProduct(BMatrices[1], BMatrices[1]);

      SparseRealMatrix FullLeftOverlapMatrix = TensorProduct(StringMatrix, StringMatrix);

      for (int i = 0; i <= MaxNbrFluxQuantaA; ++i)
	{
	  SparseRealMatrix* TmpMatrices = new SparseRealMatrix[NbrBMatrices];
	  for (int j = 0; j < NbrBMatrices; ++j)
	    {
// 	      TmpMatrices[j] = Conjugate(&(ConjugateBMatrices[j]), &FullLeftOverlapMatrix, &(BMatrices[j]), 
// 					 TmpMatrixElements, TmpColumnIndices, MaxTmpMatrixElements, Architecture.GetArchitecture()); 
	      if ((WeightAOrbitals != 0) && (j > 0))
		{
		  cout << "WeightAOrbitals[" << i << "]=" << WeightAOrbitals[i] << endl;
		  for (int k = 1; k <= j; ++k)
		    TmpMatrices[j] *= WeightAOrbitals[i] * WeightAOrbitals[i];
		}
	    }
	  FullLeftOverlapMatrix = TmpMatrices[0];
	  for (int j = 1; j < NbrBMatrices; ++j)
	    {
	      FullLeftOverlapMatrix = FullLeftOverlapMatrix + TmpMatrices[j];
	    }
	  delete[] TmpMatrices;
	}


      int NbrMPSSumIndices = 0;
      int* MPSSumIndices = MPSMatrix->GetTopologicalSectorIndices(Manager.GetInteger("topological-sector"), NbrMPSSumIndices);
      int* TmpNbrElementPerRow = new int [BMatrices[0].GetNbrRow()];
      for (int i = 0; i < BMatrices[0].GetNbrRow(); ++i)
	TmpNbrElementPerRow[i] = 1;

//       SparseRealMatrix FullLeftOverlapMatrix (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow(), TmpNbrElementPerRow);
      SparseRealMatrix FullRightOverlapMatrix (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());
      for (int i = 0; i < BMatrices[0].GetNbrRow(); ++i)
	FullLeftOverlapMatrix.SetMatrixElement(i, i, 1.0);
      
      
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
      
      for (int i = 0; i <= MaxNbrFluxQuantaA; ++i)
	{
	  SparseRealMatrix* TmpMatrices = new SparseRealMatrix[NbrBMatrices];
	  for (int j = 0; j < NbrBMatrices; ++j)
	    {
	      TmpMatrices[j] = Conjugate(&(ConjugateBMatrices[j]), &FullLeftOverlapMatrix, &(BMatrices[j]), 
					 TmpMatrixElements, TmpColumnIndices, MaxTmpMatrixElements, Architecture.GetArchitecture()); 
	      if ((WeightAOrbitals != 0) && (j > 0))
		{
		  cout << "WeightAOrbitals[" << i << "]=" << WeightAOrbitals[i] << endl;
		  for (int k = 1; k <= j; ++k)
		    TmpMatrices[j] *= WeightAOrbitals[i] * WeightAOrbitals[i];
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
		  cout << "WeightBOrbitals[" << i << "]=" << WeightBOrbitals[i] << endl;
		  for (int k = 1; k <= j; ++k)
		    TmpMatrices[j] *= WeightBOrbitals[i] * WeightBOrbitals[i];
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
      int MinQValue;
      int MaxQValue;
      MPSMatrix->GetChargeIndexRange(0, MinQValue, MaxQValue);
      for (int CurrentPLevel = 1; CurrentPLevel <= PLevel; ++CurrentPLevel)
	{
	  int TmpMinQValue;
	  int TmpMaxQValue;
	  MPSMatrix->GetChargeIndexRange(CurrentPLevel, TmpMinQValue, TmpMaxQValue);
	}
      

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
