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
#include "GeneralTools/FilenameTools.h"

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
  
  OptionManager Manager ("FQHEMPSEMatrix" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  FQHEMPSMatrixManager MPSMatrixManager(true);

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  OptionGroup* ArnoldiGroup  = new OptionGroup ("Arnoldi options");
  Architecture.AddOptionGroup(&Manager);
  Manager += ArnoldiGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new BooleanOption  ('\n', "diagonal-block", "consider only the block diagonal in P, CFT sector and Q");
  (*SystemGroup) += new BooleanOption  ('\n', "right-eigenstates", "compute the right eigenstates");
  (*SystemGroup) += new BooleanOption  ('\n', "left-eigenstates", "compute the left eigenstates");
  
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "memory", "amount of memory that can used for precalculations (in Mb)", 500);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "ematrix-memory", "amount of memory that can used for precalculations of the E matrix (in Mb)", 500);
  (*OutputGroup) += new SingleStringOption  ('o', "output-file", "output file name");
  (*OutputGroup) += new BooleanOption  ('\n', "show-ematrix", "show the transfer matrix");
  (*OutputGroup) += new BooleanOption  ('\n', "ematrix-dimonly", "only compute the dimension of the transfer matrix");
  (*ArnoldiGroup) += new SingleIntegerOption  ('\n', "full-diag", 
					       "maximum Hilbert space dimension for which full diagonalization is applied", 1000);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "disk", "enable disk storage for the Arnoldi algorithm", false);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "resume", "resume from disk datas", false);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Arnoldi iteration", false); 
  (*ArnoldiGroup) += new BooleanOption  ('\n', "power-method", "use the power method instead of the Arnoldi algorithm. A single eigenvalue is computed (the one with the largest real part)", false); 
  (*ArnoldiGroup) += new  SingleIntegerOption ('\n', "arnoldi-memory", "amount of memory when using the Arnoldi algorithm (in Mb)", 500); 
  (*ArnoldiGroup) += new BooleanOption  ('\n', "implicitly-restarted", "use the implicitly restarted Arnoldi algorithm", false); 
  (*ArnoldiGroup) += new  SingleIntegerOption ('\n', "nbr-excited", "number of eigenvalues to compute above the groundstate", 0);
  (*ArnoldiGroup) += new  SingleIntegerOption ('\n', "nbr-eigenstates", "number of eigenvalues to compute (if set to zero, this number will be deduced from the state and nbr-excited)", 0);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "sort-real", "sort the eigenvalues with respect to their real part instead of their norm", false); 
  (*ArnoldiGroup) += new SingleDoubleOption ('\n', "arnoldi-precision", "define Arnoldi precision for eigenvalues (0 if automatically defined by the program)", 0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEMPSEMatrix -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0; 
  int NbrFluxQuanta = 1;
  int TotalLz = 0;

  bool CylinderFlag = Manager.GetBoolean("normalize-cylinder");

  int LandauLevel = 0;

  AbstractFQHEMPSMatrix* MPSLeftMatrix = MPSMatrixManager.GetLeftMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 
  AbstractFQHEMPSMatrix* MPSRightMatrix = MPSMatrixManager.GetRightMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 

  
  int MPSRowIndex = 0;
  int MPSColumnIndex = 0;

  SparseRealMatrix* TmpSparseBMatrices = MPSLeftMatrix->GetMatrices();
  SparseRealMatrix* TmpSparseRightBMatrices = MPSRightMatrix->GetMatrices();
  int NbrBMatrices = 0;

  for (int i = 0; i < MPSLeftMatrix->GetNbrMatrices(); ++i)
    {
      if (SearchInUnsortedArray<unsigned long>(MPSLeftMatrix->GetPhysicalIndices()[i], MPSRightMatrix->GetPhysicalIndices(), MPSRightMatrix->GetNbrMatrices()) >= 0)
	++NbrBMatrices;
    }
  SparseRealMatrix* SparseBMatrices = new SparseRealMatrix[NbrBMatrices];
  SparseRealMatrix* SparseRightBMatrices = new SparseRealMatrix[NbrBMatrices];
  NbrBMatrices = 0;
  for (int i = 0; i < MPSLeftMatrix->GetNbrMatrices(); ++i)
    {
      int TmpIndex = SearchInUnsortedArray<unsigned long>(MPSLeftMatrix->GetPhysicalIndices()[i], MPSRightMatrix->GetPhysicalIndices(), MPSRightMatrix->GetNbrMatrices());
      if (TmpIndex >= 0)
	{
	  SparseBMatrices[NbrBMatrices] = TmpSparseBMatrices[i];
	  SparseRightBMatrices[NbrBMatrices] = TmpSparseRightBMatrices[TmpIndex];
	  ++NbrBMatrices;
	}
    }

  cout << "Left B matrix size = " << SparseBMatrices[0].GetNbrRow() << "x" << SparseBMatrices[0].GetNbrColumn() << endl;
  cout << "Right B matrix size = " << SparseRightBMatrices[0].GetNbrRow() << "x" << SparseRightBMatrices[0].GetNbrColumn() << endl;

  char* StateName = 0;
  if (strcmp(MPSLeftMatrix->GetName(), MPSRightMatrix->GetName()) == 0)
    {
      StateName = new char [strlen(MPSLeftMatrix->GetName()) + 1];
      strcpy(StateName, MPSLeftMatrix->GetName());
    }
  else
    {
      StateName = new char [strlen(MPSLeftMatrix->GetName()) + strlen(MPSRightMatrix->GetName()) + 2];
      sprintf (StateName, "%s_%s", MPSLeftMatrix->GetName(), MPSRightMatrix->GetName());
    }

  int NbrEigenstates = MPSLeftMatrix->GetTransferMatrixLargestEigenvalueDegeneracy();
  NbrEigenstates += Manager.GetInteger("nbr-excited");
  if (Manager.GetInteger("nbr-eigenstates") > 0)
    {
      NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
    }
  
  double EnergyShift = 0.0;
  if (Manager.GetBoolean("power-method") == true)
    {
      EnergyShift = 10.0;
      NbrEigenstates = 1;
    }

  char* PrefixOutputFileName = 0;
  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
      if (GetExtensionFromFileName(OutputFileName) != 0)
	{
	  int TmpLength = strlen(OutputFileName) - strlen(GetExtensionFromFileName(OutputFileName));
	  PrefixOutputFileName = new char[TmpLength + 1];
	  strncpy(PrefixOutputFileName, OutputFileName, TmpLength);
	  PrefixOutputFileName[TmpLength] = '\0';
	}
    }
  else
    {
      PrefixOutputFileName  = new char [1024];
      if (CylinderFlag == true)
	{
	  if (Manager.GetBoolean("diagonal-block"))
	    {
	      if (Manager.GetBoolean("boson") == true)
		{
		  sprintf(PrefixOutputFileName, "ematrix_diagblock_cylinder_%s_perimeter_%f_plevel_%ld_maxocc_%ld", StateName,
			  MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"), 
			  Manager.GetInteger("boson-truncation"));
		}
	      else
		{
		  sprintf(PrefixOutputFileName, "ematrix_diagblock_cylinder_%s_perimeter_%f_plevel_%ld", StateName,
			  MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"));
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("boson") == true)
		{
		  sprintf(PrefixOutputFileName, "ematrix_cylinder_%s_perimeter_%f_plevel_%ld_maxocc_%ld", StateName,
			  MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"), 
			  Manager.GetInteger("boson-truncation"));
		}
	      else
		{
		  sprintf(PrefixOutputFileName, "ematrix_cylinder_%s_perimeter_%f_plevel_%ld", StateName,
			  MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"));
		}
	    }
	}
      else
	{
	  if (Manager.GetBoolean("diagonal-block"))
	    {
	      if (Manager.GetBoolean("boson") == true)
		{
		  sprintf(PrefixOutputFileName, "ematrix_diagblock_%s_plevel_%ld_maxocc_%ld", StateName,
			  Manager.GetInteger("p-truncation"), Manager.GetInteger("boson-truncation"));
		}
	      else
		{
		  sprintf(PrefixOutputFileName, "ematrix_diagblock_%s_plevel_%ld", StateName,
			  Manager.GetInteger("p-truncation"));
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("boson") == true)
		{
		  sprintf(PrefixOutputFileName, "ematrix_%s_plevel_%ld_maxocc_%ld", StateName,
			  Manager.GetInteger("p-truncation"), Manager.GetInteger("boson-truncation"));
		}
	      else
		{
		  sprintf(PrefixOutputFileName, "ematrix_%s_plevel_%ld", StateName,
			  Manager.GetInteger("p-truncation"));
		}
	    }
	}
      OutputFileName = new char[strlen(PrefixOutputFileName) + 64];
      sprintf(OutputFileName, "%s.dat", PrefixOutputFileName);
   }


  double Error = 1e-13;
  
  double* Coefficients = new double[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      Coefficients[i] = 1.0;
    }

  int TmpBMatrixDimension = SparseBMatrices[0].GetNbrRow();
  int TmpRightBMatrixDimension = SparseRightBMatrices[0].GetNbrRow();
  
  TensorProductSparseMatrixHamiltonian* ETransposeHamiltonian = 0;
  TensorProductSparseMatrixHamiltonian* EHamiltonian = 0;
  int PTruncationLevel = MPSLeftMatrix->GetTruncationLevel();
  int NbrCFTSectors = MPSLeftMatrix->GetNbrCFTSectors();
  char** EMatrixIndexString = 0;
  if (Manager.GetBoolean("diagonal-block"))
    {
      long EffectiveDimension = 0l;
      for (int PLevel = 0; PLevel <= PTruncationLevel; ++PLevel)
	{
	  for (int i = 0; i < NbrCFTSectors; ++i)
	    {
	      int MinQValue = 0;
	      int MaxQValue = 0;
	      MPSLeftMatrix->GetChargeIndexRange(PLevel, i, MinQValue, MaxQValue);
	      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
		{
		  long Tmp = MPSLeftMatrix->GetBondIndexRange(PLevel, QValue, i);
		  EffectiveDimension += Tmp * Tmp;
		}
	    }
	}
      if (Manager.GetBoolean("ematrix-dimonly") == true)
	{
	  cout << "E matrix effective dimension = " << EffectiveDimension << "( vs " << (SparseBMatrices[0].GetNbrRow() * SparseBMatrices[0].GetNbrRow()) << ")" << endl;
	  return 0;
	}
      cout << "computing effective E matrix indices " << endl;
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
      for (int PLevel = 0; PLevel <= PTruncationLevel; ++PLevel)
	{
	  for (int l = 0; l < NbrCFTSectors; ++l)
	    {
	      int MinQValue = 0;
	      int MaxQValue = 0;
	      MPSLeftMatrix->GetChargeIndexRange(PLevel, l, MinQValue, MaxQValue);
	      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
		{
		  long Tmp = MPSLeftMatrix->GetBondIndexRange(PLevel, QValue, l);
		  for (int i = 0; i < Tmp; ++i)
		    {
		      long Tmp2 = ((long) MPSLeftMatrix->GetBondIndexWithFixedChargePLevelCFTSector(i, PLevel, QValue, l));
		      BlockIndexProductTableNbrElements[Tmp2] = Tmp;
		      BlockIndexProductTableShift[Tmp2] = EffectiveDimension;
		      BlockIndexProductTable[Tmp2] = new long[Tmp];
		      long* TmpBlockIndexProductTable = BlockIndexProductTable[Tmp2];
		      Tmp2 *= TmpBMatrixDimension;
		      for (int j = 0; j < Tmp; ++j)
			{
			  TmpBlockIndexProductTable[j] = Tmp2 + MPSLeftMatrix->GetBondIndexWithFixedChargePLevelCFTSector(j, PLevel, QValue, l);
			  EffectiveBlockIndices[EffectiveDimension] = TmpBlockIndexProductTable[j];		      
			  ++EffectiveDimension;
			}
		    }
		}
	    }
	}
      for (long i = 1l; i < EffectiveDimension; ++i)
	if (EffectiveBlockIndices[i] < EffectiveBlockIndices[i - 1l])
	  cout << "error, unsorted indices" << endl;
      cout << "E matrix effective dimension = " << EffectiveDimension << "( vs " << (SparseBMatrices[0].GetNbrRow() * SparseBMatrices[0].GetNbrRow()) << ")" << endl;
      Architecture.GetArchitecture()->SetDimension(EffectiveDimension);
      long Memory = Manager.GetInteger("ematrix-memory") << 20;
       if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      if ((Manager.GetBoolean("right-eigenstates") == true) || (Manager.GetBoolean("left-eigenstates") == false))
	{
	  ETransposeHamiltonian = new TensorProductSparseMatrixSelectedBlockHamiltonian(NbrBMatrices, SparseBMatrices, SparseBMatrices, Coefficients, EffectiveDimension, EffectiveBlockIndices, BlockIndexProductTable, BlockIndexProductTableNbrElements, BlockIndexProductTableShift, Architecture.GetArchitecture(), Manager.GetInteger("ematrix-memory") << 20);    
	}
      if (Manager.GetBoolean("left-eigenstates") == true)
	{
	  SparseRealMatrix* ConjugateSparseBMatrices = new SparseRealMatrix[NbrBMatrices];
	  for (int i = 0; i < NbrBMatrices; ++i)
	    ConjugateSparseBMatrices[i] = SparseBMatrices[i].Transpose();
	  EHamiltonian = new TensorProductSparseMatrixSelectedBlockHamiltonian(NbrBMatrices, ConjugateSparseBMatrices, ConjugateSparseBMatrices, Coefficients, 
									       EffectiveDimension, EffectiveBlockIndices, 
									       BlockIndexProductTable, BlockIndexProductTableNbrElements, BlockIndexProductTableShift, 
									       Architecture.GetArchitecture(), Manager.GetInteger("ematrix-memory") << 20);    
	}
    }
  else
    {
       if (Manager.GetBoolean("ematrix-dimonly") == true)
	{
	  cout << "E matrix dimension = " << (((long) TmpBMatrixDimension) * ((long) TmpRightBMatrixDimension)) << ")" << endl;
	  return 0;
	}
       Architecture.GetArchitecture()->SetDimension(((long) TmpBMatrixDimension) * ((long) TmpRightBMatrixDimension));
       if ((Manager.GetBoolean("right-eigenstates") == true) || (Manager.GetBoolean("left-eigenstates") == false))
	ETransposeHamiltonian = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, SparseBMatrices, SparseRightBMatrices, Coefficients,
									 Architecture.GetArchitecture()); 
      if (Manager.GetBoolean("left-eigenstates") == true)
	{
	  SparseRealMatrix* ConjugateSparseBMatrices = new SparseRealMatrix[NbrBMatrices];
	  SparseRealMatrix* ConjugateSparseRightBMatrices = new SparseRealMatrix[NbrBMatrices];
	  for (int i = 0; i < NbrBMatrices; ++i)
	    {
	      ConjugateSparseBMatrices[i] = SparseBMatrices[i].Transpose();
	      ConjugateSparseRightBMatrices[i] = SparseRightBMatrices[i].Transpose();
	    }
	  EHamiltonian = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, ConjugateSparseBMatrices, ConjugateSparseRightBMatrices, Coefficients, Architecture.GetArchitecture()); 
	}
      if (Manager.GetBoolean("show-fulllabel"))
	{
	  EMatrixIndexString = new char* [TmpBMatrixDimension * TmpRightBMatrixDimension];
	  int TmpPLevelLeft;
	  int TmpQLeft;
	  int TmpCFTSectorLeft;
	  int TmpPLevelRight;
	  int TmpQRight;
	  int TmpCFTSectorRight;	  
	  for (int i = 0; i < TmpBMatrixDimension; ++i)
	    {
	      MPSLeftMatrix->GetCFTSectorChargeAndPLevelFromMatrixIndex(i, TmpCFTSectorLeft, TmpPLevelLeft, TmpQLeft);
	      for (int j = 0; j < TmpRightBMatrixDimension; ++j)
		{
		  MPSRightMatrix->GetCFTSectorChargeAndPLevelFromMatrixIndex(j, TmpCFTSectorRight, TmpPLevelRight, TmpQRight);
		  EMatrixIndexString[(i * TmpRightBMatrixDimension) + j] = new char [128];
		  sprintf (EMatrixIndexString[(i * TmpRightBMatrixDimension) + j], 
			   "(x=%d, Q=%d, P=%d, i=%d)x(x=%d, Q=%d, P=%d, j=%d)", TmpCFTSectorLeft, TmpQLeft, TmpPLevelLeft, i,
			   TmpCFTSectorRight, TmpQRight, TmpPLevelRight, j);
		}
	    }
	}
    }

  if (Manager.GetBoolean("power-method") == true)
    {
      if (ETransposeHamiltonian != 0)
	ETransposeHamiltonian->ShiftHamiltonian(EnergyShift);
      if (EHamiltonian != 0)
	EHamiltonian->ShiftHamiltonian(EnergyShift);
    }
  
  if (Manager.GetBoolean("show-ematrix"))
    {
      if (ETransposeHamiltonian != 0)
	{
	  ComplexMatrix EMatrix (ETransposeHamiltonian->GetHilbertSpaceDimension(), ETransposeHamiltonian->GetHilbertSpaceDimension());
	  ETransposeHamiltonian->GetHamiltonian(EMatrix);
	  if (Manager.GetBoolean("show-fulllabel"))
	    {
	      EMatrix.PrintNonZero(cout, EMatrixIndexString, EMatrixIndexString);
	    }
	  else
	    {
	      EMatrix.PrintNonZero(cout);
	    }
	}
      if (EHamiltonian != 0)
	{
	  ComplexMatrix EMatrix (EHamiltonian->GetHilbertSpaceDimension(), EHamiltonian->GetHilbertSpaceDimension());
	  EHamiltonian->GetHamiltonian(EMatrix);
	  if (Manager.GetBoolean("show-fulllabel"))
	    {
	      EMatrix.PrintNonZero(cout, EMatrixIndexString, EMatrixIndexString);
	    }
	  else
	    {
	      EMatrix.PrintNonZero(cout);
	    }
	}      
    }

  if ((Manager.GetBoolean("right-eigenstates") == false) && (Manager.GetBoolean("left-eigenstates") == false))
    {
      FQHEMPSEMatrixMainTask TaskLeft(&Manager, ETransposeHamiltonian, NbrEigenstates, false, true, 1e-10, EnergyShift, OutputFileName);
      MainTaskOperation TaskOperationLeft (&TaskLeft);
      TaskOperationLeft.ApplyOperation(Architecture.GetArchitecture());
    }
  if (Manager.GetBoolean("right-eigenstates") == true)
    {
      cout << "computing right eigenstates" << endl;
      char* EigenvectorFileName = new char [strlen(PrefixOutputFileName) + 128];
      sprintf(EigenvectorFileName, "%s_right", PrefixOutputFileName);
      FQHEMPSEMatrixMainTask TaskLeft(&Manager, ETransposeHamiltonian, NbrEigenstates, true, true, 1e-10, EnergyShift, OutputFileName, 0, EigenvectorFileName);
      MainTaskOperation TaskOperationLeft (&TaskLeft);
      TaskOperationLeft.ApplyOperation(Architecture.GetArchitecture());
    }
  if (Manager.GetBoolean("left-eigenstates") == true)
    {
      cout << "computing left eigenstates" << endl;
      char* EigenvectorFileName = new char [strlen(PrefixOutputFileName) + 128];
      sprintf(EigenvectorFileName, "%s_left", PrefixOutputFileName);
      FQHEMPSEMatrixMainTask TaskLeft(&Manager, EHamiltonian, NbrEigenstates, true, false, 1e-10, EnergyShift, OutputFileName, 0, EigenvectorFileName);
      MainTaskOperation TaskOperationLeft (&TaskLeft);
      TaskOperationLeft.ApplyOperation(Architecture.GetArchitecture());
    }
  return 0;
}

