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
#include "GeneralTools/StringTools.h"

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
  
  OptionManager Manager ("FQHEMPSConvertEMatrixEigenstate" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  FQHEMPSMatrixManager MPSMatrixManager(true);

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleStringOption ('i', "input-state", "name of the file that contains the eigenvector to convert");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEMPSConvertEMatrixEigenstate -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("input-state") == 0)
    {
      cout << "error, an input file should be provided" << endl;
      return 0;
    }

  int NbrFluxQuanta = 1;

  AbstractFQHEMPSMatrix* MPSLeftMatrixWithFixedQSector = MPSMatrixManager.GetLeftMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 
  AbstractFQHEMPSMatrix* MPSRightMatrixWithFixedQSector = MPSMatrixManager.GetRightMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 
  AbstractFQHEMPSMatrix* MPSLeftMatrix = MPSLeftMatrixWithFixedQSector->GetParentMPSMatrices();
  AbstractFQHEMPSMatrix* MPSRightMatrix = MPSRightMatrixWithFixedQSector->GetParentMPSMatrices(); 

  SparseRealMatrix* SparseLeftBMatricesWithFixedQSector = MPSLeftMatrixWithFixedQSector->GetMatrices();
  SparseRealMatrix* SparseRightBMatricesWithFixedQSector = MPSRightMatrixWithFixedQSector->GetMatrices();
  SparseRealMatrix* SparseLeftBMatrices = MPSLeftMatrix->GetMatrices();
  SparseRealMatrix* SparseRightBMatrices = MPSRightMatrix->GetMatrices();

  cout << "Left B matrix size with fixed Q sector = " << SparseLeftBMatricesWithFixedQSector[0].GetNbrRow() << "x" << SparseLeftBMatricesWithFixedQSector[0].GetNbrColumn() << endl;
  cout << "Right B matrix size with fixed Q sector = " << SparseRightBMatricesWithFixedQSector[0].GetNbrRow() << "x" << SparseRightBMatricesWithFixedQSector[0].GetNbrColumn() << endl;

  cout << "Left B matrix size without fixed Q sector = " << SparseLeftBMatrices[0].GetNbrRow() << "x" << SparseLeftBMatrices[0].GetNbrColumn() << endl;
  cout << "Right B matrix size without fixed Q sector = " << SparseRightBMatrices[0].GetNbrRow() << "x" << SparseRightBMatrices[0].GetNbrColumn() << endl;

  long InitialEMatrixEigenstateDimension = ((long) SparseLeftBMatricesWithFixedQSector[0].GetNbrRow()) * ((long) SparseRightBMatricesWithFixedQSector[0].GetNbrRow());
  long TargetEMatrixEigenstateDimension = ((long) SparseLeftBMatrices[0].GetNbrRow()) * ((long) SparseRightBMatrices[0].GetNbrRow());
  
  ComplexVector InitialEMatrixEigenstate;
  if (InitialEMatrixEigenstate.ReadVector(Manager.GetString("input-state")) == false)
    {
      cout << "can't read " << Manager.GetString("input-state") << endl;
      return 0;
    }      
  if (InitialEMatrixEigenstate.GetLargeVectorDimension() != InitialEMatrixEigenstateDimension)
    {
      cout << "error, the dimension of " << Manager.GetString("input-state") << " is " 
	   << InitialEMatrixEigenstate.GetLargeVectorDimension() << ", it should be " << InitialEMatrixEigenstateDimension << endl;
      return 0;
    }
  ComplexVector TargetEMatrixEigenstate (TargetEMatrixEigenstateDimension, true);

  char* OutputFileName = ReplaceString(Manager.GetString("input-state"), "fixedq", "convertedfixedq");
  int LeftBMatricesWithFixedQSectorSize = SparseLeftBMatricesWithFixedQSector[0].GetNbrRow();
  int RightBMatricesWithFixedQSectorSize = SparseRightBMatricesWithFixedQSector[0].GetNbrRow();
  int LeftBMatricesSize = SparseLeftBMatrices[0].GetNbrRow();
  int RightBMatricesSize = SparseRightBMatrices[0].GetNbrRow();
  int* LeftConversionArray = MPSLeftMatrixWithFixedQSector->GetIndexMappingArray();
  int* RightConversionArray = MPSRightMatrixWithFixedQSector->GetIndexMappingArray();
  
  for (int i = 0; i < LeftBMatricesWithFixedQSectorSize; ++i)
    {
      if (LeftConversionArray[i] >= 0)
	{
	  int InitialShift = (i * RightBMatricesWithFixedQSectorSize);
	  int TargetShift = (LeftConversionArray[i] * RightBMatricesSize);
	  for (int j = 0; j < RightBMatricesWithFixedQSectorSize; ++j)
	    {
	      if (RightConversionArray[j] >= 0)
		{
		  TargetEMatrixEigenstate[TargetShift + RightConversionArray[j]] = InitialEMatrixEigenstate[InitialShift + j];
		}
	    }
	}
    }
  if (TargetEMatrixEigenstate.WriteVector(OutputFileName) == false)
    {
      cout << "can't write " << OutputFileName << endl;
      return 0;
    }      
  return 0;
}
