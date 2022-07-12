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
  (*SystemGroup) += new SingleDoubleOption  ('\n', "theta", "angle between the components");
  (*SystemGroup) += new BooleanOption  ('\n', "use-fillingfactor", "use the filling factor to set the angle between the components");
  (*SystemGroup) += new SingleStringOption ('\n', "use-productstate", "use a prefined product state defined through a single column formatted ASCII file, instead of finding the optimal one");
  (*SystemGroup) += new BooleanOption  ('\n', "no-normalization", "do not compute the MPS model state normalization");
  (*SystemGroup) += new BooleanOption  ('\n', "left-eigenstates", "compute the left eigenstates");
  (*SystemGroup) += new BooleanOption  ('\n', "fixed-parity", "compute the left eigenstates");

  
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
      cout << "see man page for option syntax or type FQHEMPSMixedEMatrix -h" << endl;
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
  bool ParityFlag = Manager.GetBoolean("fixed-parity");
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
  double Theta =  Manager.GetDouble("theta");
  if (Manager.GetBoolean("use-fillingfactor") == true)
    {
      int Numerator;
      int Denominator;
      MPSLeftMatrix->GetFillingFactor(Numerator, Denominator);
      Theta = acos(sqrt(((double) Numerator) / ((double) Denominator)));
    }
  char * AngleString = new char [50];
  sprintf(AngleString,"theta_%f",Theta);
  char * ParityString = new char [50];
  if (ParityFlag)
    {
      sprintf(ParityString,"%s","fixedparity_");
    }
  else
    sprintf(ParityString,"%s","");


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
		  sprintf(PrefixOutputFileName, "ematrix_diagblock_cylinder_%s_%s%s_perimeter_%f_plevel_%ld_maxocc_%ld", StateName,ParityString,AngleString,  MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"), 
			  Manager.GetInteger("boson-truncation"));
		}
	      else
		{
		  sprintf(PrefixOutputFileName, "ematrix_diagblock_cylinder_%s_%s%s_perimeter_%f_plevel_%ld", StateName,ParityString,AngleString,  MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"));
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("boson") == true)
		{
		  sprintf(PrefixOutputFileName, "ematrix_cylinder_%s_%s%s_perimeter_%f_plevel_%ld_maxocc_%ld", StateName,ParityString,AngleString,
			  MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"), 
			  Manager.GetInteger("boson-truncation"));
		}
	      else
		{
		  sprintf(PrefixOutputFileName, "ematrix_cylinder_%s_%s%s_perimeter_%f_plevel_%ld", StateName,ParityString,AngleString,
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
		  sprintf(PrefixOutputFileName, "ematrix_diagblock_%s_%s%s_plevel_%ld_maxocc_%ld", StateName,ParityString,AngleString,
			  Manager.GetInteger("p-truncation"), Manager.GetInteger("boson-truncation"));
		}
	      else
		{
		  sprintf(PrefixOutputFileName, "ematrix_diagblock_%s_%s%s_plevel_%ld", StateName,ParityString,AngleString,
			  Manager.GetInteger("p-truncation"));
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("boson") == true)
		{
		  sprintf(PrefixOutputFileName, "ematrix_%s_%s%s_plevel_%ld_maxocc_%ld", StateName,ParityString,AngleString,
			  Manager.GetInteger("p-truncation"), Manager.GetInteger("boson-truncation"));
		}
	      else
		{
		  sprintf(PrefixOutputFileName, "ematrix_%s_%s%s_plevel_%ld", StateName,ParityString,AngleString,
			  Manager.GetInteger("p-truncation"));
		}
	    }
	}
      OutputFileName = new char[strlen(PrefixOutputFileName) + 64];
      sprintf(OutputFileName, "%s.dat", PrefixOutputFileName);
   }


  double Error = 1e-13;
  cout <<"NbrBMatrices = "<< NbrBMatrices << endl;
  SparseRealMatrix  MixedEMatrix;
  SparseRealMatrix*  RightMatrices = new SparseRealMatrix[NbrBMatrices];

   //SparseComplexMatrix  MixedEMatrix;
   //SparseComplexMatrix*  RightMatrices = new SparseComplexMatrix[NbrBMatrices];

  double * Coefficients= new double[NbrBMatrices];

  for (int i = 0; i < NbrBMatrices;i++)
    {
      Coefficients[i] = 1;
      RightMatrices[i] = SparseRealMatrix(1,1);
      //RightMatrices[i] = SparseComplexMatrix(1,1);
    }
  
  if (Manager.GetString("use-productstate") == 0)
    {
      RightMatrices[0].SetMatrixElement(0, 0, sin(Theta));
      RightMatrices[1].SetMatrixElement(0, 0, cos(Theta));
    }
  else
    {
      MultiColumnASCIIFile ProductStateFile;
      if (ProductStateFile.Parse(Manager.GetString("use-productstate")) == false)
	{
	  ProductStateFile.DumpErrors(cout);
	  return -1;
	}
      if (ProductStateFile.GetNbrLines() != NbrBMatrices)
	{
	  cout << "error, " << Manager.GetString("use-productstate") << " does not have the proper number of components (" << NbrBMatrices << " vs " << ProductStateFile.GetNbrLines() << ")" << endl;
	  return -1;
	}
      double* TmpCoefficients = ProductStateFile.GetAsDoubleArray(0);
      for (int i = 0; i < NbrBMatrices;i++)
	{
	  RightMatrices[i].SetMatrixElement(0, 0, TmpCoefficients[i]);
	}
    }
  int NbrOrbitals = MPSLeftMatrix->GetNbrOrbitals();
  int NbrStatesPerOrbital = MPSLeftMatrix->GetMaximumOccupation() + 1;
  int NbrStatesPerBlock =  1;
  int NbrRMatrices = 2;
  for (int i = 0; i < NbrOrbitals ; i++)
    NbrStatesPerBlock *= NbrStatesPerOrbital;
  
  
  cout << "NbrOrbitals = " << NbrOrbitals << " NbrStatesPerOrbital = " <<NbrStatesPerOrbital<<" NbrStatesPerBlock =" <<NbrStatesPerBlock<<endl; 
  unsigned long * ArrayPhysicalIndice = MPSLeftMatrix->GetPhysicalIndices();
  
  
  SparseRealMatrix* FusedRMatrices = new SparseRealMatrix [NbrStatesPerBlock];
  //SparseComplexMatrix* FusedRMatrices = new SparseComplexMatrix [NbrStatesPerBlock];

  if (Manager.GetString("use-productstate") == 0)
    {
      int NbrUn =0 ;
      int TmpI;
      for(int i =0 ; i < NbrStatesPerBlock; i++)
	{ 
	  int NbrUn =0 ;
	  TmpI = i;
	  int Index = SearchInUnsortedArray( (unsigned long)( TmpI %  NbrStatesPerOrbital) , ArrayPhysicalIndice,  NbrRMatrices);
	  if (Index <0)
	    {
	      FusedRMatrices[i] = SparseRealMatrix(RightMatrices[0].GetNbrRow(),RightMatrices[0].GetNbrColumn());
	    }
	  else
	    {
	      FusedRMatrices[i].Copy(RightMatrices[Index]);
	    }
	  NbrUn+=TmpI %  NbrStatesPerOrbital;
	  TmpI /= NbrStatesPerOrbital;
	  for(int p = 1; p < NbrOrbitals ; p++)
	    {
	      int Index = SearchInUnsortedArray( (unsigned long)( TmpI %   NbrStatesPerOrbital) , ArrayPhysicalIndice,  NbrRMatrices);
	      if (Index <0)
		{
		  FusedRMatrices[i].ClearMatrix ();
		}
	      else
		{
		  FusedRMatrices[i].Multiply(RightMatrices[Index]);
		}
	      NbrUn+=TmpI %  NbrStatesPerOrbital;
	      TmpI /= NbrStatesPerOrbital;
	    }
	  if (( NbrUn %2 ==0)&&(ParityFlag))
	    Coefficients[i] =0;	
	}
    }
  else
    {
      for (int i = 0; i < NbrBMatrices; ++i)
	{
	  FusedRMatrices[i] = RightMatrices[i];
	}
    }

  TensorProductSparseMatrixHamiltonian* MixedETransposeHamiltonian =0;
  MixedETransposeHamiltonian = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, SparseBMatrices, FusedRMatrices, Coefficients, Architecture.GetArchitecture());
  Architecture.GetArchitecture()->SetDimension(SparseBMatrices[0].GetNbrRow());
  
  

  cout << "Computing the mixed transfer matrix" << endl;
  FQHEMPSEMatrixMainTask TaskMixedLeft(&Manager, MixedETransposeHamiltonian, 2 * NbrEigenstates, false, true, 1e-10, EnergyShift, OutputFileName);
  MainTaskOperation TaskOperationMixedLeft (&TaskMixedLeft);
  TaskOperationMixedLeft.ApplyOperation(Architecture.GetArchitecture());

  if (Manager.GetBoolean("no-normalization") == false)
    {
      Complex* MixedEigenvalues = TaskMixedLeft.GetEigenvalues();
      
      cout << "Computing the model state transfer matrix" << endl;
      TensorProductSparseMatrixHamiltonian* ETransposeHamiltonian =0;
      ETransposeHamiltonian = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, SparseBMatrices, SparseBMatrices, Coefficients, Architecture.GetArchitecture());
      Architecture.GetArchitecture()->SetDimension(SparseBMatrices[0].GetNbrRow());
      
      FQHEMPSEMatrixMainTask TaskLeft(&Manager, ETransposeHamiltonian, NbrEigenstates, false, true, 1e-10, EnergyShift, OutputFileName);
      MainTaskOperation TaskOperationLeft (&TaskLeft);
      TaskOperationLeft.ApplyOperation(Architecture.GetArchitecture());
      
      Complex* Eigenvalues = TaskLeft.GetEigenvalues();
      
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
      char* StateName = new char [256];      
      strcpy(StateName, MPSLeftMatrix->GetName());
      sprintf(TmpFileName, "%s_infinite_cylinder_%s_perimeter_%f_%s_n_0_2s_0_lz_0.0.geoent", StatisticPrefix, StateName,
	      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), TruncationName);
      ofstream File;
      File.open(TmpFileName, ios::binary | ios::out);
      File.precision(14);
      double GeometricalEntropy = -(log(SqrNorm(MixedEigenvalues[0]) / Norm(Eigenvalues[0])) / ((double) NbrOrbitals)) * (MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta) / (2.0 * M_PI));
      cout << "Geometrical entropy = " << GeometricalEntropy << endl;
      File << "# Geo.Ent. Norm(Mixed Eigenvalue)  SqrNorm(Mixed Eigenvalue)" << endl;
      File << GeometricalEntropy << " " << MixedEigenvalues[0] << " " << Eigenvalues[0] << endl;
      File.close();
    }

  return 0;
}

