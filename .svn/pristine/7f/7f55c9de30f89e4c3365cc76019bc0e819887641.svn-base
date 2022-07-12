#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "HilbertSpace/FermionOnTorus.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETorusFermionsEntanglementEntropy" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  Manager += SystemGroup;
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "kymax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-ky", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-la", "minimum size of the subsystem whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-la", "maximum size of the subsystem whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new BooleanOption  ('c', "complex", "states of the density matrix are complex");
  (*ToolsGroup) += new BooleanOption  ('\n', "use-svd", "use singular value decomposition instead of diagonalization to compute the entropy");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "kya-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Ky value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusFermionsEntanglementEntropy -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHETorusEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if ((Manager.GetString("ground-file") != 0) && 
      (IsFile(Manager.GetString("ground-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("ground-file") << endl;
      return -1;
    }
  if ((Manager.GetString("degenerated-groundstate") != 0) && 
      (IsFile(Manager.GetString("degenerated-groundstate")) == false))
    {
      cout << "can't open file " << Manager.GetString("degenerated-groundstate") << endl;
      return -1;
    }


  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int KyMax = Manager.GetInteger("kymax"); 
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  bool EigenstateFlag = Manager.GetBoolean("density-eigenstate");
  int FilterKya = Manager.GetInteger("kya-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
  int* TotalKy = 0;
  bool Statistics = true;
  int NbrSpaces = 1;
  FermionOnTorus** Spaces = 0;
  RealVector* GroundStates = 0;
  ComplexVector * ComplexGroundStates = 0;
  char** GroundStateFiles = 0;
  bool SVDFlag = Manager.GetBoolean("use-svd");
  double* Weights =0;
  Complex* ComplexWeights =0;
  bool WeightFlag = false;
  bool ComplexFlag = Manager.GetBoolean("complex");

  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalKy = new int[1];
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));      
      if (ComplexFlag == false)  
	{
	  Weights = new double[1];
	  Weights[0] = 1.0;
	}
      else
	{
	  ComplexWeights = new Complex[1];
	  ComplexWeights[0] = 1.0;
	}
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("degenerated-groundstate")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
       NbrSpaces = DegeneratedFile.GetNbrLines();
       GroundStateFiles = new char* [NbrSpaces];
       TotalKy = new int[NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	 }
      if (DegeneratedFile.GetNbrColumns() > 1)
	{
	  if (ComplexFlag == false)  
	    {
	      Weights = DegeneratedFile.GetAsDoubleArray(1);
	      WeightFlag = true;
	    }
	  else
	    {
	      ComplexWeights = DegeneratedFile.GetAsComplexArray(1);
	      WeightFlag = true;
	    }
	}
      else
	{
	  if (ComplexFlag == false)  
	    {
	      Weights = new double[NbrSpaces];
	      for (int i = 0; i < NbrSpaces; ++i)
		Weights[i] = 1.0;
	    }
	  else
	    {
	      ComplexWeights = new Complex[NbrSpaces];
	      for (int i = 0; i < NbrSpaces; ++i)
		ComplexWeights[i] = 1.0;
	    }
	}
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalKy[i] = 0;
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(GroundStateFiles[i],
						      NbrParticles, KyMax, TotalKy[i], Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
      if (Statistics == false)
	{
	  cout << GroundStateFiles[i] << " is not a fermionic state" << endl;
	  return -1;
	}
    }


  if (ComplexFlag == false)
    {
      GroundStates = new RealVector [NbrSpaces];  
      for (int i = 0; i < NbrSpaces; ++i)
	if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	  {
	    cout << "can't open vector file " << GroundStateFiles[i] << endl;
	    return -1;      
	  }
    }
  else
    {
      ComplexGroundStates = new ComplexVector [NbrSpaces];  
      for (int i = 0; i < NbrSpaces; ++i)
	if (ComplexGroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	  {
	    cout << "can't open vector file " << GroundStateFiles[i] << endl;
	    return -1;      
	  }
    }

  Spaces = new FermionOnTorus* [KyMax];
  for (int i = 0; i < KyMax; ++i)
    {
      Spaces[i] = 0;
    }
  int NbrKySectorGroundStates = 0;
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (Spaces[TotalKy[i]] == 0)
	{
	  Spaces[TotalKy[i]] = new FermionOnTorus (NbrParticles, KyMax, TotalKy[i]);
	  cout << NbrParticles << " " <<  KyMax << " " << TotalKy[i] << endl;
	  if (ComplexFlag == false)
	    {
	      if (Spaces[TotalKy[i]]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
		{
		  cout << "dimension mismatch between Hilbert space (" << Spaces[TotalKy[i]]->GetLargeHilbertSpaceDimension() << ") and ground state (" << GroundStates[i].GetLargeVectorDimension() << ") " << endl;
		  return 0;
		}
	    }
	  else
	    {
	      if (Spaces[TotalKy[i]]->GetLargeHilbertSpaceDimension() != ComplexGroundStates[i].GetLargeVectorDimension())
		{
		  cout << "dimension mismatch between Hilbert space (" << Spaces[TotalKy[i]]->GetLargeHilbertSpaceDimension() << ") and ground state (" << ComplexGroundStates[i].GetLargeVectorDimension() << ") " << endl;
		  return 0;
		}
	    }
	  ++NbrKySectorGroundStates;
	}
      else
	{
	  if (ComplexFlag == false)
	    {
	      if (Spaces[TotalKy[i]]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
		{
		  cout << "dimension mismatch between Hilbert space (" << Spaces[TotalKy[i]]->GetLargeHilbertSpaceDimension() << ") and ground state (" << GroundStates[i].GetLargeVectorDimension() << ") " << endl;
		  return 0;
		}
	    }
	  else
	    {
	      if (Spaces[TotalKy[i]]->GetLargeHilbertSpaceDimension() != ComplexGroundStates[i].GetLargeVectorDimension())
		{
		  cout << "dimension mismatch between Hilbert space (" << Spaces[TotalKy[i]]->GetLargeHilbertSpaceDimension() << ") and ground state (" << ComplexGroundStates[i].GetLargeVectorDimension() << ") " << endl;
		  return 0;
		}
	    }
	}
    }

  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "#  l_a    N    Ky    lambda";
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }

  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "ent");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  File.precision(14);
  cout.precision(14);

  int MeanSubsystemSize = KyMax >> 1;
  if ((KyMax & 1) != 0)
    ++MeanSubsystemSize;
  if (((SingleIntegerOption*) Manager["max-la"])->GetInteger() > 0)
    {
      MeanSubsystemSize = ((SingleIntegerOption*) Manager["max-la"])->GetInteger();
      if (MeanSubsystemSize > KyMax)
	MeanSubsystemSize = KyMax;
    }
  int SubsystemSize = ((SingleIntegerOption*) Manager["min-la"])->GetInteger();
  if (SubsystemSize < 1)
    SubsystemSize = 1;

  for (; SubsystemSize <= MeanSubsystemSize; ++SubsystemSize)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      long TmpDensityMatrixEigenvaluePosition = 0;
      for (int SubsystemNbrParticles = 0; SubsystemNbrParticles <= SubsystemSize; ++SubsystemNbrParticles)
	{
	  int SubsystemTotalKy = (SubsystemNbrParticles * (SubsystemNbrParticles - 1)) / 2; 
	  int SubsystemMaxTotalKy = (SubsystemSize - 1) * SubsystemNbrParticles - SubsystemTotalKy;
	  for (; SubsystemTotalKy <= SubsystemMaxTotalKy; ++SubsystemTotalKy)
	    {
	      cout << "processing subsystem size=" << SubsystemSize << "  subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Ky=" << SubsystemTotalKy << endl;
	      RealSymmetricMatrix PartialDensityMatrix;
	      RealMatrix PartialEntanglementMatrix;
	      HermitianMatrix PartialComplexDensityMatrix;
	      ComplexMatrix PartialComplexEntanglementMatrix;
	      if (SVDFlag == false)
		{
		  if (ComplexFlag == false)
		    {
		      PartialDensityMatrix = Spaces[TotalKy[0]]->EvaluatePartialDensityMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalKy, GroundStates[0]);
		      if (WeightFlag == true)
			PartialDensityMatrix *= Weights[0];
		    }
		  else
		    {
		      PartialComplexDensityMatrix = Spaces[TotalKy[0]]->EvaluatePartialDensityMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalKy, ComplexGroundStates[0]);
		      if (WeightFlag == true)
			PartialComplexDensityMatrix *= Norm(ComplexWeights[0]);
		    }
		}
	      else
		{
		  if (ComplexFlag == false)
		    {
		      if (NbrKySectorGroundStates > 1)
			PartialEntanglementMatrix = Spaces[TotalKy[0]]->EvaluatePartialEntanglementMatrixFullKyPartB(SubsystemSize, SubsystemNbrParticles, SubsystemTotalKy, GroundStates[0]);
		      else
			PartialEntanglementMatrix = Spaces[TotalKy[0]]->EvaluatePartialEntanglementMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalKy, GroundStates[0]);
		      if (WeightFlag == true)
			PartialEntanglementMatrix *= Weights[0];
		    }
		  else
		    {
		      if (NbrKySectorGroundStates > 1)
			PartialComplexEntanglementMatrix = Spaces[TotalKy[0]]->EvaluatePartialEntanglementMatrixFullKyPartB(SubsystemSize, SubsystemNbrParticles, SubsystemTotalKy, ComplexGroundStates[0]);
		      else
			PartialComplexEntanglementMatrix = Spaces[TotalKy[0]]->EvaluatePartialEntanglementMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalKy, ComplexGroundStates[0]);
		      if (WeightFlag == true)
			PartialComplexEntanglementMatrix *= ComplexWeights[0];
		    }
		}
	      for (int i = 1; i < NbrSpaces; ++i)
		{
		  RealSymmetricMatrix TmpMatrix;
		  RealMatrix TmpEntanglementMatrix;
		  HermitianMatrix TmpComplexMatrix;
		  ComplexMatrix TmpComplexEntanglementMatrix;
		  if (SVDFlag == false)
		    {
		      if (ComplexFlag == false)
			{
			  TmpMatrix = Spaces[TotalKy[i]]->EvaluatePartialDensityMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalKy, GroundStates[i]);
			  if (WeightFlag == true)
			    TmpMatrix *= Weights[i];
			  if (PartialDensityMatrix.GetNbrRow() == 0)
			    PartialDensityMatrix = TmpMatrix;
			  else
			    PartialDensityMatrix += TmpMatrix;
			}
		      else
			{
			  TmpComplexMatrix = Spaces[TotalKy[i]]->EvaluatePartialDensityMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalKy, ComplexGroundStates[i]);
			  if (WeightFlag == true)
			    TmpComplexMatrix *= Norm(ComplexWeights[i]);
			  if (PartialComplexDensityMatrix.GetNbrRow() == 0)
			    PartialComplexDensityMatrix = TmpComplexMatrix;
			  else
			    PartialComplexDensityMatrix += TmpComplexMatrix;
			}
		    }
		  else
		    {
		      if (ComplexFlag == false)
			{
			  if (NbrKySectorGroundStates > 1)
			    TmpEntanglementMatrix = Spaces[TotalKy[i]]->EvaluatePartialEntanglementMatrixFullKyPartB(SubsystemSize, SubsystemNbrParticles, SubsystemTotalKy, GroundStates[i]);
			  else
			    TmpEntanglementMatrix = Spaces[TotalKy[i]]->EvaluatePartialEntanglementMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalKy, GroundStates[i]);
			  if (WeightFlag == true)
			    TmpEntanglementMatrix *= Weights[i];
			  if (PartialEntanglementMatrix.GetNbrRow() == 0)
			    PartialEntanglementMatrix = TmpEntanglementMatrix;
			  else
			    PartialEntanglementMatrix += TmpEntanglementMatrix;
			}
		      else
			{
			  if (NbrKySectorGroundStates > 1)
			    TmpComplexEntanglementMatrix = Spaces[TotalKy[i]]->EvaluatePartialEntanglementMatrixFullKyPartB(SubsystemSize, SubsystemNbrParticles, SubsystemTotalKy, ComplexGroundStates[i]);
			  else
			    TmpComplexEntanglementMatrix = Spaces[TotalKy[i]]->EvaluatePartialEntanglementMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalKy, ComplexGroundStates[i]);
			  if (WeightFlag == true)
			    TmpComplexEntanglementMatrix *= ComplexWeights[i];
			  if (PartialComplexEntanglementMatrix.GetNbrRow() == 0)
			    PartialComplexEntanglementMatrix = TmpComplexEntanglementMatrix;
			  else
			    PartialComplexEntanglementMatrix += TmpComplexEntanglementMatrix;
			}
		    }
		}
	      if (SVDFlag == false)
		{
		  if ((NbrSpaces > 1) && (WeightFlag == false))
		    {
		      if (ComplexFlag == false)
			{
			  PartialDensityMatrix /= ((double) NbrSpaces);
			}
		      else
			{
			  PartialComplexDensityMatrix /= ((double) NbrSpaces);
			}
		    }
		}
	      else
		{
		  if ((NbrSpaces > 1) && (WeightFlag == false))
		    {
		      if (ComplexFlag == false)
			{
			  PartialEntanglementMatrix /= sqrt((double) NbrSpaces);
			}
		      else
			{
			  PartialComplexEntanglementMatrix /= sqrt((double) NbrSpaces);
			}
		    }
		}
	      
	      if (((ComplexFlag == false) && ((PartialDensityMatrix.GetNbrRow() > 1) || 
					      ((PartialEntanglementMatrix.GetNbrRow() >= 1) && (PartialEntanglementMatrix.GetNbrColumn() >= 1)))) ||
		  ((ComplexFlag == true) && ((PartialComplexDensityMatrix.GetNbrRow() > 1) || 
					     ((PartialComplexEntanglementMatrix.GetNbrRow() >= 1) && (PartialComplexEntanglementMatrix.GetNbrColumn() >= 1)))))
		{
		  RealDiagonalMatrix TmpDiag;
		  if (SVDFlag == false)
		    {
		      if (ComplexFlag == false)
			{
			  TmpDiag = RealDiagonalMatrix(PartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
			  if (LapackFlag == true)
			    PartialDensityMatrix.LapackDiagonalize(TmpDiag);
			  else
			    PartialDensityMatrix.Diagonalize(TmpDiag);
#else
			  PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
			}
		      else
			{
			  TmpDiag = RealDiagonalMatrix(PartialComplexDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
			  if (LapackFlag == true)
			    PartialComplexDensityMatrix.LapackDiagonalize(TmpDiag);
			  else
			    PartialComplexDensityMatrix.Diagonalize(TmpDiag);
#else
			  PartialComplexDensityMatrix.Diagonalize(TmpDiag);
#endif		  
			}
		    }
		  else
		    {
		      if (((ComplexFlag == false) && (PartialEntanglementMatrix.GetNbrRow() > 1) && (PartialEntanglementMatrix.GetNbrColumn() > 1)) ||
			  ((ComplexFlag == true) && (PartialComplexEntanglementMatrix.GetNbrRow() > 1) && (PartialComplexEntanglementMatrix.GetNbrColumn() > 1)))
			{
			  double* TmpValues;
			  int TmpDimension;
			  if (ComplexFlag == false)
			    {
			      TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
			      TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
			      if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
				{
				  TmpDimension = PartialEntanglementMatrix.GetNbrRow();
				}
			    }
			  else
			    {
			      TmpValues = PartialComplexEntanglementMatrix.SingularValueDecomposition();
			      TmpDimension = PartialComplexEntanglementMatrix.GetNbrColumn();
			      if (TmpDimension > PartialComplexEntanglementMatrix.GetNbrRow())
				{
				  TmpDimension = PartialComplexEntanglementMatrix.GetNbrRow();
				}
			    }
			  for (int i = 0; i < TmpDimension; ++i)
			    TmpValues[i] *= TmpValues[i];
			  TmpDiag = RealDiagonalMatrix(TmpValues, TmpDimension);
			}
		      else
			{
			  double TmpValue = 0.0;
			  if (ComplexFlag == false)
			    {
			      if (PartialEntanglementMatrix.GetNbrRow() == 1)
				{
				  for (int i = 0; i < PartialEntanglementMatrix.GetNbrColumn(); ++i)
				    TmpValue += PartialEntanglementMatrix[i][0] * PartialEntanglementMatrix[i][0];
				}
			      else
				{
				  for (int i = 0; i < PartialEntanglementMatrix.GetNbrRow(); ++i)
				    TmpValue += PartialEntanglementMatrix[0][i] * PartialEntanglementMatrix[0][i];				  
				}
			    }
			  else
			    {
			      if (PartialComplexEntanglementMatrix.GetNbrRow() == 1)
				{
				  for (int i = 0; i < PartialComplexEntanglementMatrix.GetNbrColumn(); ++i)
				    TmpValue += SqrNorm(PartialComplexEntanglementMatrix[i][0]);
				}
			      else
				{
				  for (int i = 0; i < PartialComplexEntanglementMatrix.GetNbrRow(); ++i)
				    TmpValue += SqrNorm(PartialComplexEntanglementMatrix[0][i]);				  
				}
			    }
			  TmpDiag = RealDiagonalMatrix(1, 1);
			  TmpDiag[0] = TmpValue;
			}
		    }
		  TmpDiag.SortMatrixDownOrder();
		  if (DensityMatrixFileName != 0)
		    {
		      ofstream DensityMatrixFile;
		      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		      DensityMatrixFile.precision(14);
		      for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
			DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalKy << " " << TmpDiag[i] << endl;
		      DensityMatrixFile.close();
		    }
		  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
		    {
		      if (TmpDiag[i] > 1e-14)
			{
			  EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
			  DensitySum +=TmpDiag[i];
			}
		    }
		}
	      else
		{
		  if (ComplexFlag == false)
		    {
		      if (PartialDensityMatrix.GetNbrRow() == 1)
			{
			  double TmpValue = PartialDensityMatrix(0,0);
			  if (DensityMatrixFileName != 0)
			    {
			      ofstream DensityMatrixFile;
			      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			      DensityMatrixFile.precision(14);
			      DensityMatrixFile  << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalKy << " " << TmpValue << endl;
			      DensityMatrixFile.close();
			    }		  
			  if (TmpValue > 1e-14)
			    {
			      EntanglementEntropy += TmpValue * log(TmpValue);
			      DensitySum += TmpValue;
			    }
			}
		    }
		  else
		    {
		      if (PartialComplexDensityMatrix.GetNbrRow() == 1)
			{
			  double TmpValue;
			  PartialComplexDensityMatrix.GetMatrixElement(0, 0, TmpValue);
			  if (DensityMatrixFileName != 0)
			    {
			      ofstream DensityMatrixFile;
			      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			      DensityMatrixFile.precision(14);
			      DensityMatrixFile  << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalKy << " " << TmpValue << endl;
			      DensityMatrixFile.close();
			    }		  
			  if (TmpValue > 1e-14)
			    {
			      EntanglementEntropy += TmpValue * log(TmpValue);
			      DensitySum += TmpValue;
			    }
			}
		    }
		}
	    }
	}
      File << SubsystemSize << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << endl;
    }
  File.close();
}

