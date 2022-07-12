#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "HilbertSpace/Spin1_2ChainFixedParity.h"
#include "HilbertSpace/Potts3Chain.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "Options/Options.h"


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(4); 

  // some running options and help
  OptionManager Manager ("PottsChainEntanglementEntropy" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption ('n', "zn", "Zn value of the Potts model", 3);
  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-la", "minimum size of the subsystem whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-la", "maximum size of the subsystem whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "shift-sites", "location of the leftmost site that belongs to A", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "real-vectors", "indicates  that the ground state vector is real");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension)");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with full.ent extension) to store the reduced density matrix spectrum");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
  (*ToolsGroup) += new BooleanOption  ('\n', "use-svd", "use singular value decomposition instead of diagonalization to compute the entropy");
  (*ToolsGroup) += new SingleDoubleOption  ('\n', "diag-precision", "convergence precision in non LAPACK mode", 1e-7);
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PottsChainEntanglementEntropy -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type PottsChainEntanglementEntropy -h" << endl;
      return -1;
    }

  int ZnValue = Manager.GetInteger("zn");
  int NbrSpins = 0;
  int SzValue = 0;
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  bool SVDFlag = Manager.GetBoolean("use-svd");
  bool RealVectorFlag =  Manager.GetBoolean("real-vectors");

  ComplexVector* GroundStates = 0;
  RealVector* RealGroundStates = 0;
  double* Weights =0;
  bool WeightFlag = false;
  int NbrSpaces = 1;
  char** GroundStateFiles = 0;
  AbstractSpinChain** Spaces = 0;
  int* TotalSz = 0;
  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalSz = new int[1];
      Weights = new double[1];
      Weights[0] = 1.0;
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));      
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
      TotalSz = new int[NbrSpaces];
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	  strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	}
      if (DegeneratedFile.GetNbrColumns() > 1)
	{
	  Weights = DegeneratedFile.GetAsDoubleArray(1);
	  WeightFlag = true;
	}
      else
	{
	  Weights = new double[NbrSpaces];
	  for (int i = 0; i < NbrSpaces; ++i)
	    Weights[i] = 1.0;
	}
    }
  
  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalSz[i] = 0;
      if (PottsFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrSpins, TotalSz[i]) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
    }
  
  if (RealVectorFlag == false)
    {
      GroundStates = new ComplexVector [NbrSpaces];  
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	    {
	      cout << "can't open vector file " << GroundStateFiles[i] << endl;
	      return -1;      
	    }
	}
    }
  else
    {
      RealGroundStates = new RealVector [NbrSpaces];  
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  if (RealGroundStates[i].ReadVector(GroundStateFiles[i]) == false)
	    {
	      cout << "can't open vector file " << GroundStateFiles[i] << endl;
	      return -1;      
	    }
	}
    }
  
  
  ZnValue = Manager.GetInteger("zn");
  
  for (int i = 0; i < NbrSpaces; ++i)
    {
      cout << "Filename: " << GroundStateFiles[i] << " N= " << NbrSpins << " Q= " << TotalSz[i] << endl;
    }
  
  
  Spaces = new AbstractSpinChain* [NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
      switch (ZnValue)
	{
 	case 2 :
 	  Spaces[i] = new Spin1_2ChainFixedParity (NbrSpins, TotalSz[i]);
 	  break;
	case 3 :
	  Spaces[i] = new Potts3Chain (NbrSpins, TotalSz[i], 1000000);
	  break;
	default :
	  {
	    cout << "Zn with n=" << ZnValue  << " is not available" << endl;
	    return -1;
	  }
       }
    }
  
  ofstream File;
  char* EntropyFileName;
  if (Manager.GetString("output-file") != 0)
    {
      EntropyFileName = new char[strlen(Manager.GetString("output-file")) + 1];
      strcpy (EntropyFileName, Manager.GetString("output-file"));
    }
  else
    {
      if (Manager.GetInteger("shift-sites") == 0)
	{
	  EntropyFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "ent");
	  if (EntropyFileName == 0)
	    {
	      cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	      return 0;
	    }
	}
      else
	{
	  char* TmpExtenstion = new char[64];
	  sprintf (TmpExtenstion, "shift_%d.ent", (int) Manager.GetInteger("shift-sites"));
	  EntropyFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", TmpExtenstion);
	  if (EntropyFileName == 0)
	    {
	      cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	      return 0;
	    }
	}
    }
  File.open(EntropyFileName, ios::binary | ios::out);
  char* DensityMatrixFileName = 0;
  if (Manager.GetString("density-matrix") != 0)
    {
      File.open(Manager.GetString("density-matrix"), ios::binary | ios::out);
    }
  else
    {
      DensityMatrixFileName = ReplaceExtensionToFileName(EntropyFileName, "ent", "full.ent");
      if (DensityMatrixFileName == 0)
	{
	  cout << "no ent extension was find in " << EntropyFileName << " file name" << endl;
	  return 0;
	}
   }
  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "# l_a    Q    lambda" << endl;
      DensityMatrixFile.close();
    }
  
  

  File.precision(14);
  cout.precision(14);
  File << "# N_A S_A Tr(\\rho_A)" << endl;

  int SubsystemSize = Manager.GetInteger("min-la");
  if (SubsystemSize < 1)
    SubsystemSize = 1;
  int MeanSubsystemSize = NbrSpins >> 1;
  if (Manager.GetInteger("max-la") > 0)
    {
      MeanSubsystemSize = Manager.GetInteger("max-la");
      if (MeanSubsystemSize > NbrSpins)
	MeanSubsystemSize = NbrSpins;
    }
  for (; SubsystemSize <= MeanSubsystemSize; ++SubsystemSize)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      int MinSzA = 0;
      for (; MinSzA < ZnValue; ++MinSzA)
	{
	  RealDiagonalMatrix TmpDiag;
	  if (RealVectorFlag == false)
	    {
	      if (Manager.GetInteger("shift-sites") != 0)
		{
		  cout << "error, shift-sites not implemented for  complex vectors" << endl;
		}
	      HermitianMatrix PartialDensityMatrix;
	      ComplexMatrix PartialEntanglementMatrix;
	      for (int i = 0; i < NbrSpaces; ++i)
		{
		  cout << "processing subsytem size " << SubsystemSize << " QA=" << MinSzA << endl;
		  if (SVDFlag == false)
		    {
		      HermitianMatrix TmpPartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrix(SubsystemSize, MinSzA, GroundStates[0]);;
		      if (WeightFlag == true)
			TmpPartialDensityMatrix *= Weights[i];
		      if (PartialDensityMatrix.GetNbrRow() == 0)
			PartialDensityMatrix = TmpPartialDensityMatrix;
		      else
			PartialDensityMatrix += TmpPartialDensityMatrix;
		    }
		  else
		    {
		      ComplexMatrix TmpPartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrix(SubsystemSize, MinSzA, GroundStates[0]);
		      if (WeightFlag == true)
			TmpPartialEntanglementMatrix *= sqrt(Weights[i]);
		      if (PartialEntanglementMatrix.GetNbrRow() == 0)
			PartialEntanglementMatrix = TmpPartialEntanglementMatrix;
		      else
			PartialEntanglementMatrix += TmpPartialEntanglementMatrix;
		    }  
		}	
	      if(SVDFlag == false)
		{
		  if ((NbrSpaces > 1) && (WeightFlag == false))
		    PartialDensityMatrix /= ((double) NbrSpaces);
		  
		}
	      else
		{
		  if ((NbrSpaces > 1) && (WeightFlag == false))
		    PartialEntanglementMatrix /= sqrt(((double) NbrSpaces));
		}
	      
	      if ((PartialDensityMatrix.GetNbrRow() >= 1) || ((PartialEntanglementMatrix.GetNbrRow() >= 1) && (PartialEntanglementMatrix.GetNbrColumn() >= 1)))
		{
		  TmpDiag = RealDiagonalMatrix(PartialDensityMatrix.GetNbrRow());
		  if (SVDFlag == false)
		    {
#ifdef __LAPACK__
		      if (LapackFlag == true)
			{
			  PartialDensityMatrix.LapackDiagonalize(TmpDiag);
			  TmpDiag.SortMatrixDownOrder();
			}
		      else
			{
			  PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
			  TmpDiag.SortMatrixDownOrder();
			}
#else
		      PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
		      TmpDiag.SortMatrixDownOrder();
#endif		  		
		    }
		  else 
		    {
		      if ((PartialEntanglementMatrix.GetNbrRow() > 1) && (PartialEntanglementMatrix.GetNbrColumn() > 1))
			{
			  double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
			  int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
			  if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
			    {
			      TmpDimension = PartialEntanglementMatrix.GetNbrRow();
			    }
			  for (int i = 0; i < TmpDimension; ++i)
			    TmpValues[i] *= TmpValues[i];
			  TmpDiag = RealDiagonalMatrix(TmpValues, TmpDimension);
			  TmpDiag.SortMatrixDownOrder();
			}
		      else
			{
			  double TmpValue = 0.0;
			  if (PartialEntanglementMatrix.GetNbrRow() == 1)
			    {
			      for (int i = 0; i < PartialEntanglementMatrix.GetNbrColumn(); ++i)
				{ 
				  Complex Tmp = PartialEntanglementMatrix[i][0]; 
				  TmpValue += Tmp.Re * Tmp.Re + Tmp.Im * Tmp.Im;
				}
			    }
			  else
			    {
			      for (int i = 0; i < PartialEntanglementMatrix.GetNbrRow(); ++i)
				{
				  Complex Tmp = PartialEntanglementMatrix[0][i]; 
				  TmpValue += Tmp.Re * Tmp.Re + Tmp.Im * Tmp.Im;		
				}		  
			    }
			  TmpDiag = RealDiagonalMatrix(1, 1);
			  TmpDiag[0] = TmpValue;
			}
		    }  		
		}
	      else
		{
		  if (PartialDensityMatrix.GetNbrRow() == 1)
		    {
		      TmpDiag = RealDiagonalMatrix(PartialDensityMatrix.GetNbrRow());
		      TmpDiag[1] = PartialDensityMatrix(0,0);
		    }
		}
	    }
	  else
	    {
	      RealSymmetricMatrix PartialDensityMatrix;
	      RealMatrix PartialEntanglementMatrix;
	      for (int i = 0; i < NbrSpaces; ++i)
		{
		  cout << "processing subsytem size " << SubsystemSize << " SzA=" << MinSzA << endl;
		  if (SVDFlag == false)
		    {
		      RealSymmetricMatrix TmpPartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrix(SubsystemSize, MinSzA, Manager.GetInteger("shift-sites"), RealGroundStates[0]);;
		      if (WeightFlag == true)
			TmpPartialDensityMatrix *= Weights[i];
		      if (PartialDensityMatrix.GetNbrRow() == 0)
			PartialDensityMatrix = TmpPartialDensityMatrix;
		      else
			PartialDensityMatrix += TmpPartialDensityMatrix;
		    }
		  else 
		    {
		      RealMatrix TmpPartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrix(SubsystemSize, MinSzA, Manager.GetInteger("shift-sites"), RealGroundStates[0]);
		      if (WeightFlag == true)
			TmpPartialEntanglementMatrix *= sqrt(Weights[i]);
		      if (PartialEntanglementMatrix.GetNbrRow() == 0)
			PartialEntanglementMatrix = TmpPartialEntanglementMatrix;
		      else
			PartialEntanglementMatrix += TmpPartialEntanglementMatrix;
		    }  
		}
	      if(SVDFlag == false)
		{
		  if ((NbrSpaces > 1) && (WeightFlag == false))
		    PartialDensityMatrix /= ((double) NbrSpaces);
		}
	      else
		{
		  if ((NbrSpaces > 1) && (WeightFlag == false))
		    PartialEntanglementMatrix /= sqrt(((double) NbrSpaces));
		}
	      
	      if ((PartialDensityMatrix.GetNbrRow() > 1) || ((PartialEntanglementMatrix.GetNbrRow() >= 1) && (PartialEntanglementMatrix.GetNbrColumn() >= 1)))
		{
		  TmpDiag = RealDiagonalMatrix(PartialDensityMatrix.GetNbrRow());
		  if (SVDFlag == false)
		    {
#ifdef __LAPACK__
		      if (LapackFlag == true)
			{
			  PartialDensityMatrix.LapackDiagonalize(TmpDiag);
			  TmpDiag.SortMatrixDownOrder();
			}
		      else
			{
			  PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
			  TmpDiag.SortMatrixDownOrder();
			}
#else
		      PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
		      TmpDiag.SortMatrixDownOrder();
#endif		  		
		    }
		  else 
		    {
		      if ((PartialEntanglementMatrix.GetNbrRow() > 1) && (PartialEntanglementMatrix.GetNbrColumn() > 1))
			{
			  double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
			  int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
			  if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
			    {
			      TmpDimension = PartialEntanglementMatrix.GetNbrRow();
			    }
			  for (int i = 0; i < TmpDimension; ++i)
			    TmpValues[i] *= TmpValues[i];
			  TmpDiag = RealDiagonalMatrix(TmpValues, TmpDimension);
			  TmpDiag.SortMatrixDownOrder();
			}
		      else
			{
			  double TmpValue = 0.0;
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
			  TmpDiag = RealDiagonalMatrix(1, 1);
			  TmpDiag[0] = TmpValue;
			}
		    }
		}  
	      else
		{
		  if (PartialDensityMatrix.GetNbrRow() == 1)
		    {
		      TmpDiag = RealDiagonalMatrix(PartialDensityMatrix.GetNbrRow());
		      TmpDiag[0] = PartialDensityMatrix(0,0);
		    }
		}
	    }
	  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
	    {
	      if (TmpDiag[i] > 1e-14)
		{
		  EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
		  DensitySum += TmpDiag[i];
		}
	    }
	  if (DensityMatrixFileName != 0)
	    {
	      ofstream DensityMatrixFile;
	      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
	      DensityMatrixFile.precision(14);
	      for (int i = 0; i <TmpDiag.GetNbrRow(); ++i)
		DensityMatrixFile << SubsystemSize << " " << MinSzA << " " << TmpDiag[i] << endl;
	      DensityMatrixFile.close();
	    }	 
	}
      File << SubsystemSize << " " << (-EntanglementEntropy) << " " << DensitySum << endl;
    }
  File.close();
  delete[] EntropyFileName;
  return 0;
}
