#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Hamiltonian/SpinChainHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainFull.h"
#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin2Chain.h"
#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzInversionSymmetries.h"
#include "HilbertSpace/Spin2ChainWithTranslations.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndSzInversionSymmetries.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

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
  OptionManager Manager ("SpinChainEntanglementEntropy" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 1);
  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");  
  (*SystemGroup) += new BooleanOption  ('\n', "no-sz", "ground state does not have well-defined Sz quantum number (e.g., Ising or XYZ models)");  
  (*SystemGroup) += new BooleanOption  ('c', "complex", "consider complex wave function");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-la", "minimum size of the subsystem whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-la", "maximum size of the subsystem whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*OutputGroup) += new BooleanOption ('\n', "disable-densitymatrix", "do not save the eigenvalues of the reduced density matrix");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
  (*ToolsGroup) += new BooleanOption  ('\n', "use-svd", "use singular value decomposition instead of diagonalization to compute the entropy");
  (*ToolsGroup) += new SingleDoubleOption  ('\n', "diag-precision", "convergence precision in non LAPACK mode", 1e-7);
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinChainEntanglementEntropy -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type SpinChainEntanglementEntropy -h" << endl;
      return -1;
    }

  int SpinValue = 0;
  int NbrSpins = 0;
  int SzValue = 0;
  bool NoSzFlag = Manager.GetBoolean("no-sz");
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  bool SVDFlag = Manager.GetBoolean("use-svd");
  char* DensityMatrixFileName = 0;
  double* Weights =0;
  bool WeightFlag = false;
  int NbrSpaces = 1;
  char** GroundStateFiles = 0;
  AbstractSpinChain** Spaces = 0;
  int* TotalSz = 0;
  int* Momenta = 0;
  int* InversionSectors = 0;
  int* SzSymmetrySectors = 0;

  ofstream File;

  if (NoSzFlag)
  {
    //no Sz symmetry
    cout << "Assume no Sz conservation "<<endl;

    GroundStateFiles = new char* [1];
    Momenta = new int[1];
    Weights = new double[1];
    Weights[0] = 1.0;
    GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
    strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));      

    NbrSpaces= 1;
    if (SpinAllSzFindSystemInfoFromVectorFileName(GroundStateFiles[0], NbrSpins, SpinValue, Momenta[0]) == false)
		{
			cout << "error while retrieving system parameters from file name " << GroundStateFiles[0] << endl;
		  	return -1;
	    }
	cout << "Read in file " << GroundStateFiles[0] << " NbrSpins= "<<NbrSpins<<" SpinValue= "<<SpinValue<<" Momenta= "<<Momenta[0]<<endl;	

  if (Manager.GetString("output-file") != 0)
    {
      File.open(Manager.GetString("output-file"), ios::binary | ios::out);
      if (Manager.GetBoolean("disable-densitymatrix") == false)
	{
	  DensityMatrixFileName  = ReplaceExtensionToFileName(Manager.GetString("output-file"), "ent", "full.ent");
	  if (DensityMatrixFileName == 0)
	    {
	      cout << "no ent extension was find in " << Manager.GetString("output-file") << " file name" << endl;
	      return 0;
	    }
	}
    }
  else
    {
      char* TmpFileName;
      TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "ent");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      if (Manager.GetBoolean("disable-densitymatrix") == false)
	{
	  DensityMatrixFileName  = ReplaceExtensionToFileName(TmpFileName, "ent", "full.ent");
	  if (DensityMatrixFileName == 0)
	    {
	      cout << "no ent extension was find in " <<  TmpFileName << " file name" << endl;
	      return 0;
	    }
	}
      delete[] TmpFileName;
    }


  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
	  DensityMatrixFile << "# l_a    lambda" << endl;
      DensityMatrixFile.close();
    }
  
  File.precision(14);
  cout.precision(14);

  cout << "Complex problem "<<endl;
  ComplexVector* GroundStates = 0;
     
  GroundStates = new ComplexVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    {
	 if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	   {
	     cout << "can't open vector file " << GroundStateFiles[i] << endl;
	     return -1;      
	   }
    }
 
  Spaces = new AbstractSpinChain* [NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
	{
	   switch (SpinValue)
		{
		case 1 :
		  Spaces[i] = new Spin1_2ChainWithTranslations (NbrSpins, Momenta[i], 1, 1000000, 1000000);
		  break;
		default :
		  {
		    if ((SpinValue & 1) == 0)
		      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		    else 
		      cout << "spin " << SpinValue << "/2 are not available" << endl;
		    return -1;
		  }
		}
	}

  for (int i = 0; i < NbrSpaces; ++i)
	{
	  if (Spaces[i]->GetHilbertSpaceDimension() != GroundStates[i].GetVectorDimension())
	    {
	      cout << "error, dimension mismatch for " << GroundStateFiles[i] << " (dimnension is " << GroundStates[i].GetVectorDimension() 
		   << ", should be " << Spaces[i]->GetHilbertSpaceDimension() << ")" << endl;
	      return -1;
	    }
	}

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

	 ComplexMatrix PartialEntanglementMatrix;
	 for (int i = 0; i < NbrSpaces; ++i)
	   {
		ComplexMatrix TmpPartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrix(SubsystemSize, GroundStates[0]);
	    if (WeightFlag == true)
			   TmpPartialEntanglementMatrix *= sqrt(Weights[i]);
		if (PartialEntanglementMatrix.GetNbrRow() == 0)
			   	PartialEntanglementMatrix = TmpPartialEntanglementMatrix;
			 else
			 	PartialEntanglementMatrix += TmpPartialEntanglementMatrix;
		}	 
	  if ((NbrSpaces > 1) && (WeightFlag == false))
		   PartialEntanglementMatrix /= sqrt(((double) NbrSpaces));
	    
	     
	  if (((PartialEntanglementMatrix.GetNbrRow() >= 1) && (PartialEntanglementMatrix.GetNbrColumn() >= 1)))
	       {
	       	  int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
			  if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
			   	{
			      TmpDimension = PartialEntanglementMatrix.GetNbrRow();
			   	}
		    RealDiagonalMatrix TmpDiag (TmpDimension);
	
		    PartialEntanglementMatrix.RemoveZeroColumns();
		    PartialEntanglementMatrix.RemoveZeroRows();
		    if ((PartialEntanglementMatrix.GetNbrRow() > 1) && (PartialEntanglementMatrix.GetNbrColumn() > 1))
		       {
			 	cout << "PartialEntanglementMatrix = " << PartialEntanglementMatrix.GetNbrRow() << " x " << PartialEntanglementMatrix.GetNbrColumn() << endl;
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
			   DensityMatrixFile << SubsystemSize << " " << TmpDiag[i] << endl;
		     DensityMatrixFile.close();
		   }
	     }
	     else
	       {
		 if ((PartialEntanglementMatrix.GetNbrRow() == 1) || (PartialEntanglementMatrix.GetNbrColumn() == 1))
		   {
		     double TmpValue = 0; //PartialDensityMatrix(0,0);
		     if (TmpValue > 1e-14)
		       {
			 EntanglementEntropy += TmpValue * log(TmpValue);
			 DensitySum += TmpValue;
		       }
		     if (DensityMatrixFileName != 0)
		       {
			 ofstream DensityMatrixFile;
			 DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			 DensityMatrixFile.precision(14);
			 DensityMatrixFile << SubsystemSize << " "<< TmpValue << endl;
			 DensityMatrixFile.close();
		       }		  
		   }
	       }
	 File << SubsystemSize << " " << (-EntanglementEntropy) << " " << DensitySum << endl;
     }
  File.close();
  return 0;
 }


  
  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalSz = new int[1];
      Momenta = new int[1];
      InversionSectors = new int[1];
      SzSymmetrySectors = new int[1];
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
      Momenta = new int[NbrSpaces];
      InversionSectors = new int[NbrSpaces];
      SzSymmetrySectors = new int[NbrSpaces];
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


  bool SzFlag = true;
  bool MomentumFlag = true;
  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalSz[i] = 0;
      Momenta[i] = 0;
      InversionSectors[i] = 0;
      SzSymmetrySectors[i] = 0;
      if (SpinFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrSpins, TotalSz[i], SpinValue, Momenta[i],
					       InversionSectors[i], SzSymmetrySectors[i]) == false)
	{
	  MomentumFlag = false;
	  if (SpinFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrSpins, TotalSz[i], SpinValue) == false)
	    {
	      SzFlag = false;
	      if (SpinFindSystemInfoFromFileName(GroundStateFiles[i], NbrSpins, SpinValue) == false)
		{	     
		  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
		  return -1;
		}
	    }
	}
    }
  
  //ofstream File;
  if (Manager.GetString("output-file") != 0)
    {
      File.open(Manager.GetString("output-file"), ios::binary | ios::out);
      if (Manager.GetBoolean("disable-densitymatrix") == false)
	{
	  DensityMatrixFileName  = ReplaceExtensionToFileName(Manager.GetString("output-file"), "ent", "full.ent");
	  if (DensityMatrixFileName == 0)
	    {
	      cout << "no ent extension was find in " << Manager.GetString("output-file") << " file name" << endl;
	      return 0;
	    }
	}
    }
  else
    {
      char* TmpFileName;
      TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "ent");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      if (Manager.GetBoolean("disable-densitymatrix") == false)
	{
	  DensityMatrixFileName  = ReplaceExtensionToFileName(TmpFileName, "ent", "full.ent");
	  if (DensityMatrixFileName == 0)
	    {
	      cout << "no ent extension was find in " <<  TmpFileName << " file name" << endl;
	      return 0;
	    }
	}
      delete[] TmpFileName;
    }
  
  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      if (SzFlag == true)
	DensityMatrixFile << "# l_a    Sz    lambda" << endl;
      else
	DensityMatrixFile << "# l_a    lambda" << endl;
      DensityMatrixFile.close();
    }
  
  File.precision(14);
  cout.precision(14);
      
  if (Manager.GetBoolean("complex") == false)
    {
      RealVector* GroundStates = 0;
      
      GroundStates = new RealVector [NbrSpaces];  
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	    {
	      cout << "can't open vector file " << GroundStateFiles[i] << endl;
	      return -1;      
	    }
	}
      
      
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  cout << "Filename: " << GroundStateFiles[i] << " N= " << NbrSpins << " Sz= " << TotalSz[i] << " spin= " << SpinValue << endl;
	}
      
      Spaces = new AbstractSpinChain* [NbrSpaces];
      if (MomentumFlag == false)
	{
	  for (int i = 0; i < NbrSpaces; ++i)
	    {
	      switch (SpinValue)
		{
		case 1 :
		  Spaces[i] = new Spin1_2Chain (NbrSpins, TotalSz[i], 1000000);
		  break;
		case 2 :
		  Spaces[i] = new Spin1Chain (NbrSpins, TotalSz[i], 1000000);
		  break;
		case 4 :
		  Spaces[i] = new Spin2Chain (NbrSpins, TotalSz[i], 1000000);
		  break;
		default :
		  {
		    if ((SpinValue & 1) == 0)
		      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		    else 
		      cout << "spin " << SpinValue << "/2 are not available" << endl;
		    return -1;
		  }
		}
	    }
	}
      else
	{
	  for (int i = 0; i < NbrSpaces; ++i)
	    {
	      switch (SpinValue)
		{
		case 1 :
		  Spaces[i] = new Spin1_2ChainWithTranslations (NbrSpins, Momenta[i], 1, TotalSz[i], 1000000, 1000000);
		  break;
		case 2 :
		  {
		    if (InversionSectors[i] != 0)
		      {
			if (SzSymmetrySectors[i] != 0)
			  {
			    Spaces[i] = new Spin1ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, Momenta[i], InversionSectors[i], SzSymmetrySectors[i], TotalSz[i]);
			  }
			else
			  {
			    Spaces[i] = new Spin1ChainWithTranslationsAndInversionSymmetry (NbrSpins, Momenta[i], InversionSectors[i], TotalSz[i]);
			  }
		      }
		    else
		      {
			if (SzSymmetrySectors[i] != 0)
			  {
			    Spaces[i] = new Spin1ChainWithTranslationsAndSzSymmetry (NbrSpins, Momenta[i], SzSymmetrySectors[i], TotalSz[i]);
			  }
			else
			  {
			    Spaces[i] = new Spin1ChainWithTranslations (NbrSpins, Momenta[i], TotalSz[i]);
			  }
		      }
		  }
		  break;
		case 4 :
		  {
		    if (InversionSectors[i] != 0)
		      {
			if (SzSymmetrySectors[i] != 0)
			  {
			    Spaces[i] = new Spin2ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, Momenta[i], InversionSectors[i], SzSymmetrySectors[i], TotalSz[i]);
			  }
			else
			  {
			    Spaces[i] = new Spin2ChainWithTranslationsAndInversionSymmetry (NbrSpins, Momenta[i], InversionSectors[i], TotalSz[i]);
			  }
		      }
		    else
		      {
			if (SzSymmetrySectors[i] != 0)
			  {
			    Spaces[i] = new Spin2ChainWithTranslationsAndSzSymmetry (NbrSpins, Momenta[i], SzSymmetrySectors[i], TotalSz[i]);
			  }
			else
			  {
			    Spaces[i] = new Spin2ChainWithTranslations (NbrSpins, Momenta[i], TotalSz[i]);
			  }
		      }
		  }
		  break;
		default :
		  {
		    if ((SpinValue & 1) == 0)
		      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		    else 
		      cout << "spin " << SpinValue << "/2 are not available" << endl;
		    return -1;
		  }
		}
	    }
	}
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  if (Spaces[i]->GetHilbertSpaceDimension() != GroundStates[i].GetVectorDimension())
	    {
	      cout << "error, dimension mismatch for " << GroundStateFiles[i] << " (dimnension is " << GroundStates[i].GetVectorDimension() 
		   << ", should be " << Spaces[i]->GetHilbertSpaceDimension() << ")" << endl;
	      return -1;
	    }
	}
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
	  int MaxSzA = (SubsystemSize * SpinValue);
	  int MinSzA = -MaxSzA;
	  int MaxSzB = ((NbrSpins - SubsystemSize) * SpinValue);
	  int MinSzB = -MaxSzB;
	  for (; MinSzA <= MaxSzA; MinSzA += 2)
	    {
	      RealSymmetricMatrix PartialDensityMatrix;
	      RealMatrix PartialEntanglementMatrix;
	      for (int i = 0; i < NbrSpaces; ++i)
		{
		  int SzB = SzValue - MinSzA;
		  if ((SzB <= MaxSzB) && (SzB >= MinSzB))
		    {
		      cout << "processing subsytem size " << SubsystemSize << " SzA=" << MinSzA << endl;
		      if (SVDFlag == false)
			{
			  RealSymmetricMatrix TmpPartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrix(SubsystemSize, MinSzA, GroundStates[0]);
			  if (WeightFlag == true)
			    TmpPartialDensityMatrix *= Weights[i];
			  if (PartialDensityMatrix.GetNbrRow() == 0)
			    PartialDensityMatrix = TmpPartialDensityMatrix;
			  else
			    PartialDensityMatrix += TmpPartialDensityMatrix;
			}
		      else //use SVD
			{
			  RealMatrix TmpPartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrix(SubsystemSize, MinSzA, GroundStates[0]);
			  if (WeightFlag == true)
			    TmpPartialEntanglementMatrix *= sqrt(Weights[i]);
			  if (PartialEntanglementMatrix.GetNbrRow() == 0)
			    PartialEntanglementMatrix = TmpPartialEntanglementMatrix;
			  else
			    PartialEntanglementMatrix += TmpPartialEntanglementMatrix;
                   }  
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
		  RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
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
		  else //SVD...
		    {
		      if ((PartialEntanglementMatrix.GetNbrRow() > 1) && (PartialEntanglementMatrix.GetNbrColumn() > 1))
			{
			  cout << "PartialEntanglementMatrix = " << PartialEntanglementMatrix.GetNbrRow() << " x " << PartialEntanglementMatrix.GetNbrColumn() << endl;
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
	      else
		{
		  if (PartialDensityMatrix.GetNbrRow() == 1)
		    {
		      double TmpValue = PartialDensityMatrix(0,0);
		      if (TmpValue > 1e-14)
			{
			  EntanglementEntropy += TmpValue * log(TmpValue);
			  DensitySum += TmpValue;
			}
		      if (DensityMatrixFileName != 0)
			{
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			  DensityMatrixFile << SubsystemSize << " " << MinSzA << " " << TmpValue << endl;
			  DensityMatrixFile.close();
			}		  
		    }
		}
	    }
	  File << SubsystemSize << " " << (-EntanglementEntropy) << " " << DensitySum << endl;
	}
      File.close();
    }
  else //Complex problem
    {
     cout << "Complex problem "<<endl;
     ComplexVector* GroundStates = 0;
     
     GroundStates = new ComplexVector [NbrSpaces];  
     for (int i = 0; i < NbrSpaces; ++i)
       {
	 if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	   {
	     cout << "can't open vector file " << GroundStateFiles[i] << endl;
	     return -1;      
	   }
       }
     
     
     if (SzFlag == true)
       {
	 for (int i = 0; i < NbrSpaces; ++i)
	   {
	     cout << "Filename: " << GroundStateFiles[i] << " N= " << NbrSpins << " Sz= " << TotalSz[i] << " spin= " << SpinValue << endl;
	   }
       }
     else
       {
	 for (int i = 0; i < NbrSpaces; ++i)
	   {
	     cout << "Filename: " << GroundStateFiles[i] << " N= " << NbrSpins << " spin= " << SpinValue << endl;
	   }
       }
     
     
     Spaces = new AbstractSpinChain* [NbrSpaces];
     
      if (MomentumFlag == false)
	{
	  if (SzFlag == true)
	    {
	      for (int i = 0; i < NbrSpaces; ++i)
		{
		  switch (SpinValue)
		    {
		    case 1 :
		      Spaces[i] = new Spin1_2Chain (NbrSpins, TotalSz[i], 1000000);
		      break;
		    case 2 :
		      Spaces[i] = new Spin1Chain (NbrSpins, TotalSz[i], 1000000);
		      break;
		    default :
		      {
			if ((SpinValue & 1) == 0)
			  cout << "spin " << (SpinValue / 2) << " are not available" << endl;
			else 
			  cout << "spin " << SpinValue << "/2 are not available" << endl;
			return -1;
		      }
		    }
		}
	    }
	  else
	    {
	      for (int i = 0; i < NbrSpaces; ++i)
		{
		  switch (SpinValue)
		    {
		    case 1 :
		      Spaces[i] = new Spin1_2ChainFull (NbrSpins);
		      break;
		    }
		}
	    }
	}
      else
	{
	  for (int i = 0; i < NbrSpaces; ++i)
	    {
	      switch (SpinValue)
		{
		case 1 :
		  Spaces[i] = new Spin1_2ChainWithTranslations (NbrSpins, Momenta[i], 1, TotalSz[i], 1000000, 1000000);
		  break;
		case 2 :
		  {
		    if (InversionSectors[i] != 0)
		      {
			if (SzSymmetrySectors[i] != 0)
			  {
			    Spaces[i] = new Spin1ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, Momenta[i], InversionSectors[i], SzSymmetrySectors[i], TotalSz[i]);
			  }
			else
			  {
			    Spaces[i] = new Spin1ChainWithTranslationsAndInversionSymmetry (NbrSpins, Momenta[i], InversionSectors[i], TotalSz[i]);
			  }
		      }
		    else
		      {
			if (SzSymmetrySectors[i] != 0)
			  {
			    Spaces[i] = new Spin1ChainWithTranslationsAndSzSymmetry (NbrSpins, Momenta[i], SzSymmetrySectors[i], TotalSz[i]);
			  }
			else
			  {
			    Spaces[i] = new Spin1ChainWithTranslations (NbrSpins, Momenta[i], TotalSz[i]);
			  }
		      }
		  }
		  break;
		case 4 :
		  {
		    if (InversionSectors[i] != 0)
		      {
			if (SzSymmetrySectors[i] != 0)
			  {
			    Spaces[i] = new Spin2ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, Momenta[i], InversionSectors[i], SzSymmetrySectors[i], TotalSz[i]);
			  }
			else
			  {
			    Spaces[i] = new Spin2ChainWithTranslationsAndInversionSymmetry (NbrSpins, Momenta[i], InversionSectors[i], TotalSz[i]);
			  }
		      }
		    else
		      {
			if (SzSymmetrySectors[i] != 0)
			  {
			    Spaces[i] = new Spin2ChainWithTranslationsAndSzSymmetry (NbrSpins, Momenta[i], SzSymmetrySectors[i], TotalSz[i]);
			  }
			else
			  {
			    Spaces[i] = new Spin2ChainWithTranslations (NbrSpins, Momenta[i], TotalSz[i]);
			  }
		      }
		  }
		  break;
		default :
		  {
		    if ((SpinValue & 1) == 0)
		      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		    else 
		      cout << "spin " << SpinValue << "/2 are not available" << endl;
		    return -1;
		  }
		}
	    }
	}
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  if (Spaces[i]->GetHilbertSpaceDimension() != GroundStates[i].GetVectorDimension())
	    {
	      cout << "error, dimension mismatch for " << GroundStateFiles[i] << " (dimnension is " << GroundStates[i].GetVectorDimension() 
		   << ", should be " << Spaces[i]->GetHilbertSpaceDimension() << ")" << endl;
	      return -1;
	    }
	}


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
	 int MaxSzA = (SubsystemSize * SpinValue);
	 int MinSzA = -MaxSzA;
	 int MaxSzB = ((NbrSpins - SubsystemSize) * SpinValue);
	 int MinSzB = -MaxSzB;
	 if (SzFlag == false)
	   {
	     MaxSzA = 0;
	     MinSzA = 0;
	     MaxSzB = 0;
	     MinSzB = 0;
	     SzValue = 0;
	   }
	 for (; MinSzA <= MaxSzA; MinSzA += 2)
	   {
	     HermitianMatrix PartialDensityMatrix;
	     ComplexMatrix PartialEntanglementMatrix;
	     for (int i = 0; i < NbrSpaces; ++i)
	       {
		 int SzB = SzValue - MinSzA;
		 if ((SzB <= MaxSzB) && (SzB >= MinSzB))
		   {
		     if (SzFlag == true)
		       {
			 cout << "processing subsytem size " << SubsystemSize << " SzA=" << MinSzA << endl;
		       }
		     else
		       {
			 cout << "processing subsytem size " << SubsystemSize << endl;
		       }
		     if (SVDFlag == false)
		       {
			 cout<<"Use SVD "<<endl;
			 HermitianMatrix TmpPartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrix(SubsystemSize, MinSzA, GroundStates[0]);;
			 if (WeightFlag == true)
			   TmpPartialDensityMatrix *= Weights[i];
			 if (PartialDensityMatrix.GetNbrRow() == 0)
			   PartialDensityMatrix = TmpPartialDensityMatrix;
			 else
			   PartialDensityMatrix += TmpPartialDensityMatrix;
		       }
		     else //use SVD
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
		 RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
		 if (SVDFlag == false)
		   {
		     cout<<"Use SVD "<<endl;
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
		 else //SVD...
		   {
		     PartialEntanglementMatrix.RemoveZeroColumns();
		     PartialEntanglementMatrix.RemoveZeroRows();
		     if ((PartialEntanglementMatrix.GetNbrRow() > 1) && (PartialEntanglementMatrix.GetNbrColumn() > 1))
		       {
			 cout << "PartialEntanglementMatrix = " << PartialEntanglementMatrix.GetNbrRow() << " x " << PartialEntanglementMatrix.GetNbrColumn() << endl;
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
		     if (SzFlag == true)
		       {
			 for (int i = 0; i <TmpDiag.GetNbrRow(); ++i)
			   DensityMatrixFile << SubsystemSize << " " << MinSzA << " " << TmpDiag[i] << endl;
		       }
		     else
		       {
			 for (int i = 0; i <TmpDiag.GetNbrRow(); ++i)
			   DensityMatrixFile << SubsystemSize << " " << TmpDiag[i] << endl;
		       }
		     DensityMatrixFile.close();
		   }
	       }
	     else
	       {
		 if (PartialDensityMatrix.GetNbrRow() == 1)
		   {
		     double TmpValue = PartialDensityMatrix(0,0);
		     if (TmpValue > 1e-14)
		       {
			 EntanglementEntropy += TmpValue * log(TmpValue);
			 DensitySum += TmpValue;
		       }
		     if (DensityMatrixFileName != 0)
		       {
			 ofstream DensityMatrixFile;
			 DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			 DensityMatrixFile.precision(14);
			 DensityMatrixFile << SubsystemSize << " " << MinSzA << " " << TmpValue << endl;
			 DensityMatrixFile.close();
		       }		  
		   }
	       }
	   }
	 File << SubsystemSize << " " << (-EntanglementEntropy) << " " << DensitySum << endl;
       }
  File.close();
  }
  return 0;
}
