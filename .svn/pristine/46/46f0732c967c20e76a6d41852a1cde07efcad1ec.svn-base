#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

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

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"

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
  OptionManager Manager ("FQHETopInsulatorEntanglementEntropyMomentumSpace" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "size-x", "number of sites along the x direction that belong to the rectangular subsytem A (0 if equal to the total system size along x)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "size-y", "number of sites along the y direction that belong to the rectangular subsytem A (0 if equal to the total system size along y)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "start-x", "x momentum marking the beginning of the rectangluar subsystem A", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "start-y", "y momentum marking the beginning of the rectangluar subsystem A", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with moment extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETopInsulatorEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSpaces = 1;
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  FermionOnSquareLatticeMomentumSpace** Spaces = 0;
  ComplexVector* GroundStates = 0;
  char** GroundStateFiles = 0;
  int* TotalKx = 0;
  int* TotalKy = 0;
  int NbrParticles = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  bool Statistics = true;
  double* Coefficients = 0;
  bool ShowTimeFlag = Manager.GetBoolean("show-time");

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHETopInsulatorEntanglementEntropyParticlePartition -h" << endl;
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


  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalKx = new int[1];
      TotalKy = new int[1];
      Coefficients = new double[1];
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));
      Coefficients[0] = 1.0;
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
       TotalKx = new int[NbrSpaces];
       TotalKy = new int[NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	 }
       if (DegeneratedFile.GetNbrColumns() == 1)
	 {
	   Coefficients = new double[NbrSpaces];
	   for (int i = 0; i < NbrSpaces; ++i)
	     Coefficients[i] = 1.0 / ((double) NbrSpaces);
	 }
       else
	 {
	   double TmpSum = 0.0;
	   Coefficients = DegeneratedFile.GetAsDoubleArray(1);
	   for (int i = 0; i < NbrSpaces; ++i)
	     TmpSum += Coefficients[i];
	   TmpSum = 1.0 / TmpSum;
	   for (int i = 0; i < NbrSpaces; ++i)
	     Coefficients[i] *= TmpSum;
	 }
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalKx[i] = 0;
      TotalKy[i] = 0;
      double Mass = 0.0;
      if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[i],
							      NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], Mass, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
    }


  GroundStates = new ComplexVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
      {
	cout << "can't open vector file " << GroundStateFiles[i] << endl;
	return -1;      
      }

  Spaces = new FermionOnSquareLatticeMomentumSpace*[NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
      Spaces[i] = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
      
      if (Spaces[i]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	{
	  cout << "dimension mismatch between Hilbert space and ground state" << endl;
	  return 0;
	}
    }


  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "#  N    Kx    Ky    lambda";
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }

  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "moment");
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

  int SubsystemSizeX = Manager.GetInteger("size-x");
  int SubsystemSizeY = Manager.GetInteger("size-y");
  int SubsystemStartX = Manager.GetInteger("start-x");
  int SubsystemStartY = Manager.GetInteger("start-y");

  double EntanglementEntropy = 0.0;
  double DensitySum = 0.0;
  for (int SubsystemNbrParticles = 0; SubsystemNbrParticles <= NbrParticles; ++SubsystemNbrParticles)
    {
      for (int SubsystemTotalKx = 0; SubsystemTotalKx < NbrSiteX; ++SubsystemTotalKx)
	{
	  for (int SubsystemTotalKy = 0; SubsystemTotalKy < NbrSiteY; ++SubsystemTotalKy)
	    {
	  
	      cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Kx=" << SubsystemTotalKx << " Ky=" << SubsystemTotalKy << endl;
	      timeval TotalStartingTime;
	      timeval TotalEndingTime;
	      if (ShowTimeFlag == true)
		{
		  gettimeofday (&(TotalStartingTime), 0);
		}
	      HermitianMatrix PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixMomentumSpace(SubsystemSizeX, SubsystemSizeY, SubsystemStartX, SubsystemStartY, SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStates[0]);
	      PartialDensityMatrix *= Coefficients[0];
	      for (int i = 1; i < NbrSpaces; ++i)
		{
		  HermitianMatrix TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixMomentumSpace(SubsystemSizeX, SubsystemSizeY, SubsystemStartX, SubsystemStartY, SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStates[i]);
		  TmpMatrix *= Coefficients[i];
		  PartialDensityMatrix += TmpMatrix;
		}
	      if (ShowTimeFlag == true)
		{
		  gettimeofday (&(TotalEndingTime), 0);
		  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		  cout << "reduced density matrix evaluated in " << Dt << "s" << endl;
		}
	      if (PartialDensityMatrix.GetNbrRow() > 1)
		{
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalStartingTime), 0);
		    }
		  RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
		  if (LapackFlag == true)
		    {
		      PartialDensityMatrix.LapackDiagonalize(TmpDiag);
		    }
		  else
		    {
		      PartialDensityMatrix.Diagonalize(TmpDiag);
		    }
#else
		  PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
		  TmpDiag.SortMatrixDownOrder();
		  if (DensityMatrixFileName != 0)
		    {
		      ofstream DensityMatrixFile;
		      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		      DensityMatrixFile.precision(14);
		      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << TmpDiag[i] << endl;
		      DensityMatrixFile.close();
		    }
		  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		    {
		      if (TmpDiag[i] > 1e-14)
			{
			  EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
			  DensitySum +=TmpDiag[i];
			}
		    }
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalEndingTime), 0);
		      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		      cout << "diagonalization done in " << Dt << "s" << endl;
		    }
		}
	      else
		if (PartialDensityMatrix.GetNbrRow() == 1)
		  {
		    double TmpValue = PartialDensityMatrix(0,0);
		    if (DensityMatrixFileName != 0)
		      {
			ofstream DensityMatrixFile;
			DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			DensityMatrixFile.precision(14);
			DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << TmpValue << endl;
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
  File << SubsystemSizeX << " " << SubsystemSizeY << " " << SubsystemStartX << " " << SubsystemStartY << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << endl;
  File.close();
  return 0;
}
