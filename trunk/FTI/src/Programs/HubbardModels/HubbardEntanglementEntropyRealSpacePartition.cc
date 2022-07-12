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

#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"


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
  OptionManager Manager ("HubbardEntanglementEntropyRealSpacePartition" , "0.01");
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
  (*SystemGroup) += new SingleStringOption  ('\n', "kept-sites", "column-based tesxt file that list sites that have to be kept");
  (*SystemGroup) += new BooleanOption  ('\n', "no-spin", "use this option for spinless particles");
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*OutputGroup) += new BooleanOption ('\n', "disable-densitymatrix", "do not save the eigenvalues of the reduced density matrix");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardEntanglementEntropyRealSpacePartition -h" << endl;
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
  char* EntropyFileName = 0;
  char* DensityMatrixFileName = 0;

  ComplexVector* GroundStates = 0;
  char** GroundStateFiles = 0;
  int NbrParticles = 0;
  int NbrSites = 0;
  bool Statistics = true;
  bool GutzwillerFlag = false;
  bool ShowTimeFlag = Manager.GetBoolean("show-time");
  double* Coefficients = 0;
  bool SpinlessFlag = Manager.GetBoolean("no-spin");
  bool FixedTotalSzFlag = false;
  bool TwoDimensionTranslationFlag = false;
  int TotalSz = 0;
  int XMomentum = 0;
  int YMomentum = 0;
  int XPeriodicity = 0;
  int YPeriodicity = 0;

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type HubbardEntanglementEntropyRealSpacePartition -h" << endl;
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

  if (Manager.GetString("kept-sites") == 0)
    {
      cout << "error, a file describing the sites to keep has to be provided" << endl;
      return -1;     
    }

  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
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
      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrParticles, NbrSites, TotalSz, 
									    XMomentum, YMomentum, XPeriodicity, YPeriodicity, Statistics, GutzwillerFlag) == false)
	{
	  if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrParticles, NbrSites,
										XMomentum, YMomentum, XPeriodicity, YPeriodicity, Statistics, GutzwillerFlag) == false)
	    {
	      if (FTIHubbardModelFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrParticles, NbrSites, TotalSz, Statistics, GutzwillerFlag) == false)
		{
		  if (FTIHubbardModelFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrParticles, NbrSites, Statistics, GutzwillerFlag) == false)
		    {
		      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
		      return -1;
		    }
		}
	      else
		{
		  FixedTotalSzFlag = true;
		}
	    }
	  else
	    {
	      TwoDimensionTranslationFlag = true;
	    }
	}
      else
	{
	  TwoDimensionTranslationFlag = true;
	  FixedTotalSzFlag = true;
	}
    }

  MultiColumnASCIIFile KeptOrbitalFile;
  if (KeptOrbitalFile.Parse(Manager.GetString("kept-sites")) == false)
    {
      KeptOrbitalFile.DumpErrors(cout);
      return -1;
    }
  int NbrKeptOrtbitals = KeptOrbitalFile.GetNbrLines();
  int* KeptOrbitals = 0;
  if (KeptOrbitalFile.GetNbrColumns() == 1)
    {
      KeptOrbitals = KeptOrbitalFile.GetAsIntegerArray(0);
    }
  else
    {
      int* TmpKeptOrbitalX = KeptOrbitalFile.GetAsIntegerArray(0);
      int* TmpKeptOrbitalY = KeptOrbitalFile.GetAsIntegerArray(1);
      KeptOrbitals = new int[NbrKeptOrtbitals];
      for (int i = 0; i < NbrKeptOrtbitals; ++i)
	{
	  KeptOrbitals[i] = (TmpKeptOrbitalX[i] * YPeriodicity) + TmpKeptOrbitalY[i];
	}
    }
  cout << "number of kept sites = " << NbrKeptOrtbitals << endl;
  SortArrayUpOrdering(KeptOrbitals, NbrKeptOrtbitals);


  GroundStates = new ComplexVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << GroundStateFiles[i] << endl;
	  return -1;      
	}
    }

  if (Manager.GetString("output-file") != 0)
    {
      EntropyFileName = new char[strlen(Manager.GetString("output-file")) + 1];
      strcpy(Manager.GetString("output-file"), EntropyFileName);
    }
  else
    {
      EntropyFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "ent");
      if (EntropyFileName == 0)
	{
	  cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
    }
  if (Manager.GetBoolean("disable-densitymatrix") == false)
    {
      DensityMatrixFileName  = ReplaceExtensionToFileName(EntropyFileName, "ent", "full.ent");
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
      if (FixedTotalSzFlag == false)
	{
	  DensityMatrixFile << "#  N    lambda";
	}
      else
	{
	  DensityMatrixFile << "#  N    Sz    lambda";
	}
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }

  ParticleOnSphere* Space = 0;
  if (Statistics == true)
    {
      if (FixedTotalSzFlag == false)
	{
	  if (SpinlessFlag)
	    {
	      Space = new FermionOnLatticeRealSpace (NbrParticles, NbrSites);
	      cout <<"creating FermionOnLatticeRealSpace = " <<NbrParticles<<" "<< NbrSites<<endl;
	    }
	  else
	    {  
	      if (GutzwillerFlag == false)
		Space = new FermionOnLatticeWithSpinRealSpace (NbrParticles, NbrSites);
	      else
		Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
	    }
	}
      else
	{
	  if (GutzwillerFlag == false)
	    Space = new FermionOnLatticeWithSpinRealSpace (NbrParticles, TotalSz, NbrSites);
	  else
	    Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, TotalSz, NbrSites);
	}
      if (TwoDimensionTranslationFlag == true)
	{
	  ParticleOnSphere* TmpSpace = 0;
	      if (FixedTotalSzFlag == false)
	    {
	      if (GutzwillerFlag == false)
		TmpSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum, XPeriodicity, YMomentum, YPeriodicity);
	      else
		TmpSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum, XPeriodicity, YMomentum, YPeriodicity);
	    }
	  else
	    {
	      if (GutzwillerFlag == false)
		TmpSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, TotalSz, NbrSites, XMomentum, XPeriodicity, YMomentum, YPeriodicity);
	      else
		TmpSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, TotalSz, NbrSites, XMomentum, XPeriodicity, YMomentum, YPeriodicity);
	    }	  
	  for (int i = 0; i < NbrSpaces; ++i)
	    {
	      if (GroundStates[i].GetLargeVectorDimension() != TmpSpace->GetLargeHilbertSpaceDimension())
		{
		  cout << "error, " << GroundStateFiles[i] << " dimension does not match the Hilbert space dimension (" 
		       << GroundStates[i].GetLargeVectorDimension() << " vs (" << Space->GetLargeHilbertSpaceDimension() << ")" << endl;
		  return 0;
		}
	      ComplexVector TmpVector = TmpSpace->ConvertFromKxKyBasis(GroundStates[i], Space);
	      GroundStates[i] = TmpVector;
	      GroundStates[i] /= GroundStates[i].Norm();
	    }
	  delete TmpSpace;
	}
    }
  else
    {
      cout << "not available for bosons" << endl;
      return -1;
    }
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (GroundStates[i].GetLargeVectorDimension() != Space->GetLargeHilbertSpaceDimension())
	{
	  cout << "error, " << GroundStateFiles[i] << " dimension does not match the Hilbert space dimension (" 
	       << GroundStates[i].GetLargeVectorDimension() << " vs (" << Space->GetLargeHilbertSpaceDimension() << ")" << endl;
	  return 0;
	}
    }

  ofstream File;
  File.open(EntropyFileName, ios::binary | ios::out);
  File.precision(14);
  cout.precision(14);
  

  int MaxSubsystemNbrParticles = 2 * NbrKeptOrtbitals;
  if ((GutzwillerFlag == true)||( SpinlessFlag))
    {
      MaxSubsystemNbrParticles = NbrKeptOrtbitals;
    }
  if (MaxSubsystemNbrParticles > NbrParticles)
    MaxSubsystemNbrParticles = NbrParticles;
  int SubsystemNbrParticles = 0;
  int NbrRemovedOrtbitals = NbrSites - NbrKeptOrtbitals;

  double EntanglementEntropy = 0.0;
  double DensitySum = 0.0;
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  if (FixedTotalSzFlag == false)
    {
      for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
	{
	  int NbrRemovedParticles = NbrParticles - SubsystemNbrParticles;
	  if ((NbrRemovedParticles >= 0) && (NbrRemovedParticles <= (2 * NbrRemovedOrtbitals)))
	    {
	      cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << endl;      
	      if (ShowTimeFlag == true)
		{
		  gettimeofday (&(TotalStartingTime), 0);
		}
	      ComplexMatrix PartialEntanglementMatrix;
	      if( SpinlessFlag == false)
		{
		  PartialEntanglementMatrix = ((FermionOnLatticeWithSpinRealSpace*) Space)->EvaluatePartialEntanglementMatrix(SubsystemNbrParticles, NbrKeptOrtbitals, KeptOrbitals, GroundStates[0], Architecture.GetArchitecture());
		  PartialEntanglementMatrix *= Coefficients[0];
		  for (int i = 1; i < NbrSpaces; ++i)
		    {
		      ComplexMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpace*) Space)->EvaluatePartialEntanglementMatrix(SubsystemNbrParticles, NbrKeptOrtbitals, KeptOrbitals, GroundStates[i], Architecture.GetArchitecture());
		      TmpMatrix *= Coefficients[i];
		      PartialEntanglementMatrix += TmpMatrix;
		    }
		}
	      else
		{
		  PartialEntanglementMatrix = ((FermionOnLatticeRealSpace*) Space)->EvaluatePartialEntanglementMatrix(SubsystemNbrParticles, NbrKeptOrtbitals, KeptOrbitals, GroundStates[0], Architecture.GetArchitecture());
		  PartialEntanglementMatrix *= Coefficients[0];
		  for (int i = 1; i < NbrSpaces; ++i)
		    {
		      ComplexMatrix TmpMatrix = ((FermionOnLatticeRealSpace*) Space)->EvaluatePartialEntanglementMatrix(SubsystemNbrParticles, NbrKeptOrtbitals, KeptOrbitals, GroundStates[i], Architecture.GetArchitecture());
		      TmpMatrix *= Coefficients[i];
		      PartialEntanglementMatrix += TmpMatrix;
		    }
		}
	      if (ShowTimeFlag == true)
		{
		  gettimeofday (&(TotalEndingTime), 0);
		  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		  cout << "reduced density matrix evaluated in " << Dt << "s" << endl;
		}
	      PartialEntanglementMatrix.RemoveZeroColumns();
	      PartialEntanglementMatrix.RemoveZeroRows();
	      if ((PartialEntanglementMatrix.GetNbrColumn() > 0) && (PartialEntanglementMatrix.GetNbrRow() > 0))
		{
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalStartingTime), 0);
		    }
		  double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
		  int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
		  if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
		    {
		      TmpDimension = PartialEntanglementMatrix.GetNbrRow();
		    }
		  for (int i = 0; i < TmpDimension; ++i)
		    TmpValues[i] *= TmpValues[i];
		  RealDiagonalMatrix TmpDiag = RealDiagonalMatrix(TmpValues, TmpDimension);
		  TmpDiag.SortMatrixDownOrder();
		  if (DensityMatrixFileName != 0)
		    {
		      ofstream DensityMatrixFile;
		      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		      DensityMatrixFile.precision(14);
		      for (int i = 0; i < TmpDimension; ++i)
			DensityMatrixFile << SubsystemNbrParticles << " " << TmpDiag[i] << endl;
		      DensityMatrixFile.close();
		    }
		  for (int i = 0; i < TmpDimension; ++i)
		    {
		      EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
		      DensitySum +=TmpDiag[i];
		    }
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalEndingTime), 0);
		      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		      cout << "diagonalization done in " << Dt << "s" << endl;
		    }
		}
	    }
	}
    }
  else
    {
      for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
	{
	  int MaxSubsystemTotalSz = SubsystemNbrParticles;
	  if (MaxSubsystemTotalSz > NbrKeptOrtbitals)
	    {
	      MaxSubsystemTotalSz = (2 * NbrKeptOrtbitals) - SubsystemNbrParticles;
	    }
	  for (int SubsystemTotalSz = -MaxSubsystemTotalSz; SubsystemTotalSz <= MaxSubsystemTotalSz; SubsystemTotalSz += 2)
	    {
	      int NbrRemovedParticles = NbrParticles - SubsystemNbrParticles;
	      int RemovedParticleTotalSz = TotalSz - SubsystemTotalSz;
 	      if ((NbrRemovedParticles >= 0) && (NbrRemovedParticles <= (2 * NbrRemovedOrtbitals)) 
 		  && (RemovedParticleTotalSz >= -NbrRemovedParticles) && (RemovedParticleTotalSz <= NbrRemovedParticles))
		{
		  cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " SzA=" << SubsystemTotalSz << endl;      
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalStartingTime), 0);
		    }
		  ComplexMatrix PartialEntanglementMatrix;
		  PartialEntanglementMatrix = ((FermionOnLatticeWithSpinRealSpace*) Space)->EvaluatePartialEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalSz, NbrKeptOrtbitals, KeptOrbitals, GroundStates[0], Architecture.GetArchitecture());
		  PartialEntanglementMatrix *= Coefficients[0];
		  for (int i = 1; i < NbrSpaces; ++i)
		    {
		      ComplexMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpace*) Space)->EvaluatePartialEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalSz, NbrKeptOrtbitals, KeptOrbitals, GroundStates[i], Architecture.GetArchitecture());
		      TmpMatrix *= Coefficients[i];
		      PartialEntanglementMatrix += TmpMatrix;
		    }
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalEndingTime), 0);
		      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		      cout << "reduced density matrix evaluated in " << Dt << "s" << endl;
		    }
		  PartialEntanglementMatrix.RemoveZeroColumns();
		  PartialEntanglementMatrix.RemoveZeroRows();
		  if ((PartialEntanglementMatrix.GetNbrColumn() > 0) && (PartialEntanglementMatrix.GetNbrRow() > 0))
		    {
		      if (ShowTimeFlag == true)
			{
			  gettimeofday (&(TotalStartingTime), 0);
			}
		      double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
		      int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
		      if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
			{
			  TmpDimension = PartialEntanglementMatrix.GetNbrRow();
			}
		      for (int i = 0; i < TmpDimension; ++i)
			TmpValues[i] *= TmpValues[i];
		      RealDiagonalMatrix TmpDiag = RealDiagonalMatrix(TmpValues, TmpDimension);
		      TmpDiag.SortMatrixDownOrder();
		      if (DensityMatrixFileName != 0)
			{
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			  for (int i = 0; i < TmpDimension; ++i)
			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << TmpDiag[i] << endl;
			  DensityMatrixFile.close();
			}
		      for (int i = 0; i < TmpDimension; ++i)
			{
			  if (TmpDiag[i] > 0.0)
			    {
			      EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
			      DensitySum += TmpDiag[i];
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
		}
	    }
	}
    }
  File << NbrKeptOrtbitals << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << endl;
  File.close();
  return 0;
}
