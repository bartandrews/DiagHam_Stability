#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MatrixFullDiagonalizeOperation.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"
#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong.h"
#include "HilbertSpace/FermionOnSquareLatticeRealSpaceNNExclusion.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion.h"


#include "HilbertSpace/BosonOnLatticeRealSpace.h"
#include "HilbertSpace/BosonOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation.h"


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
  OptionManager Manager ("FTIRealSpaceEntanglementEntropyParticlePartition" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-na", "minimum size of the particles whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-na", "maximum size of the particles whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-sza", "minimum value of sz in a subsystem", -1000);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-sza", "maximum value of sz in a subsystem", 1000);
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*SystemGroup) += new BooleanOption  ('\n', "decoupled", "assume that the total spin is a good quantum number of the problem");
  (*SystemGroup) += new BooleanOption  ('\n', "hofstadter", "the file name uses Hofstadter model convention");
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "particles have a SU(2) spin");
  (*SystemGroup) += new SingleStringOption ('\n', "selected-blocks", "provide a column formatted ascii file that indicates which block of the reduced density matrix should be computed");
  (*SystemGroup) += new SingleStringOption ('\n', "exclusion-file", "name of the file that describes an exclusion principle (beyond the gutzwiller projection)"); 
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with partent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");

  (*OutputGroup) += new SingleIntegerOption  ('\n', "na-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Na value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "sza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Sza value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
  (*OutputGroup) += new SingleStringOption ('\n', "import-densitymatrix", "read a single block of the reduced density matrix from a file and process it");
  (*OutputGroup) += new BooleanOption ('\n', "export-densitymatrix", "write a single block of the reduced density matrix from a file and exit");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIRealSpaceEntanglementEntropyParticlePartition -h" << endl;
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
  ComplexVector* GroundStates = 0;
  char** GroundStateFiles = 0;
  int* TotalKx = 0;
  int* TotalKy = 0;
  int NbrParticles = 0;
  int NbrSites = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  int UnitCellX = 0;
  int UnitCellY = 0;
  bool Statistics = true;
  double* Coefficients = 0;
  bool ShowTimeFlag = Manager.GetBoolean("show-time");
  int TotalSpin = 0;
  int SzSymmetrySector = 0;
  bool TwoDTranslationFlag = false;
  bool SU2SpinFlag = Manager.GetBoolean("su2-spin");
  bool GutzwillerFlag = false;
  bool HofstadterFlag = Manager.GetBoolean("hofstadter");

  bool EigenstateFlag = Manager.GetBoolean("density-eigenstate");
  int FilterNa = Manager.GetInteger("na-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
  int FilterSza  = Manager.GetInteger("sza-eigenstate");
  int MinNbrSinglets = 0;

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FTIEntanglementEntropyParticlePartition -h" << endl;
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
  bool TotalSpinConservedFlag;
  if (HofstadterFlag == false)
    {     
      if (FTIHubbardModelFindSystemInfoFromVectorFileName(GroundStateFiles[0], NbrParticles, NbrSites, Statistics, GutzwillerFlag) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[0] << endl;
	  return -1;
	  
	}
      TwoDTranslationFlag = FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(GroundStateFiles[0], NbrParticles, NbrSites, TotalKx[0], TotalKy[0], 
											     NbrSiteX, NbrSiteY, Statistics, GutzwillerFlag);
    
	
      TotalSpinConservedFlag = FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(GroundStateFiles[0], NbrParticles, NbrSites, TotalSpin, Statistics, GutzwillerFlag);



      if (TwoDTranslationFlag == true)
	{ 
	  for (int i = 0; i < NbrSpaces; ++i)
	    {
	      TotalKx[i] = 0;
	      TotalKy[i] = 0;
	      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrParticles, NbrSites, TotalKx[i], TotalKy[i],
										   NbrSiteX, NbrSiteY, Statistics, GutzwillerFlag) == false)
		{
		  cout << "error while retrieving 2D translation parameters from file name " << GroundStateFiles[i] << endl;
		  return -1;
		}
	  } 
	}
    }
  else
    {
      if (FTIHofstadterdModelFindSystemInfoFromVectorFileName(GroundStateFiles[0], NbrParticles, NbrSiteX, NbrSiteY, UnitCellX, UnitCellY, Statistics, GutzwillerFlag) == false) 
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[0] << endl;
	  return -1;
	}
      NbrSites = NbrSiteX * NbrSiteY * UnitCellX * UnitCellY;
      TwoDTranslationFlag = FTIHofstadterdModelWith2DTranslationFindSystemInfoFromVectorFileName(GroundStateFiles[0], NbrParticles, TotalKx[0], TotalKy[0],
												 NbrSiteX, NbrSiteY, UnitCellX, UnitCellY, Statistics, GutzwillerFlag);
      TotalSpinConservedFlag = FTIHofstadterModelWithSzFindSystemInfoFromVectorFileName(GroundStateFiles[0], TotalSpin, SzSymmetrySector, MinNbrSinglets);
      cout << "MinNbrSinglets=" << MinNbrSinglets << endl;
      if (TotalSpinConservedFlag == true) 
	cout << "Detecting spin" << endl;
      cout << TotalSpin <<endl;
	
      if (TwoDTranslationFlag == true)
	{ 
	  cout <<"Hofstadter with translation"<<endl;
	  for (int i = 0; i < NbrSpaces; ++i)
	    {
	      TotalKx[i] = 0;
	      TotalKy[i] = 0;
	      if (FTIHofstadterdModelWith2DTranslationFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrParticles, TotalKx[i], TotalKy[i],NbrSiteX, NbrSiteY, UnitCellX, UnitCellY, Statistics, GutzwillerFlag) == false)
		{
		  cout << "error while retrieving 2D translation parameters from file name " << GroundStateFiles[i] << endl;
		  return -1;
		}
	    } 
	}
    }

  
  GroundStates = new ComplexVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << GroundStateFiles[i] << endl;
	  return -1;      
	}
    }
    
  if ((DensityMatrixFileName != 0) && (Architecture.GetArchitecture()->CanWriteOnDisk()))
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      if (TwoDTranslationFlag == false)
	{
	  if ((Manager.GetBoolean("decoupled") == false) || (SU2SpinFlag == false))
	    {
	      DensityMatrixFile << "#  N    lambda";
	    }
	  else
	    {
	      DensityMatrixFile << "#  N    Sz    lambda";
	    }
	}
      else
	{
	  if ((Manager.GetBoolean("decoupled") == false) || (SU2SpinFlag == false))
	    {
	      DensityMatrixFile << "#  N    Kx    Ky    lambda";
	    }
	  else
	    {	      
	      if (SzSymmetrySector == 0)
		{
		  DensityMatrixFile << "#  N    Kx    Ky    Sz    lambda";
		}
	      else
		{
		  DensityMatrixFile << "#  N    Kx    Ky    Sz    SzParity    lambda";
		}
	    }
      }
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }
  cout << "NbrParticles = " << NbrParticles << " NbrSites = "  << NbrSites << "TotalSz = "<<TotalSpin <<endl;

   
  
  int TotalNbrSites = NbrSiteX * NbrSiteY;
  if (TwoDTranslationFlag == false)
    TotalNbrSites = 1;
  int* NbrGroundStatePerMomentumSector = new int[TotalNbrSites];
  ComplexVector** GroundStatePerMomentumSector = new ComplexVector*[TotalNbrSites];
  double** CoefficientPerMomentumSector = new double*[TotalNbrSites];
  for (int i = 0; i < TotalNbrSites; ++i)
    {
      NbrGroundStatePerMomentumSector[i] = 0;
      GroundStatePerMomentumSector[i] = 0;
      CoefficientPerMomentumSector[i] = 0;
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << GroundStateFiles[i] << endl;
	  return -1;      
	}
      int TmpIndex;
      if(TwoDTranslationFlag == true)
	TmpIndex = (TotalKx[i] * NbrSiteY) + TotalKy[i];
      else
	TmpIndex = 0;
      NbrGroundStatePerMomentumSector[TmpIndex]++; 
    }
  for (int i = 0; i < TotalNbrSites; ++i)
    {
      if (NbrGroundStatePerMomentumSector[i] > 0)
	{
	  GroundStatePerMomentumSector[i] = new ComplexVector[NbrGroundStatePerMomentumSector[i]];
	  CoefficientPerMomentumSector[i] = new double[NbrGroundStatePerMomentumSector[i]];
	}
      NbrGroundStatePerMomentumSector[i] = 0;
    }
  for (int i = 0; i < NbrSpaces; ++i)
    {
      int TmpIndex;
      if(TwoDTranslationFlag == true )
	TmpIndex = (TotalKx[i] * NbrSiteY) + TotalKy[i];
      else
	TmpIndex = 0;
      GroundStatePerMomentumSector[TmpIndex][NbrGroundStatePerMomentumSector[TmpIndex]] = GroundStates[i];
      CoefficientPerMomentumSector[TmpIndex][NbrGroundStatePerMomentumSector[TmpIndex]] = Coefficients[i];
      NbrGroundStatePerMomentumSector[TmpIndex]++;
    }  


  int MaxNbrSpaces;
  if (TwoDTranslationFlag == false)
    MaxNbrSpaces = 1;
  else
    MaxNbrSpaces = NbrSiteX * NbrSiteY;
  ParticleOnSphere** Spaces = new ParticleOnSphere*[MaxNbrSpaces];
  for (int i = 0; i < MaxNbrSpaces; ++i)
    {
      Spaces[i] = 0;
    }
  for (int i = 0; i < NbrSpaces; ++i)
    {
      int TmpIndex;
      if (TwoDTranslationFlag == false)
	TmpIndex = 0;
      else
	TmpIndex = TotalKx[i] * NbrSiteY + TotalKy[i];
      
      if (Spaces[TmpIndex] == 0)
	{
	  if (Statistics == true)
	    {
	      if (SU2SpinFlag == false)
		{
		  if (GutzwillerFlag == false)
		    {
		      if (TwoDTranslationFlag == false)
			{
			  Spaces[TmpIndex] = new FermionOnLatticeRealSpace (NbrParticles, NbrSites);
			}
		      else
			{
			  Spaces[TmpIndex] = new FermionOnLatticeRealSpaceAnd2DTranslation (NbrParticles, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
			}
		    }
		  else
		    {
		      if (TwoDTranslationFlag == false)
			{
			  Spaces[TmpIndex] = new FermionOnSquareLatticeRealSpaceNNExclusion (NbrParticles, NbrSites, NbrSiteY);
			}
		      else
			{
			  MultiColumnASCIIFile ExclusionFile;
			  if (ExclusionFile.Parse(Manager.GetString("exclusion-file")) == false)
			    {
			      ExclusionFile.DumpErrors(cout);
			      return -1;
			    }
			  int NbrExclusionRules = ExclusionFile.GetNbrLines();
			  int* TmpOrbitalIndices = ExclusionFile.GetAsIntegerArray(3);
			  int* TmpLinearizedIndices = ExclusionFile.GetAsIntegerArray(4);
			  int TmpNbrSitesPerUnitCell = NbrSites / (NbrSiteX * NbrSiteY);
			  int* NbrExcludedSites = new int[TmpNbrSitesPerUnitCell];
			  int** ExcludedSites = new int*[TmpNbrSitesPerUnitCell];
			  for (int k = 0; k < TmpNbrSitesPerUnitCell; ++k)
			    {
			      NbrExcludedSites[k] = 0;
			    }
			  for (int k = 0; k < NbrExclusionRules; ++k)
			    {
			      NbrExcludedSites[TmpOrbitalIndices[k]]++;
			    }
			  for (int k = 0; k < TmpNbrSitesPerUnitCell; ++k)
			    {
			      ExcludedSites[k] = new int[NbrExcludedSites[k]];
			      NbrExcludedSites[k] = 0;
			    }
			  for (int k = 0; k < NbrExclusionRules; ++k)
			    {
			      ExcludedSites[TmpOrbitalIndices[k]][NbrExcludedSites[TmpOrbitalIndices[k]]] = TmpLinearizedIndices[k];
			      NbrExcludedSites[TmpOrbitalIndices[k]]++;
			    }			  
			  Spaces[TmpIndex] = new FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion (NbrParticles, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY, ExcludedSites, NbrExcludedSites);
			}		      
		    }
		}
	      else
		{
		  if (GutzwillerFlag == false)
		    {
		      if (TotalSpinConservedFlag == false)
			{
			  if (TwoDTranslationFlag == false)
			    { 
			      if (NbrSites <= 31)
				{
				  Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpace (NbrParticles, NbrSites);
				}
			      else
				{
				  Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceLong (NbrParticles, NbrSites);
				}
			    }
			  else
			    {
			      if (SzSymmetrySector != 0)
				{
				  if (MinNbrSinglets == 0)
				    {
				      if (NbrSites <= 31)
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (NbrParticles, NbrSites, (SzSymmetrySector == -1), 
															      TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
					}
				      else
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong (NbrParticles, NbrSites, (SzSymmetrySector == -1), 
																  TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
					}
				    }
				  else
				    {
				      if (NbrSites <= 31)
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets (NbrParticles, MinNbrSinglets, NbrSites, 
																	    (SzSymmetrySector == -1), 
																	    TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
					}
				      else
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong (NbrParticles, MinNbrSinglets, NbrSites, 
																		(SzSymmetrySector == -1), 
																		TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
					}
				    }
				}
			      else
				{
				  if (MinNbrSinglets == 0)
				    {
				      if (NbrSites <= 31)
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
					}
				      else
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong (NbrParticles, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
					}
				    }
				  else
				    {
				      if (NbrSites <= 31)
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets (NbrParticles, MinNbrSinglets, NbrSites, TotalKx[i], 
																  NbrSiteX, TotalKy[i], NbrSiteY);
					}
				      else
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong (NbrParticles, MinNbrSinglets, NbrSites, TotalKx[i], 
																  NbrSiteX, TotalKy[i], NbrSiteY);
					}
				    }
				}
			    }
			}
		      else
			{
			  if (TwoDTranslationFlag == false)
			    {
			      if (NbrSites <= 31)
				{
				  Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpace (NbrParticles, TotalSpin, NbrSites, 10000000);
				}
			      else
				{
				  Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceLong (NbrParticles, TotalSpin, NbrSites, 10000000);
				}
			    }
			  else
			    {
			      if (SzSymmetrySector != 0)
				{
				  if (MinNbrSinglets == 0)
				    {
				      if (NbrSites <= 31)
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (NbrParticles, TotalSpin, NbrSites, (SzSymmetrySector == -1),
															      TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY, 10000000);
					}
				      else
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong (NbrParticles, TotalSpin, NbrSites, 
																  (SzSymmetrySector == -1),
																  TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY, 10000000);
					}
				    }
				  else
				    {
				      if (NbrSites <= 31)
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets (NbrParticles, MinNbrSinglets, TotalSpin, 
																	    NbrSites, (SzSymmetrySector == -1),
																	    TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
					}
				      else
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong (NbrParticles, MinNbrSinglets, TotalSpin, 
																	    NbrSites, (SzSymmetrySector == -1),
																	    TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
					}
				    }
				}
			      else
				{
				  if (MinNbrSinglets == 0)
				    {
				      if (NbrSites <= 31)
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, TotalSpin, NbrSites, 
														    TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY, 10000000);
					}
				      else
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong (NbrParticles, TotalSpin, NbrSites, 
															TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY, 10000000);
					}
				    }
				  else
				    {
				      if (NbrSites <= 31)
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets (NbrParticles, MinNbrSinglets, TotalSpin, NbrSites, 
																  TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY, 10000000);
					}
				      else
					{
					  Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong (NbrParticles, MinNbrSinglets, TotalSpin, NbrSites, 
																      TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY, 10000000);
					}
				    }
				}
			    }
			}
		    }
		  else
		    {
		      if (TotalSpinConservedFlag == false)
			{
			  if (TwoDTranslationFlag == false)
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
			  else
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
			}
		      else
			{
			  if (TwoDTranslationFlag == false)
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, TotalSpin, NbrSites, 10000000);
			  else
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, TotalSpin, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY, 10000000);
			}
		      
		    }
		}
	    }
	  else
	    {
	      if (SU2SpinFlag == false)
		{
		  if (GutzwillerFlag == false)
		    {
		      if (TwoDTranslationFlag == false)
			Spaces[TmpIndex] = new BosonOnLatticeRealSpace (NbrParticles, NbrSites);
		      else
			Spaces[TmpIndex] = new BosonOnLatticeRealSpaceAnd2DTranslation (NbrParticles, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
		    }
		  else
		    {
		      if (TwoDTranslationFlag == false)
			Spaces[TmpIndex] = new BosonOnLatticeGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
		      else
			Spaces[TmpIndex] = new BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
		    }
		}
	      else
		{
		  cout << "Error: Bosonic statistics not implemented" << endl;
		  return -1;		  
		}
	    }
	}
      
      if (Spaces[TmpIndex]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	{
	  cout << Spaces[TmpIndex]->GetLargeHilbertSpaceDimension() << " " << GroundStates[i].GetLargeVectorDimension() << endl;
	  cout << "dimension mismatch between Hilbert space and ground state" << endl;
	  return 0;
	}
    }

  char* OutputFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "partent");
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy(OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "partent");
      if (OutputFileName == 0)
	{
	  cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
    }

  if (Architecture.GetArchitecture()->CanWriteOnDisk())  
    {
      ofstream File;
      File.open(OutputFileName, ios::binary | ios::out);
      File.precision(14);
      File.close();
    }
  cout.precision(14);
  
  int MaxSubsystemNbrParticles = (NbrParticles >> 1) + (NbrParticles & 1);
  if (Manager.GetInteger("max-na") > 0)
    MaxSubsystemNbrParticles = Manager.GetInteger("max-na");
  int MinSubsystemNbrParticles = Manager.GetInteger("min-na");


  int TotalNbrReducedDensityMatrixBlocks = 0;
  int* SubsystemNbrParticleSectors = 0;
  int* SubsystemTotalSzSectors = 0;
  int* SubsystemSzSymmetrySectors = 0;
  int* SubsystemKxSectors = 0; 
  int* SubsystemKySectors = 0; 
  if (Manager.GetString("selected-blocks") == 0)
    {
      for (int SubsystemNbrParticles = MinSubsystemNbrParticles; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
	{
	  int ComplementaryNbrParticles = NbrParticles - SubsystemNbrParticles;
	  int MinSz = -SubsystemNbrParticles;
	  if (Manager.GetInteger("min-sza") > MinSz)
	    MinSz = Manager.GetInteger("min-sza");
	  if ((TotalSpin - ComplementaryNbrParticles) > MinSz)
	    MinSz = (TotalSpin - ComplementaryNbrParticles);
	  int MaxSz = SubsystemNbrParticles;
	  if (Manager.GetInteger("max-sza") < MaxSz)
	    MaxSz = Manager.GetInteger("max-sza");
	  if ((TotalSpin + ComplementaryNbrParticles) < MaxSz)
	    MaxSz = (TotalSpin + ComplementaryNbrParticles);
	  if ((Manager.GetBoolean("decoupled") == false) || (Manager.GetBoolean("su2-spin")) == false)
	    {
	      MinSz = 0;
	      MaxSz = 0;
	    }
	  int SubsystemTotalKxMin = 0;
	  int SubsystemTotalKyMin = 0;
	  int SubsystemTotalKxMax = 1;
	  int SubsystemTotalKyMax = 1;
	  if(TwoDTranslationFlag == true)
	    {
	      SubsystemTotalKxMax = NbrSiteX;
	      SubsystemTotalKyMax = NbrSiteY;
	    }
	  
	  for (int SubsystemTotalSz = MinSz; SubsystemTotalSz <= MaxSz; SubsystemTotalSz += 2)
	    {
	      int SubsystemSzSymmetrySector = 0;
	      int MaxSzSymmetrySector = 0;
	      if ((SubsystemTotalSz == 0) & (SzSymmetrySector != 0))
		{
		  SubsystemSzSymmetrySector = -1;
		  MaxSzSymmetrySector = 1;
		}
	      for (; SubsystemSzSymmetrySector <= MaxSzSymmetrySector; SubsystemSzSymmetrySector += 2)
		{
		  for (int SubsystemTotalKx = SubsystemTotalKxMin; SubsystemTotalKx < SubsystemTotalKxMax; ++SubsystemTotalKx)
		    {
		      for (int SubsystemTotalKy = SubsystemTotalKyMin; SubsystemTotalKy < SubsystemTotalKyMax; ++SubsystemTotalKy)
			{
			  ++TotalNbrReducedDensityMatrixBlocks;
			}
		    }
		}
	    }
	}
      SubsystemNbrParticleSectors = new int [TotalNbrReducedDensityMatrixBlocks];
      SubsystemTotalSzSectors = new int [TotalNbrReducedDensityMatrixBlocks];
      SubsystemSzSymmetrySectors = new int [TotalNbrReducedDensityMatrixBlocks];
      SubsystemKxSectors = new int [TotalNbrReducedDensityMatrixBlocks]; 
      SubsystemKySectors = new int [TotalNbrReducedDensityMatrixBlocks]; 
      TotalNbrReducedDensityMatrixBlocks = 0;  
      for (int SubsystemNbrParticles = MinSubsystemNbrParticles; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
	{
	  int ComplementaryNbrParticles = NbrParticles - SubsystemNbrParticles;
	  int MinSz = -SubsystemNbrParticles;
	  if (Manager.GetInteger("min-sza") > MinSz)
	    MinSz = Manager.GetInteger("min-sza");
	  if ((TotalSpin - ComplementaryNbrParticles) > MinSz)
	    MinSz = (TotalSpin - ComplementaryNbrParticles);
	  int MaxSz = SubsystemNbrParticles;
	  if (Manager.GetInteger("max-sza") < MaxSz)
	    MaxSz = Manager.GetInteger("max-sza");
	  if ((TotalSpin + ComplementaryNbrParticles) < MaxSz)
	    MaxSz = (TotalSpin + ComplementaryNbrParticles);
	  if ((Manager.GetBoolean("decoupled") == false) || (Manager.GetBoolean("su2-spin")) == false)
	    {
	      MinSz = 0;
	      MaxSz = 0;
	    }
	  int SubsystemTotalKxMin = 0;
	  int SubsystemTotalKyMin = 0;
	  int SubsystemTotalKxMax = 1;
	  int SubsystemTotalKyMax = 1;
	  if(TwoDTranslationFlag == true)
	    {
	      SubsystemTotalKxMax = NbrSiteX;
	      SubsystemTotalKyMax = NbrSiteY;
	    }
	  
	  for (int SubsystemTotalSz = MinSz; SubsystemTotalSz <= MaxSz; SubsystemTotalSz += 2)
	    {
	      int SubsystemSzSymmetrySector = 0;
	      int MaxSzSymmetrySector = 0;
	      if ((SubsystemTotalSz == 0) & (SzSymmetrySector != 0))
		{
		  SubsystemSzSymmetrySector = -1;
		  MaxSzSymmetrySector = 1;
		}
	      for (; SubsystemSzSymmetrySector <= MaxSzSymmetrySector; SubsystemSzSymmetrySector += 2)
		{
		  for (int SubsystemTotalKx = SubsystemTotalKxMin; SubsystemTotalKx < SubsystemTotalKxMax; ++SubsystemTotalKx)
		    {
		      for (int SubsystemTotalKy = SubsystemTotalKyMin; SubsystemTotalKy < SubsystemTotalKyMax; ++SubsystemTotalKy)
			{
			  SubsystemNbrParticleSectors[TotalNbrReducedDensityMatrixBlocks] = SubsystemNbrParticles;
			  SubsystemTotalSzSectors[TotalNbrReducedDensityMatrixBlocks] = SubsystemTotalSz;
			  SubsystemSzSymmetrySectors[TotalNbrReducedDensityMatrixBlocks] = SubsystemSzSymmetrySector;
			  SubsystemKxSectors[TotalNbrReducedDensityMatrixBlocks] = SubsystemTotalKx;
			  SubsystemKySectors[TotalNbrReducedDensityMatrixBlocks] = SubsystemTotalKy;
			  ++TotalNbrReducedDensityMatrixBlocks;
			}
		    }
		}
	    }
	}
    }
  else
    {
      MultiColumnASCIIFile BlockFile;
      if (BlockFile.Parse(Manager.GetString("selected-blocks")) == false)
	{
	  BlockFile.DumpErrors(cout);
	  return -1;
	}
      TotalNbrReducedDensityMatrixBlocks = BlockFile.GetNbrLines();
      SubsystemNbrParticleSectors = BlockFile.GetAsIntegerArray(0);
      MinSubsystemNbrParticles = 100000;
      MaxSubsystemNbrParticles = 0;
      for (int i = 0; i < TotalNbrReducedDensityMatrixBlocks; ++i)
	{
	  if (SubsystemNbrParticleSectors[i] < MinSubsystemNbrParticles)
	    MinSubsystemNbrParticles = SubsystemNbrParticleSectors[i];
	  if (SubsystemNbrParticleSectors[i] > MaxSubsystemNbrParticles)
	    MaxSubsystemNbrParticles = SubsystemNbrParticleSectors[i];
	}
      if(TwoDTranslationFlag == true)
	{
	  SubsystemKxSectors = BlockFile.GetAsIntegerArray(1);
	  SubsystemKySectors = BlockFile.GetAsIntegerArray(2);
	  if ((Manager.GetBoolean("decoupled") == false) || (Manager.GetBoolean("su2-spin")) == false)
	    {
	      for (int i = 0; i < TotalNbrReducedDensityMatrixBlocks; ++i)
		{
		  SubsystemTotalSzSectors[i] = 0;
		  SubsystemSzSymmetrySectors[i] = 0;
		}
	    }
	  else
	    {
	      SubsystemTotalSzSectors = BlockFile.GetAsIntegerArray(3);
	      if (BlockFile.GetNbrColumns() >= 5)
		{
		  SubsystemSzSymmetrySectors = BlockFile.GetAsIntegerArray(4);
		}
	      else
		{
		  for (int i = 0; i < TotalNbrReducedDensityMatrixBlocks; ++i)
		    {
		      SubsystemSzSymmetrySectors[i] = 0;
		    }
		}
	    }
	}
      else
	{
	  for (int i = 0; i < TotalNbrReducedDensityMatrixBlocks; ++i)
	    {
	      SubsystemKxSectors[i] = 0;
	      SubsystemKySectors[i] = 0;
	    }
	  if ((Manager.GetBoolean("decoupled") == false) || (Manager.GetBoolean("su2-spin")) == false)
	    {
	      for (int i = 0; i < TotalNbrReducedDensityMatrixBlocks; ++i)
		{
		  SubsystemTotalSzSectors[i] = 0;
		  SubsystemSzSymmetrySectors[i] = 0;
		}
	    }
	  else
	    {
	      SubsystemTotalSzSectors = BlockFile.GetAsIntegerArray(1);
	      if (BlockFile.GetNbrColumns() >= 3)
		{
		  SubsystemSzSymmetrySectors = BlockFile.GetAsIntegerArray(2);
		}
	      else
		{
		  for (int i = 0; i < TotalNbrReducedDensityMatrixBlocks; ++i)
		    {
		      SubsystemSzSymmetrySectors[i] = 0;
		    }
		}
	    }
	}
    }

  double* EntanglementEntropies = new double[MaxSubsystemNbrParticles - MinSubsystemNbrParticles + 1];
  double* DensitySums = new double[MaxSubsystemNbrParticles - MinSubsystemNbrParticles + 1];
  for (int SubsystemNbrParticles = MinSubsystemNbrParticles; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
    {
      EntanglementEntropies[SubsystemNbrParticles - MinSubsystemNbrParticles] = 0.0;
      DensitySums[SubsystemNbrParticles - MinSubsystemNbrParticles] = 0.0;
    }

  for (int BlockIndex = 0; BlockIndex < TotalNbrReducedDensityMatrixBlocks; ++BlockIndex)
    {
      int SubsystemNbrParticles = SubsystemNbrParticleSectors[BlockIndex];
      int SubsystemTotalSz = SubsystemTotalSzSectors[BlockIndex];
      int SubsystemSzSymmetrySector = SubsystemSzSymmetrySectors[BlockIndex];
      int SubsystemTotalKx = SubsystemKxSectors[BlockIndex];
      int SubsystemTotalKy = SubsystemKySectors[BlockIndex];
      if ((SU2SpinFlag == false) || (Manager.GetBoolean("decoupled") == false))
	{      
	  if(TwoDTranslationFlag == false)
	    cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << endl;
	  else
	    cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Kx=" << SubsystemTotalKx << " Ky=" << SubsystemTotalKy << endl;
	}
      else
	{   
	  cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles;   
	  if(TwoDTranslationFlag == true)
	    cout << " subsystem total Kx=" << SubsystemTotalKx << " Ky=" << SubsystemTotalKy ;
	  cout << " Sz = " << SubsystemTotalSz;
	  if ((SubsystemSzSymmetrySector != 0) && (SzSymmetrySector != 0))
	    cout << " Sz parity sector=" << SubsystemSzSymmetrySector;
	  cout << endl;
	}
      timeval TotalStartingTime;
      timeval TotalEndingTime;
      if (ShowTimeFlag == true)
	{
	  gettimeofday (&(TotalStartingTime), 0);
	}
      int TmpIndex = 0;
      while (NbrGroundStatePerMomentumSector[TmpIndex] == 0)
	++TmpIndex;
      HermitianMatrix PartialDensityMatrix;
      if (Manager.GetString("import-densitymatrix") == 0)
	{
	  if (Statistics == true)
	    {
	      if (SU2SpinFlag == false)
		{      
		  if (TwoDTranslationFlag == false)
		    {
		      PartialDensityMatrix = ((FermionOnLatticeRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
		      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
		      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
			{
			  HermitianMatrix TmpMatrix = ((FermionOnLatticeRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
			  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
			  PartialDensityMatrix += TmpMatrix;
			}
		    }
		  else
		    {
		      PartialDensityMatrix = ((FermionOnLatticeRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
		      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
		      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
			{
			  HermitianMatrix TmpMatrix = ((FermionOnLatticeRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
			  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
			  PartialDensityMatrix += TmpMatrix;
			}
		    }
		}
	      else
		{
		  if (GutzwillerFlag == false)
		    {
		      if (Manager.GetBoolean("decoupled") == false)
			{
			  if (TwoDTranslationFlag == false)
			    {
			      PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
			      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
			      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				{
				  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
				  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
				  PartialDensityMatrix += TmpMatrix;
				}
			    }
			  else
			    {
			      if (NbrSites <= 31)
				{
				  if (MinNbrSinglets == 0)
				    {
				      PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				  else
				    {
				      PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				}
			      else
				{
				  if (MinNbrSinglets == 0)
				    {
				      PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				  else
				    {
				      PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				}
			    }
			}
		      else
			{
			  if (TwoDTranslationFlag == false)
			    {
			      if (NbrSites <= 31)
				{
				  PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				  for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				    {
				      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
				      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
				      PartialDensityMatrix += TmpMatrix;
				    }
				}
			      else
				{
				  PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpaceLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				  for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				    {
				      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
				      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
				      PartialDensityMatrix += TmpMatrix;
				    }
				}
			    }
			  else
			    {
			      if (SubsystemSzSymmetrySector == 0)
				{
				  if (NbrSites <= 31)
				    {
				      // 					      PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles,  SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex] , Architecture.GetArchitecture());
				      if (MinNbrSinglets == 0)
					{
					  PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles,  SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
					  for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					    {
					      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					      PartialDensityMatrix += TmpMatrix;
					    }
					}
				      else
					{
					  
					  PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles,  SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
					  for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					    {
					      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					      PartialDensityMatrix += TmpMatrix;
					    }
					  
					}
				    }
				  else
				    {
				      if (MinNbrSinglets == 0)
					{
					  PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles,  SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
					  for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					    {
					      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					      PartialDensityMatrix += TmpMatrix;
					    }
					}
				      else
					{
					  PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles,  SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
					  for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					    {
					      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					      PartialDensityMatrix += TmpMatrix;
					    }
					}
				    }
				}
			      else
				{
				  if (NbrSites <= 31)
				    {
				      if (MinNbrSinglets == 0)
					{
					  PartialDensityMatrix = ((FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles,  SubsystemTotalSz, SubsystemSzSymmetrySector, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
					  for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					    {
					      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemSzSymmetrySector, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					      PartialDensityMatrix += TmpMatrix;
					    }
					}
				      else
					{
					  PartialDensityMatrix = ((FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets *) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles,  SubsystemTotalSz, SubsystemSzSymmetrySector, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
					  for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					    {
					      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets *) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemSzSymmetrySector, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					      PartialDensityMatrix += TmpMatrix;
					  
					    }
					}
				    }
				  else
				    {
				      if (MinNbrSinglets == 0)
					{
					  PartialDensityMatrix = ((FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles,  SubsystemTotalSz, SubsystemSzSymmetrySector, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
					  for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					    {
					      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemSzSymmetrySector, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					      PartialDensityMatrix += TmpMatrix;
					    } 
					}
				      else
					{
					  PartialDensityMatrix = ((FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles,  SubsystemTotalSz, SubsystemSzSymmetrySector, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
					  for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					    {
					      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemSzSymmetrySector, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					      PartialDensityMatrix += TmpMatrix;
					    } 
					  
					}
				    }
				}					  
			    }
			}
		    }
		  else
		    {
		      if (Manager.GetBoolean("decoupled") == false)
			{
			  if (TwoDTranslationFlag == false)
			    {
			      PartialDensityMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
			      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
			      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				{
				  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
				  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
				  PartialDensityMatrix += TmpMatrix;
				}
			    }
			  else
			    {
			      PartialDensityMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
			      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
			      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				{
				  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
				  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
				  PartialDensityMatrix += TmpMatrix;
				}
			    }
			}
		      else
			{
			  if (TwoDTranslationFlag == false)
			    {
			      PartialDensityMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
			      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
			      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				{
				  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
				  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
				  PartialDensityMatrix += TmpMatrix;
				}
			    }
			  else
			    {
			      PartialDensityMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
			      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
			      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				{
				  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
				  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
				  PartialDensityMatrix += TmpMatrix;
				}
			    }
			}
		    }
		}
	    }
	  else
	    {
	      if (SU2SpinFlag == false)
		{      
		  if (GutzwillerFlag == false)
		    {
		      if (TwoDTranslationFlag == false)
			{
			  PartialDensityMatrix = ((BosonOnLatticeRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
			  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
			  for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
			    {
			      HermitianMatrix TmpMatrix = ((BosonOnLatticeRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
			      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
			      PartialDensityMatrix += TmpMatrix;
			    }
			}
		      else
			{
			  PartialDensityMatrix = ((BosonOnLatticeRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
			  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
			  for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
			    {
			      HermitianMatrix TmpMatrix = ((BosonOnLatticeRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
			      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
			      PartialDensityMatrix += TmpMatrix;
			    }
			}
		    }
		  else
		    {
		      if (TwoDTranslationFlag == false)
			{
			  PartialDensityMatrix = ((BosonOnLatticeGutzwillerProjectionRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
			  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
			  for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
			    {
			      HermitianMatrix TmpMatrix = ((BosonOnLatticeGutzwillerProjectionRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
			      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
			      PartialDensityMatrix += TmpMatrix;
			    }
			}
		      else
			{
			  PartialDensityMatrix = ((BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
			  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
			  for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
			    {
			      HermitianMatrix TmpMatrix = ((BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
			      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
			      PartialDensityMatrix += TmpMatrix;
			    }
			}
		    }
		}
	      else
		{
		  cout << "Error: Bosonic statistics not implemented" << endl;
		  return -1;
		}
	    }
	  ++TmpIndex;
	  
	  for (; TmpIndex < TotalNbrSites; ++TmpIndex)
	    {
	      if (NbrGroundStatePerMomentumSector[TmpIndex] != 0)
		{
		  if (Statistics == true)
		    {
		      if (SU2SpinFlag == false)
			{      
			  if (TwoDTranslationFlag == false)
			    {
			      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				{
				  HermitianMatrix TmpMatrix = ((FermionOnLatticeRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
				  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
				  PartialDensityMatrix += TmpMatrix;
				}
			    }
			  else
			    {				      
			      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				{
				  HermitianMatrix TmpMatrix = ((FermionOnLatticeRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
			      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
			      PartialDensityMatrix += TmpMatrix;
				}
			    }
			}
		      else
			{
			  if (GutzwillerFlag == false)
			    {
			      if (Manager.GetBoolean("decoupled") == false)
				{
				  if (TwoDTranslationFlag == false)
				    {
				      if (NbrSites <= 31)
					{
					  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					    {
					      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					      PartialDensityMatrix += TmpMatrix;
					    }
					}
				      else
					{
					  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					    {
					      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					      PartialDensityMatrix += TmpMatrix;
					    }
					}
				    }
				  else
				    {
				      if (NbrSites <= 31)
					{
					  if (MinNbrSinglets == 0)
					    {
					      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						{
						  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles,  SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						  PartialDensityMatrix += TmpMatrix;
						}
					    }
					}
				      else
					{
					  if (MinNbrSinglets == 0)
					    {
					      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						{
						  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles,  SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						  PartialDensityMatrix += TmpMatrix;
						}
					    }
					  else
					    {
					      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						{
						  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles,  SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						  PartialDensityMatrix += TmpMatrix;
						}
					      
					    }
					  
					}
				    }
				}
			      else
				{
				  if (TwoDTranslationFlag == false)
				    {
				      if (NbrSites <= 31)
					{
					  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					    {
					      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					      PartialDensityMatrix += TmpMatrix;
					    }
					}
				      else
					{
					  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					    {
					      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					      PartialDensityMatrix += TmpMatrix;
					    } 
					}
				    }
				  else
				    {
				      if (SubsystemSzSymmetrySector == 0)
					{
					  if(NbrSites <= 31)
					    {
					      if (MinNbrSinglets == 0)
						{
						  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						    {
						      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						      PartialDensityMatrix += TmpMatrix;
						    }
						}
					      else
						{
						  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						    {
						      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						      PartialDensityMatrix += TmpMatrix;
						    }
						  
						}
					    }
					  else
					    {
					      if (MinNbrSinglets == 0)
						{
						  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						    {
						      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong *) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						      PartialDensityMatrix += TmpMatrix;
						    }
						}
					      else
						{
						  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						    {
						      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong *) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						      PartialDensityMatrix += TmpMatrix;
						    }
						}
					    }
					}
				      else
					{
					  
					  
					  if(NbrSites <= 31)
					    {
					      if (MinNbrSinglets == 0)
						{
						  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						    {
						      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemSzSymmetrySector, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						      PartialDensityMatrix += TmpMatrix;
						    }
						}
					      else
						{
						  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						    {
						      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets *) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemSzSymmetrySector, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						      PartialDensityMatrix += TmpMatrix;
						    }
						}
					    }
					  else
					    {
					      
					      if (MinNbrSinglets == 0)
						{
						  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						    {
						      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemSzSymmetrySector, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						      PartialDensityMatrix += TmpMatrix;
						    }
						}
					      else
						{
						  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						    {
						      HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemSzSymmetrySector, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						      PartialDensityMatrix += TmpMatrix;
						    }
						  
						}
					    }
					}
				    }
				}
			    }
			  else
			    {
			      if (Manager.GetBoolean("decoupled") == false)
				{
				  if (TwoDTranslationFlag == false)
				    {
				      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				  else
				    {
				      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles,  SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				}
			      else
				{
				  if (TwoDTranslationFlag == false)
				    {
				      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				  else
				    {
				      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				}
			    }
			}
		    }
		  else
		    {
		      if (SU2SpinFlag == false)
			{      
			  if (GutzwillerFlag == false)
			    {
			      if (TwoDTranslationFlag == false)
				{
				  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				    {
				      HermitianMatrix TmpMatrix = ((BosonOnLatticeRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
				      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
				      PartialDensityMatrix += TmpMatrix;
				    }
				}
			      else
				{
				  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				    {
				      HermitianMatrix TmpMatrix = ((BosonOnLatticeRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
				      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
				      PartialDensityMatrix += TmpMatrix;
				    }
				}
			    }
			  else
			    {
			      if (TwoDTranslationFlag == false)
				{
				  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				    {
				      HermitianMatrix TmpMatrix = ((BosonOnLatticeGutzwillerProjectionRealSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
				      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
				      PartialDensityMatrix += TmpMatrix;
				    }
				}
			      else
				{
				  for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				    {
				      HermitianMatrix TmpMatrix = ((BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
				      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
				      PartialDensityMatrix += TmpMatrix;
				    }
				}
			    }
			}
		      else
			{
			  cout << "Error: Bosonic statistics not implemented" << endl;
			  return -1;
			}
		    }			  
		}
	    }
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalEndingTime), 0);
	      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "reduced density matrix evaluated in " << Dt << "s" << endl;
	      
	    }
	  if (Manager.GetBoolean("export-densitymatrix") == true)
	    {
	      if (Architecture.GetArchitecture()->CanWriteOnDisk())
		{
		  char* TmpBlockFileName = 0;
		  if (DensityMatrixFileName == 0)
		    {
		      TmpBlockFileName = new char[256];
		      sprintf (TmpBlockFileName, "densitymatrix_na_%d_kxa_%d_kya_%d_sza_%d_szsyma_%d.mat", SubsystemNbrParticles, SubsystemTotalKx, 
			       SubsystemTotalKy, SubsystemTotalSz, SubsystemSzSymmetrySector);		  
		    }
		  else
		    {
		      TmpBlockFileName = new char[strlen(DensityMatrixFileName) + 256];
		      sprintf (TmpBlockFileName, "%s_na_%d_kxa_%d_kya_%d_sza_%d_szsyma_%d.mat", DensityMatrixFileName, SubsystemNbrParticles, SubsystemTotalKx, 
			       SubsystemTotalKy, SubsystemTotalSz, SubsystemSzSymmetrySector);
		    }
		  if (PartialDensityMatrix.WriteMatrix(TmpBlockFileName) == false)
		    {
		      cout << "error, can't write the reduced density matrix block " << TmpBlockFileName << endl;
		      return -1;
		    }
		  delete[] TmpBlockFileName;
		}
	      return 0;
	    }
	}
      else
	{
	  if (Manager.GetBoolean("use-scalapack") == true)
	    {
	      if (Architecture.GetArchitecture()->CanWriteOnDisk())
		{
		  if (PartialDensityMatrix.ReadMatrix(Manager.GetString("import-densitymatrix")) == false)
		    {
		      cout << "error, can't read the reduced density matrix block " << Manager.GetString("import-densitymatrix") << endl;
		      return -1;
		    }
		}
	      else
		{
		  PartialDensityMatrix = HermitianMatrix(2, true);
		  PartialDensityMatrix.SetMatrixElement(0, 0, 0.5);
		  PartialDensityMatrix.SetMatrixElement(1, 1, 0.5);
		}
	    }
	  else
	    {
	      if (Architecture.GetArchitecture()->CanWriteOnDisk())
		{
		  if (PartialDensityMatrix.ReadMatrix(Manager.GetString("import-densitymatrix")) == false)
		    {
		      cout << "error, can't read the reduced density matrix block " << Manager.GetString("import-densitymatrix") << endl;
		      return -1;
		    }
		}	      
	    }
	}
      if (PartialDensityMatrix.GetNbrRow() > 1)
	{
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalStartingTime), 0);
	      
	    }
	  RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
#ifdef __SCALAPACK__
	  if (Manager.GetBoolean("use-scalapack") == true)
	    {
	      if ((Architecture.GetArchitecture()->GetArchitectureID() & AbstractArchitecture::WithCommunicator) == 0)
		{
		  cout << "error : SCALAPACK requires a MPI enable architecture" << endl;
		  return 1;
		}	  
	      MatrixFullDiagonalizeOperation TmpOperation(&PartialDensityMatrix, true);
	      TmpOperation.ApplyOperation(Architecture.GetArchitecture());
	      TmpDiag = TmpOperation.GetDiagonalizedMatrix();
	    }
	  else
	    {
#endif		  
#ifdef __LAPACK__
	      if (LapackFlag == true)
		{
		  if (((EigenstateFlag == true) && (FilterNa ==  SubsystemNbrParticles) )&& (((SU2SpinFlag == false) || (Manager.GetBoolean("decoupled") == false)) || (FilterSza ==  SubsystemTotalSz)  ))
		    {
		      ComplexMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(), PartialDensityMatrix.GetNbrRow(), true);
		      TmpEigenstates.SetToIdentity();
		      PartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
		      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
		      char* TmpEigenstateName;
		      char* TmpSuffix = new char [512];
		      int MaxNbrEigenstates = NbrEigenstates;
		      if (NbrEigenstates == 0)
			MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
		      for (int i = 0; i < MaxNbrEigenstates; ++i)
			{
			  if (TmpDiag[i] > 1e-14)
			    {
			      if ((SU2SpinFlag == false) || (Manager.GetBoolean("decoupled") == false))
				{				  
				  sprintf(TmpSuffix, "partent_na_%d_kxa_%d_kya_%d.%d.vec", SubsystemNbrParticles, SubsystemTotalKx,SubsystemTotalKy, i);
				}
			      else
				{
				  if (SubsystemSzSymmetrySector != 0 )
				    {
				      sprintf(TmpSuffix, "partent_na_%d_sza_%d_szsyma_%d_kxa_%d_kya_%d.%d.vec", SubsystemNbrParticles, SubsystemTotalSz, SubsystemSzSymmetrySector ,SubsystemTotalKx,SubsystemTotalKy, i);
				    }
				  else
				    {
				      sprintf(TmpSuffix, "partent_na_%d_sza_%d_kxa_%d_kya_%d.%d.vec", SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx,SubsystemTotalKy, i);
				    }
				}
			      TmpEigenstateName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", TmpSuffix);
			      TmpEigenstates[i].WriteVector(TmpEigenstateName);
			    }
			}
		      delete[] TmpEigenstateName;
		      delete[] TmpSuffix;
		    }
		  else
		    PartialDensityMatrix.LapackDiagonalize(TmpDiag);
		}
	      else
		{
		  PartialDensityMatrix.Diagonalize(TmpDiag);
		}
#else
	      PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
#ifdef __SCALAPACK__
	    }
#endif		  
	  TmpDiag.SortMatrixDownOrder();
	  
	  if ((DensityMatrixFileName != 0) &&(Architecture.GetArchitecture()->CanWriteOnDisk()))
	    { 
	      ofstream DensityMatrixFile;
	      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 		      
	      DensityMatrixFile.precision(14);
	      double Trace = 0.0;
	      if (TwoDTranslationFlag == false)
		{
		  if ((SU2SpinFlag == false) || (Manager.GetBoolean("decoupled") == false))
		    {
		      
		      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			{
			  DensityMatrixFile << SubsystemNbrParticles << " " << TmpDiag[i] << endl;
			  Trace += TmpDiag[i];
			}
		      cout << "Trace = " << Trace << endl;
		    }
		  else
		    {
		      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			{
			  DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << TmpDiag[i] << endl;
			  Trace += TmpDiag[i];
			}
		      cout << "Trace = " << Trace << endl;
		    }
		}
	      else
		{
		  if ((SU2SpinFlag == false) || (Manager.GetBoolean("decoupled") == false))
		    {				  
		      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << TmpDiag[i] << endl;
		    }
		  else
		    {
		      if (SzSymmetrySector == 0)
			{
			  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << SubsystemTotalSz 
					      << " " << TmpDiag[i] << endl;
			}
		      else
			{
			  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << SubsystemTotalSz 
					      << " " <<  SubsystemSzSymmetrySector  << " " << TmpDiag[i] << endl;
			}
		    }
		}
	      DensityMatrixFile.close();
	    }
	  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
	    {
	      if (TmpDiag[i] > 1e-14)
		{
		  EntanglementEntropies[SubsystemNbrParticles - MinSubsystemNbrParticles] += TmpDiag[i] * log(TmpDiag[i]);
		  DensitySums[SubsystemNbrParticles - MinSubsystemNbrParticles] +=TmpDiag[i];
		}
	    }
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalEndingTime), 0);
	      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "diagonalization done in " << Dt << "s" << endl;
	    }
	}
      else
	{
	  if (PartialDensityMatrix.GetNbrRow() == 1)
	    {
	      double TmpValue = PartialDensityMatrix(0,0);
	      if ((DensityMatrixFileName != 0) &&(Architecture.GetArchitecture()->CanWriteOnDisk()))
		{
		  ofstream DensityMatrixFile;
		  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 		DensityMatrixFile.precision(14);
		  if (TwoDTranslationFlag == false)
		    {
		      if ((SU2SpinFlag == false) || (Manager.GetBoolean("decoupled") == false))
			{
			  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemNbrParticles << " " << TmpValue << endl;
					}
		      else
			{
			  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << TmpValue << endl;
			}
		    }
		  else
		    {
		      if ((SU2SpinFlag == false) || (Manager.GetBoolean("decoupled") == false))
			{
			  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << TmpValue << endl;
			}
		      else
			{
			  if (SzSymmetrySector == 0)
			    {
			      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << SubsystemTotalSz 
						  << " " <<  TmpValue << endl;
			    }
			  else
			    {
			      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << SubsystemTotalSz 
						  << " " <<  SubsystemSzSymmetrySector << " " << TmpValue << endl;
			    }
			  
			}
		    }
		  DensityMatrixFile.close();
		}
	      if (TmpValue > 1e-14)
		{
		  EntanglementEntropies[SubsystemNbrParticles - MinSubsystemNbrParticles] += TmpValue * log(TmpValue);
		  DensitySums[SubsystemNbrParticles - MinSubsystemNbrParticles] += TmpValue;
		}
	    }
	}
    }
  if (Architecture.GetArchitecture()->CanWriteOnDisk())  
    {
      ofstream File;
      File.open(OutputFileName, ios::binary | ios::out | ios::app);
      File.precision(14);
      for (int SubsystemNbrParticles = MinSubsystemNbrParticles; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
	{
	  File << SubsystemNbrParticles << " " << (-EntanglementEntropies[SubsystemNbrParticles - MinSubsystemNbrParticles]) 
	       << " " << DensitySums[SubsystemNbrParticles - MinSubsystemNbrParticles] << " " << (1.0 - DensitySums[SubsystemNbrParticles - MinSubsystemNbrParticles]) << endl;
	}
      File.close();
    }
  
  return 0;
}
