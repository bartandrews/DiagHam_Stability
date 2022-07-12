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

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/BosonOnCubicLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnCubicLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWannierSpace.h"


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
  OptionManager Manager ("FTIEntanglementEntropyParticlePartition" , "0.01");
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
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-subbands", "number of subbands", 1);
  (*SystemGroup) += new MultipleIntegerOption  ('\n', "enhance", "symmetry-enhanced unit cell with respect to file-name (for OFL with symmetry model)", ',', ',', "1,1");
  (*SystemGroup) += new BooleanOption ('\n', "decoupled", "assume that the FTI states are made of two decoupled FCI copies");
  (*SystemGroup) += new BooleanOption  ('\n', "3d", "consider a 3d model instead of a 2d model");
  (*SystemGroup) += new BooleanOption  ('\n', "Wannier", "Wannier basis");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with partent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "kya-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Ky value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIEntanglementEntropyParticlePartition -h" << endl;
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
  int* TotalKz = 0;
  int NbrParticles = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  int NbrSiteZ = 0;
  bool Statistics = true;
  double* Coefficients = 0;
  bool ShowTimeFlag = Manager.GetBoolean("show-time");
  bool Flag3d = Manager.GetBoolean("3d");
  bool FlagDecoupled = Manager.GetBoolean("decoupled");
  int TotalSpin = 0;
  int NbrBands = Manager.GetInteger("nbr-subbands");
  bool FlagWannier = Manager.GetBoolean("Wannier");
  bool EigenstateFlag = Manager.GetBoolean("density-eigenstate");
  int FilterKya = Manager.GetInteger("kya-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
  

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
      TotalKz = new int[1];
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
       TotalKz = new int[NbrSpaces];
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

  if (Flag3d == false)
    {
      NbrSiteZ = 1;
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  TotalKx[i] = 0;
	  TotalKy[i] = 0;
	  TotalKz[i] = 0;
	  double Mass = 0.0;
	  if (FlagDecoupled == false)
	    {
	      if(FlagWannier == false) 
		{
		  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[i],
									  NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], Mass, Statistics) == false)
		    {
		      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
		      return -1;
		    }
		}
	      else
		{
		  if (FQHEOnSquareLatticeWannierFindSystemInfoFromVectorFileName(GroundStateFiles[i],
										 NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], Statistics) == false)
		    {
		      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
		      return -1;
		    }
		}
	    }
	  else
	    {
	      if (FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(GroundStateFiles[i],
									      NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], TotalSpin, Statistics) == false)
		{
		  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
		  return -1;
		}
		cout << GroundStateFiles[i] << " " << NbrParticles << " " << NbrSiteX << " " << NbrSiteY << " " << NbrSiteZ << " " << TotalKx[i] << " " << TotalKy[i] << " " << TotalKz[i] << " " << TotalSpin << " " << Statistics << endl;
	    }
	}
    }
  else
    {
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  TotalKx[i] = 0;
	  TotalKy[i] = 0;
	  TotalKz[i] = 0;
	  if (FQHEOnCubicLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[i],
								 NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i], Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
// 	  cout << GroundStateFiles[i] << " " << NbrParticles << " " << NbrSiteX << " " << NbrSiteY << " " << NbrSiteZ << " " << TotalKx[i] << " " << TotalKy[i] << " " << TotalKz[i]  << endl;
	}
    }

  // decide whether to enhance the representation by symmetry factors
  int l; 
  int *SymmetryFactors = Manager.GetIntegers("enhance",l);
  if (l!=2)
    {
      cout << "Error: symmetry factors need to be given for the x- and y-direction using --enhance sx,sy"<<endl;
      exit(1);
    }
  NbrSiteX *= SymmetryFactors[0];
  NbrSiteY *= SymmetryFactors[1];

  GroundStates = new ComplexVector [NbrSpaces];  
  int TotalNbrSites;
  if(FlagWannier == false || (FlagWannier == true &&  TotalKx[0]>-1))
    TotalNbrSites = NbrSiteX * NbrSiteY * NbrSiteZ;
  else
    TotalNbrSites = NbrSiteY * NbrSiteZ;
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
      if(FlagWannier == false  || (FlagWannier == true &&  TotalKx[0]>-1) )
	TmpIndex = (((TotalKx[i] * NbrSiteY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
      else
	TmpIndex = TotalKy[i] * NbrSiteZ + TotalKz[i];
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
      if(FlagWannier == false  || (FlagWannier == true &&  TotalKx[0]>-1) )
	TmpIndex = (((TotalKx[i] * NbrSiteY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
      else
	TmpIndex = TotalKy[i] * NbrSiteZ + TotalKz[i];
      GroundStatePerMomentumSector[TmpIndex][NbrGroundStatePerMomentumSector[TmpIndex]] = GroundStates[i];
      CoefficientPerMomentumSector[TmpIndex][NbrGroundStatePerMomentumSector[TmpIndex]] = Coefficients[i];
      NbrGroundStatePerMomentumSector[TmpIndex]++;
    }  


  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      if (Flag3d == false)
	{
	  if (FlagDecoupled == false)
	    {
	      if(FlagWannier == false  || (FlagWannier == true &&  TotalKx[0]>-1) )
		DensityMatrixFile << "#  N    Kx    Ky    lambda";
	      else
		DensityMatrixFile << "#  N    Ky    lambda";
	    }
	  else
	    {
	      DensityMatrixFile << "#  N    Kx    Ky    Sz    lambda";
	    }
	}
      else
	{
	  DensityMatrixFile << "#  N    Kx    Ky    Kz    lambda";
	}
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }

  int MaxNbrSpaces; 
  if(FlagWannier == false  || (FlagWannier == true &&  TotalKx[0]>-1) )
    MaxNbrSpaces = NbrSiteX * NbrSiteY * NbrSiteZ;
  else
    MaxNbrSpaces = NbrSiteY * NbrSiteZ;
  ParticleOnSphere** Spaces = new ParticleOnSphere*[MaxNbrSpaces];
  for (int i = 0; i < MaxNbrSpaces; ++i)
    {
      Spaces[i] = 0;
    }
  for (int i = 0; i < NbrSpaces; ++i)
    {
      int TmpIndex;
      if(FlagWannier == false  || (FlagWannier == true &&  TotalKx[0]>-1) )
	TmpIndex = (((TotalKx[i] * NbrSiteY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
      else
	TmpIndex = TotalKy[i] * NbrSiteZ + TotalKz[i];
      if (Spaces[TmpIndex] == 0)
	{
	  if (Flag3d == false)
	    {
	      if (NbrBands == 1)
		{
		  if (Statistics == true)
		    Spaces[TmpIndex] = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
		  else
		    {
		      if(FlagWannier == false)
			Spaces[TmpIndex] = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
		      else
			Spaces[TmpIndex] = new BosonOnSquareLatticeWannierSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKy[i], TotalKx[i]);
		    }
		}
	      else
		{
		  if (FlagDecoupled == false)
		    {
		      if (Statistics == true)
			Spaces[TmpIndex] = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
		      else
			Spaces[TmpIndex] = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
		    }
		  else
		    {
		      if (Statistics == true)
			Spaces[TmpIndex] = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, (TotalSpin + NbrParticles) >> 1, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
		      else
			Spaces[TmpIndex] = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, (TotalSpin + NbrParticles) >> 1, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
		    }
		}
	    }
	  else
	    {
	      if (NbrBands == 1)
		{
                  if (Statistics == true)
                    Spaces[TmpIndex] = new FermionOnCubicLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i]);
                  else
                    Spaces[TmpIndex] = new BosonOnCubicLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i]);
                }
              else
                {
                  if (Statistics == true)
                    Spaces[TmpIndex] = new FermionOnCubicLatticeWithSpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i]);
                  else
                    Spaces[TmpIndex] = new BosonOnCubicLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i]);
                }
	    }
	}
      if (Spaces[TmpIndex]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	{
	      cout << "dimension mismatch between Hilbert space and ground state" << endl;
	      return 0;
	}
    }
  
  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "partent");
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
  
  int MaxSubsystemNbrParticles = (NbrParticles >> 1) + (NbrParticles & 1);
  if (Manager.GetInteger("max-na") > 0)
    MaxSubsystemNbrParticles = Manager.GetInteger("max-na");
  int SubsystemNbrParticles = Manager.GetInteger("min-na");
  
  for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      int ComplementaryNbrParticles = NbrParticles - SubsystemNbrParticles;
      int MinSz = -SubsystemNbrParticles;
      if ((TotalSpin - ComplementaryNbrParticles) > MinSz)
	MinSz = (TotalSpin - ComplementaryNbrParticles);
      int MaxSz = SubsystemNbrParticles;
      if ((TotalSpin + ComplementaryNbrParticles) < MaxSz)
	MaxSz = (TotalSpin + ComplementaryNbrParticles);
      if (FlagDecoupled == false)
	{
	  MinSz = 0;
	  MaxSz = 0;
	}
      int SubsystemTotalKxMin;
      int SubsystemTotalKxMax;
      if(FlagWannier == false  || (FlagWannier == true &&  TotalKx[0]>-1) )
	{
	  SubsystemTotalKxMin = 0;
	  SubsystemTotalKxMax = NbrSiteX;
	}
      else
	{
	  SubsystemTotalKxMin = -1;
	  SubsystemTotalKxMax = 0;
	}

      for (int SubsystemTotalSz = MinSz; SubsystemTotalSz <= MaxSz; SubsystemTotalSz += 2)
	{
	  for (int SubsystemTotalKx = SubsystemTotalKxMin; SubsystemTotalKx < SubsystemTotalKxMax; ++SubsystemTotalKx)
	    {
	      for (int SubsystemTotalKy = 0; SubsystemTotalKy < NbrSiteY; ++SubsystemTotalKy)
		{
		  for (int SubsystemTotalKz = 0; SubsystemTotalKz < NbrSiteZ; ++SubsystemTotalKz)
		    {
		      if (Flag3d == false)
			{
			  if (FlagDecoupled == false)
			  {
			    if(FlagWannier == false || (FlagWannier == true &&  TotalKx[0]>-1))
			      cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Kx=" << SubsystemTotalKx << " Ky=" << SubsystemTotalKy << endl;
			    else
			      cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " Ky=" << SubsystemTotalKy << endl;
			  }
			  else
			  {
			    cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Kx=" << SubsystemTotalKx << " Ky=" << SubsystemTotalKy << " Sz=" << SubsystemTotalSz <<endl;
			  }
			}
		      else
			cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Kx=" << SubsystemTotalKx << " Ky=" << SubsystemTotalKy << " Kz=" << SubsystemTotalKz << endl;
		      
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
		      if (Flag3d == false)
			{
			  if (NbrBands == 1)
			    {
			      if (Statistics == true)
				{
				  if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
				    {
				      PartialDensityMatrix = ((FermionOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				    }
				  else
				    {
				      PartialDensityMatrix = ((FermionOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
				    }
				}
			      else
				{
				  if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
				    {
				      if(FlagWannier == false)
					PartialDensityMatrix = ((BosonOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				      else
					PartialDensityMatrix = ((BosonOnSquareLatticeWannierSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				    }
				  else
				    {
				      if(FlagWannier == false)
					PartialDensityMatrix = ((BosonOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
				      else
					PartialDensityMatrix = ((BosonOnSquareLatticeWannierSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
				    }
				}
			    }
			  else
			    {
			      if (FlagDecoupled == false)
				{
				  if (Statistics == true)
				    {
				      if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
				      {
				      	PartialDensityMatrix = ((FermionOnSquareLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				      }
				      else
				      {
					PartialDensityMatrix = ((FermionOnSquareLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
				      }
				    }
				  else
				    {
				      if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
					{
					  PartialDensityMatrix = ((BosonOnSquareLatticeWithSU2SpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
					}
				      else
					{
					  PartialDensityMatrix = ((BosonOnSquareLatticeWithSU2SpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
					}
				    }
				}
			      else
				{
				  if (Statistics == true)
				    {
				      if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
					{
					  PartialDensityMatrix = ((FermionOnSquareLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
					}
				      else
					{
					  PartialDensityMatrix = ((FermionOnSquareLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalSz, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
					}
				    }
				  else
				    {
				      if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
					{
					  PartialDensityMatrix = ((BosonOnSquareLatticeWithSU2SpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
					}
				      else
					{
					  PartialDensityMatrix = ((BosonOnSquareLatticeWithSU2SpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalSz, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
					}
				    }
				}
			    }
			}
		      else
			{
	                  if (NbrBands == 1)
	                    {
                              if (Statistics == true)
                                {
                                  if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
                                    {
                                      PartialDensityMatrix = ((FermionOnCubicLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
                                      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
                                    }
                                  else
                                    {
                                      PartialDensityMatrix = ((FermionOnCubicLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
                                    }
                                }
                              else
                                {
                                  if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
                                    {
                                      PartialDensityMatrix = ((BosonOnCubicLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
                                      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
                                    }
                                  else
                                    {
                                      PartialDensityMatrix = ((BosonOnCubicLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
                                    }
                                }
                            }
                          else
                            {
                              if (Statistics == true)
                                {
                                  if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
                                    {
                                      PartialDensityMatrix = ((FermionOnCubicLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
                                      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
                                    }
                                  else
                                    {
                                      PartialDensityMatrix = ((FermionOnCubicLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
                                    }
                                }
                              else
                                {
                                  if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
                                    {
                                      PartialDensityMatrix = ((BosonOnCubicLatticeWithSU2SpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
                                      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
                                    }
                                  else
                                    {
                                      PartialDensityMatrix = ((BosonOnCubicLatticeWithSU2SpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
                                    }
                                }
                            }
			}
		      
		      ++TmpIndex;
		      
		      for (; TmpIndex < TotalNbrSites; ++TmpIndex)
			{
			  if (NbrGroundStatePerMomentumSector[TmpIndex] != 0)
			    {
			      if (Flag3d == false)
				{
				  if (NbrBands == 1)
				    {
				      if (Statistics == true)
					{
					  if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
					    {
					      HermitianMatrix TmpMatrix = ((FermionOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
					      PartialDensityMatrix += TmpMatrix;
					    }
					  else
					    {
					      HermitianMatrix TmpMatrix = ((FermionOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
					      PartialDensityMatrix += TmpMatrix;
					    }
					}
				      else
					{
					  if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
					    {
					      HermitianMatrix TmpMatrix;
					      if(FlagWannier == false)
						 TmpMatrix = ((BosonOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					      else
						TmpMatrix = ((BosonOnSquareLatticeWannierSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
					      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
					      PartialDensityMatrix += TmpMatrix;
					    }
					  else
					    {
					      HermitianMatrix TmpMatrix;
					      if(FlagWannier == false)
						TmpMatrix = ((BosonOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
					      else
						TmpMatrix = ((BosonOnSquareLatticeWannierSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
					      PartialDensityMatrix += TmpMatrix;
					    }
					}
				    }
				  else
				    {
				      if (FlagDecoupled == false)
					{
					  if (Statistics == true)
					    {
					      if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
						{
						  HermitianMatrix TmpMatrix = ((FermionOnSquareLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
						  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
						  PartialDensityMatrix += TmpMatrix;				      			
						}
						else
						{
						  HermitianMatrix TmpMatrix = ((FermionOnSquareLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
						  PartialDensityMatrix += TmpMatrix;				      			
						}
					    }
					  else
					    {
					      if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
						{
						  HermitianMatrix TmpMatrix = ((BosonOnSquareLatticeWithSU2SpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
						  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
						  PartialDensityMatrix += TmpMatrix;
						}
					      else
						{
						  HermitianMatrix TmpMatrix = ((BosonOnSquareLatticeWithSU2SpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
						  PartialDensityMatrix += TmpMatrix;
						}
					    }
					}
				      else
					{
					  if (Statistics == true)
					    {
					      if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
						{
						  HermitianMatrix TmpMatrix = ((FermionOnSquareLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
						  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
						  PartialDensityMatrix += TmpMatrix;
						}
					      else
						{
						  HermitianMatrix TmpMatrix = ((FermionOnSquareLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalSz, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
						  PartialDensityMatrix += TmpMatrix;
						}
					    }
					  else
					    {
					      if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
						{
						  HermitianMatrix TmpMatrix = ((BosonOnSquareLatticeWithSU2SpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
						  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
						  PartialDensityMatrix += TmpMatrix;
						}
					      else
						{
						  HermitianMatrix TmpMatrix = ((BosonOnSquareLatticeWithSU2SpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalSz, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
						  PartialDensityMatrix += TmpMatrix;
						}
					    }
					}
				    }
				}
			      else
				{
	                          if (NbrBands == 1)
	                            {
                                      if (Statistics == true)
                                        {
                                          if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
                                            {
                                              HermitianMatrix TmpMatrix = ((FermionOnCubicLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
                                              TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
                                              PartialDensityMatrix += TmpMatrix;
                                            }
                                          else
                                            {
                                              HermitianMatrix TmpMatrix = ((FermionOnCubicLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
                                              PartialDensityMatrix += TmpMatrix;
                                            }
                                        }
                                      else
                                        {
                                          if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
                                            {
                                              HermitianMatrix TmpMatrix = ((BosonOnCubicLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
                                              TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
                                              PartialDensityMatrix += TmpMatrix;
                                            }
                                          else
                                            {
                                              HermitianMatrix TmpMatrix = ((BosonOnCubicLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
                                              PartialDensityMatrix += TmpMatrix;
                                            }
                                        }
                                    }
                                  else
                                    {
                                      if (Statistics == true)
                                        {
                                          if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
                                            {
                                              HermitianMatrix TmpMatrix = ((FermionOnCubicLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
                                              TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
                                              PartialDensityMatrix += TmpMatrix;
                                            }
                                          else
                                            {
                                              HermitianMatrix TmpMatrix = ((FermionOnCubicLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
                                              PartialDensityMatrix += TmpMatrix;
                                            }
                                        }
                                      else
                                        {
                                          if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
                                            {
                                              HermitianMatrix TmpMatrix = ((BosonOnCubicLatticeWithSU2SpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
                                              TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
                                              PartialDensityMatrix += TmpMatrix;
                                            }
                                          else
                                            {
                                              HermitianMatrix TmpMatrix = ((BosonOnCubicLatticeWithSU2SpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
                                              PartialDensityMatrix += TmpMatrix;
                                            }
                                        }
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
			      if ((EigenstateFlag == true) && (FilterKya == SubsystemTotalKy))
			      {
				ComplexMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(), PartialDensityMatrix.GetNbrRow(), true);
// 				for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
// 				  TmpEigenstates[i][i] = 1.0;
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
				      sprintf(TmpSuffix, "partent_na_%d_kxa_%d_kya_%d.%d.vec", SubsystemNbrParticles, SubsystemTotalKy, SubsystemTotalKx, i);
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
			  TmpDiag.SortMatrixDownOrder();
			  if (DensityMatrixFileName != 0)
			    {
			      ofstream DensityMatrixFile;
			      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			      DensityMatrixFile.precision(14);
			      if (Flag3d == false)
				{
				  if (FlagDecoupled == false)
				    {
				      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
					DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << TmpDiag[i] << endl;
				    }
				  else
				    {
				      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
					DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << SubsystemTotalSz << " " << TmpDiag[i] << endl;
				    }
				}
			      else
				{
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << SubsystemTotalKz << " " << TmpDiag[i] << endl;
				}
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
				if (Flag3d == false)
				  {
				    if (FlagDecoupled == false)
				      {
					DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << TmpValue << endl;
				      }
				    else
				      {
					DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << SubsystemTotalSz << " " << TmpValue << endl;
				      }
				  }
				else
				  {
				    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << SubsystemTotalKz << " " << TmpValue << endl;
				  }
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
      File << SubsystemNbrParticles << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << endl;
    }
  File.close();

  return 0;
}
