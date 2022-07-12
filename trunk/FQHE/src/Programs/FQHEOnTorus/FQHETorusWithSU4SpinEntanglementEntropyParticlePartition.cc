#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include "HilbertSpace/BosonOnTorusWithSU4Spin.h"

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

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

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
  OptionManager Manager ("FQHETorusWithSU4SpinEntanglementEntropyParticlePartition" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "kymax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-ky", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-na", "minimum size of the particles whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-na", "maximum size of the particles whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with partent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "kya-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Ky value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "sza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Sz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "iza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Iz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "pza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Pz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusWithSU4SpinEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHETorusWithSU4SpinEntanglementEntropyParticlePartition -h" << endl;
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
  int FilterSza = Manager.GetInteger("sza-eigenstate");
  int FilterIza = Manager.GetInteger("iza-eigenstate");
  int FilterPza = Manager.GetInteger("pza-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
  bool ShowTimeFlag = Manager.GetBoolean("show-time");
  int* TotalKy = 0;
  int* TotalSz = 0;
  int* TotalIz = 0;
  int* TotalPz = 0;
  bool Statistics = true;
  int NbrSpaces = 1;
  ParticleOnSphereWithSU4Spin** Spaces = 0;
  RealVector* GroundStates = 0;
  char** GroundStateFiles = 0;

  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalKy = new int[1];
      TotalSz = new int[1];
      TotalIz = new int[1];
      TotalPz = new int[1];
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
       TotalKy = new int[NbrSpaces];
       TotalSz = new int[NbrSpaces];
       TotalIz = new int[NbrSpaces];
       TotalPz = new int[NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	 }
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalKy[i] = 0;
      TotalSz[i] = 0;
      TotalIz[i] = 0;
      TotalPz[i] = 0;
      if (FQHEOnTorusWithSU4SpinFindSystemInfoFromVectorFileName(GroundStateFiles[i],
								 NbrParticles, KyMax, TotalKy[i], TotalSz[i], TotalIz[i], TotalPz[i], Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
     }


  GroundStates = new RealVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
      {
	cout << "can't open vector file " << GroundStateFiles[i] << endl;
	return -1;      
      }


  Spaces = new ParticleOnSphereWithSU4Spin* [NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (Statistics == true)
	{
	  cout << "error : fermions are not yet supported" << endl;
	  return 0;
	  Spaces[i] = 0;
	}
      else
	{
	  Spaces[i] = new BosonOnTorusWithSU4Spin (NbrParticles, TotalSz[i], TotalIz[i], TotalPz[i], KyMax, TotalKy[i]);
	}
      
      if (Spaces[i]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	{
	  cout << "dimension mismatch between Hilbert space and ground state" << endl;
	  return 0;
	}
    }

  int* NbrNUpPlus = new int[NbrSpaces];
  int* NbrNUpMinus = new int[NbrSpaces];
  int* NbrNDownPlus = new int[NbrSpaces];
  int* NbrNDownMinus = new int[NbrSpaces];
  int MaxNbrNUpPlus = 0;
  int MaxNbrNUpMinus = 0;
  int MaxNbrNDownPlus = 0;
  int MaxNbrNDownMinus = 0;
  for (int i = 0; i < NbrSpaces; ++i)
    { 
      NbrNUpPlus[i] = (NbrParticles + TotalSz[i] + TotalIz[i] + TotalPz[i]);
      NbrNUpMinus[i] = (NbrParticles + TotalSz[i] - TotalIz[i] - TotalPz[i]);
      NbrNDownPlus[i] = (NbrParticles - TotalSz[i] + TotalIz[i] - TotalPz[i]);
      NbrNDownMinus[i] = (NbrParticles - TotalSz[i] - TotalIz[i] + TotalPz[i]);
      NbrNUpPlus[i] >>= 2;
      NbrNUpMinus[i] >>= 2;
      NbrNDownPlus[i] >>= 2;
      NbrNDownMinus[i] >>= 2;
      if (NbrNUpPlus[i] > MaxNbrNUpPlus)
	MaxNbrNUpPlus = NbrNUpPlus[i];
      if (NbrNUpMinus[i] > MaxNbrNUpMinus)
	MaxNbrNUpMinus = NbrNUpMinus[i];
      if (NbrNDownPlus[i] > MaxNbrNDownPlus)
	MaxNbrNDownPlus = NbrNDownPlus[i];
      if (NbrNDownMinus[i] > MaxNbrNDownMinus)
	MaxNbrNDownMinus = NbrNDownMinus[i];
    }

  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "#  N    Tz    Y    Nup    Num    Ndp    Ndm    Ky    lambda";
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
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
      for (int SubsystemNbrNUpPlus = 0; SubsystemNbrNUpPlus <= MaxNbrNUpPlus; ++SubsystemNbrNUpPlus)
	for (int SubsystemNbrNUpMinus = 0; SubsystemNbrNUpMinus <= MaxNbrNUpMinus; ++SubsystemNbrNUpMinus)
	  for (int SubsystemNbrNDownPlus = 0; SubsystemNbrNDownPlus <= MaxNbrNDownPlus; ++SubsystemNbrNDownPlus)
	    {
	      int SubsystemNbrNDownMinus = SubsystemNbrParticles - SubsystemNbrNUpPlus - SubsystemNbrNUpMinus - SubsystemNbrNDownPlus;
	      if ((SubsystemNbrNDownMinus >= 0) && (SubsystemNbrNDownMinus <= MaxNbrNDownMinus))
		{
		  int SubsystemTotalSz = (SubsystemNbrNUpPlus + SubsystemNbrNUpMinus) - (SubsystemNbrNDownPlus + SubsystemNbrNDownMinus);
		  int SubsystemTotalIz = (SubsystemNbrNUpPlus + SubsystemNbrNDownPlus) - (SubsystemNbrNUpMinus + SubsystemNbrNDownMinus);
		  int SubsystemTotalPz = (SubsystemNbrNUpPlus + SubsystemNbrNDownMinus) - (SubsystemNbrNDownPlus + SubsystemNbrNUpMinus);
		  
		  int SubsystemMaxTotalKy = KyMax - 1;
		  
		  int SubsystemTotalKy = 0; 
		  for (; SubsystemTotalKy <= SubsystemMaxTotalKy; ++SubsystemTotalKy)
		    {
		      timeval TotalStartingTime;
		      timeval TotalEndingTime;
		      if (ShowTimeFlag == true)
			{
			  gettimeofday (&(TotalStartingTime), 0);
			}
		      cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Nup=" << SubsystemNbrNUpPlus 
			   << " subsystem total Num=" << SubsystemNbrNUpMinus  << " subsystem total Ndp=" << SubsystemNbrNDownPlus <<  " subsystem total Ndm=" << SubsystemNbrNDownMinus << " subsystem total Ky=" << SubsystemTotalKy << endl;
		      RealSymmetricMatrix PartialDensityMatrix;
		      for (int i = 0; i < NbrSpaces; ++i)
			{
			  if ((SubsystemNbrNUpPlus <= NbrNUpPlus[i]) && (SubsystemNbrNUpMinus <= NbrNUpMinus[i]) && 
			      (SubsystemNbrNDownPlus <= NbrNDownPlus[i]) && (SubsystemNbrNDownMinus <= NbrNDownMinus[i]))
			    {
			      RealSymmetricMatrix TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKy, 
														       SubsystemNbrNUpPlus, SubsystemNbrNUpMinus, SubsystemNbrNDownPlus, SubsystemNbrNDownMinus, GroundStates[i], Architecture.GetArchitecture());
			    if (PartialDensityMatrix.GetNbrRow() > 0)
			      PartialDensityMatrix += TmpMatrix;
			    else
			      PartialDensityMatrix = TmpMatrix;
			    }
			}
		      if (NbrSpaces > 1)
			PartialDensityMatrix /= ((double) NbrSpaces);
		      if (ShowTimeFlag == true)
			{
			  gettimeofday (&(TotalEndingTime), 0);
			  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
						((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
			  cout << "reduced density matrix evaluated in " << Dt << "s" << endl;
			}
		      if (PartialDensityMatrix.GetNbrRow() > 1)
			{
			  RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
			  if (LapackFlag == true)
			    {
			      if ((EigenstateFlag == true) && (FilterKya == SubsystemTotalKy) && (FilterSza == SubsystemTotalSz) && 
				  (FilterIza == SubsystemTotalIz) && (FilterPza == SubsystemTotalPz))
				{
				  RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
							    PartialDensityMatrix.GetNbrRow(), true);
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    TmpEigenstates[i][i] = 1.0;
				  PartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
				  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
				  char* TmpEigenstateName = new char[512];
				  int MaxNbrEigenstates = NbrEigenstates;
				  if (NbrEigenstates == 0)
				    MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
				  for (int i = 0; i < MaxNbrEigenstates; ++i)
				    {
				      if (TmpDiag[i] > 1e-14)
					{
					  if (Statistics == true)
					    {
					      sprintf (TmpEigenstateName,
						       "fermions_torus_kysym_density_n_%d_2s_%d_ky_%d_na_%d_sza_%d_iza_%d_pza_%d_kya_%d.%d.vec",
						       NbrParticles, KyMax, TotalKy[0], 
						       SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalIz, SubsystemTotalPz, SubsystemTotalKy, i);
					    }
					  else
					    {
					      sprintf (TmpEigenstateName,
						       "bosons_torus_kysym_density_n_%d_2s_%d_ky_%d_na_%d_sza_%d_iza_%d_pza_%d_kya_%d.%d.vec",
						       NbrParticles, KyMax, TotalKy[0], 
						       SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalIz, SubsystemTotalPz, SubsystemTotalKy, i);
					    }
					  TmpEigenstates[i].WriteVector(TmpEigenstateName);
					}
				    }
				  delete[] TmpEigenstateName;
				}
			      else
				{
				  PartialDensityMatrix.LapackDiagonalize(TmpDiag);
				}
			    }
			  else
			    {
			      if ((EigenstateFlag == true) && (FilterKya == SubsystemTotalKy) && (FilterSza == SubsystemTotalSz) && 
				  (FilterIza == SubsystemTotalIz) && (FilterPza == SubsystemTotalPz))
				{
				  RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
							    PartialDensityMatrix.GetNbrRow(), true);
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    TmpEigenstates[i][i] = 1.0;
				  PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
				  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
				  char* TmpEigenstateName = new char[512];
				  int MaxNbrEigenstates = NbrEigenstates;
				  if (NbrEigenstates == 0)
				    MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
				  for (int i = 0; i < MaxNbrEigenstates; ++i)
				    {
				      if (TmpDiag[i] > 1e-14)
					{
					  if (Statistics == true)
					    {
					      sprintf (TmpEigenstateName,
						       "fermions_torus_kysym_density_n_%d_2s_%d_ky_%d_na_%d_sza_%d_iza_%d_pza_%d_kya_%d.%d.vec",
						       NbrParticles, KyMax, TotalKy[0], 
						       SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalIz, SubsystemTotalPz, SubsystemTotalKy, i);
					    }
					  else
					    {
					      sprintf (TmpEigenstateName,
						       "bosons_torus_kysym_density_n_%d_2s_%d_ky_%d_na_%d_sza_%d_iza_%d_pza_%d_kya_%d.%d.vec",
						       NbrParticles, KyMax, TotalKy[0], 
						       SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalIz, SubsystemTotalPz, SubsystemTotalKy, i);
					    }
					  TmpEigenstates[i].WriteVector(TmpEigenstateName);
					}
				    }
				  delete[] TmpEigenstateName;
				}
			      else
				{
				  PartialDensityMatrix.Diagonalize(TmpDiag);
				}
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
				DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " <<  SubsystemTotalIz << " " << SubsystemTotalPz << " " << SubsystemNbrNUpPlus << " " <<  SubsystemNbrNUpMinus << " " <<  SubsystemNbrNDownPlus << " " << SubsystemNbrNDownMinus << " " << SubsystemTotalKy << " " << TmpDiag[i] << endl;
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
				DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << SubsystemTotalIz << " " << SubsystemTotalPz << " " << SubsystemNbrNUpPlus << " " <<  SubsystemNbrNUpMinus << " " <<  SubsystemNbrNDownPlus << " " << SubsystemNbrNDownMinus << " " << SubsystemTotalKy << " " << TmpValue << endl;
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
      File << SubsystemNbrParticles << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << endl;
    }
  File.close();
}

